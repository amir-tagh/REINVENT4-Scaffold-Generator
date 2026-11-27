#!/usr/bin/env python3

"""Offline checkpoint and recovery tests for V7."""

import sys
import tempfile
import unittest

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(
    0,
    str(ROOT / "scripts" / "commercial")
)

from commercial_config import CommercialConfig
from commercial_cache import CommercialCache
from commercial_screening import CommercialScreeningEngine
from swisssimilarity_client import SwissSimilarityClient


ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"

LIBRARY = "Enamine"
METHOD = "ECFP"

RAW_RESULTS = (
    ROOT
    / "tests"
    / "commercial_screening"
    / "aspirin_enamine"
    / "raw_results.txt"
).read_text()


class MockSwissSimilarityClient:

    def __init__(
        self,
        fail_submit=False,
        fail_wait=False,
        fail_retrieve=False,
    ):

        self.fail_submit = fail_submit
        self.fail_wait = fail_wait
        self.fail_retrieve = fail_retrieve

        self.submit_calls = 0
        self.wait_calls = 0
        self.retrieve_calls = 0

        self.session_id = "999000001"

    def submit(self, smiles, library, method):

        self.submit_calls += 1

        if self.fail_submit:
            raise TimeoutError(
                "Simulated ambiguous submission timeout"
            )

        return self.session_id

    def wait(self, session_id):

        self.wait_calls += 1

        if self.fail_wait:
            raise TimeoutError(
                "Simulated interrupted screening"
            )

        return "Calculation is finished"

    def retrieve(self, session_id):

        self.retrieve_calls += 1

        if self.fail_retrieve:
            raise ConnectionError(
                "Simulated retrieval failure"
            )

        return RAW_RESULTS

    @staticmethod
    def parse_results(raw, library):

        return SwissSimilarityClient.parse_results(
            raw,
            library,
        )


class TestCheckpointRecovery(unittest.TestCase):

    def setUp(self):

        self.temp = tempfile.TemporaryDirectory()

        self.database = (
            Path(self.temp.name) / "test_cache.sqlite"
        )

        self.config = CommercialConfig()

        self.cache = CommercialCache(self.database)

    def tearDown(self):

        self.cache.close()
        self.temp.cleanup()

    def make_engine(
        self,
        client,
        allow_remote=True,
        maximum_new_jobs=1,
    ):

        return CommercialScreeningEngine(
            config=self.config,
            cache=self.cache,
            client=client,
            allow_remote=allow_remote,
            maximum_new_jobs=maximum_new_jobs,
        )

    def test_01_successful_submission(self):

        client = MockSwissSimilarityClient()

        engine = self.make_engine(client)

        raw, status, session, timestamp, error = (
            engine.retrieve_screening(
                ASPIRIN,
                LIBRARY,
            )
        )

        self.assertEqual(status, "REMOTE_SCREENING")
        self.assertEqual(client.submit_calls, 1)
        self.assertEqual(client.wait_calls, 1)
        self.assertEqual(client.retrieve_calls, 1)
        self.assertEqual(engine.new_jobs, 1)

        record = self.cache.lookup(
            ASPIRIN, LIBRARY, METHOD
        )

        self.assertEqual(record["status"], "COMPLETED")
        self.assertEqual(record["hit_count"], 400)
        self.assertEqual(record["session_id"], session)

    def test_02_completed_cache_reuse(self):

        self.cache.save_results(
            ASPIRIN,
            LIBRARY,
            METHOD,
            "999000001",
            RAW_RESULTS,
        )

        client = MockSwissSimilarityClient()
        engine = self.make_engine(client)

        result = engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(result[1], "CACHE_HIT")
        self.assertEqual(client.submit_calls, 0)
        self.assertEqual(client.wait_calls, 0)
        self.assertEqual(client.retrieve_calls, 0)
        self.assertEqual(engine.cache_hits, 1)

    def test_03_interrupted_session_recovery(self):

        client = MockSwissSimilarityClient(
            fail_wait=True
        )

        engine = self.make_engine(client)

        result = engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(
            result[1], "SCREENING_ERROR"
        )

        record = self.cache.lookup(
            ASPIRIN, LIBRARY, METHOD
        )

        self.assertEqual(record["status"], "ERROR")
        self.assertEqual(
            record["session_id"], client.session_id
        )

        # Simulate restarting the application.
        self.cache.close()
        self.cache = CommercialCache(self.database)

        recovery_client = MockSwissSimilarityClient()

        recovery_engine = self.make_engine(
            recovery_client
        )

        result = recovery_engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(
            result[1], "RESUMED_SESSION"
        )

        self.assertEqual(
            recovery_client.submit_calls, 0
        )

        self.assertEqual(
            recovery_client.retrieve_calls, 1
        )

        self.assertEqual(
            recovery_engine.resumed_jobs, 1
        )

        record = self.cache.lookup(
            ASPIRIN, LIBRARY, METHOD
        )

        self.assertEqual(record["status"], "COMPLETED")

    def test_04_ambiguous_submission_failure(self):

        client = MockSwissSimilarityClient(
            fail_submit=True
        )

        engine = self.make_engine(client)

        result = engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(
            result[1], "SUBMISSION_UNCERTAIN"
        )

        record = self.cache.lookup(
            ASPIRIN, LIBRARY, METHOD
        )

        self.assertEqual(record["status"], "ERROR")
        self.assertFalse(record["session_id"])

        recovery_client = MockSwissSimilarityClient()

        recovery_engine = self.make_engine(
            recovery_client
        )

        result = recovery_engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(
            result[1], "REVIEW_REQUIRED"
        )

        self.assertEqual(
            recovery_client.submit_calls, 0
        )

    def test_05_remote_job_limit(self):

        client = MockSwissSimilarityClient()

        engine = self.make_engine(
            client,
            maximum_new_jobs=1,
        )

        first = engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        second = engine.retrieve_screening(
            IBUPROFEN, LIBRARY
        )

        self.assertEqual(
            first[1], "REMOTE_SCREENING"
        )

        self.assertEqual(
            second[1], "JOB_LIMIT_REACHED"
        )

        self.assertEqual(client.submit_calls, 1)

        self.assertIsNone(
            self.cache.lookup(
                IBUPROFEN, LIBRARY, METHOD
            )
        )

    def test_06_offline_enforcement(self):

        client = MockSwissSimilarityClient()

        engine = self.make_engine(
            client,
            allow_remote=False,
        )

        result = engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(
            result[1], "NOT_SCREENED"
        )

        self.assertEqual(client.submit_calls, 0)
        self.assertEqual(client.wait_calls, 0)
        self.assertEqual(client.retrieve_calls, 0)

    def test_07_retrieval_failure_recovery(self):

        client = MockSwissSimilarityClient(
            fail_retrieve=True
        )

        engine = self.make_engine(client)

        result = engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(
            result[1], "SCREENING_ERROR"
        )

        record = self.cache.lookup(
            ASPIRIN, LIBRARY, METHOD
        )

        self.assertEqual(record["status"], "ERROR")
        self.assertEqual(
            record["session_id"], client.session_id
        )

        recovery_client = MockSwissSimilarityClient()

        recovery_engine = self.make_engine(
            recovery_client
        )

        result = recovery_engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(
            result[1], "RESUMED_SESSION"
        )

        self.assertEqual(
            recovery_client.submit_calls, 0
        )

        self.assertEqual(
            recovery_client.retrieve_calls, 1
        )

    def test_08_resume_disabled(self):

        self.cache.save_session(
            ASPIRIN,
            LIBRARY,
            METHOD,
            "999000001",
        )

        self.config.data["cache"][
            "resume_sessions"
        ] = False

        client = MockSwissSimilarityClient()
        engine = self.make_engine(client)

        result = engine.retrieve_screening(
            ASPIRIN, LIBRARY
        )

        self.assertEqual(
            result[1], "RESUME_DISABLED"
        )

        self.assertEqual(client.submit_calls, 0)
        self.assertEqual(client.wait_calls, 0)
        self.assertEqual(client.retrieve_calls, 0)


if __name__ == "__main__":

    unittest.main(verbosity=2)
