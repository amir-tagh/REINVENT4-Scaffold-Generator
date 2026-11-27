#!/usr/bin/env python3

"""SwissSimilarity API client for commercial molecular screening."""

import argparse
import csv
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path

import requests


BASE_URL = "https://swisssimilarity.ch:1234"


class SwissSimilarityError(RuntimeError):
    """SwissSimilarity API or response error."""


@dataclass(frozen=True)
class ScreeningHit:
    catalog_id: str
    similarity: float
    smiles: str
    library: str


class SwissSimilarityClient:

    def __init__(
        self,
        base_url=BASE_URL,
        request_timeout=60,
        poll_interval=15,
        max_wait=900,
        request_interval=3,
    ):
        self.base_url = base_url.rstrip("/")
        self.request_timeout = request_timeout
        self.poll_interval = poll_interval
        self.max_wait = max_wait
        self.request_interval = request_interval

        self.session = requests.Session()
        self.last_request = 0.0

        self.session.headers.update({
            "User-Agent": "MolecularCommercialScreening/0.1"
        })

    def close(self):
        self.session.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _request(self, method, endpoint, **kwargs):

        elapsed = time.monotonic() - self.last_request

        if elapsed < self.request_interval:
            time.sleep(self.request_interval - elapsed)

        try:
            response = self.session.request(
                method,
                f"{self.base_url}/{endpoint}",
                timeout=(10, self.request_timeout),
                **kwargs,
            )

            response.raise_for_status()

        except requests.RequestException as error:
            raise SwissSimilarityError(
                f"{endpoint}: {error}"
            ) from error

        finally:
            self.last_request = time.monotonic()

        return response.text.strip()

    def hello(self):
        return self._request("GET", "hello")

    def submit(self, smiles, library="Enamine", method="ECFP"):

        smiles = smiles.strip()

        if not smiles or any(c.isspace() for c in smiles):
            raise ValueError("Expected one nonempty SMILES string")

        response = self._request(
            "POST",
            "startscreen",
            params={
                "library": library,
                "method": method,
            },
            files={
                "mySMILES": (
                    "query.smi",
                    smiles + "\n",
                    "text/plain",
                )
            },
        )

        if not re.fullmatch(r"\d+", response):
            raise SwissSimilarityError(
                f"Invalid session response: {response[:500]}"
            )

        return response

    def status(self, session_id):

        return self._request(
            "GET",
            "checksession",
            params={"sessionNumber": session_id},
        )

    def wait(self, session_id):

        start = time.monotonic()

        while time.monotonic() - start < self.max_wait:

            status = self.status(session_id)
            normalized = status.lower()

            print(
                f"Session {session_id}: {status}",
                flush=True,
            )

            if any(word in normalized for word in (
                "error",
                "failed",
                "cancelled",
                "canceled",
                "not found",
                "does not exist",
            )):
                raise SwissSimilarityError(status)

            if (
                "calculation is finished" in normalized
                or normalized == "finished"
            ):
                return status

            time.sleep(self.poll_interval)

        raise TimeoutError(
            f"Session {session_id} exceeded "
            f"{self.max_wait} seconds. "
            "The remote job may still be running."
        )

    def retrieve(self, session_id):

        return self._request(
            "GET",
            "retrievesession",
            params={"sessionNumber": session_id},
        )

    @staticmethod
    def parse_results(raw_text, library):

        hits = []

        if not raw_text.strip():
            raise SwissSimilarityError(
                "Empty result response; inspect session status"
            )

        for line_number, line in enumerate(
            raw_text.splitlines(), start=1
        ):

            line = line.strip()

            if not line:
                continue

            fields = line.split()

            if len(fields) != 3:
                raise SwissSimilarityError(
                    f"Unexpected result at line {line_number}: "
                    f"{line[:200]}"
                )

            catalog_id, score, smiles = fields

            try:
                similarity = float(score)
            except ValueError as error:
                raise SwissSimilarityError(
                    f"Invalid similarity score: {score}"
                ) from error

            if (
                not math.isfinite(similarity)
                or not 0.0 <= similarity <= 1.0
            ):
                raise SwissSimilarityError(
                    f"Invalid similarity value: {similarity}"
                )

            hits.append(
                ScreeningHit(
                    catalog_id=catalog_id,
                    similarity=similarity,
                    smiles=smiles,
                    library=library,
                )
            )

        return sorted(
            hits,
            key=lambda hit: (
                -hit.similarity,
                hit.catalog_id,
            ),
        )

    @staticmethod
    def save_csv(hits, output_file):

        output_file = Path(output_file)
        output_file.parent.mkdir(
            parents=True, exist_ok=True
        )

        with output_file.open("w", newline="") as handle:

            writer = csv.writer(handle)

            writer.writerow([
                "Library",
                "Catalog_ID",
                "Similarity",
                "Matched_SMILES",
            ])

            for hit in hits:
                writer.writerow([
                    hit.library,
                    hit.catalog_id,
                    hit.similarity,
                    hit.smiles,
                ])


def main():

    parser = argparse.ArgumentParser(
        description="SwissSimilarity commercial screening client"
    )

    query = parser.add_mutually_exclusive_group(required=True)

    query.add_argument("--smiles")
    query.add_argument("--session-id")

    parser.add_argument("--library", default="Enamine")
    parser.add_argument("--method", default="ECFP")
    parser.add_argument("--output", required=True)
    parser.add_argument("--poll-interval", type=float, default=15)
    parser.add_argument("--max-wait", type=float, default=900)

    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    with SwissSimilarityClient(
        poll_interval=args.poll_interval,
        max_wait=args.max_wait,
    ) as client:

        if args.session_id:

            session_id = args.session_id

            print("Reusing session:", session_id)

        else:

            print("Submitting molecule:", args.smiles)
            print("Library:", args.library)
            print("Method:", args.method)

            session_id = client.submit(
                args.smiles,
                library=args.library,
                method=args.method,
            )

        (output / "session_id.txt").write_text(
            session_id + "\n"
        )

        print("Session ID:", session_id, flush=True)

        client.wait(session_id)

        raw = client.retrieve(session_id)

        (output / "raw_results.txt").write_text(raw + "\n")

        hits = client.parse_results(raw, args.library)

        client.save_csv(
            hits,
            output / "screening_results.csv",
        )

        perfect = sum(
            hit.similarity == 1.0 for hit in hits
        )

        print()
        print("=" * 60)
        print("SWISSSIMILARITY SCREENING COMPLETE")
        print("=" * 60)
        print("Library:", args.library)
        print("Session:", session_id)
        print("Returned compounds:", len(hits))
        print("Perfect fingerprint matches:", perfect)
        print("Output directory:", output)
        print("=" * 60)


if __name__ == "__main__":
    main()
