#!/usr/bin/env python3

import argparse
import csv
import json
import subprocess
import sys

from datetime import datetime, timezone
from pathlib import Path

from commercial_cache import CommercialCache
from commercial_config import CommercialConfig
from molecular_identity import get_identity
from swisssimilarity_client import SwissSimilarityClient


ROOT = Path(__file__).resolve().parents[2]

DEFAULT_CONFIG = ROOT / "config/commercial_sources.json"

DEFAULT_OUTPUT = (
    ROOT / "tests/commercial_screening/library_validation"
)

SCREENING_SCRIPT = (
    ROOT / "scripts/commercial/commercial_screening.py"
)

REFERENCE_MOLECULES = {
    "ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, data):
    path.write_text(
        json.dumps(data, indent=4) + "\n"
    )


def read_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def validate_library(
    library,
    smiles,
    original_config,
    database,
    output,
    allow_remote,
    remaining_jobs,
):

    library_dir = output / library
    library_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(json.dumps(original_config))

    config["libraries"][library]["enabled"] = True

    config["profiles"]["validation_test"] = [library]

    config["api"]["allow_remote_submission"] = allow_remote
    config["api"]["maximum_new_jobs_per_run"] = min(
        remaining_jobs, 1
    )
    config["api"]["parallel_jobs"] = 1

    config["cache"]["database"] = str(database)

    config_path = library_dir / "validation_config.json"

    write_json(config_path, config)

    result_dir = library_dir / "results"

    result = {
        "Library": library,
        "Reference_SMILES": smiles,
        "Validation_Status": "NOT_VALIDATED",
        "Screening_Status": "",
        "Session_ID": "",
        "Returned_Hits": 0,
        "Exact_Matches": 0,
        "Exact_Reported_Structure": 0,
        "Reported_Analogs": 0,
        "Invalid_Structures": 0,
        "New_Submission_Attempts": 0,
        "Cache_Verified": False,
        "Error": "",
        "Timestamp": utc_now(),
    }

    command = [
        sys.executable,
        "-u",
        str(SCREENING_SCRIPT),
        "--smiles", smiles,
        "--config", str(config_path),
        "--profile", "validation_test",
        "--mode", "both",
        "--minimum-similarity", "0.5",
        "--maximum-analogs", "20",
        "--output", str(result_dir),
    ]

    if allow_remote:
        command.extend([
            "--allow-remote",
            "--maximum-new-jobs",
            str(min(remaining_jobs, 1)),
        ])

    with CommercialCache(database) as cache:
        before = cache.lookup(
            smiles, library, "ECFP"
        )

    try:

        process = subprocess.run(
            command,
            cwd=ROOT,
            check=False,
        )

        if process.returncode != 0:
            raise RuntimeError(
                f"Screening exited with code "
                f"{process.returncode}"
            )

        summary = json.loads(
            (result_dir / "screening_summary.json").read_text()
        )

        result["New_Submission_Attempts"] = int(
            summary["New_Submission_Attempts"]
        )

        jobs = read_csv(
            result_dir / "screening_jobs.csv"
        )

        if len(jobs) != 1:
            raise RuntimeError(
                f"Expected one job, found {len(jobs)}"
            )

        job = jobs[0]

        if job["Library"] != library:
            raise RuntimeError("Incorrect library")

        expected_identity = get_identity(smiles)

        observed_identity = get_identity(
            job["Query_SMILES"]
        )

        if expected_identity != observed_identity:
            raise RuntimeError("Incorrect query identity")

        result["Screening_Status"] = job["Status"]
        result["Session_ID"] = job["Session_ID"]

        for key in [
            "Returned_Hits",
            "Exact_Matches",
            "Exact_Reported_Structure",
            "Reported_Analogs",
            "Invalid_Structures",
        ]:
            result[key] = int(job[key])

        result["Error"] = job["Error"]

        successful = {
            "CACHE_HIT",
            "REMOTE_SCREENING",
            "RESUMED_SESSION",
        }

        if job["Status"] not in successful:
            return result

        if not job["Session_ID"]:
            raise RuntimeError("Missing session ID")

        if result["Returned_Hits"] <= 0:
            raise RuntimeError("No retrieved catalog entries")

        if result["Invalid_Structures"] != 0:
            raise RuntimeError("Invalid catalog structures")

        if result["Reported_Analogs"] > 20:
            raise RuntimeError("Analog limit exceeded")

        with CommercialCache(database) as cache:

            record = cache.get_completed(
                smiles, library, "ECFP"
            )

            if record is None:
                raise RuntimeError(
                    "Completed cache record missing"
                )

            if record["session_id"] != job["Session_ID"]:
                raise RuntimeError(
                    "Cache session ID mismatch"
                )

            hits = SwissSimilarityClient.parse_results(
                record["raw_results"], library
            )

            if len(hits) != result["Returned_Hits"]:
                raise RuntimeError(
                    "Cached hit count mismatch"
                )

            if record["hit_count"] != len(hits):
                raise RuntimeError(
                    "Stored hit count mismatch"
                )

        result["Cache_Verified"] = True
        result["Validation_Status"] = "PASSED"

    except Exception as exc:

        result["Validation_Status"] = "FAILED"
        result["Error"] = str(exc)

    finally:

        # Conservatively account for a possible submission
        # if the screening subprocess failed before writing
        # its summary.

        with CommercialCache(database) as cache:
            after = cache.lookup(
                smiles, library, "ECFP"
            )

        if (
            allow_remote
            and remaining_jobs > 0
            and before is None
            and after is not None
        ):
            result["New_Submission_Attempts"] = max(
                result["New_Submission_Attempts"], 1
            )

        write_json(
            library_dir / "validation_result.json",
            result,
        )

    return result


def main():

    parser = argparse.ArgumentParser(
        description="Validate commercial screening libraries"
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
    )

    parser.add_argument(
        "--libraries",
        nargs="+",
    )

    parser.add_argument(
        "--reference",
        choices=REFERENCE_MOLECULES,
        default="ibuprofen",
    )

    parser.add_argument(
        "--allow-remote",
        action="store_true",
    )

    parser.add_argument(
        "--maximum-new-jobs",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
    )

    args = parser.parse_args()

    if args.maximum_new_jobs < 0:
        parser.error("Job limit cannot be negative")

    if args.maximum_new_jobs and not args.allow_remote:
        parser.error(
            "Remote submission requires --allow-remote"
        )

    if args.allow_remote and args.maximum_new_jobs == 0:
        parser.error(
            "Specify a positive --maximum-new-jobs limit"
        )

    config_path = args.config.resolve()

    original_bytes = config_path.read_bytes()

    config = CommercialConfig(config_path)

    original_config = json.loads(original_bytes)

    database = config.database_path()

    libraries = args.libraries or [
        name
        for name, settings in config.libraries.items()
        if not settings["validated"]
    ]

    libraries = list(dict.fromkeys(libraries))

    unknown = [
        name for name in libraries
        if name not in config.libraries
    ]

    if unknown:
        parser.error(f"Unknown libraries: {unknown}")

    output = args.output.resolve()

    output.mkdir(parents=True, exist_ok=True)

    smiles = REFERENCE_MOLECULES[args.reference]

    print("\nCOMMERCIAL LIBRARY VALIDATION")
    print("Reference:", args.reference)
    print("Libraries:", libraries)
    print("Remote allowed:", args.allow_remote)
    print("Maximum new jobs:", args.maximum_new_jobs)
    print("Cache:", database)

    results = []

    remaining = args.maximum_new_jobs

    for index, library in enumerate(libraries, 1):

        print(
            f"\n[{index}/{len(libraries)}] "
            f"Validating {library}"
        )

        result = validate_library(
            library=library,
            smiles=smiles,
            original_config=original_config,
            database=database,
            output=output,
            allow_remote=(
                args.allow_remote and remaining > 0
            ),
            remaining_jobs=remaining,
        )

        remaining -= result["New_Submission_Attempts"]

        results.append(result)

        print(
            library,
            result["Validation_Status"],
            result["Screening_Status"],
            f"hits={result['Returned_Hits']}",
        )

        if result["Validation_Status"] == "FAILED":
            print("Stopping after validation failure")
            break

        if result["Screening_Status"] in {
            "SUBMISSION_UNCERTAIN",
            "SCREENING_ERROR",
            "REVIEW_REQUIRED",
            "CLASSIFICATION_ERROR",
        }:
            print("Stopping for manual review")
            break

    if config_path.read_bytes() != original_bytes:
        raise RuntimeError(
            "Original configuration was modified"
        )

    passed = [
        row["Library"] for row in results
        if row["Validation_Status"] == "PASSED"
    ]

    summary = {
        "Reference": args.reference,
        "Requested_Libraries": libraries,
        "Attempted_Libraries": len(results),
        "Passed_Libraries": passed,
        "Passed_Count": len(passed),
        "New_Submission_Attempts": sum(
            row["New_Submission_Attempts"]
            for row in results
        ),
        "Original_Configuration_Preserved": True,
        "Production_Configuration_Updated": False,
        "Stock_Status": "UNVERIFIED",
        "Results": results,
    }

    write_json(
        output / "validation_summary.json",
        summary,
    )

    if results:

        with (
            output / "validation_results.csv"
        ).open("w", newline="") as handle:

            writer = csv.DictWriter(
                handle,
                fieldnames=list(results[0]),
            )

            writer.writeheader()
            writer.writerows(results)

    print("\nVALIDATION SUMMARY")
    print(json.dumps(summary, indent=4))


if __name__ == "__main__":
    main()
