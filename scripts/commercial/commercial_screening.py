#!/usr/bin/env python3

"""Integrated commercial molecular screening engine."""

import argparse
import csv
import json

from pathlib import Path
from collections import Counter

from commercial_config import CommercialConfig
from commercial_cache import CommercialCache
from molecular_identity import get_identity, classify_match
from swisssimilarity_client import SwissSimilarityClient


def write_csv(filename, records, fields):

    with Path(filename).open("w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(records)


def load_queries(smiles=None, input_file=None,
                 smiles_column="SMILES"):

    records = []

    if smiles is not None:

        records.append(("Query_000001", smiles))

    else:

        path = Path(input_file)

        if path.suffix.lower() == ".csv":

            with path.open(newline="") as handle:

                reader = csv.DictReader(handle)

                if smiles_column not in (
                    reader.fieldnames or []
                ):
                    raise ValueError(
                        f"Missing column: {smiles_column}"
                    )

                for i, row in enumerate(reader, 1):

                    records.append((
                        f"Query_{i:06d}",
                        row[smiles_column],
                    ))

        else:

            with path.open() as handle:

                for i, line in enumerate(handle, 1):

                    line = line.strip()

                    if not line or line.startswith("#"):
                        continue

                    fields = line.split(maxsplit=1)

                    query_id = (
                        fields[1]
                        if len(fields) > 1
                        else f"Query_{i:06d}"
                    )

                    records.append((
                        query_id,
                        fields[0],
                    ))

    unique = {}
    invalid = []

    for query_id, raw_smiles in records:

        try:

            identity = get_identity(raw_smiles)

        except (ValueError, TypeError) as error:

            invalid.append({
                "Query_ID": query_id,
                "SMILES": raw_smiles,
                "Error": str(error),
            })

            continue

        key = identity.canonical_smiles

        if key not in unique:

            unique[key] = {
                "Query_ID": query_id,
                "SMILES": key,
                "InChIKey": identity.inchikey,
                "Identity": identity,
                "Input_Count": 1,
                "Input_IDs": [query_id],
            }

        else:

            unique[key]["Input_Count"] += 1
            unique[key]["Input_IDs"].append(query_id)

    return list(unique.values()), invalid, len(records)


class CommercialScreeningEngine:

    def __init__(
        self,
        config,
        cache,
        client,
        mode=None,
        minimum_similarity=None,
        maximum_analogs=None,
        allow_remote=False,
        maximum_new_jobs=1,
    ):

        self.config = config
        self.cache = cache
        self.client = client

        self.method = config.screening["method"]

        self.mode = (
            mode or config.screening["mode"]
        )

        self.minimum_similarity = (
            minimum_similarity
            if minimum_similarity is not None
            else config.screening[
                "minimum_analog_similarity"
            ]
        )

        self.maximum_analogs = (
            maximum_analogs
            if maximum_analogs is not None
            else config.screening[
                "maximum_analogs_per_library"
            ]
        )

        self.allow_remote = allow_remote
        self.maximum_new_jobs = maximum_new_jobs

        self.new_jobs = 0
        self.cache_hits = 0
        self.resumed_jobs = 0

        if not 0 <= self.minimum_similarity <= 1:
            raise ValueError("Invalid similarity threshold")

        if self.maximum_analogs < 0:
            raise ValueError("Invalid analog limit")

        if self.maximum_new_jobs < 0:
            raise ValueError("Invalid new-job limit")

    def retrieve_screening(self, smiles, library):

        record = self.cache.lookup(
            smiles,
            library,
            self.method,
        )

        if record and record["status"] == "COMPLETED":

            self.cache_hits += 1

            return (
                record["raw_results"],
                "CACHE_HIT",
                record["session_id"],
                record["updated_at"],
                "",
            )

        if not self.allow_remote:

            return (
                None,
                "NOT_SCREENED",
                "",
                "",
                "Remote access disabled",
            )

        session_id = None

        if record:

            session_id = record["session_id"]

            if not session_id:

                return (
                    None,
                    "REVIEW_REQUIRED",
                    "",
                    "",
                    "Existing record has no recoverable session",
                )

            if not self.config.cache["resume_sessions"]:

                return (
                    None,
                    "RESUME_DISABLED",
                    session_id,
                    "",
                    "",
                )

            source = "RESUMED_SESSION"
            self.resumed_jobs += 1

        else:

            if self.new_jobs >= self.maximum_new_jobs:

                return (
                    None,
                    "JOB_LIMIT_REACHED",
                    "",
                    "",
                    "",
                )

            self.new_jobs += 1
            source = "REMOTE_SCREENING"

            # Record submission intent before contacting
            # the remote service. This prevents automatic
            # resubmission after an ambiguous failure.
            self.cache.save_session(
                smiles,
                library,
                self.method,
                "",
            )

            try:

                session_id = self.client.submit(
                    smiles,
                    library,
                    self.method,
                )

                self.cache.save_session(
                    smiles,
                    library,
                    self.method,
                    session_id,
                )

            except Exception as error:

                self.cache.save_error(
                    smiles,
                    library,
                    self.method,
                    str(error),
                )

                return (
                    None,
                    "SUBMISSION_UNCERTAIN",
                    session_id or "",
                    "",
                    str(error),
                )

        try:

            self.client.wait(session_id)

            raw = self.client.retrieve(session_id)

            self.cache.save_results(
                smiles,
                library,
                self.method,
                session_id,
                raw,
            )

            completed = self.cache.get_completed(
                smiles,
                library,
                self.method,
            )

            return (
                raw,
                source,
                session_id,
                completed["updated_at"],
                "",
            )

        except Exception as error:

            self.cache.save_error(
                smiles,
                library,
                self.method,
                str(error),
            )

            return (
                None,
                "SCREENING_ERROR",
                session_id,
                "",
                str(error),
            )

    def classify_results(self, query, library, raw):

        hits = self.client.parse_results(raw, library)

        records = []
        counts = Counter()

        for hit in hits:

            try:

                candidate = get_identity(hit.smiles)

                match_type = classify_match(
                    query["Identity"],
                    candidate,
                    hit.similarity,
                )

            except ValueError:

                counts["INVALID_STRUCTURE"] += 1
                continue

            counts[match_type] += 1

            records.append({
                "Query_ID": query["Query_ID"],
                "Query_SMILES": query["SMILES"],
                "Query_InChIKey": query["InChIKey"],
                "Library": library,
                "Catalog_ID": hit.catalog_id,
                "Similarity": hit.similarity,
                "Matched_SMILES": hit.smiles,
                "Matched_Canonical_SMILES":
                    candidate.canonical_smiles,
                "Matched_InChIKey": candidate.inchikey,
                "Match_Type": match_type,
                "Availability_Status": "UNVERIFIED",
            })

        exact = [
            row for row in records
            if row["Match_Type"] in {
                "EXACT_MATCH",
                "EXACT_REPORTED_STRUCTURE",
            }
        ]

        analogs = [
            row for row in records
            if row["Match_Type"] == "STRUCTURAL_ANALOG"
            and row["Similarity"] >= self.minimum_similarity
        ]

        analogs.sort(
            key=lambda row: (
                -row["Similarity"],
                row["Catalog_ID"],
            )
        )

        analogs = analogs[:self.maximum_analogs]

        special = [
            row for row in records
            if row["Match_Type"] in {
                "STEREOCHEMISTRY_MISMATCH",
                "FINGERPRINT_ONLY_MATCH",
            }
        ]

        return exact, analogs, special, counts, len(hits)

    def screen(self, queries, libraries, output_dir):

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        exact_records = []
        analog_records = []
        special_records = []
        job_records = []
        molecule_records = []

        for index, query in enumerate(queries, 1):

            print(
                f"\n[{index}/{len(queries)}] "
                f"{query['Query_ID']} {query['SMILES']}",
                flush=True,
            )

            query_exact = []
            query_analogs = []
            query_special = []

            completed_libraries = 0

            for library in libraries:

                raw, status, session_id, timestamp, error = (
                    self.retrieve_screening(
                        query["SMILES"],
                        library,
                    )
                )

                job = {
                    "Query_ID": query["Query_ID"],
                    "Query_SMILES": query["SMILES"],
                    "Library": library,
                    "Method": self.method,
                    "Status": status,
                    "Session_ID": session_id,
                    "Cache_Updated_At": timestamp,
                    "Returned_Hits": 0,
                    "Exact_Matches": 0,
                    "Exact_Reported_Structure": 0,
                    "Reported_Analogs": 0,
                    "Invalid_Structures": 0,
                    "Error": error,
                }

                if raw is not None:

                    try:

                        exact, analogs, special, counts, total = (
                            self.classify_results(
                                query,
                                library,
                                raw,
                            )
                        )

                        completed_libraries += 1

                        job.update({
                            "Returned_Hits": total,
                            "Exact_Matches":
                                counts["EXACT_MATCH"],
                            "Exact_Reported_Structure":
                                counts["EXACT_REPORTED_STRUCTURE"],
                            "Reported_Analogs": len(analogs),
                            "Invalid_Structures":
                                counts["INVALID_STRUCTURE"],
                        })

                        query_exact.extend(exact)
                        query_analogs.extend(analogs)
                        query_special.extend(special)

                    except Exception as exc:

                        job["Status"] = "CLASSIFICATION_ERROR"
                        job["Error"] = str(exc)

                job_records.append(job)

                print(
                    f"  {library:20s} "
                    f"{job['Status']:24s} "
                    f"hits={job['Returned_Hits']} "
                    f"exact={job['Exact_Matches']} "
                    f"reported_exact="
                    f"{job['Exact_Reported_Structure']}",
                    flush=True,
                )

            strict_exact = [
                row for row in query_exact
                if row["Match_Type"] == "EXACT_MATCH"
            ]

            reported_exact = [
                row for row in query_exact
                if row["Match_Type"]
                == "EXACT_REPORTED_STRUCTURE"
            ]

            if strict_exact:

                molecular_status = "EXACT_CATALOG_MATCH"

            elif reported_exact:

                molecular_status = "EXACT_REPORTED_STRUCTURE"

            elif completed_libraries == 0:

                molecular_status = "NOT_SCREENED"

            elif completed_libraries < len(libraries):

                molecular_status = "PARTIALLY_SCREENED"

            else:

                molecular_status = "NO_EXACT_MATCH_RETURNED"

            molecule_records.append({
                "Query_ID": query["Query_ID"],
                "SMILES": query["SMILES"],
                "InChIKey": query["InChIKey"],
                "Input_Count": query["Input_Count"],
                "Input_IDs": json.dumps(query["Input_IDs"]),
                "Requested_Libraries": len(libraries),
                "Completed_Libraries": completed_libraries,
                "Screening_Status": molecular_status,
                "Exact_Matches": len(strict_exact),
                "Exact_Reported_Structure": len(reported_exact),
                "Reported_Analogs": len(query_analogs),
                "Special_Matches": len(query_special),
                "Exact_Match_Libraries": ";".join(sorted({
                    row["Library"]
                    for row in query_exact
                })),
                "Best_Analog_Similarity": max(
                    (
                        row["Similarity"]
                        for row in query_analogs
                    ),
                    default="",
                ),
                "Availability_Status": "UNVERIFIED",
            })

            exact_records.extend(query_exact)
            analog_records.extend(query_analogs)
            special_records.extend(query_special)

        match_fields = [
            "Query_ID",
            "Query_SMILES",
            "Query_InChIKey",
            "Library",
            "Catalog_ID",
            "Similarity",
            "Matched_SMILES",
            "Matched_Canonical_SMILES",
            "Matched_InChIKey",
            "Match_Type",
            "Availability_Status",
        ]

        write_csv(
            output_dir / "exact_matches.csv",
            exact_records
            if self.mode in {"exact", "both"} else [],
            match_fields,
        )

        write_csv(
            output_dir / "commercial_analogs.csv",
            analog_records
            if self.mode in {"similarity", "both"} else [],
            match_fields,
        )

        write_csv(
            output_dir / "special_matches.csv",
            special_records,
            match_fields,
        )

        write_csv(
            output_dir / "Commercial_Availability.csv",
            molecule_records,
            list(molecule_records[0]),
        )

        write_csv(
            output_dir / "screening_jobs.csv",
            job_records,
            list(job_records[0]),
        )

        summary = {
            "Unique_Query_Molecules": len(queries),
            "Libraries": libraries,
            "Method": self.method,
            "Mode": self.mode,
            "Minimum_Analog_Similarity":
                self.minimum_similarity,
            "Maximum_Analogs_Per_Library":
                self.maximum_analogs,
            "Requested_Screening_Jobs": len(job_records),
            "Cache_Hits": self.cache_hits,
            "New_Submission_Attempts": self.new_jobs,
            "Resumed_Jobs": self.resumed_jobs,
            "Exact_Matches": sum(
                row["Match_Type"] == "EXACT_MATCH"
                for row in exact_records
            ),
            "Exact_Reported_Structure": sum(
                row["Match_Type"]
                == "EXACT_REPORTED_STRUCTURE"
                for row in exact_records
            ),
            "Reported_Analogs": len(analog_records),
            "Special_Matches": len(special_records),
            "Job_Status_Counts": dict(Counter(
                row["Status"] for row in job_records
            )),
            "Molecular_Status_Counts": dict(Counter(
                row["Screening_Status"]
                for row in molecule_records
            )),
            "Availability_Status": "UNVERIFIED",
        }

        return summary


def main():

    parser = argparse.ArgumentParser(
        description="Commercial molecular screening engine"
    )

    inputs = parser.add_mutually_exclusive_group(
        required=True
    )

    inputs.add_argument("--smiles")
    inputs.add_argument("--input")

    parser.add_argument(
        "--smiles-column",
        default="SMILES",
    )

    parser.add_argument("--config")
    parser.add_argument("--profile")

    parser.add_argument(
        "--libraries",
        nargs="+",
    )

    parser.add_argument(
        "--mode",
        choices=["exact", "similarity", "both"],
    )

    parser.add_argument(
        "--minimum-similarity",
        type=float,
    )

    parser.add_argument(
        "--maximum-analogs",
        type=int,
    )

    parser.add_argument(
        "--allow-remote",
        action="store_true",
    )

    parser.add_argument(
        "--maximum-new-jobs",
        type=int,
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    args = parser.parse_args()

    config = (
        CommercialConfig(args.config)
        if args.config
        else CommercialConfig()
    )

    libraries = (
        args.libraries
        or config.get_libraries(args.profile)
    )

    libraries = list(dict.fromkeys(libraries))

    if not libraries:
        parser.error("No enabled libraries selected")

    for library in libraries:

        if library not in config.libraries:
            parser.error(f"Unknown library: {library}")

        if not config.libraries[library]["enabled"]:
            parser.error(f"Disabled library: {library}")

    allow_remote = args.allow_remote

    if allow_remote and not config.api[
        "allow_remote_submission"
    ]:
        parser.error(
            "Remote submission disabled in configuration"
        )

    maximum_new_jobs = (
        args.maximum_new_jobs
        if args.maximum_new_jobs is not None
        else config.api["maximum_new_jobs_per_run"]
    )

    if maximum_new_jobs > config.api[
        "maximum_new_jobs_per_run"
    ]:
        parser.error("Requested job limit exceeds configuration")

    if config.cache["refresh"]:
        parser.error("Cache refresh is not yet supported")

    if not config.cache["enabled"]:
        parser.error("Persistent cache is required")

    if not config.cache["reuse_completed"]:
        parser.error("Completed-cache reuse is required")

    queries, invalid, input_count = load_queries(
        smiles=args.smiles,
        input_file=args.input,
        smiles_column=args.smiles_column,
    )

    if not queries:
        parser.error("No valid input molecules")

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    with CommercialCache(
        config.database_path()
    ) as cache:

        with SwissSimilarityClient(
            base_url=config.base_url,
            request_timeout=config.api["request_timeout"],
            poll_interval=config.api["poll_interval"],
            max_wait=config.api["maximum_wait"],
            request_interval=config.api["request_interval"],
        ) as client:

            engine = CommercialScreeningEngine(
                config=config,
                cache=cache,
                client=client,
                mode=args.mode,
                minimum_similarity=args.minimum_similarity,
                maximum_analogs=args.maximum_analogs,
                allow_remote=allow_remote,
                maximum_new_jobs=maximum_new_jobs,
            )

            summary = engine.screen(
                queries,
                libraries,
                output,
            )

    write_csv(
        output / "invalid_queries.csv",
        invalid,
        ["Query_ID", "SMILES", "Error"],
    )

    summary["Input_Molecules"] = input_count
    summary["Invalid_Queries"] = len(invalid)
    summary["Duplicate_Queries_Removed"] = (
        input_count - len(invalid) - len(queries)
    )

    (output / "screening_summary.json").write_text(
        json.dumps(summary, indent=4) + "\n"
    )

    print("\n" + "=" * 60)
    print("COMMERCIAL SCREENING SUMMARY")
    print("=" * 60)
    print(json.dumps(summary, indent=4))
    print("\nOutput directory:", output)


if __name__ == "__main__":
    main()
