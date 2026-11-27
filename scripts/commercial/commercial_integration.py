#!/usr/bin/env python3

"""Integrate commercial screening results with CMOS V2."""

import argparse
import json
from pathlib import Path
from collections import Counter

import pandas as pd

from molecular_identity import get_identity


def read_csv(path):

    return pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
    )


def canonicalize(smiles):

    try:
        return get_identity(smiles).canonical_smiles

    except (ValueError, TypeError):
        return None


def integrate(input_file, screening_dir, output_dir):

    input_file = Path(input_file)
    screening_dir = Path(screening_dir)
    output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    original = read_csv(input_file)

    required = {
        "ParentID",
        "GeneratedID",
        "SMILES",
        "CMOS_V2_Intrinsic",
        "CMOS_V2_Cumulative",
        "CMOS_V2_SelectionRank",
    }

    missing = required - set(original.columns)

    if missing:
        raise ValueError(
            f"Missing CMOS columns: {missing}"
        )

    availability = read_csv(
        screening_dir / "Commercial_Availability.csv"
    )

    exact = read_csv(
        screening_dir / "exact_matches.csv"
    )

    analogs = read_csv(
        screening_dir / "commercial_analogs.csv"
    )

    jobs = read_csv(
        screening_dir / "screening_jobs.csv"
    )

    screening_summary = json.loads(
        (
            screening_dir / "screening_summary.json"
        ).read_text()
    )

    if availability["SMILES"].duplicated().any():
        raise ValueError(
            "Duplicate commercial summary identities"
        )

    if screening_summary["Mode"] != "both":
        raise ValueError(
            "Integration requires screening mode 'both'"
        )

    lookup = availability.set_index(
        "SMILES"
    ).to_dict("index")

    exact_groups = {
        key: group.to_dict("records")
        for key, group in exact.groupby("Query_SMILES")
    }

    analog_groups = {
        key: group.to_dict("records")
        for key, group in analogs.groupby("Query_SMILES")
    }

    job_groups = {
        key: group.to_dict("records")
        for key, group in jobs.groupby("Query_SMILES")
    }

    annotations = []

    for _, row in original.iterrows():

        canonical = canonicalize(row["SMILES"])

        record = lookup.get(canonical)

        matches = exact_groups.get(canonical, [])
        related = analog_groups.get(canonical, [])
        query_jobs = job_groups.get(canonical, [])

        strict = [
            item for item in matches
            if item["Match_Type"] == "EXACT_MATCH"
        ]

        reported = [
            item for item in matches
            if item["Match_Type"]
            == "EXACT_REPORTED_STRUCTURE"
        ]

        vendor_ids = sorted({
            f"{item['Library']}:{item['Catalog_ID']}"
            for item in matches
        })

        related.sort(
            key=lambda item: (
                -float(item["Similarity"]),
                item["Library"],
                item["Catalog_ID"],
            )
        )

        best = related[0] if related else {}

        if canonical is None:
            status = "INVALID_STRUCTURE"

        elif record is None:
            status = "NOT_SCREENED"

        else:
            status = record["Screening_Status"]

        requested = int(
            record["Requested_Libraries"]
        ) if record else 0

        completed = int(
            record["Completed_Libraries"]
        ) if record else 0

        if completed > requested:
            raise ValueError(
                f"Invalid screening coverage: {canonical}"
            )

        if record:

            if len(strict) != int(record["Exact_Matches"]):
                raise ValueError(
                    f"Strict match count mismatch: {canonical}"
                )

            if len(reported) != int(
                record["Exact_Reported_Structure"]
            ):
                raise ValueError(
                    f"Reported match count mismatch: {canonical}"
                )

            if len(related) != int(record["Reported_Analogs"]):
                raise ValueError(
                    f"Analog count mismatch: {canonical}"
                )

        coverage = (
            completed / requested
            if requested else 0.0
        )

        if requested == 0 or completed == 0:
            coverage_status = "NOT_SCREENED"

        elif completed < requested:
            coverage_status = "PARTIAL"

        else:
            coverage_status = "COMPLETE"

        annotations.append({

            "Commercial_CanonicalSMILES":
                canonical or "",

            "Commercial_InChIKey":
                record["InChIKey"] if record else "",

            "Commercial_Screening_Status":
                status,

            "Commercial_Coverage_Status":
                coverage_status,

            "Commercial_Requested_Libraries":
                requested,

            "Commercial_Completed_Libraries":
                completed,

            "Commercial_Screening_Coverage":
                coverage,

            "Commercial_Exact_Matches":
                len(strict),

            "Commercial_Exact_Reported_Structure":
                len(reported),

            "Commercial_Total_Catalog_Matches":
                len(matches),

            "Commercial_Exact_Match_Libraries":
                ";".join(sorted({
                    item["Library"] for item in matches
                })),

            "Commercial_Catalog_IDs":
                ";".join(vendor_ids),

            "Commercial_Analog_Count":
                len(related),

            "Commercial_Best_Analog_Similarity":
                best.get("Similarity", ""),

            "Commercial_Best_Analog_Library":
                best.get("Library", ""),

            "Commercial_Best_Analog_Catalog_ID":
                best.get("Catalog_ID", ""),

            "Commercial_Best_Analog_SMILES":
                best.get("Matched_Canonical_SMILES", ""),

            "Commercial_Special_Matches":
                int(record["Special_Matches"])
                if record else 0,

            "Commercial_Screening_Method":
                screening_summary["Method"],

            "Commercial_Cache_Timestamps":
                json.dumps({
                    item["Library"]: item["Cache_Updated_At"]
                    for item in query_jobs
                    if item["Cache_Updated_At"]
                }),

            "Commercial_Stock_Status":
                "UNVERIFIED",
        })

    commercial = pd.DataFrame(
        annotations,
        index=original.index,
    )

    overlapping = set(original.columns) & set(
        commercial.columns
    )

    if overlapping:
        raise ValueError(
            f"Existing commercial columns: {overlapping}"
        )

    annotated = pd.concat(
        [original, commercial],
        axis=1,
    )

    # Verify preservation of all original data.
    pd.testing.assert_frame_equal(
        annotated[original.columns],
        original,
    )

    assert len(annotated) == len(original)

    exact_candidates = annotated[
        annotated["Commercial_Total_Catalog_Matches"] > 0
    ]

    strict_candidates = annotated[
        annotated["Commercial_Exact_Matches"] > 0
    ]

    analog_candidates = annotated[
        annotated["Commercial_Analog_Count"] > 0
    ]

    unscreened = annotated[
        annotated["Commercial_Completed_Libraries"] == 0
    ]

    annotated.to_csv(
        output_dir / "CMOS_V2_Commercial_Annotated.csv",
        index=False,
    )

    exact_candidates.to_csv(
        output_dir / "CMOS_V2_Commercial_Exact.csv",
        index=False,
    )

    strict_candidates.to_csv(
        output_dir / "CMOS_V2_Commercial_Strict_Exact.csv",
        index=False,
    )

    analog_candidates.to_csv(
        output_dir / "CMOS_V2_Commercial_Analog_Queries.csv",
        index=False,
    )

    unscreened.to_csv(
        output_dir / "CMOS_V2_Commercial_Unscreened.csv",
        index=False,
    )

    parent_summary = annotated.groupby(
        "ParentID",
        sort=False,
    ).agg(
        Molecules=("GeneratedID", "size"),

        Strict_Exact_Matches=(
            "Commercial_Exact_Matches",
            lambda x: int((x > 0).sum()),
        ),

        Catalog_Matched_Molecules=(
            "Commercial_Total_Catalog_Matches",
            lambda x: int((x > 0).sum()),
        ),

        Molecules_With_Analogs=(
            "Commercial_Analog_Count",
            lambda x: int((x > 0).sum()),
        ),

        Screened_Molecules=(
            "Commercial_Completed_Libraries",
            lambda x: int((x > 0).sum()),
        ),

        Mean_Screening_Coverage=(
            "Commercial_Screening_Coverage",
            "mean",
        ),
    ).reset_index()

    parent_summary.to_csv(
        output_dir / "Commercial_Per_Parent_Summary.csv",
        index=False,
    )

    summary = {
        "Input_File": str(input_file),
        "Input_Molecules": len(original),
        "Original_Columns": len(original.columns),
        "Commercial_Columns": len(commercial.columns),
        "Output_Columns": len(annotated.columns),
        "Original_Data_Preserved": True,
        "Catalog_Matched_Molecules": len(exact_candidates),
        "Strict_Exact_Matched_Molecules":
            len(strict_candidates),
        "Molecules_With_Analogs": len(analog_candidates),
        "Unscreened_Molecules": len(unscreened),
        "Screening_Status_Counts": dict(Counter(
            annotated["Commercial_Screening_Status"]
        )),
        "Coverage_Status_Counts": dict(Counter(
            annotated["Commercial_Coverage_Status"]
        )),
        "Mean_Screening_Coverage": float(
            annotated[
                "Commercial_Screening_Coverage"
            ].mean()
        ),
        "Stock_Status": "UNVERIFIED",
        "Screening_Summary": screening_summary,
    }

    (
        output_dir / "Commercial_Screening_Summary.json"
    ).write_text(
        json.dumps(summary, indent=4) + "\n"
    )

    print("\n" + "=" * 60)
    print("CMOS COMMERCIAL INTEGRATION COMPLETE")
    print("=" * 60)
    print(json.dumps(summary, indent=4))

    print("\nPer-parent summary:")
    print(parent_summary.to_string(index=False))

    return annotated, summary


def main():

    parser = argparse.ArgumentParser(
        description="Integrate commercial screening with CMOS"
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--screening-dir", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    integrate(
        args.input,
        args.screening_dir,
        args.output,
    )


if __name__ == "__main__":
    main()
