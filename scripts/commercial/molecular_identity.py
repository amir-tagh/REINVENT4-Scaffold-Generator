#!/usr/bin/env python3

"""Molecular identity verification for commercial screening."""

import argparse
import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path

from rdkit import Chem


@dataclass(frozen=True)
class MolecularIdentity:

    canonical_smiles: str
    connectivity_smiles: str
    inchikey: str
    connectivity_key: str
    specified_stereo: int
    unspecified_stereo: int


def get_identity(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"Invalid SMILES: {smiles}")

    canonical = Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True,
    )

    connectivity_mol = Chem.Mol(mol)
    Chem.RemoveStereochemistry(connectivity_mol)

    connectivity = Chem.MolToSmiles(
        connectivity_mol,
        canonical=True,
        isomericSmiles=True,
    )

    inchikey = Chem.MolToInchiKey(mol)

    if not inchikey:
        raise ValueError(
            f"InChIKey generation failed: {smiles}"
        )

    stereo = Chem.FindPotentialStereo(mol)

    specified = sum(
        item.specified == Chem.StereoSpecified.Specified
        for item in stereo
    )

    unspecified = sum(
        item.specified != Chem.StereoSpecified.Specified
        for item in stereo
    )

    return MolecularIdentity(
        canonical_smiles=canonical,
        connectivity_smiles=connectivity,
        inchikey=inchikey,
        connectivity_key=inchikey.split("-")[0],
        specified_stereo=specified,
        unspecified_stereo=unspecified,
    )


def classify_match(query, candidate, similarity):

    exact = (
        query.canonical_smiles
        == candidate.canonical_smiles
        and query.inchikey == candidate.inchikey
    )

    connectivity_match = (
        query.connectivity_smiles
        == candidate.connectivity_smiles
    )

    if exact:

        if (
            query.unspecified_stereo > 0
            or candidate.unspecified_stereo > 0
        ):
            return "EXACT_REPORTED_STRUCTURE"

        return "EXACT_MATCH"

    if connectivity_match:
        return "STEREOCHEMISTRY_MISMATCH"

    if similarity == 1.0:
        return "FINGERPRINT_ONLY_MATCH"

    return "STRUCTURAL_ANALOG"


def verify_results(query_smiles, input_csv, output_dir):

    query = get_identity(query_smiles)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = []

    summary = {
        "Total": 0,
        "EXACT_MATCH": 0,
        "EXACT_REPORTED_STRUCTURE": 0,
        "STEREOCHEMISTRY_MISMATCH": 0,
        "FINGERPRINT_ONLY_MATCH": 0,
        "STRUCTURAL_ANALOG": 0,
        "INVALID_STRUCTURE": 0,
    }

    with open(input_csv, newline="") as handle:

        reader = csv.DictReader(handle)

        required = {
            "Library",
            "Catalog_ID",
            "Similarity",
            "Matched_SMILES",
        }

        if not required.issubset(
            reader.fieldnames or []
        ):
            raise ValueError(
                "Missing required screening columns"
            )

        for row in reader:

            summary["Total"] += 1

            record = dict(row)

            record.update({
                "Query_SMILES": query.canonical_smiles,
                "Query_InChIKey": query.inchikey,
                "Query_Unspecified_Stereo":
                    query.unspecified_stereo,
                "Matched_Canonical_SMILES": "",
                "Matched_InChIKey": "",
                "Matched_Unspecified_Stereo": "",
                "Match_Type": "",
                "Exact_Match": False,
                "Connectivity_Match": False,
                "Availability_Status": "UNVERIFIED",
                "Error": "",
            })

            try:

                candidate = get_identity(
                    row["Matched_SMILES"]
                )

                similarity = float(row["Similarity"])

                match_type = classify_match(
                    query,
                    candidate,
                    similarity,
                )

                record.update({
                    "Matched_Canonical_SMILES":
                        candidate.canonical_smiles,
                    "Matched_InChIKey":
                        candidate.inchikey,
                    "Matched_Unspecified_Stereo":
                        candidate.unspecified_stereo,
                    "Match_Type": match_type,
                    "Exact_Match": (
                        match_type in {
                            "EXACT_MATCH",
                            "EXACT_REPORTED_STRUCTURE",
                        }
                    ),
                    "Connectivity_Match": (
                        query.connectivity_smiles
                        == candidate.connectivity_smiles
                    ),
                })

            except ValueError as error:

                match_type = "INVALID_STRUCTURE"

                record["Match_Type"] = match_type
                record["Error"] = str(error)

            summary[match_type] += 1
            records.append(record)

    if not records:
        raise ValueError("No screening records found")

    fields = list(records[0].keys())

    def write_csv(filename, rows):

        with (output_dir / filename).open(
            "w", newline=""
        ) as handle:

            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
            )

            writer.writeheader()
            writer.writerows(rows)

    exact = [
        row for row in records
        if row["Exact_Match"]
    ]

    analogs = [
        row for row in records
        if row["Match_Type"] == "STRUCTURAL_ANALOG"
    ]

    stereo_matches = [
        row for row in records
        if row["Match_Type"]
        == "STEREOCHEMISTRY_MISMATCH"
    ]

    write_csv("verified_results.csv", records)
    write_csv("exact_matches.csv", exact)
    write_csv("structural_analogs.csv", analogs)
    write_csv("stereochemistry_matches.csv", stereo_matches)

    report = {
        "Query": asdict(query),
        "Summary": summary,
        "Availability_Status": "UNVERIFIED",
    }

    (output_dir / "identity_summary.json").write_text(
        json.dumps(report, indent=4) + "\n"
    )

    print("\n" + "=" * 60)
    print("MOLECULAR IDENTITY VERIFICATION")
    print("=" * 60)

    print("Query SMILES:", query.canonical_smiles)
    print("Query InChIKey:", query.inchikey)
    print()

    for key, value in summary.items():
        print(f"{key:30s}: {value}")

    print("\nVerified catalog structure matches:")

    for row in exact:

        print(
            row["Library"],
            row["Catalog_ID"],
            row["Similarity"],
            row["Match_Type"],
        )

    print("\nOutput directory:", output_dir)
    print("=" * 60)

    return report


def main():

    parser = argparse.ArgumentParser(
        description="Verify commercial molecular identity"
    )

    parser.add_argument("--smiles", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    verify_results(
        args.smiles,
        args.input,
        args.output,
    )


if __name__ == "__main__":
    main()
