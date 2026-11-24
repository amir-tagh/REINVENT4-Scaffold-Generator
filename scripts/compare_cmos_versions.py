#!/usr/bin/env python

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold


FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)


PROPERTY_COLUMNS = [
    "QED",
    "SA",
    "SCScore",
    "MW",
    "cLogP",
    "TPSA",
    "HBD",
    "HBA",
    "RotB",
    "FractionCSP3",
    "ParentSimilarity",
    "LogS_ESOL",
    "LogS_ADMETAI",
    "LogS_Consensus",
    "LogD_7_4_ADMETAI",
    "CNS_MPO",
    "Score_CNS_Partial",
    "Score_BBB",
    "Score_ADME",
    "Score_Toxicity",
    "Toxicity_Alert_Count"
]


def canonicalize(smiles):

    if pd.isna(smiles):
        raise ValueError("Missing SMILES")

    mol = Chem.MolFromSmiles(str(smiles))

    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"Invalid SMILES: {smiles}")

    return Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True
    )


def prepare_library(df, parent_column):

    df = df.copy()

    if df[parent_column].isna().any():
        raise ValueError("Missing parent identifiers")

    df["BenchmarkParent"] = (
        df[parent_column].astype(str)
    )

    cache = {
        s: canonicalize(s)
        for s in df["SMILES"].unique()
    }

    df["BenchmarkSMILES"] = (
        df["SMILES"].map(cache)
    )

    keys = [
        "BenchmarkParent",
        "BenchmarkSMILES"
    ]

    if df.duplicated(keys).any():
        raise ValueError(
            "Duplicate parent-molecule identities detected"
        )

    return df


def molecular_keys(df):

    return set(
        zip(
            df["BenchmarkParent"],
            df["BenchmarkSMILES"]
        )
    )


def selected_mask(values):

    values = (
        values.astype(str)
        .str.strip()
        .str.lower()
    )

    valid = values.isin([
        "true", "false", "1", "0"
    ])

    if not valid.all():
        raise ValueError(
            "Unexpected selection status values"
        )

    return values.isin(["true", "1"])


def compare_properties(groups):

    records = []

    for prop in PROPERTY_COLUMNS:

        row = {"Property": prop}

        for name, df in groups.items():

            values = pd.to_numeric(
                df[prop],
                errors="coerce"
            )

            row[f"{name}_Mean"] = values.mean()
            row[f"{name}_SD"] = values.std()
            row[f"{name}_Median"] = values.median()
            row[f"{name}_N"] = values.count()

        row["V1_vs_Library"] = (
            row["CMOS_V1_Mean"]
            - row["Full_Library_Mean"]
        )

        row["V2_vs_Library"] = (
            row["CMOS_V2_Mean"]
            - row["Full_Library_Mean"]
        )

        row["V2_vs_V1"] = (
            row["CMOS_V2_Mean"]
            - row["CMOS_V1_Mean"]
        )

        records.append(row)

    return pd.DataFrame(records)


def threshold_analysis(groups):

    records = []

    log100 = np.log10(100e-6)
    log300 = np.log10(300e-6)

    for name, df in groups.items():

        conditions = {
            "QED_GE_0.7":
                df["QED"].ge(0.7),

            "SA_LE_3":
                df["SA"].le(3),

            "SCScore_LE_3.5":
                df["SCScore"].le(3.5),

            "Consensus_Solubility_GE_100uM":
                df["LogS_Consensus"].ge(log100),

            "Consensus_Solubility_GE_300uM":
                df["LogS_Consensus"].ge(log300),

            "Both_Models_GE_100uM":
                (
                    df["LogS_ESOL"].ge(log100)
                    & df["LogS_ADMETAI"].ge(log100)
                ),

            "Both_Models_GE_300uM":
                (
                    df["LogS_ESOL"].ge(log300)
                    & df["LogS_ADMETAI"].ge(log300)
                ),

            "Partial_CNS_GE_0.8":
                df["Score_CNS_Partial"].ge(0.8),

            "Full_CNS_MPO_GE_4":
                df["CNS_MPO"].ge(4),

            "BBB_GE_0.5":
                df["Score_BBB"].ge(0.5),

            "ADME_GE_0.7":
                df["Score_ADME"].ge(0.7),

            "Toxicity_Desirability_GE_0.7":
                df["Score_Toxicity"].ge(0.7),

            "Zero_Toxicity_Alerts":
                df["Toxicity_Alert_Count"].eq(0)
        }

        for metric, mask in conditions.items():

            count = int(mask.sum())

            records.append({
                "Group": name,
                "Metric": metric,
                "Count": count,
                "Total": len(df),
                "Fraction": count / len(df),
                "Percentage": 100 * count / len(df)
            })

    return pd.DataFrame(records)


def diversity_analysis(df):

    # Use unique molecular structures so that
    # cross-parent duplicates do not inflate similarity.
    smiles = sorted(
        df["BenchmarkSMILES"].unique()
    )

    mols = [
        Chem.MolFromSmiles(s)
        for s in smiles
    ]

    fps = [
        FP_GENERATOR.GetFingerprint(m)
        for m in mols
    ]

    scaffolds = [
        MurckoScaffold.MurckoScaffoldSmiles(
            mol=m,
            includeChirality=False
        )
        for m in mols
    ]

    n = len(fps)

    nearest = np.full(n, np.nan)

    if n > 1:
        nearest[:] = 0.0

    similarity_sum = 0.0
    similarity_squared_sum = 0.0
    pair_count = 0

    for i in range(1, n):

        similarities = np.asarray(
            DataStructs.BulkTanimotoSimilarity(
                fps[i],
                fps[:i]
            ),
            dtype=float
        )

        similarity_sum += similarities.sum()

        similarity_squared_sum += (
            similarities ** 2
        ).sum()

        pair_count += len(similarities)

        nearest[i] = similarities.max()

        nearest[:i] = np.maximum(
            nearest[:i],
            similarities
        )

    mean_similarity = (
        similarity_sum / pair_count
        if pair_count else np.nan
    )

    similarity_sd = (
        np.sqrt(
            max(
                0.0,
                similarity_squared_sum / pair_count
                - mean_similarity ** 2
            )
        )
        if pair_count else np.nan
    )

    nonempty_scaffolds = [
        s for s in scaffolds if s
    ]

    return {
        "Total_Rows": len(df),
        "Unique_Molecules": n,
        "Cross_Parent_Duplicate_Rows": len(df) - n,
        "Unique_Murcko_Scaffolds":
            len(set(nonempty_scaffolds)),
        "Acyclic_Molecules":
            scaffolds.count(""),
        "Pairwise_Comparisons": pair_count,
        "Mean_Pairwise_Tanimoto": mean_similarity,
        "SD_Pairwise_Tanimoto": similarity_sd,
        "Mean_Pairwise_Diversity":
            1.0 - mean_similarity,
        "Mean_Nearest_Neighbor_Tanimoto":
            float(nearest.mean()) if n > 1 else np.nan,
        "Median_Nearest_Neighbor_Tanimoto":
            float(np.median(nearest)) if n > 1 else np.nan
    }


def overlap_analysis(v1, v2):

    records = []

    parents = sorted(
        set(v1["BenchmarkParent"])
        | set(v2["BenchmarkParent"])
    )

    for parent in ["ALL"] + parents:

        if parent == "ALL":
            a = molecular_keys(v1)
            b = molecular_keys(v2)
        else:
            a = molecular_keys(
                v1[v1["BenchmarkParent"].eq(parent)]
            )
            b = molecular_keys(
                v2[v2["BenchmarkParent"].eq(parent)]
            )

        intersection = a & b
        union = a | b

        records.append({
            "Parent": parent,
            "V1_Selected": len(a),
            "V2_Selected": len(b),
            "Shared": len(intersection),
            "V1_Exclusive": len(a - b),
            "V2_Exclusive": len(b - a),
            "Union": len(union),
            "Jaccard": (
                len(intersection) / len(union)
                if union else np.nan
            )
        })

    return pd.DataFrame(records)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--v1", required=True)
    parser.add_argument("--v2", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument(
        "--select-per-hit",
        type=int,
        default=30
    )

    parser.add_argument(
        "--full-diversity",
        action="store_true"
    )

    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    v1 = prepare_library(
        pd.read_csv(args.v1),
        "Hit"
    )

    v2 = prepare_library(
        pd.read_csv(args.v2),
        "ParentID"
    )

    keys_v1 = molecular_keys(v1)
    keys_v2 = molecular_keys(v2)

    print("CMOS V1 molecules:", len(v1))
    print("CMOS V2 molecules:", len(v2))

    print(
        "Shared molecular identities:",
        len(keys_v1 & keys_v2)
    )

    print("V1-only identities:", len(keys_v1 - keys_v2))
    print("V2-only identities:", len(keys_v2 - keys_v1))

    if keys_v1 != keys_v2:
        raise ValueError(
            "Input molecular libraries are not identical"
        )

    print("\nMolecular identity verification: PASSED")

    missing = [
        c for c in PROPERTY_COLUMNS
        if c not in v2.columns
    ]

    if missing:
        raise ValueError(
            f"Missing evaluation properties: {missing}"
        )

    v2[PROPERTY_COLUMNS] = v2[
        PROPERTY_COLUMNS
    ].apply(pd.to_numeric, errors="coerce")

    if not np.isfinite(
        v2[PROPERTY_COLUMNS].to_numpy(dtype=float)
    ).all():
        raise ValueError(
            "Missing or non-finite evaluation properties"
        )

    v1_selected = v1.loc[
        selected_mask(v1["Selected"])
    ].copy()

    v2_selected = v2.loc[
        selected_mask(v2["CMOS_V2_Selected"])
    ].copy()

    # Confirm equal selection budget for every parent.
    for name, selected in [
        ("V1", v1_selected),
        ("V2", v2_selected)
    ]:

        counts = (
            selected["BenchmarkParent"]
            .value_counts()
            .reindex(
                sorted(v2["BenchmarkParent"].unique()),
                fill_value=0
            )
        )

        if not counts.eq(args.select_per_hit).all():
            raise ValueError(
                f"{name}: unexpected selection counts"
            )

    keys = [
        "BenchmarkParent",
        "BenchmarkSMILES"
    ]

    # Attach identical evaluation properties to V1
    # selections using the existing V2/ADMET library.
    v1_selected = v1_selected[
        keys + [
            "IntrinsicScore",
            "CumulativeScore",
            "SelectionRank",
            "Score_Diversity"
        ]
    ].merge(
        v2,
        on=keys,
        how="left",
        validate="one_to_one"
    )

    groups = {
        "Full_Library": v2,
        "CMOS_V1": v1_selected,
        "CMOS_V2": v2_selected
    }

    print("\nSelected molecules:")
    print("CMOS V1:", len(v1_selected))
    print("CMOS V2:", len(v2_selected))

    # Property comparison
    properties = compare_properties(groups)

    properties.to_csv(
        output / "property_comparison.csv",
        index=False
    )

    print("\nProperty comparison:")

    print(
        properties[[
            "Property",
            "Full_Library_Mean",
            "CMOS_V1_Mean",
            "CMOS_V2_Mean",
            "V2_vs_V1"
        ]].round(4).to_string(index=False)
    )

    # Threshold comparison
    thresholds = threshold_analysis(groups)

    thresholds.to_csv(
        output / "threshold_comparison.csv",
        index=False
    )

    threshold_table = thresholds.pivot(
        index="Metric",
        columns="Group",
        values="Percentage"
    )

    print("\nThreshold comparison (%):")
    print(threshold_table.round(2).to_string())

    # Solubility geometric means
    solubility_records = []

    for name, df in groups.items():

        for prop in [
            "LogS_ESOL",
            "LogS_ADMETAI",
            "LogS_Consensus"
        ]:

            mean_log = df[prop].mean()

            solubility_records.append({
                "Group": name,
                "Model": prop,
                "Mean_LogS": mean_log,
                "Geometric_Mean_uM":
                    1e6 * 10 ** mean_log
            })

    solubility = pd.DataFrame(solubility_records)

    solubility.to_csv(
        output / "solubility_comparison.csv",
        index=False
    )

    print("\nSolubility geometric means (uM):")

    print(
        solubility.pivot(
            index="Model",
            columns="Group",
            values="Geometric_Mean_uM"
        ).round(2).to_string()
    )

    # Selection overlap
    overlap = overlap_analysis(
        v1_selected,
        v2_selected
    )

    overlap.to_csv(
        output / "selection_overlap.csv",
        index=False
    )

    print("\nSelection overlap:")
    print(overlap.round(4).to_string(index=False))

    # Structural diversity
    diversity_records = []

    for name, df in groups.items():

        if name == "Full_Library" and not args.full_diversity:
            continue

        print(f"\nCalculating diversity: {name}")

        record = diversity_analysis(df)
        record["Group"] = name

        diversity_records.append(record)

    diversity = pd.DataFrame(diversity_records)

    diversity.to_csv(
        output / "diversity_comparison.csv",
        index=False
    )

    print("\nStructural diversity:")

    print(
        diversity.set_index("Group")
        .T.round(4).to_string()
    )

    # Parent-specific property and diversity analysis
    parent_records = []
    parent_diversity = []

    for parent in sorted(v2["BenchmarkParent"].unique()):

        parent_groups = {
            name: df[
                df["BenchmarkParent"].eq(parent)
            ]
            for name, df in groups.items()
        }

        table = compare_properties(parent_groups)
        table.insert(0, "ParentID", parent)

        parent_records.append(table)

        for name in ["CMOS_V1", "CMOS_V2"]:

            record = diversity_analysis(
                parent_groups[name]
            )

            record["ParentID"] = parent
            record["Group"] = name

            parent_diversity.append(record)

    pd.concat(
        parent_records,
        ignore_index=True
    ).to_csv(
        output / "per_parent_comparison.csv",
        index=False
    )

    pd.DataFrame(parent_diversity).to_csv(
        output / "per_parent_diversity.csv",
        index=False
    )

    # Export selection membership for all molecules.
    selected_keys_v1 = molecular_keys(v1_selected)
    selected_keys_v2 = molecular_keys(v2_selected)

    membership = v2.copy()

    row_keys = list(
        zip(
            membership["BenchmarkParent"],
            membership["BenchmarkSMILES"]
        )
    )

    membership["Benchmark_Selected_V1"] = [
        k in selected_keys_v1 for k in row_keys
    ]

    membership["Benchmark_Selected_V2"] = [
        k in selected_keys_v2 for k in row_keys
    ]

    a = membership["Benchmark_Selected_V1"]
    b = membership["Benchmark_Selected_V2"]

    membership["Benchmark_Selection_Group"] = np.select(
        [a & b, a & ~b, ~a & b],
        ["Shared", "V1_Only", "V2_Only"],
        default="Neither"
    )

    membership.to_csv(
        output / "library_selection_membership.csv",
        index=False
    )

    v1_selected.to_csv(
        output / "CMOS_V1_Selected_Annotated.csv",
        index=False
    )

    v2_selected.to_csv(
        output / "CMOS_V2_Selected_Annotated.csv",
        index=False
    )

    summary = {
        "identical_input_libraries": True,
        "total_molecules": len(v2),
        "v1_selected": len(v1_selected),
        "v2_selected": len(v2_selected),
        "shared_selected": len(
            selected_keys_v1 & selected_keys_v2
        ),
        "v1_unique_selected_structures":
            int(v1_selected["BenchmarkSMILES"].nunique()),
        "v2_unique_selected_structures":
            int(v2_selected["BenchmarkSMILES"].nunique()),
        "select_per_hit": args.select_per_hit,
        "fingerprint_radius": 2,
        "fingerprint_bits": 2048,
        "full_library_diversity": args.full_diversity
    }

    with open(
        output / "benchmark_summary.json", "w"
    ) as handle:

        json.dump(summary, handle, indent=4)

    print("\nBENCHMARK COMPLETE")
    print("Output directory:", output)


if __name__ == "__main__":
    main()
