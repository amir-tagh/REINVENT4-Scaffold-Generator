#!/usr/bin/env python3

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs

RDLogger.DisableLog("rdApp.warning")

PROPERTIES = [
    "MW", "cLogP", "TPSA", "HBD", "HBA",
    "QED", "SA", "SCScore", "RotB",
    "FractionCSP3", "ParentSimilarity"
]


def analyze(root, output):
    output.mkdir(parents=True, exist_ok=True)

    records = []
    fingerprints = {}
    scaffold_counts = {}
    morgan = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048
    )

    for hit_dir in sorted(root.glob("Hit_*")):
        csv_file = hit_dir / "generated.csv"
        if not csv_file.exists():
            continue

        df = pd.read_csv(csv_file)
        df["Hit"] = hit_dir.name

        for col in PROPERTIES:
            if col in df:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        mols = [Chem.MolFromSmiles(str(s)) for s in df["SMILES"]]
        valid = [i for i, mol in enumerate(mols) if mol is not None]
        df = df.iloc[valid].copy()
        mols = [mols[i] for i in valid]

        df["CanonicalSMILES"] = [
            Chem.MolToSmiles(mol, isomericSmiles=True)
            for mol in mols
        ]

        df["MurckoScaffold"] = [
            MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
            for mol in mols
        ]

        fps = [morgan.GetFingerprint(mol) for mol in mols]

        fingerprints[hit_dir.name] = fps
        scaffold_counts[hit_dir.name] = df["MurckoScaffold"].nunique()

        records.append(df)

        print(
            f"{hit_dir.name}: {len(df)} valid molecules, "
            f"{scaffold_counts[hit_dir.name]} Murcko scaffolds"
        )

    if not records:
        raise RuntimeError("No generated.csv files found")

    combined = pd.concat(records, ignore_index=True)
    combined.to_csv(output / "all_molecules.csv", index=False)

    # Property statistics
    stats = (
        combined.groupby("Hit")[PROPERTIES]
        .agg(["mean", "std", "median", "min", "max"])
    )
    stats.to_csv(output / "property_statistics.csv")

    # Scaffold retention
    if "ScaffoldRetained" in combined:
        retention = (
            combined.assign(
                ScaffoldRetainedNumeric=combined[
                    "ScaffoldRetained"
                ].astype(str).str.lower().map({
                    "true": 1, "false": 0,
                    "1": 1, "0": 0
                })
            )
            .groupby("Hit")["ScaffoldRetainedNumeric"]
            .mean()
            .mul(100)
        )
    else:
        retention = pd.Series(dtype=float)

    # Diversity: exact mean pairwise Tanimoto dissimilarity
    diversity = {}

    for hit, fps in fingerprints.items():
        similarity_sum = 0.0
        pairs = 0

        for i in range(1, len(fps)):
            similarities = DataStructs.BulkTanimotoSimilarity(
                fps[i], fps[:i]
            )
            similarity_sum += sum(similarities)
            pairs += len(similarities)

        diversity[hit] = (
            1.0 - similarity_sum / pairs if pairs else np.nan
        )

    # Summary
    summary = combined.groupby("Hit").agg(
        Molecules=("SMILES", "size"),
        UniqueMolecules=("CanonicalSMILES", "nunique"),
        MeanQED=("QED", "mean"),
        MeanSA=("SA", "mean"),
        MeanSCScore=("SCScore", "mean"),
        MeanParentSimilarity=("ParentSimilarity", "mean")
    )

    summary["UniqueMurckoScaffolds"] = pd.Series(scaffold_counts)
    summary["MeanPairwiseDissimilarity"] = pd.Series(diversity)
    summary["ScaffoldRetentionPercent"] = retention

    summary.to_csv(output / "hit_summary.csv")

    # Cross-hit overlap
    overlap = (
        combined.groupby("CanonicalSMILES")["Hit"]
        .nunique()
    )

    duplicated = set(overlap[overlap > 1].index)

    combined[
        combined["CanonicalSMILES"].isin(duplicated)
    ].to_csv(output / "cross_hit_duplicates.csv", index=False)

    # Property distributions
    for prop in ["QED", "SA", "SCScore", "MW",
                 "cLogP", "TPSA", "ParentSimilarity"]:

        if prop not in combined:
            continue

        fig, ax = plt.subplots(figsize=(9, 5.5))

        data = [
            group[prop].dropna().values
            for _, group in combined.groupby("Hit")
        ]

        labels = [
            hit.replace("Hit_", "Hit ")
            for hit in sorted(combined["Hit"].unique())
        ]

        ax.boxplot(data, tick_labels=labels, showfliers=False)
        ax.set_ylabel(prop)
        ax.set_xlabel("Parent hit")
        ax.set_title(f"{prop} distribution across generated libraries")

        fig.tight_layout()
        fig.savefig(output / f"{prop}_boxplot.png", dpi=300)
        fig.savefig(output / f"{prop}_boxplot.pdf")
        plt.close(fig)

    print("\nAnalysis completed.")
    print(f"Total records: {len(combined):,}")
    print(f"Unique canonical molecules: {combined['CanonicalSMILES'].nunique():,}")
    print(f"Cross-hit shared molecules: {len(duplicated):,}")
    print(f"Results: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="results/v5/REINVENT_V5_All_8_Hits_Relaxed"
    )
    parser.add_argument(
        "--output",
        default="results/v5/REINVENT_V5_Analysis"
    )
    args = parser.parse_args()

    analyze(Path(args.input), Path(args.output))
