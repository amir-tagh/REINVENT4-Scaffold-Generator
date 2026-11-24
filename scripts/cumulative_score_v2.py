#!/usr/bin/env python

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


COMPONENTS = [
    "QED",
    "SA",
    "SCScore",
    "ParentSimilarity",
    "Diversity",
    "Solubility",
    "CNS",
    "BBB",
    "ADME",
    "Toxicity"
]


SOURCE_COLUMNS = {
    "QED": "QED",
    "SA": "SA",
    "SCScore": "SCScore",
    "ParentSimilarity": "ParentSimilarity",
    "Solubility": "Score_Solubility_Consensus",
    "CNS": "Score_CNS_Partial",
    "BBB": "Score_BBB",
    "ADME": "Score_ADME",
    "Toxicity": "Score_Toxicity"
}


FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)


def load_profile(path, profile):

    with open(path) as handle:
        profiles = json.load(handle)

    if profile not in profiles:
        raise ValueError(
            f"Unknown profile: {profile}. "
            f"Available: {list(profiles)}"
        )

    weights = profiles[profile]

    if set(weights) != set(COMPONENTS):
        raise ValueError(
            "Profile contains missing or unknown components"
        )

    weights = {
        k: float(v)
        for k, v in weights.items()
    }

    if not all(
        np.isfinite(v) and v >= 0
        for v in weights.values()
    ):
        raise ValueError("Invalid profile weights")

    if not np.isclose(sum(weights.values()), 1.0):
        raise ValueError("Weights must sum to 1")

    if weights["Diversity"] >= 1.0:
        raise ValueError(
            "Diversity weight must be less than 1"
        )

    return weights


def calculate_component_scores(
    df,
    weights,
    parent_center=0.55,
    parent_sigma=0.20
):

    df = df.copy()

    for component, source in SOURCE_COLUMNS.items():

        target = f"CMOS_V2_{component}"

        if source not in df.columns:

            if weights[component] > 0:
                raise ValueError(
                    f"Missing required column: {source}"
                )

            df[target] = np.nan
            continue

        values = pd.to_numeric(
            df[source],
            errors="coerce"
        ).replace([np.inf, -np.inf], np.nan)

        if component == "SA":

            values = values.where(
                values.between(1, 10)
            )

            score = (10.0 - values) / 9.0

        elif component == "SCScore":

            values = values.where(
                values.between(1, 5)
            )

            score = (5.0 - values) / 4.0

        elif component == "ParentSimilarity":

            values = values.where(
                values.between(0, 1)
            )

            score = np.exp(
                -0.5 * (
                    (values - parent_center)
                    / parent_sigma
                ) ** 2
            )

        else:

            score = values.where(
                values.between(0, 1)
            )

        df[target] = score.clip(0, 1)

    active = {
        f"CMOS_V2_{k}": v
        for k, v in weights.items()
        if k != "Diversity" and v > 0
    }

    columns = list(active)

    df["CMOS_V2_MissingComponents"] = (
        df[columns].isna().sum(axis=1)
    )

    weighted = df[columns].mul(
        pd.Series(active),
        axis=1
    )

    df["CMOS_V2_Intrinsic"] = (
        weighted.sum(
            axis=1,
            min_count=len(columns)
        ) / sum(active.values())
    )

    return df


def prepare_fingerprints(df):

    cache = {}
    fingerprints = {}
    canonical = {}

    for idx, smiles in df["SMILES"].items():

        if pd.isna(smiles):
            continue

        key = str(smiles)

        if key not in cache:

            mol = Chem.MolFromSmiles(key)

            if mol is None:
                cache[key] = None

            else:
                cache[key] = (
                    Chem.MolToSmiles(
                        mol,
                        canonical=True,
                        isomericSmiles=True
                    ),
                    FP_GENERATOR.GetFingerprint(mol)
                )

        result = cache[key]

        if result is not None:

            canonical[idx] = result[0]
            fingerprints[idx] = result[1]

    return fingerprints, canonical


def select_diverse_molecules(
    df,
    weights,
    select_per_hit=30
):

    df = df.copy()

    fingerprints, canonical = (
        prepare_fingerprints(df)
    )

    df["CMOS_V2_ValidSMILES"] = (
        df.index.isin(fingerprints)
    )

    df["CMOS_V2_CanonicalSMILES"] = (
        pd.Series(canonical, dtype="object")
    )

    df["CMOS_V2_Diversity"] = np.nan
    df["CMOS_V2_Cumulative"] = np.nan
    df["CMOS_V2_SelectionRank"] = pd.Series(
        pd.NA,
        index=df.index,
        dtype="Int64"
    )

    df["CMOS_V2_Selected"] = False

    df["CMOS_V2_Eligible"] = (
        df["CMOS_V2_ValidSMILES"]
        & df["CMOS_V2_Intrinsic"].notna()
        & df["ParentID"].notna()
    )

    diversity_weight = weights["Diversity"]

    for parent, group in df.groupby(
        "ParentID",
        sort=True
    ):

        candidates = group.loc[
            group["CMOS_V2_Eligible"]
        ].copy()

        # Deterministic ordering and duplicate removal
        candidates = candidates.sort_values(
            [
                "CMOS_V2_Intrinsic",
                "CMOS_V2_CanonicalSMILES"
            ],
            ascending=[False, True],
            kind="stable"
        )

        candidates = candidates.drop_duplicates(
            subset="CMOS_V2_CanonicalSMILES",
            keep="first"
        )

        indices = candidates.index.to_numpy()

        if len(indices) == 0:
            continue

        fps = [
            fingerprints[idx]
            for idx in indices
        ]

        intrinsic = candidates[
            "CMOS_V2_Intrinsic"
        ].to_numpy(dtype=float)

        base_scores = (
            (1.0 - diversity_weight)
            * intrinsic
        )

        max_similarity = np.zeros(
            len(indices),
            dtype=float
        )

        available = np.ones(
            len(indices),
            dtype=bool
        )

        n_select = min(
            select_per_hit,
            len(indices)
        )

        for rank in range(1, n_select + 1):

            diversity = 1.0 - max_similarity

            cumulative = (
                base_scores
                + diversity_weight * diversity
            )

            cumulative[~available] = -np.inf

            best = int(np.argmax(cumulative))

            idx = indices[best]

            df.at[
                idx, "CMOS_V2_Diversity"
            ] = diversity[best]

            df.at[
                idx, "CMOS_V2_Cumulative"
            ] = cumulative[best]

            df.at[
                idx, "CMOS_V2_SelectionRank"
            ] = rank

            df.at[
                idx, "CMOS_V2_Selected"
            ] = True

            available[best] = False

            similarities = np.asarray(
                DataStructs.BulkTanimotoSimilarity(
                    fps[best],
                    fps
                ),
                dtype=float
            )

            max_similarity = np.maximum(
                max_similarity,
                similarities
            )

        print(
            f"{parent}: "
            f"{len(group)} input, "
            f"{len(indices)} eligible unique, "
            f"{n_select} selected"
        )

    return df


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument(
        "--profiles",
        default="config/optimization_profiles.json"
    )

    parser.add_argument(
        "--profile",
        default="cns"
    )

    parser.add_argument(
        "--select-per-hit",
        type=int,
        default=30
    )

    parser.add_argument(
        "--parent-center",
        type=float,
        default=0.55
    )

    parser.add_argument(
        "--parent-sigma",
        type=float,
        default=0.20
    )

    args = parser.parse_args()

    if args.select_per_hit < 1:
        parser.error("--select-per-hit must be positive")

    if not 0 <= args.parent_center <= 1:
        parser.error("--parent-center must be in [0,1]")

    if (
        not np.isfinite(args.parent_sigma)
        or args.parent_sigma <= 0
    ):
        parser.error("--parent-sigma must be positive")

    weights = load_profile(
        args.profiles,
        args.profile
    )

    df = pd.read_csv(args.input)

    required = ["ParentID", "SMILES"]

    missing = [
        c for c in required
        if c not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing columns: {missing}"
        )

    print("Input molecules:", len(df))
    print("Optimization profile:", args.profile)

    print("\nProfile weights:")
    for key, value in weights.items():
        print(f"{key:20s} {value:.3f}")

    df = calculate_component_scores(
        df,
        weights,
        parent_center=args.parent_center,
        parent_sigma=args.parent_sigma
    )

    print("\nDiversity-aware selection:")

    df = select_diverse_molecules(
        df,
        weights,
        select_per_hit=args.select_per_hit
    )

    df["CMOS_V2_Profile"] = args.profile

    selected = df.loc[
        df["CMOS_V2_Selected"]
    ].copy()

    selected = selected.sort_values(
        ["ParentID", "CMOS_V2_SelectionRank"]
    )

    path = Path(args.output)

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    df.to_csv(path, index=False)

    selected_path = path.with_name(
        path.stem + "_selected.csv"
    )

    selected.to_csv(
        selected_path,
        index=False
    )

    smiles_path = path.with_name(
        path.stem + "_selected.smi"
    )

    with open(smiles_path, "w") as handle:

        for _, row in selected.iterrows():

            name = (
                f"{row['ParentID']}_"
                f"Rank_{int(row['CMOS_V2_SelectionRank']):03d}"
            )

            handle.write(
                f"{row['CMOS_V2_CanonicalSMILES']}"
                f"\t{name}\n"
            )

    columns = [
        "CMOS_V2_Intrinsic",
        "CMOS_V2_Diversity",
        "CMOS_V2_Cumulative"
    ]

    print("\nTotal molecules:", len(df))

    print(
        "Eligible molecules:",
        int(df["CMOS_V2_Eligible"].sum())
    )

    print(
        "Selected molecules:",
        len(selected)
    )

    print("\nFull-library intrinsic score:")
    print(
        df["CMOS_V2_Intrinsic"]
        .describe()
        .round(4)
    )

    print("\nSelected molecule scores:")
    print(
        selected[columns]
        .describe()
        .round(4)
    )

    summary = {
        "profile": args.profile,
        "weights": weights,
        "input": args.input,
        "total_molecules": len(df),
        "eligible_molecules": int(
            df["CMOS_V2_Eligible"].sum()
        ),
        "selected_molecules": len(selected),
        "select_per_hit": args.select_per_hit,
        "parent_similarity_center": args.parent_center,
        "parent_similarity_sigma": args.parent_sigma,
        "pka_used_for_optimization": False,
        "fingerprint": {
            "type": "Morgan",
            "radius": 2,
            "bits": 2048
        },
        "selected_per_parent": {
            str(k): int(v)
            for k, v in selected[
                "ParentID"
            ].value_counts().items()
        }
    }

    summary_path = path.with_suffix(
        ".summary.json"
    )

    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=4)

    print("\nFull output:", path)
    print("Selected CSV:", selected_path)
    print("Selected SMILES:", smiles_path)
    print("Summary:", summary_path)


if __name__ == "__main__":
    main()
