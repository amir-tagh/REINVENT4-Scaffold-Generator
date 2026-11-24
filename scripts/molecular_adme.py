#!/usr/bin/env python

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ADME_WEIGHTS = {
    "Score_ADME_HIA": 0.20,
    "Score_ADME_Bioavailability": 0.20,
    "Score_ADME_PAMPA": 0.15,
    "Score_ADME_Pgp": 0.20,
    "Score_ADME_CYP": 0.25
}


CYP_ENDPOINTS = [
    "CYP1A2_Veith",
    "CYP2C19_Veith",
    "CYP2C9_Veith",
    "CYP2D6_Veith",
    "CYP3A4_Veith"
]


def probability(df, column):

    if column not in df.columns:
        raise ValueError(
            f"Missing required column: {column}"
        )

    values = pd.to_numeric(
        df[column], errors="coerce"
    )

    invalid = (
        values.notna()
        & ~values.between(0, 1)
    )

    if invalid.any():
        raise ValueError(
            f"{column}: {invalid.sum()} "
            "values outside [0, 1]"
        )

    return values


def calculate_adme_scores(df):

    df = df.copy()

    # Absorption
    df["Score_ADME_HIA"] = probability(
        df, "ADMETAI_HIA_Hou"
    )

    df["Score_ADME_Bioavailability"] = probability(
        df, "ADMETAI_Bioavailability_Ma"
    )

    df["Score_ADME_PAMPA"] = probability(
        df, "ADMETAI_PAMPA_NCATS"
    )

    # Efflux liability
    df["Score_ADME_Pgp"] = 1.0 - probability(
        df, "ADMETAI_Pgp_Broccatelli"
    )

    # Individual CYP inhibition scores
    cyp_columns = []

    for endpoint in CYP_ENDPOINTS:

        enzyme = endpoint.split("_")[0]

        score_column = f"Score_ADME_{enzyme}"

        df[score_column] = 1.0 - probability(
            df, f"ADMETAI_{endpoint}"
        )

        cyp_columns.append(score_column)

    # Mean CYP desirability
    df["Score_ADME_CYP"] = df[
        cyp_columns
    ].sum(
        axis=1,
        min_count=len(cyp_columns)
    ) / len(cyp_columns)

    # Worst individual CYP liability
    df["Score_ADME_CYP_Worst"] = df[
        cyp_columns
    ].min(axis=1, skipna=False)

    # BBB is an independent CMOS component
    df["Score_BBB"] = probability(
        df, "ADMETAI_BBB_Martins"
    )

    # Weighted ADME score
    components = list(ADME_WEIGHTS)

    weighted = df[components].mul(
        pd.Series(ADME_WEIGHTS),
        axis=1
    )

    df["Score_ADME"] = weighted.sum(
        axis=1,
        min_count=len(components)
    )

    df["ADME_MissingComponents"] = df[
        components
    ].isna().sum(axis=1)

    df["ADME_Status"] = np.where(
        df["ADME_MissingComponents"].eq(0),
        "OK",
        "Incomplete"
    )

    return df


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    df = pd.read_csv(args.input)

    print("Input molecules:", len(df))

    output = calculate_adme_scores(df)

    if len(output) != len(df):
        raise RuntimeError("Unexpected row count change")

    path = Path(args.output)

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    output.to_csv(path, index=False)

    columns = list(ADME_WEIGHTS) + [
        "Score_ADME_CYP_Worst",
        "Score_BBB",
        "Score_ADME"
    ]

    print("\nADME scoring coverage:")
    print(output[columns].count())

    print("\nADME score statistics:")
    print(output[columns].describe().round(3))

    print("\nADME status:")
    print(output["ADME_Status"].value_counts())

    summary = {
        "total_molecules": len(output),
        "weights": ADME_WEIGHTS,
        "complete_scores": int(
            output["Score_ADME"].count()
        ),
        "mean_adme_score": float(
            output["Score_ADME"].mean()
        ),
        "mean_bbb_score": float(
            output["Score_BBB"].mean()
        ),
        "adme_ge_0_7": int(
            output["Score_ADME"].ge(0.7).sum()
        ),
        "bbb_ge_0_5": int(
            output["Score_BBB"].ge(0.5).sum()
        )
    }

    summary_path = path.with_suffix(
        ".summary.json"
    )

    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=4)

    print("\nOutput:", path)
    print("Summary:", summary_path)


if __name__ == "__main__":
    main()
