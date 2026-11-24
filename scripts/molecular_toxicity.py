#!/usr/bin/env python

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


TOXICITY_WEIGHTS = {
    "Score_Tox_AMES": 0.25,
    "Score_Tox_hERG": 0.25,
    "Score_Tox_DILI": 0.20,
    "Score_Tox_ClinTox": 0.20,
    "Score_Tox_Carcinogens": 0.10
}


TOXICITY_ENDPOINTS = {
    "Score_Tox_AMES": "ADMETAI_AMES",
    "Score_Tox_hERG": "ADMETAI_hERG",
    "Score_Tox_DILI": "ADMETAI_DILI",
    "Score_Tox_ClinTox": "ADMETAI_ClinTox",
    "Score_Tox_Carcinogens": "ADMETAI_Carcinogens_Lagunin"
}


def calculate_toxicity_scores(df):

    df = df.copy()

    for score, endpoint in TOXICITY_ENDPOINTS.items():

        if endpoint not in df.columns:
            raise ValueError(
                f"Missing required column: {endpoint}"
            )

        values = pd.to_numeric(
            df[endpoint], errors="coerce"
        )

        invalid = (
            values.notna()
            & ~values.between(0, 1)
        )

        if invalid.any():
            raise ValueError(
                f"{endpoint}: {invalid.sum()} "
                "values outside [0, 1]"
            )

        df[score] = 1.0 - values

    components = list(TOXICITY_WEIGHTS)

    weighted = df[components].mul(
        pd.Series(TOXICITY_WEIGHTS),
        axis=1
    )

    df["Score_Toxicity"] = weighted.sum(
        axis=1,
        min_count=len(components)
    )

    # Worst individual toxicity desirability
    df["Score_Toxicity_Worst"] = df[
        components
    ].min(axis=1, skipna=False)

    # Maximum predicted toxicity liability
    df["Toxicity_Max_Probability"] = (
        1.0 - df["Score_Toxicity_Worst"]
    )

    # Number of endpoints with predicted
    # positive-class probability >= 0.5
    risk_columns = list(
        TOXICITY_ENDPOINTS.values()
    )

    risk_values = df[risk_columns].apply(
        pd.to_numeric, errors="coerce"
    )

    complete = risk_values.notna().all(axis=1)

    df["Toxicity_Alert_Count"] = (
        risk_values.ge(0.5)
        .sum(axis=1)
        .where(complete)
    )

    df["Toxicity_MissingComponents"] = df[
        components
    ].isna().sum(axis=1)

    df["Toxicity_Status"] = np.where(
        df["Toxicity_MissingComponents"].eq(0),
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

    output = calculate_toxicity_scores(df)

    if len(output) != len(df):
        raise RuntimeError("Unexpected row count change")

    path = Path(args.output)

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    output.to_csv(path, index=False)

    columns = list(TOXICITY_WEIGHTS) + [
        "Score_Toxicity",
        "Score_Toxicity_Worst",
        "Toxicity_Max_Probability",
        "Toxicity_Alert_Count"
    ]

    print("\nToxicity scoring coverage:")
    print(output[columns].count())

    print("\nToxicity score statistics:")
    print(output[columns].describe().round(3))

    print("\nToxicity status:")
    print(output["Toxicity_Status"].value_counts())

    print("\nToxicity alert distribution:")
    print(
        output["Toxicity_Alert_Count"]
        .value_counts(dropna=False)
        .sort_index()
    )

    summary = {
        "total_molecules": len(output),
        "weights": TOXICITY_WEIGHTS,
        "complete_scores": int(
            output["Score_Toxicity"].count()
        ),
        "mean_toxicity_score": float(
            output["Score_Toxicity"].mean()
        ),
        "toxicity_score_ge_0_7": int(
            output["Score_Toxicity"].ge(0.7).sum()
        ),
        "zero_toxicity_alerts": int(
            output["Toxicity_Alert_Count"].eq(0).sum()
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
