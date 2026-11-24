#!/usr/bin/env python

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def decreasing_desirability(values, optimal, worst):

    values = pd.to_numeric(values, errors="coerce")
    values = values.replace([np.inf, -np.inf], np.nan)

    return ((worst - values) / (worst - optimal)).clip(0, 1)


def tpsa_desirability(values):

    values = pd.to_numeric(values, errors="coerce")
    values = values.replace([np.inf, -np.inf], np.nan)

    lower = ((values - 20.0) / 20.0).clip(0, 1)
    upper = ((120.0 - values) / 30.0).clip(0, 1)

    return pd.Series(
        np.minimum(lower, upper),
        index=values.index
    )


def calculate_cns_scores(df):

    df = df.copy()

    required = [
        "MW",
        "cLogP",
        "LogD_7_4_ADMETAI",
        "TPSA",
        "HBD",
        "pKa_Basic",
        "Num_Basic_Sites",
        "pKa_Status"
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    # Molecular weight
    df["Score_CNS_MW"] = decreasing_desirability(
        df["MW"], 360.0, 500.0
    )

    # Lipophilicity
    df["Score_CNS_LogP"] = decreasing_desirability(
        df["cLogP"], 3.0, 5.0
    )

    # Distribution coefficient at pH 7.4
    df["Score_CNS_LogD"] = decreasing_desirability(
        df["LogD_7_4_ADMETAI"], 2.0, 4.0
    )

    # Polar surface area
    df["Score_CNS_TPSA"] = tpsa_desirability(
        df["TPSA"]
    )

    # Hydrogen bond donors
    df["Score_CNS_HBD"] = decreasing_desirability(
        df["HBD"], 0.5, 3.5
    )

    # Basic pKa
    df["Score_CNS_pKa"] = decreasing_desirability(
        df["pKa_Basic"], 8.0, 10.0
    )

    # Distinguish absent detected basic sites
    # from actual prediction failures
    no_basic = (
        df["Num_Basic_Sites"].eq(0)
        & df["pKa_Status"].isin([
            "OK",
            "No_Ionizable_Sites"
        ])
    )

    df.loc[no_basic, "Score_CNS_pKa"] = 1.0

    failed = ~df["pKa_Status"].isin([
        "OK",
        "No_Ionizable_Sites"
    ])

    df.loc[failed, "Score_CNS_pKa"] = np.nan

    # Five-component CNS score
    partial_components = [
        "Score_CNS_MW",
        "Score_CNS_LogP",
        "Score_CNS_LogD",
        "Score_CNS_TPSA",
        "Score_CNS_HBD"
    ]

    df["CNS_Partial_5"] = df[
        partial_components
    ].sum(axis=1, min_count=5)

    df["Score_CNS_Partial"] = (
        df["CNS_Partial_5"] / 5.0
    )

    # Full six-component CNS MPO
    full_components = partial_components + [
        "Score_CNS_pKa"
    ]

    df["CNS_MPO"] = df[
        full_components
    ].sum(axis=1, min_count=6)

    df["Score_CNS_MPO"] = (
        df["CNS_MPO"] / 6.0
    )

    # Explicit quality and usage labels
    df["CNS_MPO_pKa_Quality"] = "Unvalidated_MolGpKa"

    df["CNS_Optimization_Score"] = (
        df["Score_CNS_Partial"]
    )

    df["CNS_Optimization_Method"] = (
        "Partial_5_Components_No_pKa"
    )

    # Conventional CNS MPO threshold
    df["CNS_MPO_GE4"] = (
        df["CNS_MPO"].ge(4.0)
        .where(df["CNS_MPO"].notna())
    )

    # Heuristic partial CNS threshold
    df["CNS_Partial_GE_0_8"] = (
        df["Score_CNS_Partial"].ge(0.8)
        .where(df["Score_CNS_Partial"].notna())
    )

    return df


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    df = pd.read_csv(args.input)

    print("Input molecules:", len(df))

    output = calculate_cns_scores(df)

    if len(output) != len(df):
        raise RuntimeError("Unexpected row count change")

    output_path = Path(args.output)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    output.to_csv(output_path, index=False)

    columns = [
        "Score_CNS_MW",
        "Score_CNS_LogP",
        "Score_CNS_LogD",
        "Score_CNS_TPSA",
        "Score_CNS_HBD",
        "Score_CNS_pKa",
        "CNS_MPO",
        "CNS_Partial_5",
        "Score_CNS_Partial"
    ]

    print("\nCNS prediction coverage:")
    print(output[columns].count())

    print("\nCNS score statistics:")
    print(output[columns].describe().round(3))

    full_pass = int(
        output["CNS_MPO"].ge(4.0).sum()
    )

    partial_pass = int(
        output["Score_CNS_Partial"].ge(0.8).sum()
    )

    print("\nCNS MPO >= 4:", full_pass)

    print(
        "Normalized partial CNS score >= 0.8:",
        partial_pass
    )

    summary = {
        "input": args.input,
        "output": args.output,
        "total_molecules": len(output),
        "full_cns_mpo_coverage": int(
            output["CNS_MPO"].count()
        ),
        "partial_cns_coverage": int(
            output["Score_CNS_Partial"].count()
        ),
        "mean_cns_mpo": float(
            output["CNS_MPO"].mean()
        ),
        "mean_partial_cns_score": float(
            output["Score_CNS_Partial"].mean()
        ),
        "cns_mpo_ge4": full_pass,
        "partial_cns_ge_0_8": partial_pass,
        "optimization_method": "Partial_5_Components",
        "pka_used_for_optimization": False,
        "pka_quality": "Unvalidated_MolGpKa"
    }

    summary_path = output_path.with_suffix(
        ".summary.json"
    )

    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=4)

    print("\nOutput:", output_path)
    print("Summary:", summary_path)


if __name__ == "__main__":
    main()
