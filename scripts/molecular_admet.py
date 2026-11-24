
#!/usr/bin/env python3

"""
Molecular ADMET prediction module for REINVENT4.

Backends:
    ESOL: RDKit-based aqueous solubility prediction.
    ADMET-AI: pretrained neural network predictions.
    Both: combined predictions and comparison.

All solubility values are expressed as log10(mol/L).
"""

import argparse
import json
from pathlib import Path
from importlib.metadata import version

import numpy as np
import pandas as pd

from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors


def predict_esol(mol):
    """RDKit-adapted ESOL regression."""

    if mol is None or mol.GetNumHeavyAtoms() == 0:
        return np.nan

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)

    aromatic_atoms = sum(
        atom.GetIsAromatic()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() > 1
    )

    ap = aromatic_atoms / mol.GetNumHeavyAtoms()

    return (
        0.26121066137801696
        - 0.7416739523408995 * logp
        - 0.0066138847738667125 * mw
        + 0.003451545565957996 * rb
        - 0.42624840441316975 * ap
    )


def solubility_score(logs):
    """Provisional continuous desirability score."""

    return np.clip((logs + 6.0) / 3.0, 0.0, 1.0)


def canonicalize(smiles):

    if not isinstance(smiles, str):
        return None

    mol = Chem.MolFromSmiles(smiles)

    if mol is None or mol.GetNumHeavyAtoms() == 0:
        return None

    return Chem.MolToSmiles(mol, isomericSmiles=True)


def calculate_esol(df):

    records = {}

    unique_smiles = df["ADMET_SMILES"].dropna().unique()

    for smiles in unique_smiles:

        mol = Chem.MolFromSmiles(smiles)

        logs = predict_esol(mol)
        mw = Descriptors.MolWt(mol)

        records[smiles] = {
            "LogS_ESOL": logs,
            "ESOL_Solubility_uM": 10.0 ** logs * 1e6,
            "ESOL_Solubility_mg_mL": 10.0 ** logs * mw,
            "Score_Solubility_ESOL": solubility_score(logs),
        }

    predictions = pd.DataFrame.from_dict(
        records, orient="index"
    )

    predictions = predictions.reindex(columns=[
        "LogS_ESOL",
        "ESOL_Solubility_uM",
        "ESOL_Solubility_mg_mL",
        "Score_Solubility_ESOL",
    ])

    df = df.join(predictions, on="ADMET_SMILES")

    print(f"ESOL: predicted {len(records)} unique molecules")

    return df


def calculate_admet_ai(df, batch_size=512):

    from admet_ai import ADMETModel

    smiles = (
        df["ADMET_SMILES"]
        .dropna()
        .unique()
        .tolist()
    )

    if not smiles:
        raise ValueError("No valid molecules for ADMET-AI")

    print("Loading ADMET-AI models...", flush=True)

    model = ADMETModel()

    batches = []

    for start in range(0, len(smiles), batch_size):

        batch = smiles[start:start + batch_size]

        preds = model.predict(smiles=batch)

        if not isinstance(preds, pd.DataFrame):
            raise TypeError("Unexpected ADMET-AI output type")

        if not preds.index.is_unique:
            raise ValueError("Duplicate prediction indices")

        missing = set(batch) - set(preds.index)

        if missing:
            raise ValueError(
                f"Missing ADMET-AI predictions: {len(missing)}"
            )

        batches.append(preds)

        print(
            f"ADMET-AI: "
            f"{min(start + batch_size, len(smiles))}"
            f"/{len(smiles)}",
            flush=True,
        )

    predictions = pd.concat(batches)

    required = [
        "Solubility_AqSolDB",
        "Lipophilicity_AstraZeneca",
    ]

    for col in required:

        if col not in predictions.columns:
            raise KeyError(
                f"Required ADMET-AI endpoint missing: {col}"
            )

    # Preserve every ADMET-AI endpoint.
    predictions = predictions.add_prefix("ADMETAI_")

    df = df.join(predictions, on="ADMET_SMILES")

    df["LogS_ADMETAI"] = pd.to_numeric(
        df["ADMETAI_Solubility_AqSolDB"],
        errors="coerce",
    )

    df["LogD_7_4_ADMETAI"] = pd.to_numeric(
        df["ADMETAI_Lipophilicity_AstraZeneca"],
        errors="coerce",
    )

    df["ADMETAI_Solubility_uM"] = (
        10.0 ** df["LogS_ADMETAI"] * 1e6
    )

    mw_lookup = {
        s: Descriptors.MolWt(Chem.MolFromSmiles(s))
        for s in smiles
    }

    mw = df["ADMET_SMILES"].map(mw_lookup)

    df["ADMETAI_Solubility_mg_mL"] = (
        10.0 ** df["LogS_ADMETAI"] * mw
    )

    df["Score_Solubility_ADMETAI"] = (
        solubility_score(df["LogS_ADMETAI"])
    )

    return df


def compare_predictions(df):

    df["LogS_Difference"] = (
        df["LogS_ESOL"] - df["LogS_ADMETAI"]
    )

    df["LogS_Absolute_Difference"] = (
        df["LogS_Difference"].abs()
    )

    df["LogS_Consensus"] = (
        df["LogS_ESOL"] + df["LogS_ADMETAI"]
    ) / 2.0

    df["Score_Solubility_Consensus"] = (
        solubility_score(df["LogS_Consensus"])
    )

    difference = df["LogS_Absolute_Difference"]

    df["LogS_Model_Agreement"] = np.select(
        [
            difference.isna(),
            difference <= 0.5,
            difference <= 1.0,
        ],
        [
            "Missing",
            "High",
            "Moderate",
        ],
        default="Low",
    )

    return df


def main():

    parser = argparse.ArgumentParser(
        description="Molecular ADMET and solubility prediction"
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument(
        "--backend",
        choices=["esol", "admetai", "both"],
        default="esol",
    )

    parser.add_argument(
        "--smiles-column",
        default="SMILES",
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=512,
    )

    args = parser.parse_args()

    if args.batch_size <= 0:
        parser.error("--batch-size must be positive")

    df = pd.read_csv(args.input)

    if args.smiles_column not in df.columns:
        raise ValueError(
            f"Missing SMILES column: {args.smiles_column}"
        )

    df["ADMET_SMILES"] = (
        df[args.smiles_column].map(canonicalize)
    )

    valid = df["ADMET_SMILES"].notna().sum()

    print(f"Input molecules: {len(df)}")
    print(f"Valid SMILES: {valid}")
    print(f"Backend: {args.backend}")

    if valid == 0:
        raise ValueError("No valid input molecules")

    if args.backend in ("esol", "both"):
        df = calculate_esol(df)

    if args.backend in ("admetai", "both"):
        df = calculate_admet_ai(
            df, batch_size=args.batch_size
        )

    if args.backend == "both":
        df = compare_predictions(df)

    output = Path(args.output)

    output.parent.mkdir(
        parents=True, exist_ok=True
    )

    df.to_csv(output, index=False)

    summary = {
        "input": str(args.input),
        "output": str(output),
        "backend": args.backend,
        "total_molecules": len(df),
        "valid_molecules": int(valid),
        "rdkit_version": rdBase.rdkitVersion,
        "pandas_version": pd.__version__,
        "admet_ai_version": (
            version("admet-ai")
            if args.backend in ("admetai", "both")
            else None
        ),
    }

    for col in [
        "LogS_ESOL",
        "LogS_ADMETAI",
        "LogD_7_4_ADMETAI",
        "LogS_Consensus",
        "LogS_Absolute_Difference",
    ]:

        if col in df.columns:

            values = df[col].dropna()

            if not values.empty:

                summary[col] = {
                    "count": int(len(values)),
                    "mean": float(values.mean()),
                    "median": float(values.median()),
                    "min": float(values.min()),
                    "max": float(values.max()),
                }

    if args.backend == "both":

        pair = df[
            ["LogS_ESOL", "LogS_ADMETAI"]
        ].dropna()

        pearson = pair.corr().iloc[0, 1]
        spearman = pair.corr(
            method="spearman"
        ).iloc[0, 1]

        summary["model_comparison"] = {
            "pearson": (
                float(pearson)
                if np.isfinite(pearson) else None
            ),
            "spearman": (
                float(spearman)
                if np.isfinite(spearman) else None
            ),
            "mean_absolute_difference": (
                float(
                    df["LogS_Absolute_Difference"].mean()
                )
                if not pair.empty else None
            ),
        }

    summary_path = output.with_suffix(".summary.json")

    summary_path.write_text(
        json.dumps(summary, indent=4)
    )

    print(json.dumps(summary, indent=4))
    print(f"COMPLETE: {output}")


if __name__ == "__main__":
    main()
