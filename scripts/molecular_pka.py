#!/usr/bin/env python

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

import prop_profiler

from prop_profiler.models.molgpka_wrapper import MolGpKaWrapper
from prop_profiler.external.MolGpKa.predict_pka import model_pred
from prop_profiler.external.MolGpKa.utils.ionization_group import (
    get_ionization_aid
)


class MolecularPKA:

    def __init__(self, device="cpu"):

        weights = (
            Path(prop_profiler.__file__).parent
            / "model_weights"
        )

        self.device = device

        self.model = MolGpKaWrapper(
            acid_model_path=weights / "weight_acid.pth",
            base_model_path=weights / "weight_base.pth",
            device=device,
            verbose=False
        )

    def predict(self, smiles):

        result = {
            "pKa_Basic": np.nan,
            "pKa_Acidic_Min": np.nan,
            "pKa_Acidic_Max": np.nan,
            "Num_Basic_Sites": np.nan,
            "Num_Acidic_Sites": np.nan,
            "pKa_Basic_All": "",
            "pKa_Acidic_All": "",
            "pKa_Status": "OK",
            "pKa_Error": ""
        }

        try:
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                raise ValueError("Invalid SMILES")

            mol = self.model.standardize(mol)

            base_idxs = list(
                get_ionization_aid(
                    mol,
                    acid_or_base="base"
                )
            )

            acid_idxs = list(
                get_ionization_aid(
                    mol,
                    acid_or_base="acid"
                )
            )

            result["Num_Basic_Sites"] = len(base_idxs)
            result["Num_Acidic_Sites"] = len(acid_idxs)

            base_vals = [
                float(model_pred(
                    mol,
                    aid,
                    self.model.base_model,
                    device=self.device
                ))
                for aid in base_idxs
            ]

            acid_vals = [
                float(model_pred(
                    mol,
                    aid,
                    self.model.acid_model,
                    device=self.device
                ))
                for aid in acid_idxs
            ]

            if not all(np.isfinite(base_vals + acid_vals)):
                raise ValueError("Non-finite pKa prediction")

            if base_vals:
                result["pKa_Basic"] = max(base_vals)

            if acid_vals:
                result["pKa_Acidic_Min"] = min(acid_vals)
                result["pKa_Acidic_Max"] = max(acid_vals)

            result["pKa_Basic_All"] = ";".join(
                f"{v:.4f}" for v in base_vals
            )

            result["pKa_Acidic_All"] = ";".join(
                f"{v:.4f}" for v in acid_vals
            )

            if not base_vals and not acid_vals:
                result["pKa_Status"] = "No_Ionizable_Sites"

        except Exception as exc:
            result["pKa_Status"] = "Failed"
            result["pKa_Error"] = str(exc)

        return result


def basic_pka_desirability(row):

    if row["pKa_Status"] == "Failed":
        return np.nan

    if row["Num_Basic_Sites"] == 0:
        return 1.0

    pka = row["pKa_Basic"]

    if pd.isna(pka):
        return np.nan

    return float(
        np.clip((10.0 - pka) / 2.0, 0.0, 1.0)
    )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--smiles-column", default="SMILES")
    parser.add_argument(
        "--device",
        choices=["cpu", "cuda"],
        default="cpu"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.input)

    if args.smiles_column not in df.columns:
        raise ValueError("SMILES column not found")

    print("Input molecules:", len(df))
    print("Device:", args.device)

    predictor = MolecularPKA(device=args.device)

    cache = {}
    records = []

    for smiles in tqdm(
        df[args.smiles_column],
        total=len(df),
        desc="Predicting pKa"
    ):

        key = str(smiles)

        if key not in cache:
            cache[key] = predictor.predict(key)

        records.append(cache[key])

    predictions = pd.DataFrame(records)

    predictions["Score_CNS_pKa"] = predictions.apply(
        basic_pka_desirability,
        axis=1
    )

    output = pd.concat(
        [
            df.reset_index(drop=True),
            predictions.reset_index(drop=True)
        ],
        axis=1
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    output.to_csv(output_path, index=False)

    print("\nPrediction status:")
    print(output["pKa_Status"].value_counts())

    print("\nUnique input SMILES:", len(cache))
    print("Output:", output_path)


if __name__ == "__main__":
    main()
