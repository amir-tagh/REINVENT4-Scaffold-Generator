#!/usr/bin/env python3

"""Synchronize computational optimization profiles with the website."""

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

SOURCE = ROOT / "config" / "optimization_profiles.json"
TARGET = ROOT / "docs" / "assets" / "data" / "optimization_profiles.json"

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
    "Toxicity",
]

with SOURCE.open() as handle:
    profiles = json.load(handle)

for name, weights in profiles.items():

    if set(weights) != set(COMPONENTS):
        raise ValueError(
            "Invalid components in profile: {}".format(name)
        )

    for component, weight in weights.items():
        if not isinstance(weight, (int, float)):
            raise ValueError(
                "Non-numeric weight: {} / {}".format(name, component)
            )

        if not math.isfinite(weight) or not 0 <= weight <= 1:
            raise ValueError(
                "Invalid weight: {} / {}".format(name, component)
            )

    total = sum(weights.values())

    if not math.isclose(total, 1.0, abs_tol=1e-8):
        raise ValueError(
            "{} weights sum to {}, expected 1.0".format(name, total)
        )

    print("[PASS] {:20s} Total: {:.2f}%".format(name, total * 100))

TARGET.parent.mkdir(parents=True, exist_ok=True)

with TARGET.open("w") as handle:
    json.dump(profiles, handle, indent=4)
    handle.write("\n")

print("\nSynchronized {} optimization profiles.".format(len(profiles)))
print("Website configuration: {}".format(TARGET))
