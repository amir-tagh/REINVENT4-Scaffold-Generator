#!/bin/bash
#SBATCH --job-name=REINVENT_V6
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=reinvent_v6_%j.out
#SBATCH --error=reinvent_v6_%j.err

set -euo pipefail

# Repository location
# Run from the repository root when submitting through Slurm.
if [[ -n "${REINVENT_ROOT:-}" ]]; then
    ROOT="$REINVENT_ROOT"
elif [[ -n "${SLURM_JOB_ID:-}" ]]; then
    ROOT="${SLURM_SUBMIT_DIR:?Missing Slurm submission directory}"
else
    ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

cd "$ROOT"

# Load runtime configuration
RUNTIME_CONFIG="${RUNTIME_CONFIG:-config/runtime.conf}"

if [[ ! -f "$RUNTIME_CONFIG" ]]; then
    RUNTIME_CONFIG="config/runtime.example.conf"
fi

source "$RUNTIME_CONFIG"

resolve_python() {
    local env_name="$1"

    conda run --no-capture-output -n "$env_name" \
        python -c 'import sys; print(sys.executable)'
}

REINVENT_PY="${REINVENT_PY:-$(resolve_python "$REINVENT_ENV")}"
ADMET_PY="${ADMET_PY:-$(resolve_python "$ADMET_ENV")}"
PKA_PY="${PKA_PY:-$(resolve_python "$PKA_ENV")}"

# Input arguments
INPUT=${1:?Provide input SMILES file}
N=${2:-1000}
OUT=${3:-results/v6/REINVENT_V6_Results}
PROFILE=${4:-cns}
SELECT=${5:-30}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-${CPU_THREADS:-4}}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-${CPU_THREADS:-4}}
export PYTHONUNBUFFERED=1

mkdir -p "$OUT"
OUT=$(realpath "$OUT")

GEN="$OUT/generation"
ADMET="$OUT/admet"
CMOS="$OUT/cmos"

mkdir -p "$ADMET" "$CMOS"

stage() {
    echo
    echo "=================================================="
    echo "$1"
    echo "Time: $(date)"
    echo "=================================================="
}

trap 'echo "PIPELINE FAILED at line $LINENO" >&2' ERR

stage "REINVENT V6 PIPELINE INITIALIZATION"

echo "Input: $INPUT"
echo "Molecules per parent: $N"
echo "Output: $OUT"
echo "Optimization profile: $PROFILE"
echo "Selected per parent: $SELECT"

if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi || true
fi

# Validate dependencies
for PYTHON in "$REINVENT_PY" "$ADMET_PY" "$PKA_PY"
do
    test -x "$PYTHON"
done

test -s "$INPUT"
test -s "$LIBINVENT_PRIOR"
test -d "$SCSCORE_REPO"
test -s "$SCSCORE_MODEL"

"$REINVENT_PY" -c '
import json, sys
with open("config/optimization_profiles.json") as f:
    profiles = json.load(f)
weights = profiles[sys.argv[1]]
assert abs(sum(weights.values()) - 1) < 1e-8
print("Optimization weights:", weights)
' "$PROFILE"

# Save configuration
cp config/optimization_profiles.json \
    "$OUT/optimization_profiles.json"

# --------------------------------------------------
# Stage 1: Molecular generation
# --------------------------------------------------

stage "STAGE 1: REINVENT4 MOLECULAR GENERATION"

"$REINVENT_PY" -u \
    scripts/scaffold_libinvent_v5_devicefix.py \
    --input "$INPUT" \
    --num-molecules "$N" \
    --out "$GEN" \
    --model "$LIBINVENT_PRIOR" \
    --device "$REINVENT_DEVICE" \
    --qed-min 0.5 \
    --sa-max 4.0 \
    --sc-max 5.0 \
    --scscore-repo "$SCSCORE_REPO" \
    --scscore-model "$SCSCORE_MODEL" \
    --max-rounds 30 \
    --batch-size 512 \
    --min-samples 128

test -s "$GEN/combined_generated.csv"

# --------------------------------------------------
# Stage 2: Solubility and ADMET prediction
# --------------------------------------------------

stage "STAGE 2: ESOL AND ADMET-AI PREDICTION"

"$ADMET_PY" -u \
    scripts/molecular_admet.py \
    --input "$GEN/combined_generated.csv" \
    --output "$ADMET/ADMET_Predictions.csv" \
    --backend both \
    --batch-size 512

test -s "$ADMET/ADMET_Predictions.csv"

# --------------------------------------------------
# Stage 3: pKa prediction
# --------------------------------------------------

stage "STAGE 3: GPU-ACCELERATED PKA PREDICTION"

"$PKA_PY" -u \
    scripts/molecular_pka.py \
    --input "$ADMET/ADMET_Predictions.csv" \
    --output "$ADMET/ADMET_pKa.csv" \
    --device "$PKA_DEVICE"

test -s "$ADMET/ADMET_pKa.csv"

# --------------------------------------------------
# Stage 4: CNS scoring
# --------------------------------------------------

stage "STAGE 4: CNS DESIRABILITY SCORING"

"$REINVENT_PY" -u \
    scripts/molecular_cns.py \
    --input "$ADMET/ADMET_pKa.csv" \
    --output "$ADMET/ADMET_CNS.csv"

test -s "$ADMET/ADMET_CNS.csv"

# --------------------------------------------------
# Stage 5: ADME scoring
# --------------------------------------------------

stage "STAGE 5: ADME SCORING"

"$REINVENT_PY" -u \
    scripts/molecular_adme.py \
    --input "$ADMET/ADMET_CNS.csv" \
    --output "$ADMET/ADMET_CNS_ADME.csv"

test -s "$ADMET/ADMET_CNS_ADME.csv"

# --------------------------------------------------
# Stage 6: Toxicity scoring
# --------------------------------------------------

stage "STAGE 6: TOXICITY SCORING"

"$REINVENT_PY" -u \
    scripts/molecular_toxicity.py \
    --input "$ADMET/ADMET_CNS_ADME.csv" \
    --output "$ADMET/ADMET_Complete.csv"

test -s "$ADMET/ADMET_Complete.csv"

# --------------------------------------------------
# Stage 7: CMOS V2 optimization
# --------------------------------------------------

stage "STAGE 7: CMOS V2 MULTIPARAMETER OPTIMIZATION"

"$REINVENT_PY" -u \
    scripts/cumulative_score_v2.py \
    --input "$ADMET/ADMET_Complete.csv" \
    --output "$CMOS/CMOS_V2.csv" \
    --profiles "$OUT/optimization_profiles.json" \
    --profile "$PROFILE" \
    --select-per-hit "$SELECT"

test -s "$CMOS/CMOS_V2_selected.csv"

# --------------------------------------------------
# Stage 8: Final results
# --------------------------------------------------

stage "STAGE 8: FINAL LIBRARY SUMMARY"

"$REINVENT_PY" - "$OUT" "$SELECT" <<'PY'

import sys
from pathlib import Path
import pandas as pd
import numpy as np

root = Path(sys.argv[1])
target = int(sys.argv[2])

full = pd.read_csv(root / "cmos/CMOS_V2.csv")
selected = pd.read_csv(root / "cmos/CMOS_V2_selected.csv")

print("\nGenerated molecules:", len(full))
print("Selected molecules:", len(selected))
print("Unique selected structures:",
      selected["CMOS_V2_CanonicalSMILES"].nunique())

counts = selected.groupby("ParentID").size().reindex(
    sorted(full["ParentID"].unique()), fill_value=0
)

print("\nSelected molecules per parent:")
print(counts.to_string())

if not counts.eq(target).all():
    print("\nWARNING: Some parents did not reach selection target")

columns = [
    "QED",
    "SA",
    "SCScore",
    "MW",
    "cLogP",
    "TPSA",
    "HBD",
    "HBA",
    "ParentSimilarity",
    "LogS_ESOL",
    "LogS_ADMETAI",
    "LogS_Consensus",
    "LogD_7_4_ADMETAI",
    "Score_CNS_Partial",
    "Score_BBB",
    "Score_ADME",
    "Score_Toxicity",
    "Toxicity_Alert_Count",
    "CMOS_V2_Intrinsic"
]

summary = pd.DataFrame({
    "Full_Library": full[columns].mean(),
    "Selected": selected[columns].mean()
})

summary["Difference"] = (
    summary["Selected"] - summary["Full_Library"]
)

print("\nProperty comparison:")
print(summary.round(4).to_string())

summary.to_csv(root / "Final_Property_Comparison.csv")

log300 = np.log10(300e-6)

print("\nConsensus solubility geometric mean (uM):")

for name, df in [
    ("Full library", full),
    ("Selected", selected)
]:

    solubility = 1e6 * 10 ** df["LogS_Consensus"].mean()

    both300 = (
        df["LogS_ESOL"].ge(log300)
        & df["LogS_ADMETAI"].ge(log300)
    ).sum()

    print(f"{name}: {solubility:.2f} uM")
    print(f"Both models >= 300 uM: {both300}/{len(df)}")

print("\nZero predicted toxicity alerts:",
      selected["Toxicity_Alert_Count"].eq(0).sum())

print("\nTop-ranked compound per parent:")

top = selected[
    selected["CMOS_V2_SelectionRank"].eq(1)
]

print(top[[
    "ParentID",
    "GeneratedID",
    "QED",
    "SA",
    "SCScore",
    "LogS_Consensus",
    "Score_CNS_Partial",
    "Score_ADME",
    "Score_Toxicity",
    "CMOS_V2_Cumulative"
]].round(4).to_string(index=False))

PY

stage "REINVENT V6 PIPELINE COMPLETE"

echo "Results: $OUT"
echo "Full library: $CMOS/CMOS_V2.csv"
echo "Selected library: $CMOS/CMOS_V2_selected.csv"
echo "Selected SMILES: $CMOS/CMOS_V2_selected.smi"
