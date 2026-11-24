#!/bin/bash
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"
INPUT="${1:-inputs/Tau-hits.txt}"
N="${2:-1000}"
OUT="${3:-results/v5/REINVENT_V5_Relaxed_$(date +%Y%m%d_%H%M%S)}"
[[ -f "$INPUT" ]] || { echo "Missing input: $INPUT" >&2; exit 1; }
[[ -f scripts/run_multiple_hits_v4_standardized.sh ]] || { echo 'Missing generation launcher' >&2; exit 1; }
[[ -f scripts/run_reinvent_analysis.sh ]] || { echo 'Missing analysis launcher' >&2; exit 1; }
[[ ! -e "$OUT" ]] || { echo "Output already exists: $OUT" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
JOB=$(sbatch --parsable scripts/run_multiple_hits_v4_standardized.sh \
    "$INPUT" "$N" "$OUT" \
    --qed-min 0.5 --sa-max 4.0 --sc-max 5.0 \
    --max-rounds 15 --batch-size 512 --min-samples 128)
JOB="${JOB%%;*}"
ANALYSIS_JOB=$(sbatch --parsable --dependency="afterok:${JOB}" \
    scripts/run_reinvent_analysis.sh "$OUT" "${OUT}_Analysis")
echo "Generation job: $JOB"
echo "Analysis job (starts after successful generation): ${ANALYSIS_JOB%%;*}"
echo "Generated molecules: $OUT"
echo "Analysis report: ${OUT}_Analysis/report.html"
