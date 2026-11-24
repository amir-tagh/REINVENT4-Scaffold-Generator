#!/bin/bash
#SBATCH --job-name=REINVENT_analysis
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --partition=hpg-default
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=reinvent_analysis_%j.out
#SBATCH --error=reinvent_analysis_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
INPUT="${1:?Usage: sbatch scripts/run_reinvent_analysis.sh input_dir output_dir}"
OUTPUT="${2:?Missing output_dir}"
/blue/mdisney/a.taghavi/reinvent4_env/bin/python -u analysis/analyze_reinvent_comprehensive.py --input "$INPUT" --output "$OUTPUT"
