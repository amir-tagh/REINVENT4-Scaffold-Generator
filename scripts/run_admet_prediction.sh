#!/bin/bash
#SBATCH --job-name=REINVENT_ADMET
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --partition=hpg-turin
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=admet_%j.out
#SBATCH --error=admet_%j.err

set -euo pipefail

cd /blue/mdisney/a.taghavi/REINVENT4

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

INPUT=${1:?Provide input CSV}
OUTPUT=${2:?Provide output CSV}
BACKEND=${3:-both}

nvidia-smi

srun /blue/mdisney/a.taghavi/admet_ai_env/bin/python -u \
    scripts/molecular_admet.py \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --backend "$BACKEND" \
    --batch-size 512
