#!/bin/bash
#SBATCH --job-name=REINVENT_pKa
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --partition=hpg-turin
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=pka_%j.out
#SBATCH --error=pka_%j.err

set -euo pipefail

cd /blue/mdisney/a.taghavi/REINVENT4

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

INPUT=${1:?Provide input CSV}
OUTPUT=${2:?Provide output CSV}

nvidia-smi

time srun /blue/mdisney/share/.conda/envs/cns_mpo/bin/python -u \
    scripts/molecular_pka.py \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --device cuda
