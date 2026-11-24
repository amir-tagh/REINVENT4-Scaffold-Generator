#!/bin/bash
#SBATCH --job-name=REINVENT4_batch
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --partition=hpg-turin
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=reinvent_v4_%j.out
#SBATCH --error=reinvent_v4_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
INPUT=${1:?Usage: sbatch run_multiple_hits_v4.sh input.txt [count] [output] [extra Python options...]}
N=${2:-1000}
OUT=${3:-REINVENT_V4_Results}
if (( $# >= 3 )); then shift 3; else shift "$#"; fi
/blue/mdisney/a.taghavi/reinvent4_env/bin/python -u scripts/scaffold_libinvent_v5_devicefix.py \
    --input "$INPUT" --num-molecules "$N" --out "$OUT" \
    --model priors/libinvent.prior --device cuda:0 \
    --qed-min 0.7 --sa-max 4.0 --sc-max 3.5 \
    --scscore-repo ./scscore \
    --scscore-model scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz \
    --max-rounds 30 --batch-size 512 --min-samples 128 "$@"
