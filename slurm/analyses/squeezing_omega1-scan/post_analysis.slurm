#!/bin/bash
#SBATCH --job-name=mcwf-squeeze-analysis
#SBATCH --cpus-per-task=64
#SBATCH --mem=16G
#SBATCH --time=02:30:00
#SBATCH --output=slurm/logs/%j_%x.out
#SBATCH --error=slurm/logs/%j_%x.err

# Submit this job with an afterok dependency on the simulation array job:
# sbatch --dependency=afterok:ARRAY_JOB_ID \
#   slurm/analyses/squeezing_omega1-scan/post_analysis.slurm ARRAY_JOB_ID

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch post_analysis.slurm ARRAY_JOB_ID" >&2
  exit 2
fi

source ~/.venvs/cavity-oat_env/bin/activate

RUN_ID="$1"
INPUT_DIR="/scratch-shared/${USER}/cavity-OAT/squeezing_omega1-scan/${RUN_ID}"
OUTPUT_DIR="slurm/analyses/squeezing_omega1-scan/output"
NUM_FILES="4"
OMEGA_1_LIST=("1.0" "0.9" "0.8" "0.7")

for OMEGA_1 in "${OMEGA_1_LIST[@]}"; do
  python slurm/analyses/squeezing_omega1-scan/post_analysis.py \
    --filename "v1_omega_${OMEGA_1}" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --num-files "${NUM_FILES}" \
    --n-processes "${SLURM_CPUS_PER_TASK}"
done
