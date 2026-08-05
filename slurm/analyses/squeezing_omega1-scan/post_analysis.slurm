#!/bin/bash
#SBATCH --job-name=pa
#SBATCH --array=0-3
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G

#SBATCH --time=00:10:00
#SBATCH --output=../../logs/%A_%a_%x.out
#SBATCH --error=../../logs/%A_%a_%x.err

# Submit this job with an afterok dependency on the simulation array job:
# sbatch --dependency=afterok:ARRAY_JOB_ID \
#   post_analysis.slurm ARRAY_JOB_ID

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch post_analysis.slurm ARRAY_JOB_ID" >&2
  exit 2
fi

source ~/.venvs/cavity-oat_env/bin/activate

ANALYSIS_DIR="${SLURM_SUBMIT_DIR}"
RUN_ID="$1"
INPUT_DIR="/scratch-shared/${USER}/cavity-OAT/squeezing_omega1-scan/${RUN_ID}"
OUTPUT_DIR="${ANALYSIS_DIR}/output"
NUM_FILES="4"
OMEGA_1_LIST=("1.0" "0.9" "0.8" "0.7")
OMEGA_1="${OMEGA_1_LIST[${SLURM_ARRAY_TASK_ID}]}"

python "${ANALYSIS_DIR}/post_analysis.py" \
  --filename "v1_omega_${OMEGA_1}" \
  --input-dir "${INPUT_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --num-files "${NUM_FILES}" \
  --n-processes "${SLURM_CPUS_PER_TASK}"
