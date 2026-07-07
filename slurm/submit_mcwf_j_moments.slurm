#!/bin/bash
#SBATCH --job-name=mcwf-j
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --time=00:10:00
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err

set -euo pipefail

source ~/.venvs/cavity-oat_env/bin/activate

FILENAME="mcwf_run"
NTRAJ="100"

python -m slurm.mcwf_j_moments_run \
  --n-processes "${SLURM_CPUS_PER_TASK}" \
  --ntraj "${NTRAJ}" \
  --filename "${FILENAME}" \
  --array-index "${SLURM_ARRAY_TASK_ID}"
