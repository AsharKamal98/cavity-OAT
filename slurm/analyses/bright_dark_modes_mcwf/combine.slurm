#!/bin/bash
#SBATCH --job-name=bright-dark-combine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=../../logs/%j_%x.out
#SBATCH --error=../../logs/%j_%x.err

# submit.sh schedules this only after all four array tasks succeed.
# The memory request covers 4 * 50 trajectories and merge copies.
# For 4 * 250 trajectories at this saved-time resolution, request 192G.
set -euo pipefail

source ~/.venvs/cavity-oat_env/bin/activate

ANALYSIS_DIR="${SLURM_SUBMIT_DIR}"
REPO_ROOT="$(cd "${ANALYSIS_DIR}/../../.." && pwd)"
ARRAY_JOB_ID="${1:?Pass the simulation array job ID.}"

python "${REPO_ROOT}/slurm/combine_trajectory_ensembles.py" \
  --input-dir "${ANALYSIS_DIR}/output/${ARRAY_JOB_ID}" \
  --num-files 4
