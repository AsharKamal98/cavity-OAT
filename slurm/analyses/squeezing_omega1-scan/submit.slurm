#!/bin/bash
#SBATCH --job-name=v4
#SBATCH --array=0-3

#SBATCH --cpus-per-task=64
#SBATCH --time=00:05:00
#SBATCH --output=../../logs/%A_%a_%x.out
#SBATCH --error=../../logs/%A_%a_%x.err

# Submit from this analysis directory with: sbatch submit.slurm
# The SLURM job name prefixes every saved array file.

set -euo pipefail

source ~/.venvs/cavity-oat_env/bin/activate

ANALYSIS_DIR="${SLURM_SUBMIT_DIR}"
REPO_ROOT="$(cd "${ANALYSIS_DIR}/../../.." && pwd)"
OUTPUT_DIR="${ANALYSIS_DIR}/output"
mkdir -p "${OUTPUT_DIR}"

NTRAJ="1000"
GAMMA="1.0"
NI="25,25"
DN="2"
NUM_SNAPSHOTS="500"
OMEGA_1_LIST=("1.0" "0.9" "0.8" "0.7")
OMEGA_FACTOR="0.1"
DELTA_FACTOR="0.02"
T1="10.0"
T2="100.0"
T3="10.0"

OMEGA_1="${OMEGA_1_LIST[${SLURM_ARRAY_TASK_ID}]}"
FILENAME="${SLURM_JOB_NAME}_omega_${OMEGA_1}"

python "${REPO_ROOT}/slurm/run_solvers/run_mcwf.py" \
  --n-processes "${SLURM_CPUS_PER_TASK}" \
  --ntraj "${NTRAJ}" \
  --filename "${FILENAME}" \
  --output-dir "${OUTPUT_DIR}" \
  --array-index "${SLURM_ARRAY_TASK_ID}" \
  --Gamma "${GAMMA}" \
  --Ni "${NI}" \
  --dN "${DN}" \
  --num-snapshots "${NUM_SNAPSHOTS}" \
  --omega-i "${OMEGA_1}" \
  --Omega-factor "${OMEGA_FACTOR}" \
  --delta-factor "${DELTA_FACTOR}" \
  --T1 "${T1}" \
  --T2 "${T2}" \
  --T3 "${T3}" \
  --save-ensemble
