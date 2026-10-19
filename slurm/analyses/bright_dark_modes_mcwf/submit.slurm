#!/bin/bash
#SBATCH --job-name=bright-dark-v1

#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00
#SBATCH --output=../../logs/%j_%x.out
#SBATCH --error=../../logs/%j_%x.err

# Submit from this analysis directory with: sbatch submit.slurm
# One ensemble using the bright_dark_modes_mcwf notebook parameters.

set -euo pipefail

source ~/.venvs/cavity-oat_env/bin/activate

ANALYSIS_DIR="${SLURM_SUBMIT_DIR}"
REPO_ROOT="$(cd "${ANALYSIS_DIR}/../../.." && pwd)"
OUTPUT_DIR="${ANALYSIS_DIR}/output"
mkdir -p "${OUTPUT_DIR}"

NTRAJ="500"
GAMMA="1.0"
NI="5,5"
DN="0"
NUM_SNAPSHOTS="500"
OMEGA_1="0.7"
OMEGA_FACTOR="0.1"
DELTA_FACTOR="0.02"
T1="10.0"
T2="10000.0"
T3="10.0"

FILENAME="${NI}${SLURM_JOB_NAME}"

python "${REPO_ROOT}/slurm/run_solvers/run_mcwf.py" \
  --n-processes "${SLURM_CPUS_PER_TASK}" \
  --ntraj "${NTRAJ}" \
  --filename "${FILENAME}" \
  --output-dir "${OUTPUT_DIR}" \
  --array-index "0" \
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
