#!/bin/bash
#SBATCH --job-name=bright-dark-v2

#SBATCH --cpus-per-task=16
## SBATCH --mem=128G
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

NTRAJ="50"
GAMMA="1.0"
NI="5,5"
# dN=None: omit --dN to disable the total-NJ window.
SAVED_TIME_STEP="0.1"
OMEGA_1="0.6666666666666666"
OMEGA_FACTOR="0.119"
DELTA_FACTOR="0.283"
JOINT_PROBABILITY="0.40"
T1="100.0"
T2="10000.0"
T3="100.0"

# Match the notebook's saved-time spacing of 0.1 / Gamma.
NUM_SNAPSHOTS="$(python - "${T1}" "${T2}" "${T3}" "${SAVED_TIME_STEP}" "${GAMMA}" <<'PY'
import math
import sys

t1, t2, t3, saved_time_step, gamma = map(float, sys.argv[1:])
print(math.ceil((t1 + t2 + t3) / (saved_time_step / gamma)) + 1)
PY
)"

FILENAME="${NI}${SLURM_JOB_NAME}"

python "${REPO_ROOT}/slurm/run_solvers/run_mcwf.py" \
  --n-processes "${SLURM_CPUS_PER_TASK}" \
  --ntraj "${NTRAJ}" \
  --filename "${FILENAME}" \
  --output-dir "${OUTPUT_DIR}" \
  --array-index "0" \
  --Gamma "${GAMMA}" \
  --Ni "${NI}" \
  --joint-probability "${JOINT_PROBABILITY}" \
  --num-snapshots "${NUM_SNAPSHOTS}" \
  --omega-i "${OMEGA_1}" \
  --Omega-factor "${OMEGA_FACTOR}" \
  --delta-factor "${DELTA_FACTOR}" \
  --T1 "${T1}" \
  --T2 "${T2}" \
  --T3 "${T3}" \
  --save-ensemble
