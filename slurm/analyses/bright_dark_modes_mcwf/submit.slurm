#!/bin/bash
#SBATCH --job-name=bright-dark-v2
#SBATCH --array=0-3

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
## SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=../../logs/%A_%a_%x.out
#SBATCH --error=../../logs/%A_%a_%x.err

# Submit the array and its dependent combine job with: bash submit.sh
# All four tasks use the same physical parameters.

set -euo pipefail

source ~/.venvs/cavity-oat_env/bin/activate

ANALYSIS_DIR="${SLURM_SUBMIT_DIR}"
REPO_ROOT="$(cd "${ANALYSIS_DIR}/../../.." && pwd)"
OUTPUT_DIR="${ANALYSIS_DIR}/output/${SLURM_ARRAY_JOB_ID}"
mkdir -p "${OUTPUT_DIR}"

NTRAJ="50"  # Trajectories per array task; the combined total is 4 * NTRAJ.
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

python "${REPO_ROOT}/slurm/run_solvers/run_mcwf.py" \
  --n-processes "${SLURM_CPUS_PER_TASK}" \
  --ntraj "${NTRAJ}" \
  --filename "ensemble" \
  --output-dir "${OUTPUT_DIR}" \
  --array-index "${SLURM_ARRAY_TASK_ID}" \
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
