#!/bin/bash
#SBATCH --job-name=nj-ratio-v1
#SBATCH --array=0-4

#SBATCH --cpus-per-task=64
#SBATCH --time=00:10:00
#SBATCH --output=../../logs/%A_%a_%x.out
#SBATCH --error=../../logs/%A_%a_%x.err

# Submit from this analysis directory with: sbatch submit.slurm
# Each array task changes only the constraint NJ1 = A * NJ2.

set -euo pipefail

source ~/.venvs/cavity-oat_env/bin/activate

ANALYSIS_DIR="${SLURM_SUBMIT_DIR}"
REPO_ROOT="$(cd "${ANALYSIS_DIR}/../../.." && pwd)"
OUTPUT_DIR="${ANALYSIS_DIR}/output"
mkdir -p "${OUTPUT_DIR}"

NTRAJ="1000"
GAMMA="1.0"
NI="50,50"
MAX_OMEGA_OVER_OMEGA_C="0.8"
JOINT_PROBABILITY="0.6"
NUM_SNAPSHOTS="500"
OMEGA_1="0.7"
OMEGA_FACTOR="0.1"
DELTA_FACTOR="0.02"
T1="10.0"
T2="100.0"
T3="10.0"
A_VALUES=("1" "1.5" "2" "3" "4")
CONSTRAINT_ROWS=("1,-1" "2,-3" "1,-2" "1,-3" "1,-4")

A="${A_VALUES[${SLURM_ARRAY_TASK_ID}]}"
LINEAR_CONSTRAINT="${CONSTRAINT_ROWS[${SLURM_ARRAY_TASK_ID}]}=0"
FILENAME="${SLURM_JOB_NAME}_A_${A}"

python "${REPO_ROOT}/slurm/run_solvers/run_mcwf.py" \
  --n-processes "${SLURM_CPUS_PER_TASK}" \
  --ntraj "${NTRAJ}" \
  --filename "${FILENAME}" \
  --output-dir "${OUTPUT_DIR}" \
  --array-index "${SLURM_ARRAY_TASK_ID}" \
  --Gamma "${GAMMA}" \
  --Ni "${NI}" \
  --max-omega-over-omega-c "${MAX_OMEGA_OVER_OMEGA_C}" \
  --linear-constraints "${LINEAR_CONSTRAINT}" \
  --joint-probability "${JOINT_PROBABILITY}" \
  --num-snapshots "${NUM_SNAPSHOTS}" \
  --omega-i "${OMEGA_1}" \
  --Omega-factor "${OMEGA_FACTOR}" \
  --delta-factor "${DELTA_FACTOR}" \
  --T1 "${T1}" \
  --T2 "${T2}" \
  --T3 "${T3}" \
  --save-ensemble
