#!/bin/bash
#SBATCH --job-name=mcwf-j
#SBATCH --array=0-3
#SBATCH --cpus-per-task=64
#SBATCH --time=02:30:00
#SBATCH --output=slurm/logs/%A_%a_%x.out
#SBATCH --error=slurm/logs/%A_%a_%x.err

## TIME GUIDELINE
## N1=N2=10, 3 min / processs / trajectory
## One two-group run, T1=T2=50 000 NTRAJ = 250, 64 CPUs for above took 12 minutes

## Vary N, scale drive parameters as N Gamma (decay rate scaling). How well do solutions agree with mean-field.
## how large systems can we simulate

set -euo pipefail

source ~/.venvs/cavity-oat_env/bin/activate

OUTPUT_DIR="slurm/finite_size_effects/outputs"

NTRAJ="250"
GAMMA="1.0"
NI_LIST=("10,10" "15,15" "20,20")
dN="0"
DT="0.01"
NUM_SNAPSHOTS="500"
OMEGA_1="0.8"
OMEGA_FACTOR="0.1"
DELTA_FACTOR="0.01"
T1="50000.0"
T2="50000.0"
T3="1000.0"


for NI in "${NI_LIST[@]}"; do
  FILENAME=two_group_Ni_${NI//,/x}

  python slurm/run_solvers/run_mcwf_two_group.py \
    --n-processes "${SLURM_CPUS_PER_TASK}" \
    --ntraj "${NTRAJ}" \
    --filename "${FILENAME}" \
    --output-dir "${OUTPUT_DIR}" \
    --array-index "${SLURM_ARRAY_TASK_ID}" \
    --Gamma "${GAMMA}" \
    --Ni "${NI}" \
    --dN "${dN}" \
    --dt "${DT}" \
    --num-snapshots "${NUM_SNAPSHOTS}" \
    --omega-i "${OMEGA_1}" \
    --Omega-factor "${OMEGA_FACTOR}" \
    --delta-factor "${DELTA_FACTOR}" \
    --T1 "${T1}" \
    --T2 "${T2}" \
    --T3 "${T3}"

  FILENAME=single_group_Ni_${NI//,/x}
  python slurm/run_solvers/run_mcwf_single_group.py \
    --n-processes "${SLURM_CPUS_PER_TASK}" \
    --ntraj "${NTRAJ}" \
    --filename "${FILENAME}" \
    --output-dir "${OUTPUT_DIR}" \
    --array-index "${SLURM_ARRAY_TASK_ID}" \
    --Gamma "${GAMMA}" \
    --Ni "${NI}" \
    --dN "${dN}" \
    --dt "${DT}" \
    --num-snapshots "${NUM_SNAPSHOTS}" \
    --Omega-factor "${OMEGA_FACTOR}" \
    --delta-factor "${DELTA_FACTOR}" \
    --T1 "${T1}" \
    --T2 "${T2}" \
    --T3 "${T3}"
done
