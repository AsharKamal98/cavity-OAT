#!/bin/bash
# Submit from this analysis directory with: bash submit.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p ../../logs

ARRAY_JOB_ID="$(sbatch --parsable submit.slurm)"
# Federated clusters may append ';cluster' to the numeric job ID.
ARRAY_JOB_ID="${ARRAY_JOB_ID%%;*}"
sbatch --dependency="afterok:${ARRAY_JOB_ID}" combine.slurm "${ARRAY_JOB_ID}"
echo "Array job: ${ARRAY_JOB_ID}"
echo "Combined artifact: output/${ARRAY_JOB_ID}/ensemble.pkl"
