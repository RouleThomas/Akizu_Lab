#!/bin/bash
#SBATCH --job-name=filtered_exp
#SBATCH --output=logs/exp_%A_%a.out
#SBATCH --error=logs/exp_%A_%a.err
#SBATCH --array=REPLACE_ME%228  # Placeholder; real array passed by submit script
#SBATCH --cpus-per-task=2
#SBATCH --mem=36G
#SBATCH --time=01:00:00
#SBATCH --partition=general

# Activate environment
source /etc/profile
module load anaconda3/2023.09.0
source activate mutsim

# Signature and parameter setup (⛳ change 1)
SIGNATURE_LIST=( $(</rs01/home/fontbuj/mutsim/signatures/all_experimental_sigs.txt) )
N_LIST=(500 1000 2000 4000 8000 16000 32000)
REPS=10

TOTAL_SIGNATURES=${#SIGNATURE_LIST[@]}
TOTAL_NS=${#N_LIST[@]}
TOTAL_JOBS=$((TOTAL_SIGNATURES * TOTAL_NS * REPS))

# Compute effective job index
BASE_OFFSET=${BASE_OFFSET:-0}
IDX=$((SLURM_ARRAY_TASK_ID + BASE_OFFSET))

SIGNATURE_IDX=$((IDX / (REPS * TOTAL_NS)))
N_IDX=$(( (IDX / REPS) % TOTAL_NS ))
REP=$((IDX % REPS + 1))

SIGNATURE=${SIGNATURE_LIST[$SIGNATURE_IDX]}
MUTATION_COUNT=${N_LIST[$N_IDX]}

echo "SLURM ID: $SLURM_ARRAY_TASK_ID (offset: $BASE_OFFSET → global ID: $IDX)"
echo "Signature: $SIGNATURE"
echo "Mutations: $MUTATION_COUNT"
echo "Replicate: $REP"
echo "Seed: $IDX"

# Run simulation (⛳ change 2: different sigfile, ⛳ change 3: new output path)
python /rs01/home/fontbuj/mutsim/scripts/simulate_sbs5_array.py \
  --signature "$SIGNATURE" \
  --n "$MUTATION_COUNT" \
  --rep "$REP" \
  --seed "$IDX" \
  --sigfile /rs01/home/fontbuj/mutsim/signatures/experimental_with_flat.txt \
  --outdir /rs01/home/fontbuj/mutsim/results_experimental

