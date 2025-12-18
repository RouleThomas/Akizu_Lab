#!/bin/bash
#SBATCH --job-name=filtered
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --array=0-5459%50  # 5460 total jobs= "$TOTAL_SIGNATURES * $TOTAL_NS * $REPS"; will run 50 at same time # 78*7*10
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --partition=defq # default long run time check name with `sinfo`



# Signature and parameter setup 
SIGNATURE_LIST=( $(<signatures/context_signature_list_DBS.txt) )
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

python scripts/simulate_array_DBS_v2.py \
  --signature "$SIGNATURE" \
  --n "$MUTATION_COUNT" \
  --rep "$REP" \
  --seed "$IDX" \
  --sigfile signatures/context_sigs_fixed_DBS.txt \
  --contexts_list signatures/context_signature_list_DBS.txt \
  --context indices/context78.pkl \
  --exome parquet \
  --outdir results_contexts_DBS



