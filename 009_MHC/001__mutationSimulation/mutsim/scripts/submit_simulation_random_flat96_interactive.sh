#!/bin/bash

# Parameters
N_LIST=(500 1000 2000 4000 8000 16000 32000 64000)
REPS=10

# Loop through all combinations of N and REP
for ((n_idx=0; n_idx<${#N_LIST[@]}; n_idx++)); do
  for ((rep=1; rep<=REPS; rep++)); do

    MUTATION_COUNT=${N_LIST[$n_idx]}
    IDX=$((n_idx * REPS + rep - 1))

    # Set output file path
    OUTFILE=results/random_interactive/n_${MUTATION_COUNT}/rep_$(printf "%02d" $rep).annot.parquet
    mkdir -p $(dirname "$OUTFILE")

    echo "Running simulation: n=$MUTATION_COUNT, rep=$rep, seed=$IDX"

    python scripts/simulate_random_flat96.py \
      --n $MUTATION_COUNT \
      --rep $rep \
      --seed $IDX \
      --exome-dir parquet \
      --context-index indices/context96.pkl \
      --out $OUTFILE

  done
done
