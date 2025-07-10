#!/bin/bash

SLURM_SCRIPT="run_filtered_contexts.slurm"
CONCURRENT_LIMIT=228
MAX_BATCH_SIZE=1001

START=0
END=6719  # 96 contexts × 7 Ns × 10 replicates

for ((current=START; current<=END; current+=MAX_BATCH_SIZE)); do
  batch_start=$current
  batch_end=$((current + MAX_BATCH_SIZE - 1))
  if [ $batch_end -gt $END ]; then
    batch_end=$END
  fi
  array_size=$((batch_end - batch_start))

  echo "🚀 Submitting jobs $batch_start-$batch_end as array 0-$array_size"
  BASE_OFFSET=$batch_start sbatch --export=ALL,BASE_OFFSET=$batch_start \
    --array=0-$array_size%$CONCURRENT_LIMIT "$SLURM_SCRIPT"
done

echo "🎉 All batches submitted!"

