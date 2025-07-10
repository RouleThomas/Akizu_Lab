#!/bin/bash

SLURM_SCRIPT="run_filtered_experimental.slurm"  # 🔄 Experimental-specific SLURM script
CONCURRENT_LIMIT=228

START=0
END=7139  # 🔄 Total number of jobs for experimental run
MAX_BATCH_SIZE=1001

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

