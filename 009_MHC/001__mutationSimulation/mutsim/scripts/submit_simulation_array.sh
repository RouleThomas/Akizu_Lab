#!/bin/bash
#SBATCH --job-name=sim_sbs90
#SBATCH --output=logs/sbs90_%A_%a.out
#SBATCH --error=logs/sbs90_%A_%a.err
#SBATCH --array=0-79
#SBATCH --cpus-per-task=2
#SBATCH --mem=36G
#SBATCH --time=01:00:00

# Activate environment
source /etc/profile
module load anaconda3/2023.09.0
source activate mutsim

# Parameters
SIGNATURE="SBS90"
N_LIST=(500 1000 2000 4000 8000 16000 32000 64000)
REPS=10

# Map SLURM_ARRAY_TASK_ID to mutation count and replicate
IDX=$SLURM_ARRAY_TASK_ID
N_IDX=$((IDX / REPS))
REP=$((IDX % REPS + 1))
MUTATION_COUNT=${N_LIST[$N_IDX]}

# Run the job from scripts/ directory context
python simulate_sbs5_array.py \
  --signature $SIGNATURE \
  --n $MUTATION_COUNT \
  --rep $REP \
  --seed $IDX
