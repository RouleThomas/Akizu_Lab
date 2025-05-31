#!/bin/bash
#SBATCH --job-name=1_2_dimethylhydrazine_9cffdc696a32
#SBATCH --array=0-79
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --output=/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/logs/sbs90_%A_%a.out
#SBATCH --error=/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/logs/sbs90_%A_%a.err




# Parameters
SIGNATURE="1_2_dimethylhydrazine_9cffdc696a32"
N_LIST=(500 1000 2000 4000 8000 16000 32000 64000)
REPS=10

# Map SLURM_ARRAY_TASK_ID to mutation count and replicate
IDX=$SLURM_ARRAY_TASK_ID
N_IDX=$((IDX / REPS))
REP=$((IDX % REPS + 1))
MUTATION_COUNT=${N_LIST[$N_IDX]}

# Run the job from scripts/ directory context
python scripts/simulate_array_experimental.py \
  --signature $SIGNATURE \
  --n $MUTATION_COUNT \
  --rep $REP \
  --seed $IDX


