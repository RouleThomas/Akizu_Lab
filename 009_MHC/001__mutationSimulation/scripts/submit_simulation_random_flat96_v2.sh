#!/bin/bash
#SBATCH --job-name=sim_random_flat96
#SBATCH --array=0-79
#SBATCH --cpus-per-task=2
#SBATCH --mem=250G
#SBATCH --time=01:00:00
#SBATCH --output=/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/logs/random_flat96_v2_%A_%a.out
#SBATCH --error=/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/logs/random_flat96_v2_%A_%a.err

# Parameters
N_LIST=(500 1000 2000 4000 8000 16000 32000 64000)
REPS=10

# Map SLURM_ARRAY_TASK_ID to mutation count and replicate
IDX=$SLURM_ARRAY_TASK_ID
N_IDX=$((IDX / REPS))
REP=$((IDX % REPS + 1))
MUTATION_COUNT=${N_LIST[$N_IDX]}

# Create output directory
OUTDIR=results/random_v2/n_${MUTATION_COUNT}
mkdir -p $OUTDIR

# Set output file path
OUTFILE=results/random_v2/n_${MUTATION_COUNT}/rep_$(printf "%02d" $REP).annot.parquet

# Run the simulation
python scripts/simulate_random_flat96_v2.py \
  --n $MUTATION_COUNT \
  --rep $REP \
  --seed $IDX \
  --exome-dir parquet \
  --context-index indices/context96.pkl \
  --out $OUTFILE \
  --dbnsfp ref/dbNSFP5.1a_grch38.gz 

