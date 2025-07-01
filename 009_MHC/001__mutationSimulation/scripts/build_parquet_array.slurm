#!/bin/bash
#SBATCH --job-name=exome_tab
#SBATCH --array=1-22
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/logs/tab_%A_%a.out
#SBATCH --error=/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/logs/tab_%A_%a.err


MUTSIM_HOME=/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/
CHR=${SLURM_ARRAY_TASK_ID}

python scripts/build_exome_table.py \
       --cds-fa  $MUTSIM_HOME/mutsim/gtf/cds.fa \
       --chrom   chr${CHR} \
       --output  $MUTSIM_HOME/mutsim/parquet/chr${CHR}.parquet

