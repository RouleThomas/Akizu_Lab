#!/bin/bash
#SBATCH --job-name=exome_tab
#SBATCH --array=1-22
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=24:00:00


MUTSIM_HOME=/scr1/users/roulet/Akizu_Lab/009_MHC/007__mutationSimulation
CHR=${SLURM_ARRAY_TASK_ID}

python scripts/build_exome_table_gtf_genome.py \
  --gtf      $MUTSIM_HOME/gtf/gencode.v49.annotation.gtf.gz \
  --genome   $MUTSIM_HOME/ref/GRCh38.primary_assembly.genome.fa \
  --chrom    chr${CHR} \
  --output   $MUTSIM_HOME/parquet/chr${CHR}.parquet
