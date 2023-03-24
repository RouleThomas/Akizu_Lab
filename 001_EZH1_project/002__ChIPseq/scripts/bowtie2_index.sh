#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=72:00:00


bowtie2-build ../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta ../../Master/meta/bowtie2_genome_dir/GRCh38