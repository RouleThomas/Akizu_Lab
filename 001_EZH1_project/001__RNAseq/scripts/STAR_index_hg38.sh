#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=250G

module load STAR/2.7.3a-GCC-9.3.0

nthreads=12

STAR --runThreadN 12 \
	--runMode genomeGenerate \
	--genomeDir /scr1/users/roulet/Akizu_Lab/Master/meta/STAR_hg38 \
	--genomeFastaFiles /scr1/users/roulet/Akizu_Lab/Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	--sjdbGTFfile /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf

    
