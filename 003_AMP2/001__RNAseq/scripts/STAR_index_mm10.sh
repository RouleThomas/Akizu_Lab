#!/bin/bash
#SBATCH --mem=500G

module load STAR/2.7.3a-GCC-9.3.0

nthreads=12

STAR --runThreadN 12 \
	--runMode genomeGenerate \
	--genomeDir /scr1/users/roulet/Akizu_Lab/Master/meta_mice/STAR_mm10 \
	--genomeFastaFiles /scr1/users/roulet/Akizu_Lab/Master/meta_mice/mm10_no_alt_analysis_set_ENCODE.fasta \
	--sjdbGTFfile /scr1/users/roulet/Akizu_Lab/Master/meta_mice/ENCFF871VGR.gtf

    
