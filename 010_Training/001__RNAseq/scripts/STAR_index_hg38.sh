#!/bin/bash
#SBATCH --mem=250G

module load STAR/2.7.3a-GCC-9.3.0

nthreads=12


# --> UPDATE EACH PATH ACCORDINGLY!!!!!!!!!!!!!!!!!!! 


STAR --runThreadN 12 \
	--runMode genomeGenerate \
	--genomeDir /Master/meta/STAR_hg38 \
	--genomeFastaFiles /Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	--sjdbGTFfile /Master/meta/ENCFF159KBI.gtf

    
