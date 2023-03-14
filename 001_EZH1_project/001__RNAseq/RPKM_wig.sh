#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=72:00:00

bamCoverage --bam output/STAR/fastp/4wN_WT_R1_Aligned.sortedByCoord.out.bam --outFileName output/temp/4wN_WT_R1_Aligned.sortedByCoord.out.bigwig --outFileFormat bigwig --normalizeUsing RPKM --binSize 10  
bamCoverage --bam output/STAR/fastp/4wN_KO_R1_Aligned.sortedByCoord.out.bam --outFileName output/temp/4wN_KO_R1_Aligned.sortedByCoord.out.bigwig --outFileFormat bigwig --normalizeUsing RPKM --binSize 10  
bamCoverage --bam output/STAR/fastp/4wN_HET_R1_Aligned.sortedByCoord.out.bam --outFileName output/temp/4wN_HET_R1_Aligned.sortedByCoord.out.bigwig --outFileFormat bigwig --normalizeUsing RPKM --binSize 10 