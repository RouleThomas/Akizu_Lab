#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


bamCoverage --bam output/bowtie2/CPC_untreated_YAP1_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_YAP1_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 127 \
    --scaleFactor 0.5

bamCoverage --bam output/bowtie2/CPC_RA_YAP1_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_YAP1_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_untreated_YAP1_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_YAP1_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 136 \
    --scaleFactor 0.5



bamCoverage --bam output/bowtie2/CPC_RA_YAP1_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_YAP1_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_untreated_TEAD4_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_TEAD4_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 122 \
    --scaleFactor 0.5

bamCoverage --bam output/bowtie2/CPC_RA_TEAD4_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_TEAD4_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 123 \
    --scaleFactor 0.5

bamCoverage --bam output/bowtie2/CPC_untreated_TEAD4_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_TEAD4_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 132 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_RA_TEAD4_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_TEAD4_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 130 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_untreated_NR2F2_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_NR2F2_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 130 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_RA_NR2F2_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_NR2F2_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 153 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_untreated_NR2F2_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_NR2F2_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 136 \
    --scaleFactor 0.5

bamCoverage --bam output/bowtie2/CPC_RA_NR2F2_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_NR2F2_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 136 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_RA_NR2F2_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_NR2F2_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 130 \
    --scaleFactor 0.5

bamCoverage --bam output/bowtie2/CPC_untreated_NR2F2_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_NR2F2_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --scaleFactor 0.5

bamCoverage --bam output/bowtie2/CPC_untreated_input_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_input_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 81 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_RA_input_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_input_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 81 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_untreated_YAP1_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_YAP1_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 84 \
    --scaleFactor 0.5

bamCoverage --bam output/bowtie2/CPC_RA_YAP1_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_YAP1_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 142 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/CPC_untreated_TEAD4_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_untreated_TEAD4_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 143 \
    --scaleFactor 0.5

bamCoverage --bam output/bowtie2/CPC_RA_TEAD4_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/CPC_RA_TEAD4_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 124 \
    --scaleFactor 0.5


