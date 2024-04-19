#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


bamCoverage --bam output/bowtie2/CPC_untreated_YAP1_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_YAP1_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 127 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/CPC_RA_YAP1_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_YAP1_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_untreated_YAP1_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_YAP1_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 136 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398



bamCoverage --bam output/bowtie2/CPC_RA_YAP1_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_YAP1_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_untreated_TEAD4_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_TEAD4_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 122 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/CPC_RA_TEAD4_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_TEAD4_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 123 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/CPC_untreated_TEAD4_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_TEAD4_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 132 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_RA_TEAD4_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_TEAD4_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 130 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_untreated_NR2F2_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_NR2F2_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 130 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_RA_NR2F2_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_NR2F2_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 153 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_untreated_NR2F2_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_NR2F2_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 136 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/CPC_RA_NR2F2_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_NR2F2_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 136 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_RA_NR2F2_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_NR2F2_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 130 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/CPC_untreated_NR2F2_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_NR2F2_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/CPC_untreated_input_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_input_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 81 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_RA_input_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_input_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 81 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_untreated_YAP1_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_YAP1_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 84 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/CPC_RA_YAP1_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_YAP1_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 142 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/CPC_untreated_TEAD4_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_untreated_TEAD4_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 143 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/CPC_RA_TEAD4_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_RPGC/CPC_RA_TEAD4_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 124 \
    --normalizeUsing RPGC \
    --ignoreForNormalization chrX chrY chrM \
    --effectiveGenomeSize 2913022398


