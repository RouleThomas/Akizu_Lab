#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "8wN_HET_IGG_R1" 1.3153942428035
  "8wN_HET_H3K27me3_R1" 1.50803188228081
  "8wN_HET_IGG_R2" 1.07884856070088
  "8wN_HET_H3K27me3_R2" 1.79362354383814
  "8wN_KO_IGG_R1" 1.3768115942029
  "8wN_KO_H3K27me3_R1" 1.30657216494845
  "8wN_KO_IGG_R2" 1
  "8wN_KO_H3K27me3_R2" 1
  "8wN_iPSCpatient_IGG_R1" 1.36049107142857
  "8wN_iPSCpatient_H3K27me3_R1" 1.5370709981661
  "8wN_WT_IGG_R1" 1.82217343578485
  "8wN_WT_H3K27me3_R1" 1
  "8wN_WT_IGG_R2" 2.33479692645445
  "8wN_WT_H3K27me3_R2" 1.73020349058761
  "8wN_HET_IGG_R3" 1
  "8wN_HET_H3K27me3_R3" 1
  "8wN_HET_IGG_R4" 2.43178973717146
  "8wN_HET_H3K27me3_R4" 2.71894543225015
  "8wN_KO_IGG_R3" 1.42555994729908
  "8wN_KO_H3K27me3_R3" 3.80953608247423
  "8wN_KO_IGG_R4" 1.02766798418972
  "8wN_KO_H3K27me3_R4" 1.12899484536082
  "8wN_iPSCpatient_IGG_R2" 1
  "8wN_iPSCpatient_H3K27me3_R2" 1
  "8wN_WT_IGG_R3" 1
  "8wN_WT_H3K27me3_R3" 1.21332215532353
  "8wN_WT_IGG_R4" 4.34577387486279
  "8wN_WT_H3K27me3_R4" 1.82966237329472
)

for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"

  libSize=$(samtools view -c -F 260 output/bowtie2/${sample}.dupmark.sorted.bam)
  scale=$(echo "15000000/($libSize*$SF)" | bc -l)

  bamCoverage --bam output/bowtie2/${sample}.dupmark.sorted.bam \
      --outFileName output/bigwig_histone_lib/${sample}.dupmark.sorted.bw \
      --outFileFormat bigwig \
      --binSize 10 \
      --numberOfProcessors 7 \
      --extendReads \
      --scaleFactor $scale
done