#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "2dN_HET_H3K27me3_R1" 2.27
  "2dN_HET_H3K27me3_R2" 2.03
  "2dN_KO_H3K27me3_R1" 1.69
  "2dN_KO_H3K27me3_R2" 1.14
  "2dN_WT_H3K27me3_R1" 1.4
  "2dN_WT_H3K27me3_R2" 1.84
  "ESC_HET_H3K27me3_R1" 10.12
  "ESC_HET_H3K27me3_R2" 23.95
  "ESC_KO_H3K27me3_R1" 9.7
  "ESC_KO_H3K27me3_R2" 15.53
  "ESC_WT_H3K27me3_R1" 1.98
  "ESC_WT_H3K27me3_R2" 4.65
  "ESC_WT_H3K27me3_R3" 1
  "NPC_HET_H3K27me3_R1" 1.24
  "NPC_HET_H3K27me3_R2" 1.61
  "NPC_KO_H3K27me3_R1" 1.7
  "NPC_KO_H3K27me3_R2" 3.14
  "NPC_WT_H3K27me3_R1" 1.7
  "NPC_WT_H3K27me3_R2" 1.66
)

for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"

  libSize=$(samtools view -c -F 260 output/bowtie2_endtoend/${sample}.dupmark.sorted.bam)
  scale=$(echo "15000000/($libSize*$SF)" | bc -l)

  bamCoverage --bam output/bowtie2_endtoend/${sample}.dupmark.sorted.bam \
      --outFileName output/bigwig_ChIPseqSpikeInFree/${sample}.dupmark.sorted.bw \
      --outFileFormat bigwig \
      --binSize 10 \
      --numberOfProcessors 7 \
      --extendReads \
      --scaleFactor $scale
done