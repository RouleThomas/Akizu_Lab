#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=100:00:00


module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
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