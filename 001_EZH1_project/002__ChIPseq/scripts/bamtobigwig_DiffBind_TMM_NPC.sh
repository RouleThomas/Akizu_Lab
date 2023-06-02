#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "NPC_HET_H3K27me3_R1" 3.242710387
  "NPC_HET_H3K27me3_R2" 1.310448705
  "NPC_KO_H3K27me3_R1" 1.829832193
  "NPC_KO_H3K27me3_R2" 1.2635567
  "NPC_WT_H3K27me3_R1" 1.425558409
  "NPC_WT_H3K27me3_R2" 1.556845489
)


for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"

  bamCoverage --bam output/bowtie2_endtoend/${sample}.dupmark.sorted.bam \
      --outFileName output/bigwig_DiffBind_TMM/${sample}.bw \
      --outFileFormat bigwig \
      --binSize 1 \
      --numberOfProcessors 7 \
      --extendReads \
      --scaleFactor $SF
done

