#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "NPC_HET_H3K27me3_R1" 1.864883589304147
  "NPC_HET_H3K27me3_R2" 0.754880092204074
  "NPC_KO_H3K27me3_R1" 1.04970114483556
  "NPC_KO_H3K27me3_R2" 0.7278916734142225
  "NPC_WT_H3K27me3_R1" 0.8212428590880295
  "NPC_WT_H3K27me3_R2" 0.8952950632893033
)


for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"

  bamCoverage --bam output/bowtie2_endtoend/${sample}.dupmark.sorted.bam \
      --outFileName output/bigwig_DiffBind_TMM_ONLY/${sample}.bw \
      --outFileFormat bigwig \
      --binSize 1 \
      --numberOfProcessors 7 \
      --extendReads \
      --scaleFactor $SF
done

