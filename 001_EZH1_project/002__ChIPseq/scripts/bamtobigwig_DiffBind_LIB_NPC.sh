#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "NPC_HET_H3K27me3_R1" 7.950440136365949
  "NPC_HET_H3K27me3_R2" 3.070502420170008
  "NPC_KO_H3K27me3_R1" 3.450437739783858
  "NPC_KO_H3K27me3_R2" 1.625662843703735
  "NPC_WT_H3K27me3_R1" 3.223457486947415
  "NPC_WT_H3K27me3_R2" 3.334499296587373
)


for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"

  bamCoverage --bam output/bowtie2_endtoend/${sample}.dupmark.sorted.bam \
      --outFileName output/bigwig_DiffBind_LIB/${sample}.bw \
      --outFileFormat bigwig \
      --binSize 1 \
      --numberOfProcessors 7 \
      --extendReads \
      --scaleFactor $SF
done

