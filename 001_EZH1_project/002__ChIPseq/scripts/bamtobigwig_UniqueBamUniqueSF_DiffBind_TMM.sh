#!/bin/bash
#SBATCH --mem=400G

module load SAMtools/1.16.1-GCC-11.3.0

samples_and_scaling_factors=(
      "2dN_HET_H3K27me3_R1" 1.328505368
  "2dN_HET_H3K27me3_R2" 1.386290557
  "2dN_KO_H3K27me3_R1" 1.431311572
  "2dN_KO_H3K27me3_R2" 1.38381988
  "2dN_WT_H3K27me3_R1" 1.281454092
  "2dN_WT_H3K27me3_R2" 1.073748366
    "ESC_HET_H3K27me3_R1" 2.299537264
  "ESC_HET_H3K27me3_R2" 3.45113078
  "ESC_KO_H3K27me3_R1" 3.703952006
  "ESC_KO_H3K27me3_R2" 4.454690232
  "ESC_WT_H3K27me3_R1" 1.010075401
  "ESC_WT_H3K27me3_R2" 0.905242149
  "NPC_HET_H3K27me3_R1" 2.22126017
  "NPC_HET_H3K27me3_R2" 0.986039553
  "NPC_KO_H3K27me3_R1" 1.46400002
  "NPC_KO_H3K27me3_R2" 1.106896051
  "NPC_WT_H3K27me3_R1" 1.029216259
  "NPC_WT_H3K27me3_R2" 1.047211652
)


for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"

  bamCoverage --bam output/bowtie2_endtoend/${sample}.dupmark.sorted.bam \
      --outFileName output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/${sample}.bw \
      --outFileFormat bigwig \
      --binSize 1 \
      --numberOfProcessors 7 \
      --extendReads \
      --scaleFactor $SF
done

