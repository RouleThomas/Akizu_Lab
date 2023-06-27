#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00



module load SAMtools/1.16.1-GCC-11.3.0

samples_and_scaling_factors=(
  "ESC_WT_H3K27me3_R1" 1.010075401
  "ESC_WT_H3K27me3_R2" 0.905242149
  "NPC_HET_H3K27me3_R1" 2.22126017
  "NPC_HET_H3K27me3_R2" 0.986039553
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





