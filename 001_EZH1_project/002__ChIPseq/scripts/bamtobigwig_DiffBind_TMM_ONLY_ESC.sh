#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "ESC_HET_H3K27me3_R1" 1.24207679214292
  "ESC_HET_H3K27me3_R2" 1.603887309594518
  "ESC_KO_H3K27me3_R1" 1.280243840362811
  "ESC_KO_H3K27me3_R2" 1.833315336287782
  "ESC_WT_H3K27me3_R1" 0.7227935624246804
  "ESC_WT_H3K27me3_R2" 0.7840824980241121
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

