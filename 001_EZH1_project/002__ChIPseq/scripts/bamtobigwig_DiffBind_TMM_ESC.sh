#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "ESC_HET_H3K27me3_R1" 2.158245122
  "ESC_HET_H3K27me3_R2" 2.796412762
  "ESC_KO_H3K27me3_R1" 2.237804803
  "ESC_KO_H3K27me3_R2" 3.1922803
  "ESC_WT_H3K27me3_R1" 1.25660203
  "ESC_WT_H3K27me3_R2" 1.363597191
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

