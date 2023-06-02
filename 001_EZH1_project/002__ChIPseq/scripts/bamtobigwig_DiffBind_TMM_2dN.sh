#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "2dN_HET_H3K27me3_R1" 1.515964164
  "2dN_HET_H3K27me3_R2" 1.591067619
  "2dN_KO_H3K27me3_R1" 1.831466946
  "2dN_KO_H3K27me3_R2" 1.793907639
  "2dN_WT_H3K27me3_R1" 1.548264434
  "2dN_WT_H3K27me3_R2" 1.691816077
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

