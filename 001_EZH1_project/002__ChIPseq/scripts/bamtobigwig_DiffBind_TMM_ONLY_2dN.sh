#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "2dN_HET_H3K27me3_R1" 0.8731719490367953
  "2dN_HET_H3K27me3_R2" 0.9161719185611117
  "2dN_KO_H3K27me3_R1" 1.055764755341827
  "2dN_KO_H3K27me3_R2" 1.026240244688514
  "2dN_WT_H3K27me3_R1" 0.887629323154237
  "2dN_WT_H3K27me3_R2" 0.9702459524977284
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

