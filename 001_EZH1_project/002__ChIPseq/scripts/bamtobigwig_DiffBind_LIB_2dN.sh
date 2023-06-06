#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "2dN_HET_H3K27me3_R1" 2.600552305298599
  "2dN_HET_H3K27me3_R2" 3.081437145307326
  "2dN_KO_H3K27me3_R1" 4.261940572353047
  "2dN_KO_H3K27me3_R2" 5.046432222882744
  "2dN_WT_H3K27me3_R1" 3.927336421528873
  "2dN_WT_H3K27me3_R2" 3.5348142083978
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

