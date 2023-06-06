#!/bin/bash
#SBATCH --mem=500G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "ESC_HET_H3K27me3_R1" 0.4021031765386951
  "ESC_HET_H3K27me3_R2" 0.2335863377032899
  "ESC_KO_H3K27me3_R1" 0.4196930910333201
  "ESC_KO_H3K27me3_R2" 0.3992484068890632
  "ESC_WT_H3K27me3_R1" 0.5703153484578302
  "ESC_WT_H3K27me3_R2" 1.057889624451299
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

