#!/bin/bash
#SBATCH --mem=250G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "2dN_HET_input_R1" 0.8731719490367953
  "2dN_HET_input_R2" 0.9161719185611117
  "2dN_KO_input_R1" 1.055764755341827
  "2dN_KO_input_R2" 1.026240244688514
  "2dN_WT_input_R1" 0.887629323154237
  "2dN_WT_input_R2" 0.9702459524977284
  "NPC_HET_input_R1" 1.864883589304147
  "NPC_HET_input_R2" 0.754880092204074
  "NPC_KO_input_R1" 1.04970114483556
  "NPC_KO_input_R2" 0.7278916734142225
  "NPC_WT_input_R1" 0.8212428590880295
  "NPC_WT_input_R2" 0.8952950632893033
  "ESC_HET_input_R1" 1.24207679214292
  "ESC_HET_input_R2" 1.603887309594518
  "ESC_KO_input_R1" 1.280243840362811
  "ESC_KO_input_R2" 1.833315336287782
  "ESC_WT_input_R1" 0.7227935624246804
  "ESC_WT_input_R2" 0.7840824980241121
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

