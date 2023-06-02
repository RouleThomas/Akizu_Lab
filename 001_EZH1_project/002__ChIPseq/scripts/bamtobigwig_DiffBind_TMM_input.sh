#!/bin/bash
#SBATCH --mem=250G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "2dN_HET_input_R1" 1.515964164
  "2dN_HET_input_R2" 1.591067619
  "2dN_KO_input_R1" 1.831466946
  "2dN_KO_input_R2" 1.793907639
  "2dN_WT_input_R1" 1.548264434
  "2dN_WT_input_R2" 1.691816077
  "NPC_HET_input_R1" 3.242710387
  "NPC_HET_input_R2" 1.310448705
  "NPC_KO_input_R1" 1.829832193
  "NPC_KO_input_R2" 1.2635567
  "NPC_WT_input_R1" 1.425558409
  "NPC_WT_input_R2" 1.556845489
  "ESC_HET_input_R1" 2.158245122
  "ESC_HET_input_R2" 2.796412762
  "ESC_KO_input_R1" 2.237804803
  "ESC_KO_input_R2" 3.1922803
  "ESC_WT_input_R1" 1.25660203
  "ESC_WT_input_R2" 1.363597191
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

