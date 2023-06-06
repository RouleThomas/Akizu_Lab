#!/bin/bash
#SBATCH --mem=250G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "2dN_HET_input_R1" 2.600552305298599
  "2dN_HET_input_R2" 3.081437145307326
  "2dN_KO_input_R1" 4.261940572353047
  "2dN_KO_input_R2" 5.046432222882744
  "2dN_WT_input_R1" 3.927336421528873
  "2dN_WT_input_R2" 3.5348142083978
  "NPC_HET_input_R1" 7.950440136365949
  "NPC_HET_input_R2" 3.070502420170008
  "NPC_KO_input_R1" 3.450437739783858
  "NPC_KO_input_R2" 1.625662843703735
  "NPC_WT_input_R1" 3.223457486947415
  "NPC_WT_input_R2" 3.334499296587373
  "ESC_HET_input_R1" 0.4021031765386951
  "ESC_HET_input_R2" 0.2335863377032899
  "ESC_KO_input_R1" 0.4196930910333201
  "ESC_KO_input_R2" 0.3992484068890632
  "ESC_WT_input_R1" 0.5703153484578302
  "ESC_WT_input_R2" 1.057889624451299
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

