#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=250G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "8wN_HET_H3K27me3_R1"	381.9585298
  "8wN_HET_H3K27me3_R2"	281.4296848
  "8wN_HET_H3K27me3_R3"	480.4419375
  "8wN_HET_H3K27me3_R4"	244.5464754
  "8wN_HET_IGG_R1"	4293.590866
  "8wN_HET_IGG_R2"	3660.308585
  "8wN_HET_IGG_R3"	4432.289111
  "8wN_HET_IGG_R4"	2385.665466
  "8wN_KO_H3K27me3_R1"	396.4809153
  "8wN_KO_H3K27me3_R2"	619.4858247
  "8wN_KO_H3K27me3_R3"	497.4274068
  "8wN_KO_H3K27me3_R4"	432.1408515
  "8wN_KO_IGG_R1"	3117.745455
  "8wN_KO_IGG_R2"	5938.901186
  "8wN_KO_IGG_R3"	2843.05268
  "8wN_KO_IGG_R4"	5404.483333
  "8wN_WT_H3K27me3_R1"	413.840942
  "8wN_WT_H3K27me3_R2"	201.0010131
  "8wN_WT_H3K27me3_R3"	275.8128141
  "8wN_WT_H3K27me3_R4"	174.170367
  "8wN_WT_IGG_R1"	2539.238554
  "8wN_WT_IGG_R2"	1388.843912
  "8wN_WT_IGG_R3"	4787.347969
  "8wN_WT_IGG_R4"	941.211417
  "8wN_iPSCpatient_H3K27me3_R1"	256.474348
  "8wN_iPSCpatient_H3K27me3_R2"	428.7848223
  "8wN_iPSCpatient_IGG_R1"	2979.26087
  "8wN_iPSCpatient_IGG_R2"	4398.535714
)

for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"


  bamCoverage --bam output/bowtie2/${sample}.dupmark.sorted.bam \
      --outFileName output/bigwig_Epicypher/${sample}.dupmark.sorted.bw \
      --outFileFormat bigwig \
      --binSize 50 \
      --numberOfProcessors 7 \
      --extendReads \
      --scaleFactor ${SF}
done