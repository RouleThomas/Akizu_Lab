#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=250G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
    "8wN_HET_H3K27me3_R1"	468.1179867
    "8wN_HET_H3K27me3_R2"	339.3440213
    "8wN_HET_H3K27me3_R3"	591.3682403
    "8wN_HET_H3K27me3_R4"	296.4537501
    "8wN_HET_IGG_R1"	5409.246432
    "8wN_HET_IGG_R2"	4555.041763
    "8wN_HET_IGG_R3"	5595.981227
    "8wN_HET_IGG_R4"	2938.219249
    "8wN_KO_H3K27me3_R1"	476.3562482
    "8wN_KO_H3K27me3_R2"	742.989433
    "8wN_KO_H3K27me3_R3"	598.9701982
    "8wN_KO_H3K27me3_R4"	520.665335
    "8wN_KO_IGG_R1"	3867.438278
    "8wN_KO_IGG_R2"	7371.197628
    "8wN_KO_IGG_R3"	3534.970425
    "8wN_KO_IGG_R4"	6650.710256
    "8wN_WT_H3K27me3_R1"	507.0964865
    "8wN_WT_H3K27me3_R2"	241.3607612
    "8wN_WT_H3K27me3_R3"	335.3427764
    "8wN_WT_H3K27me3_R4"	209.3100346
    "8wN_WT_IGG_R1"	3202.795783
    "8wN_WT_IGG_R2"	1721.155148
    "8wN_WT_IGG_R3"	5905.055982
    "8wN_WT_IGG_R4"	1166.006567
    "8wN_iPSCpatient_H3K27me3_R1"	324.2998693
    "8wN_iPSCpatient_H3K27me3_R2"	531.3842459
    "8wN_iPSCpatient_IGG_R1"	3797.317473
    "8wN_iPSCpatient_IGG_R2"	5487.300223
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