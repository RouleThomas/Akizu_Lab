#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=100G

module load sam-bcf-tools/1.6

samples_and_scaling_factors=(
  "8wN_HET_H3K27me3_R1" 1.7585726939137
  "8wN_HET_H3K27me3_R2" 1.62938124762748
  "8wN_HET_H3K27me3_R3" 1
  "8wN_HET_H3K27me3_R4" 1.76768315829432
  "8wN_HET_IGG_R1" 2.0094484954298
  "8wN_HET_IGG_R2" 2.39088014788949
  "8wN_HET_IGG_R3" 1.15040566909726
  "8wN_HET_IGG_R4" 1
  "8wN_KO_H3K27me3_R1" 3.22396768402154
  "8wN_KO_H3K27me3_R2" 1.44389587073609
  "8wN_KO_H3K27me3_R3" 7.20436864153202
  "8wN_KO_H3K27me3_R4" 1
  "8wN_KO_IGG_R1" 3.07833619282755
  "8wN_KO_IGG_R2" 1.62484830731891
  "8wN_KO_IGG_R3" 2.72301962589446
  "8wN_KO_IGG_R4" 1
  "8wN_WT_H3K27me3_R1" 1
  "8wN_WT_H3K27me3_R2" 1.98686838426203
  "8wN_WT_H3K27me3_R3" 1.00838825644098
  "8wN_WT_H3K27me3_R4" 2.09002396644697
  "8wN_WT_IGG_R1" 1.76358031548545
  "8wN_WT_IGG_R2" 3.90124416796268
  "8wN_WT_IGG_R3" 1
  "8wN_WT_IGG_R4" 5.3987447233948
  "8wN_iPSCpatient_H3K27me3_R2" 1
  "8wN_iPSCpatient_IGG_R2" 1
  "8wN_iPSCpatient_H3K27me3_R1" 2.08938964333682
  "8wN_iPSCpatient_IGG_R1" 1.37815014000622
)

for ((i=0; i<${#samples_and_scaling_factors[@]}; i+=2)); do
  sample="${samples_and_scaling_factors[i]}"
  SF="${samples_and_scaling_factors[i+1]}"

  libSize=$(samtools view -c -F 260 output/bowtie2/${sample}.dupmark.sorted.bam)
  scale=$(echo "15000000/($libSize*$SF)" | bc -l)

  bamCoverage --bam output/bowtie2/${sample}.dupmark.sorted.bam \
      --outFileName output/bigwig_MG1655_lib/${sample}.dupmark.sorted.bw \
      --outFileFormat bigwig \
      --binSize 50 \
      --numberOfProcessors 7 \
      --extendReads \
      --scaleFactor $scale
done