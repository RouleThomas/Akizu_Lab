#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8


x=(
    "PSC_WT_R1"
    "PSC_WT_R2"
    "PSC_WT_R3"
    "PSC_KO_R1"
    "PSC_KO_R2"
    "PSC_KO_R3"
    "PSC_KOEF1aEZH1_R1"
    "PSC_KOEF1aEZH1_R2"
    "PSC_KOEF1aEZH1_R3"
    )
        
for x in "${x[@]}"; do
echo "Processing sample ${x}"
salmon quant -i ../../Master/meta/salmon/Homo_sapiens -l A \
         -1 output/fastp/${x}_1.fq.gz \
         -2 output/fastp/${x}_2.fq.gz \
         -p 8 --validateMappings -o salmon_gtf_allchr/${x}_quant -g ../../Master/meta/gencode.v47.chr_patch_hapl_scaff.annotation.gtf
done







