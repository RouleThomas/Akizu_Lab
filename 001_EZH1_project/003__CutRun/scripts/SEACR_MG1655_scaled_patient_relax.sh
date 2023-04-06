#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=75G
#SBATCH --time=100:00:00

module load R/4.2.2



bash scripts/SEACR_1.3.sh output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph non relaxed output/SEACR/8wN_iPSCpatient_H3K27me3_R1_MG1655_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph non relaxed output/SEACR/8wN_iPSCpatient_H3K27me3_R2_MG1655_scaled
