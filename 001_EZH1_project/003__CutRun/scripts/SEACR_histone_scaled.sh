#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=75G
#SBATCH --time=100:00:00

module load R/4.2.2



bash scripts/SEACR_1.3.sh output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_WT_H3K27me3_R1_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_WT_H3K27me3_R2_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_WT_H3K27me3_R3_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_WT_H3K27me3_R4_histone_scaled
   
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_KO_H3K27me3_R1_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_KO_H3K27me3_R2_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_KO_H3K27me3_R3_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_KO_H3K27me3_R4_histone_scaled
  
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_HET_H3K27me3_R1_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_HET_H3K27me3_R2_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_HET_H3K27me3_R3_histone_scaled
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.fragments.histonescaled.bedgraph output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.fragments.histonescaled.bedgraph non stringent output/SEACR/8wN_HET_H3K27me3_R4_histone_scaled
   
  

