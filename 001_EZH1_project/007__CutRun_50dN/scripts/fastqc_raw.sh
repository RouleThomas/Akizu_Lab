#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"50dN_KOEF1aEZH1_EZH1cs_R1"
"50dN_KOEF1aEZH1_EZH1cs_R2"
"50dN_KOEF1aEZH1_EZH2_R1"
"50dN_KOEF1aEZH1_EZH2_R2"
"50dN_KOEF1aEZH1_H3K27me3_R1"
"50dN_KOEF1aEZH1_H3K27me3_R2"
"50dN_KOEF1aEZH1_IGG_R1"
"50dN_KOEF1aEZH1_SUZ12_R1"
"50dN_KOEF1aEZH1_SUZ12_R2"
"50dN_KO_EZH1cs_R1"
"50dN_KO_EZH1cs_R2"
"50dN_KO_EZH2_R1"
"50dN_KO_EZH2_R2"
"50dN_KO_H3K27me3_R1"
"50dN_KO_H3K27me3_R2"
"50dN_KO_IGG_R1"
"50dN_KO_IGG_R2"
"50dN_KO_SUZ12_R1"
"50dN_KO_SUZ12_R2"
"50dN_WTQ731E_EZH1cs_R1"
"50dN_WTQ731E_EZH1cs_R2"
"50dN_WTQ731E_EZH2_R1"
"50dN_WTQ731E_EZH2_R2"
"50dN_WTQ731E_H3K27me3_R3"
"50dN_WTQ731E_H3K27me3_R1"
"50dN_WTQ731E_H3K27me3_R2"
"50dN_WTQ731E_IGG_R1"
"50dN_WTQ731E_IGG_R2"
"50dN_WTQ731E_SUZ12_R1"
"50dN_WTQ731E_SUZ12_R2"
"PSC_WT_EZH1cs_01FA"
"PSC_WT_EZH1cs_1FA"
"PSC_WT_H3K27me1_01FA"
"PSC_WT_H3K27me1_1FA"
)


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/raw input/${x}_1.fq.gz
    fastqc -o output/fastqc/raw input/${x}_2.fq.gz
done