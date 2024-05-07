#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/020/SRR14005120/SRR14005120_1.fastq.gz -o SRR14005120_GSM5183585_Methyl_caputre_seq_1_on_H1_WT_hESCs_Homo_sapiens_Bisulfite-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/020/SRR14005120/SRR14005120_2.fastq.gz -o SRR14005120_GSM5183585_Methyl_caputre_seq_1_on_H1_WT_hESCs_Homo_sapiens_Bisulfite-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/021/SRR14005121/SRR14005121_1.fastq.gz -o SRR14005121_GSM5183586_Methyl_caputre_seq_2_on_H1_WT_hESCs_Homo_sapiens_Bisulfite-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/021/SRR14005121/SRR14005121_2.fastq.gz -o SRR14005121_GSM5183586_Methyl_caputre_seq_2_on_H1_WT_hESCs_Homo_sapiens_Bisulfite-Seq_2.fastq.gz


