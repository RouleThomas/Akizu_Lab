#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00


curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/083/SRR14005183/SRR14005183_1.fastq.gz -o SRR14005183_GSM5183529_QSER1-FLAG_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/083/SRR14005183/SRR14005183_2.fastq.gz -o SRR14005183_GSM5183529_QSER1-FLAG_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/037/SRR14005237/SRR14005237_1.fastq.gz -o SRR14005237_GSM5183583_Input_sample_from_QSER1-FLAG_TET1-V5_QSER1_KO_H1_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/037/SRR14005237/SRR14005237_2.fastq.gz -o SRR14005237_GSM5183583_Input_sample_from_QSER1-FLAG_TET1-V5_QSER1_KO_H1_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/082/SRR14005182/SRR14005182_1.fastq.gz -o SRR14005182_GSM5183528_QSER1-FLAG_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/082/SRR14005182/SRR14005182_2.fastq.gz -o SRR14005182_GSM5183528_QSER1-FLAG_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/030/SRR14005230/SRR14005230_1.fastq.gz -o SRR14005230_GSM5183576_Input_sample_from_QSER1-FLAG_TET1-V5_WT_H1_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/030/SRR14005230/SRR14005230_2.fastq.gz -o SRR14005230_GSM5183576_Input_sample_from_QSER1-FLAG_TET1-V5_WT_H1_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz