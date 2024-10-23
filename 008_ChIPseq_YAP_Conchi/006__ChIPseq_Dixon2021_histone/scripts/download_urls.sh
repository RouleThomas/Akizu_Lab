#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00


curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/024/SRR14005224/SRR14005224_1.fastq.gz -o SRR14005224_GSM5183570_EZH2_ChIP-seq_on_H1_QSER1_KO_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/024/SRR14005224/SRR14005224_2.fastq.gz -o SRR14005224_GSM5183570_EZH2_ChIP-seq_on_H1_QSER1_KO_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/002/SRR14005202/SRR14005202_1.fastq.gz -o SRR14005202_GSM5183548_EZH2_ChIP-seq_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/002/SRR14005202/SRR14005202_2.fastq.gz -o SRR14005202_GSM5183548_EZH2_ChIP-seq_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz