#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00


curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/096/SRR14005196/SRR14005196_1.fastq.gz -o SRR14005196_GSM5183542_H3K27me3_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/096/SRR14005196/SRR14005196_2.fastq.gz -o SRR14005196_GSM5183542_H3K27me3_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/089/SRR14005189/SRR14005189_1.fastq.gz -o SRR14005189_GSM5183535_DNMT3B_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/089/SRR14005189/SRR14005189_2.fastq.gz -o SRR14005189_GSM5183535_DNMT3B_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/031/SRR14005231/SRR14005231_1.fastq.gz -o SRR14005231_GSM5183577_Input_sample_from_H1_WT_hESCs_1_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/031/SRR14005231/SRR14005231_2.fastq.gz -o SRR14005231_GSM5183577_Input_sample_from_H1_WT_hESCs_1_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/097/SRR14005197/SRR14005197_1.fastq.gz -o SRR14005197_GSM5183543_H3K27me3_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/097/SRR14005197/SRR14005197_2.fastq.gz -o SRR14005197_GSM5183543_H3K27me3_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/032/SRR14005232/SRR14005232_1.fastq.gz -o SRR14005232_GSM5183578_Input_sample_from_H1_WT_hESCs_2_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/032/SRR14005232/SRR14005232_2.fastq.gz -o SRR14005232_GSM5183578_Input_sample_from_H1_WT_hESCs_2_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/094/SRR14005194/SRR14005194_1.fastq.gz -o SRR14005194_GSM5183540_H3K4me3_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/094/SRR14005194/SRR14005194_2.fastq.gz -o SRR14005194_GSM5183540_H3K4me3_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/088/SRR14005188/SRR14005188_1.fastq.gz -o SRR14005188_GSM5183534_DNMT3B_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/088/SRR14005188/SRR14005188_2.fastq.gz -o SRR14005188_GSM5183534_DNMT3B_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/095/SRR14005195/SRR14005195_1.fastq.gz -o SRR14005195_GSM5183541_H3K4me3_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/095/SRR14005195/SRR14005195_2.fastq.gz -o SRR14005195_GSM5183541_H3K4me3_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/087/SRR14005187/SRR14005187_1.fastq.gz -o SRR14005187_GSM5183533_DNMT3A_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/087/SRR14005187/SRR14005187_2.fastq.gz -o SRR14005187_GSM5183533_DNMT3A_ChIP-seq_2_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/086/SRR14005186/SRR14005186_1.fastq.gz -o SRR14005186_GSM5183532_DNMT3A_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR140/086/SRR14005186/SRR14005186_2.fastq.gz -o SRR14005186_GSM5183532_DNMT3A_ChIP-seq_1_on_H1_WT_hESCs_Homo_sapiens_ChIP-Seq_2.fastq.gz


