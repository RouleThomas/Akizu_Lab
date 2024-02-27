#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00



curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/050/SRR17882750/SRR17882750_1.fastq.gz -o SRR17882750_GSM5860172_CnR_IgG_NPC_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/050/SRR17882750/SRR17882750_2.fastq.gz -o SRR17882750_GSM5860172_CnR_IgG_NPC_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/059/SRR17882759/SRR17882759_1.fastq.gz -o SRR17882759_GSM5860163_CnR_H3K4me3_NPC_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/059/SRR17882759/SRR17882759_2.fastq.gz -o SRR17882759_GSM5860163_CnR_H3K4me3_NPC_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/051/SRR17882751/SRR17882751_1.fastq.gz -o SRR17882751_GSM5860171_CnR_IgG_NPC_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/051/SRR17882751/SRR17882751_2.fastq.gz -o SRR17882751_GSM5860171_CnR_IgG_NPC_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/069/SRR17882769/SRR17882769_1.fastq.gz -o SRR17882769_GSM5860153_CnR_H3K27ac_D53_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/069/SRR17882769/SRR17882769_2.fastq.gz -o SRR17882769_GSM5860153_CnR_H3K27ac_D53_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/062/SRR17882762/SRR17882762_1.fastq.gz -o SRR17882762_GSM5860160_CnR_H3K27me3_NPC_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/062/SRR17882762/SRR17882762_2.fastq.gz -o SRR17882762_GSM5860160_CnR_H3K27me3_NPC_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/061/SRR17882761/SRR17882761_1.fastq.gz -o SRR17882761_GSM5860161_CnR_H3K4me3_D53_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/061/SRR17882761/SRR17882761_2.fastq.gz -o SRR17882761_GSM5860161_CnR_H3K4me3_D53_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/068/SRR17882768/SRR17882768_1.fastq.gz -o SRR17882768_GSM5860154_CnR_H3K27ac_D53_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/068/SRR17882768/SRR17882768_2.fastq.gz -o SRR17882768_GSM5860154_CnR_H3K27ac_D53_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/052/SRR17882752/SRR17882752_1.fastq.gz -o SRR17882752_GSM5860170_CnR_IgG_D53_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/052/SRR17882752/SRR17882752_2.fastq.gz -o SRR17882752_GSM5860170_CnR_IgG_D53_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/063/SRR17882763/SRR17882763_1.fastq.gz -o SRR17882763_GSM5860159_CnR_H3K27me3_NPC_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/063/SRR17882763/SRR17882763_2.fastq.gz -o SRR17882763_GSM5860159_CnR_H3K27me3_NPC_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/053/SRR17882753/SRR17882753_1.fastq.gz -o SRR17882753_GSM5860169_CnR_IgG_D53_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/053/SRR17882753/SRR17882753_2.fastq.gz -o SRR17882753_GSM5860169_CnR_IgG_D53_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/067/SRR17882767/SRR17882767_1.fastq.gz -o SRR17882767_GSM5860155_CnR_H3K27ac_NPC_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/067/SRR17882767/SRR17882767_2.fastq.gz -o SRR17882767_GSM5860155_CnR_H3K27ac_NPC_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/054/SRR17882754/SRR17882754_1.fastq.gz -o SRR17882754_GSM5860168_CnR_H3K9me3_NPC_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/054/SRR17882754/SRR17882754_2.fastq.gz -o SRR17882754_GSM5860168_CnR_H3K9me3_NPC_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/055/SRR17882755/SRR17882755_1.fastq.gz -o SRR17882755_GSM5860167_CnR_H3K9me3_NPC_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/055/SRR17882755/SRR17882755_2.fastq.gz -o SRR17882755_GSM5860167_CnR_H3K9me3_NPC_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/065/SRR17882765/SRR17882765_1.fastq.gz -o SRR17882765_GSM5860157_CnR_H3K27me3_D53_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/065/SRR17882765/SRR17882765_2.fastq.gz -o SRR17882765_GSM5860157_CnR_H3K27me3_D53_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/056/SRR17882756/SRR17882756_1.fastq.gz -o SRR17882756_GSM5860166_CnR_H3K9me3_D53_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/056/SRR17882756/SRR17882756_2.fastq.gz -o SRR17882756_GSM5860166_CnR_H3K9me3_D53_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/066/SRR17882766/SRR17882766_1.fastq.gz -o SRR17882766_GSM5860156_CnR_H3K27ac_NPC_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/066/SRR17882766/SRR17882766_2.fastq.gz -o SRR17882766_GSM5860156_CnR_H3K27ac_NPC_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/064/SRR17882764/SRR17882764_1.fastq.gz -o SRR17882764_GSM5860158_CnR_H3K27me3_D53_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/064/SRR17882764/SRR17882764_2.fastq.gz -o SRR17882764_GSM5860158_CnR_H3K27me3_D53_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/057/SRR17882757/SRR17882757_1.fastq.gz -o SRR17882757_GSM5860165_CnR_H3K9me3_D53_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/057/SRR17882757/SRR17882757_2.fastq.gz -o SRR17882757_GSM5860165_CnR_H3K9me3_D53_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/060/SRR17882760/SRR17882760_1.fastq.gz -o SRR17882760_GSM5860162_CnR_H3K4me3_D53_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/060/SRR17882760/SRR17882760_2.fastq.gz -o SRR17882760_GSM5860162_CnR_H3K4me3_D53_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/058/SRR17882758/SRR17882758_1.fastq.gz -o SRR17882758_GSM5860164_CnR_H3K4me3_NPC_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/058/SRR17882758/SRR17882758_2.fastq.gz -o SRR17882758_GSM5860164_CnR_H3K4me3_NPC_rep2_Homo_sapiens_OTHER_2.fastq.gz


