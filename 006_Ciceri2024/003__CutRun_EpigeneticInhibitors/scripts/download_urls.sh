#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00




curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/059/SRR26553059/SRR26553059_1.fastq.gz -o SRR26553059_GSM7869114_CnR_EZH2inh_noab_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/059/SRR26553059/SRR26553059_2.fastq.gz -o SRR26553059_GSM7869114_CnR_EZH2inh_noab_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/017/SRR26553017/SRR26553017_1.fastq.gz -o SRR26553017_GSM7869101_CnR_DOT1Linh_noab_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/017/SRR26553017/SRR26553017_2.fastq.gz -o SRR26553017_GSM7869101_CnR_DOT1Linh_noab_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/018/SRR26553018/SRR26553018_1.fastq.gz -o SRR26553018_GSM7869100_CnR_DOT1Linh_H3K4me3_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/018/SRR26553018/SRR26553018_2.fastq.gz -o SRR26553018_GSM7869100_CnR_DOT1Linh_H3K4me3_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/060/SRR26553060/SRR26553060_1.fastq.gz -o SRR26553060_GSM7869113_CnR_EZH2inh_noab_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/060/SRR26553060/SRR26553060_2.fastq.gz -o SRR26553060_GSM7869113_CnR_EZH2inh_noab_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/023/SRR26553023/SRR26553023_1.fastq.gz -o SRR26553023_GSM7869095_CnR_DMSO_H3K4me3_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/023/SRR26553023/SRR26553023_2.fastq.gz -o SRR26553023_GSM7869095_CnR_DMSO_H3K4me3_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/065/SRR26553065/SRR26553065_1.fastq.gz -o SRR26553065_GSM7869092_CnR_DMSO_H3K27me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/065/SRR26553065/SRR26553065_2.fastq.gz -o SRR26553065_GSM7869092_CnR_DMSO_H3K27me3_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/064/SRR26553064/SRR26553064_1.fastq.gz -o SRR26553064_GSM7869093_CnR_DMSO_H3K27me3_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/064/SRR26553064/SRR26553064_2.fastq.gz -o SRR26553064_GSM7869093_CnR_DMSO_H3K27me3_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/019/SRR26553019/SRR26553019_1.fastq.gz -o SRR26553019_GSM7869099_CnR_DOT1Linh_H3K4me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/019/SRR26553019/SRR26553019_2.fastq.gz -o SRR26553019_GSM7869099_CnR_DOT1Linh_H3K4me3_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/061/SRR26553061/SRR26553061_1.fastq.gz -o SRR26553061_GSM7869112_CnR_EZH2inh_H3K4me3_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/061/SRR26553061/SRR26553061_2.fastq.gz -o SRR26553061_GSM7869112_CnR_EZH2inh_H3K4me3_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/024/SRR26553024/SRR26553024_1.fastq.gz -o SRR26553024_GSM7869094_CnR_DMSO_H3K4me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/024/SRR26553024/SRR26553024_2.fastq.gz -o SRR26553024_GSM7869094_CnR_DMSO_H3K4me3_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/021/SRR26553021/SRR26553021_1.fastq.gz -o SRR26553021_GSM7869097_CnR_DOT1Linh_H3K27me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/021/SRR26553021/SRR26553021_2.fastq.gz -o SRR26553021_GSM7869097_CnR_DOT1Linh_H3K27me3_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/063/SRR26553063/SRR26553063_1.fastq.gz -o SRR26553063_GSM7869110_CnR_EZH2inh_H3K27me3_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/063/SRR26553063/SRR26553063_2.fastq.gz -o SRR26553063_GSM7869110_CnR_EZH2inh_H3K27me3_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/062/SRR26553062/SRR26553062_1.fastq.gz -o SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/062/SRR26553062/SRR26553062_2.fastq.gz -o SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/022/SRR26553022/SRR26553022_1.fastq.gz -o SRR26553022_GSM7869096_CnR_DMSO_noab_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/022/SRR26553022/SRR26553022_2.fastq.gz -o SRR26553022_GSM7869096_CnR_DMSO_noab_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/020/SRR26553020/SRR26553020_1.fastq.gz -o SRR26553020_GSM7869098_CnR_DOT1Linh_H3K27me3_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/020/SRR26553020/SRR26553020_2.fastq.gz -o SRR26553020_GSM7869098_CnR_DOT1Linh_H3K27me3_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/009/SRR26553009/SRR26553009_1.fastq.gz -o SRR26553009_GSM7869109_CnR_EZH2inh_H3K27me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/009/SRR26553009/SRR26553009_2.fastq.gz -o SRR26553009_GSM7869109_CnR_EZH2inh_H3K27me3_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/010/SRR26553010/SRR26553010_1.fastq.gz -o SRR26553010_GSM7869108_CnR_EHMTinh_noab_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/010/SRR26553010/SRR26553010_2.fastq.gz -o SRR26553010_GSM7869108_CnR_EHMTinh_noab_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/016/SRR26553016/SRR26553016_1.fastq.gz -o SRR26553016_GSM7869102_CnR_DOT1Linh_noab_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/016/SRR26553016/SRR26553016_2.fastq.gz -o SRR26553016_GSM7869102_CnR_DOT1Linh_noab_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/011/SRR26553011/SRR26553011_1.fastq.gz -o SRR26553011_GSM7869107_CnR_EHMTinh_noab_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/011/SRR26553011/SRR26553011_2.fastq.gz -o SRR26553011_GSM7869107_CnR_EHMTinh_noab_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/014/SRR26553014/SRR26553014_1.fastq.gz -o SRR26553014_GSM7869104_CnR_EHMTinh_H3K27me3_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/014/SRR26553014/SRR26553014_2.fastq.gz -o SRR26553014_GSM7869104_CnR_EHMTinh_H3K27me3_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/013/SRR26553013/SRR26553013_1.fastq.gz -o SRR26553013_GSM7869105_CnR_EHMTinh_H3K4me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/013/SRR26553013/SRR26553013_2.fastq.gz -o SRR26553013_GSM7869105_CnR_EHMTinh_H3K4me3_rep1_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/012/SRR26553012/SRR26553012_1.fastq.gz -o SRR26553012_GSM7869106_CnR_EHMTinh_H3K4me3_rep2_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/012/SRR26553012/SRR26553012_2.fastq.gz -o SRR26553012_GSM7869106_CnR_EHMTinh_H3K4me3_rep2_Homo_sapiens_OTHER_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/015/SRR26553015/SRR26553015_1.fastq.gz -o SRR26553015_GSM7869103_CnR_EHMTinh_H3K27me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/015/SRR26553015/SRR26553015_2.fastq.gz -o SRR26553015_GSM7869103_CnR_EHMTinh_H3K27me3_rep1_Homo_sapiens_OTHER_2.fastq.gz