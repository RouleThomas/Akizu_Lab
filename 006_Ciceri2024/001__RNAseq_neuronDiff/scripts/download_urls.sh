#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00


curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/051/SRR17873851/SRR17873851_1.fastq.gz -o SRR17873851_GSM5859029_CN3_d25_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/051/SRR17873851/SRR17873851_2.fastq.gz -o SRR17873851_GSM5859029_CN3_d25_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/053/SRR17873853/SRR17873853_1.fastq.gz -o SRR17873853_GSM5859027_CN3_ESC_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/053/SRR17873853/SRR17873853_2.fastq.gz -o SRR17873853_GSM5859027_CN3_ESC_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/057/SRR17873857/SRR17873857_1.fastq.gz -o SRR17873857_GSM5859023_CN1_d25_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/057/SRR17873857/SRR17873857_2.fastq.gz -o SRR17873857_GSM5859023_CN1_d25_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/059/SRR17873859/SRR17873859_1.fastq.gz -o SRR17873859_GSM5859021_CN1_ESC_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/059/SRR17873859/SRR17873859_2.fastq.gz -o SRR17873859_GSM5859021_CN1_ESC_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/056/SRR17873856/SRR17873856_1.fastq.gz -o SRR17873856_GSM5859024_CN1_d50_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/056/SRR17873856/SRR17873856_2.fastq.gz -o SRR17873856_GSM5859024_CN1_d50_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/045/SRR17873845/SRR17873845_1.fastq.gz -o SRR17873845_GSM5859035_CN5_d25_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/045/SRR17873845/SRR17873845_2.fastq.gz -o SRR17873845_GSM5859035_CN5_d25_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/055/SRR17873855/SRR17873855_1.fastq.gz -o SRR17873855_GSM5859025_CN1_d75_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/055/SRR17873855/SRR17873855_2.fastq.gz -o SRR17873855_GSM5859025_CN1_d75_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/046/SRR17873846/SRR17873846_1.fastq.gz -o SRR17873846_GSM5859034_CN5_NPC_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/046/SRR17873846/SRR17873846_2.fastq.gz -o SRR17873846_GSM5859034_CN5_NPC_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/054/SRR17873854/SRR17873854_1.fastq.gz -o SRR17873854_GSM5859026_CN1_d100_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/054/SRR17873854/SRR17873854_2.fastq.gz -o SRR17873854_GSM5859026_CN1_d100_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/047/SRR17873847/SRR17873847_1.fastq.gz -o SRR17873847_GSM5859033_CN5_ESC_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/047/SRR17873847/SRR17873847_2.fastq.gz -o SRR17873847_GSM5859033_CN5_ESC_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/050/SRR17873850/SRR17873850_1.fastq.gz -o SRR17873850_GSM5859030_CN3_d50_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/050/SRR17873850/SRR17873850_2.fastq.gz -o SRR17873850_GSM5859030_CN3_d50_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/048/SRR17873848/SRR17873848_1.fastq.gz -o SRR17873848_GSM5859032_CN3_d100_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/048/SRR17873848/SRR17873848_2.fastq.gz -o SRR17873848_GSM5859032_CN3_d100_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/049/SRR17873849/SRR17873849_1.fastq.gz -o SRR17873849_GSM5859031_CN3_d75_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/049/SRR17873849/SRR17873849_2.fastq.gz -o SRR17873849_GSM5859031_CN3_d75_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/052/SRR17873852/SRR17873852_1.fastq.gz -o SRR17873852_GSM5859028_CN3_NPC_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/052/SRR17873852/SRR17873852_2.fastq.gz -o SRR17873852_GSM5859028_CN3_NPC_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/043/SRR17873843/SRR17873843_1.fastq.gz -o SRR17873843_GSM5859037_CN5_d75_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/043/SRR17873843/SRR17873843_2.fastq.gz -o SRR17873843_GSM5859037_CN5_d75_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/042/SRR17873842/SRR17873842_1.fastq.gz -o SRR17873842_GSM5859038_CN5_d100_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/042/SRR17873842/SRR17873842_2.fastq.gz -o SRR17873842_GSM5859038_CN5_d100_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/058/SRR17873858/SRR17873858_1.fastq.gz -o SRR17873858_GSM5859022_CN1_NPC_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/058/SRR17873858/SRR17873858_2.fastq.gz -o SRR17873858_GSM5859022_CN1_NPC_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/044/SRR17873844/SRR17873844_1.fastq.gz -o SRR17873844_GSM5859036_CN5_d50_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/044/SRR17873844/SRR17873844_2.fastq.gz -o SRR17873844_GSM5859036_CN5_d50_Homo_sapiens_RNA-Seq_2.fastq.gz