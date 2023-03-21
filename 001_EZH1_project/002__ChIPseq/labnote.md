# Import files from Google drive to the cluster
#### 20230309

I cannot use a bash script to import files as if disconnection the transfer will fail. So cp the old-fashion way, let computer running o/n.
```bash
cp -r /home/roulet/tsclient/roule/Google\ Drive\ Streaming/Shared\ drives/akizulab/Primary\ Data/ChIPseqs/Fastq\ files/ESC\ NPC\ and\ Neurons/ Dec/ 2020/LopezN_FASTQ-216179970/* \
/scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq/input
``` 

# Concatenate fastq files
Here, each sample sequenced from 2 Illumina lanes, need concatenate them into one.

**Move all fastq files within the input folder** (now [file] are in input/folder1/folder2/[file]) using a loop:
```bash
for file in input/*/*/;
    do mv "$file"* input/;
done
```
- `input/*/*/` this go 2 folder downstream; whatever their names
- we move file to `input/` directory 

## Light test of concatenation with 1 sample to make sure we keep file integrity
```bash
# copy files as backup
## Read1 (lane1 lane2)
cp input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz input/backup/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz
cp input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz input/backup/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz
## Read2 (lane1 lane2)
cp input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz input/backup/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz
cp input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz input/backup/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz

# concatenate lane 
## Read1 (lane1 lane2)
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz > input/ESC_KO_input_R2_1.fq.gz
## Read2
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz > input/ESC_KO_input_R2_2.fq.gz
# fastqc
## Read1
fastqc -o output/fastqc input/ESC_KO_input_R2_1.fq.gz
fastqc -o output/fastqc input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz
## Read2
fastqc -o output/fastqc input/ESC_KO_input_R2_2.fq.gz
fastqc -o output/fastqc input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz
```
--> The number of reads is similar in the 2 (paired)-read; and reduced in one of the lane; so concatenation keep data integrity and combine all reads

--> In case, backup has been performed from all files: **To be deleted after mapping (XXX)**


## Concatenate the 2 lanes of all samples


```bash
cat 10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq 10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq > ESC_KO_input_R2_1
cat 10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq 10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq > ESC_KO_input_R2_2

cat 10-H3K27me3-Neurons--DID2--Het2A-2_S20_L001_R1_001.fastq 10-H3K27me3-Neurons--DID2--Het2A-2_S20_L002_R1_001.fastq > 2dN_HET_H3K27me3_R2_1
cat 10-H3K27me3-Neurons--DID2--Het2A-2_S20_L001_R2_001.fastq 10-H3K27me3-Neurons--DID2--Het2A-2_S20_L002_R2_001.fastq > 2dN_HET_H3K27me3_R2_2

cat 11-ESCs-Pan-Het-2A-1-INPUT_S9_L001_R1_001.fastq 11-ESCs-Pan-Het-2A-1-INPUT_S9_L002_R1_001.fastq > ESC_HET_input_R1_1
cat 11-ESCs-Pan-Het-2A-1-INPUT_S9_L001_R2_001.fastq 11-ESCs-Pan-Het-2A-1-INPUT_S9_L002_R2_001.fastq > ESC_HET_input_R1_2

cat 11-H3K27me3-Neurons--DID2--EZH1-KO-1_S21_L001_R1_001.fastq 11-H3K27me3-Neurons--DID2--EZH1-KO-1_S21_L002_R1_001.fastq > 2dN_KO_H3K27me3_R1_1
cat 11-H3K27me3-Neurons--DID2--EZH1-KO-1_S21_L001_R2_001.fastq 11-H3K27me3-Neurons--DID2--EZH1-KO-1_S21_L002_R2_001.fastq > 2dN_KO_H3K27me3_R1_2

cat 12-ESCs-Pan-Het-2A-2-INPUT_S10_L001_R1_001.fastq 12-ESCs-Pan-Het-2A-2-INPUT_S10_L002_R1_001.fastq > ESC_HET_input_R2_1
cat 12-ESCs-Pan-Het-2A-2-INPUT_S10_L001_R2_001.fastq 12-ESCs-Pan-Het-2A-2-INPUT_S10_L002_R2_001.fastq > ESC_HET_input_R2_2

cat 12-H3K27me3-Neurons--DID2--EZH1-KO-2_S22_L001_R1_001.fastq 12-H3K27me3-Neurons--DID2--EZH1-KO-2_S22_L002_R1_001.fastq > 2dN_KO_H3K27me3_R2_1
cat 12-H3K27me3-Neurons--DID2--EZH1-KO-2_S22_L001_R2_001.fastq 12-H3K27me3-Neurons--DID2--EZH1-KO-2_S22_L002_R2_001.fastq > 2dN_KO_H3K27me3_R2_2

cat 13-INPUT-NPCs-WTS2-1_S23_L001_R1_001.fastq 13-INPUT-NPCs-WTS2-1_S23_L002_R1_001.fastq > NPC_WT_input_R1_1
cat 13-INPUT-NPCs-WTS2-1_S23_L001_R2_001.fastq 13-INPUT-NPCs-WTS2-1_S23_L002_R2_001.fastq > NPC_WT_input_R1_2

cat 14-INPUT-NPCs-WTS2-2_S24_L001_R1_001.fastq 14-INPUT-NPCs-WTS2-2_S24_L002_R1_001.fastq > NPC_WT_input_R2_1
cat 14-INPUT-NPCs-WTS2-2_S24_L001_R2_001.fastq 14-INPUT-NPCs-WTS2-2_S24_L002_R2_001.fastq > NPC_WT_input_R2_2

cat 15-INPUT-NPCs-Het2A-1_S25_L001_R1_001.fastq 15-INPUT-NPCs-Het2A-1_S25_L002_R1_001.fastq > NPC_HET_input_R1_1
cat 15-INPUT-NPCs-Het2A-1_S25_L001_R2_001.fastq 15-INPUT-NPCs-Het2A-1_S25_L002_R2_001.fastq > NPC_HET_input_R1_2

cat 16-INPUT-NPCs-Het2A-2_S26_L001_R1_001.fastq 16-INPUT-NPCs-Het2A-2_S26_L002_R1_001.fastq > NPC_HET_input_R2_1
cat 16-INPUT-NPCs-Het2A-2_S26_L001_R2_001.fastq 16-INPUT-NPCs-Het2A-2_S26_L002_R2_001.fastq > NPC_HET_input_R2_2

cat 17-INPUT-NPCs-EZH1-KO-1_S27_L001_R1_001.fastq 17-INPUT-NPCs-EZH1-KO-1_S27_L002_R1_001.fastq > NPC_KO_input_R1_1
cat 17-INPUT-NPCs-EZH1-KO-1_S27_L001_R2_001.fastq 17-INPUT-NPCs-EZH1-KO-1_S27_L002_R2_001.fastq > NPC_KO_input_R1_2

cat 18-INPUT-NPCs-EZH1-KO-2_S28_L001_R1_001.fastq 18-INPUT-NPCs-EZH1-KO-2_S28_L002_R1_001.fastq > NPC_KO_input_R2_1
cat 18-INPUT-NPCs-EZH1-KO-2_S28_L001_R2_001.fastq 18-INPUT-NPCs-EZH1-KO-2_S28_L002_R2_001.fastq > NPC_KO_input_R2_2

cat 19-INPUT-Neurons--DID2--WTS2-1_S29_L001_R1_001.fastq 19-INPUT-Neurons--DID2--WTS2-1_S29_L002_R1_001.fastq > 2dN_WT_input_R1_1
cat 19-INPUT-Neurons--DID2--WTS2-1_S29_L001_R2_001.fastq 19-INPUT-Neurons--DID2--WTS2-1_S29_L002_R2_001.fastq > 2dN_WT_input_R1_2

cat 1-ESCs-WT-S2-1-H3k27me3_S1_L001_R1_001.fastq 1-ESCs-WT-S2-1-H3k27me3_S1_L002_R1_001.fastq > ESC_WT_H3K27me3_R3_1
cat 1-ESCs-WT-S2-1-H3k27me3_S1_L001_R2_001.fastq 1-ESCs-WT-S2-1-H3k27me3_S1_L002_R2_001.fastq > ESC_WT_H3K27me3_R3_2

cat 1-H3K27me3-ESCs-H9-WT-1_S35_L001_R1_001.fastq 1-H3K27me3-ESCs-H9-WT-1_S35_L002_R1_001.fastq > ESC_WT_H3K27me3_R1_1
cat 1-H3K27me3-ESCs-H9-WT-1_S35_L001_R2_001.fastq 1-H3K27me3-ESCs-H9-WT-1_S35_L002_R2_001.fastq > ESC_WT_H3K27me3_R1_2

cat 1-H3K27me3-NPCs-WTS2-1_S11_L001_R1_001.fastq 1-H3K27me3-NPCs-WTS2-1_S11_L002_R1_001.fastq > NPC_WT_H3K27me3_R1_1
cat 1-H3K27me3-NPCs-WTS2-1_S11_L001_R2_001.fastq 1-H3K27me3-NPCs-WTS2-1_S11_L002_R2_001.fastq > NPC_WT_H3K27me3_R1_2

cat 20-INPUT-Neurons--DID2--WTS2-2_S30_L001_R1_001.fastq 20-INPUT-Neurons--DID2--WTS2-2_S30_L002_R1_001.fastq > 2dN_WT_input_R2_1
cat 20-INPUT-Neurons--DID2--WTS2-2_S30_L001_R2_001.fastq 20-INPUT-Neurons--DID2--WTS2-2_S30_L002_R2_001.fastq > 2dN_WT_input_R2_2

cat 21-INPUT-Neurons--DID2--Het2A-1_S31_L001_R1_001.fastq 21-INPUT-Neurons--DID2--Het2A-1_S31_L002_R1_001.fastq > 2dN_HET_input_R1_1
cat 21-INPUT-Neurons--DID2--Het2A-1_S31_L001_R2_001.fastq 21-INPUT-Neurons--DID2--Het2A-1_S31_L002_R2_001.fastq > 2dN_HET_input_R1_2

cat 22-INPUT-Neurons--DID2--Het2A-2_S32_L001_R1_001.fastq 22-INPUT-Neurons--DID2--Het2A-2_S32_L002_R1_001.fastq > 2dN_HET_input_R2_1
cat 22-INPUT-Neurons--DID2--Het2A-2_S32_L001_R2_001.fastq 22-INPUT-Neurons--DID2--Het2A-2_S32_L002_R2_001.fastq > 2dN_HET_input_R2_2

cat 23-INPUT-Neurons--DID2--EZH1-KO-1_S33_L001_R1_001.fastq 23-INPUT-Neurons--DID2--EZH1-KO-1_S33_L002_R1_001.fastq > 2dN_KO_input_R1_1
cat 23-INPUT-Neurons--DID2--EZH1-KO-1_S33_L001_R2_001.fastq 23-INPUT-Neurons--DID2--EZH1-KO-1_S33_L002_R2_001.fastq > 2dN_KO_input_R1_2

cat 24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L001_R1_001.fastq 24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L002_R1_001.fastq > 2dN_KO_input_R2_1
cat 24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L001_R2_001.fastq 24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L002_R2_001.fastq > 2dN_KO_input_R2_1

cat 2-H3K27me3-ESCs-H9-WT-2_S36_L001_R1_001.fastq 2-H3K27me3-ESCs-H9-WT-2_S36_L002_R1_001.fastq > ESC_WT_H3K27me3_R2_1
cat 2-H3K27me3-ESCs-H9-WT-2_S36_L001_R2_001.fastq 2-H3K27me3-ESCs-H9-WT-2_S36_L002_R2_001.fastq > ESC_WT_H3K27me3_R2_2

cat 2-H3K27me3-NPCs-WTS2-2_S12_L001_R1_001.fastq 2-H3K27me3-NPCs-WTS2-2_S12_L002_R1_001.fastq > NPC_WT_H3K27me3_R2_1
cat 2-H3K27me3-NPCs-WTS2-2_S12_L001_R2_001.fastq 2-H3K27me3-NPCs-WTS2-2_S12_L002_R2_001.fastq > NPC_WT_H3K27me3_R2_2

cat 3-ESCs-EZH1-KO-1-H3k27me3_S2_L001_R1_001.fastq 3-ESCs-EZH1-KO-1-H3k27me3_S2_L002_R1_001.fastq > ESC_KO_H3K27me3_R1_1
cat 3-ESCs-EZH1-KO-1-H3k27me3_S2_L001_R2_001.fastq 3-ESCs-EZH1-KO-1-H3k27me3_S2_L002_R2_001.fastq > ESC_KO_H3K27me3_R1_2

cat 3-H3K27me3-NPCs-Het2A-1_S13_L001_R1_001.fastq 3-H3K27me3-NPCs-Het2A-1_S13_L002_R1_001.fastq > NPC_HET_H3K27me3_R1_1
cat 3-H3K27me3-NPCs-Het2A-1_S13_L001_R2_001.fastq 3-H3K27me3-NPCs-Het2A-1_S13_L002_R2_001.fastq > NPC_HET_H3K27me3_R1_2

cat 3-INPUT-ESCs-H9-WT-1_S37_L001_R1_001.fastq 3-INPUT-ESCs-H9-WT-1_S37_L002_R1_001.fastq > ESC_WT_input_R1_1
cat 3-INPUT-ESCs-H9-WT-1_S37_L001_R2_001.fastq 3-INPUT-ESCs-H9-WT-1_S37_L002_R2_001.fastq > ESC_WT_input_R1_2

cat 4-ESCs-EZH1-KO-2-H3k27me3_S3_L001_R1_001.fastq 4-ESCs-EZH1-KO-2-H3k27me3_S3_L002_R1_001.fastq > ESC_KO_H3K27me3_R2_1
cat 4-ESCs-EZH1-KO-2-H3k27me3_S3_L001_R1_001.fastq 4-ESCs-EZH1-KO-2-H3k27me3_S3_L002_R1_001.fastq > ESC_KO_H3K27me3_R2_2

cat 4-H3K27me3-NPCs-Het2A-2_S14_L001_R1_001.fastq 4-H3K27me3-NPCs-Het2A-2_S14_L002_R1_001.fastq > NPC_HET_H3K27me3_R2_1
cat 4-H3K27me3-NPCs-Het2A-2_S14_L001_R2_001.fastq 4-H3K27me3-NPCs-Het2A-2_S14_L002_R2_001.fastq > NPC_HET_H3K27me3_R2_2

cat 4-INPUT-ESCs-H9-WT-2_S38_L001_R1_001.fastq 4-INPUT-ESCs-H9-WT-2_S38_L002_R1_001.fastq > ESC_WT_input_R2_1
cat 4-INPUT-ESCs-H9-WT-2_S38_L001_R2_001.fastq 4-INPUT-ESCs-H9-WT-2_S38_L002_R2_001.fastq > ESC_WT_input_R2_2

cat 5-ESCs-Pan-Het-2A-1-H3k27me3_S4_L001_R1_001.fastq 5-ESCs-Pan-Het-2A-1-H3k27me3_S4_L002_R1_001.fastq > ESC_HET_H3K27me3_R1_1
cat 5-ESCs-Pan-Het-2A-1-H3k27me3_S4_L001_R2_001.fastq 5-ESCs-Pan-Het-2A-1-H3k27me3_S4_L002_R2_001.fastq > ESC_HET_H3K27me3_R1_2

cat 5-H3K27me3-NPCs-EZH1-KO-1_S15_L001_R1_001.fastq 5-H3K27me3-NPCs-EZH1-KO-1_S15_L002_R1_001.fastq > NPC_KO_H3K27me3_R1_1
cat 5-H3K27me3-NPCs-EZH1-KO-1_S15_L001_R2_001.fastq 5-H3K27me3-NPCs-EZH1-KO-1_S15_L002_R2_001.fastq > NPC_KO_H3K27me3_R1_2

cat 6-ESCs-Pan-Het-2A-2-H3k27me3_S5_L001_R1_001.fastq 6-ESCs-Pan-Het-2A-2-H3k27me3_S5_L002_R1_001.fastq > ESC_HET_H3K27me3_R2_1
cat 6-ESCs-Pan-Het-2A-2-H3k27me3_S5_L001_R2_001.fastq 6-ESCs-Pan-Het-2A-2-H3k27me3_S5_L002_R2_001.fastq > ESC_HET_H3K27me3_R2_2

cat 6-H3K27me3-NPCs-EZH1-KO-2_S16_L001_R1_001.fastq 6-H3K27me3-NPCs-EZH1-KO-2_S16_L002_R1_001.fastq > NPC_KO_H3K27me3_R2_1
cat 6-H3K27me3-NPCs-EZH1-KO-2_S16_L001_R2_001.fastq 6-H3K27me3-NPCs-EZH1-KO-2_S16_L002_R2_001.fastq > NPC_KO_H3K27me3_R2_2

cat 7-ESCs-WT-S2-1-INPUT_S6_L001_R1_001.fastq 7-ESCs-WT-S2-1-INPUT_S6_L002_R1_001.fastq > ESC_WT_input_R3_1
cat 7-ESCs-WT-S2-1-INPUT_S6_L001_R2_001.fastq 7-ESCs-WT-S2-1-INPUT_S6_L002_R2_001.fastq > ESC_WT_input_R3_2

cat 7-H3K27me3-Neurons--DID2--WTS2-1_S17_L001_R1_001.fastq 7-H3K27me3-Neurons--DID2--WTS2-1_S17_L002_R1_001.fastq > 2dN_WT_H3K27me3_R1_1
cat 7-H3K27me3-Neurons--DID2--WTS2-1_S17_L001_R2_001.fastq 7-H3K27me3-Neurons--DID2--WTS2-1_S17_L002_R2_001.fastq > 2dN_WT_H3K27me3_R1_2

cat 8-H3K27me3-Neurons--DID2--WTS2-2_S18_L001_R1_001.fastq 8-H3K27me3-Neurons--DID2--WTS2-2_S18_L002_R1_001.fastq > 2dN_WT_H3K27me3_R2_1
cat 8-H3K27me3-Neurons--DID2--WTS2-2_S18_L001_R2_001.fastq 8-H3K27me3-Neurons--DID2--WTS2-2_S18_L002_R2_001.fastq > 2dN_WT_H3K27me3_R2_2

cat 9-ESCs-EZH1-KO-1-INPUT_S7_L001_R1_001.fastq 9-ESCs-EZH1-KO-1-INPUT_S7_L002_R1_001.fastq > ESC_KO_input_R1_1
cat 9-ESCs-EZH1-KO-1-INPUT_S7_L001_R2_001.fastq 9-ESCs-EZH1-KO-1-INPUT_S7_L002_R2_001.fastq > ESC_KO_input_R1_2

cat 9-H3K27me3-Neurons--DID2--Het2A-1_S19_L001_R1_001.fastq 9-H3K27me3-Neurons--DID2--Het2A-1_S19_L002_R1_001.fastq > 2dN_HET_H3K27me3_R1_1
cat 9-H3K27me3-Neurons--DID2--Het2A-1_S19_L001_R2_001.fastq 9-H3K27me3-Neurons--DID2--Het2A-1_S19_L002_R2_001.fastq > 2dN_HET_H3K27me3_R1_2
```

















