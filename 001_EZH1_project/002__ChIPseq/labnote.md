# Import files from Google drive to the cluster


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
--> The number of reads is similar in the 2 (paired)-read; and reduced/2 when looking at only 1 of the 2 lanes; so concatenation keep data integrity and combine all reads. Just in case, backup has been performed from all files: **To be deleted after mapping (XXX)**


## Concatenate the 2 lanes of all samples
```bash
# example for 1 file:
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz > input/ESC_KO_input_R2_1.fq.gz
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz > input/ESC_KO_input_R2_2.fq.gz

# time per time:
sbatch scripts/concat_1.sh # 11256971; fail with sample17, repeat it manually
sbatch scripts/concat_2.sh # 11256972
sbatch scripts/concat_3.sh # 11256973
```
--> `2dN_KO_input_R2_2.fq.gz` is missing... So I look for it `grep -rn '2dN_KO_input_R2' scripts/`, mistake corrected in concat_2, and re-run manually for both 2dN_KO_input_R1 and 2dN_KO_input_R2. Will double check all is good at the fastqc step



# Trimming and Fastqc
```
2dN_HET_H3K27me3_R1_1
2dN_HET_H3K27me3_R1_2
2dN_HET_H3K27me3_R2_1
2dN_HET_H3K27me3_R2_2
2dN_HET_input_R1_1
2dN_HET_input_R1_2
2dN_HET_input_R2_1
2dN_HET_input_R2_2
2dN_KO_H3K27me3_R1_1
2dN_KO_H3K27me3_R1_2
2dN_KO_H3K27me3_R2_1
2dN_KO_H3K27me3_R2_2
2dN_KO_input_R1_1
2dN_KO_input_R1_2
2dN_KO_input_R2_1
2dN_KO_input_R2_2
2dN_WT_H3K27me3_R1_1
2dN_WT_H3K27me3_R1_2
2dN_WT_H3K27me3_R2_1
2dN_WT_H3K27me3_R2_2
2dN_WT_input_R1_1
2dN_WT_input_R1_2
2dN_WT_input_R2_1
2dN_WT_input_R2_2
ESC_HET_H3K27me3_R1_1
ESC_HET_H3K27me3_R1_2
ESC_HET_H3K27me3_R2_1
ESC_HET_H3K27me3_R2_2
ESC_HET_input_R1_1
ESC_HET_input_R1_2
ESC_HET_input_R2_1
ESC_HET_input_R2_2
ESC_KO_H3K27me3_R1_1
ESC_KO_H3K27me3_R1_2
ESC_KO_H3K27me3_R2_1
ESC_KO_H3K27me3_R2_2
ESC_KO_input_R1_1
ESC_KO_input_R1_2
ESC_KO_input_R2_1
ESC_KO_input_R2_2
ESC_WT_H3K27me3_R1_1
ESC_WT_H3K27me3_R1_2
ESC_WT_H3K27me3_R2_1
ESC_WT_H3K27me3_R2_2
ESC_WT_H3K27me3_R3_1
ESC_WT_H3K27me3_R3_2
ESC_WT_input_R1_1
ESC_WT_input_R1_2
ESC_WT_input_R2_1
ESC_WT_input_R2_2
ESC_WT_input_R3_1
ESC_WT_input_R3_2
NPC_HET_H3K27me3_R1_1
NPC_HET_H3K27me3_R1_2
NPC_HET_H3K27me3_R2_1
NPC_HET_H3K27me3_R2_2
NPC_HET_input_R1_1
NPC_HET_input_R1_2
NPC_HET_input_R2_1
NPC_HET_input_R2_2
NPC_KO_H3K27me3_R1_1
NPC_KO_H3K27me3_R1_2
NPC_KO_H3K27me3_R2_1
NPC_KO_H3K27me3_R2_2
NPC_KO_input_R1_1
NPC_KO_input_R1_2
NPC_KO_input_R2_1
NPC_KO_input_R2_2
NPC_WT_H3K27me3_R1_1
NPC_WT_H3K27me3_R1_2
NPC_WT_H3K27me3_R2_1
NPC_WT_H3K27me3_R2_2
NPC_WT_input_R1_1
NPC_WT_input_R1_2
NPC_WT_input_R2_1
NPC_WT_input_R2_2
```

## Fastqc on raw reads
```bash
# example for 1 file:
fastqc -o output/fastqc input/ESC_KO_input_R2_1.fq.gz

# time per time:
sbatch 
```
--> Double check files 2dN_KO_input_R1 and 2dN_KO_input_R2 are not the same!





## Fastqc on fastp-trimmed-reads
### Trimming with Fastp
XXX

### Fastqc trimmed reads
XXX












