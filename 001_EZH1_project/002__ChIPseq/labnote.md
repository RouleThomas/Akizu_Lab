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
sbatch scripts/concat_2.sh # 11256972 ok
sbatch scripts/concat_3.sh # 11256973 ok
```
--> `2dN_KO_input_R2_2.fq.gz` is missing... So I look for it `grep -rn '2dN_KO_input_R2' scripts/`, mistake corrected in concat_2, and re-run manually for both 2dN_KO_input_R1 and 2dN_KO_input_R2. Will double check all is good at the fastqc step



# Trimming and Fastqc
```



```

## Fastqc on raw reads
```bash
# example for 1 file:
fastqc -o output/fastqc input/2dN_HET_H3K27me3_R1_1.fq.gz

# grouped:
sbatch scripts/fastqc_1.sh # 11262560 ok 
sbatch scripts/fastqc_2.sh # 11262458 2dN_KO_input_R2_2 fastqc fail (file was empty)
sbatch scripts/fastqc_3.sh # 11262515 ok
```

--> Re-generate the 2dN_KO_input_R2:
```bash
# Copy backup
cp input/backup/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L001_R2_001.fastq.gz input/
cp input/backup/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L002_R2_001.fastq.gz input/

# Concatenate backup
cat input/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L001_R2_001.fastq.gz input/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L002_R2_001.fastq.gz > input/2dN_KO_input_R2_2.fq.gz

# fastqc
fastqc -o output/fastqc input/2dN_KO_input_R2_2.fq.gz
```
--> 2dN_KO_input_R1 and 2dN_KO_input_R2 are not the same, so we're good

ESC_KO_H3K27me3_R2_2 and ESC_KO_H3K27me3_R2_1 are the exact same one according to fastqc... I did `grep -rn 'ESC_KO_H3K27me3_R2' scripts` and found the cat command was wrong for ESC_KO_H3K27me3_R2_2 file. Let's regenerate ESC_KO_H3K27me3_R2_2:
```bash
# Copy backup
cp input/backup/4-ESCs-EZH1-KO-2-H3k27me3_S3_L001_R2_001.fastq.gz input/
cp input/backup/4-ESCs-EZH1-KO-2-H3k27me3_S3_L002_R2_001.fastq.gz input/

# Concatenate backup
cat input/4-ESCs-EZH1-KO-2-H3k27me3_S3_L001_R2_001.fastq.gz input/4-ESCs-EZH1-KO-2-H3k27me3_S3_L002_R2_001.fastq.gz > input/ESC_KO_H3K27me3_R2_2.fq.gz


# fastqc
fastqc -o output/fastqc input/ESC_KO_H3K27me3_R2_2.fq.gz
```
--> All good now



## Fastqc on fastp-trimmed-reads
### Trimming with Fastp
```bash
# example for 1 file:
fastp -i input/2dN_HET_H3K27me3_R1_1.fq.gz -I input/2dN_HET_H3K27me3_R1_2.fq.gz \
      -o output/fastp/2dN_HET_H3K27me3_R1_1.fq.gz -O output/fastp/2dN_HET_H3K27me3_R1_2.fq.gz \
	  -h output/fastp/2dN_HET_H3K27me3_R1 -j output/fastp/2dN_HET_H3K27me3_R1

# grouped:
sbatch fastp_ESC.sh # 11262881 ok
sbatch fastp_NPC.sh # 11262882 ok 
sbatch fastp_2dN.sh # 11262874, 2dN_KO_input_R2 fail
```
Let's repeat fastp-triming for ESC_KO_H3K27me3_R2_2:
```bash
fastp -i input/ESC_KO_H3K27me3_R2_1.fq.gz -I input/ESC_KO_H3K27me3_R2_2.fq.gz \
    -o output/fastp/ESC_KO_H3K27me3_R2_1.fq.gz -O output/fastp/ESC_KO_H3K27me3_R2_2.fq.gz \
	-h output/fastp/ESC_KO_H3K27me3_R2 -j output/fastp/ESC_KO_H3K27me3_R2
```


### Fastqc trimmed reads
Scripts adapted so that they can only start when there respective `fastp.sh` job is finish:
```bash
# grouped:
sbatch --dependency=afterany:11262881 scripts/fastqc_fastp_ESC.sh # 11264158 ok
sbatch --dependency=afterany:11262882 scripts/fastqc_fastp_NPC.sh # 11264163 ok
sbatch --dependency=afterany:11262874 scripts/fastqc_fastp_2dN.sh # 11264174, 2dN_KO_input_R2
```
--> Fail, I took raw fastq input for fastqc, not the fastp-trim reads, script corrected and re-launch:
```bash
sbatch scripts/fastqc_fastp_ESC.sh # 11313011
sbatch scripts/fastqc_fastp_NPC.sh # 11313006
sbatch scripts/fastqc_fastp_2dN.sh # 11313014
```

--> Re-generate the 2dN_KO_input_R2 fastp and fastqc: 
```bash
sbatch fastp_fastqc_2dN_KO_input_R2.sh # 11312794
```
Let's repeat fastqc for fastp-triming ESC_KO_H3K27me3_R2_2 XXX:
```bash
fastqc -o output/fastqc/fastp output/fastp/ESC_KO_H3K27me3_R2_2.fq.gz
```


--> All is XXX










