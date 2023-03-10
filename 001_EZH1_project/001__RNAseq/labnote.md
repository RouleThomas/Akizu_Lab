# Import meta (genome) files
Files format to follow (ENCODE):\
**chr**:
`chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrM,chrX,chrY`\
**gene**:
`ENSG00000228253.1`

Files download from [ENCODE](https://www.encodeproject.org/data-standards/reference-sequences/) and cp to cluster: 
```bash
cp /home/roulet/tsclient/roule/Google\ Drive\ Streaming/Shared\ drives/akizulab/Personal\ Folders/Thomas/meta/* \
/scr1/users/roulet/Akizu_Lab/Master/meta`
```
ENCODE genome files in `/scr1/users/roulet/Akizu_Lab/Master/meta`:
- **GRCh38 fasta genome** = GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
- **GRCh38 gtf file** = ENCFF159KBI.gtf
- **hg19 fasta genome** = male.hg19.fasta
- **hg19 gtf file** = gencode.v19.annotation.gtf

# pipeline used from previous analyses
Infos from `Method_RNAseq, DEG, volcano, GSEA, heatmap.SZ.docx`:
- Raw fastq cleaned with *fastp* (adaptor removed, >10% poly-N sequences removed, low quality removed)
	- Q20, Q30, and GC content of the clean data were calculated.
    - Downstream analyses on the clean data with high quality.
		- Seems to be the default *fastp* parameter
 - Reads map on hg19 with *STAR* 
 - Reads count on gene feature with *featureCounts*
 - DEG with *DESEq2* (padj<0.05)
 - GSEA performed and plotted with *ClusterProfiler*


# Import files from Google drive to the cluster
##### 20230308, 20230309, 20230310
I cannot use a bash script to import files as if disconnection the transfer will fail. So cp the old-fashion way, let computer running o/n.\
**ESC, NPC, 2 days-neurons**
```bash
cp /home/roulet/tsclient/roule/Google\ Drive\ Streaming/Shared\ drives/akizulab/Primary\ Data/RNAseqs/EZH1\ RNAseq/1\ and\ 2\ month\ neuron\ RNAseq\ Aug2022/01.RawData/* \
/scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq/input
``` 
**1, 2 months-neurons**
```bash
cp -r /home/roulet/tsclient/roule/Google\ Drive\ Streaming/Shared\ drives/akizulab/Primary\ Data/RNAseqs/EZH1\ RNAseq/1\ and\ 2\ month\ neuron\ RNAseq\ Aug2022/01.RawData/ \
/scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq/input
``` 
# File renaiming
##### 20230309, 20230310
Made a custom bash script to rename each files (files are already compressed so I modified script `organize_raw.sh` to keep only renaming function). 
```bash
# example for 1 file:
outdir="input"

x="NPC_WT_R1_1"
raw_f="P_WT_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"

if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi
# Run command time-per-time (ESC, then 2dN, then PNC):
sbatch rename_raw_ESC.sh
sbatch rename_raw_2dN.sh
sbatch rename_raw_NPC.sh
```

# Quality control with FASTQC on raw fastq
##### 20230310
FASTQC is not an available module. Let's download it:
```bash
# Download in Master/Software/
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
# Unzip
unzip fastqc_v0.12.1.zip
# Add execution right to the file
chmod +x FastQC/fastqc
# Add shortcut so that we only use *fastqc* to run it
## backup .bashrc file in case
cp ~/.bashrc ~/.bashrc.backup
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software/FastQC
# Restart terminal
```
Run fastqc
```bash
# example for 1 file:
fastqc -o output/fastqc input/ESC_HET_R1_1.fq.gz

# Run time-per-time (ESC, then 2dN, then PNC):
sbatch fastqc_raw_ESC.sh # 10979372 complete
sbatch fastqc_raw_2dN.sh # 10980363 complete
sbatch fastqc_raw_NPC.sh # 10982864
sbatch fastqc_raw_4wN.sh # 10983894
sbatch fastqc_raw_4wN.sh # 10984036
```
Copy report to google drive
```bash
XXX wait help
```

# Quality control with FASTP (trim)
##### 20210310
Install [fastp](https://github.com/OpenGene/fastp).
```bash
# Download in Master/Software/
wget https://opengene.org/fastp/fastp
chmod a+x ./fastp
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software
# Restart terminal
```
Run fastp
```bash
# example for 1 file:
fastp -i input/ESC_WT_R1_1.fq.gz -I input/ESC_WT_R1_2.fq.gz \
      -o output/fastp/ESC_WT_R1_1.fq.gz -O output/fastp/ESC_WT_R1_2.fq.gz \
	  -h output/fastp/ESC_WT_R1 -j output/fastp/ESC_WT_R1

# Run time-per-time (ESC, then 2dN, then PNC):
sbatch scripts/fastp_raw_ESC.sh # 10980247
sbatch scripts/fastp_raw_NPC2dN.sh # 10984692
sbatch scripts/fastp_raw_4wN8wN.sh # 10984925
```
Run fastqc on fastp-trimmed files
```
XXX
```

# Mapping with STAR
##### 20230310

## Index the genome
*NOTE: theorically optimal size for `--sjdbOverhang` is [max read length]-1, thus need create specific index for specific read size. But the effect is marginal according to the [creator](https://github.com/alexdobin/STAR/issues/931). So let's keep it default.*

hg19 genome with 12CPU and 50G mem (time=XXX)
```bash
module load STAR/2.7.3a-GCC-9.3.0
# command
STAR --runThreadN 12 \
	--runMode genomeGenerate \
	--genomeDir /scr1/users/roulet/Akizu_Lab/Master/meta/STAR_hg19 \
	--genomeFastaFiles /scr1/users/roulet/Akizu_Lab/Master/meta/male.hg19.fasta \
	--sjdbGTFfile /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf 

# Run in slurm
sbatch STAR_index_hg19.sh # 10982789
```


### Untrimmed fastq
XXX https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html https://biocorecrg.github.io/RNAseq_course_2019/alnpractical.html XXX
```bash
module load STAR/2.7.3a-GCC-9.3.0
# example for 1 file:

# Run time-per-time (ESC, then 2dN, then PNC):

```



### Fastp-trimmed fastq







