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
sbatch fastqc_raw_NPC.sh # 10982864 complete
sbatch fastqc_raw_4wN.sh # 10983894 complete
sbatch fastqc_raw_8wN.sh # 10984036 complete
# By mistake, I replace output/fastqc by the fastp-fastqc. 
# Let's re-generate fastqc for raw files in output/fastqc/raw
sbatch fastqc_raw_ESC.sh # 11062004
sbatch fastqc_raw_2dN.sh # 11062006
sbatch fastqc_raw_NPC.sh # 11062007
sbatch fastqc_raw_4wN.sh # 11062008
sbatch fastqc_raw_8wN.sh # 11062009

### XXX Double check date for the new fastqc-fastp is correct, 20230313 for all XXX
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
sbatch scripts/fastp_raw_ESC.sh # 10980247 complete (important infos)
sbatch scripts/fastp_raw_NPC2dN.sh # 10984692 complete (important infos)
sbatch scripts/fastp_raw_4wN8wN.sh # 10984925 complete (important infos)
sbatch scripts/fastp_raw_8wN_miss.sh # 11033843 complete (important infos)
```
Run fastqc on fastp-trimmed files
```bash
sbatch scripts/fastqc_fastp.sh # 11034677 (cancelled "due to time limit")
sbatch scripts/fastqc_fastp_1.sh # 11060976 complete
```

# Mapping with STAR
##### 20230310, 20230313

## Index the genome
*NOTE: theorically optimal size for `--sjdbOverhang` is [max read length]-1, thus need create specific index for specific read size. But the effect is marginal according to the [creator](https://github.com/alexdobin/STAR/issues/931). So let's keep it default.*

hg19 genome with 12CPU and 50G mem (time=<1.5day)
```bash
module load STAR/2.7.3a-GCC-9.3.0
# command
STAR --runThreadN 12 \
	--runMode genomeGenerate \
	--genomeDir /scr1/users/roulet/Akizu_Lab/Master/meta/STAR_hg19 \
	--genomeFastaFiles /scr1/users/roulet/Akizu_Lab/Master/meta/male.hg19.fasta \
	--sjdbGTFfile /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf 

# Run in slurm
sbatch STAR_index_hg19.sh # 10982789 complete
```


### Untrimmed fastq
Keep standard parameter as adapted for mammalian genome. Some examples [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html) and [here](https://biocorecrg.github.io/RNAseq_course_2019/alnpractical.html).

Mapping
```bash
module load STAR/2.7.3a-GCC-9.3.0
# example for 1 file:
STAR --genomeDir ../../Master/meta/STAR_hg19/ \
	--runThreadN 12 \
	--readFilesCommand zcat \
	--readFilesIn input/NPC_WT_R1_1.fq.gz input/NPC_WT_R1_2.fq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix output/STAR/NPC_WT_

# Run time-per-time:
sbatch scripts/STAR_raw_NPC.sh # 11034673, slight test, output ok
sbatch scripts/STAR_raw_NPC_1.sh # 11061777, 11065401
sbatch scripts/STAR_raw_ESC.sh # 11061825, 11065404
sbatch scripts/STAR_raw_2dN.sh # 11061828, 11065418
sbatch scripts/STAR_raw_4wN.sh # 11061917, 11065420
sbatch scripts/STAR_raw_8wN.sh # 11061920, 11065422
```
Mapping indexation
```bash
# example for 1 file:
XXX
# Run time-per-time:
XXX
```
### Fastp-trimmed fastq
Mapping
```bash
module load STAR/2.7.3a-GCC-9.3.0
# example for 1 file:
STAR --genomeDir ../../Master/meta/STAR_hg19/ \
	--runThreadN 12 \
	--readFilesCommand zcat \
	--readFilesIn output/fastp/NPC_WT_R1_1.fq.gz output/fastp/NPC_WT_R1_2.fq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix output/STAR/fastp/NPC_WT_R1_

# Run time-per-time:
sbatch scripts/STAR_fastp_ESC.sh # 11062736, 11065467
sbatch scripts/STAR_fastp_NPC.sh # 11062742, 11065493
sbatch scripts/STAR_fastp_2dN.sh # 11062743, 11065515
sbatch scripts/STAR_fastp_4wN.sh # 11062744, 11065547
sbatch scripts/STAR_fastp_8wN.sh # 11062745, 11065561
```
Mapping indexation
```
XXX
```
*NOTE: next time do indexation after mapping; same script*
# Install IGV for vizualization
Go to [IGV](https://software.broadinstitute.org/software/igv/download) and copy link for Linux download.\
Go to `Master/software`
```bash
# download and install
wget https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_2.16.0_WithJava.zip
unzip IGV*
# Add shortcut so that we only use *igv.sh* to run it
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software/IGV_Linux_2.16.0
# Restart terminal
```




# Generate coverage (wiggle) files

XXX:
Coverage file:
FIND PARAMETER FOR HUMAN


--outWigType wiggle --outWigStrand Unstranded
Let's use bamCoverage from Deeptools suite within conda activate CondaGS to generate the bigwig files, that is better than using STAR as STAR cannot normalize based on TPM(=BPM) but RPM or raw reads...:
bamCoverage --bam data/mapping/raw/${x}_Aligned.sortedByCoord.out.bam --outFileName data/mapping/raw/${x}_Aligned.sortedByCoord.out.bigwig --outFileFormat bigwig --normalizeUsing BPM --binSize 10  



