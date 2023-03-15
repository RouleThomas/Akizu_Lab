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
##### 20230310, 20230313, 20230314
Data are unstranded.
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
```bash
# example for 1 file:
module load sam-bcf-tools
samtools index output/STAR/NPC_WT_Aligned.sortedByCoord.out.bam

# time per time (raw and fastp together):
sbatch STAR_index_NPC.sh # 11075019
sbatch STAR_index_ESC.sh # 11075041
sbatch STAR_index_2dN.sh # 11075067
sbatch STAR_index_4wN.sh # 11075079
sbatch STAR_index_8wN.sh # 11075108
```
*NOTE: next time do indexation after mapping; same script*

--> Let's compil the number of uniquely mapped reads for all files (add it in the Google Drive `RNAseq_infos.xlsx` file)
```bash
# Print nb of uniq map reads for raw mapping
for file in output/STAR/raw/*Log.final.out; do
    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $file | awk '{print $NF}')
    echo "$file: Number of uniquely mapped reads: $uniquely_mapped_reads"
done > output/STAR/raw/uniq_map_reads_counts_raw.txt

# Print nb of uniq map reads for fastp mapping
for file in output/STAR/fastp/*Log.final.out; do
    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $file | awk '{print $NF}')
    echo "$file: Number of uniquely mapped reads: $uniquely_mapped_reads"
done > output/STAR/fastp/uniq_map_reads_counts_fastp.txt
```
**More number of uniquely mapped reads for the fastp trimmed reads**, thus from now on, data analyses with the fastp-trimmed.
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

# Install Anaconda
##### 20230313
```bash
module load Python/3.9.6*
```
Go to `Master/software` 
```bash
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh`
bash Anaconda3-2022.10-Linux-x86_64.sh
# Follow installation; accept all
```
Anaconda3 installed to `/home/roulet/anaconda3`; no need to module load it.
# Generate coverage (wiggle) files
Create a **deeptools; conda environment**
```bash
conda create -n deeptools -c bioconda deeptools
```
Generate few files RPKM-normalized for comparison with novogene analyses:
- 4wN_WT_R1 (H9)
- 4wN_KO_R1 (8del)
- 4wN_HET_R1 (het5)
```bash
conda activate deeptools
# example for 1 file:
bamCoverage --bam output/STAR/fastp/4wN_WT_R1_Aligned.sortedByCoord.out.bam \
	--outFileName output/temp/4wN_WT_R1_Aligned.sortedByCoord.out.bigwig \
	--outFileFormat bigwig \ 
	--normalizeUsing RPKM \ 
	--binSize 10
# Run for all 3 files:
sbatch RPKM_wig.sh # 11075632
```
Generate coverage files **TPM-normalized** (=BPM) for all samples (fastp-trimmmed reads)

*NOTE: RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb)). BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions)*

```bash
conda activate deeptools
# example for 1 file with correct BPM-TPM normalization:
bamCoverage --bam output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam \
	--outFileName output/bigwig/${x}_Aligned.sortedByCoord.out.bw \
	--outFileFormat bigwig \
	--normalizeUsing BPM \
	--binSize 10   
# run time-per-time:
sbatch scripts/TPM_bw_NPC.sh # 11075894
sbatch scripts/TPM_bw_ESC.sh # 11075901
sbatch scripts/TPM_bw_2dN.sh # 11075962
sbatch scripts/TPM_bw_4wN.sh # 11075965
sbatch scripts/TPM_bw_8wN.sh # 11075967
```

# Count the reads on gene feature
##### 20230314
Create a **featurecounts; conda environment**
```bash
conda create -c bioconda -n featurecounts subread
conda activate featurecounts
```

Let's first do a comparison raw vs fastp to confirm we end up with more counts using fastp
```bash
# ESC_WT_R1 raw:
featureCounts -p -C -O \
	-P -B -d 30 -D 1000 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/featurecounts/ESC_WT_R1.txt output/STAR/raw/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 53%
# ESC_WT_R1 fastp:
featureCounts -p -C -O \
	-P -B -d 30 -D 1000 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/temp/ESC_WT_R1.txt output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 52%
```
--> fastp is again better.

Many fragment are unassigned for Fragment_Lenght. Samples not enough sonicated maybe? Let's try some fine-tuning:

```bash
# Both ends mapped required but fragment short
featureCounts -p -C -O \
	-P -B -d 50 -D 600\
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/temp/ESC_WT_R1.txt output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 46.9%
# Both ends mapped required but fragment long
featureCounts -p -C -O \
	-P -B -d 30 -D 2500 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/temp/ESC_WT_R1.txt output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 66.6%
# Both ends mapped required but fragment even long
featureCounts -p -C -O \
	-P -B -d 30 -D 10000 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/temp/ESC_WT_R1.txt output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 81.4%
# Both ends mapped not required 
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	o output/featurecounts/ESC_WT_R1.txt output/STAR/raw/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 87.2%
```
--> Seems that the longer fragment length allowed the better it is. (expected as no sonication in RNAseq; even though library prep focus for 300bp fragment)
- `-O` to count on meta feature (gene)
- `-p` (paired-end) with `-C` (not count paired-reads on 2 diff. chr.)
- NOT THIS: `-P -B -d 30 -D 1000` (set min and max paired reads to 30 to 1000bp)

Count on gene features with parameter
```bash
conda activate featurecounts

# example for 1 file:
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/featurecounts/${x}.txt output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam

# time per time:
sbatch scripts/featurecounts_ESC.sh # 11084182
sbatch scripts/featurecounts_NPC.sh # 11084181
sbatch scripts/featurecounts_2dN.sh # 11084178
sbatch scripts/featurecounts_4wN.sh # 11084176
sbatch scripts/featurecounts_8wN.sh # 11084168
```
# Quality control metrics
##### 20230314
Print number of succesfully assigned alignments for each sample (add to drive `RNAseq_infos.xlsx`)
```bash
for file in output/featurecounts/*.summary; do
    assigned_reads=$(grep "Assigned" $file | awk '{print $NF}')
    echo "$file: Assigned: $assigned_reads"
done > output/featurecounts/assigned_reads_counts.tsv
```
Print the total number of reads
```bash
for file in output/STAR/fastp/*.final.out; do
    input_reads=$(grep "Number of input reads" $file | awk '{print $NF}')
    echo "$file: Number of input reads: $input_reads"
done > output/STAR/fastp/input_reads_counts.txt
```
Add these values to the `RNAseq_infos.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >80% input reads as been assigned to (gene) features

# Install Bioconductor
##### 20230315
Create a **deseq2; conda environment**
```bash
conda create -n deseq2 -c bioconda bioconductor-deseq2
conda create -n deseq2 -c "bioconda/label/broken" bioconductor-deseq2
conda create -n deseq2 -c "bioconda/label/cf201901" bioconductor-deseq2
conda create -n deseq2 -c "bioconda/label/gcc7" bioconductor-deseq2
```
All these commands failed. So create environment with R/4.2.2 and install it within R
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
# PCA and clustering with deseq2
##### 20230315
*NOTE: Tons of deseq2 ressource [here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).*

Open R/4.2.2 with ressource:
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
In R; Import all counts and combine into one matrix, deseq2 dataframe (dds)
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "4wN_WT_R1", "4wN_WT_R2", "4wN_KO_R1",
   "4wN_KO_R2", "4wN_HET_R1", "4wN_HET_R2",
   "4wN_HET_R3", "4wN_HET_R4", "4wN_iPSCWT_R1",
   "4wN_iPSCWT_R2", "4wN_iPSCpatient_R1", "4wN_iPSCpatient_R2",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4", "8wN_iPSCpatient_R1", "8wN_iPSCpatient_R2",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype + time)

# Data normalization
vsd <- vst(dds, blind=TRUE)
rld <- rlog(dds, blind=TRUE) # last 5min

# Vizualization for quality metrics
## Heatmap of the sample-to-sample distances
### vsd 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$time, vsd$genotype, vsd$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("output/deseq2/heatmap_cluster_vsd.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

### rld
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$time, rld$genotype, rld$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("output/deseq2/heatmap_cluster_rld.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

## PCA
### vsd 
pdf("output/deseq2/PCA_vsd.pdf", width=10, height=10)
pcaData <- plotPCA(vsd, intgroup=c("time", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()
dev.off()

### rld 
pdf("output/deseq2/PCA_rld.pdf", width=10, height=10)
pcaData <- plotPCA(rld, intgroup=c("time", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()
dev.off()
```
*NOTE data normalization: There is two ways of normalizing data; **vsd** or **rld**. **rld** better as it accounts for the differences in library sizes and normalizes the data accordingly. However **vsd** is better (provide better variance stabilization) if working with very different library size or super-different expression levels. **Always good to test both.***

*NOTE blind parameter: Set to TRUE for unbiased quality assessment*

--> Here **rld** show better clustering

# DEGs with deseq2
##### 20230315
## 'one-by-one' comparison
Comparison WT vs mutant (KO or HET) for each time-points:
- NPC KO vs WT
- NPC HET vs WT
- ESC KO vs WT
- ESC HET vs WT
- 2dN KO vs WT
- 2dN HET vs WT
- 4wN KO vs WT
- 4wN HET vs WT
- 8wN KO vs WT
- 8wN HET vs WT
- 4wN iPSCpatient vs iPSCWT
### NPC KO vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")

## Export result as 'raw_NPC_KO_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_NPC_KO_vs_NPC_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_NPC_KO_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_NPC_KO_vs_NPC_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()


```
--> Use these data for comparison with novogene analyses

--> Save output and do venn diagram in the [venn_webtool](https://www.biovenn.nl/index.php). Files for comparison is `Google Drive///output/deseq2/Table_for_comparison_Novogene.xlsx`


### NPC HET vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_HET_vs_WT", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_NPC_HET_vs_NPC_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_NPC_HET_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_NPC_HET_vs_NPC_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
### ESC KO vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1" ,"ESC_KO_R2" ,"ESC_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")

## Export result as 'raw_ESC_KO_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_ESC_KO_vs_NPC_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_KO_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_ESC_KO_vs_NPC_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
### ESC HET vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_HET_R1" ,"ESC_HET_R2" ,"ESC_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_HET_vs_WT", type="apeglm")

## Export result as 'raw_ESC_HET_vs_ESC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_ESC_HET_vs_ESC_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_ESC_HET_vs_ESC_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
### 2dN KO vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2" ,"2dN_WT_R3",
   "2dN_KO_R1" ,"2dN_KO_R2" ,"2dN_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")

## Export result as 'raw_2dN_KO_vs_2dN_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_2dN_KO_vs_2dN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_2dN_KO_vs_2dN_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
### 2dN HET vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2" ,"2dN_WT_R3",
   "2dN_HET_R1" ,"2dN_HET_R2" ,"2dN_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_HET_vs_WT", type="apeglm")

## Export result as 'raw_2dN_HET_vs_2dN_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_2dN_HET_vs_2dN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_2dN_HET_vs_2dN_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
### 4wN KO vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_KO_R1",
   "4wN_KO_R2" )

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")

## Export result as 'raw_4wN_KO_vs_4wN_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_4wN_KO_vs_4wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_4wN_KO_vs_4wN_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
### 4wN HET vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_HET_R1", "4wN_HET_R2",
   "4wN_HET_R3" ,"4wN_HET_R4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_HET_vs_WT", type="apeglm")

## Export result as 'raw_4wN_HET_vs_4wN_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_4wN_HET_vs_4wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_4wN_HET_vs_4wN_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
--> This comparison result in very few DEGs, in agreement with clustering that shows they are similar... I suspect mis-annotation and this sample is not WT but HET.
### 8wN KO vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("8wN_WT_R1", "8wN_WT_R2" ,"8wN_WT_R3" ,"8wN_WT_R4" ,"8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")

## Export result as 'raw_8wN_KO_vs_8wN_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_8wN_KO_vs_8wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_8wN_KO_vs_8wN_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
--> XXX IT FAIL WITH WEIRD ERROR to repeat and solve later
### 8wN HET vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("8wN_WT_R1", "8wN_WT_R2" ,"8wN_WT_R3" ,"8wN_WT_R4" ,"8wN_HET_R1",
   "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_HET_vs_WT", type="apeglm")

## Export result as 'raw_8wN_HET_vs_8wN_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_8wN_HET_vs_8wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_8wN_HET_vs_8wN_WT.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
### 4wN iPSCpatient vs iPSCWT
Take ressource
```bash
module load R/4.2.2
srun --mem=50g --pty bash -l
R
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("4wN_iPSCWT_R1",
   "4wN_iPSCWT_R2" ,"4wN_iPSCpatient_R1", "4wN_iPSCpatient_R2")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_iPSCWT_vs_iPSCpatient", type="apeglm")

## Export result as 'raw_4wN_iPSCWT_vs_4wN_iPSCpatient.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_4wN_iPSCWT_vs_4wN_iPSCpatient.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_4wN_iPSCWT_vs_4wN_iPSCpatient.pdf", width=5, height=4)
### Identify DEGs and count them
res_df <- res %>% as.data.frame() %>% select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < 0 & res_df$padj == TRUE, na.rm = TRUE)
### Plot and add DEGs counts
plotMA(res_df, ylim=c(-2,2))
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = 1.75, labels = paste(n_upregulated), adj = c(0, 0.5), cex = 0.9)
text(x = max(res_df$baseMean, na.rm = TRUE) * 0.25, y = -1.75, labels = paste(n_downregulated), adj = c(0, 0.5), cex = 0.9)

dev.off()
```
