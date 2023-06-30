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

Re-installation new cluster, failed with `FindBin is missing`; need to install it with Perl
```bash
module load Perl/5.34.1-GCCcore-11.3.0
# Install a Perl module
cpanm FindBin


find /home/roulet/perl5/ -name \*FindBin\*



```
*NOTE: if error is link to perl and software missing use the same method to install Perl-dependent modules*

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
Data are unstranded.
## Index the genome
*NOTE: theorically optimal size for `--sjdbOverhang` is [max read length]-1, thus need create specific index for specific read size. But the effect is marginal according to the [creator](https://github.com/alexdobin/STAR/issues/931). So let's keep it default.*

hg19 genome with 12CPU and 50G mem (time=<1.5day)
```bash
module load STAR/2.7.3a-GCC-9.3.0
module load STAR/2.7.10b # New cluster
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
module load STAR/2.7.10b # New cluster
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
module load STAR/2.7.10b # New cluster
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
conda create -n deeptools -c bioconda deeptools # command can be launch from anywhere (directory and node)
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

**Re-install on the new cluster RES-RHEL-RH9HPC:**
Create a **deseq2; conda environment**
```bash
conda create -n deseq2 -c conda-forge r-base=4.2.2
conda install -c conda-forge r-ragg # needed as install.package(tidyverse) failed with this dependency
conda install -c conda-forge r-xml
annotate
geneplotter
```
Go in R and install deseq2:
```R
install.packages("tidyverse")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packges("pheatmap")
BiocManager::install("apeglm")
install.packges("factoextra")
BiocManager::install("rtracklayer")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("org.Hs.eg.db")
```
Works all good! I even added ChIPseeker in it!

# PCA and clustering with deseq2
*NOTE: Tons of deseq2 ressource [here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and [here](https://f1000research.com/articles/4-1070/v2).*

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

samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
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
pdf("output/deseq2/heatmap_cluster_vsd_no4wN.pdf", width=5, height=6)
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
pdf("output/deseq2/heatmap_cluster_rld_no4wN.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

## PCA
### vsd 
pdf("output/deseq2/PCA_vsd.pdf", width=10, height=10)
pdf("output/deseq2/PCA_vsd_no4wN.pdf", width=10, height=10)
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
pdf("output/deseq2/PCA_rld_no4wN.pdf", width=10, height=10)
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
- 4wN WT vs 4wN iPSCpatient
- 8wN WT vs 8wN iPSCpatient
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

## somthing weird with our samples, try different comparison
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_HET_R1", "4wN_HET_R2")
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_HET_R3" ,"4wN_HET_R4")
samples <- c("4wN_HET_R1", "4wN_HET_R2" ,"4wN_HET_R3" ,"4wN_HET_R4")

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

write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_4wN_HET_R3R4_vs_4wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_4wN_HET_R3R4_vs_4wN_WT.pdf", width=5, height=4)
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

--> Taking Replicate 3 and 4 for HET lead to a higher nb of DEGs. Let's check how similar are list of DEGs in 4wN vs 8wN with venn diagram in the [venn_webtool](https://www.biovenn.nl/index.php). Files for comparison is `Google Drive///output/deseq2/Table_for_comparison_4wN.xlsx`

### 8wN KO vs WT
Take ressource
```bash
srun --mem=50g --pty bash -l
module load R/4.2.2
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
--> IT FAILED at the shrinkage normalization
```bash
 # res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")
 error in optimHess(par=init,fn+nbinomFn,gr=nbinomGr,x=x,y=y,:non-finite value supplied by optim
 ```
--> Re-run the pipeline with >=10 minimum reads works, then re-run it again with >=5 works... 
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
--> The 4wN samples are weird, few DE, bad clustering, there may have been a sample missanotation? Let's put aside 4wN from now on

### 4wN WT vs 4wN iPSCpatient
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
samples <- c("4wN_WT_R1", "4wN_WT_R2", "4wN_iPSCpatient_R1",
   "4wN_iPSCpatient_R2")

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
res <- lfcShrink(dds, coef="genotype_iPSCpatient_vs_WT", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_4wN_iPSCpatient_vs_4wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_NPC_HET_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_4wN_iPSCpatient_vs_4wN_WT.pdf", width=5, height=4)
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

### 8wN WT vs 8wN iPSCpatient
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
samples <- c("8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_iPSCpatient_R1",
   "8wN_iPSCpatient_R2")

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
res <- lfcShrink(dds, coef="genotype_iPSCpatient_vs_WT", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/raw_8wN_iPSCpatient_vs_8wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_NPC_HET_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2/plotMA_res_8wN_iPSCpatient_vs_8wN_WT.pdf", width=5, height=4)
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
--> Compare WT vs iPSC patient within the xcell file in Google `output/deseq2/Table_for_comparison_iPSC.xlsx`







# Time-course with deseq2
## ESC, NPC, 2dN, 8wN KO and HET vs WT
Test whether there are genotype-specific differences over time (ESC, NPC, 2dN, 8wN); general trend of difference over time. 

### Identification of the significant genes in Time-Course analyses across genotypes 
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
library("factoextra")

# import all except 4wN
## collect all samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
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
write.csv(counts_all, file="output/deseq2/counts_all.txt")
### If need to import: counts_all <- read_csv("output/deseq2/counts_all.txt") #To import


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

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet Time-Course 
### desgin = full-model
ddsTC <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                                colData = coldata,
                                design = ~ genotype + time + genotype:time)

### Define the reduced model (not including interaction term for comparison with full model)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ genotype + time)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)

## Save the time-course analyses
## Export result as 'raw_4wN_iPSCWT_vs_4wN_iPSCpatient.txt'
write.csv(resTC %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2/resTC.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import
```
*NOTE: Here we the full model is` ~ genotype + time + genotype:time` and the reduced model do not include the interaction term `genotype:time`. Then we compare by LRT the full vs reduced model to assess the interaction effect, thus identify genes that behave differentially between genotypes over the time-course.*

**Gene expression profile** using normalized counts (deseq2); more adapted than TPM here as directly related to the deseq2-time-course analyses for **identification of the 'optimal' qvalue treshold**

Here we will test different qvalue treshold and look at the 20 last genes. We will select a qvalue treshold that show clear separation of genotypes.

```R
# Normalized counts with deseq2 and tidy it
normalized_counts <- as_tibble(counts(ddsTC, normalized = TRUE), rownames = "gene") %>%
  gather(key = "sample", value = "norm_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

normalized_counts$time <-
  factor(normalized_counts$time,
         c("ESC", "NPC", "2dN", "8wN"))
normalized_counts$genotype <-
  factor(normalized_counts$genotype,
         c("WT", "KO", "HET"))

# Gather the 10 first significant deseq2-TC genes
significant_deseq2TC_genes <- as_tibble(resTC, rownames = "gene") %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  select(gene) %>%
  left_join(normalized_counts)

# Stat
stat_significant_deseq2TC_genes <- significant_deseq2TC_genes %>%
  select(-replicate) %>%
  group_by(gene, time, genotype) %>% summarise(mean=mean(norm_counts), median= median(norm_counts), SD=sd(norm_counts), n=n(), SE=SD/sqrt(n)) 	



pdf("output/deseq2/deseq2_TC_Top10genes.pdf", width=11, height=6)
stat_significant_deseq2TC_genes %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color=genotype), size=0.75) +
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width=.2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(~gene, nrow = 2, scale = "free")  +	
  xlab(label = "deseq2 normalized counts") +
  ggtitle("Top 10 Significant Genes in Time-Course Analysis Across Genotypes") +
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()


# Find optimal qvalue treshold
# Gather the 20 LAST p0.05 significant deseq2-TC genes
significant_deseq2TC_genes <- as_tibble(resTC, rownames = "gene") %>%
  filter(padj <= 0.0005) %>% # !!! Here we change this number and re-run the script !!!
  arrange(desc(padj)) %>%
  slice_head(n = 20) %>%
  select(gene) %>%
  left_join(normalized_counts)


  
# Stat
stat_significant_deseq2TC_genes <- significant_deseq2TC_genes %>%
  select(-replicate) %>%
  group_by(gene, time, genotype) %>% summarise(mean=mean(norm_counts), median= median(norm_counts), SD=sd(norm_counts), n=n(), SE=SD/sqrt(n)) 	



pdf("output/deseq2/deseq2_TC_Last20genes_p0.0005.pdf", width=11, height=6) # !!! Here we change the file name accordingly !!!
stat_significant_deseq2TC_genes %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color=genotype), size=0.75) +
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width=.2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(~gene, nrow = 4, scale = "free")  +	
  xlab(label = "deseq2 normalized counts") +
  ggtitle("Last 20 p0.0005 Significant Genes in Time-Course Analysis Across Genotypes") + # !!! Here we change the file name accordingly !!!
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()




```
--> We tested different qvalue treshold and look at the less significant ones:
- p0.05 (n = 14,322 genes): A lot so not show really different pattern (but some do)
- p0.01 (n = 11,578 genes): A lot so not show really different pattern (but some do)
- p0.005 (n = 10,673 genes): A lot so not show really different pattern (but some do)
- p0.001 (n = 8,974 genes): Here we majoritaly have different patterns 
- p0.0005 (n = 8,404 genes): Here we majoritaly have different patterns

--> p0.001 looks good



### Clustering of the significant genes in Time-Course analyses across genotypes using deseq2 vst norm-counts
Prepare the data, use deseq2 norm counts vst or rlog transformed (test both)
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")

# Import list of TC-DEGs
resTC <- read_csv("output/deseq2/resTC.txt") 

# Normalized counts with deseq2 and tidy it
normalized_counts <- as_tibble(counts(ddsTC, normalized = TRUE), rownames = "gene") %>%
  gather(key = "sample", value = "norm_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

normalized_counts$time <-
  factor(normalized_counts$time,
         c("ESC", "NPC", "2dN", "8wN"))
normalized_counts$genotype <-
  factor(normalized_counts$genotype,
         c("WT", "KO", "HET"))

# Clustering
## transform the norm counts
vst_counts <- vst(ddsTC)
rlog_counts <- rlog(ddsTC) # take 5 min

## Extract the transformed counts
vst_counts_matrix <- assay(vst_counts)
rlog_counts_matrix <- assay(rlog_counts)

## Isolate the signifificant genes by gene id
signif_TC_genes = as_tibble(resTC, rownames = "gene") %>%
  filter(padj<= 0.001) %>% # !!! Here can play with the qvalue !!!
  select(gene) %>%
  unique() 

## Convert into vector 
signif_TC_genes_vector <- signif_TC_genes$gene

## Filter our matrix with the significant genes
vst_counts_matrix_sig <- vst_counts_matrix[rownames(vst_counts_matrix) %in% signif_TC_genes_vector, ]
nrow(vst_counts_matrix_sig) # double-check the nb of genes is same s in signif_TC_genes
rlog_counts_matrix_sig <- rlog_counts_matrix[rownames(rlog_counts_matrix) %in% signif_TC_genes_vector, ]

## Reorder columns
ordered_columns <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3", "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3", "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3", "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3", "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3", "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3", "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4")

vst_counts_matrix_sig_ordered <- vst_counts_matrix_sig[, ordered_columns]

# Raw heatmap
pdf("output/deseq2/raw_pheatmap_vst_p0.001.pdf", width=8, height=10) # !!! Here change qvalue accordingly !!!
pheatmap(vst_counts_matrix_sig_ordered, 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         show_rownames = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

# Clustered heatmap with 10 clusters
pdf("output/deseq2/clustered_pheatmap_vst_p0.001_20cl.pdf", width=8, height=10)  # !!! Here change qvalue accordingly !!!
pheatmap(vst_counts_matrix_sig_ordered, 
                    scale = "row", 
                    clustering_distance_rows = "euclidean", 
                    clustering_distance_cols = "euclidean", 
                    clustering_method = "complete", 
                    show_rownames = FALSE,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    cutree_rows = 20,                                   # !!! Here change tree accordingly !!!
                    annotation_colors = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Set1"))(10))

dev.off()
```
*NOTE: Here for better vizualisation, stabilize the variance accross the multiple conditions, we transform our norm counts with  vst or rlog.*

It seems that is not the best way to vizualize the data, I get trouble comparing the data; instead let's use **line plots** and not heatmap.

```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")

# Import files to generete the DESeq2Dataset
## samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")
## Import counts_all and transform to matrix
counts_all <- read_csv("output/deseq2/counts_all.txt") %>% 
  select(-1)

### Transform merged_data into a matrix
#### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
#### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet Time-Course 
### desgin = full-model
ddsTC <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                                colData = coldata,
                                design = ~ genotype + time + genotype:time)

### Define the reduced model (not including interaction term for comparison with full model)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ genotype + time)
resTC <- results(ddsTC)


# Normalize the counts
## Normalized counts with deseq2 and tidy it
normalized_counts <- as_tibble(counts(ddsTC, normalized = TRUE), rownames = "gene") %>%
  gather(key = "sample", value = "norm_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

normalized_counts$time <-
  factor(normalized_counts$time,
         c("ESC", "NPC", "2dN", "8wN"))
normalized_counts$genotype <-
  factor(normalized_counts$genotype,
         c("WT", "KO", "HET"))

# Clustering
## transform the norm counts
vst_counts <- vst(ddsTC)
rlog_counts <- rlog(ddsTC) # take 5 min

## Extract the transformed counts
vst_counts_matrix <- assay(vst_counts)
rlog_counts_matrix <- assay(rlog_counts) # take ???

## Isolate the significant genes by gene id
signif_TC_genes = as_tibble(resTC, rownames = "gene") %>%
  filter(padj<= 0.05) %>%               # !!! Here change qvalue accordingly !!!
  select(gene) %>%
  unique() 

## Convert into vector 
signif_TC_genes_vector <- signif_TC_genes$gene

## Filter our matrix with the significant genes
vst_counts_matrix_sig <- vst_counts_matrix[rownames(vst_counts_matrix) %in% signif_TC_genes_vector, ]
nrow(vst_counts_matrix_sig) # double-check the nb of genes is same s in signif_TC_genes
rlog_counts_matrix_sig <- rlog_counts_matrix[rownames(rlog_counts_matrix) %in% signif_TC_genes_vector, ]

## Reorder columns
ordered_columns <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3", "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3", "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3", "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3", "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3", "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3", "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4")

vst_counts_matrix_sig_ordered <- vst_counts_matrix_sig[, ordered_columns]



# Define computationaly the optimal nb of cluster
## Transpose the matrix
vst_counts_matrix_sig_ordered_t <- t(vst_counts_matrix_sig_ordered)

## ELBOW method
## Determine the optimal number of clusters using the elbow method
pdf("output/deseq2/Elbow_cluster_p0.05.pdf", width=14, height=20)       
fviz_nbclust(vst_counts_matrix_sig_ordered_t, FUNcluster = hcut, method = "wss") +
  geom_vline(xintercept = 50, linetype = 2, color = "red") +
  labs(subtitle = "Elbow method")
dev.off()

optimal_clusters <- 5
### 5 is too low!!!

## silhouette method
pdf("output/deseq2/Silhouette_cluster_p0.05.pdf", width=14, height=20)       
fviz_nbclust(vst_counts_matrix_sig_ordered_t, FUNcluster = hcut, method = "silhouette") +
  geom_vline(xintercept = 50, linetype = 2, color = "red") +
  labs(subtitle = "Average Silhouette method")
dev.off()

## Set the optimal number of clusters (manually based on the elbow method plot)
optimal_clusters <- 3
### 3 is too low!!!




# Make a clean table with significant deseq2-TC genes
vst_counts_tidy <- as_tibble(vst_counts_matrix_sig_ordered, rownames = "gene") %>%
  gather(key = "sample", value = "vst_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(vst_counts_matrix_sig_ordered, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into 10 clusters
row_clusters <- cutree(row_hclust, k = 25)                   # !!! Here change tree nb accordingly !!!

## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(vst_counts_matrix_sig_ordered),
                           cluster = row_clusters)

## Compil both
vst_counts_tidy <- vst_counts_tidy %>%
  left_join(cluster_gene, by = "gene")

vst_counts_tidy$time <-
  factor(vst_counts_tidy$time,
         c("ESC", "NPC", "2dN", "8wN"))
vst_counts_tidy$genotype <-
  factor(vst_counts_tidy$genotype,
         c("WT", "KO", "HET"))

# Plot vst_transform norm deseq2 count with 'loess' method
pdf("output/deseq2/line_vst_p0.05_cl20.pdf", width=14, height=20)           # !!! Here change title accordingly !!!
ggplot(vst_counts_tidy, aes(x = time, y = vst_counts, color = genotype, group = genotype)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.8) +
  facet_wrap(~cluster, scale = "free")
dev.off()



# Plot vst_transform norm deseq2 count with 'loess' method pretty


## Calculate the number of genes per cluster
genes_per_cluster <- vst_counts_tidy %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2/line_vst_p0.05_cl25_pretty.pdf", width=20, height=14)
ggplot(vst_counts_tidy, aes(x = time, y = vst_counts)) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  # Add mean value for each genotype at each time point
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  # Add standard error bars around the mean points
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  # Add number of genes per cluster to the facet_wrap panels
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()



# Display some genes expression profile cluster-per-cluster

## Gather the 10 first significant deseq2-TC genes of specified cluster
significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
  inner_join(cluster_gene) %>%
  filter(cluster == 22) %>% # !!! change cluster nb !!!
  arrange(padj) %>%
  select(gene,padj) %>%
  left_join(normalized_counts) 

significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
select(gene) %>%
unique() %>%
slice_head(n=4)
  
significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
  inner_join(significant_deseq2TC_genes_cluster_genes)


## Stat
stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
  select(-replicate) %>%
  group_by(gene, time, genotype) %>% summarise(mean=mean(norm_counts), median= median(norm_counts), SD=sd(norm_counts), n=n(), SE=SD/sqrt(n)) 	



pdf("output/deseq2/deseq2_TC_Top4genes_cluster22.pdf", width=11, height=6)  # !!! change cluster nb !!!
stat_significant_deseq2TC_genes %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color=genotype), size=0.75) +
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width=.2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(~gene, nrow = 1, scale = "free")  +	
  xlab(label = "deseq2 normalized counts") +
  ggtitle("Top 4 significant genes in cluster22") +   # !!! change cluster nb !!!
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()



# Display genes expression profile for all cluster
## Non-vst norm counts
## Create a function to generate data for the top 2 significant genes for a given cluster
get_top2_genes_data <- function(cluster_number) {
  significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
    inner_join(cluster_gene) %>%
    filter(cluster == cluster_number) %>%
    arrange(padj) %>%
    select(gene, padj, cluster) %>%
    left_join(normalized_counts)
  
  significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
    select(gene) %>%
    unique() %>%
    slice_head(n = 2)
  
  significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
    inner_join(significant_deseq2TC_genes_cluster_genes)
  
  stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
    select(-replicate) %>%
    group_by(gene, time, genotype, cluster) %>% summarise(mean = mean(norm_counts), median = median(norm_counts), SD = sd(norm_counts), n = n(), SE = SD / sqrt(n))
  
  return(stat_significant_deseq2TC_genes)
}

# Generate data for all 25 clusters
all_clusters_data <- map_dfr(1:25, get_top2_genes_data)

# Plot top 2 significant genes for each of the 25 clusters
pdf("output/deseq2/deseq2_TC_Top2genes_AllClusters.pdf", width = 30, height = 14)
all_clusters_data %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color = genotype), size = 0.75) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(cluster ~ gene, nrow = 3, labeller = labeller(gene = function(x) paste("Cluster", unlist(all_clusters_data[1, "cluster"]), x)), scales = "free") +
  xlab(label = "deseq2 normalized counts") +
  ggtitle("Top 2 significant genes in each cluster") +
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()




## vst norm counts
## Create a function to generate data for the top 2 significant genes for a given cluster
get_top2_genes_data <- function(cluster_number) {
  significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
    inner_join(cluster_gene) %>%
    filter(cluster == cluster_number) %>%
    arrange(padj) %>%
    select(gene, padj, cluster) %>%
    left_join(vst_counts_tidy) # here vst counts
  
  significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
    select(gene) %>%
    unique() %>%
    slice_head(n = 2)
  
  significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
    inner_join(significant_deseq2TC_genes_cluster_genes)
  
  stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
    select(-replicate) %>%
    group_by(gene, time, genotype, cluster) %>% summarise(mean = mean(vst_counts), median = median(vst_counts), SD = sd(vst_counts), n = n(), SE = SD / sqrt(n))
  
  return(stat_significant_deseq2TC_genes)
}

# Generate data for all 25 clusters
all_clusters_data <- map_dfr(1:25, get_top2_genes_data)

# Plot top 2 significant genes for each of the 25 clusters
pdf("output/deseq2/deseq2_TC_Top2genes_AllClusters_vst_counts.pdf", width = 30, height = 14)
all_clusters_data %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color = genotype), size = 0.75) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(cluster ~ gene, nrow = 3, labeller = labeller(gene = function(x) paste("Cluster", unlist(all_clusters_data[1, "cluster"]), x)), scales = "free") +
  xlab(label = "vst-norm deseq2 normalized counts") +
  ggtitle("Top 2 significant genes in each cluster") +
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()





## tpm counts 
### Load tpm
tpm_all <- read_csv("output/tpm/tpm_all.txt") %>% 
  select(-1) #To import


## Transform tpm tibble into matrix
tpm_all_matrix = make_matrix(select(tpm_all, -Geneid), pull(tpm_all, Geneid)) 

tpm_all_raw_tidy <- as_tibble(tpm_all_matrix, rownames = "gene") %>%
  gather(key = "sample", value = "tpm", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

tpm_all_raw_tidy$time <-
  factor(tpm_all_raw_tidy$time,
         c("ESC", "NPC", "2dN", "8wN"))
tpm_all_raw_tidy$genotype <-
  factor(tpm_all_raw_tidy$genotype,
         c("WT", "KO", "HET"))

## Create a function to generate data for the top 2 significant genes for a given cluster
get_top2_genes_data <- function(cluster_number) {
  significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
    inner_join(cluster_gene) %>%
    filter(cluster == cluster_number) %>%
    arrange(padj) %>%
    select(gene, padj, cluster) %>%
    left_join(tpm_all_raw_tidy) # here tpm
  
  significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
    select(gene) %>%
    unique() %>%
    slice_head(n = 2)
  
  significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
    inner_join(significant_deseq2TC_genes_cluster_genes)
  
  stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
    select(-replicate) %>%
    group_by(gene, time, genotype, cluster) %>% summarise(mean = mean(tpm), median = median(tpm), SD = sd(tpm), n = n(), SE = SD / sqrt(n))
  
  return(stat_significant_deseq2TC_genes)
}

# Generate data for all 25 clusters
all_clusters_data <- map_dfr(1:25, get_top2_genes_data)

# Plot top 2 significant genes for each of the 25 clusters
pdf("output/deseq2/deseq2_TC_Top2genes_AllClusters_tpm.pdf", width = 30, height = 14)
all_clusters_data %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color = genotype), size = 0.75) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(cluster ~ gene, nrow = 3, labeller = labeller(gene = function(x) paste("Cluster", unlist(all_clusters_data[1, "cluster"]), x)), scales = "free") +
  xlab(label = "TPM") +
  ggtitle("Top 2 significant genes in each cluster") +
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()






# Plot vst_transform norm deseq2 count scale between -1 and +1

## Calculate the min and max vst_counts values
min_vst_counts <- min(vst_counts_tidy$vst_counts)
max_vst_counts <- max(vst_counts_tidy$vst_counts)
## Normalize the vst_counts values between -1 and +1
vst_counts_tidy <- vst_counts_tidy %>%
  mutate(vst_counts_normalized = 2 * ((vst_counts - min_vst_counts) / (max_vst_counts - min_vst_counts)) - 1)
pdf("output/deseq2/line_vst_p0.001_cl50_scale.pdf", width=14, height=20)  # !!! Here change title accordingly !!!
ggplot(vst_counts_tidy, aes(x = time, y = vst_counts_normalized, color = genotype, group = genotype)) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~cluster)+
  scale_y_continuous(limits = c(-1.5, 1.5))
dev.off()
## --> Scaling looks like shit, maybe too much samples and variability so that scaling get rid of a lot of information
```
*NOTE: I only tried vst normalization but I could try rlog normalization too.*

--> I tried computationaly defined (Elbow and Silhouette method) the optimal nb of cluster but propose 3 or 5 LOL

--> Overall 25 clusters with a 0.001qval looks optimal (all cluster display with unique biologically relevant profile)

--> I play with the span parameter. If reduce it less 'approximate' the trend. .8 look optimal. Even though I added boxplot so that we really see where are the data

--> *loess* is a good method as not linear and do not require numeric x-axis (*gam* method required it, and the code is buggy)

--> Looking at individual genes, we sometime do not really have the pattern show in geom_smooth. 

Lets; try to be closer to our data and work with **TPM for clustering, and vizualization**


### Clustering of the significant genes in Time-Course analyses across genotypes using TPM
#### Generate TPM-normalized counts

Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
# Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch featurecounts_TPM.sh # 11318014; 11318576 ok
# mv all output to output/tpm or rpkm folder
mv output/featurecounts/*tpm* output/tpm/
mv output/featurecounts/*rpkm* output/rpkm/
```

#### Clustering with TPM
```bash
module load R/4.2.2
```

```R
# Import TPM reads
## collect all samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/tpm/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR.fastp.")) %>%
    rename(!!sample := starts_with("output.STAR.fastp."))
}

## Merge all dataframe into a single one
tpm_all <- reduce(sample_data, full_join, by = "Geneid")
write.csv(tpm_all, file="output/tpm/tpm_all.txt")
### If need to import: tpm_all <- read_csv("output/tpm/tpm_all.txt") #To import


## Transform tpm tibble into matrix
tpm_all_matrix = make_matrix(select(tpm_all, -Geneid), pull(tpm_all, Geneid)) 

tpm_all_raw_tidy <- as_tibble(tpm_all_matrix, rownames = "gene") %>%
  gather(key = "sample", value = "tpm", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

tpm_all_raw_tidy$time <-
  factor(tpm_all_raw_tidy$time,
         c("ESC", "NPC", "2dN", "8wN"))
tpm_all_raw_tidy$genotype <-
  factor(tpm_all_raw_tidy$genotype,
         c("WT", "KO", "HET"))


## Isolate the significant genes by gene id
signif_TC_genes = as_tibble(resTC, rownames = "gene") %>%
  filter(padj<= 0.05) %>%               # !!! Here change qvalue accordingly !!!
  select(gene) %>%
  unique() 

## Convert into vector 
signif_TC_genes_vector <- signif_TC_genes$gene

## Filter our matrix with the significant genes
tpm_all_matrix_sig <- tpm_all_matrix[rownames(tpm_all_matrix) %in% signif_TC_genes_vector, ]


## Reorder columns
ordered_columns <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3", "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3", "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3", "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3", "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3", "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3", "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4")

tpm_all_matrix_sig_ordered <- tpm_all_matrix_sig[, ordered_columns]



# Make a clean table with significant deseq2-TC genes
tpm_all_tidy <- as_tibble(tpm_all_matrix_sig_ordered, rownames = "gene") %>%
  gather(key = "sample", value = "tpm", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")


# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(tpm_all_matrix_sig_ordered, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into 10 clusters
row_clusters <- cutree(row_hclust, k = 25)                   # !!! Here change tree nb accordingly !!!

## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(tpm_all_matrix_sig_ordered),
                           cluster = row_clusters)

## Compil both
tpm_all_tidy <- tpm_all_tidy %>%
  left_join(cluster_gene, by = "gene")

tpm_all_tidy$time <-
  factor(tpm_all_tidy$time,
         c("ESC", "NPC", "2dN", "8wN"))
tpm_all_tidy$genotype <-
  factor(tpm_all_tidy$genotype,
         c("WT", "KO", "HET"))

# Plot tpm count with 'loess' method
pdf("output/deseq2/line_tpm_p0.05_cl10.pdf", width=14, height=20)           # !!! Here change title accordingly !!!
ggplot(tpm_all_tidy, aes(x = time, y = tpm, color = genotype, group = genotype)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.8) +
  facet_wrap(~cluster, scale = "free")
dev.off()

# Plot tpm count with 'loess' method pretty
# Calculate the mean expression for each gene at each time point
mean_tpm_all_raw_tidy <- tpm_all_tidy %>%
  group_by(gene, time, genotype, cluster) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup()

genes_per_cluster <- tpm_all_tidy %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2/line_tpm_p0.05_cl25_pretty.pdf", width=20, height=14)
ggplot(mean_tpm_all_raw_tidy, aes(x = time, y = mean_tpm)) +
#  geom_line(aes(color = genotype, group = interaction(gene, genotype)), alpha = 0.1) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()


```

--> The tpm clustering poorly separate the pattern. Bad method here.

--> From now on, the better is to do vst-norm for clustering and present tpm values gene-per-gene

Let's try rlog normalization for clustering, to see if it better represent our data (show tpm for gene verification)



### Clustering of the significant genes in Time-Course analyses across genotypes using deseq2 rlog norm-counts
```bash
module load R/4.2.2
```
The below code is clean to generate clustering plot (qvalue and number of cluster can be adjusted), vizualize them and look at invividual gene expressin (rlog counts and tpm)
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")

# Import files to generete the DESeq2Dataset
## samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")
## Import counts_all and transform to matrix
counts_all <- read_csv("output/deseq2/counts_all.txt") %>% 
  select(-1)

### Transform merged_data into a matrix
#### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
#### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet Time-Course 
### desgin = full-model
ddsTC <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                                colData = coldata,
                                design = ~ genotype + time + genotype:time)

### Define the reduced model (not including interaction term for comparison with full model)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ genotype + time)
resTC <- results(ddsTC)


# Normalize the counts
## Normalized counts with deseq2 and tidy it
normalized_counts <- as_tibble(counts(ddsTC, normalized = TRUE), rownames = "gene") %>%
  gather(key = "sample", value = "norm_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

normalized_counts$time <-
  factor(normalized_counts$time,
         c("ESC", "NPC", "2dN", "8wN"))
normalized_counts$genotype <-
  factor(normalized_counts$genotype,
         c("WT", "KO", "HET"))

# Clustering

rlog_counts <- rlog(ddsTC) # take 5 min

## Extract the transformed counts
rlog_counts_matrix <- assay(rlog_counts) 

## Isolate the significant genes by gene id
signif_TC_genes = as_tibble(resTC, rownames = "gene") %>%
  filter(padj<= 0.05) %>%               # !!! Here change qvalue accordingly !!!
  select(gene) %>%
  unique() 

## Convert into vector 
signif_TC_genes_vector <- signif_TC_genes$gene

## Filter our matrix with the significant genes
rlog_counts_matrix_sig <- rlog_counts_matrix[rownames(rlog_counts_matrix) %in% signif_TC_genes_vector, ]
 # double-check the nb of genes is same s in signif_TC_genes
nrow(rlog_counts_matrix_sig)

## Reorder columns
ordered_columns <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3", "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3", "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3", "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3", "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3", "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3", "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4")

rlog_counts_matrix_sig_ordered <- rlog_counts_matrix_sig[, ordered_columns]


# Make a clean table with significant deseq2-TC genes
rlog_counts_tidy <- as_tibble(rlog_counts_matrix_sig_ordered, rownames = "gene") %>%
  gather(key = "sample", value = "rlog_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(rlog_counts_matrix_sig_ordered, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into k clusters
row_clusters <- cutree(row_hclust, k = 25)                   # !!! Here change tree nb accordingly !!!

## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(rlog_counts_matrix_sig_ordered),
                           cluster = row_clusters)


### Save dataframe
write.csv(cluster_gene, file="output/deseq2/cluster_gene_rlog_25cl.txt")


## Compil both
rlog_counts_tidy <- rlog_counts_tidy %>%
  left_join(cluster_gene, by = "gene")

rlog_counts_tidy$time <-
  factor(rlog_counts_tidy$time,
         c("ESC", "NPC", "2dN", "8wN"))
rlog_counts_tidy$genotype <-
  factor(rlog_counts_tidy$genotype,
         c("WT", "KO", "HET"))



# Plot rlog_transform norm deseq2 count with 'loess' method pretty
## Calculate the number of genes per cluster
genes_per_cluster <- rlog_counts_tidy %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2/line_rlog_p0.05_cl25_pretty_span07.pdf", width=20, height=14)      # !!! Here change title accordingly !!!
ggplot(rlog_counts_tidy, aes(x = time, y = rlog_counts)) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.7) +
  # Add mean value for each genotype at each time point
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  # Add standard error bars around the mean points
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  # Add number of genes per cluster to the facet_wrap panels
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()

# Plot rlog_transform norm deseq2 count with 'loess' method pretty with a light transparent color
# Calculate the mean expression for each gene at each time point
mean_rlog_counts_tidy <- rlog_counts_tidy %>%
  group_by(gene, time, genotype, cluster) %>%
  summarise(mean_rlog_counts = mean(rlog_counts)) %>%
  ungroup()

genes_per_cluster <- rlog_counts_tidy %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2/line_rlog_p0.05_cl25_pretty_span08_genes_bg.pdf", width=20, height=14)
ggplot(mean_rlog_counts_tidy, aes(x = time, y = mean_rlog_counts)) +
  geom_line(aes(color = genotype, group = interaction(gene, genotype)), alpha = 0.1) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()


pdf("output/deseq2/line_rlog_p0.05_cl25_pretty_span08_genes_WT-KO.pdf", width=20, height=14)
mean_rlog_counts_tidy %>%
  filter(genotype %in% c("WT","KO")) %>%
ggplot(., aes(x = time, y = mean_rlog_counts)) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()

# Plot TPM count with 'loess' method pretty

# Calculate the mean expression for each gene at each time point
mean_tpm_all_raw_tidy <- tpm_all_raw_tidy %>%
  inner_join(cluster_gene, by = "gene") %>%
  group_by(gene, time, genotype, cluster) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup()

genes_per_cluster <- tpm_all_raw_tidy %>%
  inner_join(cluster_gene, by = "gene") %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2/line_rlog_p0.05_cl25_pretty_span08_genes_bg_tpm.pdf", width=20, height=14)
ggplot(mean_tpm_all_raw_tidy, aes(x = time, y = mean_tpm)) +
  geom_line(aes(color = genotype, group = interaction(gene, genotype)), alpha = 0.1) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()



# Display some genes expression profile cluster-per-cluster

## Gather the 10 first significant deseq2-TC genes of specified cluster
significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
  inner_join(cluster_gene) %>%
  filter(cluster == 4) %>% # !!! change cluster nb !!!
  arrange(padj) %>%
  select(gene,padj) %>%
  left_join(rlog_counts_tidy) 

significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
select(gene) %>%
unique() %>%
slice_head(n=4)
  
significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
  inner_join(significant_deseq2TC_genes_cluster_genes)


## Stat
stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
  select(-replicate) %>%
  group_by(gene, time, genotype) %>% summarise(mean=mean(rlog_counts), median= median(rlog_counts), SD=sd(rlog_counts), n=n(), SE=SD/sqrt(n)) 	



pdf("output/deseq2/deseq2_TC_Top4genes_cluster4_rlog.pdf", width=11, height=6)  # !!! change cluster nb !!!
stat_significant_deseq2TC_genes %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color=genotype), size=0.75) +
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width=.2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(~gene, nrow = 1, scale = "free")  +	
  xlab(label = "deseq2 normalized counts") +
  ggtitle("Top 4 significant genes in cluster4_rlog") +   # !!! change cluster nb !!!
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()



# Display genes expression profile for all cluster
## rlog norm counts
## Create a function to generate data for the top 2 significant genes for a given cluster
get_top2_genes_data <- function(cluster_number) {
  significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
    inner_join(cluster_gene) %>%
    filter(cluster == cluster_number) %>%
    arrange(padj) %>%
    select(gene, padj, cluster) %>%
    left_join(rlog_counts_tidy)
  
  significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
    select(gene) %>%
    unique() %>%
    slice_head(n = 2)
  
  significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
    inner_join(significant_deseq2TC_genes_cluster_genes)
  
  stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
    select(-replicate) %>%
    group_by(gene, time, genotype, cluster) %>% summarise(mean = mean(rlog_counts), median = median(rlog_counts), SD = sd(rlog_counts), n = n(), SE = SD / sqrt(n))
  
  return(stat_significant_deseq2TC_genes)
}

# Generate data for all 25 clusters
all_clusters_data <- map_dfr(1:25, get_top2_genes_data)

# Plot top 2 significant genes for each of the 25 clusters
pdf("output/deseq2/deseq2_TC_Top2genes_AllClusters_rlog.pdf", width = 30, height = 14)
all_clusters_data %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color = genotype), size = 0.75) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(cluster ~ gene, nrow = 3, labeller = labeller(gene = function(x) paste("Cluster", unlist(all_clusters_data[1, "cluster"]), x)), scales = "free") +
  xlab(label = "rlog counts") +
  ggtitle("Top 2 significant genes in each cluster") +
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()





## repersent tpm counts 
### Load tpm
tpm_all <- read_csv("output/tpm/tpm_all.txt") %>% 
  select(-1) #To import


## Transform tpm tibble into matrix
tpm_all_matrix = make_matrix(select(tpm_all, -Geneid), pull(tpm_all, Geneid)) 

tpm_all_raw_tidy <- as_tibble(tpm_all_matrix, rownames = "gene") %>%
  gather(key = "sample", value = "tpm", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")



tpm_all_raw_tidy$time <-
  factor(tpm_all_raw_tidy$time,
         c("ESC", "NPC", "2dN", "8wN"))
tpm_all_raw_tidy$genotype <-
  factor(tpm_all_raw_tidy$genotype,
         c("WT", "KO", "HET"))

## Create a function to generate data for the top 2 significant genes for a given cluster
get_top2_genes_data <- function(cluster_number) {
  significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
    inner_join(cluster_gene) %>%
    filter(cluster == cluster_number) %>%
    arrange(padj) %>%
    select(gene, padj, cluster) %>%
    left_join(tpm_all_raw_tidy) # here tpm
  
  significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
    select(gene) %>%
    unique() %>%
    slice_head(n = 2)
  
  significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
    inner_join(significant_deseq2TC_genes_cluster_genes)
  
  stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
    select(-replicate) %>%
    group_by(gene, time, genotype, cluster) %>% summarise(mean = mean(tpm), median = median(tpm), SD = sd(tpm), n = n(), SE = SD / sqrt(n))
  
  return(stat_significant_deseq2TC_genes)
}

# Generate data for all 25 clusters
all_clusters_data <- map_dfr(1:25, get_top2_genes_data)

# Plot top 2 significant genes for each of the 25 clusters
pdf("output/deseq2/deseq2_TC_Top2genes_AllClusters_tpm_rlog.pdf", width = 30, height = 14)
all_clusters_data %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color = genotype), size = 0.75) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(cluster ~ gene, nrow = 3, labeller = labeller(gene = function(x) paste("Cluster", unlist(all_clusters_data[1, "cluster"]), x)), scales = "free") +
  xlab(label = "TPM") +
  ggtitle("Top 2 significant genes in each cluster") +
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()
```





### Clustering of the significant genes in Time-Course analyses across genotypes using deseq2 rlog norm-counts
```bash
module load R/4.2.2
```
The below code is clean to generate clustering plot (qvalue and number of cluster can be adjusted), vizualize them and look at invividual gene expressin (rlog counts and tpm)
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")

# Import files to generete the DESeq2Dataset
## samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")
## Import counts_all and transform to matrix
counts_all <- read_csv("output/deseq2/counts_all.txt") %>% 
  select(-1)

### Transform merged_data into a matrix
#### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
#### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet Time-Course 
### desgin = full-model
ddsTC <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                                colData = coldata,
                                design = ~ genotype + time + genotype:time)

### Define the reduced model (not including interaction term for comparison with full model)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ genotype + time)
resTC <- results(ddsTC)


# Normalize the counts
## Normalized counts with deseq2 and tidy it
normalized_counts <- as_tibble(counts(ddsTC, normalized = TRUE), rownames = "gene") %>%
  gather(key = "sample", value = "norm_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

normalized_counts$time <-
  factor(normalized_counts$time,
         c("ESC", "NPC", "2dN", "8wN"))
normalized_counts$genotype <-
  factor(normalized_counts$genotype,
         c("WT", "KO", "HET"))

# Clustering

rlog_counts <- rlog(ddsTC) # take 5 min

## Extract the transformed counts
rlog_counts_matrix <- assay(rlog_counts) 

## Isolate the significant genes by gene id
signif_TC_genes = as_tibble(resTC, rownames = "gene") %>%
  filter(padj<= 0.05) %>%               # !!! Here change qvalue accordingly !!!
  select(gene) %>%
  unique() 

## Convert into vector 
signif_TC_genes_vector <- signif_TC_genes$gene

## Filter our matrix with the significant genes
rlog_counts_matrix_sig <- rlog_counts_matrix[rownames(rlog_counts_matrix) %in% signif_TC_genes_vector, ]
 # double-check the nb of genes is same s in signif_TC_genes
nrow(rlog_counts_matrix_sig)

## Reorder columns
ordered_columns <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3", "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3", "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3", "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3", "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3", "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3", "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4")

rlog_counts_matrix_sig_ordered <- rlog_counts_matrix_sig[, ordered_columns]


# Make a clean table with significant deseq2-TC genes
rlog_counts_tidy <- as_tibble(rlog_counts_matrix_sig_ordered, rownames = "gene") %>%
  gather(key = "sample", value = "rlog_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(rlog_counts_matrix_sig_ordered, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into k clusters
row_clusters <- cutree(row_hclust, k = 25)                   # !!! Here change tree nb accordingly !!!

## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(rlog_counts_matrix_sig_ordered),
                           cluster = row_clusters)


### Save dataframe
write.csv(cluster_gene, file="output/deseq2/cluster_gene_rlog_25cl.txt")


## Compil both
rlog_counts_tidy <- rlog_counts_tidy %>%
  left_join(cluster_gene, by = "gene")

rlog_counts_tidy$time <-
  factor(rlog_counts_tidy$time,
         c("ESC", "NPC", "2dN", "8wN"))
rlog_counts_tidy$genotype <-
  factor(rlog_counts_tidy$genotype,
         c("WT", "KO", "HET"))



# Plot rlog_transform norm deseq2 count with 'loess' method pretty
## Calculate the number of genes per cluster
genes_per_cluster <- rlog_counts_tidy %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2/line_rlog_p0.05_cl25_pretty_span07.pdf", width=20, height=14)      # !!! Here change title accordingly !!!
ggplot(rlog_counts_tidy, aes(x = time, y = rlog_counts)) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.7) +
  # Add mean value for each genotype at each time point
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  # Add standard error bars around the mean points
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  # Add number of genes per cluster to the facet_wrap panels
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()

# Plot rlog_transform norm deseq2 count with 'loess' method pretty with a light transparent color
# Calculate the mean expression for each gene at each time point
mean_rlog_counts_tidy <- rlog_counts_tidy %>%
  group_by(gene, time, genotype, cluster) %>%
  summarise(mean_rlog_counts = mean(rlog_counts)) %>%
  ungroup()

genes_per_cluster <- rlog_counts_tidy %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2/line_rlog_p0.05_cl25_pretty_span08_genes_bg.pdf", width=20, height=14)
ggplot(mean_rlog_counts_tidy, aes(x = time, y = mean_rlog_counts)) +
  geom_line(aes(color = genotype, group = interaction(gene, genotype)), alpha = 0.1) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()


pdf("output/deseq2/line_rlog_p0.05_cl25_pretty_span08_genes_WT-KO.pdf", width=20, height=14)
mean_rlog_counts_tidy %>%
  filter(genotype %in% c("WT","KO")) %>%
ggplot(., aes(x = time, y = mean_rlog_counts)) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()

# Plot TPM count with 'loess' method pretty

# Calculate the mean expression for each gene at each time point
mean_tpm_all_raw_tidy <- tpm_all_raw_tidy %>%
  inner_join(cluster_gene, by = "gene") %>%
  group_by(gene, time, genotype, cluster) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup()

genes_per_cluster <- tpm_all_raw_tidy %>%
  inner_join(cluster_gene, by = "gene") %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2/line_rlog_p0.05_cl25_pretty_span08_genes_bg_tpm.pdf", width=20, height=14)
ggplot(mean_tpm_all_raw_tidy, aes(x = time, y = mean_tpm)) +
  geom_line(aes(color = genotype, group = interaction(gene, genotype)), alpha = 0.1) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()



# Display some genes expression profile cluster-per-cluster

## Gather the 10 first significant deseq2-TC genes of specified cluster
significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
  inner_join(cluster_gene) %>%
  filter(cluster == 4) %>% # !!! change cluster nb !!!
  arrange(padj) %>%
  select(gene,padj) %>%
  left_join(rlog_counts_tidy) 

significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
select(gene) %>%
unique() %>%
slice_head(n=4)
  
significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
  inner_join(significant_deseq2TC_genes_cluster_genes)


## Stat
stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
  select(-replicate) %>%
  group_by(gene, time, genotype) %>% summarise(mean=mean(rlog_counts), median= median(rlog_counts), SD=sd(rlog_counts), n=n(), SE=SD/sqrt(n)) 	



pdf("output/deseq2/deseq2_TC_Top4genes_cluster4_rlog.pdf", width=11, height=6)  # !!! change cluster nb !!!
stat_significant_deseq2TC_genes %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color=genotype), size=0.75) +
  geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE), width=.2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(~gene, nrow = 1, scale = "free")  +	
  xlab(label = "deseq2 normalized counts") +
  ggtitle("Top 4 significant genes in cluster4_rlog") +   # !!! change cluster nb !!!
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()



# Display genes expression profile for all cluster
## rlog norm counts
## Create a function to generate data for the top 2 significant genes for a given cluster
get_top2_genes_data <- function(cluster_number) {
  significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
    inner_join(cluster_gene) %>%
    filter(cluster == cluster_number) %>%
    arrange(padj) %>%
    select(gene, padj, cluster) %>%
    left_join(rlog_counts_tidy)
  
  significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
    select(gene) %>%
    unique() %>%
    slice_head(n = 2)
  
  significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
    inner_join(significant_deseq2TC_genes_cluster_genes)
  
  stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
    select(-replicate) %>%
    group_by(gene, time, genotype, cluster) %>% summarise(mean = mean(rlog_counts), median = median(rlog_counts), SD = sd(rlog_counts), n = n(), SE = SD / sqrt(n))
  
  return(stat_significant_deseq2TC_genes)
}

# Generate data for all 25 clusters
all_clusters_data <- map_dfr(1:25, get_top2_genes_data)

# Plot top 2 significant genes for each of the 25 clusters
pdf("output/deseq2/deseq2_TC_Top2genes_AllClusters_rlog.pdf", width = 30, height = 14)
all_clusters_data %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color = genotype), size = 0.75) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(cluster ~ gene, nrow = 3, labeller = labeller(gene = function(x) paste("Cluster", unlist(all_clusters_data[1, "cluster"]), x)), scales = "free") +
  xlab(label = "rlog counts") +
  ggtitle("Top 2 significant genes in each cluster") +
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()





## repersent tpm counts 
### Load tpm
tpm_all <- read_csv("output/tpm/tpm_all.txt") %>% 
  select(-1) #To import


## Transform tpm tibble into matrix
tpm_all_matrix = make_matrix(select(tpm_all, -Geneid), pull(tpm_all, Geneid)) 

tpm_all_raw_tidy <- as_tibble(tpm_all_matrix, rownames = "gene") %>%
  gather(key = "sample", value = "tpm", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")



tpm_all_raw_tidy$time <-
  factor(tpm_all_raw_tidy$time,
         c("ESC", "NPC", "2dN", "8wN"))
tpm_all_raw_tidy$genotype <-
  factor(tpm_all_raw_tidy$genotype,
         c("WT", "KO", "HET"))

## Create a function to generate data for the top 2 significant genes for a given cluster
get_top2_genes_data <- function(cluster_number) {
  significant_deseq2TC_genes_cluster_all <- as_tibble(resTC, rownames = "gene") %>%
    inner_join(cluster_gene) %>%
    filter(cluster == cluster_number) %>%
    arrange(padj) %>%
    select(gene, padj, cluster) %>%
    left_join(tpm_all_raw_tidy) # here tpm
  
  significant_deseq2TC_genes_cluster_genes <- significant_deseq2TC_genes_cluster_all %>%
    select(gene) %>%
    unique() %>%
    slice_head(n = 2)
  
  significant_deseq2TC_genes_cluster_all_genes <- significant_deseq2TC_genes_cluster_all %>%
    inner_join(significant_deseq2TC_genes_cluster_genes)
  
  stat_significant_deseq2TC_genes <- significant_deseq2TC_genes_cluster_all_genes %>%
    select(-replicate) %>%
    group_by(gene, time, genotype, cluster) %>% summarise(mean = mean(tpm), median = median(tpm), SD = sd(tpm), n = n(), SE = SD / sqrt(n))
  
  return(stat_significant_deseq2TC_genes)
}

# Generate data for all 25 clusters
all_clusters_data <- map_dfr(1:25, get_top2_genes_data)

# Plot top 2 significant genes for each of the 25 clusters
pdf("output/deseq2/deseq2_TC_Top2genes_AllClusters_tpm_rlog.pdf", width = 30, height = 14)
all_clusters_data %>%
  ggplot(., aes(x = time, y = mean, group = genotype)) +
  geom_line(aes(color = genotype), size = 0.75) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
  geom_point(aes(y = mean), size = .75, shape = 15) +
  theme_bw() +
  facet_wrap(cluster ~ gene, nrow = 3, labeller = labeller(gene = function(x) paste("Cluster", unlist(all_clusters_data[1, "cluster"]), x)), scales = "free") +
  xlab(label = "TPM") +
  ggtitle("Top 2 significant genes in each cluster") +
  scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
dev.off()
```


## Check expression of genes involved in neuronal functionality

Genes collected from this [preprint](https://www.biorxiv.org/content/10.1101/2022.06.02.490114v1.full).

```R
# Package
library(tidyverse)
library(rtracklayer)
library(readxl)
library(ggpubr)

## Create table with gene ID and gene name
### Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)

### Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

### Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()
  

## import neuron-related genes 
neurons_genes <- read_excel("output/temp/GeneList_Neural_Function.xlsx") %>%
  inner_join(gene_id_name)


## Plot the neuron-related genes over the time course
### rlog counts


neurons_genes_stat_rlog = rlog_counts_tidy %>%
  dplyr::rename(gene_id = gene) %>%
  inner_join(neurons_genes) %>%
  group_by(gene_id, time, genotype, gene_name, `function`) %>%
  summarise(mean = mean(rlog_counts), median = median(rlog_counts), SD = sd(rlog_counts), n = n(), SE = SD / sqrt(n))
  
neurons_genes_stat_rlog$time <-
  factor(neurons_genes_stat_rlog$time,
         c("ESC", "NPC", "2dN", "8wN"))


  
#### Plot per gene functions
plot_nav_channels = neurons_genes_stat_rlog %>%
  filter(`function` == "Nav Channels") %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color = genotype), size = 0.75) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
    geom_point(aes(y = mean), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, scales = "free", nrow = 1) +
    ylab(label = "rlog_counts") +
    ggtitle("Nav Channels") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
plot_kv_channels = neurons_genes_stat_rlog %>%
  filter(`function` == "Kv Channels") %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color = genotype), size = 0.75) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
    geom_point(aes(y = mean), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, scales = "free", nrow = 1) +
    ylab(label = "rlog_counts") +
    ggtitle("Kv Channels") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
plot_cav_channels = neurons_genes_stat_rlog %>%
  filter(`function` == "Cav Channels") %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color = genotype), size = 0.75) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
    geom_point(aes(y = mean), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, scales = "free", nrow = 1) +
    ylab(label = "rlog_counts") +
    ggtitle("Cav Channels") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
plot_kcl_transporters = neurons_genes_stat_rlog %>%
  filter(`function` == "K+/Cl- transporters") %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color = genotype), size = 0.75) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
    geom_point(aes(y = mean), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, scales = "free", nrow = 1) +
    ylab(label = "rlog_counts") +
    ggtitle("K+/Cl- transporters") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
plot_nak_atpases = neurons_genes_stat_rlog %>%
  filter(`function` == "Na+/K+ ATPases") %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color = genotype), size = 0.75) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
    geom_point(aes(y = mean), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, scales = "free", nrow = 1) +
    ylab(label = "rlog_counts") +
    ggtitle("Na+/K+ ATPases") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
plot_camks_snares_pre_synaptic = neurons_genes_stat_rlog %>%
  filter(`function` == "CaMKs SNAREs and pre-synaptic") %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color = genotype), size = 0.75) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
    geom_point(aes(y = mean), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, scales = "free", nrow = 2) +
    ylab(label = "rlog_counts") +
    ggtitle("CaMKs SNAREs and pre-synaptic") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
plot_post_synaptic = neurons_genes_stat_rlog %>%
  filter(`function` == "Post-synaptic") %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color = genotype), size = 0.75) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
    geom_point(aes(y = mean), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, scales = "free", nrow = 2) +
    ylab(label = "rlog_counts") +
    ggtitle("Post-synaptic") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))
plot_receptors_neurotransmitters = neurons_genes_stat_rlog %>%
  filter(`function` == "Receptors for Neurotransmitters") %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color = genotype), size = 0.75) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = .2) +
    geom_point(aes(y = mean), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, scales = "free", nrow = 2) +
    ylab(label = "rlog_counts") +
    ggtitle("Receptors for Neurotransmitters") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green"))


all_plots <- ggarrange(plot_nav_channels, plot_kv_channels, plot_cav_channels,
                       plot_kcl_transporters, plot_nak_atpases, plot_camks_snares_pre_synaptic,
                       plot_post_synaptic, plot_receptors_neurotransmitters,
                       ncol = 1, nrow = 8)

pdf("output/temp/neurons_genes.pdf", width = 11, height = 40)
print(all_plots)
dev.off()

```
Generate heatmap with vst counts:
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")

# Import files to generete the DESeq2Dataset
## samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", 
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3")
## Import counts_all and transform to matrix
counts_all <- read_csv("output/deseq2/counts_all.txt") %>% 
  select(-1) %>%
  select(Geneid, samples)

### Transform merged_data into a matrix
#### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
#### execute function
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate) %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet Time-Course 
### desgin = full-model
ddsTC <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                                colData = coldata,
                                design = ~ genotype + time + genotype:time)

### Define the reduced model (not including interaction term for comparison with full model)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ genotype + time)

### XXX --> Filter out the genes!!!


#### Normalize the counts and keep the WT one only
normalized_counts <- as_tibble(counts(ddsTC, normalized = TRUE), rownames = "gene") %>%
  inner_join(gene_id_name) %>%
  gather(key = "sample", value = "norm_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

normalized_counts$time <-
  factor(normalized_counts$time,
         c("ESC", "NPC", "2dN", "8wN"))


# Clustering
## transform the norm counts
vst_counts <- vst(ddsTC)


### vst
sampleDists <- dist(t(assay(vst_counts)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$time, rld$genotype, rld$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("output/temp/heatmap_cluster_vst_neurons_genes.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
```

# Check Maturation state of our samples

## Maturity state
Here let's figure out the **maturity state of all our samples**, using 'maturity marker genes' from [EZH1 paper](https://www.medrxiv.org/content/10.1101/2022.08.09.22278430v1.full-text):

### GSEA-related method from the paper
- Perform DEGs from KO vs WT and HET vs WT (padj < 0.05) = *input gene list*
- Extract ranked *cell-type-specific gene lists* from the literature (Uzquiano and Kedaigle et al.), containing the top 500 significantly enriched genes in aRG, CFuPN, and CPN scRNA clusters. These gene lists represent NSC, early-born neuron, and late-born neuron gene sets, respectively.
- Perform GSEA using ClusterProfiler with the prepared input lists and the cell-type-specific gene sets.
- Generate core enrichment lists based on the GSEA results.
- Create heatmaps to visualize the expression patterns of the core enrichment lists using baseR.

Generate input gene list manually from the one-by-one analyses --> Generate input gene list in Google Drive `output/deseq2/DEGs_gene_list.xlsx` --> cp to `/output/age`

```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")
library("readxl")
library(rtracklayer)
library(ggpubr)
library(dendextend)


# Create table with gene ID and gene name
## Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()
  
# import input gene list and marker genes
neurons_genes_maturation <- read_excel("output/age/GeneList_Neural_MaturationNeurons_FromGraciaDiazPreprint.xlsx") %>%
  inner_join(gene_id_name) %>%
  rename(Geneid = gene_id)

input_genes_neurons_genes_maturation <- read_excel("output/age/DEGs_gene_list.xlsx") %>%
  inner_join(neurons_genes_maturation %>% filter(type == "aRG_up"))      ## HERE CHANGE (aRG_up, CFuPN_up, CPN_Up) !!!!!!!!!!!!!!!!!!!!

# Some stat to see how many genes per category:

input_genes_neurons_genes_maturation %>%
  select(Geneid,type) %>%
  unique() %>%
  group_by(type) %>%
  summarise(n=n())


# Show TPM
input_genes_neurons_genes_maturation



### Load tpm
tpm_all_input_genes_neurons_genes_maturation <- read_csv("output/tpm/tpm_all_sample.txt") %>% 
  select(-1) %>%
  inner_join(input_genes_neurons_genes_maturation) %>%
  select(-Geneid, -DEGs,-time,-type,-comparison)

                            


## Transform tpm tibble into matrix
tpm_all_input_genes_neurons_genes_maturation_matrix = make_matrix(select(tpm_all_input_genes_neurons_genes_maturation, -gene_name), pull(tpm_all_input_genes_neurons_genes_maturation, gene_name)) 

# LOG2 !!!!!!! To change !!!!!!!!
tpm_all_input_genes_neurons_genes_maturation_matrix = log2(tpm_all_input_genes_neurons_genes_maturation_matrix + 1) 


# Reorder
col_order <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3", "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3", "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3", "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3", "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3", "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3", "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4", "4wN_iPSCWT_R1", "4wN_iPSCWT_R2", "4wN_iPSCpatient_R1", "4wN_iPSCpatient_R2", "8wN_iPSCpatient_R1", "8wN_iPSCpatient_R2")

  
tpm_all_input_genes_neurons_genes_maturation_matrix_ordered <- tpm_all_input_genes_neurons_genes_maturation_matrix[, col_order]


# Raw heatmap                    ## HERE CHANGE TITLE !!!!!!!!!!!!!!!!!!!!
pdf("output/age/heatmap_GeneList_Neural_MaturationNeurons_aRG_up_DEGsOneByOne_raw_log2tpm.pdf", width=5, height=6)
pheatmap(tpm_all_input_genes_neurons_genes_maturation_matrix_ordered, cluster_rows=F, cluster_cols=F, color= colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
```

--> It is a mess, too many genes are displayed, let's instead of selecting the one-by-one DEGs, selecte the TC-DEgs


```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")
library("readxl")
library(rtracklayer)
library(ggpubr)
library(dendextend)


# Create table with gene ID and gene name
## Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()
  
# import input gene list and marker genes
neurons_genes_maturation <- read_excel("output/age/GeneList_Neural_MaturationNeurons_FromGraciaDiazPreprint.xlsx") %>%
  inner_join(gene_id_name) %>%
  rename(Geneid = gene_id)


input_TC_genes_neurons_genes_maturation <- read_csv("output/deseq2/resTC.txt") %>%
  select(-"...1") %>%
  filter(padj <= 0.05, log2FoldChange > 1 | log2FoldChange < -1 ) %>%
  rename(Geneid = gene) %>%
  inner_join(neurons_genes_maturation %>% filter(type == "aRG_up"))      ## HERE CHANGE (aRG_up, CFuPN_up, CPN_Up) !!!!!!!!!!!!!!!!!!!!

# Some stat to see how many genes per category:

input_TC_genes_neurons_genes_maturation %>%
  select(Geneid,type) %>%
  unique() %>%
  group_by(type) %>%
  summarise(n=n())


# Show TPM

### Load tpm
tpm_all_input_genes_neurons_genes_maturation <- read_csv("output/tpm/tpm_all_sample.txt") %>% 
  select(-1) %>%
  inner_join(input_TC_genes_neurons_genes_maturation) %>%
  select(-Geneid, -baseMean, -log2FoldChange, -lfcSE, -stat, -pvalue,-padj,-type)

                            


## Transform tpm tibble into matrix
tpm_all_input_genes_neurons_genes_maturation_matrix = make_matrix(select(tpm_all_input_genes_neurons_genes_maturation, -gene_name), pull(tpm_all_input_genes_neurons_genes_maturation, gene_name)) 

# LOG2 !!!!!!! To change !!!!!!!!
tpm_all_input_genes_neurons_genes_maturation_matrix = log2(tpm_all_input_genes_neurons_genes_maturation_matrix + 1) 


# Reorder
col_order <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3", "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3", "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3", "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3", "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3", "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3", "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4", "4wN_iPSCWT_R1", "4wN_iPSCWT_R2", "4wN_iPSCpatient_R1", "4wN_iPSCpatient_R2", "8wN_iPSCpatient_R1", "8wN_iPSCpatient_R2")

  
tpm_all_input_genes_neurons_genes_maturation_matrix_ordered <- tpm_all_input_genes_neurons_genes_maturation_matrix[, col_order]


# Raw heatmap                    ## HERE CHANGE TITLE !!!!!!!!!!!!!!!!!!!!
pdf("output/age/heatmap_GeneList_Neural_MaturationNeurons_aRG_up_DEGsTC_raw_log2tpm_p0.05_log2FC1.pdf", width=5, height=6)
pheatmap(tpm_all_input_genes_neurons_genes_maturation_matrix_ordered, cluster_rows=F, cluster_cols=F, color= colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
```

--> Unclear what is our goal here

Here is alternative approach:
- Collect DEGs WT vs HET Rep 3 and 4 at 4wN (up-regulated only)
- Isolate from these genes the cell-type-specific gene lists (top 500 Up; aRG, CFuPN, and CPN)
- Heatmap to look at their expression in the WT and HET at 4week neurons?


```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")
library("readxl")
library(rtracklayer)
library(ggpubr)
library(dendextend)


# Create table with gene ID and gene name
## Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()
  
# import input gene list and marker genes
neurons_genes_maturation <- read_excel("output/age/GeneList_Neural_MaturationNeurons_FromGraciaDiazPreprint.xlsx") %>%
  inner_join(gene_id_name) %>%
  rename(Geneid = gene_id)

input_genes_neurons_genes_maturation <- read_excel("output/age/DEGs_gene_list.xlsx") %>%
  inner_join(neurons_genes_maturation %>% filter(type == "aRG_up"))   %>%   ## HERE CHANGE (aRG_up, CFuPN_up, CPN_Up) !!!!!!!!!!!!!!!!!!!!
  filter(time == "4wN", comparison == "HETr3r4vsWT")

# Some stat to see how many genes per category:

input_genes_neurons_genes_maturation %>%
  select(Geneid,type) %>%
  unique() %>%
  group_by(type) %>%
  summarise(n=n())


# Show TPM
input_genes_neurons_genes_maturation


### Load tpm
tpm_all_input_genes_neurons_genes_maturation <- read_csv("output/tpm/tpm_all_sample.txt") %>% 
  select(-1) %>%
  inner_join(input_genes_neurons_genes_maturation) %>%
  select(gene_name, "4wN_WT_R1", "4wN_WT_R2", "4wN_HET_R1", "4wN_HET_R2", "4wN_HET_R3", "4wN_HET_R4")

make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
       
## Transform tpm tibble into matrix
tpm_all_input_genes_neurons_genes_maturation_matrix = make_matrix(select(tpm_all_input_genes_neurons_genes_maturation, -gene_name), pull(tpm_all_input_genes_neurons_genes_maturation, gene_name)) 

# LOG2 !!!!!!! To change !!!!!!!!
tpm_all_input_genes_neurons_genes_maturation_matrix = log2(tpm_all_input_genes_neurons_genes_maturation_matrix + 1) 


# Reorder
col_order <- c("4wN_WT_R1", "4wN_WT_R2", "4wN_HET_R1", "4wN_HET_R2", "4wN_HET_R3", "4wN_HET_R4")

  
tpm_all_input_genes_neurons_genes_maturation_matrix_ordered <- tpm_all_input_genes_neurons_genes_maturation_matrix[, col_order]


# Raw heatmap                    ## HERE CHANGE TITLE !!!!!!!!!!!!!!!!!!!!
pdf("output/age/heatmap_GeneList_Neural_MaturationNeurons_aRG_up_DEGsHET_log2tpm.pdf", width=5, height=6)
pheatmap(tpm_all_input_genes_neurons_genes_maturation_matrix_ordered, cluster_rows=F, cluster_cols=F, color= colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
```

--> HET R1 and R2 4wN shows same tpm as WT; however Rep3 and Rep4 HET are always higher. 


### Alternative method using genes from other preprint
- Collect the 67 *cell-type-specific gene lists* as input gene list = cluster of our tree (3 cluster)
- Look at vst-counts for these genes for all time-course point and per genotype
- Represent data with heatmap (3 clusters = marker genes and our sampling (ESC, NPC, 2dN,...) at the top); 1 heatmap per genotypes
  - Do not pull replicate so that we can observed individual samples and maybe identify issue

Do this in `/output/age` (transfer gene input lists)

```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")
library("readxl")
library(rtracklayer)
library(ggpubr)
library(dendextend)


# Create table with gene ID and gene name
## Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()
  
# import input gene list and add gene name # !!!! CHOSE GENE LIST HERE !!!!!
neurons_genes_maturation <- read_excel("output/age/GeneList_Neural_MaturationNeurons_FromGraciaDiazPreprint.xlsx") %>%
  inner_join(gene_id_name) %>%
  rename(Geneid = gene_id)

neurons_genes_maturation <- read_excel("output/age/GeneList_Neural_Function_FromCiceripreprint.xlsx") %>%
  inner_join(gene_id_name) %>%
  rename(Geneid = gene_id)

# First WT sample
## samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", 
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3")

## Import counts_all and join with the list of marker genes
counts_all <- read_csv("output/deseq2/counts_all.txt") %>% 
  select(-1) %>%
  select(Geneid, samples) %>%
  inner_join(neurons_genes_maturation) %>%
  select(gene_name, everything(), -Geneid, -type) # put gene_name first column %>%


### Transform merged_data into a matrix
#### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
#### execute function
counts_all_matrix = make_matrix(select(counts_all, -gene_name), pull(counts_all, gene_name)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  select(-replicate, -genotype) %>%
  bind_cols(data.frame(samples))

### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_")%>%
  select(-genotype) %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))


## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet Time-Course for WT only

ddsWT <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                                colData = coldata,
                                design = ~ time)


# Run DESeq with default settings (Wald test)
ddsWT <- DESeq(ddsWT)  # we do not use LRT model as for time course


# vst normalization of our counts !!!! CHOSE NORMALIZATION HERE !!!!!
ddsWT_vst <- vst(ddsWT)
ddsWT_vst <- rlog(ddsWT) #

# make into matrix
ddsWT_vst_matrix = make_matrix(data.frame(assay(ddsWT_vst)))
ddsWT_vst_df = data.frame(assay(ddsWT_vst))

# Reorder
col_order <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
  "X2dN_WT_R1", "X2dN_WT_R2", "X2dN_WT_R3", "X8wN_WT_R1", "X8wN_WT_R2", "X8wN_WT_R3", "X8wN_WT_R4")
  
ddsWT_vst_matrix_ordered <- ddsWT_vst_matrix[, col_order]


# Raw heatmap !!!! CHOSE GENE LIST HERE !!!!!
pdf("output/age/heatmap_GeneList_Neural_MaturationNeurons_WT_raw.pdf", width=5, height=6)
pdf("output/age/heatmap_GeneList_Neural_Function_WT_raw.pdf", width=5, height=6) # 
pheatmap(ddsWT_vst_matrix_ordered, cluster_rows=F, cluster_cols=F)
dev.off()



# Cluster heatmap
## make the gene type matrix
neurons_genes_maturation_df <- neurons_genes_maturation %>%
  select(gene_name, type) %>%
  as.data.frame()

rownames(neurons_genes_maturation_df) <- neurons_genes_maturation_df$gene_name
neurons_genes_maturation_df$gene_name <- NULL

# Sort the rows of the heatmap based on gene types
sorted_rows <- order(neurons_genes_maturation_df$type)
sorted_ddsWT_vst_matrix <- ddsWT_vst_matrix[sorted_rows, ]
sorted_neuron_genes_maturation_df <- neurons_genes_maturation_df[sorted_rows, ]

# Create the annotation dataframe
# annotation_df <- data.frame(type = sorted_neuron_genes_maturation_df$type) This shit was working at a time...
annotation_df = data.frame(sorted_neuron_genes_maturation_df)



# Reorder
col_order <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
  "X2dN_WT_R1", "X2dN_WT_R2", "X2dN_WT_R3", "X8wN_WT_R1", "X8wN_WT_R2", "X8wN_WT_R3", "X8wN_WT_R4")
  
sorted_ddsWT_vst_matrix_ordered <- sorted_ddsWT_vst_matrix[, col_order]


# Plot the sorted heatmap without dendrogram !!!! CHOSE GENE LIST HERE !!!!!
pdf("output/age/heatmap_GeneList_Neural_MaturationNeurons_WT_cluster.pdf", width=5, height=6)
pdf("output/age/heatmap_GeneList_Neural_Function_WT_cluster.pdf", width=5, height=6)
pheatmap(sorted_ddsWT_vst_matrix_ordered, annotation_row = annotation_df, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


# Lets display TPM
# import input gene list and add gene name # !!!! CHOSE GENE LIST HERE !!!!!
neurons_genes_maturation <- read_excel("output/age/GeneList_Neural_MaturationNeurons_FromGraciaDiazPreprint.xlsx") %>%
  inner_join(gene_id_name) %>%
  rename(Geneid = gene_id)

neurons_genes_maturation <- read_excel("output/age/GeneList_Neural_Function_FromCiceripreprint.xlsx") %>%
  inner_join(gene_id_name) %>%
  rename(Geneid = gene_id)


### Load tpm
tpm_all <- read_csv("output/tpm/tpm_all.txt") %>% 
  select(-1) %>%
  inner_join(neurons_genes_maturation) %>%
  select(gene_name, everything(), -Geneid, -type) %>%
  select(gene_name, "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
  "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4")

## Transform tpm tibble into matrix
tpm_all_matrix = make_matrix(select(tpm_all, -gene_name), pull(tpm_all, gene_name)) 

# Reorder
col_order <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
  "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", )

  
tpm_all_matrix_ordered <- tpm_all_matrix[, col_order]


# Raw heatmap !!!! CHOSE GENE LIST HERE !!!!!
pdf("output/age/heatmap_GeneList_Neural_MaturationNeurons_WT_raw_tpm.pdf", width=5, height=6)
pdf("output/age/heatmap_GeneList_Neural_Function_WT_raw_tpm.pdf", width=5, height=6) # 
pheatmap(tpm_all_matrix_ordered, cluster_rows=F, cluster_cols=F, color= colorRampPalette(c("red", "yellow", "green", "cyan", "blue"))(100))
dev.off()



# Cluster heatmap
## make the gene type matrix
neurons_genes_maturation_df <- neurons_genes_maturation %>%
  select(gene_name, type) %>%
  as.data.frame()

rownames(neurons_genes_maturation_df) <- neurons_genes_maturation_df$gene_name
neurons_genes_maturation_df$gene_name <- NULL

# Sort the rows of the heatmap based on gene types
sorted_rows <- order(neurons_genes_maturation_df$type)
sorted_tpm_all_matrix <- tpm_all_matrix[sorted_rows, ]
sorted_neuron_genes_maturation_df <- neurons_genes_maturation_df[sorted_rows, ]

# Create the annotation dataframe
# annotation_df <- data.frame(type = sorted_neuron_genes_maturation_df$type) This shit was working at a time...
annotation_df = data.frame(sorted_neuron_genes_maturation_df)



# Reorder
col_order <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
  "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4")
  
tpm_all_matrix_ordered <- tpm_all_matrix[, col_order]


# Plot the sorted heatmap without dendrogram !!!! CHOSE GENE LIST HERE !!!!!
pdf("output/age/heatmap_GeneList_Neural_MaturationNeurons_WT_cluster_tpm.pdf", width=5, height=6)
pdf("output/age/heatmap_GeneList_Neural_Function_WT_cluster_tpm.pdf", width=5, height=6)
pheatmap(tpm_all_matrix_ordered, annotation_row = annotation_df, cluster_rows = FALSE, cluster_cols = FALSE, color= colorRampPalette(c("red", "yellow", "green", "cyan", "blue"))(100))
dev.off()
```
*NOTE: I used the wald model here to construct the DEseq2 table, not the LRT nested model as for the TimeCourse experiment as LRT is used to compare two nested models: the full model, which includes the interaction term (genotype:time), and the reduced model, which does not include the interaction term.*

--> Using TPM is working pretty well, notably with the genes identified from Ciceri preprint. Now lets load all the tpm for all samples and see if we observe changes in transcript abundance for key marker maturity genes

```R
# Generate TPM for ALL samples
## collect all samples ID
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
                     "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
                     "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
                     "4wN_WT_R1", "4wN_WT_R2", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4",
                     "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
                     "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3",
                     "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
                     "4wN_HET_R1", "4wN_HET_R2", "4wN_HET_R3", "4wN_HET_R4",
                     "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4",
                     "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
                     "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
                     "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
                     "4wN_KO_R1", "4wN_KO_R2", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4",
                     "4wN_iPSCWT_R1", "4wN_iPSCWT_R2", "4wN_iPSCpatient_R1", "4wN_iPSCpatient_R2", "8wN_iPSCpatient_R1", "8wN_iPSCpatient_R2")

## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/tpm/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR.fastp.")) %>%
    rename(!!sample := starts_with("output.STAR.fastp."))
}

## Merge all dataframe into a single one
tpm_all_sample <- reduce(sample_data, full_join, by = "Geneid")
write.csv(tpm_all_sample, file="output/tpm/tpm_all_sample.txt")
### If need to import: tpm_all_sample <- read_csv("output/tpm/tpm_all_sample.txt") #To import


# import input gene list and add gene name # !!!! CHOSE GENE LIST HERE !!!!!
neurons_genes_maturation <- read_excel("output/age/GeneList_Neural_Function_FromCiceripreprint.xlsx") %>%
  inner_join(gene_id_name) %>%
  rename(Geneid = gene_id)


tpm_all_sample_neurons = tpm_all_sample %>%
  inner_join(neurons_genes_maturation) %>%
  select(gene_name, everything(), -Geneid, -type)


## Transform tpm tibble into matrix
tpm_all_sample_neurons_matrix = make_matrix(select(tpm_all_sample_neurons, -gene_name), pull(tpm_all_sample_neurons, gene_name)) 


# Raw heatmap !!!! CHOSE GENE LIST HERE !!!!!
pdf("output/age/heatmap_GeneList_Neural_Function_WT_raw_tpm_all.pdf", width=10, height=8) # 
pheatmap(tpm_all_sample_neurons_matrix, cluster_rows=F, cluster_cols=F, color= colorRampPalette(c("red", "yellow", "green", "cyan", "blue"))(100))
dev.off()

## Log2
tpm_all_sample_neurons_matrix_log2 = log2(tpm_all_sample_neurons_matrix+1)

# Raw heatmap 
pdf("output/age/heatmap_GeneList_Neural_Function_WT_raw_log2tpm_all.pdf", width=10, height=8) # 
pheatmap(tpm_all_sample_neurons_matrix_log2, cluster_rows=F, cluster_cols=F, color= colorRampPalette(c("blue", "white", "red"))(50))
dev.off()


# Represent with geom_smooth instead
tpm_all_sample_neurons_tidy = as_tibble(tpm_all_sample_neurons_matrix, rownames = "gene_name") %>%
  gather(key = "sample", value = "tpm", -gene_name) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")


tpm_all_sample_neurons_tidy$time <-
  factor(tpm_all_sample_neurons_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
tpm_all_sample_neurons_tidy$genotype <-
  factor(tpm_all_sample_neurons_tidy$genotype,
         c("WT", "KO", "HET", "iPSCWT","iPSCpatient"))


tpm_all_sample_neurons_tidy_stat = tpm_all_sample_neurons_tidy %>%
  group_by(gene_name, time, genotype) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup()


pdf("output/age/line_GeneList_Neural_Function_WT_raw_tpm_all.pdf", width=20, height=14)
ggplot(tpm_all_sample_neurons_tidy_stat, aes(x = time, y = mean_tpm)) +
  geom_line(aes(color = genotype, group = interaction(gene_name, genotype)), alpha = 0.1) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()

pdf("output/age/line_GeneList_Neural_Function_WT_raw_tpm_all.pdf", width=20, height=14)
ggplot(tpm_all_sample_neurons_tidy_stat, aes(x = time, y = mean_tpm)) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()


# LOG

tpm_all_sample_neurons_tidy = as_tibble(tpm_all_sample_neurons_matrix, rownames = "gene_name") %>%
  gather(key = "sample", value = "tpm", -gene_name) %>%
  mutate(tpm = log2(tpm+1)) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")



tpm_all_sample_neurons_tidy$time <-
  factor(tpm_all_sample_neurons_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
tpm_all_sample_neurons_tidy$genotype <-
  factor(tpm_all_sample_neurons_tidy$genotype,
         c("WT", "KO", "HET", "iPSCWT","iPSCpatient"))


tpm_all_sample_neurons_tidy_stat = tpm_all_sample_neurons_tidy %>%
  group_by(gene_name, time, genotype) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup()


pdf("output/age/line_GeneList_Neural_Function_WT_raw_log2tpm_all.pdf", width=20, height=14)
ggplot(tpm_all_sample_neurons_tidy_stat, aes(x = time, y = mean_tpm)) +
  geom_line(aes(color = genotype, group = interaction(gene_name, genotype)), alpha = 0.1) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()

pdf("output/age/line_GeneList_Neural_Function_WT_raw_log2tpm_all.pdf", width=20, height=14)
ggplot(tpm_all_sample_neurons_tidy_stat, aes(x = time, y = mean_tpm)) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()

## Comparison WT vs HET replicates

tpm_all_sample_neurons_tidy = as_tibble(tpm_all_sample_neurons_matrix, rownames = "gene_name") %>%
  gather(key = "sample", value = "tpm", -gene_name) %>%
  mutate(tpm = log2(tpm+1)) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_") %>%
  filter(genotype %in% c("WT", "HET")) %>%
  filter(!(genotype == "HET" & !(replicate %in% c("R1", "R2"))))   # !!!! Change HERE !!!! to keep replicate to show



tpm_all_sample_neurons_tidy$time <-
  factor(tpm_all_sample_neurons_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
tpm_all_sample_neurons_tidy$genotype <-
  factor(tpm_all_sample_neurons_tidy$genotype,
         c("WT", "KO", "HET", "iPSCWT","iPSCpatient"))


tpm_all_sample_neurons_tidy_stat = tpm_all_sample_neurons_tidy %>%
  group_by(gene_name, time, genotype) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup()



pdf("output/age/line_GeneList_Neural_Function_raw_log2tpm_HETR1R2.pdf", width=20, height=14) # !!!! Change HERE NAME !!!!
ggplot(tpm_all_sample_neurons_tidy_stat, aes(x = time, y = mean_tpm)) + 
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()


## Comparison WT vs KO genotype

tpm_all_sample_neurons_tidy = as_tibble(tpm_all_sample_neurons_matrix, rownames = "gene_name") %>%
  gather(key = "sample", value = "tpm", -gene_name) %>%
  mutate(tpm = log2(tpm+1)) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_") %>%
  filter(genotype %in% c("WT", "KO")) 
  

tpm_all_sample_neurons_tidy$time <-
  factor(tpm_all_sample_neurons_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
tpm_all_sample_neurons_tidy$genotype <-
  factor(tpm_all_sample_neurons_tidy$genotype,
         c("WT", "KO", "HET", "iPSCWT","iPSCpatient"))


tpm_all_sample_neurons_tidy_stat = tpm_all_sample_neurons_tidy %>%
  group_by(gene_name, time, genotype) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup()



pdf("output/age/line_GeneList_Neural_Function_raw_log2tpm_KO.pdf", width=20, height=14) # !!!! Change HERE NAME !!!!
ggplot(tpm_all_sample_neurons_tidy_stat, aes(x = time, y = mean_tpm)) + 
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()
```
--> It seems that replicate 1/2 have similar TPM as WT as 4wN, however the replicate 3/4 have lower TPM. That may explain the DEGs only observed for replicate 3/4

--> We expected the KO to be less mature, and that is the case; however, the HET does not mature faster than WT, **at least when looking at these specific marker genes**


# Clean code to check log2TPM value for heatmap (gene expression)_hg19

*NOTE: Carefull if looking at mutants need filtering! (iPSC suck and 4wN HET R1 and R2 sucks)*
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("factoextra")
library("gridExtra")
library("readxl")
library(rtracklayer)
library(ggpubr)
library(dendextend)


### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}

# Import tpm df

tpm_all_sample <- read_csv("output/tpm/tpm_all_sample.txt") %>%
  select(-"...1")

# Join with gene names for convenience
## Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


# Tidy df
tpm_all_sample_tidy = as_tibble(tpm_all_sample) %>%
  gather(key = "sample", value = "tpm", -Geneid) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")
  

tpm_all_sample_tidy$time <-
  factor(tpm_all_sample_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
tpm_all_sample_tidy$genotype <-
  factor(tpm_all_sample_tidy$genotype,
         c("WT", "KO", "HET", "iPSCWT","iPSCpatient"))


tpm_all_sample_tidy_gene_name = tpm_all_sample_tidy %>%
  dplyr::rename(gene_id=Geneid) %>%
  left_join(gene_id_name) %>%
  unique()

# stat and tidy data for matrix

## PRC2subunits
tpm_all_sample_tidy_gene_name_stat = tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("EZH2", "EZH1", "EED", "SUZ12", "RBBP7", "RBBP4", "AEBP2", "JARID2", "PHF1", "MTF2", "PHF19"),
         genotype == "WT") %>%  
  select(-replicate, -genotype, -gene_id) %>%
  group_by(gene_name, time) %>%
  summarise(mean=mean(tpm)) %>% 
  pivot_wider(names_from = time, values_from = mean) %>%
  ungroup()
## Transform tpm tibble into matrix
tpm_all_sample_tidy_gene_name_stat_matrix = make_matrix(select(tpm_all_sample_tidy_gene_name_stat, -gene_name), pull(tpm_all_sample_tidy_gene_name_stat, gene_name)) 

tpm_all_sample_tidy_gene_name_stat_matrix_log2 = log2(tpm_all_sample_tidy_gene_name_stat_matrix+1)
# Raw heatmap 
pdf("output/deseq2/heatmap_PRC2subunits_log2tpm.pdf", width=5, height=6) # 
pheatmap(tpm_all_sample_tidy_gene_name_stat_matrix_log2, cluster_rows=F, cluster_cols=F, color= colorRampPalette(c("blue", "white", "red"))(50))
dev.off()


## Genes involved in chromatin remodeling and transcription regulation causative for Neurodevelopmental disorders NDD (Gabriele et al 2018)
tpm_all_sample_tidy_gene_name_stat = tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("ADNP", "ARID1A", "ARID1B", "ARX", "ATRX", "AUTS2", "BAZ1B", "BCOR", "BRD1", "BRD2", "BRD3", "BRWD3", "CDKL5", "CECR2", "CHD2", "CHD7", "CHD8", "CREBBP", "CTCF", "DNMT1", "DNMT3A", "DNMT3B", "EBF3", "EED", "EHMT1", "EP300", "EZH2", "KAT2A", "GTF2I", "GTF2IRD1", "HDAC8", "PHC1", "KDM5C", "KANSL1", "KDM5C", "KDM6A", "KMT2A", "KMT2C", "KMT2D", "SETD1A", "MBD5", "MECP2", "MED12", "NBEA", "NSD1", "WHSC1", "PHF21A", "PHF6", "PHF8", "POLR1C", "POLR1D", "RNASEH2C", "RPS6KA3", "RPS6KA6", "SATB2", "SETD2", "SMARCA2", "SMARCA4", "SMARCB1", "SMARCE1", "SOX11", "SOX3", "TAF1", "YY1", "ZBTB20", "ZNF41", "ZNF674", "ZNF711", "ZNF81"),
         genotype == "WT") %>%  
  select(-replicate, -genotype, -gene_id) %>%
  group_by(gene_name, time) %>%
  summarise(mean=mean(tpm)) %>% 
  pivot_wider(names_from = time, values_from = mean) %>%
  ungroup()
## Transform tpm tibble into matrix
tpm_all_sample_tidy_gene_name_stat_matrix = make_matrix(select(tpm_all_sample_tidy_gene_name_stat, -gene_name), pull(tpm_all_sample_tidy_gene_name_stat, gene_name)) 

tpm_all_sample_tidy_gene_name_stat_matrix_log2 = log2(tpm_all_sample_tidy_gene_name_stat_matrix+1)
# Raw heatmap 
pdf("output/deseq2/heatmap_CR_NDD_log2tpm.pdf", width=3, height=8) # 
pheatmap(tpm_all_sample_tidy_gene_name_stat_matrix_log2, cluster_rows=F, cluster_cols=F, color= colorRampPalette(c("blue", "white", "red"))(50))
dev.off()


# As line: FAIL NEED TROUBLSHOOT IF NEEDED
tpm_all_sample_tidy_gene_name_stat_matrix_log2_tidy = as_tibble(tpm_all_sample_tidy_gene_name_stat_matrix_log2, rownames = "gene_name") %>%
  gather(key = "time", value = "log2tpm", -gene_name)

tpm_all_sample_tidy_gene_name_stat_matrix_log2_tidy$time <-
  factor(tpm_all_sample_tidy_gene_name_stat_matrix_log2_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))


pdf("output/age/line_CR_NDD_log2tpm_all.pdf", width=20, height=14)
ggplot(tpm_all_sample_tidy_gene_name_stat_matrix_log2_tidy, aes(x = time, y = log2tpm)) +
  geom_line(aes(color = gene_name), alpha = 0.1) +
  geom_smooth(aes(color = gene_name, group = gene_name), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = gene_name, group = gene_name), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = gene_name, group = gene_name), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()





```


# Clean code to check TPM of individual genes (gene expression)_hg19
```R
# Load packages
library("tidyverse")
library("RColorBrewer")
library(rtracklayer)
library(ggpubr)

# Import tpm df

tpm_all_sample <- read_csv("output/tpm/tpm_all_sample.txt") %>%
  select(-"...1")

# Join with gene names for convenience
## Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


# Tidy df
tpm_all_sample_tidy = as_tibble(tpm_all_sample) %>%
  gather(key = "sample", value = "tpm", -Geneid) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")
  

tpm_all_sample_tidy$time <-
  factor(tpm_all_sample_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
tpm_all_sample_tidy$genotype <-
  factor(tpm_all_sample_tidy$genotype,
         c("WT", "KO", "HET", "iPSCWT","iPSCpatient"))


tpm_all_sample_tidy_gene_name = tpm_all_sample_tidy %>%
  dplyr::rename(gene_id=Geneid) %>%
  left_join(gene_id_name) %>%
  unique()






# Display plot per gene
# Stat
tpm_all_sample_tidy_gene_name_stat <- tpm_all_sample_tidy_gene_name %>%
  select(-replicate) %>%
  group_by(gene_id, gene_name, time, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	



pdf("output/deseq2/genes_EZH.pdf", width=8, height=5)
tpm_all_sample_tidy_gene_name_stat %>%
  filter(gene_name %in% c("EZH1", "EZH2")) %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color=genotype), size=0.75) +
    geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE,color=genotype), width=.2) +
    geom_point(aes(y = mean,color=genotype), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, nrow = 1, scale = "free")  +	
    ylab(label = "tpm") +
    ggtitle("") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()


tpm_all_sample_tidy_gene_name_stat <- tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("JARID2", "ASTN2")) %>%
  select(-replicate) %>%
  group_by(gene_id, gene_name, time, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	
pdf("output/deseq2/genes_JARID2_ASTN2.pdf", width=8, height=5)
tpm_all_sample_tidy_gene_name_stat %>%
  filter(gene_name %in% c("JARID2", "ASTN2"),
         genotype %in% c("WT")) %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color=genotype), size=0.75) +
    geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE,color=genotype), width=.2) +
    geom_point(aes(y = mean,color=genotype), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, nrow = 1, scale = "free")  +	
    ylab(label = "tpm") +
    ggtitle("") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()


pdf("output/deseq2/genes_JARID2_ASTN2_allgenotypes.pdf", width=8, height=5)
tpm_all_sample_tidy_gene_name_stat %>%
  filter(gene_name %in% c("JARID2", "ASTN2"),
         genotype %in% c("WT", "KO", "HET")) %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color=genotype), size=0.75) +
    geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE,color=genotype), width=.2) +
    geom_point(aes(y = mean,color=genotype), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, nrow = 1, scale = "free")  +	
    ylab(label = "tpm") +
    ggtitle("") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()




tpm_all_sample_tidy_gene_name_stat <- tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("EZH2", "EZH1", "EED", "SUZ12", "RBBP7", "RBBP4", "AEBP2", "JARID2", "PHF1", "MTF2", "PHF19")) %>%
  select(-replicate) %>%
  group_by(gene_id, gene_name, time, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	
pdf("output/deseq2/genes_PRC2subunits.pdf", width=15, height=10)
tpm_all_sample_tidy_gene_name_stat %>%
  filter(gene_name %in% c("EZH2", "EZH1", "EED", "SUZ12", "RBBP7", "RBBP4", "AEBP2", "JARID2", "PHF1", "MTF2", "PHF19"),
         genotype %in% c("WT")) %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color=genotype), size=0.75) +
    geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE,color=genotype), width=.2) +
    geom_point(aes(y = mean,color=genotype), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, nrow = 2, scale = "free")  +	
    ylab(label = "tpm") +
    ggtitle("PRC2 subunits") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()



tpm_all_sample_tidy_gene_name_stat <- tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("EZH2", "EZH1", "EED", "SUZ12", "RBBP7", "RBBP4", "AEBP2", "JARID2", "PHF1", "MTF2", "PHF19")) %>%
  select(-replicate) %>%
  group_by(gene_id, gene_name, time, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	
pdf("output/deseq2/genes_PRC2subunits_allGenotypes.pdf", width=15, height=10)
tpm_all_sample_tidy_gene_name_stat %>%
  filter(gene_name %in% c("EZH2", "EZH1", "EED", "SUZ12", "RBBP7", "RBBP4", "AEBP2", "JARID2", "PHF1", "MTF2", "PHF19"),
         genotype %in% c("WT", "KO", "HET")) %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color=genotype), size=0.75) +
    geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE,color=genotype), width=.2) +
    geom_point(aes(y = mean,color=genotype), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, nrow = 2, scale = "free")  +	
    ylab(label = "tpm") +
    ggtitle("PRC2 subunits") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()




tpm_all_sample_tidy_gene_name_stat <- tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("NEUROG1", "NEUROG2")) %>%
  select(-replicate) %>%
  group_by(gene_id, gene_name, time, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	
pdf("output/deseq2/genes_NEUROG.pdf", width=8, height=5)
tpm_all_sample_tidy_gene_name_stat %>%
  filter(gene_name %in% c("NEUROG1", "NEUROG2"),
         genotype %in% c("WT", "KO", "HET")) %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color=genotype), size=0.75) +
    geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE,color=genotype), width=.2) +
    geom_point(aes(y = mean,color=genotype), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, nrow = 1, scale = "free")  +	
    ylab(label = "tpm") +
    ggtitle("") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()



tpm_all_sample_tidy_gene_name_stat <- tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("NEUROG2")) %>%
  select(-replicate) %>%
  group_by(gene_id, gene_name, time, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	
pdf("output/deseq2/gene_NEUROG2_WT-KO.pdf", width=8, height=5)
tpm_all_sample_tidy_gene_name_stat %>%
  filter(gene_name %in% c("NEUROG2"),
         genotype %in% c("WT", "KO")) %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color=genotype), size=0.75) +
    geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE,color=genotype), width=.2) +
    geom_point(aes(y = mean,color=genotype), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, nrow = 1, scale = "free")  +	
    ylab(label = "tpm") +
    ggtitle("") +
    scale_color_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()




# Filtered some replicate

tpm_all_sample_tidy_gene_name_selected = tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("EZH1", "EZH2")) %>%
  mutate(keep = ifelse(genotype == "HET" & (replicate %in% c("R1", "R2")) & time == "4wN", TRUE, genotype != "HET")) %>%
  filter(keep) %>%
  filter(time == "4wN", genotype %in% c("WT","HET")) %>%
  group_by(gene_name, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	
pdf("output/deseq2/genes_EZH_HET4wNRep12.pdf", width=4, height=3)
tpm_all_sample_tidy_gene_name_selected %>%
  ggplot(., aes(x = genotype, y = mean, fill = genotype)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.2) +
    ylab("tpm") +
    facet_wrap(~gene_name, nrow = 1, scale = "free")  +	
    ggtitle("4 week neurons; Het replicate 1 and 2") +
    theme_bw()+
    scale_fill_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()



tpm_all_sample_tidy_gene_name_selected = tpm_all_sample_tidy_gene_name %>%
  filter(gene_name %in% c("EZH1", "EZH2")) %>%
  mutate(keep = ifelse(genotype == "HET" & (replicate %in% c("R3", "R4")) & time == "4wN", TRUE, genotype != "HET")) %>%
  filter(keep) %>%
  filter(time == "4wN", genotype %in% c("WT","HET")) %>%
  group_by(gene_name, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	
pdf("output/deseq2/genes_EZH_HET4wNRep34.pdf", width=4, height=3)
tpm_all_sample_tidy_gene_name_selected %>%
  ggplot(., aes(x = genotype, y = mean, fill = genotype)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.2) +
    ylab("tpm") +
    facet_wrap(~gene_name, nrow = 1, scale = "free")  +	
    ggtitle("4 week neurons; Het replicate 3 and 4") +
    theme_bw()+
    scale_fill_manual(values = c("WT" = "grey", "KO" = "red", "HET" = "green", "iPSCWT" = "black", "iPSCpatient" = "orange"))
dev.off()
```


# Gene ontology
We will use clusterProfile package. Tutorial [here](https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html).

## Light test on 2 clusters 

Let's do a test of the pipeline with genes from cluster4 amd cluster14 from the rlog counts. Our background list will be all genes tested for differential expression.

```R
# packages
library(clusterProfiler)
library(pathview)
library(DOSE)
library(org.Hs.eg.db)
library(enrichplot)
library(rtracklayer)

## Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


# Import genes_cluster list and background list
cluster_4 = read_csv("output/deseq2/cluster_gene_rlog_25cl.txt") %>%
  filter(cluster == 4) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(cluster_4$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Save GO analyses
GO_summary <- data.frame(ego)

write.csv(GO_summary, "output/GO/cluster4_BP.csv")
write.csv(GO_summary, "output/GO/cluster4_MF.csv")
write.csv(GO_summary, "output/GO/cluster4_CC.csv")



# Vizualization

pdf("output/GO/dotplot_BP_cluster_4.pdf", width=8, height=11)
dotplot(ego, showCategory=50)
dev.off()

pdf("output/GO/emapplot_BP_cluster_4.pdf", width=8, height=11)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()



# Vizualization of some additonal clusters (2/4/5/8/9/14/18/21/22):
## Cluster 2
# Import genes_cluster list and background list
cluster = read_csv("output/deseq2/cluster_gene_rlog_25cl.txt") %>%
  filter(cluster == 2) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(cluster$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Vizualization
pdf("output/GO/dotplot_BP_cluster_2.pdf", width=8, height=11)
dotplot(ego, showCategory=50)
dev.off()

pdf("output/GO/emapplot_BP_cluster_2.pdf", width=8, height=11)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()


## Cluster 8
# Import genes_cluster list and background list
cluster = read_csv("output/deseq2/cluster_gene_rlog_25cl.txt") %>%
  filter(cluster == 8) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(cluster$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Vizualization
pdf("output/GO/dotplot_BP_cluster_8.pdf", width=8, height=11)
dotplot(ego, showCategory=50)
dev.off()

pdf("output/GO/emapplot_BP_cluster_8.pdf", width=8, height=11)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()







## Cluster 14
# Import genes_cluster list and background list
cluster = read_csv("output/deseq2/cluster_gene_rlog_25cl.txt") %>%
  filter(cluster == 14) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(cluster$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Vizualization
pdf("output/GO/dotplot_BP_cluster_14.pdf", width=8, height=11)
dotplot(ego, showCategory=50)
dev.off()

pdf("output/GO/emapplot_BP_cluster_14.pdf", width=8, height=11)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()



## Cluster 16
# Import genes_cluster list and background list
cluster = read_csv("output/deseq2/cluster_gene_rlog_25cl.txt") %>%
  filter(cluster == 16) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(cluster$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Vizualization
pdf("output/GO/dotplot_BP_cluster_16.pdf", width=8, height=11)
dotplot(ego, showCategory=50)
dev.off()

pdf("output/GO/emapplot_BP_cluster_16.pdf", width=8, height=11)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()





## Cluster 18
# Import genes_cluster list and background list
cluster = read_csv("output/deseq2/cluster_gene_rlog_25cl.txt") %>%
  filter(cluster == 18) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(cluster$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Vizualization
pdf("output/GO/dotplot_BP_cluster_18.pdf", width=8, height=11)
dotplot(ego, showCategory=50)
dev.off()

pdf("output/GO/emapplot_BP_cluster_18.pdf", width=8, height=11)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()




## Cluster 21
# Import genes_cluster list and background list
cluster = read_csv("output/deseq2/cluster_gene_rlog_25cl.txt") %>%
  filter(cluster == 21) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(cluster$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Vizualization
pdf("output/GO/dotplot_BP_cluster_21.pdf", width=8, height=20)
dotplot(ego, showCategory=50)
dev.off()

pdf("output/GO/emapplot_BP_cluster_21.pdf", width=15, height=20)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()


## Cluster 22
# Import genes_cluster list and background list
cluster = read_csv("output/deseq2/cluster_gene_rlog_25cl.txt") %>%
  filter(cluster == 22) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


# Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(cluster$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Vizualization
pdf("output/GO/dotplot_BP_cluster_22.pdf", width=8, height=20)
dotplot(ego, showCategory=50)
dev.off()
pdf("output/GO/dotplot_BP_cluster_22_short.pdf", width=8, height=10)
dotplot(ego, showCategory=25)
dev.off()

pdf("output/GO/emapplot_BP_cluster_22.pdf", width=15, height=15)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()
```

--> The pipeline works GREAT !!

Now let's re-analyze everything with the more up to date genome


# Re-analyzes using the GRCh38/hg38 genome

Input files are still our fastp-clean-trimmed data

## Mapping with STAR
### Index the genome
hg19 genome with 12CPU and 50G mem (time=<1.5day). So let's go for: 12CPU and 250G mem (time=<3hrs)

```bash
sbatch scripts/STAR_index_hg38.sh # 12323950 ok
```

### Mapp the reads to features (fastp-clean)
Mapping and indexation:
```bash
sbatch scripts/STAR_mapping_hg38_1.sh # 12325921 ok
sbatch scripts/STAR_mapping_hg38_2.sh # 12325922 ok
```

--> Let's compil the number of uniquely mapped reads for all files (add it in the Google Drive `RNAseq_infos.xlsx` file)

```bash
# Print nb of uniq map reads for raw mapping
for file in output/STAR_hg38/*Log.final.out; do
    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $file | awk '{print $NF}')
    echo "$file: Number of uniquely mapped reads: $uniquely_mapped_reads"
done > output/STAR_hg38/uniq_map_reads_counts.txt
```

## Count with featurecounts


Count on gene features with parameter
```bash
conda activate featurecounts

# example for 1 file:
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/featurecounts/${x}.txt output/STAR_hg38/${x}_Aligned.sortedByCoord.out.bam

# all samples:
sbatch scripts/featurecounts_hg38_1.sh # 12342501
sbatch scripts/featurecounts_hg38_1.sh # 12342500
```
# Quality control metrics
Print number of succesfully assigned alignments for each sample (add to drive `RNAseq_infos.xlsx`)
```bash
for file in output/featurecounts_hg38/*.summary; do
    assigned_reads=$(grep "Assigned" $file | awk '{print $NF}')
    echo "$file: Assigned: $assigned_reads"
done > output/featurecounts_hg38/assigned_reads_counts.tsv
```
Print the total number of reads
```bash
for file in output/STAR_hg38/*.final.out; do
    input_reads=$(grep "Number of input reads" $file | awk '{print $NF}')
    echo "$file: Number of input reads: $input_reads"
done > output/STAR_hg38/input_reads_counts.txt
```
Add these values to the `RNAseq_infos.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.



# Generate Bigwig coverage files

Let's generate **TPM coverage**:

```bash
conda activate deeptools
# run time-per-time:
sbatch scripts/TPM_bw_hg38_1.sh # 12377715 ok
sbatch scripts/TPM_bw_hg38_2.sh # 12377716 ok
```


Let's merge the bigwig into 1 file with wiggletools (will do average of bigwig signal and not sum, many options see [github](https://github.com/Ensembl/WiggleTools)):


**Run wiggletools:**
```bash
conda activate BedToBigwig
sbatch scripts/bigwigmerge_TPM.sh # 12452262 ok 
```
*NOTE: bigwig are merge into 1 bedgraph which is then converted into 1 bigwig (wiggletools cannot output bigwig directly so need to pass by bedgraph or wiggle in between)*



# Calculate TPM and RPKM

Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
# Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_TPM_hg38.sh # 12377786 ok
# mv all output to output/tpm or rpkm folder
mv output/featurecounts_hg38/*tpm* output/tpm_hg38/
mv output/featurecounts_hg38/*rpkm* output/rpkm_hg38/
```

All good. 

If needed to display gene with TPM:

```R
# Plot with TPM instead of baseMean (Naiara plot)
## import tpm
#### Generate TPM for ALL samples
#### collect all samples ID
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
                     "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
                     "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
                     "4wN_WT_R1", "4wN_WT_R2", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4",
                     "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
                     "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3",
                     "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
                     "4wN_HET_R1", "4wN_HET_R2", "4wN_HET_R3", "4wN_HET_R4",
                     "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4",
                     "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
                     "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
                     "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
                     "4wN_KO_R1", "4wN_KO_R2", "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4",
                     "4wN_iPSCWT_R1", "4wN_iPSCWT_R2", "4wN_iPSCpatient_R1", "4wN_iPSCpatient_R2", "8wN_iPSCpatient_R1", "8wN_iPSCpatient_R2")

## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/tpm_hg38/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR_hg38.")) %>%
    rename(!!sample := starts_with("output.STAR_hg38."))
}

## Merge all dataframe into a single one
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")
write.csv(tpm_all_sample, file="../001__RNAseq/output/tpm_hg38/tpm_all_sample.txt")
### If need to import: tpm_all_sample <- read_csv("../001__RNAseq/output/tpm_hg38/tpm_all_sample.txt") #To import


XXX look for gene of interest XXX
```


## DEGs with deseq2

### 'one-by-one' comparison
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
- 4wN WT vs 4wN iPSCpatient
- 8wN WT vs 8wN iPSCpatient
- WT ESC vs NPC
- HET ESC vs NPC
- KO ESC vs NPC


### NPC KO vs WT
Take ressource
```bash
module load R/4.2.2
srun --mem=100g --pty bash -l
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_NPC_KO_vs_NPC_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_NPC_KO_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_NPC_KO_vs_NPC_WT.pdf", width=5, height=4)
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_NPC_HET_vs_NPC_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_NPC_HET_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_NPC_HET_vs_NPC_WT.pdf", width=5, height=4)
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_ESC_KO_vs_ESC_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_KO_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_ESC_KO_vs_ESC_WT.pdf", width=5, height=4)
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_ESC_HET_vs_ESC_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_ESC_HET_vs_ESC_WT.pdf", width=5, height=4)
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_2dN_KO_vs_2dN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_2dN_KO_vs_2dN_WT.pdf", width=5, height=4)
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_2dN_HET_vs_2dN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_2dN_HET_vs_2dN_WT.pdf", width=5, height=4)
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_4wN_KO_vs_4wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_4wN_KO_vs_4wN_WT.pdf", width=5, height=4)
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

## somthing weird with our samples, try different comparison
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_HET_R1", "4wN_HET_R2")
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_HET_R3" ,"4wN_HET_R4") # Comparison to choose
samples <- c("4wN_HET_R1", "4wN_HET_R2" ,"4wN_HET_R3" ,"4wN_HET_R4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_4wN_HET_vs_4wN_WT.txt")

write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_4wN_HET_R3R4_vs_4wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_4wN_HET_vs_4wN_WT.pdf", width=5, height=4)

pdf("output/deseq2_hg38/plotMA_res_4wN_HET_R3R4_vs_4wN_WT.pdf", width=5, height=4)
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
--> Here HET R3 and R4 is preferable to choose for diff analyses, it seems that R1 and R2 present same maturation level as the WT, thus no diff expr genes

### 8wN KO vs WT
Take ressource
```bash
srun --mem=50g --pty bash -l
module load R/4.2.2
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_8wN_KO_vs_8wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_8wN_KO_vs_8wN_WT.pdf", width=5, height=4)
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_8wN_HET_vs_8wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_8wN_HET_vs_8wN_WT.pdf", width=5, height=4)
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


### 8wN HET vs KO
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
samples <- c("8wN_HET_R1",
   "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
dds$genotype <- relevel(dds$genotype, ref = "HET")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_HET", type="apeglm")

## Export result as 'raw_8wN_KO_vs_8wN_HET.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_8wN_KO_vs_8wN_HET.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_8wN_KO_vs_8wN_HET.pdf", width=5, height=4)
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_4wN_iPSCWT_vs_4wN_iPSCpatient.txt")
### If need to import: res <- read_csv("output/deseq2/raw_ESC_HET_vs_ESC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_4wN_iPSCWT_vs_4wN_iPSCpatient.pdf", width=5, height=4)
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
--> The 4wN samples are weird because the WT is weird, we decided to use the WT to compare with iPSCpatient, and not the iPSCWT

### 4wN WT vs 4wN iPSCpatient
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
samples <- c("4wN_WT_R1", "4wN_WT_R2", "4wN_iPSCpatient_R1",
   "4wN_iPSCpatient_R2")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
res <- lfcShrink(dds, coef="genotype_iPSCpatient_vs_WT", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_4wN_iPSCpatient_vs_4wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_NPC_HET_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_4wN_iPSCpatient_vs_4wN_WT.pdf", width=5, height=4)
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

### 8wN WT vs 8wN iPSCpatient
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
samples <- c("8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_iPSCpatient_R1",
   "8wN_iPSCpatient_R2")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
res <- lfcShrink(dds, coef="genotype_iPSCpatient_vs_WT", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_8wN_iPSCpatient_vs_8wN_WT.txt")
### If need to import: res <- read_csv("output/deseq2/raw_NPC_HET_vs_NPC_WT.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_8wN_iPSCpatient_vs_8wN_WT.pdf", width=5, height=4)
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


### WT ESC vs NPC
Take ressource
```bash
srun --mem=200g --pty bash -l
conda activate deseq2
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
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3", "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
                              design= ~ time)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$time <- relevel(dds$time, ref = "ESC")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="time_NPC_vs_ESC", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_WT_ESC_vs_WT_NPC.txt")
### If need to import: res <- read_csv("output/deseq2/raw_WT_ESC_vs_WT_NPC.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_WT_ESC_vs_WT_NPC.pdf", width=5, height=4)
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



### HET ESC vs NPC
Take ressource
```bash
srun --mem=50g --pty bash -l
conda activate deseq2
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
samples <- c("ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3", "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
                              design= ~ time)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$time <- relevel(dds$time, ref = "ESC")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="time_NPC_vs_ESC", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_HET_ESC_vs_HET_NPC.txt")
### If need to import: res <- read_csv("output/deseq2/raw_WT_ESC_vs_WT_NPC.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_HET_ESC_vs_HET_NPC.pdf", width=5, height=4)
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




### KO ESC vs NPC
Take ressource
```bash
srun --mem=50g --pty bash -l
conda activate deseq2
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
samples <- c("ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3", "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
                              design= ~ time)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$time <- relevel(dds$time, ref = "ESC")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="time_NPC_vs_ESC", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_KO_ESC_vs_KO_NPC.txt")
### If need to import: res <- read_csv("output/deseq2/raw_WT_ESC_vs_WT_NPC.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_KO_ESC_vs_KO_NPC.pdf", width=5, height=4)
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


### WT NPC vs 2dN

```bash
srun --mem=50g --pty bash -l
conda activate deseq2
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
samples <- c("NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3", "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
                              design= ~ time)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$time <- relevel(dds$time, ref = "NPC")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="time_2dN_vs_NPC", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_WT_NPC_vs_WT_2dN.txt")
### If need to import: res <- read_csv("output/deseq2/raw_WT_ESC_vs_WT_NPC.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_WT_NPC_vs_WT_2dN.pdf", width=5, height=4)
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


### HET NPC vs 2dN

```bash
srun --mem=50g --pty bash -l
conda activate deseq2
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
samples <- c("NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3", "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
                              design= ~ time)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$time <- relevel(dds$time, ref = "NPC")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="time_2dN_vs_NPC", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_HET_NPC_vs_HET_2dN.txt")
### If need to import: res <- read_csv("output/deseq2/raw_WT_ESC_vs_WT_NPC.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_HET_NPC_vs_HET_2dN.pdf", width=5, height=4)
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



### KO NPC vs 2dN

```bash
srun --mem=50g --pty bash -l
conda activate deseq2
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
samples <- c("NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3", "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
                              design= ~ time)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$time <- relevel(dds$time, ref = "NPC")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="time_2dN_vs_NPC", type="apeglm")

## Export result as 'raw_NPC_HET_vs_NPC_WT.txt'
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_KO_NPC_vs_KO_2dN.txt")
### If need to import: res <- read_csv("output/deseq2/raw_WT_ESC_vs_WT_NPC.txt") #To import

## Plot-MA
pdf("output/deseq2_hg38/plotMA_res_KO_NPC_vs_KO_2dN.pdf", width=5, height=4)
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





## PCA and clustering with deseq2 hg38

Changes according to previous analyses with hg19:
- iPSCWT removed (weird)
- 4wN HET Rep3 and4 only (maturation state)


*NOTE: Tons of deseq2 ressource [here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) and [here](https://f1000research.com/articles/4-1070/v2).*

Open R/4.2.2 with ressource:
```bash
module load R/4.2.2
srun --mem=200g --pty bash -l
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

# Filter out 4wN HET R1 R2 and iPSCWT
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "4wN_WT_R1", "4wN_WT_R2", "4wN_KO_R1",
   "4wN_KO_R2",
   "4wN_HET_R3", "4wN_HET_R4", "4wN_iPSCpatient_R1", "4wN_iPSCpatient_R2",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4", "8wN_iPSCpatient_R1", "8wN_iPSCpatient_R2",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

# Only WT and KO
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "4wN_WT_R1", "4wN_WT_R2", "4wN_KO_R1",
   "4wN_KO_R2",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3")



## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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
pdf("output/deseq2_hg38/heatmap_cluster_vsd.pdf", width=5, height=6)
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
pdf("output/deseq2_hg38/heatmap_cluster_rld.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
pdf("output/deseq2_hg38/heatmap_cluster_rld_WT-KO.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

## PCA
### vsd 
pdf("output/deseq2_hg38/PCA_vsd.pdf", width=10, height=10)
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
pdf("output/deseq2_hg38/PCA_rld.pdf", width=10, height=10)
pcaData <- plotPCA(rld, intgroup=c("time", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()
dev.off()

pdf("output/deseq2_hg38/PCA_rld_WT-KO.pdf", width=10, height=10)
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

--> Clustering rld is slightly better (2dN vs NPC better clustered)


## time-course analyses with deseq2 hg38

Changes according to previous analyses with hg19:
- iPSCWT all removed (even though initially only iPSCWT was weird)
- 4wN HET Rep3 and4 only (maturation state)

```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")

# import featurecounts output and keep only gene ID and counts
# Filtered 4wN HET R1 R2 and iPSCWT and iPSCpatient
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "4wN_WT_R1", "4wN_WT_R2", "4wN_KO_R1",
   "4wN_KO_R2",
   "4wN_HET_R3", "4wN_HET_R4",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

### Additional filtering if needed
#### Only WT and HET
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
   "4wN_WT_R1", "4wN_WT_R2", 
   "4wN_HET_R3", "4wN_HET_R4",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")


#### Only WT and KO
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
   "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3",
   "4wN_WT_R1", "4wN_WT_R2", "4wN_KO_R1",
   "4wN_KO_R2",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4",
   "ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3")

#### collect all samples ID
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
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
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

## Construct the DESeqDataSet Time-Course 
### desgin = full-model
ddsTC <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                                colData = coldata,
                                design = ~ genotype + time + genotype:time)


### Define the reduced model (not including interaction term for comparison with full model)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ genotype + time)
resTC <- results(ddsTC)
write.csv(resTC %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/resTC.txt")
# resTC <- read.csv("output/deseq2_hg38/resTC.txt", row.names = 1)

# Data normalization
vst_counts <- vst(ddsTC, blind=FALSE)
save(vst_counts, file = "output/deseq2_hg38/ddsTC_vsd_filter.RData")

rlog_counts <- rlog(ddsTC, blind=FALSE) # last 5min
save(rlog_counts, file = "output/deseq2_hg38/ddsTC_rld_filter.RData")
load("output/deseq2_hg38/ddsTC_rld_filter.RData")


## Extract the transformed counts
vst_counts_matrix <- assay(vst_counts)
rlog_counts_matrix <- assay(rlog_counts) 

## Isolate the significant genes by gene id
signif_TC_genes = as_tibble(resTC, rownames = "gene") %>%
  filter(padj<= 0.05) %>%               # !!! Here change qvalue accordingly !!!
  select(gene) %>%
  unique() 

## Convert into vector 
signif_TC_genes_vector <- signif_TC_genes$gene

# VST counts

## Filter our matrix with the significant genes
vst_counts_matrix_sig <- vst_counts_matrix[rownames(vst_counts_matrix) %in% signif_TC_genes_vector, ]
nrow(vst_counts_matrix_sig) # double-check the nb of genes is same s in signif_TC_genes

# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(vst_counts_matrix_sig, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into k clusters
row_clusters <- cutree(row_hclust, k = 20)                   # !!! Here change tree nb accordingly !!!

## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(vst_counts_matrix_sig),
                           cluster = row_clusters)


### Save dataframe
write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_vst_20cl.txt")

# Make a clean table with significant deseq2-TC genes
vst_counts_tidy <- as_tibble(vst_counts_matrix_sig, rownames = "gene") %>%
  gather(key = "sample", value = "vst_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

## Compil both
vst_counts_tidy <- vst_counts_tidy %>%
  left_join(cluster_gene, by = "gene")

vst_counts_tidy$time <-
  factor(vst_counts_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
vst_counts_tidy$genotype <-
  factor(vst_counts_tidy$genotype,
         c("WT", "KO", "HET"))

## Calculate the number of genes per cluster
genes_per_cluster <- vst_counts_tidy %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2_hg38/line_vst_p0.05_cl25_pretty.pdf", width=20, height=14)
ggplot(vst_counts_tidy, aes(x = time, y = vst_counts)) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  # Add mean value for each genotype at each time point
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  # Add standard error bars around the mean points
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  # Add number of genes per cluster to the facet_wrap panels
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()


# RLOG counts

## Filter our matrix with the significant genes
rlog_counts_matrix_sig <- rlog_counts_matrix[rownames(rlog_counts_matrix) %in% signif_TC_genes_vector, ]
nrow(rlog_counts_matrix_sig) # double-check the nb of genes is same s in signif_TC_genes

# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(rlog_counts_matrix_sig, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into k clusters
row_clusters <- cutree(row_hclust, k = 8)                   # !!! Here change tree nb accordingly !!!
## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(rlog_counts_matrix_sig),
                           cluster = row_clusters)
### Save dataframe
write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_30cl.txt") # !!! Here change tree nb accordingly !!
write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_8cl_WTvsHET.txt")

# Make a clean table with significant deseq2-TC genes
rlog_counts_tidy <- as_tibble(rlog_counts_matrix_sig, rownames = "gene") %>%
  gather(key = "sample", value = "rlog_counts", -gene) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")

## Compil both
rlog_counts_tidy <- rlog_counts_tidy %>%
  left_join(cluster_gene, by = "gene")

rlog_counts_tidy$time <-
  factor(rlog_counts_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
rlog_counts_tidy$genotype <-
  factor(rlog_counts_tidy$genotype,
         c("WT", "KO", "HET"))

## Calculate the number of genes per cluster
genes_per_cluster <- rlog_counts_tidy %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

pdf("output/deseq2_hg38/line_rlog_p0.05_cl30_pretty_noSmooth.pdf", width=20, height=14)   # !!! Here change tree nb accordingly !!
ggplot(rlog_counts_tidy, aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
  geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, inherit.aes = FALSE) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()




# WT vs HET (NaiaraPlot)
pdf("output/deseq2_hg38/line_rlog_p0.05_cl8_pretty_noSmooth_WTvsHET.pdf", width=20, height=10)   # !!! Here change tree nb accordingly !!
ggplot(rlog_counts_tidy, aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
  geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, inherit.aes = FALSE) +
  facet_wrap(~cluster, scale = "free", nrow = 2) +
  scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16),       # Increase x-axis legend text size
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),       # Increase x-axis legend text size
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16) 
  )
dev.off()



# Plot rlog_transform norm deseq2 count with 'loess' method pretty with a light transparent color
pdf("output/deseq2_hg38/line_rlog_p0.05_cl25_pretty_span08_genes_bg.pdf", width=20, height=14)
ggplot(mean_rlog_counts_tidy, aes(x = time, y = mean_rlog_counts)) +
  geom_line(aes(color = genotype, group = interaction(gene, genotype)), alpha = 0.1) +
  geom_smooth(aes(color = genotype, group = genotype), method = "loess", se = TRUE, span = 0.8) +
  stat_summary(aes(color = genotype, group = genotype), fun = mean, geom = "point", shape = 18, size = 3, stroke = 1.5) +
  stat_summary(aes(color = genotype, group = genotype), fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5) +
  facet_wrap(~cluster, scale = "free", nrow = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),
    axis.title.x = element_text(size = 16)
  )
dev.off()

```
- *NOTE: count are saved under `output/deseq2_hg38/`: `vsd_filter.RData` and `rld_filter.RData` --> filter because HET 4wN R1/R2 and iPSCpatient removed*
- ***NOTE: I used now `blind=FALSE` for vsd/rld normalization !! Better for downstream analyses!!***


--> vst/rlog 25 clusters is very bad, lot of clusters with very few genes 

--> Let's prioritize rlog, since it perform better for unbiased (blind=TRUE) clustering

--> Using noSmooth is MUCH better!!!

--> The more cluster the better, as we are inform about all possible profiles, we can then focus on some clusters



# Gene ontology_hg38
We will use clusterProfile package. Tutorial [here](https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html).

Let's do a test of the pipeline with genes from cluster4 amd cluster14 from the rlog counts. Our background list will be all genes tested for differential expression.

**IMPORTANT NOTE: When doing GO, do NOT set a universe (background list of genes) it perform better!**

```R
# packages
library(clusterProfiler)
library(pathview)
library(DOSE)
library(org.Hs.eg.db)
library(enrichplot)
library(rtracklayer)
library(tidyverse)

## Read GTF file
gtf_file <- "../../Master/meta/ENCFF159KBI.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()



# WT vs HET (NairaPlot)
## Import genes_cluster list and background list
cluster_6 = read_csv("output/deseq2_hg38/cluster_gene_rlog_8cl_WTvsHET.txt") %>%
  filter(cluster == 6) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name) 


background = read_csv("output/deseq2_hg38/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


## Run GO enrichment analysis 
### with contrast
ego <- enrichGO(gene = as.character(cluster_6$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                qvalueCutoff = 0.05, 
                readable = TRUE)

### more relaxed parameter
ego <- enrichGO(gene = as.character(cluster_6$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
###

## Save GO analyses
GO_summary <- data.frame(ego)

write.csv(GO_summary, "output/GO_hg38/cluster4_BP.csv")
write.csv(GO_summary, "output/GO_hg38/cluster4_MF.csv")
write.csv(GO_summary, "output/GO_hg38/cluster4_CC.csv")



## Vizualization

pdf("output/GO_hg38/dotplot_BP_cluster_6.pdf", width=7, height=22)
dotplot(ego, showCategory=50)
dev.off()
pdf("output/GO_hg38/emapplot_BP_cluster_6.pdf", width=12, height=14)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()


pdf("output/GO_hg38/dotplot_BP_cluster_6_small15.pdf", width=5, height=8)
dotplot(ego, showCategory=15)
dev.off()
pdf("output/GO_hg38/emapplot_BP_cluster_6_small.pdf", width=8, height=9)
emapplot(pairwise_termsim(ego), showCategory = 20)
dev.off()



## Cluster 3/4/6/7/8 together (NaiaraPlot)
cluster_3 = read_csv("output/deseq2_hg38/cluster_gene_rlog_8cl_WTvsHET.txt") %>%
  filter(cluster == 3) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()

cluster_4 = read_csv("output/deseq2_hg38/cluster_gene_rlog_8cl_WTvsHET.txt") %>%
  filter(cluster == 4)  %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
  
cluster_6 = read_csv("output/deseq2_hg38/cluster_gene_rlog_8cl_WTvsHET.txt") %>%
  filter(cluster == 6)  %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
  
cluster_7 = read_csv("output/deseq2_hg38/cluster_gene_rlog_8cl_WTvsHET.txt") %>%
  filter(cluster == 7)  %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
  
cluster_8 = read_csv("output/deseq2_hg38/cluster_gene_rlog_8cl_WTvsHET.txt") %>%
  filter(cluster == 8)  %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
  

# Extract the 'entrez_id' column as a vector
cluster_3_vector <- pull(cluster_3, entrez_id)
cluster_4_vector <- pull(cluster_4, entrez_id)
cluster_6_vector <- pull(cluster_6, entrez_id)
cluster_7_vector <- pull(cluster_7, entrez_id)
cluster_8_vector <- pull(cluster_8, entrez_id)

# Combine the vectors into a list
entrez_list <- list(cluster_3 = cluster_3_vector, 
                    cluster_4 = cluster_4_vector, 
                    cluster_6 = cluster_6_vector, 
                    cluster_7 = cluster_7_vector, 
                    cluster_8 = cluster_8_vector)




## GO
compGO <- compareCluster(geneCluster   = entrez_list,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
pdf("output/GO_hg38/dotplot_BP_cluster34678.pdf", width=15, height=15)
dotplot(compGO, showCategory = 15, title = "GO_Biological Process Enrichment Analysis")
dev.off()



# Other comparionsv

```






