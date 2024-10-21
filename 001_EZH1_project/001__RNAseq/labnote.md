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

samples <- c("8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4")



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
### Genotype only
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype )

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
pdf("output/deseq2/PCA_rld_8wN.pdf", width=10, height=10)
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


pdf("output/deseq2/genes_HET_H3K27me3atESC_allgenotypes.pdf", width=10, height=6)
tpm_all_sample_tidy_gene_name_stat %>%
  filter(gene_name %in% c("WNT10B","SKOR2","ROR2","IRX3","ISL2","PITX2","PTPRM","ATOH1","BMP4","LYN",     "NTF3","IRX6","ATP8B1","EMX1","LHX8","OTP","NRN1"),
         genotype %in% c("WT", "KO", "HET"),
         time %in% c("ESC","NPC")) %>%
    ggplot(., aes(x = time, y = mean, group = genotype)) +
    geom_line(aes(color=genotype), size=0.75) +
    geom_errorbar(aes(ymin = mean-SE, ymax = mean+SE,color=genotype), width=.2) +
    geom_point(aes(y = mean,color=genotype), size = .75, shape = 15) +
    theme_bw() +
    facet_wrap(~gene_name, nrow = 3, scale = "free")  +	
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

If needed to **display gene with TPM**:

```R
# library
library("biomaRt")

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




# Display some genes in TPM: 
# ---> The code below is not perfect; issue at the geneSymbol conversin; to troubleshoot later; but it work
#### Generate TPM for ALL samples
#### collect all samples ID
samples <- c("8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4",
                     "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4","8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4")

samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
                     "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
                     "2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
                     "4wN_WT_R1", "4wN_WT_R2", "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4")

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

# plot some genes
tpm_all_sample_tidy <- tpm_all_sample %>%
  gather(key = 'variable', value = 'tpm', -Geneid) %>%
  separate(variable, into = c('time', 'genotype', 'replicate'), sep = "_") %>%
  rename(gene = Geneid)

tpm_all_sample_tidy$gene <- gsub("\\..*", "", tpm_all_sample_tidy$gene)

## convert gene Ensembl to symbol 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensembl gene IDs to gene symbols
genesymbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = tpm_all_sample_tidy$gene,
                     mart = ensembl)

# Merge gene symbols to your dataframe
tpm_all_sample_tidy <- left_join(tpm_all_sample_tidy, genesymbols, 
                                 by = c("gene" = "ensembl_gene_id"))

# autism-related genes 
c("SOX2", "GLI2", "EZH2", "PTN") # opposite behavior and well known
 c("SHANK3", "FMR1", "MECP2", "TSC1", "TSC2", "NLGN3", "NLGN4", "AUTS2", "CHD8", "SYNGAP1", "PTEN", "CNTNAP2", "ADNP", "FOXP1") # well known
c("NR2F2", "RELN", "GAD1", "GAD2", "GRIN2A", "GRIN2B", "NRD2A", "NRD2B", "RXRG") # labmeeting genes
c("GABRA1", "GABRA4", "ADCY8", "CACNA1B", "GABRG2", "GABRB3", "CACNA1S", "GABRB1", "GABRA5", "GNG3", "GABRD", "ADCY5", "SLC38A2", "GABRA2") # labmeeting genes_GABRA


plot_data <- tpm_all_sample_tidy %>%
  unique() %>%
  filter(external_gene_name %in% c("EZH1")) %>%
  group_by(gene, genotype,external_gene_name, time) %>%
  summarise(mean_log2tpm = mean(log2(tpm + 1)),
            se_log2tpm = sd(log2(tpm + 1)) / sqrt(n())) %>%
  ungroup()


plot_data$genotype <-
  factor(plot_data$genotype,
         c("WT", "HET","KO"))

        
# Plot
pdf("output/tpm_hg38/autism_genes_opposite_Disease_Gene_and_Drug_Signatures_from_GEO.pdf", width=5, height=4)
pdf("output/tpm_hg38/autism_genes_opposite_Disease_Gene_and_Drug_Signatures_from_GEO_wellKnown.pdf", width=8, height=4)
pdf("output/tpm_hg38/Labmeeting_20230822_genes.pdf", width=6, height=4)
pdf("output/tpm_hg38/Labmeeting_20230822_genes_GABRA.pdf", width=8, height=4)

ggplot(plot_data, aes(x = external_gene_name, y = mean_log2tpm, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = mean_log2tpm - se_log2tpm, ymax = mean_log2tpm + se_log2tpm),
    width = 0.25,
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  theme_bw() +
  ylab("log2(TPM + 1)") +
  xlab("Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



## PLOT line with statistic

plot_data <- tpm_all_sample_tidy %>%
  unique() %>%
  filter(external_gene_name %in% c("EZH1"),
         time != "4wN") %>%
  mutate(log2tpm = log2(tpm+1))
### Summarize data to calculate mean and standard deviation for each time point
plot_data_summary <- plot_data %>%
  group_by(time) %>%
  summarize(
    mean_log2tpm = mean(log2tpm),
    sd_log2tpm = sd(log2tpm)
  )
### Perform pairwise comparisons for each time point compared to ESC
comparisons <- list(
  c("ESC", "NPC"),
  c("ESC", "2dN"),
  c("ESC", "8wN")
)
### Generate the line plot with error bars and statistical comparison
pdf("output/tpm_hg38/linePlot_WT_EZH1.pdf", width=4, height=4)
ggline(data = plot_data, x = "time", y = "log2tpm", 
       add = "mean_se", color = "black") +
  stat_compare_means(comparisons = comparisons, 
                     label = "p.signif", 
                     method = "t.test",
                     ref.group = "ESC") +
  labs(
    x = "Time",
    y = "Log2(TPM+1)",
    title = "EZH1"
  ) +
  theme_bw()+
  theme(
    axis.text = element_text(size = 14)  # Increase the size of x and y axis tick labels
  )
dev.off()




# Here is updated code to test; that may add staitsict:
# Remove gene versions (e.g., ".5" from "ENSG000.5")
tpm_all_sample_tidy$gene <- gsub("\\..*", "", tpm_all_sample_tidy$gene)

# Create a Mart object for Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensembl gene IDs to gene symbols
genesymbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = unique(tpm_all_sample_tidy$gene),
                     mart = ensembl)

# Merge gene symbols to your dataframe
tpm_all_sample_tidy <- left_join(tpm_all_sample_tidy, genesymbols, 
                                 by = c("gene" = "ensembl_gene_id"))

# Rename columns if needed
tpm_all_sample_tidy <- tpm_all_sample_tidy %>%
  rename(geneSymbol = external_gene_name) %>%
  # Remove potential duplicates
  distinct()

# Filtering the data for the genes of interest
genes_of_interest <- c("ENSG001", "ENSG002")  # Replace with your actual gene IDs
filtered_data <- tpm_all_sample_tidy %>% 
  filter(gene %in% genes_of_interest)

# Plotting
p <- ggplot(filtered_data, aes(x = geneSymbol, y = log2(tpm + 1), fill = genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(y = "log2(TPM + 1)", title = "Gene expression levels") +
  scale_fill_manual(values = c("WT" = "black", "HET" = "blue", "KO" = "red"))

# Adding statistics using ggpubr
p_stat <- stat_compare_means(aes(group = genotype), label = "p.signif", 
                             comparisons = list(c("WT", "HET"), c("WT", "KO")), 
                             method = "t.test")

p + p_stat
```


# DEGs with deseq2 (and CutRun integration)

**IMPORTANT NOTE: Here it is advisable to REMOVE all genes from chromosome X and Y BEFORE doing the DEGs analysis (X chromosome re-activation occurs in some samples, notably these with more cell passage; in our case, the HET and KO)**
--> It is good to do this on the count matrix see [here](https://support.bioconductor.org/p/119932/)
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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  dplyr::select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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
write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_NPC_KO_vs_NPC_WT_filtered.txt")
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



## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols




# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'




pdf("output/deseq2_hg38/plotVolcano_res_q05_NPC_KO_vs_NPC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, NPC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_q05_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, NPC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()


upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q05_NPC_KO_vs_NPC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_NPC_KO_vs_NPC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_NPC_KO_vs_NPC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'




pdf("output/deseq2_hg38/plotVolcano_res_q01_NPC_KO_vs_NPC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, NPC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_q01_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, NPC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()


upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q01_NPC_KO_vs_NPC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01_NPC_KO_vs_NPC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01_NPC_KO_vs_NPC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


write.table(res, file = "output/deseq2_hg38/filtered_NPC_KO_vs_NPC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) 



# Plot CutRun RNAseq integration (Jasmine CutRun 005__)

## import gene list

# H3K27me3
### GAIN



H3K27me3_qval50_Gain = read.table("../005__CutRun_NPC_PSC/output/ChIPseeker/H3K27me3_annot_gain_qval50_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = H3K27me3_qval50_Gain %>% 
  left_join(res_tibble) 


### LOST
H3K27me3_qval50_Lost = read.table("../005__CutRun_NPC_PSC/output/ChIPseeker/H3K27me3_annot_lost_qval50_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()



#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = H3K27me3_qval50_Lost %>% 
  left_join(res_tibble)

## PLOT
### GAIN
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Gain_H3K27me3_qval50_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Gain_H3K27me3_qval50_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Gain[res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, ]
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, ]
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_KO_vs_NPC_WT_Gain_H3K27me3_qval50.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_KO_vs_NPC_WT_Gain_H3K27me3_qval50.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### LOST
highlight_genes <- c("") # NA

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Lost_H3K27me3_qval50_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Lost_H3K27me3_qval50_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Lost[res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, ]
upregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Lost[res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, ]
downregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_KO_vs_NPC_WT_Lost_H3K27me3_qval50.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_KO_vs_NPC_WT_Lost_H3K27me3_qval50.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





# EZH2
### GAIN



EZH2_qval10_Gain = read.table("../005__CutRun_NPC_PSC/output/ChIPseeker/EZH2_annot_gain_qval10_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = EZH2_qval10_Gain %>% 
  left_join(res_tibble) 


### LOST
EZH2_qval10_Lost = read.table("../005__CutRun_NPC_PSC/output/ChIPseeker/EZH2_annot_lost_qval10_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()



#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = EZH2_qval10_Lost %>% 
  left_join(res_tibble)

## PLOT
### GAIN
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Gain_EZH2_qval10_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Gain_EZH2_qval10_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Gain[res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, ]
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, ]
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_KO_vs_NPC_WT_Gain_EZH2_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_KO_vs_NPC_WT_Gain_EZH2_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### LOST
highlight_genes <- c("") # NA

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Lost_EZH2_qval10_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Lost_EZH2_qval10_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Lost[res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, ]
upregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Lost[res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, ]
downregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_KO_vs_NPC_WT_Lost_EZH2_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_KO_vs_NPC_WT_Lost_EZH2_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)










# SUZ12
### GAIN

SUZ12_qval10_Gain = read.table("../005__CutRun_NPC_PSC/output/ChIPseeker/SUZ12_annot_gain_qval10_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = SUZ12_qval10_Gain %>% 
  left_join(res_tibble) 


### LOST
SUZ12_qval10_Lost = read.table("../005__CutRun_NPC_PSC/output/ChIPseeker/SUZ12_annot_lost_qval10_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()



#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = SUZ12_qval10_Lost %>% 
  left_join(res_tibble)

## PLOT
### GAIN
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Gain_SUZ12_qval10_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Gain_SUZ12_qval10_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Gain[res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, ]
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, ]
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_KO_vs_NPC_WT_Gain_SUZ12_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_KO_vs_NPC_WT_Gain_SUZ12_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### LOST
highlight_genes <- c("") # NA

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Lost_SUZ12_qval10_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Lost_SUZ12_qval10_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Lost[res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, ]
upregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Lost[res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, ]
downregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_KO_vs_NPC_WT_Lost_SUZ12_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_KO_vs_NPC_WT_Lost_SUZ12_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)







# H3K4me3
### GAIN

H3K4me3_qval10_Gain = read.table("../005__CutRun_NPC_PSC/output/ChIPseeker/H3K4me3_annot_gain_qval10_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()


#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Gain = H3K4me3_qval10_Gain %>% 
  left_join(res_tibble) 


### LOST
H3K4me3_qval10_Lost = read.table("../005__CutRun_NPC_PSC/output/ChIPseeker/H3K4me3_annot_lost_qval10_promoterAnd5_geneSymbol.txt", 
                                           header = FALSE, 
                                           col.names = "GeneSymbol") %>%
                               as_tibble()



#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene") %>%
  drop_na()   # ADDING THIS AVOID THE BUG WITH DUPPLCIATED NAME 

res_Lost = H3K4me3_qval10_Lost %>% 
  left_join(res_tibble)

## PLOT
### GAIN
highlight_genes <- c("") # 

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Gain_H3K4me3_qval10_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Gain_H3K4me3_qval10_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Gain[res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, ]
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, ]
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_KO_vs_NPC_WT_Gain_H3K4me3_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_KO_vs_NPC_WT_Gain_H3K4me3_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### LOST
highlight_genes <- c("") # NA

# FILTER ON QVALUE 0.01 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Lost_H3K4me3_qval10_NPC_KO_vs_NPC_WT.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Lost_H3K4me3_qval10_NPC_KO_vs_NPC_WT_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Lost[res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, ]
upregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Lost[res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, ]
downregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_KO_vs_NPC_WT_Lost_H3K4me3_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_KO_vs_NPC_WT_Lost_H3K4me3_qval10.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)






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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")


# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 


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
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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




## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'




pdf("output/deseq2_hg38/plotVolcano_res_q05_NPC_HET_vs_NPC_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, NPC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_NPC_HET_vs_NPC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_NPC_HET_vs_NPC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.01 FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_NPC_HET_vs_NPC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, NPC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)



### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_NPC_HET_vs_NPC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_NPC_HET_vs_NPC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Export complete table with geneSymbol ID (without X Y chromosome)
write.table(res, file = "output/deseq2_hg38/filtered_NPC_HET_vs_NPC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) 
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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")


# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1" ,"ESC_KO_R2" ,"ESC_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

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
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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




## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'




pdf("output/deseq2_hg38/plotVolcano_res_q05_ESC_KO_vs_ESC_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.01 FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_ESC_KO_vs_ESC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)



### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Export complete table with geneSymbol ID (without X Y chromosome)
write.table(res, file = "output/deseq2_hg38/filtered_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) 


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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_HET_R1" ,"ESC_HET_R2" ,"ESC_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 


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
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'




pdf("output/deseq2_hg38/plotVolcano_res_q05_ESC_HET_vs_ESC_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, ESC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_ESC_HET_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_ESC_HET_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.01 FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_ESC_HET_vs_ESC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, ESC',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)



### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_ESC_HET_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_ESC_HET_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Export complete table with geneSymbol ID (without X Y chromosome)
write.table(res, file = "output/deseq2_hg38/filtered_ESC_HET_vs_ESC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) 




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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2" ,"2dN_WT_R3",
   "2dN_KO_R1" ,"2dN_KO_R2" ,"2dN_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

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
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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






## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'




pdf("output/deseq2_hg38/plotVolcano_res_q05_2dN_KO_vs_2dN_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 2dN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_2dN_KO_vs_2dN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_2dN_KO_vs_2dN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.01 FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_2dN_KO_vs_2dN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 2dN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)



### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_2dN_KO_vs_2dN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_2dN_KO_vs_2dN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Export complete table with geneSymbol ID (without X Y chromosome)
write.table(res, file = "output/deseq2_hg38/filtered_2dN_KO_vs_2dN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) 

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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("2dN_WT_R1", "2dN_WT_R2" ,"2dN_WT_R3",
   "2dN_HET_R1" ,"2dN_HET_R2" ,"2dN_HET_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 


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
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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



## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'




pdf("output/deseq2_hg38/plotVolcano_res_q05_2dN_HET_vs_2dN_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 2dN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_2dN_HET_vs_2dN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_2dN_HET_vs_2dN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.01 FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_2dN_HET_vs_2dN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 2dN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)



### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_2dN_HET_vs_2dN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_2dN_HET_vs_2dN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Export complete table with geneSymbol ID (without X Y chromosome)
write.table(res, file = "output/deseq2_hg38/filtered_2dN_HET_vs_2dN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) 




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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")


# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_KO_R1",
   "4wN_KO_R2" )

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

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
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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



## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q05_4wN_KO_vs_2dN_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 4wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_4wN_KO_vs_4wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_4wN_KO_vs_4wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.01 FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_4wN_KO_vs_4wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 4wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)



### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_4wN_KO_vs_4wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_4wN_KO_vs_4wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Export complete table with geneSymbol ID (without X Y chromosome)
write.table(res, file = "output/deseq2_hg38/filtered_4wN_KO_vs_4wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) 


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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_HET_R1", "4wN_HET_R2",
   "4wN_HET_R3" ,"4wN_HET_R4")

## somthing weird with our samples, try different comparison
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_HET_R1", "4wN_HET_R2")
samples <- c("4wN_HET_R1", "4wN_HET_R2" ,"4wN_HET_R3" ,"4wN_HET_R4")
samples <- c("4wN_WT_R1", "4wN_WT_R2" ,"4wN_HET_R3" ,"4wN_HET_R4") # Comparison to choose. THE GOOD ONE !!


## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}


# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  dplyr::select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'




pdf("output/deseq2_hg38/plotVolcano_res_q05_4wN_HET_vs_4wN_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 4wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q05_4wN_HET_vs_4wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_4wN_HET_vs_4wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_4wN_HET_vs_4wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.01 FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_4wN_HET_vs_4wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 4wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q01FC05_4wN_HET_vs_4wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_4wN_HET_vs_4wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_4wN_HET_vs_4wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Export complete table with geneSymbol ID (without X Y chromosome)
write.table(res, file = "output/deseq2_hg38/filtered_4wN_HET_vs_4wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) 



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
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")


# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("8wN_WT_R1", "8wN_WT_R2" ,"8wN_WT_R3" ,"8wN_WT_R4" ,"8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  dplyr::select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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
# write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_8wN_KO_vs_8wN_WT.txt")
# write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/filtered_8wN_KO_vs_8wN_WT.txt") --> FILTERED IS HERE WITHOUT X AND Y GENES
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


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols


## FILTER on pvalue NOT GOOD!!!


keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 10e-6, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 10e-6, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_8wN_KO_vs_8wN_WT.pdf", width=7, height=6)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'WT vs LOF',
  pCutoff = 10e-6,
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1)  + 
  theme_bw() 
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$pvalue < 10e-6)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$pvalue < 10e-6)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$pvalue < 10e-6, ]
#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$pvalue < 10e-6, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


# Post GO analysis; volcano plot with highlighted genes
## KEGG2016__KO
# Set up gene symbols and their colors
highlight_genes <- c("EPHB6", "EPHA5", "SEMA5B", "EPHA4", "SEMA4A", "NTN4", "CXCR4", "UNC5D", "SEMA3F", "NFATC4", "RND1", "CXCL12", "SLIT1", "PLXNB2", "PLXNB1", "SRGAP2", "EPHB1", "SRGAP1", "EPHB4", "NGEF", "EPHA3")
highlight_genes <- c("PDGFRB", "MAGI1", "CSF1", "ANGPT1", "LPAR2", "PIK3R3", "ARAP3", "FGF1", "GRIN2B", "ACTB", "EGFR", "RAP1GAP", "ACTG1", "ADORA2A", "PDGFD", "P2RY1", "PLCE1", "PRKD1", "DRD2", "FGFR3", "FGF12", "FGFR2", "PFN2", "RAPGEF4")
highlight_genes <- c("NOTCH2", "LFNG", "NOTCH3", "NOTCH1", "JAG1", "MFNG", "DLL1", "HES5")
highlight_genes <- c("NOTCH2", "LFNG", "NOTCH3", "NOTCH1", "JAG1", "MFNG", "DLL1", "HES5")
highlight_genes <- c("NR2F2", "RELN", "GAD1", "GAD2", "GRIN2A", "GRIN2B", "NRD2A", "NRD2B", "RXRG") # labmeeting genes


# Create custom key-value pairs for gene highlighting
keyvals <- ifelse(res$GeneSymbol %in% highlight_genes, 'midnight blue', 'gray88')
names(keyvals)[keyvals == 'midnight blue'] <- 'Highlighted Genes'
names(keyvals)[keyvals == 'gray88'] <- 'Others'

pdf("output/deseq2_hg38/plotVolcano_res_8wN_KO_vs_8wN_WT_AxonGuidance.pdf", width=10, height=8)    
pdf("output/deseq2_hg38/plotVolcano_res_8wN_KO_vs_8wN_WT_Rap1Signaling.pdf", width=10, height=8)    
pdf("output/deseq2_hg38/plotVolcano_res_8wN_KO_vs_8wN_WT_NotchSignaling.pdf", width=10, height=8)    
pdf("output/deseq2_hg38/plotVolcano_res_8wN_KO_vs_8wN_WT_labmeetingGenes20230822.pdf", width=10, height=8)    

EnhancedVolcano(res,
    lab = res$GeneSymbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = res$GeneSymbol[which(names(keyvals) == 'Highlighted Genes')],
    pCutoff = 10e-6,
    FCcutoff = 0.5,   
    pointSize = 1.5,
    labSize = 4.5,
    shape = 20,
    colCustom = keyvals,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'black',
    max.overlaps = 100,
    arrowheads = FALSE) +theme_bw() # +xlim(-2.5,0)+ylim(0,20)
dev.off()



# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_q05_8wN_KO_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_q05_8wN_KO_vs_8wN_WT_prettyV1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q05_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.05 and FC 1
keyvals <- ifelse(
  res$log2FoldChange < -1 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 1 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 1)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -1)'

pdf("output/deseq2_hg38/plotVolcano_res_q05FC1_8wN_KO_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 1,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 1 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -1 & res$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q05FC1_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 1 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 1 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -1 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -1 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC1_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC1_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



# FILTER ON QVALUE 0.01 and FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_8wN_KO_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q01FC05_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 1 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -1 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





# FILTER ON QVALUE 0.01 and FC 1
keyvals <- ifelse(
  res$log2FoldChange < -1 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 1 & res$padj < 1e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 1)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -1)'

pdf("output/deseq2_hg38/plotVolcano_res_q01FC1_8wN_KO_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 1,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 1 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -1 & res$padj < 1e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q01FC1_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 1 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 1 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -1 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -1 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC1_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC1_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



# FILTER ON QVALUE 0.05 and FC 2 (NaiaraPlot 20240119)
keyvals <- ifelse(
  res$log2FoldChange < -2 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 2 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 2)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -2)'

pdf("output/deseq2_hg38/plotVolcano_res_q05FC2_8wN_KO_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 2,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 2 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -2 & res$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q05FC2_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 2 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 2 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -1 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -2 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC2_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC2_8wN_KO_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



# Plot CutRun RNAseq integration (PosterMidatlantic)
## import gene list
### GAIN
THOR_qval15_KO_Gain_DEG_Down = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Down.txt", 
                                           header = FALSE, 
                                           col.names = "gene") %>%
                               as_tibble()
THOR_qval15_KO_Gain_DEG_Up = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Up.txt", 
                                           header = FALSE, 
                                           col.names = "gene") %>%
                               as_tibble()

THOR_qval15_KO_Gain = THOR_qval15_KO_Gain_DEG_Down %>%
  bind_rows(THOR_qval15_KO_Gain_DEG_Up)
#### Convert tibble into character vector for GeneSymbol conversion
THOR_qval15_KO_Gain_char_vector <- as.character(THOR_qval15_KO_Gain$gene)

THOR_qval15_KO_Gain_gene_symbols <- mapIds(org.Hs.eg.db, keys = THOR_qval15_KO_Gain_char_vector,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

THOR_qval15_KO_Gain$GeneSymbol <- THOR_qval15_KO_Gain_gene_symbols
#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene")

res_Gain = THOR_qval15_KO_Gain %>% 
  left_join(res_tibble)

### LOST
THOR_qval15_KO_Lost_DEG_Down = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Down.txt", 
                                           header = FALSE, 
                                           col.names = "gene") %>%
                               as_tibble()
THOR_qval15_KO_Lost_DEG_Up = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Up.txt", 
                                           header = FALSE, 
                                           col.names = "gene") %>%
                               as_tibble()                               


THOR_qval15_KO_Lost = THOR_qval15_KO_Lost_DEG_Down %>%
  bind_rows(THOR_qval15_KO_Lost_DEG_Up)
#### Convert tibble into character vector for GeneSymbol conversion
THOR_qval15_KO_Lost_char_vector <- as.character(THOR_qval15_KO_Lost$gene)

THOR_qval15_KO_Lost_gene_symbols <- mapIds(org.Hs.eg.db, keys = THOR_qval15_KO_Lost_char_vector,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

THOR_qval15_KO_Lost$GeneSymbol <- THOR_qval15_KO_Lost_gene_symbols
#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene")

res_Lost = THOR_qval15_KO_Lost %>% 
  left_join(res_tibble)

## PLOT
### GAIN
highlight_genes <- c("GRIN1", "ADCY1", "GRM8", "GRIK3", "ADCY5") # Glutamatergic genes
highlight_genes <- c("") # none

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Gain_8wN_KO_vs_8wN_WT_glutamatergicSynapse.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Gain_8wN_KO_vs_8wN_WT_glutamatergicSynapse_prettyV1.pdf", width=8, height=8)  
pdf("output/deseq2_hg38/plotVolcano_res_Gain_8wN_KO_vs_8wN_WT_prettyV1.pdf", width=8, height=8)  

EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Gain[res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, ]
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, ]
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_8wN_KO_vs_8wN_WT_Gain_H3K27me3_poster.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_8wN_KO_vs_8wN_WT_Gain_H3K27me3_poster.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### LOST
highlight_genes <- c("") # NA

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Lost_8wN_KO_vs_8wN_WT_glutamatergicSynapse.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Lost_8wN_KO_vs_8wN_WT_glutamatergicSynapse_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Lost[res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, ]
upregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Lost[res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, ]
downregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_8wN_KO_vs_8wN_WT_Lost_H3K27me3_poster.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_8wN_KO_vs_8wN_WT_Lost_H3K27me3_poster.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





################## Naira 20240119 tasks --- Plot Cut Run integration

## Re do more properly the  CutRun expression integration 

### import the 679? 346? genes that gain H3K27me3 in promoter in KO
WTvsKO_annot_gain <- as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "gain") %>%
    dplyr::select(gene,geneSymbol,H3K27me3) %>% 
    unique() %>%
    rename(GeneSymbol = geneSymbol)
WTvsKO_annot_lost <- as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC < 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "lost") %>%
    dplyr::select(gene,geneSymbol,H3K27me3) %>% 
    unique() %>%
    rename(GeneSymbol = geneSymbol)
WTvsKO_annot_gain_lost = WTvsKO_annot_gain %>%
  bind_rows(WTvsKO_annot_lost)

#### save output gene list gain lost
write.table(as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "gain") %>%
    unique() %>%
    rename(GeneSymbol = geneSymbol), file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15_gain_Promoter5.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 
write.table(as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "gain") %>%
    dplyr::select(gene,geneSymbol,H3K27me3) %>% 
    unique() %>%
    rename(GeneSymbol = geneSymbol), file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15_gain_Promoter5_geneSymbol.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 
write.table(as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC < 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "lost") %>%
    unique() %>%
    rename(GeneSymbol = geneSymbol), file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15_lost_Promoter5.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC < 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "lost") %>%
    dplyr::select(gene,geneSymbol,H3K27me3) %>% 
    unique() %>%
    rename(GeneSymbol = geneSymbol), file = "../003__CutRun/output/ChIPseeker/annotation_WTvsKO_unique_Keepdup_qval15_lost_Promoter5_geneSymbol.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


### add gene expression information
#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene")

#### combine expression and gain lost
WTvsKO_annot_gain_lost_expression = WTvsKO_annot_gain_lost %>%
  left_join(res_tibble %>% dplyr::select(gene, log2FoldChange,pvalue,padj)) %>%
  unique() %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
         padj = ifelse(is.na(padj), 1, padj),
         pvalue = ifelse(is.na(pvalue), 1, pvalue))
##### save output: write.table(WTvsKO_annot_gain_lost_expression, file = "output/deseq2_hg38/WTvsKO_annot_gain_lost_expression.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
WTvsKO_annot_gain_expression = WTvsKO_annot_gain_lost_expression %>%
  filter(H3K27me3 == "gain")
WTvsKO_annot_lost_expression = WTvsKO_annot_gain_lost_expression %>%
  filter(H3K27me3 == "lost")

#### volcano plot ------- GAIN
###### FILTER ON QVALUE 0.05 FC 05
keyvals <- ifelse(
  WTvsKO_annot_gain_expression$log2FoldChange < -0.5 & WTvsKO_annot_gain_expression$padj < 5e-2, 'Sky Blue',
    ifelse(WTvsKO_annot_gain_expression$log2FoldChange > 0.5 & WTvsKO_annot_gain_expression$padj < 5e-2, 'Orange',
      'grey'))


keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.5)'


pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsHET_annot_gain.pdf", width=7, height=8)    
EnhancedVolcano(WTvsKO_annot_gain_expression,
  lab = WTvsKO_annot_gain_expression$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 8wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsHET_annot_gain_nolabel.pdf", width=7, height=8)    
EnhancedVolcano(WTvsKO_annot_gain_expression,
  lab = c(""),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 8wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()
upregulated_genes <- sum(WTvsKO_annot_gain_expression$log2FoldChange > 0.5 & WTvsKO_annot_gain_expression$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(WTvsKO_annot_gain_expression$log2FoldChange < -0.5 & WTvsKO_annot_gain_expression$padj < 5e-2, na.rm = TRUE)

#### mean log2fc signficiant vs non signfiicant
WTvsKO_annot_gain_expression %>%
  filter(log2FoldChange < -0.5, padj < 5e-2) %>%
  summarise(mean_log2FoldChange = mean(log2FoldChange, na.rm = TRUE))
WTvsKO_annot_gain_expression_signif = WTvsKO_annot_gain_expression %>%
  filter(log2FoldChange < -0.5, padj < 5e-2) %>%
  dplyr::select(gene,GeneSymbol)
WTvsKO_annot_gain_expression %>%
anti_join(WTvsKO_annot_gain_expression_signif) %>%
  summarise(mean_log2FoldChange = mean(log2FoldChange, na.rm = TRUE))
####

# Save as gene list for GO analysis
#### Filter for up-regulated genes
upregulated <- WTvsKO_annot_gain_expression[WTvsKO_annot_gain_expression$log2FoldChange > 0.5 & WTvsKO_annot_gain_expression$padj < 5e-2, ]
upregulated <- WTvsKO_annot_gain_expression[!is.na(WTvsKO_annot_gain_expression$log2FoldChange) & !is.na(WTvsKO_annot_gain_expression$padj) & WTvsKO_annot_gain_expression$log2FoldChange > 0.5 & WTvsKO_annot_gain_expression$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- WTvsKO_annot_gain_expression[WTvsKO_annot_gain_expression$log2FoldChange < -0.5 & WTvsKO_annot_gain_expression$padj < 5e-2, ]
downregulated <- WTvsKO_annot_gain_expression[!is.na(WTvsKO_annot_gain_expression$log2FoldChange) & !is.na(WTvsKO_annot_gain_expression$padj) & WTvsKO_annot_gain_expression$log2FoldChange < -1 & WTvsKO_annot_gain_expression$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_gain.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_gain.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


#### volcano plot ------- LOST
###### FILTER ON QVALUE 0.05 FC 05
keyvals <- ifelse(
  WTvsKO_annot_lost_expression$log2FoldChange < -0.5 & WTvsKO_annot_lost_expression$padj < 5e-2, 'Sky Blue',
    ifelse(WTvsKO_annot_lost_expression$log2FoldChange > 0.5 & WTvsKO_annot_lost_expression$padj < 5e-2, 'Orange',
      'grey'))


keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.5)'


pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsHET_annot_lost.pdf", width=7, height=8)    
EnhancedVolcano(WTvsKO_annot_lost_expression,
  lab = WTvsKO_annot_lost_expression$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 8wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsHET_annot_lost_nolabel.pdf", width=7, height=8)    
EnhancedVolcano(WTvsKO_annot_lost_expression,
  lab = c(""),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, 8wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()
upregulated_genes <- sum(WTvsKO_annot_lost_expression$log2FoldChange > 0.5 & WTvsKO_annot_lost_expression$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(WTvsKO_annot_lost_expression$log2FoldChange < -0.5 & WTvsKO_annot_lost_expression$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis
#### Filter for up-regulated genes
upregulated <- WTvsKO_annot_lost_expression[WTvsKO_annot_lost_expression$log2FoldChange > 0.5 & WTvsKO_annot_lost_expression$padj < 5e-2, ]
upregulated <- WTvsKO_annot_lost_expression[!is.na(WTvsKO_annot_lost_expression$log2FoldChange) & !is.na(WTvsKO_annot_lost_expression$padj) & WTvsKO_annot_lost_expression$log2FoldChange > 0.5 & WTvsKO_annot_lost_expression$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- WTvsKO_annot_lost_expression[WTvsKO_annot_lost_expression$log2FoldChange < -0.5 & WTvsKO_annot_lost_expression$padj < 5e-2, ]
downregulated <- WTvsKO_annot_lost_expression[!is.na(WTvsKO_annot_lost_expression$log2FoldChange) & !is.na(WTvsKO_annot_lost_expression$padj) & WTvsKO_annot_lost_expression$log2FoldChange < -1 & WTvsKO_annot_lost_expression$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_lost.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_lost.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


############################

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
library("EnhancedVolcano")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("8wN_WT_R1", "8wN_WT_R2" ,"8wN_WT_R3" ,"8wN_WT_R4" ,"8wN_HET_R1",
   "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

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
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  dplyr::select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

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
# write.csv(res %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/raw_8wN_HET_vs_8wN_WT.txt")


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



## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$pvalue < 10e-6, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$pvalue < 10e-6, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; FC < -0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_8wN_HET_vs_8wN_WT.pdf", width=7, height=6)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'WT vs LOF',
  pCutoff = 10e-6,
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1)  + 
  theme_bw() 
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$pvalue < 10e-6)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$pvalue < 10e-6)

# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$pvalue < 10e-6, ]
#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$pvalue < 10e-6, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

## KEGG2016__HET
# Set up gene symbols and their colors
highlight_genes <- c("ITGB1", "CSF1", "ITGB5", "ITGB4", "LPAR1", "TNC", "BRCA1", "LAMC1", "FGF1", "FGF2", "THBS1", "EGFR", "THBS4", "GNG10", "CCND2", "BCL2L11", "CCND1", "GNG5", "EIF4EBP1", "TNR", "MCL1", "IFNAR2", "PDGFRA", "ANGPT1", "HGF", "RPS6", "GNG12", "EFNA1", "KITLG", "COL1A2", "CDK6", "CCNE2", "COL4A2", "COL4A1", "DDIT4", "CDK2", "ITGA7", "ITGA6", "FGFR3")
highlight_genes <- c("GABRB3", "VIPR1", "GABRB1", "NPFFR1", "GIPR", "CHRM5", "OPRL1", "ADRB1", "ADRA1A", "HTR6", "GRM4", "GRM6", "P2RY1", "KISS1R", "DRD1", "DRD2", "GABRD", "S1PR4", "GABRA2", "GABRA1", "P2RY11", "NPY5R", "GPR35", "GABRA5", "GRID1", "GABRA4", "HTR1D", "NPY1R", "HTR1B", "SSTR1", "ADRA2C", "HTR5A", "SSTR3", "GABRG2", "GRIN1", "P2RX5", "GRIN3A", "MC1R", "ADORA2A", "GLRB", "F2RL1")
highlight_genes <- c("GABRB3", "GABRA2", "GABRA1", "GABRB1", "GABRA5", "GABRA4", "CACNA1B", "ADCY8", "GABRG2", "ADCY5", "GNG3", "CACNA1S", "SLC38A2", "GABRD")
highlight_genes <- c("YAP1", "CRB1", "WWC1", "FGF1", "NKD1", "GLI2", "SOX2", "CCND2", "CCND1", "CDH1", "TEAD1", "TEAD2", "TEAD3", "FZD1", "WWTR1", "SMAD1", "TGFB2", "FZD7", "WNT5A", "FZD6", "BMP7", "TGFBR2", "MOB1A", "FRMD6", "APC", "ID2", "BIRC5", "CTNNB1", "BMPR1B", "BIRC2")
highlight_genes <- c("NR2F2", "RELN", "GAD1", "GAD2", "GRIN2A", "GRIN2B", "NRD2A", "NRD2B", "RXRG") # labmeeting genes


# Create custom key-value pairs for gene highlighting # 'midnight blue'  OR 'darkorange3'
keyvals <- ifelse(res$GeneSymbol %in% highlight_genes, 'darkorange3', 'gray88')
names(keyvals)[keyvals == 'darkorange3'] <- 'Highlighted Genes'
names(keyvals)[keyvals == 'gray88'] <- 'Others'

pdf("output/deseq2_hg38/plotVolcano_res_8wN_HET_vs_8wN_WT_PI3Kpathway.pdf", width=10, height=8)    
pdf("output/deseq2_hg38/plotVolcano_res_8wN_HET_vs_8wN_WT_neuroactiveLigand.pdf", width=10, height=8)  
pdf("output/deseq2_hg38/plotVolcano_res_8wN_HET_vs_8wN_WT_GABAergicSynapse.pdf", width=10, height=8) 
pdf("output/deseq2_hg38/plotVolcano_res_8wN_HET_vs_8wN_WT_Hippopathway.pdf", width=10, height=8)        
pdf("output/deseq2_hg38/plotVolcano_res_8wN_HET_vs_8wN_WT_labmeetingGenes20230822.pdf", width=10, height=8)         

EnhancedVolcano(res,
    lab = res$GeneSymbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = res$GeneSymbol[which(names(keyvals) == 'Highlighted Genes')],
    pCutoff = 10e-6,
    FCcutoff = 0.5,   
    pointSize = 1.5,
    labSize = 4.5,
    shape = 20,
    colCustom = keyvals,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'black',
    max.overlaps = 100,
    arrowheads = FALSE) +theme_bw() # +xlim(0,3)+ylim(0,55)
dev.off()






# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q05_8wN_HET_vs_8wN_WT.pdf", width=7, height=8)    

EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_q05_8wN_HET_vs_8wN_WT_prettyV1.pdf", width=7, height=8)   
pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_HET_vs_8wN_WT.pdf", width=7, height=8)    

EnhancedVolcano(res,
  lab = "",
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q05_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





# FILTER ON QVALUE 0.01 FC 0.5
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < 0.5)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC05_8wN_HET_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.5 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -0.5 & res$padj < 1e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q01FC05_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.5 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.01 FC 1
keyvals <- ifelse(
  res$log2FoldChange < -1 & res$padj < 1e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 1 & res$padj < 1e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.01; log2FC > 1)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.01; log2FC < -1)'



pdf("output/deseq2_hg38/plotVolcano_res_q01FC1_8wN_HET_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 2wN',
  pCutoff = 1e-2,         #
  FCcutoff = 1,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 1 & res$padj < 1e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -1 & res$padj < 1e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q01FC1_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 1 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 1 & res$padj < 1e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -1 & res$padj < 1e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q01FC1_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q01FC1_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)






# FILTER ON QVALUE 0.05 FC 2
keyvals <- ifelse(
  res$log2FoldChange < -2 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 2 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 2)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -2)'



pdf("output/deseq2_hg38/plotVolcano_res_q05FC2_8wN_HET_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 2,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 2 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -2 & res$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q05FC2_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 1 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 2 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -2 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC2_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC2_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




# FILTER ON QVALUE 0.05 FC 1
keyvals <- ifelse(
  res$log2FoldChange < -1 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 1 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 1)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -1)'



pdf("output/deseq2_hg38/plotVolcano_res_q05FC1_8wN_HET_vs_8wN_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 1,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 1 & res$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < -1 & res$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2_hg38/filtered_q05FC1_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 1 & res$padj < 1e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 1 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -1 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC1_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC1_8wN_HET_vs_8wN_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)









# Plot CutRun RNAseq integration (PosterMidatlantic)
## import gene list
### GAIN
THOR_qval15_HET_Gain_DEG_Down = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Down.txt", 
                                           header = FALSE, 
                                           col.names = "gene") %>%
                               as_tibble()
THOR_qval15_HET_Gain_DEG_Up = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Up.txt", 
                                           header = FALSE, 
                                           col.names = "gene") %>%
                               as_tibble()

THOR_qval15_HET_Gain = THOR_qval15_HET_Gain_DEG_Down %>%
  bind_rows(THOR_qval15_HET_Gain_DEG_Up)
#### Convert tibble into character vector for GeneSymbol conversion
THOR_qval15_HET_Gain_char_vector <- as.character(THOR_qval15_HET_Gain$gene)

THOR_qval15_HET_Gain_gene_symbols <- mapIds(org.Hs.eg.db, keys = THOR_qval15_HET_Gain_char_vector,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

THOR_qval15_HET_Gain$GeneSymbol <- THOR_qval15_HET_Gain_gene_symbols
#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene")

res_Gain = THOR_qval15_HET_Gain %>% 
  left_join(res_tibble)

### LOST
THOR_qval15_HET_Lost_DEG_Down = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Down.txt", 
                                           header = FALSE, 
                                           col.names = "gene") %>%
                               as_tibble()
THOR_qval15_HET_Lost_DEG_Up = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Up.txt", 
                                           header = FALSE, 
                                           col.names = "gene") %>%
                               as_tibble()                               


THOR_qval15_HET_Lost = THOR_qval15_HET_Lost_DEG_Down %>%
  bind_rows(THOR_qval15_HET_Lost_DEG_Up)
#### Convert tibble into character vector for GeneSymbol conversion
THOR_qval15_HET_Lost_char_vector <- as.character(THOR_qval15_HET_Lost$gene)

THOR_qval15_HET_Lost_gene_symbols <- mapIds(org.Hs.eg.db, keys = THOR_qval15_HET_Lost_char_vector,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

THOR_qval15_HET_Lost$GeneSymbol <- THOR_qval15_HET_Lost_gene_symbols
#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene")

res_Lost = THOR_qval15_HET_Lost %>% 
  left_join(res_tibble)

## PLOT
### GAIN
highlight_genes <- c("GRIN1", "SHANK1", "GRM4", "PRKCB", "PLCB2", "GRM1") # Glutamatergic genes

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Gain_8wN_HET_vs_8wN_WT_glutamatergicSynapse.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'HET vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Gain_8wN_HET_vs_8wN_WT_glutamatergicSynapse_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Gain[res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, ]
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.5 & res_Gain$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, ]
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.5 & res_Gain$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_8wN_HET_vs_8wN_WT_Gain_H3K27me3_poster.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_8wN_HET_vs_8wN_WT_Gain_H3K27me3_poster.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



### LOST
highlight_genes <- c("") # NA

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < 0.5)'

pdf("output/deseq2_hg38/plotVolcano_res_Lost_8wN_HET_vs_8wN_WT_glutamatergicSynapse.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.5,
  labSize = 4.5,
  shape = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_Lost_8wN_HET_vs_8wN_WT_glutamatergicSynapse_prettyV1.pdf", width=8, height=8)  
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = highlight_genes,
  title = 'KO vs WT, 2wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 5,
  labSize = 9,   # gene highlight size
  shape = 20,
  axisLabSize = 25,
  captionLabSize = 20,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  colConnectors = 'black',
  max.overlaps = 100,
  arrowheads = FALSE)  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) )
dev.off()

upregulated_genes <- sum(res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:

upregulated <- res_Lost[res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, ]
upregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange > 0.5 & res_Lost$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Lost[res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, ]
downregulated <- res_Lost[!is.na(res_Lost$log2FoldChange) & !is.na(res_Lost$padj) & res_Lost$log2FoldChange < -0.5 & res_Lost$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_8wN_HET_vs_8wN_WT_Lost_H3K27me3_poster.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_8wN_HET_vs_8wN_WT_Lost_H3K27me3_poster.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)






################## Naira 20240119 tasks --- Plot Cut Run integration

## Re do more properly the  CutRun expression integration 

### import the 679? 346? genes that gain H3K27me3 in promoter in HET
WTvsHET_annot_gain <- as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "gain") %>%
    dplyr::select(gene,geneSymbol,H3K27me3) %>% 
    unique() %>%
    rename(GeneSymbol = geneSymbol)
WTvsHET_annot_lost <- as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC < 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "lost") %>%
    dplyr::select(gene,geneSymbol,H3K27me3) %>% 
    unique() %>%
    rename(GeneSymbol = geneSymbol)
WTvsHET_annot_gain_lost = WTvsHET_annot_gain %>%
  bind_rows(WTvsHET_annot_lost)

#### save output gene list gain lost
write.table(as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "gain") %>%
    unique() %>%
    rename(GeneSymbol = geneSymbol), file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15_gain_Promoter5.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 
write.table(as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "gain") %>%
    dplyr::select(gene,geneSymbol,H3K27me3) %>% 
    unique() %>%
    rename(GeneSymbol = geneSymbol), file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15_gain_Promoter5_geneSymbol.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 
write.table(as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC < 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "lost") %>%
    unique() %>%
    rename(GeneSymbol = geneSymbol), file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15_lost_Promoter5.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(as_tibble(read.csv(file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)) %>%
    filter(FC < 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "lost") %>%
    dplyr::select(gene,geneSymbol,H3K27me3) %>% 
    unique() %>%
    rename(GeneSymbol = geneSymbol), file = "../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15_lost_Promoter5_geneSymbol.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
### add gene expression information
#### Remove gene version on the res and compil with THOR diff genes
rownames(res) <- gsub("\\..*", "", rownames(res))
res_tibble <- res %>% 
  as_tibble(rownames = "gene")

#### combine expression and gain lost
WTvsHET_annot_gain_lost_expression = WTvsHET_annot_gain_lost %>%
  left_join(res_tibble %>% dplyr::select(gene, log2FoldChange,pvalue,padj)) %>%
  unique() %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
         padj = ifelse(is.na(padj), 1, padj),
         pvalue = ifelse(is.na(pvalue), 1, pvalue))
##### save output: write.table(WTvsHET_annot_gain_lost_expression, file = "output/deseq2_hg38/WTvsHET_annot_gain_lost_expression.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
WTvsHET_annot_gain_expression = WTvsHET_annot_gain_lost_expression %>%
  filter(H3K27me3 == "gain")
WTvsHET_annot_lost_expression = WTvsHET_annot_gain_lost_expression %>%
  filter(H3K27me3 == "lost")

#### volcano plot ------- GAIN
###### FILTER ON QVALUE 0.05 FC 05
keyvals <- ifelse(
  WTvsHET_annot_gain_expression$log2FoldChange < -0.5 & WTvsHET_annot_gain_expression$padj < 5e-2, 'Sky Blue',
    ifelse(WTvsHET_annot_gain_expression$log2FoldChange > 0.5 & WTvsHET_annot_gain_expression$padj < 5e-2, 'Orange',
      'grey'))


keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.5)'


pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_gain.pdf", width=7, height=8)    
EnhancedVolcano(WTvsHET_annot_gain_expression,
  lab = WTvsHET_annot_gain_expression$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 8wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_gain_nolabel.pdf", width=7, height=8)    
EnhancedVolcano(WTvsHET_annot_gain_expression,
  lab = c(""),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 8wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()
upregulated_genes <- sum(WTvsHET_annot_gain_expression$log2FoldChange > 0.5 & WTvsHET_annot_gain_expression$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(WTvsHET_annot_gain_expression$log2FoldChange < -0.5 & WTvsHET_annot_gain_expression$padj < 5e-2, na.rm = TRUE)

#### mean log2fc signficiant vs non signfiicant
WTvsHET_annot_gain_expression %>%
  filter(log2FoldChange < -0.5, padj < 5e-2) %>%
  summarise(mean_log2FoldChange = mean(log2FoldChange, na.rm = TRUE))
WTvsHET_annot_gain_expression_signif = WTvsHET_annot_gain_expression %>%
  filter(log2FoldChange < -0.5, padj < 5e-2) %>%
  dplyr::select(gene,GeneSymbol)
WTvsHET_annot_gain_expression %>%
anti_join(WTvsHET_annot_gain_expression_signif) %>%
  summarise(mean_log2FoldChange = mean(log2FoldChange, na.rm = TRUE))
####

# Save as gene list for GO analysis
#### Filter for up-regulated genes
upregulated <- WTvsHET_annot_gain_expression[WTvsHET_annot_gain_expression$log2FoldChange > 0.5 & WTvsHET_annot_gain_expression$padj < 5e-2, ]
upregulated <- WTvsHET_annot_gain_expression[!is.na(WTvsHET_annot_gain_expression$log2FoldChange) & !is.na(WTvsHET_annot_gain_expression$padj) & WTvsHET_annot_gain_expression$log2FoldChange > 0.5 & WTvsHET_annot_gain_expression$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- WTvsHET_annot_gain_expression[WTvsHET_annot_gain_expression$log2FoldChange < -0.5 & WTvsHET_annot_gain_expression$padj < 5e-2, ]
downregulated <- WTvsHET_annot_gain_expression[!is.na(WTvsHET_annot_gain_expression$log2FoldChange) & !is.na(WTvsHET_annot_gain_expression$padj) & WTvsHET_annot_gain_expression$log2FoldChange < -1 & WTvsHET_annot_gain_expression$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_gain.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_gain.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


#### volcano plot ------- LOST
###### FILTER ON QVALUE 0.05 FC 05
keyvals <- ifelse(
  WTvsHET_annot_lost_expression$log2FoldChange < -0.5 & WTvsHET_annot_lost_expression$padj < 5e-2, 'Sky Blue',
    ifelse(WTvsHET_annot_lost_expression$log2FoldChange > 0.5 & WTvsHET_annot_lost_expression$padj < 5e-2, 'Orange',
      'grey'))


keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.5)'


pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_lost.pdf", width=7, height=8)    
EnhancedVolcano(WTvsHET_annot_lost_expression,
  lab = WTvsHET_annot_lost_expression$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 8wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()

pdf("output/deseq2_hg38/plotVolcano_res_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_lost_nolabel.pdf", width=7, height=8)    
EnhancedVolcano(WTvsHET_annot_lost_expression,
  lab = c(""),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'HET vs WT, 8wN',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
  pointSize = 2.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()
upregulated_genes <- sum(WTvsHET_annot_lost_expression$log2FoldChange > 0.5 & WTvsHET_annot_lost_expression$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(WTvsHET_annot_lost_expression$log2FoldChange < -0.5 & WTvsHET_annot_lost_expression$padj < 5e-2, na.rm = TRUE)


# Save as gene list for GO analysis
#### Filter for up-regulated genes
upregulated <- WTvsHET_annot_lost_expression[WTvsHET_annot_lost_expression$log2FoldChange > 0.5 & WTvsHET_annot_lost_expression$padj < 5e-2, ]
upregulated <- WTvsHET_annot_lost_expression[!is.na(WTvsHET_annot_lost_expression$log2FoldChange) & !is.na(WTvsHET_annot_lost_expression$padj) & WTvsHET_annot_lost_expression$log2FoldChange > 0.5 & WTvsHET_annot_lost_expression$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- WTvsHET_annot_lost_expression[WTvsHET_annot_lost_expression$log2FoldChange < -0.5 & WTvsHET_annot_lost_expression$padj < 5e-2, ]
downregulated <- WTvsHET_annot_lost_expression[!is.na(WTvsHET_annot_lost_expression$log2FoldChange) & !is.na(WTvsHET_annot_lost_expression$padj) & WTvsHET_annot_lost_expression$log2FoldChange < -1 & WTvsHET_annot_lost_expression$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2_hg38/upregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_lost.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2_hg38/downregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_lost.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


############################


```

--> **Naiara task 20240119**: I re-do more propelry the CutRun RNAseq integration; there was over-filtering; because I imported already signficiant genes that gain and down expr; let's instead import list of genes that gain HEt, and then incorprorate expression; so I did:
- Import all genes that gain / lost H3K27me3 in promoter (by NOT using `read.table()` this function is FUCKED !! [Not import all rows](https://www.biostars.org/p/221983/))
- Put expression information
- Plot enhanced volcano



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

# Only 8wN
samples <- c("8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4", "8wN_KO_R1",
   "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4", "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4")

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
### Genotype only
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype )
#### Time and genotype
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

pdf("output/deseq2_hg38/PCA_rld_8wN.pdf", width=10, height=10)
pcaData <- plotPCA(rld, intgroup=c("time", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("blue", "red", "black")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
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
#### Only WT and HET at 4 and 8 weeks
samples <- c("4wN_WT_R1", "4wN_WT_R2",
   "4wN_HET_R3", "4wN_HET_R4",
   "8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4",
   "8wN_HET_R1", "8wN_HET_R2",
   "8wN_HET_R3", "8wN_HET_R4")


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

#### Only ESC and NPC
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_KO_R1", "ESC_KO_R2", "ESC_KO_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

#### Only ESC and NPC; WT vs HET
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
   "ESC_HET_R1", "ESC_HET_R2", "ESC_HET_R3",
   "NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
   "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3")

#### Only ESC and NPC; WT vs KO
samples <- c("ESC_WT_R1", "ESC_WT_R2", "ESC_WT_R3",
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
write.csv(resTC %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/resTC_ESC_NPC.txt")
write.csv(resTC %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/resTC_ESC_NPC_WT_HET.txt")
write.csv(resTC %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/resTC_ESC_NPC_WT_KO.txt")
write.csv(resTC %>% as.data.frame() %>% rownames_to_column("gene") %>% as.tibble(), file="output/deseq2_hg38/resTC_4wN_8wN_WT_HET.txt")
# resTC <- read.csv("output/deseq2_hg38/resTC.txt", row.names = 1)

# Data normalization
vst_counts <- vst(ddsTC, blind=FALSE)
save(vst_counts, file = "output/deseq2_hg38/ddsTC_vsd_filter.RData")

rlog_counts <- rlog(ddsTC, blind=FALSE) # last 5min
save(rlog_counts, file = "output/deseq2_hg38/ddsTC_rld_filter.RData")
save(rlog_counts, file = "output/deseq2_hg38/ddsTC_rld_filter_ESC_NPC.RData")
save(rlog_counts, file = "output/deseq2_hg38/ddsTC_rld_filter_ESC_NPC_WT_HET.RData")
save(rlog_counts, file = "output/deseq2_hg38/ddsTC_rld_filter_ESC_NPC_WT_KO.RData")

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
row_clusters <- cutree(row_hclust, k = 5)                   # !!! Here change tree nb accordingly !!!
## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(rlog_counts_matrix_sig),
                           cluster = row_clusters)
### Save dataframe
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_30cl.txt") # !!! Here change tree nb accordingly !!
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_8cl_WTvsHET.txt")
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_15cl_ESC_NPC.txt")
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET.txt")
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO.txt")

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
pdf("output/deseq2_hg38/line_rlog_p0.05_cl15_pretty_noSmooth_ESC_NPC.pdf", width=20, height=14)   # !!! Here change tree nb accordingly !!
pdf("output/deseq2_hg38/line_rlog_p0.05_cl15_pretty_noSmooth_ESC_NPC_WT_HET.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!
pdf("output/deseq2_hg38/line_rlog_p0.05_cl5_pretty_noSmooth_ESC_NPC_WT_KO.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!
pdf("output/deseq2_hg38/line_rlog_p0.05_cl5_pretty_noSmooth_4wN_8wN_WT_HET.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!

ggplot(rlog_counts_tidy, aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
  geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, inherit.aes = FALSE) +
  facet_wrap(~cluster, scale = "free", nrow = 1) +
  scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()




## Put together with the H3K27me3-dynamic gene (from ChIPseq ESC vs NPC in WT)__HET
### File: output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt (THOR qval25)

H3K27me3_genes <- as_tibble(read.table("../002__ChIPseq/output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt", header = TRUE, stringsAsFactors = FALSE,sep="\t")) %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
  dplyr::select(gene, qvalue,FC,annotation,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2) %>%
  unique()

H3K27me3_genes <- as_tibble(read.table("../002__ChIPseq/output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt", header = TRUE, stringsAsFactors = FALSE,sep="\t")) %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
  dplyr::select(gene) %>%
  unique()

## Calculate the number of genes per cluster
genes_per_cluster <- H3K27me3_genes %>%
  left_join(rlog_counts_tidy %>% separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE)) %>%
  na.omit() %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

rlog_counts_tidy %>%
  separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE) %>%
  left_join(H3K27me3_genes) %>%
  na.omit()
 

pdf("output/deseq2_hg38/line_rlog_p0.05_cl15_pretty_noSmooth_ESC_NPC_WT_HET_H3K27me3genes.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!


H3K27me3_genes %>%
  left_join(rlog_counts_tidy %>% separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE)) %>%
  na.omit() %>%
    ggplot(., aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
      geom_line(stat = "summary", fun = mean) +
      geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
      geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
      geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, inherit.aes = FALSE) +
      facet_wrap(~cluster, scale = "free", nrow = 1) +
      scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
        axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
      )
    dev.off()



## Put together with the H3K27me3-dynamic gene (from ChIPseq ESC vs NPC in WT)__KO
### File: output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt (THOR qval25)

H3K27me3_genes <- as_tibble(read.table("../002__ChIPseq/output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt", header = TRUE, stringsAsFactors = FALSE,sep="\t")) %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
  dplyr::select(gene, qvalue,FC,annotation,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2) %>%
  unique()

H3K27me3_genes <- as_tibble(read.table("../002__ChIPseq/output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt", header = TRUE, stringsAsFactors = FALSE,sep="\t")) %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
  dplyr::select(gene) %>%
  unique()

## Calculate the number of genes per cluster
genes_per_cluster <- H3K27me3_genes %>%
  left_join(rlog_counts_tidy %>% separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE)) %>%
  na.omit() %>%
  group_by(cluster) %>%
  summarise(num_genes = n_distinct(gene))

rlog_counts_tidy %>%
  separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE) %>%
  left_join(H3K27me3_genes) %>%
  na.omit()
 

pdf("output/deseq2_hg38/line_rlog_p0.05_cl15_pretty_noSmooth_ESC_NPC_WT_KO_H3K27me3genes.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!


H3K27me3_genes %>%
  left_join(rlog_counts_tidy %>% separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE)) %>%
  na.omit() %>%
    ggplot(., aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
      geom_line(stat = "summary", fun = mean) +
      geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
      geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
      geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, inherit.aes = FALSE) +
      facet_wrap(~cluster, scale = "free", nrow = 1) +
      scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
        axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
      )
    dev.off()




### Repeat clustering on the new list of genes H3K27me3_regulated (postClustering)__HET
#### Filter the H3K27me3-bound within our signif TC list signif_TC_genes_vector
rlog_counts_matrix_gene = rlog_counts_matrix # remove the version of the gene in the rlog matrix
rownames(rlog_counts_matrix_gene) <- gsub("\\..*", "", rownames(rlog_counts_matrix))

## Filter the significant TC and H3K27me3 genes
H3K27me3_genes <- as_tibble(read.table("../002__ChIPseq/output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt", header = TRUE, stringsAsFactors = FALSE,sep="\t")) %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
  dplyr::select(gene) %>%
  unique()

signif_TC_H3K27me3 = signif_TC_genes %>%
  separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE) %>%
  dplyr::select(-version) %>%
  inner_join(H3K27me3_genes)

signif_TC_H3K27me3_vector <- signif_TC_H3K27me3$gene



## Filter our matrix with the significant genes
rlog_counts_matrix_gene_sig <- rlog_counts_matrix_gene[rownames(rlog_counts_matrix_gene) %in% signif_TC_H3K27me3_vector, ]
nrow(rlog_counts_matrix_gene_sig) # double-check the nb of genes is same s in signif_TC_genes

# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(rlog_counts_matrix_gene_sig, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into k clusters
row_clusters <- cutree(row_hclust, k = 5)                   # !!! Here change tree nb accordingly !!!
## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(rlog_counts_matrix_gene_sig),
                           cluster = row_clusters)
### Save dataframe
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt")

# Make a clean table with significant deseq2-TC genes
rlog_counts_tidy <- as_tibble(rlog_counts_matrix_gene_sig, rownames = "gene") %>%
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


pdf("output/deseq2_hg38/line_rlog_p0.05_cl15_pretty_noSmooth_ESC_NPC_WT_HET_H3K27me3postClustering.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!

ggplot(rlog_counts_tidy, aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
  geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, inherit.aes = FALSE) +
  facet_wrap(~cluster, scale = "free", nrow = 1) +
  scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()







### Repeat clustering on the CutRun diff bound genes__HET (../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15.txt)
#### Filter the H3K27me3-bound within our signif TC list signif_TC_genes_vector
rlog_counts_matrix_gene = rlog_counts_matrix # remove the version of the gene in the rlog matrix
rownames(rlog_counts_matrix_gene) <- gsub("\\..*", "", rownames(rlog_counts_matrix))

## Filter the significant TC and H3K27me3 genes AND THAT GAIN H3K27me3 in HET
H3K27me3_genes <- as_tibble(read.table("../003__CutRun/output/ChIPseeker/annotation_WTvsHET_unique_Keepdup_qval15.txt", header = TRUE, stringsAsFactors = FALSE,sep="\t")) %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"),
         FC > 1) %>%
  dplyr::select(gene) %>%
  drop_na() %>%
  unique()

signif_TC_H3K27me3 = signif_TC_genes %>%
  separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE) %>%
  dplyr::select(-version) %>%
  inner_join(H3K27me3_genes)

signif_TC_H3K27me3_vector <- signif_TC_H3K27me3$gene



## Filter our matrix with the significant genes
rlog_counts_matrix_gene_sig <- rlog_counts_matrix_gene[rownames(rlog_counts_matrix_gene) %in% signif_TC_H3K27me3_vector, ]
nrow(rlog_counts_matrix_gene_sig) # double-check the nb of genes is same s in signif_TC_genes

# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(rlog_counts_matrix_gene_sig, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into k clusters
row_clusters <- cutree(row_hclust, k = 4)                   # !!! Here change tree nb accordingly !!!
## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(rlog_counts_matrix_gene_sig),
                           cluster = row_clusters)
### Save dataframe
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO_H3K27me3postClustering.txt")
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_4cl_4wN_8wN_WT_HET_H3K27me3postClustering.txt")

# Make a clean table with significant deseq2-TC genes
rlog_counts_tidy <- as_tibble(rlog_counts_matrix_gene_sig, rownames = "gene") %>%
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


pdf("output/deseq2_hg38/line_rlog_p0.05_cl15_pretty_noSmooth_ESC_NPC_WT_KO_H3K27me3postClustering.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!
pdf("output/deseq2_hg38/line_rlog_p0.05_cl4_pretty_noSmooth_4wN_8wN_WT_HET_H3K27me3postClustering.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!
ggplot(rlog_counts_tidy, aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
  geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, inherit.aes = FALSE) +
  facet_wrap(~cluster, scale = "free", nrow = 1) +
  scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()






### Repeat clustering on the new list of genes H3K27me3_regulated (postClustering)__KO
#### Filter the H3K27me3-bound within our signif TC list signif_TC_genes_vector
rlog_counts_matrix_gene = rlog_counts_matrix # remove the version of the gene in the rlog matrix
rownames(rlog_counts_matrix_gene) <- gsub("\\..*", "", rownames(rlog_counts_matrix))

## Filter the significant TC and H3K27me3 genes
H3K27me3_genes <- as_tibble(read.table("../002__ChIPseq/output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt", header = TRUE, stringsAsFactors = FALSE,sep="\t")) %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
  dplyr::select(gene) %>%
  unique()

signif_TC_H3K27me3 = signif_TC_genes %>%
  separate(gene, into = c("gene", "version"), sep = "\\.", remove = TRUE) %>%
  dplyr::select(-version) %>%
  inner_join(H3K27me3_genes)

signif_TC_H3K27me3_vector <- signif_TC_H3K27me3$gene



## Filter our matrix with the significant genes
rlog_counts_matrix_gene_sig <- rlog_counts_matrix_gene[rownames(rlog_counts_matrix_gene) %in% signif_TC_H3K27me3_vector, ]
nrow(rlog_counts_matrix_gene_sig) # double-check the nb of genes is same s in signif_TC_genes

# Obtain gene cluster ID from the heatmap
## Perform hierarchical clustering
row_dist <- dist(rlog_counts_matrix_gene_sig, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")

## Cut the tree into k clusters
row_clusters <- cutree(row_hclust, k = 5)                   # !!! Here change tree nb accordingly !!!
## Create a data frame with gene names and their corresponding clusters
cluster_gene <- data.frame(gene = rownames(rlog_counts_matrix_gene_sig),
                           cluster = row_clusters)
### Save dataframe
# write.csv(cluster_gene, file="output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO_H3K27me3postClustering.txt")

# Make a clean table with significant deseq2-TC genes
rlog_counts_tidy <- as_tibble(rlog_counts_matrix_gene_sig, rownames = "gene") %>%
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


pdf("output/deseq2_hg38/line_rlog_p0.05_cl15_pretty_noSmooth_ESC_NPC_WT_KO_H3K27me3postClustering.pdf", width=20, height=5)   # !!! Here change tree nb accordingly !!

ggplot(rlog_counts_tidy, aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
  geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
  geom_text(data = genes_per_cluster, aes(label = paste0("Genes: ", num_genes), x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, inherit.aes = FALSE) +
  facet_wrap(~cluster, scale = "free", nrow = 1) +
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
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")

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


# ESC_NPC Timecourse WT vs HET and WT vs KO (associated with H3K27me3-dynamics WT genes)
## Files
output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt
output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO_H3K27me3postClustering.txt

### WT HET
cluster_1 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt") %>%
  filter(cluster == 1) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_2 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt") %>%
  filter(cluster == 2) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_3 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt") %>%
  filter(cluster == 3) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_4 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt") %>%
  filter(cluster == 4) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_5 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt") %>%
  filter(cluster == 5) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()


ego <- enrichGO(gene = as.character(cluster_1$entrez_id), 
                keyType = "ENTREZID",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO_hg38/dotplot_BP_WT_HET_H3K27me3postClustering_cluster_1.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_HET_H3K27me3postClustering_cluster_2.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_HET_H3K27me3postClustering_cluster_2.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_HET_H3K27me3postClustering_cluster_3.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_HET_H3K27me3postClustering_cluster_4.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_HET_H3K27me3postClustering_cluster_5.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()


enrichKEGG <- enrichKEGG(gene   = as.character(cluster_4$entrez_id),
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/GO_hg38/dotplot_KEGG_WT_HET_H3K27me3postClustering_cluster_1.pdf", width=7, height=5)
dotplot(enrichKEGG, showCategory=20)
dev.off()


### WT KO


cluster_1 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO_H3K27me3postClustering.txt") %>%
  filter(cluster == 1) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_2 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO_H3K27me3postClustering.txt") %>%
  filter(cluster == 2) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_3 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO_H3K27me3postClustering.txt") %>%
  filter(cluster == 3) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_4 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO_H3K27me3postClustering.txt") %>%
  filter(cluster == 4) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_5 = read_csv("output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_KO_H3K27me3postClustering.txt") %>%
  filter(cluster == 5) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()


ego <- enrichGO(gene = as.character(cluster_5$entrez_id), 
                keyType = "ENTREZID",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO_hg38/dotplot_BP_WT_KO_H3K27me3postClustering_cluster_1.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_KO_H3K27me3postClustering_cluster_2.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_KO_H3K27me3postClustering_cluster_3.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_KO_H3K27me3postClustering_cluster_4.pdf", width=7, height=7)
pdf("output/GO_hg38/dotplot_BP_WT_KO_H3K27me3postClustering_cluster_5.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()


enrichKEGG <- enrichKEGG(gene   = as.character(cluster_4$entrez_id),
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/GO_hg38/dotplot_KEGG_WT_HET_H3K27me3postClustering_cluster_1.pdf", width=7, height=5)
dotplot(enrichKEGG, showCategory=20)
dev.off()


# ESC vs NPC
## Import genes_cluster list 

cluster_1 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 1) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()

cluster_2 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 2) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()

cluster_3 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 3) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()

cluster_4 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 4) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()

cluster_5 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 5) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_6 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 6) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit() 
cluster_7 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 7) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_8 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 8) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_9 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 9) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_10 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 10) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_11 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 11) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()
cluster_12 = read_csv("output/deseq2_hg38/cluster_gene_rlog_12cl_ESC_NPC.txt") %>%
  filter(cluster == 12) %>%
  mutate(gene = sub("\\..*", "", gene), # Remove version number
  entrez_id = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")) %>%
  dplyr::select(entrez_id) %>%
  na.omit()

  

# Extract the 'entrez_id' column as a vector
cluster_1_vector <- pull(cluster_1, entrez_id)
cluster_2_vector <- pull(cluster_2, entrez_id)
cluster_3_vector <- pull(cluster_3, entrez_id)
cluster_4_vector <- pull(cluster_4, entrez_id)
cluster_5_vector <- pull(cluster_5, entrez_id)
cluster_6_vector <- pull(cluster_6, entrez_id)
cluster_7_vector <- pull(cluster_7, entrez_id)
cluster_8_vector <- pull(cluster_8, entrez_id)
cluster_9_vector <- pull(cluster_9, entrez_id)
cluster_10_vector <- pull(cluster_10, entrez_id)
cluster_11_vector <- pull(cluster_11, entrez_id)
cluster_12_vector <- pull(cluster_12, entrez_id)
# Combine the vectors into a list
# Combine the vectors into a list
entrez_list <- list(
  cluster_1 = cluster_1_vector,
  cluster_2 = cluster_2_vector,
  cluster_3 = cluster_3_vector,
  cluster_4 = cluster_4_vector,
  cluster_5 = cluster_5_vector,
  cluster_6 = cluster_6_vector,
  cluster_7 = cluster_7_vector,
  cluster_8 = cluster_8_vector,
  cluster_9 = cluster_9_vector,
  cluster_10 = cluster_10_vector,
  cluster_11 = cluster_11_vector,
  cluster_12 = cluster_12_vector
)





## GO
compGO <- compareCluster(geneCluster   = entrez_list,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
pdf("output/GO_hg38/dotplot_BP_cluster_ESC_NPC.pdf", width=15, height=50)
dotplot(compGO, showCategory = 15, title = "GO_Biological Process Enrichment Analysis")
dev.off()

## KEGG
compKEGG <- compareCluster(geneCluster   = entrez_list,
                         fun           = "enrichKEGG",
                         organism = "human",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
pdf("output/GO_hg38/dotplot_KEGG_cluster_ESC_NPC.pdf", width=12, height=18)
dotplot(compKEGG, showCategory = 15, title = "KEGG pathway Enrichment Analysis")
dev.off()




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





# Gene ontology_hg38__ cleaner representation for GO up and down; enrichR

Let's follow this nice [blog](http://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html) and this [discussion](https://www.biostars.org/p/343196/) for opposite barplot.

--> There are plenty of different databases to explore; see [here](https://maayanlab.cloud/Enrichr/#libraries)

```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")
library("enrichR")
library("biomaRt")

## Files

# Define databases for enrichment
dbs <- c("KEGG_2016") # NEED TO TRY KEGG2023 !!!


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt



# IF starting with EnesbmlID
### Ensembl ID; list of signif up/down with H3K27me3 changes in each genotypes
output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Up.txt
output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Down.txt
output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Up.txt
output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Down.txt
output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Up.txt
output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Down.txt
output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Up.txt
output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Down.txt
#### Convert to GeneSymbol
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_down = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Down.txt", header = FALSE, stringsAsFactors = FALSE) 
gene_down = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Down.txt", header = FALSE, stringsAsFactors = FALSE)
gene_names_down <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                        filters = "ensembl_gene_id", 
                        values = gene_down$V1, 
                        mart = ensembl) %>%
                   dplyr::select(external_gene_name)
list_down <- unique(as.character(gene_names_down$external_gene_name))
edown <- enrichr(list_down, dbs)


gene_up = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Up.txt", header = FALSE, stringsAsFactors = FALSE)
gene_up = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Up.txt", header = FALSE, stringsAsFactors = FALSE)

gene_names_up <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                        filters = "ensembl_gene_id", 
                        values = gene_up$V1, 
                        mart = ensembl)%>%
                   dplyr::select(external_gene_name)
list_up <- unique(as.character(gene_names_up$external_gene_name))
eup <- enrichr(list_up, dbs)



# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$KEGG_2016
down <- edown$KEGG_2016
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("Homo sapiens hsa[0-9]*", "", gos$Term)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_KEGG_8wN_KO_vs_8wN_WT.pdf", width=12, height=5)
pdf("output/GO_hg38/enrichR_KEGG_8wN_HET_vs_8wN_WT.pdf", width=12, height=6)
pdf("output/GO_hg38/enrichR_KEGG_8wN_HET_vs_8wN_WT__H3K27me3_DEG.pdf", width=12, height=6)
pdf("output/GO_hg38/enrichR_KEGG_8wN_KO_vs_8wN_WT__H3K27me3_DEG.pdf", width=12, height=6) # qvalue at 1.1!

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_KEGG_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_KEGG_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_KEGG_8wN_HET_vs_8wN_WT__H3K27me3_DEG.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_KEGG_8wN_KO_vs_8wN_WT__H3K27me3_DEG.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("KEGG_2021_Human") # 
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt

# IF starting with EnesbmlID
### Ensembl ID; list of signif up/down with H3K27me3 changes in each genotypes
output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Up.txt
output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Down.txt
output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Up.txt
output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Down.txt
output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Up.txt
output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Down.txt
output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Up.txt
output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Down.txt
#### Convert to GeneSymbol
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_down = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Down.txt", header = FALSE, stringsAsFactors = FALSE) 
gene_down = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Down.txt", header = FALSE, stringsAsFactors = FALSE)
gene_names_down <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                        filters = "ensembl_gene_id", 
                        values = gene_down$V1, 
                        mart = ensembl) %>%
                   dplyr::select(external_gene_name)
list_down <- unique(as.character(gene_names_down$external_gene_name))
edown <- enrichr(list_down, dbs)

gene_up = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Up.txt", header = FALSE, stringsAsFactors = FALSE)
gene_up = read.table("../003__CutRun/output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Up.txt", header = FALSE, stringsAsFactors = FALSE)

gene_names_up <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                        filters = "ensembl_gene_id", 
                        values = gene_up$V1, 
                        mart = ensembl)%>%
                   dplyr::select(external_gene_name)
list_up <- unique(as.character(gene_names_up$external_gene_name))
eup <- enrichr(list_up, dbs)




# IF starting with geneSymbol
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.05 and FC 0.5 _ NPC
output/deseq2_hg38/downregulated_q05FC05_NPC_KO_vs_NPC_WT.txt
output/deseq2_hg38/upregulated_q05FC05_NPC_KO_vs_NPC_WT.txt
output/deseq2_hg38/downregulated_q05FC05_NPC_HET_vs_NPC_WT.txt
output/deseq2_hg38/upregulated_q05FC05_NPC_HET_vs_NPC_WT.txt


## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q05FC05_NPC_KO_vs_NPC_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q05FC05_NPC_KO_vs_NPC_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)


# Extracting KEGG data and assigning types
up <- eup$KEGG_2021_Human
down <- edown$KEGG_2021_Human
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos <- gos %>% filter(abs(logAdjP) > 2.6)
## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)



# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_KEGG_2021_8wN_KO_vs_8wN_WT_q01F05.pdf", width=12, height=5)
pdf("output/GO_hg38/enrichR_KEGG_2021_8wN_HET_vs_8wN_WT_q01F05.pdf", width=12, height=8)
pdf("output/GO_hg38/enrichR_KEGG_2021_8wN_HET_vs_8wN_WT__H3K27me3_DEG.pdf", width=12, height=6)
pdf("output/GO_hg38/enrichR_KEGG_2021_8wN_KO_vs_8wN_WT__H3K27me3_DEG.pdf", width=12, height=6) # qvalue at 1.1!
pdf("output/GO_hg38/enrichR_KEGG_2021_NPC_KO_vs_NPC_WT_q05F05.pdf", width=12, height=6) # qvalue at 1.1!

pdf("output/GO_hg38/enrichR_KEGG_2021_NPC_KO_vs_NPC_WT_q05F05_logAdjP2.6.pdf", width=12, height=6) # qvalue at 1.1!

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 9, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 25)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_KEGG_2021_NPC_KO_vs_NPC_WT_q05F05.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("GO_Biological_Process_2018")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2018
down <- edown$GO_Biological_Process_2018
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_GO_BP_8wN_KO_vs_8wN_WT.pdf", width=12, height=7)
ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()




# Define databases for enrichment
dbs <- c("GO_Biological_Process_2018")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2018
down <- edown$GO_Biological_Process_2018
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_GO_BP_8wN_KO_vs_8wN_WT.pdf", width=12, height=7)
ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()







# Define databases for enrichment
dbs <- c("RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO'
down <- edown$'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO'
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("GSE[0-9]+$", "", gos$Term)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Disease_Gene_and_Drug_Signatures_from_GEO_8wN_KO_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_Disease_Gene_and_Drug_Signatures_from_GEO_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Disease_Gene_and_Drug_Signatures_from_GEO_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Disease_Gene_and_Drug_Signatures_from_GEO_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("Jensen_DISEASES")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Jensen_DISEASES
down <- edown$Jensen_DISEASES
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Jensen_DISEASES_8wN_KO_vs_8wN_WT.pdf", width=10, height=7)
pdf("output/GO_hg38/enrichR_Jensen_DISEASES_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_Jensen_DISEASES_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=10, height=2)
pdf("output/GO_hg38/enrichR_Jensen_DISEASES_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Jensen_DISEASES_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Jensen_DISEASES_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("OMIM_Disease")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$OMIM_Disease
down <- edown$OMIM_Disease
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Jensen_DISEASES_8wN_KO_vs_8wN_WT.pdf", width=10, height=7)
pdf("output/GO_hg38/enrichR_Jensen_DISEASES_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Jensen_DISEASES_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Jensen_DISEASES_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("MAGMA_Drugs_and_Diseases")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$MAGMA_Drugs_and_Diseases
down <- edown$MAGMA_Drugs_and_Diseases
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_MAGMA_Drugs_and_Diseases_8wN_KO_vs_8wN_WT.pdf", width=8, height=4)
pdf("output/GO_hg38/enrichR_MAGMA_Drugs_and_Diseases_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_MAGMA_Drugs_and_Diseases_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_MAGMA_Drugs_and_Diseases_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("Rare_Diseases_AutoRIF_Gene_Lists")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Rare_Diseases_AutoRIF_Gene_Lists
down <- edown$Rare_Diseases_AutoRIF_Gene_Lists
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_KO_vs_8wN_WT.pdf", width=8, height=6)
pdf("output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=8, height=5)
pdf("output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=6)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("UK_Biobank_GWAS_v1")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$UK_Biobank_GWAS_v1
down <- edown$UK_Biobank_GWAS_v1
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_KO_vs_8wN_WT.pdf", width=8, height=4)
pdf("output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_MAGMA_Drugs_and_Diseases_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("ClinVar_2019")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$ClinVar_2019
down <- edown$ClinVar_2019
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_KO_vs_8wN_WT.pdf", width=8, height=4)
pdf("output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Rare_Diseases_AutoRIF_Gene_Lists_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_MAGMA_Drugs_and_Diseases_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("dbGaP")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$dbGaP
down <- edown$dbGaP
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_dbGaP_8wN_KO_vs_8wN_WT.pdf", width=8, height=2)
pdf("output/GO_hg38/enrichR_dbGaP_8wN_HET_vs_8wN_WT.pdf", width=14, height=2)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_dbGaP_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_dbGaP_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("DisGeNET")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt

### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$DisGeNET
down <- edown$DisGeNET
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 10)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 10)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_DisGeNET_8wN_KO_vs_8wN_WT.pdf", width=8, height=4)
pdf("output/GO_hg38/enrichR_DisGeNET_8wN_HET_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_DisGeNET_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=8, height=3)
pdf("output/GO_hg38/enrichR_DisGeNET_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=3)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_DisGeNET_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_DisGeNET_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("Elsevier_Pathway_Collection")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt

### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Elsevier_Pathway_Collection
down <- edown$Elsevier_Pathway_Collection
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Elsevier_Pathway_Collection_8wN_KO_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_Elsevier_Pathway_Collection_8wN_HET_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_Elsevier_Pathway_Collection_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_Elsevier_Pathway_Collection_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=6)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Elsevier_Pathway_Collection_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Elsevier_Pathway_Collection_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("ENCODE_Histone_Modifications_2015")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$ENCODE_Histone_Modifications_2015
down <- edown$ENCODE_Histone_Modifications_2015
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)




# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_ENCODE_Histone_Modifications_2015_8wN_KO_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_ENCODE_Histone_Modifications_2015_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_ENCODE_Histone_Modifications_2015_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_ENCODE_Histone_Modifications_2015_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("ENCODE_TF_ChIP-seq_2015")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$`ENCODE_TF_ChIP-seq_2015`
down <- edown$`ENCODE_TF_ChIP-seq_2015`
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_ENCODE_TF_ChIP-seq_2015_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_ENCODE_TF_ChIP-seq_2015_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_ENCODE_TF_ChIP-seq_2015_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_ENCODE_TF_ChIP-seq_2015_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)









# Define databases for enrichment
dbs <- c("ESCAPE")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$ESCAPE
down <- edown$ESCAPE
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_ESCAPE_8wN_KO_vs_8wN_WT.pdf", width=14, height=8)
pdf("output/GO_hg38/enrichR_ESCAPE_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_ESCAPE_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_ESCAPE_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment
dbs <- c("GeDiPNet_2023")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GeDiPNet_2023
down <- edown$GeDiPNet_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_GeDiPNet_2023_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_GeDiPNet_2023_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)
pdf("output/GO_hg38/enrichR_GeDiPNet_2023_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=3)
pdf("output/GO_hg38/enrichR_GeDiPNet_2023_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=5)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_GeDiPNet_2023_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_GeDiPNet_2023_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)









# Define databases for enrichment
dbs <- c("Genome_Browser_PWMs")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Genome_Browser_PWMs
down <- edown$Genome_Browser_PWMs
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Genome_Browser_PWMs_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_Genome_Browser_PWMs_8wN_HET_vs_8wN_WT.pdf", width=14, height=4)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Genome_Browser_PWMs_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Genome_Browser_PWMs_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023")

### Gene list symbol for qval 0.05 and FC 0.5 _ 8wN
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5 _ 8wN
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt


### Gene list symbol for qval 0.05 and FC 0.5 _ NPC
output/deseq2_hg38/downregulated_q05FC05_NPC_KO_vs_NPC_WT.txt
output/deseq2_hg38/upregulated_q05FC05_NPC_KO_vs_NPC_WT.txt
output/deseq2_hg38/downregulated_q05FC05_NPC_HET_vs_NPC_WT.txt
output/deseq2_hg38/upregulated_q05FC05_NPC_HET_vs_NPC_WT.txt
xxxxxxxxxxxxxxxxxxxx


# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q05FC05_NPC_HET_vs_NPC_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q05FC05_NPC_HET_vs_NPC_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_KO_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=8)

pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_NPC_KO_vs_NPC_WT_q05FC05.pdf", width=14, height=8)
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_NPC_HET_vs_NPC_WT_q05FC05.pdf", width=14, height=8)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_GO_Biological_Process_2023_NPC_HET_vs_NPC_WT_q05FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)








# Define databases for enrichment
dbs <- c("GO_Cellular_Component_2023")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Cellular_Component_2023
down <- edown$GO_Cellular_Component_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_GO_Cellular_Component_2023_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_GO_Cellular_Component_2023_8wN_HET_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_GO_Cellular_Component_2023_8wN_KO_vs_8wN_WT_q01Fc05.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_GO_Cellular_Component_2023_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_GO_Cellular_Component_2023_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_GO_Cellular_Component_2023_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("GO_Molecular_Function_2023")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Molecular_Function_2023
down <- edown$GO_Molecular_Function_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_GO_Molecular_Function_2023_8wN_KO_vs_8wN_WT.pdf", width=14, height=4)
pdf("output/GO_hg38/enrichR_GO_Molecular_Function_2023_8wN_HET_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_GO_Molecular_Function_2023_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=4)
pdf("output/GO_hg38/enrichR_GO_Molecular_Function_2023_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=4)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_GO_Molecular_Function_2023_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_GO_Molecular_Function_2023_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment
dbs <- c("GWAS_Catalog_2023")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GWAS_Catalog_2023
down <- edown$GWAS_Catalog_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_GWAS_Catalog_2023_8wN_KO_vs_8wN_WT.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_GWAS_Catalog_2023_8wN_HET_vs_8wN_WT.pdf", width=14, height=3)
pdf("output/GO_hg38/enrichR_GWAS_Catalog_2023_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=3)
pdf("output/GO_hg38/enrichR_GWAS_Catalog_2023_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=2)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_GWAS_Catalog_2023_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_GWAS_Catalog_2023_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("BioPlanet_2019")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt

### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$BioPlanet_2019
down <- edown$BioPlanet_2019
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_BioPlanet_2019_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_BioPlanet_2019_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)
pdf("output/GO_hg38/enrichR_BioPlanet_2019_8wN_KO_vs_8wN_WT_q01Fc05.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_BioPlanet_2019_8wN_HET_vs_8wN_WT_q01Fc05.pdf", width=14, height=7)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_BioPlanet_2019_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_BioPlanet_2019_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)








# Define databases for enrichment
dbs <- c("HMDB_Metabolites")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$HMDB_Metabolites
down <- edown$HMDB_Metabolites
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_HMDB_Metabolites_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_HMDB_Metabolites_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_HMDB_Metabolites_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_HMDB_Metabolites_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("Human_Phenotype_Ontology")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Human_Phenotype_Ontology
down <- edown$Human_Phenotype_Ontology
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(HP:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Human_Phenotype_Ontology_8wN_KO_vs_8wN_WT.pdf", width=14, height=4)
pdf("output/GO_hg38/enrichR_Human_Phenotype_Ontology_8wN_HET_vs_8wN_WT.pdf", width=14, height=4)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Human_Phenotype_Ontology_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Human_Phenotype_Ontology_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("HumanCyc_2016")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$HumanCyc_2016
down <- edown$HumanCyc_2016
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_HumanCyc_2016_8wN_KO_vs_8wN_WT.pdf", width=14, height=4)
pdf("output/GO_hg38/enrichR_HumanCyc_2016_8wN_HET_vs_8wN_WT.pdf", width=14, height=4)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_HumanCyc_2016_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_HumanCyc_2016_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment
dbs <- c("MGI_Mammalian_Phenotype_Level_4_2021")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$MGI_Mammalian_Phenotype_Level_4_2021
down <- edown$MGI_Mammalian_Phenotype_Level_4_2021
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("MP:[0-9]+", "", gos$Term)  # Regular expression to match the pattern "MP: [number]" and replace it with an empty

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_MGI_Mammalian_Phenotype_Level_4_2021_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_MGI_Mammalian_Phenotype_Level_4_2021_8wN_HET_vs_8wN_WT.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_MGI_Mammalian_Phenotype_Level_4_2021_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_MGI_Mammalian_Phenotype_Level_4_2021_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_MGI_Mammalian_Phenotype_Level_4_2021_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_MGI_Mammalian_Phenotype_Level_4_2021_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)



# Define databases for enrichment
dbs <- c("KOMP2_Mouse_Phenotypes_2022")


### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$KOMP2_Mouse_Phenotypes_2022
down <- edown$KOMP2_Mouse_Phenotypes_2022
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
## NO ENRICHMENT


# Define databases for enrichment
dbs <- c("OMIM_Expanded")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$OMIM_Expanded
down <- edown$OMIM_Expanded
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_OMIM_Expanded_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_OMIM_Expanded_8wN_HET_vs_8wN_WT.pdf", width=14, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_OMIM_Expanded_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_OMIM_Expanded_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("Orphanet_Augmented_2021")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Orphanet_Augmented_2021
down <- edown$Orphanet_Augmented_2021
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("ORPHA:[0-9]+", "", gos$Term)  # Regular expression to match the pattern "MP: [number]" and replace it with an empty

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Orphanet_Augmented_2021_8wN_KO_vs_8wN_WT.pdf", width=14, height=9)
pdf("output/GO_hg38/enrichR_Orphanet_Augmented_2021_8wN_HET_vs_8wN_WT.pdf", width=14, height=10)
pdf("output/GO_hg38/enrichR_Orphanet_Augmented_2021_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=9)
pdf("output/GO_hg38/enrichR_Orphanet_Augmented_2021_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=10)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Orphanet_Augmented_2021_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Orphanet_Augmented_2021_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("PFOCR_Pathways_2023")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$PFOCR_Pathways_2023
down <- edown$PFOCR_Pathways_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_PFOCR_Pathways_2023_8wN_KO_vs_8wN_WT.pdf", width=14, height=9)
pdf("output/GO_hg38/enrichR_PFOCR_Pathways_2023_8wN_HET_vs_8wN_WT.pdf", width=14, height=10)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_PFOCR_Pathways_2023_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_PFOCR_Pathways_2023_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)







# Define databases for enrichment
dbs <- c("PheWeb_2019")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$PheWeb_2019
down <- edown$PheWeb_2019
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_PheWeb_2019_8wN_KO_vs_8wN_WT.pdf", width=14, height=2)
pdf("output/GO_hg38/enrichR_PheWeb_2019_8wN_HET_vs_8wN_WT.pdf", width=14, height=10)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_PheWeb_2019_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_PheWeb_2019_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)



# Define databases for enrichment
dbs <- c("Rare_Diseases_GeneRIF_Gene_Lists")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Rare_Diseases_GeneRIF_Gene_Lists
down <- edown$Rare_Diseases_GeneRIF_Gene_Lists
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Rare_Diseases_GeneRIF_Gene_Lists_8wN_KO_vs_8wN_WT.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_Rare_Diseases_GeneRIF_Gene_Lists_8wN_HET_vs_8wN_WT.pdf", width=14, height=10)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Rare_Diseases_GeneRIF_Gene_Lists_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Rare_Diseases_GeneRIF_Gene_Lists_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("SynGO_2022")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$SynGO_2022
down <- edown$SynGO_2022
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_SynGO_2022_8wN_KO_vs_8wN_WT.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_SynGO_2022_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_SynGO_2022_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_SynGO_2022_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("TG_GATES_2020")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$TG_GATES_2020
down <- edown$TG_GATES_2020
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_TG_GATES_2020_8wN_KO_vs_8wN_WT.pdf", width=14, height=8)
pdf("output/GO_hg38/enrichR_TG_GATES_2020_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_TG_GATES_2020_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_TG_GATES_2020_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("TRRUST_Transcription_Factors_2019")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$TRRUST_Transcription_Factors_2019
down <- edown$TRRUST_Transcription_Factors_2019
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_TRRUST_Transcription_Factors_2019_8wN_KO_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_TRRUST_Transcription_Factors_2019_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_TRRUST_Transcription_Factors_2019_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_TRRUST_Transcription_Factors_2019_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("WikiPathway_2021_Human")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt

### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$WikiPathway_2021_Human
down <- edown$WikiPathway_2021_Human
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("WP[0-9]+", "", gos$Term)  # Regular expression to match the pattern "MP: [number]" and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_WikiPathway_2021_Human_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=6)
pdf("output/GO_hg38/enrichR_WikiPathway_2021_Human_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)
pdf("output/GO_hg38/enrichR_WikiPathway_2021_Human_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_WikiPathway_2021_Human_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_WikiPathway_2021_Human_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)



# Define databases for enrichment
dbs <- c("TRRUST_Transcription_Factors_2019")

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$TRRUST_Transcription_Factors_2019
down <- edown$TRRUST_Transcription_Factors_2019
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_TRRUST_Transcription_Factors_2019_8wN_KO_vs_8wN_WT.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_TRRUST_Transcription_Factors_2019_8wN_HET_vs_8wN_WT.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_TRRUST_Transcription_Factors_2019_8wN_KO_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_TRRUST_Transcription_Factors_2019_8wN_HET_vs_8wN_WT.txt", sep="\t", row.names=FALSE, quote=FALSE)



# Define databases for enrichment
dbs <- c("Descartes_Cell_Types_and_Tissue_2021")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Descartes_Cell_Types_and_Tissue_2021
down <- edown$Descartes_Cell_Types_and_Tissue_2021
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Descartes_Cell_Types_and_Tissue_2021_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=5)
pdf("output/GO_hg38/enrichR_Descartes_Cell_Types_and_Tissue_2021_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Descartes_Cell_Types_and_Tissue_2021_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Descartes_Cell_Types_and_Tissue_2021_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("CellMarker_Augmented_2021")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$CellMarker_Augmented_2021
down <- edown$CellMarker_Augmented_2021
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_CellMarker_Augmented_2021_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_CellMarker_Augmented_2021_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=9)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_CellMarker_Augmented_2021_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_CellMarker_Augmented_2021_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("Azimuth_Cell_Types_2021")

### Gene list symbol for qval 0.05 and FC 0.5
output/deseq2_hg38/downregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05_8wN_HET_vs_8wN_WT.txt
### Gene list symbol for qval 0.01 and FC 0.5
output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt

# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_KO_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q01FC05_8wN_HET_vs_8wN_WT.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$Azimuth_Cell_Types_2021
down <- edown$Azimuth_Cell_Types_2021
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("CL[0-9]+", "", gos$Term)  # Regular expression to match the pattern "MP: [number]" and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO_hg38/enrichR_Azimuth_Cell_Types_2021_8wN_KO_vs_8wN_WT_q01FC05.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_Azimuth_Cell_Types_2021_8wN_HET_vs_8wN_WT_q01FC05.pdf", width=14, height=6)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_Azimuth_Cell_Types_2021_8wN_KO_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_Azimuth_Cell_Types_2021_8wN_HET_vs_8wN_WT_q01FC05.txt", sep="\t", row.names=FALSE, quote=FALSE)



# For NaiaraPlot (20240119)

# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023")

### Gene list symbol for qval 0.05 and FC 1
output/deseq2_hg38/downregulated_q05FC1_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05FC1_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05FC1_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05FC1_8wN_KO_vs_8wN_WT.txt

### Gene list symbol for qval 0.05 and FC 2
output/deseq2_hg38/downregulated_q05FC2_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05FC2_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05FC2_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05FC2_8wN_KO_vs_8wN_WT.txt

### Gene list symbol for RNAseq CutRun integration 20240119
output/deseq2_hg38/downregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_gain.txt
output/deseq2_hg38/upregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_lost.txt

output/deseq2_hg38/downregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_gain.txt
output/deseq2_hg38/upregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_lost.txt



# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_gain.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_lost.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics 
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_HET_vs_8wN_WT_q05FC1.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_KO_vs_8wN_WT_q05FC1.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_HET_vs_8wN_WT_q05FC2.pdf", width=14, height=7) #no enricghment
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_KO_vs_8wN_WT_q05FC2.pdf", width=14, height=7) 

pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot.pdf", width=10, height=3) 
pdf("output/GO_hg38/enrichR_GO_Biological_Process_2023_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot.pdf", width=10, height=3) 

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_GO_Biological_Process_2023_8wN_KO_vs_8wN_WT_q05FC2.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_GO_Biological_Process_2023_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment
dbs <- c("KEGG_2021_Human")

### Gene list symbol for qval 0.05 and FC 1
output/deseq2_hg38/downregulated_q05FC1_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05FC1_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05FC1_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05FC1_8wN_KO_vs_8wN_WT.txt

### Gene list symbol for qval 0.05 and FC 2
output/deseq2_hg38/downregulated_q05FC2_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05FC2_8wN_HET_vs_8wN_WT.txt
output/deseq2_hg38/downregulated_q05FC2_8wN_KO_vs_8wN_WT.txt
output/deseq2_hg38/upregulated_q05FC2_8wN_KO_vs_8wN_WT.txt

### Gene list symbol for RNAseq CutRun integration 20240119
output/deseq2_hg38/downregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_gain.txt
output/deseq2_hg38/upregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_lost.txt

output/deseq2_hg38/downregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_gain.txt
output/deseq2_hg38/upregulated_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot_lost.txt




# Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2_hg38/downregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_gain.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

# Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2_hg38/upregulated_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot_lost.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$KEGG_2021_Human
down <- edown$KEGG_2021_Human
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics 
pdf("output/GO_hg38/enrichR_KEGG_2021_Human_8wN_HET_vs_8wN_WT_q05FC1.pdf", width=14, height=7)
pdf("output/GO_hg38/enrichR_KEGG_2021_Human_8wN_KO_vs_8wN_WT_q05FC1.pdf", width=14, height=3)
pdf("output/GO_hg38/enrichR_KEGG_2021_Human_8wN_KO_vs_8wN_WT_q05FC2.pdf", width=14, height=3)


pdf("output/GO_hg38/enrichR_KEGG_2021_Human_q05FC05_8wN_KO_vs_8wN_WT_20240119_WTvsKO_annot.pdf", width=10, height=6)
pdf("output/GO_hg38/enrichR_KEGG_2021_Human_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot.pdf", width=10, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()
## save output
write.table(gos, "output/GO_hg38/enrichR_KEGG_2021_Human_8wN_KO_vs_8wN_WT_q05FC1.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO_hg38/enrichR_KEGG_2021_Human_q05FC05_8wN_HET_vs_8wN_WT_20240119_WTvsHET_annot.txt", sep="\t", row.names=FALSE, quote=FALSE)







```






# RNAseq deconvolution

Let's try various deconvolution method to predict cell type population:


- [BrainDeconvShiny](https://voineagulab.shinyapps.io/BrainDeconvShiny/)

## BrainDeconvShiny

For this one; let's just generate a table of tpm-level gene expression (they said rpkm; but let's try with tpm...) with with rownames as Ensembl gene IDs and colnames as sampleIDs



```R
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

samples <- c("8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4",
                     "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4",
                     "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4")



## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/tpm_hg38/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR_hg38.")) %>%
    rename(!!sample := starts_with("output.STAR_hg38."))
}

## Merge all dataframe into a single one
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")

tpm_all_sample_tidy <- tpm_all_sample %>%
  gather(key = 'variable', value = 'tpm', -Geneid) %>%
  separate(variable, into = c('time', 'genotype', 'replicate'), sep = "_") %>%
  rename(gene = Geneid)

tpm_all_sample_tidy$gene <- gsub("\\..*", "", tpm_all_sample_tidy$gene) # remove Ensembl gene id version

## Calculate median for each sample
tpm_all_sample_tidy_median = tpm_all_sample_tidy %>%
  group_by(gene, genotype) %>%
  summarise(median = median(tpm))

## save
write.table(tpm_all_sample_tidy_median, file = "output/tpm_hg38/tpm_all_sample_tidy_median.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## pivot table for BrainDeconvShiny
tpm_all_sample_tidy_median_BrainDeconvShiny = tpm_all_sample_tidy_median %>%
  spread(key = genotype, value = median)

write.table(tpm_all_sample_tidy_median_BrainDeconvShiny, file = "output/tpm_hg38/tpm_all_sample_tidy_median_BrainDeconvShiny.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Save sample per sample
tpm_all_sample_tidy_8wN_KO_BrainDeconvShiny = tpm_all_sample_tidy %>%
  filter(time == "8wN", genotype == "KO") %>%
  dplyr::select(gene,replicate,tpm) %>%
  group_by(gene, replicate) %>%
  summarise(tpm = mean(tpm, na.rm = TRUE)) %>%  # some gene id are dupplicated as we remove version id; need do mean...
  spread(key = replicate, value = tpm)

write.table(tpm_all_sample_tidy_8wN_KO_BrainDeconvShiny, file = "output/tpm_hg38/tpm_all_sample_tidy_8wN_KO_BrainDeconvShiny.txt", sep = "\t", quote = FALSE, row.names = FALSE)


```


Let's generate rpkm too in case...



```R
#### collect all samples ID
samples <- c("8wN_WT_R1", "8wN_WT_R2", "8wN_WT_R3", "8wN_WT_R4",
                     "8wN_HET_R1", "8wN_HET_R2", "8wN_HET_R3", "8wN_HET_R4",
                     "8wN_KO_R1", "8wN_KO_R2", "8wN_KO_R3", "8wN_KO_R4")
samples <- c("4wN_WT_R1", "4wN_WT_R2",
            "4wN_HET_R3", "4wN_HET_R4",
            "4wN_KO_R1", "4wN_KO_R2")
samples <- c("NPC_WT_R1", "NPC_WT_R2", "NPC_WT_R3",
            "NPC_HET_R1", "NPC_HET_R2", "NPC_HET_R3",
            "NPC_KO_R1", "NPC_KO_R2", "NPC_KO_R3")
samples <- c("2dN_WT_R1", "2dN_WT_R2", "2dN_WT_R3",
            "2dN_HET_R1", "2dN_HET_R2", "2dN_HET_R3",
            "2dN_KO_R1", "2dN_KO_R2", "2dN_KO_R3")
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
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/rpkm_hg38/", sample, "_rpkm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR_hg38.")) %>%
    rename(!!sample := starts_with("output.STAR_hg38."))
}


## Merge all dataframe into a single one
rpkm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")
## rpkm_all_sample$Geneid <- gsub("\\..*", "", rpkm_all_sample$Geneid) # remove Ensembl gene id version
## write.csv(rpkm_all_sample %>% unique(), file="output/rpkm_hg38/rpkm_all_sample.txt")



rpkm_all_sample_tidy <- rpkm_all_sample %>%
  gather(key = 'variable', value = 'rpkm', -Geneid) %>%
  separate(variable, into = c('time', 'genotype', 'replicate'), sep = "_") %>%
  rename(gene = Geneid)

rpkm_all_sample_tidy$gene <- gsub("\\..*", "", rpkm_all_sample_tidy$gene) # remove Ensembl gene id version

## Calculate median for each sample
rpkm_all_sample_tidy_median = rpkm_all_sample_tidy %>%
  group_by(gene, genotype) %>%
  summarise(median = median(rpkm))
## save
write.table(rpkm_all_sample_tidy_median, file = "output/rpkm_hg38/rpkm_all_sample_tidy_median.txt", sep = "\t", quote = FALSE, row.names = FALSE)
## pivot table for BrainDeconvShiny
rpkm_all_sample_tidy_median_BrainDeconvShiny = rpkm_all_sample_tidy_median %>%
  spread(key = genotype, value = median)
write.table(rpkm_all_sample_tidy_median_BrainDeconvShiny, file = "output/rpkm_hg38/rpkm_all_sample_tidy_median_BrainDeconvShiny.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## Save sample per sample
rpkm_all_sample_tidy_8wN_KO_BrainDeconvShiny = rpkm_all_sample_tidy %>%
  filter(time == "2dN", genotype == "KO") %>%
  dplyr::select(gene,replicate,rpkm) %>%
  group_by(gene, replicate) %>%
  summarise(rpkm = mean(rpkm, na.rm = TRUE)) %>%  # some gene id are dupplicated as we remove version id; need do mean...
  spread(key = replicate, value = rpkm)
write.table(rpkm_all_sample_tidy_8wN_KO_BrainDeconvShiny, file = "output/rpkm_hg38/rpkm_all_sample_tidy_2dN_KO_BrainDeconvShiny.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```


--> Export file to Local and convert to csv (to avoid `Not a matrix error`, re-order the columns WT, then HET, then KO)
----> CA last ~5min to run!

*NOTE: [This](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02290-6) paper recommend to use TPM as normalized per library size! [This](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03016-6) paper mention other normalization method worth to be tested (logNormlaize, CPM, lognormnalize,...)*


--> I tested genotype per genotype but the r goodness of  fit remains the same...

--> tpm sometime perform better! Quite random; but for now db MM seems to be the best



# GSEA

Objective here is to do GSEA to **1.** assess cell type population changes (notably the changes of Astrocyte and neuronal population, as it change with the EZH1 mutations according to deconvolution analysis)

- Let's follow this [tutorial](https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html).
- Look here interesting cell type markers [here](https://www.gsea-msigdb.org/gsea/msigdb/human/)
  - DESCARTES_FETAL: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/DESCARTES_MAIN_FETAL_ENS_NEURONS.html 
  - FAN_EMBRYONIC_CTX https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/FAN_EMBRYONIC_CTX_BIG_GROUPS_EXCITATORY_NEURON.html 
  - ZHONG_PFC : https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/ZHONG_PFC_C1_DLX5_POS_INTERNEURON.html 


Let's also do GSEA to **2.** check steroids-related genes at each time points for WT vs HET/GOF

Use `conda activate deseq2` and R:



```R
# Packages
library("tidyverse")
library("clusterProfiler")
library("msigdbr") # BiocManager::install("msigdbr")
library("org.Mm.eg.db")
library("enrichplot") # for gseaplot2()
library("pheatmap")

# import DEGs
## filtered_8wN_KO_vs_8wN_WT
## filtered_8wN_HET_vs_8wN_WT

KO <- read.table("output/deseq2_hg38/filtered_8wN_KO_vs_8wN_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() 
KO_geneSymbol = KO %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene

HET <- read.table("output/deseq2_hg38/filtered_2dN_HET_vs_2dN_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()
HET_geneSymbol = HET %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene
  

  



# import msigdbr cell marker db 
hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "C8"   # From C8 cell marker category
)
hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "C2"   # From C2 pathways
)

# Order our DEG
## Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- KO_geneSymbol$log2FoldChange
names(lfc_vector) <- KO_geneSymbol$GeneSymbol
## We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
### Set the seed so our results are reproducible:
set.seed(42)


# run GSEA

## without pvalue cutoff
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    hs_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

gsea_result_df <- data.frame(gsea_results@result)



# ASTROCYTE -----------------
## From this file I look for significance in 'Astrocyte' and it return:
## For KO
DESCARTES_MAIN_FETAL_ASTROCYTES
DESCARTES_FETAL_CEREBELLUM_ASTROCYTES
DESCARTES_FETAL_CEREBRUM_ASTROCYTES
ZHONG_PFC_MAJOR_TYPES_ASTROCYTES
ZHONG_PFC_C3_ASTROCYTE
FAN_EMBRYONIC_CTX_ASTROCYTE_1
FAN_EMBRYONIC_CTX_ASTROCYTE_2
## For HET
ZHONG_PFC_MAJOR_TYPES_ASTROCYTES
FAN_EMBRYONIC_CTX_ASTROCYTE_2
FAN_EMBRYONIC_CTX_ASTROCYTE_1
ZHONG_PFC_C1_ASTROCYTE
ZHONG_PFC_C3_ASTROCYTE
DESCARTES_FETAL_CEREBELLUM_ASTROCYTES
DESCARTES_FETAL_CEREBRUM_ASTROCYTES
DESCARTES_FETAL_EYE_ASTROCYTES


# plots
DESCARTES_MAIN_FETAL_ASTROCYTES
pdf("output/gsea/FAN_EMBRYONIC_CTX_ASTROCYTE_2_8wN_WTvsKO.pdf", width=14, height=8)

pdf("output/gsea/DESCARTES_FETAL_EYE_ASTROCYTES_8wN_WTvsHET.pdf", width=14, height=8)
enrichplot::gseaplot(
  gsea_results,
  geneSetID = "DESCARTES_FETAL_EYE_ASTROCYTES",
  title = "DESCARTES_FETAL_EYE_ASTROCYTES",
  color.line = "#0d76ff"
)
dev.off()

## gseaplot2() multiple gene lists
DESCARTES_FETAL_CEREBRUM_EXCITATORY_NEURONS 
DESCARTES_FETAL_CEREBRUM_INHIBITORY_NEURONS 
DESCARTES_FETAL_CEREBRUM_ASTROCYTES 
DESCARTES_FETAL_CEREBRUM_OLIGODENDROCYTES 
DESCARTES_FETAL_CEREBRUM_MICROGLIA 

## Find corresponding gene ID number
gsea_results$ID # TO CHECK WHICH NB (=gene lists genes) TO USE
which(gsea_results$ID == "DESCARTES_FETAL_CEREBRUM_MICROGLIA")

pdf("output/gsea/DESCARTES_FETAL_CEREBRUM_8wN_WTvsKO.pdf", width=14, height=15)
gseaplot2(gsea_results, geneSetID = c(694,84,82,401,124),  pvalue_table = TRUE)
dev.off()

## heatmap Enrichment Score
gsea_result_df_KO = as_tibble(gsea_result_df) %>%
  dplyr::select(ID, enrichmentScore, qvalue) %>%
  add_column(genotype = "KO")
gsea_result_df_HET = as_tibble(gsea_result_df) %>%
  dplyr::select(ID, enrichmentScore, qvalue) %>%
  add_column(genotype = "HET")
gsea_result_df_tidy = gsea_result_df_KO %>%
  bind_rows(gsea_result_df_HET)

### Set up the heatmap
desired_ids <- c(
"DESCARTES_FETAL_CEREBRUM_ASTROCYTES",
"DESCARTES_FETAL_CEREBRUM_EXCITATORY_NEURONS",
"DESCARTES_FETAL_CEREBRUM_INHIBITORY_NEURONS",
"DESCARTES_FETAL_CEREBRUM_LIMBIC_SYSTEM_NEURONS",
"DESCARTES_FETAL_CEREBRUM_MICROGLIA",
"DESCARTES_FETAL_CEREBRUM_OLIGODENDROCYTES",
"DESCARTES_FETAL_CEREBRUM_SKOR2_NPSR1_POSITIVE_CELLS",
"DESCARTES_FETAL_CEREBRUM_VASCULAR_ENDOTHELIAL_CELLS"
)
desired_ids <- c(
"DESCARTES_FETAL_CEREBELLUM_ASTROCYTES",
"DESCARTES_FETAL_CEREBELLUM_GRANULE_NEURONS",
"DESCARTES_FETAL_CEREBELLUM_INHIBITORY_INTERNEURONS",
"DESCARTES_FETAL_CEREBELLUM_MICROGLIA",
"DESCARTES_FETAL_CEREBELLUM_OLIGODENDROCYTES",
"DESCARTES_FETAL_CEREBELLUM_PURKINJE_NEURONS",
"DESCARTES_FETAL_CEREBELLUM_SLC24A4_PEX5L_POSITIVE_CELLS",
"DESCARTES_FETAL_CEREBELLUM_UNIPOLAR_BRUSH_CELLS",
"DESCARTES_FETAL_CEREBELLUM_VASCULAR_ENDOTHELIAL_CELLS"
)
desired_ids <- c(
  "ZHONG_PFC_MAJOR_TYPES_EXCITATORY_NEURON", 
  "ZHONG_PFC_MAJOR_TYPES_INTERNEURON", 
  "ZHONG_PFC_MAJOR_TYPES_ASTROCYTES", 
  "ZHONG_PFC_MAJOR_TYPES_OPC", 
  "ZHONG_PFC_MAJOR_TYPES_MICROGLIA",
  "ZHONG_PFC_MAJOR_TYPES_NPCS"
)
desired_ids <- c(
"FAN_EMBRYONIC_CTX_BIG_GROUPS_BRAIN_ENDOTHELIAL",
"FAN_EMBRYONIC_CTX_BIG_GROUPS_BRAIN_IMMUNE",
"FAN_EMBRYONIC_CTX_BIG_GROUPS_CAJAL_RETZIUS",
"FAN_EMBRYONIC_CTX_BIG_GROUPS_GLIAL",
"FAN_EMBRYONIC_CTX_BIG_GROUPS_INHIBITORY",
"FAN_EMBRYONIC_CTX_BIG_GROUPS_MICROGLIA",
"FAN_EMBRYONIC_CTX_ASTROCYTE_1",
"FAN_EMBRYONIC_CTX_ASTROCYTE_2",
"FAN_EMBRYONIC_CTX_BRAIN_B_CELL",
"FAN_EMBRYONIC_CTX_BRAIN_EFFECTOR_T_CELL",
"FAN_EMBRYONIC_CTX_BRAIN_ENDOTHELIAL_1",
"FAN_EMBRYONIC_CTX_BRAIN_ENDOTHELIAL_2",
"FAN_EMBRYONIC_CTX_BRAIN_MYELOID",
"FAN_EMBRYONIC_CTX_BRAIN_NAIVE_LIKE_T_CELL",
"FAN_EMBRYONIC_CTX_EX_1_EXCITATORY_NEURON",
"FAN_EMBRYONIC_CTX_EX_2_EXCITATORY_NEURON",
"FAN_EMBRYONIC_CTX_EX_4_EXCITATORY_NEURON",
"FAN_EMBRYONIC_CTX_IN_1_INTERNEURON",
"FAN_EMBRYONIC_CTX_IN_2_INTERNEURON",
"FAN_EMBRYONIC_CTX_IN_3_INTERNEURON",
"FAN_EMBRYONIC_CTX_IN_4_INTERNEURON",
"FAN_EMBRYONIC_CTX_IN_5_INTERNEURON",
"FAN_EMBRYONIC_CTX_IN_6_INTERNEURON",
"FAN_EMBRYONIC_CTX_MICROGLIA_1",
"FAN_EMBRYONIC_CTX_MICROGLIA_2",
"FAN_EMBRYONIC_CTX_MICROGLIA_3",
"FAN_EMBRYONIC_CTX_NSC_1",
"FAN_EMBRYONIC_CTX_NSC_2",
"FAN_EMBRYONIC_CTX_OLIG",
"FAN_EMBRYONIC_CTX_OPC" 
)



# Filter the data for desired IDs
filtered_data <- gsea_result_df_tidy %>%
  filter(ID %in% desired_ids)

pdf("output/gsea/heatmap_DESCARTES_FETAL_CEREBRUM_8wN.pdf", width=1.5, height=5)
pdf("output/gsea/heatmap_ZHONG_PFC_MAJOR_TYPES_8wN.pdf", width=1.25, height=5)
pdf("output/gsea/heatmap_FAN_EMBRYONIC_CTX_8wN.pdf", width=4, height=4)
pdf("output/gsea/heatmap_DESCARTES_FETAL_CEREBELLUM_8wN.pdf", width=1.5, height=5)

ggplot(filtered_data, aes(x=ID, y=genotype, fill=enrichmentScore)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Enrichment\nScore") +
  geom_text(aes(label=sprintf("%.2f", enrichmentScore)), 
            color = ifelse(filtered_data$qvalue <= 0.05, "black", "grey50"), 
            size=1) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()



# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_8wN_WTvsHET_complete.tsv"
  )
)

## From this file I look for significance in 'neurons' and it return:
## For KO
FAN_EMBRYONIC_CTX_EX_4_EXCITATORY_NEURON
FAN_EMBRYONIC_CTX_BIG_GROUPS_EXCITATORY_NEURON
ZHONG_PFC_C7_ORG_UNDERGOING_NEURONAL_DIFFERENTIATION

DESCARTES_FETAL_LUNG_VISCERAL_NEURONS
DESCARTES_FETAL_CEREBELLUM_GRANULE_NEURONS
DESCARTES_FETAL_STOMACH_ENS_NEURONS

DESCARTES_FETAL_INTESTINE_ENS_NEURONS
DESCARTES_FETAL_PANCREAS_ENS_NEURONS
DESCARTES_FETAL_CEREBRUM_INHIBITORY_NEURONS
ZHONG_PFC_C8_UNKNOWN_NEUROD2_POS_INTERNEURON

## For HET
DESCARTES_FETAL_HEART_VISCERAL_NEURONS

DESCARTES_FETAL_INTESTINE_ENS_NEURONS
DESCARTES_FETAL_PANCREAS_ENS_NEURONS
DESCARTES_FETAL_CEREBRUM_INHIBITORY_NEURONS
ZHONG_PFC_C8_UNKNOWN_NEUROD2_POS_INTERNEURON





## From this file I look for significance in 'dendrocyte' and it return:
## For KO
DESCARTES_MAIN_FETAL_OLIGODENDROCYTES


## For HET
DESCARTES_FETAL_CEREBELLUM_OLIGODENDROCYTES



## From this file I look for significance in 'microglia' and it return:
## For KO
FAN_EMBRYONIC_CTX_MICROGLIA_1
ZHONG_PFC_C3_MICROGLIA
ZHONG_PFC_C1_MICROGLIA
ZHONG_PFC_MAJOR_TYPES_MICROGLIA
DESCARTES_FETAL_CEREBELLUM_MICROGLIA
DESCARTES_FETAL_CEREBRUM_MICROGLIA
HU_FETAL_RETINA_MICROGLIA

## For HET
FAN_EMBRYONIC_CTX_MICROGLIA_1
ZHONG_PFC_C1_MICROGLIA
DESCARTES_FETAL_CEREBELLUM_MICROGLIA
DESCARTES_FETAL_CEREBRUM_MICROGLIA
FAN_EMBRYONIC_CTX_BIG_GROUPS_MICROGLIA
ZHONG_PFC_MAJOR_TYPES_MICROGLIA
DESCARTES_FETAL_EYE_MICROGLIA
HU_FETAL_RETINA_MICROGLIA



## From this file I look for significance in 'endothelia' and it return:
## For KO
FAN_EMBRYONIC_CTX_BRAIN_ENDOTHELIAL_2
DESCARTES_FETAL_CEREBRUM_VASCULAR_ENDOTHELIAL_CELLS
MENON_FETAL_KIDNEY_9_ENDOTHELIAL_CELLS
LAKE_ADULT_KIDNEY_C22_ENDOTHELIAL_CELLS_GLOMERULAR_CAPILLARIES
GAUTAM_EYE_IRIS_CILIARY_BODY_CILIARY_BODY_ENDOTHELIAL_CELLS
DESCARTES_FETAL_EYE_VASCULAR_ENDOTHELIAL_CELLS
LAKE_ADULT_KIDNEY_C25_ENDOTHELIAL_CELLS_UNASSIGNED
CUI_DEVELOPING_HEART_C4_ENDOTHELIAL_CELL
CUI_DEVELOPING_HEART_VALVAR_ENDOTHELIAL_CELL
DESCARTES_FETAL_KIDNEY_VASCULAR_ENDOTHELIAL_CELLS
MURARO_PANCREAS_ENDOTHELIAL_CELL


## For HET
FAN_OVARY_CL7_ANGEIOGENIC_ENDOTHELIAL_CELL
FAN_EMBRYONIC_CTX_BIG_GROUPS_BRAIN_ENDOTHELIAL
LAKE_ADULT_KIDNEY_C23_ENDOTHELIAL_CELLS_AVR
FAN_OVARY_CL9_PUTATIVE_APOPTOTIC_ENDOTHELIAL_CELL
LAKE_ADULT_KIDNEY_C24_ENDOTHELIAL_CELLS_AEA_AND_DVR
RUBENSTEIN_SKELETAL_MUSCLE_PCV_ENDOTHELIAL_CELLS
FAN_EMBRYONIC_CTX_BRAIN_ENDOTHELIAL_2
MENON_FETAL_KIDNEY_9_ENDOTHELIAL_CELLS
MURARO_PANCREAS_ENDOTHELIAL_CELL
FAN_OVARY_CL16_LYMPHATIC_ENDOTHELIAL_CELL
FAN_EMBRYONIC_CTX_BRAIN_ENDOTHELIAL_1
RUBENSTEIN_SKELETAL_MUSCLE_ENDOTHELIAL_CELLS
CUI_DEVELOPING_HEART_VASCULAR_ENDOTHELIAL_CELL
DESCARTES_FETAL_CEREBRUM_VASCULAR_ENDOTHELIAL_CELLS
CUI_DEVELOPING_HEART_C4_ENDOTHELIAL_CELL
CUI_DEVELOPING_HEART_VALVAR_ENDOTHELIAL_CELL
LAKE_ADULT_KIDNEY_C22_ENDOTHELIAL_CELLS_GLOMERULAR_CAPILLARIES
DESCARTES_FETAL_KIDNEY_VASCULAR_ENDOTHELIAL_CELLS
GAUTAM_EYE_CHOROID_SCLERA_CHOROID_ENDOTHELIAL_CELLS
GAUTAM_EYE_IRIS_CILIARY_BODY_CILIARY_BODY_ENDOTHELIAL_CELLS
DESCARTES_FETAL_INTESTINE_VASCULAR_ENDOTHELIAL_CELLS
DESCARTES_MAIN_FETAL_LYMPHATIC_ENDOTHELIAL_CELLS




# STEROIDS -----------------
# GO THERE to find steroid pathway name : https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=C2

pdf("output/gsea/KEGG_STEROID_BIOSYNTHESIS__8wN_WTvsHET.pdf", width=14, height=8)
pdf("output/gsea/KEGG_STEROID_BIOSYNTHESIS__2dN_WTvsHET.pdf", width=14, height=8)

pdf("output/gsea/KEGG_STEROID_BIOSYNTHESIS__8wN_WTvsKO.pdf", width=14, height=8)

enrichplot::gseaplot(
  gsea_results,
  geneSetID = "KEGG_STEROID_BIOSYNTHESIS",
  title = "KEGG_STEROID_BIOSYNTHESIS",
  color.line = "#0d76ff"
)
dev.off()

# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_8wN_WTvsKO_complete_C2.tsv"
  )
)



```

--> I think useless to use the qvalue cutoff in the GSEA function; otherwise it will NOT return rows for non-significant pathway

--> `filtered*` for the DEGs means without X and Y chromosomes

--> **Astroctye enrichment for GOF, not for LOF**; strengthen that GOF maturate faster than LOF as compared to WT (the older neurons get, the more astrocyte they accumulate)

## GSEA on NSC, early/late-born neurons - ON A USER GENE LISTS

*From [EZH1 paper](https://www.nature.com/articles/s41467-023-39645-5) supp. data 3*; Top 500 significantly enriched genes in aRG, CFuPN, and CPN scRNA clusters from Uzquiano and [Kedaigle et al](https://doi.org/10.1016/j.cell.2022.09.010) to generate the NSC, early-born neuron and late-born neuron gene sets. 

```bash
conda activate deseq2
```

```R
# Packages
library("tidyverse")
library("clusterProfiler")
library("msigdbr") # BiocManager::install("msigdbr")
library("org.Mm.eg.db")
library("enrichplot") # for gseaplot2()
library("pheatmap")

# import DEGs 8wN
## filtered_8wN_KO_vs_8wN_WT
## filtered_8wN_HET_vs_8wN_WT
KO <- read.table("output/deseq2_hg38/filtered_8wN_KO_vs_8wN_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() 
KO_geneSymbol = KO %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene

HET <- read.table("output/deseq2_hg38/filtered_8wN_HET_vs_8wN_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()
HET_geneSymbol = HET %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene


# import gene set
NSC = read.table("output/gsea/aRG_NSC_gseaGeneList.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(cellName = "NSC")
early_born = read.table("output/gsea/CFuPN_earlyBorn_gseaGeneList.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(cellName = "early_born")
late_born = read.table("output/gsea/CPN_lateBorn_gseaGeneList.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(cellName = "late_born")

all_data <- NSC %>%
  bind_rows(early_born) %>%
  bind_rows(late_born)


# Order our DEG
## Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- KO_geneSymbol$log2FoldChange  ### CHAGNE HERE DATA!!!!!!!
names(lfc_vector) <- KO_geneSymbol$GeneSymbol ### CHAGNE HERE DATA!!!!!!!
## We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
### Set the seed so our results are reproducible:
set.seed(42)


# run GSEA
## without pvalue
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = all_data %>% dplyr::select(cellName,gene), # Need to be in that order...
)



gsea_result_df <- data.frame(gsea_results@result)
# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results__NSC_earlyLateNeurons_complete_WTvsHET.tsv"
  )
)

# plots
c("NSC", "early_born", "late_born")


pdf("output/gsea/NSC_earlyLateNeurons_WTvsHET.pdf", width=7, height=8)
pdf("output/gsea/NSC_earlyLateNeurons_WTvsKO.pdf", width=7, height=8)

enrichplot::gseaplot2(
  gsea_results,
  geneSetID = c("NSC", "early_born", "late_born"),
  title = "NSC_earlyLateNeurons"
)
dev.off()



# import DEGs NPC
## filtered_NPC_KO_vs_8wN_WT
## filtered_NPC_HET_vs_8wN_WT
KO <- read.table("output/deseq2_hg38/filtered_NPC_KO_vs_NPC_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() 
KO_geneSymbol = KO %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene

HET <- read.table("output/deseq2_hg38/filtered_NPC_HET_vs_NPC_WT.txt", header = TRUE, sep = "\t", row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()
HET_geneSymbol = HET %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene


# import gene set
NSC = read.table("output/gsea/aRG_NSC_gseaGeneList.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(cellName = "NSC")
early_born = read.table("output/gsea/CFuPN_earlyBorn_gseaGeneList.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(cellName = "early_born")
late_born = read.table("output/gsea/CPN_lateBorn_gseaGeneList.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(cellName = "late_born")

all_data <- NSC %>%
  bind_rows(early_born) %>%
  bind_rows(late_born)


# Order our DEG
## Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- KO_geneSymbol$log2FoldChange  ### CHAGNE HERE DATA!!!!!!!
names(lfc_vector) <- KO_geneSymbol$GeneSymbol ### CHAGNE HERE DATA!!!!!!!
## We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
### Set the seed so our results are reproducible:
set.seed(42)


# run GSEA
## without pvalue
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = all_data %>% dplyr::select(cellName,gene), # Need to be in that order...
)



gsea_result_df <- data.frame(gsea_results@result)
# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results__NSC_earlyLateNeurons_complete_WTvsHET_NPC.tsv"
  )
)

# plots
c("NSC", "early_born", "late_born")

pdf("output/gsea/NSC_earlyLateNeurons_WTvsKO_NPC.pdf", width=7, height=8)
pdf("output/gsea/NSC_earlyLateNeurons_WTvsHET_NPC.pdf", width=7, height=8)
pdf("output/gsea/NSC_earlyLateNeurons_WTvsKO_NPC_v1.pdf", width=10, height=8)

enrichplot::gseaplot2(
  gsea_results,
  geneSetID = c("NSC", "early_born", "late_born"),
  title = "NSC_earlyLateNeurons",
  base_size = 20
)
dev.off()



```

--> GOF/HET show late born neurons enrichment; LOF/KO show NPC enrichment




# Shiny app to explore RNAseq results

Let's create a shiny app for the lab members to explore the RNAseq results; (`conda activate scRNAseqV2`)

Version1= type a gene and display it's TPM value in barplot (for any genotype and time-points)

```R

# packages
library("shiny")
library("ggplot2")
library("dplyr")
library("tidyr")
library("tidyverse")
library("biomaRt")
library("rsconnect")

# Generate the log2tpm output file
## import tpm 
tpm_all <- read_csv("tpm_all_sample.txt") %>%   # tpm from all_sample.txt from tpm_hg38/
  dplyr::select(-1) #To import
## add geneSymbol
tpm_all$Geneid <- gsub("\\..*", "", tpm_all$Geneid) # remove Ensembl gene id version

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## Convert Ensembl gene IDs to gene symbols
tpm_all_genesymbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = tpm_all$Geneid,
                     mart = ensembl)

tpm_all_withGenesymbols = tpm_all %>% 
  dplyr::rename("ensembl_gene_id" = "Geneid") %>%
  left_join(tpm_all_genesymbols)


## Convert data from wide to long keep only geneSymbol
long_data <- tpm_all_withGenesymbols %>% 
  dplyr::select(-ensembl_gene_id) %>%
  drop_na() %>%
  tidyr::pivot_longer(-external_gene_name, names_to = "condition", values_to = "TPM") %>% 
  tidyr::separate(condition, into = c("Tissue", "Genotype", "Replicate"), sep = "_")

long_data_log2tpm = long_data %>%
  mutate(log2tpm = log2(TPM + 1)) %>%
  mutate(Genotype = recode(Genotype, "HET" = "GOF", "KO" = "LOF"))
## Save
write.table(long_data_log2tpm, file = c("output/shinyApp/RNAseqDataViewer_V1/long_data_log2tpm.txt"), sep = "\t", quote = FALSE, row.names = FALSE)



# ShinyApp
## 

global <- '
library(dplyr)

# import log2tpm file
long_data_log2tpm = read.table("long_data_log2tpm.txt", sep = "\t", header = TRUE)
'
writeLines(global, "output/shinyApp/RNAseqDataViewer_V1/global.R")


## WORK GREAT:

ui_code <- '
library(shiny)

fluidPage(
  titlePanel("RNAseq Data Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput("gene", "Select Gene:", unique(long_data_log2tpm$external_gene_name)),
      checkboxGroupInput("tissue", "Select Tissue:", choices = c("ESC", "NPC", "2dN", "4wN", "8wN")),
      checkboxGroupInput("genotype", "Select Genotype:", choices = c("WT", "GOF", "LOF")),
      numericInput("plot_width", label = "Plot Width", value = 500, min = 100, max = 2000, step = 10),
      numericInput("plot_height", label = "Plot Height", value = 500, min = 100, max = 2000, step = 10),
      radioButtons("plot_type", "Choose Plot Type:", choices = c("Bar Plot", "Line Plot"))
    ),
    
    mainPanel(
      uiOutput("dynamicPlotOutput")
    )
  )
)
'
writeLines(ui_code, "output/shinyApp/RNAseqDataViewer_V1/ui.R")



server_code <- '
library(shiny)
library(dplyr)
library(ggplot2)

server <- function(input, output) {
  
  plotWidth <- reactive({
    paste0(input$plot_width, "px")
  })

  plotHeight <- reactive({
    paste0(input$plot_height, "px")
  })

  output$dynamicPlotOutput <- renderUI({
    plotOutput("chosenPlot", width = plotWidth(), height = plotHeight())
  })
  
  output$chosenPlot <- renderPlot({
    subset_data <- long_data_log2tpm %>% 
      filter(external_gene_name == input$gene, Tissue %in% input$tissue, Genotype %in% input$genotype) %>% 
      group_by(Genotype, Tissue) %>% 
      summarise(
        Median = median(log2tpm),
        SEM = sd(log2tpm) / sqrt(n())
      ) %>%
      arrange(factor(Genotype, levels = c("WT", "GOF", "LOF")), factor(Tissue, levels = c("ESC", "NPC", "2dN", "4wN", "8wN")))

    if(input$plot_type == "Bar Plot") {
      return(
        ggplot(subset_data, aes(x = interaction(Tissue, Genotype, sep = "_"), y = Median, fill = Genotype)) + 
          geom_bar(stat = "identity", position = position_dodge()) + 
          geom_errorbar(aes(ymin = Median - SEM, ymax = Median + SEM), width = 0.25, position = position_dodge(0.9)) + 
          scale_fill_manual(values = c("WT" = "black", "GOF" = "blue", "LOF" = "red")) + 
          theme_bw() +
          labs(x = "Tissue_Genotype", y = "Median log2(tpm+1)") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_x_discrete(limits = unique(as.character(interaction(subset_data$Tissue, subset_data$Genotype, sep="_"))))
      )
    } else {
      return(
        ggplot(subset_data, aes(x = factor(Tissue, levels = c("ESC", "NPC", "2dN", "4wN", "8wN")), y = Median, color = Genotype, group = Genotype)) + 
          geom_line() + 
          geom_errorbar(aes(ymin = Median - SEM, ymax = Median + SEM), width = 0.2, size = 1) + 
          geom_point(shape = 18, size = 3, stroke = 1.5) +
          scale_color_manual(values = c("WT" = "black", "GOF" = "blue", "LOF" = "red")) + 
          theme_bw() +
          labs(x = "Tissue", y = "Median log2(tpm+1)") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
    }
  })
}
'
writeLines(server_code, "output/shinyApp/RNAseqDataViewer_V1/server.R")


# Deploy the app

rsconnect::deployApp("output/shinyApp/RNAseqDataViewer_V1")
                     

```


- NOTE: the `global.R` is run automaticcaly and detected by shiny app; should contain ll objects and meta


--> The version 1 works great!




## Shiny app V2; including Ciceri RNAseq neuron diff dataset


cp files to `long_data_log2tpm*.txt` files to `output/shinyApp/RNAseqDataViewer_V2`
- `output/tpm/long_data_log2tpm_Ciceri.txt` = Ciceri
- `output/shinyApp/*V1/long_data_log2tpm.txt`= our

for the line plot comparison, let's assume the following:
001 / Ciceri \ combined
ESC / ESC \ ESC
NPC / NPC \ NPC
2dN / 25dN \ 2dN_25dN
4wN / 50dN \ 4wN_50dN
8wN / 75dN \ 8wN_75dN
NA / 100dN \NA_100dN



```bash
conda activate scRNAseqV2
```


```R

# packages
library("shiny")
library("ggplot2")
library("dplyr")
library("tidyr")
library("tidyverse")
library("biomaRt")
library("rsconnect")






# prep the data

long_data_log2tpm = read.table("output/shinyApp/RNAseqDataViewer_V2/long_data_log2tpm.txt", sep = "\t", header = TRUE)
long_data_log2tpm_Ciceri = read.table("output/shinyApp/RNAseqDataViewer_V2/long_data_log2tpm_Ciceri.txt", sep = "\t", header = TRUE)

# add combined time
time_combined = tibble(
  Tissue = c("ESC", "NPC", "2dN", "4wN", "8wN", NA),
  Ciceri = c("ESC", "NPC", "25dN", "50dN", "75dN", "100dN"),
  "Akizu / Ciceri" = c("ESC / ESC", "NPC / NPC", "2dN / 25dN", "4wN / 50dN", "8wN / 75dN", "NA / 100dN")
)



# ShinyApp 

global <- '
library(dplyr)

# Assuming the data is read from the respective files for Akizu and Ciceri
long_data_log2tpm_Akizu = read.table("long_data_log2tpm.txt", sep = "\t", header = TRUE)
long_data_log2tpm_Ciceri = read.table("long_data_log2tpm_Ciceri.txt", sep = "\t", header = TRUE)
'
writeLines(global, "output/shinyApp/RNAseqDataViewer_V2/global.R")


## WORK GREAT:

ui_code <- '

library(shiny)

ui <- fluidPage(
  titlePanel("RNAseq Data Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Select Dataset:", choices = c("Akizu", "Ciceri")),
      selectizeInput("gene", "Select Gene:", choices = NULL, options = list("server" = TRUE)),
      checkboxGroupInput("tissue", "Select Tissue:", choices = NULL),
      checkboxGroupInput("genotype", "Select Genotype:", choices = NULL),
      numericInput("plot_width", label = "Plot Width", value = 500, min = 100, max = 2000, step = 10),
      numericInput("plot_height", label = "Plot Height", value = 500, min = 100, max = 2000, step = 10),
      radioButtons("plot_type", "Choose Plot Type:", choices = c("Bar Plot", "Line Plot"))
    ),
    
    mainPanel(
      plotOutput("chosenPlot", width = "100%")  # Removed the height function here
    )
  )
)

'
writeLines(ui_code, "output/shinyApp/RNAseqDataViewer_V2/ui.R")



server_code <- '
library(shiny)
library(dplyr)
library(ggplot2)

server <- function(input, output, session) {
  
  # Update gene selection based on the dataset
  observeEvent(input$dataset, {
    if (input$dataset == "Akizu") {
      genes <- unique(long_data_log2tpm_Akizu$external_gene_name)
      tissues <- c("ESC", "NPC", "2dN", "4wN", "8wN")
      genotypes <- c("WT", "GOF", "LOF")
    } else if (input$dataset == "Ciceri") {
      genes <- unique(long_data_log2tpm_Ciceri$external_gene_name)
      tissues <- c("ESC", "NPC", "25dN", "50dN", "75dN", "100dN")
      genotypes <- c("WT")
    }
    updateSelectizeInput(session, "gene", choices = genes, server = TRUE)
    updateCheckboxGroupInput(session, "tissue", choices = tissues)
    updateCheckboxGroupInput(session, "genotype", choices = genotypes)
  })

  # Reactive expression to filter data based on selections
  filtered_data <- reactive({
    if (input$dataset == "Akizu") {
      data <- long_data_log2tpm_Akizu
    } else if (input$dataset == "Ciceri") {
      data <- long_data_log2tpm_Ciceri
    }

    data <- data %>%
      filter(external_gene_name == input$gene,
             Tissue %in% input$tissue,
             Genotype %in% input$genotype)

    return(data)
  })

  # Render the plot output
  output$chosenPlot <- renderPlot({
    data <- filtered_data()
    
    if (nrow(data) == 0) {
      return(NULL)
    }

    data_summary <- data %>%
      group_by(Genotype, Tissue) %>%
      summarise(
        Median = median(log2tpm, na.rm = TRUE),
        SEM = sd(log2tpm, na.rm = TRUE) / sqrt(n())
      ) %>%
      ungroup() %>% 
      filter(Tissue %in% input$tissue)

    # Plot the data based on the selected plot type
    if(input$plot_type == "Bar Plot") {
      ggplot(data_summary, aes(x = Tissue, y = Median, fill = Genotype)) + 
        geom_bar(stat = "identity", position = position_dodge()) + 
        geom_errorbar(aes(ymin = Median - SEM, ymax = Median + SEM), width = 0.25, position = position_dodge(0.9)) + 
        scale_fill_manual(values = c("WT" = "black", "GOF" = "blue", "LOF" = "red")) + 
        theme_bw() +
        labs(x = "Tissue", y = "Median log2(tpm+1)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(limits = input$tissue)
    } else {
      ggplot(data_summary, aes(x = Tissue, y = Median, color = Genotype, group = Genotype)) + 
        geom_line() + 
        geom_errorbar(aes(ymin = Median - SEM, ymax = Median + SEM), width = 0.2, size = 1) + 
        geom_point(shape = 18, size = 3, stroke = 1.5) +
        scale_color_manual(values = c("WT" = "black", "GOF" = "blue", "LOF" = "red")) + 
        theme_bw() +
        labs(x = "Tissue", y = "Median log2(tpm+1)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(limits = input$tissue)
    }
  }, width = function() { 
    input$plot_width 
  }, height = function() { 
    if (is.numeric(input$plot_height) && input$plot_height > 0) {
      input$plot_height
    } else {
      400  # Default height if input is invalid
    }
  })
}
'


writeLines(server_code, "output/shinyApp/RNAseqDataViewer_V2/server.R")


# Deploy the app

rsconnect::deployApp("output/shinyApp/RNAseqDataViewer_V2")
                     

```

 










In version 3, let's implement other functino such as:
- looking at multiple genes? Heatmap?
- statistical test?










# MuSiC2 RNAseq bulk deconvolution

Let's try to use [MuSiC2](https://academic.oup.com/bib/article-abstract/23/6/bbac430/6751147?redirectedFrom=fulltext) for RNAseq deconvolution. Good as it use just 1 scRNAseq file for reference and we can compare different bulk RNAseq samples! Otherwise we need scRNAseq for each of our genotype...


- scRNAseq signature: [Human brain organoid](https://www.nature.com/articles/s41586-023-06473-y); files downloaded `CHOOSE_CTRL_annot_srt.rds`
- [Tutorial](https://jiaxin-fan.github.io/MuSiC2/articles/introduction.html) to follow for MuSiC2 
- Antother more recent [Tutorial](https://xuranw.github.io/MuSiC/articles/pages/MuSiC2.html)

--> Install MuSic; need devtools and Seurat so let's clone `scRNAseqV2` and install MuSiC2:

```bash
conda activate scRNAseqV2

conda create --name scRNAseqV3 --clone scRNAseqV2

conda activate scRNAseqV3
```
```R
# install MuSiC2; requires installation of dependencies separaetly
BiocManager::install("TOAST")
devtools::install_github('xuranw/MuSiC')
remotes::install_github("renozao/xbioc")
if (!"MuSiC2" %in% rownames(installed.packages())) {
  devtools::install_github('Jiaxin-Fan/MuSiC2')
}
```


Check wether the control scRNAseq seurat object is good (generate UMAP with cell type label); is so, start MuSiC with `conda activate scRNAseqV3`

--> The file looks good, but very few cells for the control `CHOOSE_CTRL_annot_srt.rds`; maybe better and to make sure to download the whole set `CHOOSE_full_dataset_srt.rds` and subtract from it the control cells.

--> For Bulk RNAseq, seems I need to use the raw read counts (from tutorial)

```R

# library
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") 
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")
library("MuSiC2")
library("MuSiC")
library("biomaRt")


# CHECK CONTROL FILE


## import seurat object
CHOOSE_CTRL_annot_srt <- readRDS(file = "output/deconv/CHOOSE_CTRL_annot_srt.rds")
DefaultAssay(CHOOSE_CTRL_annot_srt) <- "RNA" # 


pdf("output/deconv/UMAP_CHOOSE_CTRL_annot_srt.celltype_cl_coarse2.pdf", width=10, height=6)
DimPlot(CHOOSE_CTRL_annot_srt, reduction = "umap", group.by = "celltype_cl_coarse2", label=TRUE)
dev.off()

## isolate the cells of interest; control one
cells_to_keep <- WhichCells(CHOOSE_CTRL_annot_srt, expression = celltype_cl_coarse2 != "NA")

CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2 <- subset(CHOOSE_CTRL_annot_srt, cells = cells_to_keep)


pdf("output/deconv/UMAP_CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.pdf", width=10, height=6)
DimPlot(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2, reduction = "umap", group.by = "celltype_cl_coarse2", label=TRUE)
dev.off()

pdf("output/deconv/UMAP_CHOOSE_CTRL_annot_srt.subset_pseudotime_ranks.pdf", width=10, height=6)
FeaturePlot(
  object = CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2,
  features = 'pseudotime_ranks', 
  reduction = 'umap', 
  pt.size = 1, # Adjust point size if needed
  cols = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")) # Use reverse Spectral palette for color gradient
)
dev.off()

## Transform scRNAseq data in ExpressionSet Class
CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.sceset = ExpressionSet(assayData = as.matrix(GetAssayData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2)), phenoData =  new("AnnotatedDataFrame",CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2@meta.data))


## Transform scRNAseq data in SingleCellExperiment (for music2 t statistics)

CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment <- as.SingleCellExperiment(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2, assay = "RNA")







# Import raw RNAseq read counts
# code for ensembl to genesymbol conv
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
genes <- X8wN_WT_R1$Geneid
#### Get the mapping from Ensembl ID to gene symbol
genes_mapped  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = genes,
                      mart = ensembl)
####

### import featureCounts output
#### 8wN WT R1
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
X8wN_WT_R1_geneSymbol <- merge(X8wN_WT_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam')
## Calculate median as gene dupplicated with the conversino!
X8wN_WT_R1_count_summary <- X8wN_WT_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R1_count_matrix <- as.matrix(X8wN_WT_R1_count_summary$median_count)
rownames(X8wN_WT_R1_count_matrix) <- X8wN_WT_R1_count_summary$external_gene_name
colnames(X8wN_WT_R1_count_matrix) <- "8wN_WT_R1"



#### 8wN WT R2
X8wN_WT_R2 <- read.delim("output/featurecounts_hg38/8wN_WT_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R2$Geneid <- gsub("\\..*", "", X8wN_WT_R2$Geneid)
X8wN_WT_R2_geneSymbol <- merge(X8wN_WT_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R2_count_summary <- X8wN_WT_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R2_count_matrix <- as.matrix(X8wN_WT_R2_count_summary$median_count)
rownames(X8wN_WT_R2_count_matrix) <- X8wN_WT_R2_count_summary$external_gene_name
colnames(X8wN_WT_R2_count_matrix) <- "8wN_WT_R2"


#### 8wN WT R3
X8wN_WT_R3 <- read.delim("output/featurecounts_hg38/8wN_WT_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R3$Geneid <- gsub("\\..*", "", X8wN_WT_R3$Geneid)
X8wN_WT_R3_geneSymbol <- merge(X8wN_WT_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R3_count_summary <- X8wN_WT_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R3_count_matrix <- as.matrix(X8wN_WT_R3_count_summary$median_count)
rownames(X8wN_WT_R3_count_matrix) <- X8wN_WT_R3_count_summary$external_gene_name
colnames(X8wN_WT_R3_count_matrix) <- "8wN_WT_R3"


#### 8wN WT R4
X8wN_WT_R4 <- read.delim("output/featurecounts_hg38/8wN_WT_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R4$Geneid <- gsub("\\..*", "", X8wN_WT_R4$Geneid)
X8wN_WT_R4_geneSymbol <- merge(X8wN_WT_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R4_count_summary <- X8wN_WT_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R4_count_matrix <- as.matrix(X8wN_WT_R4_count_summary$median_count)
rownames(X8wN_WT_R4_count_matrix) <- X8wN_WT_R4_count_summary$external_gene_name
colnames(X8wN_WT_R4_count_matrix) <- "8wN_WT_R4"


#### 8wN HET R1
X8wN_HET_R1 <- read.delim("output/featurecounts_hg38/8wN_HET_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R1$Geneid <- gsub("\\..*", "", X8wN_HET_R1$Geneid)
X8wN_HET_R1_geneSymbol <- merge(X8wN_HET_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R1_count_summary <- X8wN_HET_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R1_count_matrix <- as.matrix(X8wN_HET_R1_count_summary$median_count)
rownames(X8wN_HET_R1_count_matrix) <- X8wN_HET_R1_count_summary$external_gene_name
colnames(X8wN_HET_R1_count_matrix) <- "8wN_HET_R1"



#### 8wN HET R2
X8wN_HET_R2 <- read.delim("output/featurecounts_hg38/8wN_HET_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R2$Geneid <- gsub("\\..*", "", X8wN_HET_R2$Geneid)
X8wN_HET_R2_geneSymbol <- merge(X8wN_HET_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R2_count_summary <- X8wN_HET_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R2_count_matrix <- as.matrix(X8wN_HET_R2_count_summary$median_count)
rownames(X8wN_HET_R2_count_matrix) <- X8wN_HET_R2_count_summary$external_gene_name
colnames(X8wN_HET_R2_count_matrix) <- "8wN_HET_R2"


#### 8wN HET R3
X8wN_HET_R3 <- read.delim("output/featurecounts_hg38/8wN_HET_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R3$Geneid <- gsub("\\..*", "", X8wN_HET_R3$Geneid)
X8wN_HET_R3_geneSymbol <- merge(X8wN_HET_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R3_count_summary <- X8wN_HET_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R3_count_matrix <- as.matrix(X8wN_HET_R3_count_summary$median_count)
rownames(X8wN_HET_R3_count_matrix) <- X8wN_HET_R3_count_summary$external_gene_name
colnames(X8wN_HET_R3_count_matrix) <- "8wN_HET_R3"


#### 8wN HET R4
X8wN_HET_R4 <- read.delim("output/featurecounts_hg38/8wN_HET_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R4$Geneid <- gsub("\\..*", "", X8wN_HET_R4$Geneid)
X8wN_HET_R4_geneSymbol <- merge(X8wN_HET_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R4_count_summary <- X8wN_HET_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R4_count_matrix <- as.matrix(X8wN_HET_R4_count_summary$median_count)
rownames(X8wN_HET_R4_count_matrix) <- X8wN_HET_R4_count_summary$external_gene_name
colnames(X8wN_HET_R4_count_matrix) <- "8wN_HET_R4"


#### 8wN KO R1
X8wN_KO_R1 <- read.delim("output/featurecounts_hg38/8wN_KO_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R1$Geneid <- gsub("\\..*", "", X8wN_KO_R1$Geneid)
X8wN_KO_R1_geneSymbol <- merge(X8wN_KO_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R1_count_summary <- X8wN_KO_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R1_count_matrix <- as.matrix(X8wN_KO_R1_count_summary$median_count)
rownames(X8wN_KO_R1_count_matrix) <- X8wN_KO_R1_count_summary$external_gene_name
colnames(X8wN_KO_R1_count_matrix) <- "8wN_KO_R1"



#### 8wN KO R2
X8wN_KO_R2 <- read.delim("output/featurecounts_hg38/8wN_KO_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R2$Geneid <- gsub("\\..*", "", X8wN_KO_R2$Geneid)
X8wN_KO_R2_geneSymbol <- merge(X8wN_KO_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R2_count_summary <- X8wN_KO_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R2_count_matrix <- as.matrix(X8wN_KO_R2_count_summary$median_count)
rownames(X8wN_KO_R2_count_summary) <- X8wN_KO_R2_count_summary$external_gene_name
colnames(X8wN_KO_R2_count_summary) <- "8wN_KO_R2"



#### 8wN KO R3
X8wN_KO_R3 <- read.delim("output/featurecounts_hg38/8wN_KO_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R3$Geneid <- gsub("\\..*", "", X8wN_KO_R3$Geneid)
X8wN_KO_R3_geneSymbol <- merge(X8wN_KO_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R3_count_summary <- X8wN_KO_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R3_count_matrix <- as.matrix(X8wN_KO_R3_count_summary$median_count)
rownames(X8wN_KO_R3_count_summary) <- X8wN_KO_R3_count_summary$external_gene_name
colnames(X8wN_KO_R3_count_summary) <- "8wN_KO_R3"



#### 8wN KO R4
X8wN_KO_R4 <- read.delim("output/featurecounts_hg38/8wN_KO_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R4$Geneid <- gsub("\\..*", "", X8wN_KO_R4$Geneid)
X8wN_KO_R4_geneSymbol <- merge(X8wN_KO_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R4_count_summary <- X8wN_KO_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R4_count_matrix <- as.matrix(X8wN_KO_R4_count_summary$median_count)
rownames(X8wN_KO_R4_count_summary) <- X8wN_KO_R4_count_summary$external_gene_name
colnames(X8wN_KO_R4_count_summary) <- "8wN_KO_R4"



## WT versus HET comparison
all_counts <- cbind(X8wN_WT_R1_count_matrix, X8wN_WT_R2_count_matrix, X8wN_WT_R3_count_matrix, X8wN_WT_R4_count_matrix, X8wN_HET_R1_count_matrix, X8wN_HET_R2_count_matrix, X8wN_HET_R3_count_matrix, X8wN_HET_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("WT", "WT","WT", "WT", "HET", "HET", "HET", "HET"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
WT_HET_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)
table(WT_HET_8wN.bulkeset$group)

bulk.control.mtx = exprs(WT_HET_8wN.bulkeset)[, WT_HET_8wN.bulkeset$group == 'WT']
bulk.case.mtx = exprs(WT_HET_8wN.bulkeset)[, WT_HET_8wN.bulkeset$group == 'HET']

# TO AVOID SUBSCRIPTY OUT OF BOUND ERROR:
### subset bulk features that are only present in scRNAseq assay
###### Assuming 'featureNames' can be used to retrieve the features from both ExpressionSets
common_features <- intersect(featureNames(WT_HET_8wN.bulkeset), featureNames(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.sceset))
# Subset the WT_HET_8wN.bulkeset to only include common features
WT_HET_8wN.bulkeset_subset <- WT_HET_8wN.bulkeset[common_features, ]


bulk.control.mtx = exprs(WT_HET_8wN.bulkeset_subset)[, WT_HET_8wN.bulkeset_subset$group == 'WT']
bulk.case.mtx = exprs(WT_HET_8wN.bulkeset_subset)[, WT_HET_8wN.bulkeset_subset$group == 'HET']

# music2 deconvolution music2_prop_t_statistics
set.seed(42)
CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment$sampleID <- rownames(colData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment))

est = music2_prop_t_statistics(bulk.control.mtx = bulk.control.mtx, bulk.case.mtx = bulk.case.mtx, sc.sce = CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment, clusters = 'celltype_cl_coarse2', samples = 'sampleID', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "IP", "INP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "oRG", "RG","vRG")
, n_resample=20, sample_prop=0.5,cutoff_c=0.05,cutoff_r=0.01)

est.prop = est$Est.prop










# music2 deconvolution TOAST
## The TOAST is buggy, create a new function modifying the source code
set.seed(42)
CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment$sampleID <- rownames(colData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment))

est = my_music2_prop_toast_V2(bulk.control.mtx = bulk.control.mtx, bulk.case.mtx = bulk.case.mtx, sc.sce = CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment, clusters = 'celltype_cl_coarse2', samples = 'sampleID', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "IP", "INP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "oRG", "RG","vRG"))







```

--> Not sure what is the correct cluster name:
- celltype_jf_refined: looks good
- state: neuron, progenitor
- lineage: excitatory, inhibitory
- lineage_jf: excitatory, inhibitory but broader
- gRNA: Control2 and NA
- induced: False (maybe false is control?)
- **celltype_cl_coarse2**: astrocyte present; and all layer of excitatory and inhibitory!!!
- **celltype_cl_refined**: same as celltype_cl_coarse2
- pseudotime_ranks: pseudotime information


--> For the control file `CHOOSE_CTRL_annot_srt`; celltype_cl_coarse2 is the file to use; only non-NA cells as been extracted. Pseudtoime is in agreement with cell types

- *NOTE: got error with subscrupt out of bound when running `music2_prop_t_statistics()`. I ended up subtract my bulkRNAseq object to only keep the features present in the sCRNAseq obect. Issue discussed [here](https://github.com/xuranw/MuSiC/issues/61)*
- *NOTE" got error `subscript contains invalid names` when running `music2_prop_t_statistics()`and solved (here)[https://github.com/Jiaxin-Fan/MuSiC2/issues/6]*



--> `music2_prop_t_statistics()` is too long to run in local; here a custom script top run it:
```bash
conda activate scRNAseqV3
# 8wN WT vs HET using Control scRNAseq data
sbatch scripts/music2_prop_t_statistics.8wN_WTvsHET_CHOOSE_CTRL.sh # 6729977 ok FAIL with if (sum(abs(p.weight.new - p.weight)) error...
sbatch scripts/my_music2_prop_toast_V2.8wN_WTvsHET_CHOOSE_CTRL.sh # 


```





New `music2_prop_toast()` funcction for Music2 TOAST, solution found [here](https://github.com/xuranw/MuSiC/issues/108).
But this version is still buggy with error:

```
Error in if (sum(abs(p.weight.new - p.weight)) < eps) { :
missing value where TRUE/FALSE needed
```

So let's try generate a `my_music2_prop_toast_V2` too

```R
my_music2_prop_toast <- function (bulk.control.mtx, bulk.case.mtx, sc.sce, clusters, 
          samples, select.ct, expr_low = 20, prop_r = 0.1, eps_c = 0.05, 
          eps_r = 0.01, cutoff_c = 10^(-3), cutoff_r = 10^(-3), cap = 0.3, 
          maxiter = 200, markers = NULL, ct.cov = FALSE, cell_size = NULL,
          centered = FALSE, normalize = FALSE) {
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if (length(gene.bulk) < 0.1 * min(nrow(bulk.control.mtx), 
                                    nrow(bulk.case.mtx))) {
    stop("Not enough genes for bulk data! Please check gene annotations.")
  }
  bulk.mtx = cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, 
  ])
  Pheno = data.frame(condition = factor(c(rep("control", ncol(bulk.control.mtx)), 
                                          rep("case", ncol(bulk.case.mtx))), levels = c("control", 
                                                                                        "case")))
  rownames(Pheno) = colnames(bulk.mtx)
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if (length(gene_all) < 0.2 * min(length(gene.bulk), nrow(sc.sce))) {
    stop("Not enough genes between bulk and single-cell data! Please check gene annotations.")
  }
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr >= expr_low])
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)]
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)]
  prop_control = music_prop(bulk.mtx = bulk.control, sc.sce = sc.sce, 
                            clusters = clusters, samples = samples, select.ct = select.ct, 
                            markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                            iter.max = 1000, nu = 1e-04, eps = 0.01, centered = centered, 
                            normalize = normalize, verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.mtx = bulk.case, sc.sce = sc.sce, 
                             clusters = clusters, samples = samples, select.ct = select.ct, 
                             markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                             iter.max = 1000, nu = 1e-04, eps = 0.01, centered = centered, 
                             normalize = normalize, verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control, prop_CASE)
  iter = 1
  ncell = length(select.ct)
  id_conv = NULL
  while (iter <= maxiter) {
    Y_raw = log1p(bulk.mtx)
    design = Pheno
    Prop <- prop_all[rownames(Pheno), ]
    Design_out <- makeDesign(design, Prop)
    fitted_model <- fitModel(Design_out, Y_raw)
    res_table <- csTest(fitted_model, coef = "condition", verbose = F)
    mex = apply(prop_all, 2, mean)
    lr = NULL
    for (celltype in select.ct) {
      m = mex[celltype]
      DE = res_table[[celltype]]
      pval = DE$fdr
      names(pval) = rownames(DE)
      pval = pval[names(pval) %in% exp_genel]
      if (m >= prop_r) {
        lr = c(lr, names(pval[pval <= cutoff_c & pval <= 
                                quantile(pval, prob = cap)]))
      }
      else {
        lr = c(lr, names(pval[pval <= cutoff_r & pval <= 
                                quantile(pval, prob = cap)]))
      }
    }
    lr = unique(lr)
    l = setdiff(gene_all, lr)
    sc.iter.sce = sc.sce[l, ]
    if (length(id_conv) > 0) {
      case_sample = bulk.case[, !colnames(bulk.case) %in% 
                                id_conv]
    }
    else {
      case_sample = bulk.case
    }
    prop_case = music_prop(bulk.mtx = case_sample, sc.sce = sc.iter.sce, 
                           clusters = clusters, samples = samples, select.ct = select.ct, 
                           verbose = F)$Est.prop.weighted
    prop_CASE = rbind(prop_case, prop_case_fix)
    if (length(id_conv) == 1) {
      rownames(prop_CASE) = c(rownames(prop_case), id_conv)
    }
    prop_all = rbind(prop_control, prop_CASE)
    prop_case = prop_case[rownames(prop_case_ini), ]
    pc = abs(prop_case - prop_case_ini)
    conv = pc
    conv[, ] = 1
    conv[prop_case_ini <= prop_r] = ifelse(pc[prop_case_ini <= 
                                                prop_r] < eps_r, 0, 1)
    pc[prop_case_ini > prop_r] = pc[prop_case_ini > prop_r]/prop_case_ini[prop_case_ini > 
                                                                            prop_r]
    conv[prop_case_ini > prop_r] = ifelse(pc[prop_case_ini > 
                                               prop_r] < eps_c, 0, 1)
    convf = apply(conv, 1, function(x) {
      all(x == 0)
    })
    all_converge = FALSE
    id_conv = c(id_conv, names(convf[convf == TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% 
                                id_conv, ]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv, 
    ]
    if (is.vector(prop_case_ini)) {
      all_converge = TRUE
      break
    }
    else if (nrow(prop_case_ini) == 0) {
      all_converge = TRUE
      break
    }
    iter = iter + 1
  }
  if (all_converge) {
    return(list(Est.prop = prop_all, convergence = TRUE, 
                n.iter = iter, DE.genes = lr))
  }
  else {
    return(list(Est.prop = prop_all, convergence = FALSE, 
                id.not.converge = rownames(prop_case_ini)))
  }
}














my_music2_prop_toast_V2 <- function (bulk.control.mtx, bulk.case.mtx, sc.sce, clusters, 
          samples, select.ct, expr_low = 20, prop_r = 0.1, eps_c = 0.05, 
          eps_r = 0.01, cutoff_c = 10^(-3), cutoff_r = 10^(-3), cap = 0.3, 
          maxiter = 200, markers = NULL, ct.cov = FALSE, cell_size = NULL,
          centered = FALSE, normalize = FALSE) {
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if (length(gene.bulk) < 0.1 * min(nrow(bulk.control.mtx), 
                                    nrow(bulk.case.mtx))) {
    stop("Not enough genes for bulk data! Please check gene annotations.")
  }
  bulk.mtx = cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, 
  ])
  Pheno = data.frame(condition = factor(c(rep("control", ncol(bulk.control.mtx)), 
                                          rep("case", ncol(bulk.case.mtx))), levels = c("control", 
                                                                                        "case")))
  rownames(Pheno) = colnames(bulk.mtx)
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if (length(gene_all) < 0.2 * min(length(gene.bulk), nrow(sc.sce))) {
    stop("Not enough genes between bulk and single-cell data! Please check gene annotations.")
  }
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr >= expr_low])
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)]
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)]
  prop_control = music_prop(bulk.mtx = bulk.control, sc.sce = sc.sce, 
                            clusters = clusters, samples = samples, select.ct = select.ct, 
                            markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                            iter.max = 1000, nu = 1e-04, eps = 0.01, centered = centered, 
                            normalize = normalize, verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.mtx = bulk.case, sc.sce = sc.sce, 
                             clusters = clusters, samples = samples, select.ct = select.ct, 
                             markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                             iter.max = 1000, nu = 1e-04, eps = 0.01, centered = centered, 
                             normalize = normalize, verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control, prop_CASE)
  iter = 1
  ncell = length(select.ct)
  id_conv = NULL
  while (iter <= maxiter) {
    Y_raw = log1p(bulk.mtx)
    design = Pheno
    Prop <- prop_all[rownames(Pheno), ]
    Design_out <- makeDesign(design, Prop)
    fitted_model <- fitModel(Design_out, Y_raw)
    res_table <- csTest(fitted_model, coef = "condition", verbose = F)
    mex = apply(prop_all, 2, mean)
    lr = NULL
    for (celltype in select.ct) {
      m = mex[celltype]
      DE = res_table[[celltype]]
      pval = DE$fdr
      names(pval) = rownames(DE)
      pval = pval[names(pval) %in% exp_genel]
      if (m >= prop_r) {
        lr = c(lr, names(pval[pval <= cutoff_c & pval <= 
                                quantile(pval, prob = cap)]))
      }
      else {
        lr = c(lr, names(pval[pval <= cutoff_r & pval <= 
                                quantile(pval, prob = cap)]))
      }
    }
    lr = unique(lr)
    l = setdiff(gene_all, lr)
    sc.iter.sce = sc.sce[l, ]
    if (length(id_conv) > 0) {
      case_sample = bulk.case[, !colnames(bulk.case) %in% 
                                id_conv]
    }
    else {
      case_sample = bulk.case
    }
    prop_case = music_prop(bulk.mtx = case_sample, sc.sce = sc.iter.sce, 
                           clusters = clusters, samples = samples, select.ct = select.ct, 
                           verbose = F)$Est.prop.weighted
    prop_CASE = rbind(prop_case, prop_case_fix)
    if (length(id_conv) == 1) {
      rownames(prop_CASE) = c(rownames(prop_case), id_conv)
    }
    prop_all = rbind(prop_control, prop_CASE)
    prop_case = prop_case[rownames(prop_case_ini), ]

    replace_problematic_values <- function(x, replace_with = 0) {
      ifelse(is.na(x) | is.nan(x) | is.infinite(x), replace_with, x)
    }

    # Ensure 'prop_case' and 'prop_case_ini' do not contain problematic values
    prop_case <- replace_problematic_values(prop_case)
    prop_case_ini <- replace_problematic_values(prop_case_ini)

    # Use 'safe_division' function to avoid division by zero
    safe_division <- function(x, y) {
      ifelse(y == 0, 0, x / y)
    }

    # Calculate absolute differences and replace problematic values with zero
    pc <- abs(prop_case - prop_case_ini)
    pc <- replace_problematic_values(pc)

    # Initialize convergence matrix with ones
    conv <- matrix(1, nrow = nrow(pc), ncol = ncol(pc))

    # Update the convergence matrix based on conditions
    conv[prop_case_ini <= prop_r] <- ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
    pc[prop_case_ini > prop_r] <- safe_division(pc[prop_case_ini > prop_r], prop_case_ini[prop_case_ini > prop_r])
    conv[prop_case_ini > prop_r] <- ifelse(pc[prop_case_ini > prop_r] < eps_c, 0, 1)


                                               
    convf = apply(conv, 1, function(x) {
      all(x == 0)
    })
    all_converge = FALSE
    id_conv = c(id_conv, names(convf[convf == TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% 
                                id_conv, ]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv, 
    ]
    if (is.vector(prop_case_ini)) {
      all_converge = TRUE
      break
    }
    else if (nrow(prop_case_ini) == 0) {
      all_converge = TRUE
      break
    }
    iter = iter + 1
  }
  if (all_converge) {
    return(list(Est.prop = prop_all, convergence = TRUE, 
                n.iter = iter, DE.genes = lr))
  }
  else {
    return(list(Est.prop = prop_all, convergence = FALSE, 
                id.not.converge = rownames(prop_case_ini)))
  }
}


```

--> MuSiC2 isw no more supported, the author is NOT responsive on github...
----> Let's try music version1

# MuSiC1 

Try MuSiC1 and do deconvolution genotype per genotype and then compare.
--> Follow this [tutorial](https://xuranw.github.io/MuSiC/articles/MuSiC.html)

```bash
conda activate scRNAseqV3
```




```R

# library
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") 
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")
library("MuSiC2")
library("MuSiC")
library("biomaRt")



# scRNASeq import all
CHOOSE_full_dataset_srt <- readRDS(file = "output/deconv/CHOOSE_full_dataset_srt.rds")
DefaultAssay(CHOOSE_full_dataset_srt) <- "RNA" # 
## Transform scRNAseq data in SingleCellExperiment
CHOOSE_full_dataset_srt.SingleCellExperiment <- as.SingleCellExperiment(CHOOSE_full_dataset_srt, assay = "RNA")
## Transform scRNAseq data in SingleCellExperiment (FOR the out of bound error)
CHOOSE_full_dataset_srt.sceset = ExpressionSet(assayData = as.matrix(GetAssayData(CHOOSE_full_dataset_srt)), phenoData =  new("AnnotatedDataFrame",CHOOSE_full_dataset_srt@meta.data))


# Import raw RNAseq read counts
### code for ensembl to genesymbol conv
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <-  useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "uswest") # uswest useast asia
ensembl <-  useEnsembl(biomart = 'ensembl', 
                       dataset = 'hsapiens_gene_ensembl',
                       version = 110)
##                       
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
genes <- X8wN_WT_R1$Geneid
#### Get the mapping from Ensembl ID to gene symbol
genes_mapped  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = genes,
                      mart = ensembl)
####

### import featureCounts output
#### 8wN WT R1
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
X8wN_WT_R1_geneSymbol <- merge(X8wN_WT_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam')
## Calculate median as gene dupplicated with the conversino!
X8wN_WT_R1_count_summary <- X8wN_WT_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R1_count_matrix <- as.matrix(X8wN_WT_R1_count_summary$median_count)
rownames(X8wN_WT_R1_count_matrix) <- X8wN_WT_R1_count_summary$external_gene_name
colnames(X8wN_WT_R1_count_matrix) <- "8wN_WT_R1"



#### 8wN WT R2
X8wN_WT_R2 <- read.delim("output/featurecounts_hg38/8wN_WT_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R2$Geneid <- gsub("\\..*", "", X8wN_WT_R2$Geneid)
X8wN_WT_R2_geneSymbol <- merge(X8wN_WT_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R2_count_summary <- X8wN_WT_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R2_count_matrix <- as.matrix(X8wN_WT_R2_count_summary$median_count)
rownames(X8wN_WT_R2_count_matrix) <- X8wN_WT_R2_count_summary$external_gene_name
colnames(X8wN_WT_R2_count_matrix) <- "8wN_WT_R2"


#### 8wN WT R3
X8wN_WT_R3 <- read.delim("output/featurecounts_hg38/8wN_WT_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R3$Geneid <- gsub("\\..*", "", X8wN_WT_R3$Geneid)
X8wN_WT_R3_geneSymbol <- merge(X8wN_WT_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R3_count_summary <- X8wN_WT_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R3_count_matrix <- as.matrix(X8wN_WT_R3_count_summary$median_count)
rownames(X8wN_WT_R3_count_matrix) <- X8wN_WT_R3_count_summary$external_gene_name
colnames(X8wN_WT_R3_count_matrix) <- "8wN_WT_R3"


#### 8wN WT R4
X8wN_WT_R4 <- read.delim("output/featurecounts_hg38/8wN_WT_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R4$Geneid <- gsub("\\..*", "", X8wN_WT_R4$Geneid)
X8wN_WT_R4_geneSymbol <- merge(X8wN_WT_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R4_count_summary <- X8wN_WT_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R4_count_matrix <- as.matrix(X8wN_WT_R4_count_summary$median_count)
rownames(X8wN_WT_R4_count_matrix) <- X8wN_WT_R4_count_summary$external_gene_name
colnames(X8wN_WT_R4_count_matrix) <- "8wN_WT_R4"


#### 8wN HET R1
X8wN_HET_R1 <- read.delim("output/featurecounts_hg38/8wN_HET_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R1$Geneid <- gsub("\\..*", "", X8wN_HET_R1$Geneid)
X8wN_HET_R1_geneSymbol <- merge(X8wN_HET_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R1_count_summary <- X8wN_HET_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R1_count_matrix <- as.matrix(X8wN_HET_R1_count_summary$median_count)
rownames(X8wN_HET_R1_count_matrix) <- X8wN_HET_R1_count_summary$external_gene_name
colnames(X8wN_HET_R1_count_matrix) <- "8wN_HET_R1"



#### 8wN HET R2
X8wN_HET_R2 <- read.delim("output/featurecounts_hg38/8wN_HET_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R2$Geneid <- gsub("\\..*", "", X8wN_HET_R2$Geneid)
X8wN_HET_R2_geneSymbol <- merge(X8wN_HET_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R2_count_summary <- X8wN_HET_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R2_count_matrix <- as.matrix(X8wN_HET_R2_count_summary$median_count)
rownames(X8wN_HET_R2_count_matrix) <- X8wN_HET_R2_count_summary$external_gene_name
colnames(X8wN_HET_R2_count_matrix) <- "8wN_HET_R2"


#### 8wN HET R3
X8wN_HET_R3 <- read.delim("output/featurecounts_hg38/8wN_HET_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R3$Geneid <- gsub("\\..*", "", X8wN_HET_R3$Geneid)
X8wN_HET_R3_geneSymbol <- merge(X8wN_HET_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R3_count_summary <- X8wN_HET_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R3_count_matrix <- as.matrix(X8wN_HET_R3_count_summary$median_count)
rownames(X8wN_HET_R3_count_matrix) <- X8wN_HET_R3_count_summary$external_gene_name
colnames(X8wN_HET_R3_count_matrix) <- "8wN_HET_R3"


#### 8wN HET R4
X8wN_HET_R4 <- read.delim("output/featurecounts_hg38/8wN_HET_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R4$Geneid <- gsub("\\..*", "", X8wN_HET_R4$Geneid)
X8wN_HET_R4_geneSymbol <- merge(X8wN_HET_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R4_count_summary <- X8wN_HET_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R4_count_matrix <- as.matrix(X8wN_HET_R4_count_summary$median_count)
rownames(X8wN_HET_R4_count_matrix) <- X8wN_HET_R4_count_summary$external_gene_name
colnames(X8wN_HET_R4_count_matrix) <- "8wN_HET_R4"


#### 8wN KO R1
X8wN_KO_R1 <- read.delim("output/featurecounts_hg38/8wN_KO_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R1$Geneid <- gsub("\\..*", "", X8wN_KO_R1$Geneid)
X8wN_KO_R1_geneSymbol <- merge(X8wN_KO_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R1_count_summary <- X8wN_KO_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R1_count_matrix <- as.matrix(X8wN_KO_R1_count_summary$median_count)
rownames(X8wN_KO_R1_count_matrix) <- X8wN_KO_R1_count_summary$external_gene_name
colnames(X8wN_KO_R1_count_matrix) <- "8wN_KO_R1"



#### 8wN KO R2
X8wN_KO_R2 <- read.delim("output/featurecounts_hg38/8wN_KO_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R2$Geneid <- gsub("\\..*", "", X8wN_KO_R2$Geneid)
X8wN_KO_R2_geneSymbol <- merge(X8wN_KO_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R2_count_summary <- X8wN_KO_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R2_count_matrix <- as.matrix(X8wN_KO_R2_count_summary$median_count)
rownames(X8wN_KO_R2_count_matrix) <- X8wN_KO_R2_count_summary$external_gene_name
colnames(X8wN_KO_R2_count_matrix) <- "8wN_KO_R2"



#### 8wN KO R3
X8wN_KO_R3 <- read.delim("output/featurecounts_hg38/8wN_KO_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R3$Geneid <- gsub("\\..*", "", X8wN_KO_R3$Geneid)
X8wN_KO_R3_geneSymbol <- merge(X8wN_KO_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R3_count_summary <- X8wN_KO_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R3_count_matrix <- as.matrix(X8wN_KO_R3_count_summary$median_count)
rownames(X8wN_KO_R3_count_matrix) <- X8wN_KO_R3_count_summary$external_gene_name
colnames(X8wN_KO_R3_count_matrix) <- "8wN_KO_R3"



#### 8wN KO R4
X8wN_KO_R4 <- read.delim("output/featurecounts_hg38/8wN_KO_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R4$Geneid <- gsub("\\..*", "", X8wN_KO_R4$Geneid)
X8wN_KO_R4_geneSymbol <- merge(X8wN_KO_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R4_count_summary <- X8wN_KO_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R4_count_matrix <- as.matrix(X8wN_KO_R4_count_summary$median_count)
rownames(X8wN_KO_R4_count_matrix) <- X8wN_KO_R4_count_summary$external_gene_name
colnames(X8wN_KO_R4_count_matrix) <- "8wN_KO_R4"



# WT expressionSet
all_counts <- cbind(X8wN_WT_R1_count_matrix, X8wN_WT_R2_count_matrix, X8wN_WT_R3_count_matrix, X8wN_WT_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("WT", "WT","WT", "WT"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
WT_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)
### Then convert to expression matrix
WT.bulk.mtx = exprs(WT_8wN.bulkeset)


# HET expressionSet
all_counts <- cbind(X8wN_HET_R1_count_matrix, X8wN_HET_R2_count_matrix, X8wN_HET_R3_count_matrix, X8wN_HET_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("HET", "HET","HET", "HET"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
HET_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)  
### Then convert to expression matrix
HET.bulk.mtx = exprs(HET_8wN.bulkeset)

# KO expressionSet
all_counts <- cbind(X8wN_KO_R1_count_matrix, X8wN_KO_R2_count_matrix, X8wN_KO_R3_count_matrix, X8wN_KO_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("KO", "KO","KO", "KO"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
KO_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)  
### Then convert to expression matrix
KO.bulk.mtx = exprs(KO_8wN.bulkeset)



# all genotype together expressionSet
all_counts <- cbind(X8wN_WT_R1_count_matrix, X8wN_WT_R2_count_matrix, X8wN_WT_R3_count_matrix, X8wN_WT_R4_count_matrix, X8wN_HET_R1_count_matrix, X8wN_HET_R2_count_matrix, X8wN_HET_R3_count_matrix, X8wN_HET_R4_count_matrix, X8wN_KO_R1_count_matrix, X8wN_KO_R2_count_matrix, X8wN_KO_R3_count_matrix, X8wN_KO_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("WT", "WT","WT", "WT", "HET", "HET","HET", "HET", "KO", "KO","KO", "KO" ),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
all_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)
### Then convert to expression matrix
all.bulk.mtx = exprs(all_8wN.bulkeset)



# MuSiC1 deconvolution
## Estimate cell type proportions
Est.prop.WT = music_prop(bulk.mtx = WT.bulk.mtx, CHOOSE_full_dataset_srt.SingleCellExperiment, clusters = 'celltype_ctrl_transfer', samples = 'gRNA', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "INP", "IP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "OPC", "oRG", "RG","vRG"), verbose = TRUE)
Est.prop.HET = music_prop(bulk.mtx = HET.bulk.mtx, CHOOSE_full_dataset_srt.SingleCellExperiment, clusters = 'celltype_ctrl_transfer', samples = 'gRNA', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "INP", "IP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "OPC", "oRG", "RG","vRG"), verbose = TRUE)


## --> KO and all bug with SUBSCRIPTY OUT OF BOUND ERROR for no reason
Est.prop.KO = music_prop(bulk.mtx = KO.bulk.mtx, CHOOSE_full_dataset_srt.SingleCellExperiment, clusters = 'celltype_ctrl_transfer', samples = 'gRNA', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "INP", "IP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "OPC", "oRG", "RG","vRG"), verbose = TRUE)
# TO AVOID SUBSCRIPTY OUT OF BOUND ERROR:
### subset bulk features that are only present in scRNAseq assay
###### Assuming 'featureNames' can be used to retrieve the features from both ExpressionSets
common_features <- intersect(featureNames(KO_8wN.bulkeset), featureNames(CHOOSE_full_dataset_srt.sceset))
# Subset the WT_HET_8wN.bulkeset to only include common features
KO_8wN.bulkeset_subset <- KO_8wN.bulkeset[common_features, ]
KO.bulk.mtx = exprs(KO_8wN.bulkeset_subset)
Est.prop.KO = music_prop(bulk.mtx = KO.bulk.mtx, CHOOSE_full_dataset_srt.SingleCellExperiment, clusters = 'celltype_ctrl_transfer', samples = 'gRNA', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "INP", "IP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "OPC", "oRG", "RG","vRG"), verbose = TRUE)

common_features <- intersect(featureNames(all_8wN.bulkeset), featureNames(CHOOSE_full_dataset_srt.sceset))
# Subset the WT_HET_8wN.bulkeset to only include common features
all_8wN.bulkeset_subset <- all_8wN.bulkeset[common_features, ]
all.bulk.mtx = exprs(all_8wN.bulkeset_subset)
Est.prop.all = music_prop(bulk.mtx = all.bulk.mtx, CHOOSE_full_dataset_srt.SingleCellExperiment, clusters = 'celltype_ctrl_transfer', samples = 'gRNA', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "INP", "IP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "OPC", "oRG", "RG","vRG"), verbose = TRUE)

## write output
console_output <- capture.output(print(Est.prop.all$Est.prop.weighted))
writeLines(console_output, "output/deconv/MuSiC-Est.prop.weighted-all.txt")




# Check diff between estimation method (MuSiC vs NNLS)
## Jitter plot of estimated cell type proportions
pdf("output/deconv/MuSiC-jitterPlot_Est.prop.weighted-WT.pdf", width=14, height=20)  
Jitter_Est(list(data.matrix(Est.prop.WT$Est.prop.weighted),
                             data.matrix(Est.prop.WT$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
dev.off()
pdf("output/deconv/MuSiC-jitterPlot_Est.prop.weighted-HET.pdf", width=14, height=20)  
Jitter_Est(list(data.matrix(Est.prop.HET$Est.prop.weighted),
                             data.matrix(Est.prop.HET$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
dev.off()
pdf("output/deconv/MuSiC-jitterPlot_Est.prop.weighted-KO.pdf", width=14, height=20)  
Jitter_Est(list(data.matrix(Est.prop.KO$Est.prop.weighted),
                             data.matrix(Est.prop.KO$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
dev.off()

XXX plopt with cvolor
```

--> It work great; notably MuSiC method; produce similar result as braindeconv with more astrocyte HET and less for KO

Let's do deconvolution for all previous RNAseq point as QC:

```R


# library
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") 
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")
library("MuSiC2")
library("MuSiC")
library("biomaRt")



# scRNASeq import all
CHOOSE_full_dataset_srt <- readRDS(file = "output/deconv/CHOOSE_full_dataset_srt.rds")
DefaultAssay(CHOOSE_full_dataset_srt) <- "RNA" # 
## Transform scRNAseq data in SingleCellExperiment
CHOOSE_full_dataset_srt.SingleCellExperiment <- as.SingleCellExperiment(CHOOSE_full_dataset_srt, assay = "RNA")
## Transform scRNAseq data in SingleCellExperiment (FOR the out of bound error)
CHOOSE_full_dataset_srt.sceset = ExpressionSet(assayData = as.matrix(GetAssayData(CHOOSE_full_dataset_srt)), phenoData =  new("AnnotatedDataFrame",CHOOSE_full_dataset_srt@meta.data))


# Import raw RNAseq read counts
### code for ensembl to genesymbol conv
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <-  useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "uswest") # uswest useast asia
ensembl <-  useEnsembl(biomart = 'ensembl', 
                       dataset = 'hsapiens_gene_ensembl',
                       version = 110)
##                       
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
genes <- X8wN_WT_R1$Geneid
#### Get the mapping from Ensembl ID to gene symbol
genes_mapped  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = genes,
                      mart = ensembl)
####

### import featureCounts output
## 2dN
#### 2dN WT R1
X2dN_WT_R1 <- read.delim("output/featurecounts_hg38/2dN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_WT_R1$Geneid <- gsub("\\..*", "", X2dN_WT_R1$Geneid)
X2dN_WT_R1_geneSymbol <- merge(X2dN_WT_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_WT_R1_Aligned.sortedByCoord.out.bam')
## Calculate median as gene dupplicated with the conversino!
X2dN_WT_R1_count_summary <- X2dN_WT_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_WT_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_WT_R1_count_matrix <- as.matrix(X2dN_WT_R1_count_summary$median_count)
rownames(X2dN_WT_R1_count_matrix) <- X2dN_WT_R1_count_summary$external_gene_name
colnames(X2dN_WT_R1_count_matrix) <- "2dN_WT_R1"



#### 2dN WT R2
X2dN_WT_R2 <- read.delim("output/featurecounts_hg38/2dN_WT_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_WT_R2$Geneid <- gsub("\\..*", "", X2dN_WT_R2$Geneid)
X2dN_WT_R2_geneSymbol <- merge(X2dN_WT_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_WT_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_WT_R2_count_summary <- X2dN_WT_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_WT_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_WT_R2_count_matrix <- as.matrix(X2dN_WT_R2_count_summary$median_count)
rownames(X2dN_WT_R2_count_matrix) <- X2dN_WT_R2_count_summary$external_gene_name
colnames(X2dN_WT_R2_count_matrix) <- "2dN_WT_R2"


#### 2dN WT R3
X2dN_WT_R3 <- read.delim("output/featurecounts_hg38/2dN_WT_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_WT_R3$Geneid <- gsub("\\..*", "", X2dN_WT_R3$Geneid)
X2dN_WT_R3_geneSymbol <- merge(X2dN_WT_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_WT_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_WT_R3_count_summary <- X2dN_WT_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_WT_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_WT_R3_count_matrix <- as.matrix(X2dN_WT_R3_count_summary$median_count)
rownames(X2dN_WT_R3_count_matrix) <- X2dN_WT_R3_count_summary$external_gene_name
colnames(X2dN_WT_R3_count_matrix) <- "2dN_WT_R3"


#### 2dN WT R4
X2dN_WT_R4 <- read.delim("output/featurecounts_hg38/2dN_WT_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_WT_R4$Geneid <- gsub("\\..*", "", X2dN_WT_R4$Geneid)
X2dN_WT_R4_geneSymbol <- merge(X2dN_WT_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_WT_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_WT_R4_count_summary <- X2dN_WT_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_WT_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_WT_R4_count_matrix <- as.matrix(X2dN_WT_R4_count_summary$median_count)
rownames(X2dN_WT_R4_count_matrix) <- X2dN_WT_R4_count_summary$external_gene_name
colnames(X2dN_WT_R4_count_matrix) <- "2dN_WT_R4"


#### 2dN HET R1
X2dN_HET_R1 <- read.delim("output/featurecounts_hg38/2dN_HET_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_HET_R1$Geneid <- gsub("\\..*", "", X2dN_HET_R1$Geneid)
X2dN_HET_R1_geneSymbol <- merge(X2dN_HET_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_HET_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_HET_R1_count_summary <- X2dN_HET_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_HET_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_HET_R1_count_matrix <- as.matrix(X2dN_HET_R1_count_summary$median_count)
rownames(X2dN_HET_R1_count_matrix) <- X2dN_HET_R1_count_summary$external_gene_name
colnames(X2dN_HET_R1_count_matrix) <- "2dN_HET_R1"



#### 2dN HET R2
X2dN_HET_R2 <- read.delim("output/featurecounts_hg38/2dN_HET_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_HET_R2$Geneid <- gsub("\\..*", "", X2dN_HET_R2$Geneid)
X2dN_HET_R2_geneSymbol <- merge(X2dN_HET_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_HET_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_HET_R2_count_summary <- X2dN_HET_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_HET_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_HET_R2_count_matrix <- as.matrix(X2dN_HET_R2_count_summary$median_count)
rownames(X2dN_HET_R2_count_matrix) <- X2dN_HET_R2_count_summary$external_gene_name
colnames(X2dN_HET_R2_count_matrix) <- "2dN_HET_R2"


#### 2dN HET R3
X2dN_HET_R3 <- read.delim("output/featurecounts_hg38/2dN_HET_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_HET_R3$Geneid <- gsub("\\..*", "", X2dN_HET_R3$Geneid)
X2dN_HET_R3_geneSymbol <- merge(X2dN_HET_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_HET_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_HET_R3_count_summary <- X2dN_HET_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_HET_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_HET_R3_count_matrix <- as.matrix(X2dN_HET_R3_count_summary$median_count)
rownames(X2dN_HET_R3_count_matrix) <- X2dN_HET_R3_count_summary$external_gene_name
colnames(X2dN_HET_R3_count_matrix) <- "2dN_HET_R3"


#### 2dN HET R4
X2dN_HET_R4 <- read.delim("output/featurecounts_hg38/2dN_HET_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_HET_R4$Geneid <- gsub("\\..*", "", X2dN_HET_R4$Geneid)
X2dN_HET_R4_geneSymbol <- merge(X2dN_HET_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_HET_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_HET_R4_count_summary <- X2dN_HET_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_HET_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_HET_R4_count_matrix <- as.matrix(X2dN_HET_R4_count_summary$median_count)
rownames(X2dN_HET_R4_count_matrix) <- X2dN_HET_R4_count_summary$external_gene_name
colnames(X2dN_HET_R4_count_matrix) <- "2dN_HET_R4"


#### 2dN KO R1
X2dN_KO_R1 <- read.delim("output/featurecounts_hg38/2dN_KO_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_KO_R1$Geneid <- gsub("\\..*", "", X2dN_KO_R1$Geneid)
X2dN_KO_R1_geneSymbol <- merge(X2dN_KO_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_KO_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_KO_R1_count_summary <- X2dN_KO_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_KO_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_KO_R1_count_matrix <- as.matrix(X2dN_KO_R1_count_summary$median_count)
rownames(X2dN_KO_R1_count_matrix) <- X2dN_KO_R1_count_summary$external_gene_name
colnames(X2dN_KO_R1_count_matrix) <- "2dN_KO_R1"



#### 2dN KO R2
X2dN_KO_R2 <- read.delim("output/featurecounts_hg38/2dN_KO_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_KO_R2$Geneid <- gsub("\\..*", "", X2dN_KO_R2$Geneid)
X2dN_KO_R2_geneSymbol <- merge(X2dN_KO_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_KO_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_KO_R2_count_summary <- X2dN_KO_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_KO_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_KO_R2_count_matrix <- as.matrix(X2dN_KO_R2_count_summary$median_count)
rownames(X2dN_KO_R2_count_matrix) <- X2dN_KO_R2_count_summary$external_gene_name
colnames(X2dN_KO_R2_count_matrix) <- "2dN_KO_R2"



#### 2dN KO R3
X2dN_KO_R3 <- read.delim("output/featurecounts_hg38/2dN_KO_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_KO_R3$Geneid <- gsub("\\..*", "", X2dN_KO_R3$Geneid)
X2dN_KO_R3_geneSymbol <- merge(X2dN_KO_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_KO_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_KO_R3_count_summary <- X2dN_KO_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_KO_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_KO_R3_count_matrix <- as.matrix(X2dN_KO_R3_count_summary$median_count)
rownames(X2dN_KO_R3_count_matrix) <- X2dN_KO_R3_count_summary$external_gene_name
colnames(X2dN_KO_R3_count_matrix) <- "2dN_KO_R3"



#### 2dN KO R4
X2dN_KO_R4 <- read.delim("output/featurecounts_hg38/2dN_KO_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X2dN_KO_R4$Geneid <- gsub("\\..*", "", X2dN_KO_R4$Geneid)
X2dN_KO_R4_geneSymbol <- merge(X2dN_KO_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.2dN_KO_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X2dN_KO_R4_count_summary <- X2dN_KO_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.2dN_KO_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X2dN_KO_R4_count_matrix <- as.matrix(X2dN_KO_R4_count_summary$median_count)
rownames(X2dN_KO_R4_count_matrix) <- X2dN_KO_R4_count_summary$external_gene_name
colnames(X2dN_KO_R4_count_matrix) <- "2dN_KO_R4"



## NPC
#### NPC WT R1
XNPC_WT_R1 <- read.delim("output/featurecounts_hg38/NPC_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_WT_R1$Geneid <- gsub("\\..*", "", XNPC_WT_R1$Geneid)
XNPC_WT_R1_geneSymbol <- merge(XNPC_WT_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_WT_R1_Aligned.sortedByCoord.out.bam')
## Calculate median as gene dupplicated with the conversino!
XNPC_WT_R1_count_summary <- XNPC_WT_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_WT_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_WT_R1_count_matrix <- as.matrix(XNPC_WT_R1_count_summary$median_count)
rownames(XNPC_WT_R1_count_matrix) <- XNPC_WT_R1_count_summary$external_gene_name
colnames(XNPC_WT_R1_count_matrix) <- "NPC_WT_R1"



#### NPC WT R2
XNPC_WT_R2 <- read.delim("output/featurecounts_hg38/NPC_WT_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_WT_R2$Geneid <- gsub("\\..*", "", XNPC_WT_R2$Geneid)
XNPC_WT_R2_geneSymbol <- merge(XNPC_WT_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_WT_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_WT_R2_count_summary <- XNPC_WT_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_WT_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_WT_R2_count_matrix <- as.matrix(XNPC_WT_R2_count_summary$median_count)
rownames(XNPC_WT_R2_count_matrix) <- XNPC_WT_R2_count_summary$external_gene_name
colnames(XNPC_WT_R2_count_matrix) <- "NPC_WT_R2"


#### NPC WT R3
XNPC_WT_R3 <- read.delim("output/featurecounts_hg38/NPC_WT_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_WT_R3$Geneid <- gsub("\\..*", "", XNPC_WT_R3$Geneid)
XNPC_WT_R3_geneSymbol <- merge(XNPC_WT_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_WT_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_WT_R3_count_summary <- XNPC_WT_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_WT_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_WT_R3_count_matrix <- as.matrix(XNPC_WT_R3_count_summary$median_count)
rownames(XNPC_WT_R3_count_matrix) <- XNPC_WT_R3_count_summary$external_gene_name
colnames(XNPC_WT_R3_count_matrix) <- "NPC_WT_R3"


#### NPC WT R4
XNPC_WT_R4 <- read.delim("output/featurecounts_hg38/NPC_WT_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_WT_R4$Geneid <- gsub("\\..*", "", XNPC_WT_R4$Geneid)
XNPC_WT_R4_geneSymbol <- merge(XNPC_WT_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_WT_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_WT_R4_count_summary <- XNPC_WT_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_WT_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_WT_R4_count_matrix <- as.matrix(XNPC_WT_R4_count_summary$median_count)
rownames(XNPC_WT_R4_count_matrix) <- XNPC_WT_R4_count_summary$external_gene_name
colnames(XNPC_WT_R4_count_matrix) <- "NPC_WT_R4"


#### NPC HET R1
XNPC_HET_R1 <- read.delim("output/featurecounts_hg38/NPC_HET_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_HET_R1$Geneid <- gsub("\\..*", "", XNPC_HET_R1$Geneid)
XNPC_HET_R1_geneSymbol <- merge(XNPC_HET_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_HET_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_HET_R1_count_summary <- XNPC_HET_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_HET_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_HET_R1_count_matrix <- as.matrix(XNPC_HET_R1_count_summary$median_count)
rownames(XNPC_HET_R1_count_matrix) <- XNPC_HET_R1_count_summary$external_gene_name
colnames(XNPC_HET_R1_count_matrix) <- "NPC_HET_R1"



#### NPC HET R2
XNPC_HET_R2 <- read.delim("output/featurecounts_hg38/NPC_HET_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_HET_R2$Geneid <- gsub("\\..*", "", XNPC_HET_R2$Geneid)
XNPC_HET_R2_geneSymbol <- merge(XNPC_HET_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_HET_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_HET_R2_count_summary <- XNPC_HET_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_HET_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_HET_R2_count_matrix <- as.matrix(XNPC_HET_R2_count_summary$median_count)
rownames(XNPC_HET_R2_count_matrix) <- XNPC_HET_R2_count_summary$external_gene_name
colnames(XNPC_HET_R2_count_matrix) <- "NPC_HET_R2"


#### NPC HET R3
XNPC_HET_R3 <- read.delim("output/featurecounts_hg38/NPC_HET_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_HET_R3$Geneid <- gsub("\\..*", "", XNPC_HET_R3$Geneid)
XNPC_HET_R3_geneSymbol <- merge(XNPC_HET_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_HET_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_HET_R3_count_summary <- XNPC_HET_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_HET_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_HET_R3_count_matrix <- as.matrix(XNPC_HET_R3_count_summary$median_count)
rownames(XNPC_HET_R3_count_matrix) <- XNPC_HET_R3_count_summary$external_gene_name
colnames(XNPC_HET_R3_count_matrix) <- "NPC_HET_R3"


#### NPC HET R4
XNPC_HET_R4 <- read.delim("output/featurecounts_hg38/NPC_HET_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_HET_R4$Geneid <- gsub("\\..*", "", XNPC_HET_R4$Geneid)
XNPC_HET_R4_geneSymbol <- merge(XNPC_HET_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_HET_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_HET_R4_count_summary <- XNPC_HET_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_HET_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_HET_R4_count_matrix <- as.matrix(XNPC_HET_R4_count_summary$median_count)
rownames(XNPC_HET_R4_count_matrix) <- XNPC_HET_R4_count_summary$external_gene_name
colnames(XNPC_HET_R4_count_matrix) <- "NPC_HET_R4"


#### NPC KO R1
XNPC_KO_R1 <- read.delim("output/featurecounts_hg38/NPC_KO_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_KO_R1$Geneid <- gsub("\\..*", "", XNPC_KO_R1$Geneid)
XNPC_KO_R1_geneSymbol <- merge(XNPC_KO_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_KO_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_KO_R1_count_summary <- XNPC_KO_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_KO_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_KO_R1_count_matrix <- as.matrix(XNPC_KO_R1_count_summary$median_count)
rownames(XNPC_KO_R1_count_matrix) <- XNPC_KO_R1_count_summary$external_gene_name
colnames(XNPC_KO_R1_count_matrix) <- "NPC_KO_R1"



#### NPC KO R2
XNPC_KO_R2 <- read.delim("output/featurecounts_hg38/NPC_KO_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_KO_R2$Geneid <- gsub("\\..*", "", XNPC_KO_R2$Geneid)
XNPC_KO_R2_geneSymbol <- merge(XNPC_KO_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_KO_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_KO_R2_count_summary <- XNPC_KO_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_KO_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_KO_R2_count_matrix <- as.matrix(XNPC_KO_R2_count_summary$median_count)
rownames(XNPC_KO_R2_count_matrix) <- XNPC_KO_R2_count_summary$external_gene_name
colnames(XNPC_KO_R2_count_matrix) <- "NPC_KO_R2"



#### NPC KO R3
XNPC_KO_R3 <- read.delim("output/featurecounts_hg38/NPC_KO_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_KO_R3$Geneid <- gsub("\\..*", "", XNPC_KO_R3$Geneid)
XNPC_KO_R3_geneSymbol <- merge(XNPC_KO_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_KO_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_KO_R3_count_summary <- XNPC_KO_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_KO_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_KO_R3_count_matrix <- as.matrix(XNPC_KO_R3_count_summary$median_count)
rownames(XNPC_KO_R3_count_matrix) <- XNPC_KO_R3_count_summary$external_gene_name
colnames(XNPC_KO_R3_count_matrix) <- "NPC_KO_R3"



#### NPC KO R4
XNPC_KO_R4 <- read.delim("output/featurecounts_hg38/NPC_KO_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
XNPC_KO_R4$Geneid <- gsub("\\..*", "", XNPC_KO_R4$Geneid)
XNPC_KO_R4_geneSymbol <- merge(XNPC_KO_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.NPC_KO_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
XNPC_KO_R4_count_summary <- XNPC_KO_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.NPC_KO_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
XNPC_KO_R4_count_matrix <- as.matrix(XNPC_KO_R4_count_summary$median_count)
rownames(XNPC_KO_R4_count_matrix) <- XNPC_KO_R4_count_summary$external_gene_name
colnames(XNPC_KO_R4_count_matrix) <- "NPC_KO_R4"



## 4wN
#### 4wN WT R1
X4wN_WT_R1 <- read.delim("output/featurecounts_hg38/4wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_WT_R1$Geneid <- gsub("\\..*", "", X4wN_WT_R1$Geneid)
X4wN_WT_R1_geneSymbol <- merge(X4wN_WT_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_WT_R1_Aligned.sortedByCoord.out.bam')
## Calculate median as gene dupplicated with the conversino!
X4wN_WT_R1_count_summary <- X4wN_WT_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_WT_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_WT_R1_count_matrix <- as.matrix(X4wN_WT_R1_count_summary$median_count)
rownames(X4wN_WT_R1_count_matrix) <- X4wN_WT_R1_count_summary$external_gene_name
colnames(X4wN_WT_R1_count_matrix) <- "4wN_WT_R1"



#### 4wN WT R2
X4wN_WT_R2 <- read.delim("output/featurecounts_hg38/4wN_WT_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_WT_R2$Geneid <- gsub("\\..*", "", X4wN_WT_R2$Geneid)
X4wN_WT_R2_geneSymbol <- merge(X4wN_WT_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_WT_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_WT_R2_count_summary <- X4wN_WT_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_WT_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_WT_R2_count_matrix <- as.matrix(X4wN_WT_R2_count_summary$median_count)
rownames(X4wN_WT_R2_count_matrix) <- X4wN_WT_R2_count_summary$external_gene_name
colnames(X4wN_WT_R2_count_matrix) <- "4wN_WT_R2"


#### 4wN WT R3
X4wN_WT_R3 <- read.delim("output/featurecounts_hg38/4wN_WT_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_WT_R3$Geneid <- gsub("\\..*", "", X4wN_WT_R3$Geneid)
X4wN_WT_R3_geneSymbol <- merge(X4wN_WT_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_WT_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_WT_R3_count_summary <- X4wN_WT_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_WT_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_WT_R3_count_matrix <- as.matrix(X4wN_WT_R3_count_summary$median_count)
rownames(X4wN_WT_R3_count_matrix) <- X4wN_WT_R3_count_summary$external_gene_name
colnames(X4wN_WT_R3_count_matrix) <- "4wN_WT_R3"


#### 4wN WT R4
X4wN_WT_R4 <- read.delim("output/featurecounts_hg38/4wN_WT_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_WT_R4$Geneid <- gsub("\\..*", "", X4wN_WT_R4$Geneid)
X4wN_WT_R4_geneSymbol <- merge(X4wN_WT_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_WT_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_WT_R4_count_summary <- X4wN_WT_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_WT_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_WT_R4_count_matrix <- as.matrix(X4wN_WT_R4_count_summary$median_count)
rownames(X4wN_WT_R4_count_matrix) <- X4wN_WT_R4_count_summary$external_gene_name
colnames(X4wN_WT_R4_count_matrix) <- "4wN_WT_R4"


#### 4wN HET R1
X4wN_HET_R1 <- read.delim("output/featurecounts_hg38/4wN_HET_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_HET_R1$Geneid <- gsub("\\..*", "", X4wN_HET_R1$Geneid)
X4wN_HET_R1_geneSymbol <- merge(X4wN_HET_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_HET_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_HET_R1_count_summary <- X4wN_HET_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_HET_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_HET_R1_count_matrix <- as.matrix(X4wN_HET_R1_count_summary$median_count)
rownames(X4wN_HET_R1_count_matrix) <- X4wN_HET_R1_count_summary$external_gene_name
colnames(X4wN_HET_R1_count_matrix) <- "4wN_HET_R1"



#### 4wN HET R2
X4wN_HET_R2 <- read.delim("output/featurecounts_hg38/4wN_HET_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_HET_R2$Geneid <- gsub("\\..*", "", X4wN_HET_R2$Geneid)
X4wN_HET_R2_geneSymbol <- merge(X4wN_HET_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_HET_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_HET_R2_count_summary <- X4wN_HET_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_HET_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_HET_R2_count_matrix <- as.matrix(X4wN_HET_R2_count_summary$median_count)
rownames(X4wN_HET_R2_count_matrix) <- X4wN_HET_R2_count_summary$external_gene_name
colnames(X4wN_HET_R2_count_matrix) <- "4wN_HET_R2"


#### 4wN HET R3
X4wN_HET_R3 <- read.delim("output/featurecounts_hg38/4wN_HET_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_HET_R3$Geneid <- gsub("\\..*", "", X4wN_HET_R3$Geneid)
X4wN_HET_R3_geneSymbol <- merge(X4wN_HET_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_HET_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_HET_R3_count_summary <- X4wN_HET_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_HET_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_HET_R3_count_matrix <- as.matrix(X4wN_HET_R3_count_summary$median_count)
rownames(X4wN_HET_R3_count_matrix) <- X4wN_HET_R3_count_summary$external_gene_name
colnames(X4wN_HET_R3_count_matrix) <- "4wN_HET_R3"


#### 4wN HET R4
X4wN_HET_R4 <- read.delim("output/featurecounts_hg38/4wN_HET_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_HET_R4$Geneid <- gsub("\\..*", "", X4wN_HET_R4$Geneid)
X4wN_HET_R4_geneSymbol <- merge(X4wN_HET_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_HET_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_HET_R4_count_summary <- X4wN_HET_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_HET_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_HET_R4_count_matrix <- as.matrix(X4wN_HET_R4_count_summary$median_count)
rownames(X4wN_HET_R4_count_matrix) <- X4wN_HET_R4_count_summary$external_gene_name
colnames(X4wN_HET_R4_count_matrix) <- "4wN_HET_R4"


#### 4wN KO R1
X4wN_KO_R1 <- read.delim("output/featurecounts_hg38/4wN_KO_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_KO_R1$Geneid <- gsub("\\..*", "", X4wN_KO_R1$Geneid)
X4wN_KO_R1_geneSymbol <- merge(X4wN_KO_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_KO_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_KO_R1_count_summary <- X4wN_KO_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_KO_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_KO_R1_count_matrix <- as.matrix(X4wN_KO_R1_count_summary$median_count)
rownames(X4wN_KO_R1_count_matrix) <- X4wN_KO_R1_count_summary$external_gene_name
colnames(X4wN_KO_R1_count_matrix) <- "4wN_KO_R1"



#### 4wN KO R2
X4wN_KO_R2 <- read.delim("output/featurecounts_hg38/4wN_KO_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_KO_R2$Geneid <- gsub("\\..*", "", X4wN_KO_R2$Geneid)
X4wN_KO_R2_geneSymbol <- merge(X4wN_KO_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_KO_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_KO_R2_count_summary <- X4wN_KO_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_KO_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_KO_R2_count_matrix <- as.matrix(X4wN_KO_R2_count_summary$median_count)
rownames(X4wN_KO_R2_count_matrix) <- X4wN_KO_R2_count_summary$external_gene_name
colnames(X4wN_KO_R2_count_matrix) <- "4wN_KO_R2"



#### 4wN KO R3
X4wN_KO_R3 <- read.delim("output/featurecounts_hg38/4wN_KO_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_KO_R3$Geneid <- gsub("\\..*", "", X4wN_KO_R3$Geneid)
X4wN_KO_R3_geneSymbol <- merge(X4wN_KO_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_KO_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_KO_R3_count_summary <- X4wN_KO_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_KO_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_KO_R3_count_matrix <- as.matrix(X4wN_KO_R3_count_summary$median_count)
rownames(X4wN_KO_R3_count_matrix) <- X4wN_KO_R3_count_summary$external_gene_name
colnames(X4wN_KO_R3_count_matrix) <- "4wN_KO_R3"



#### 4wN KO R4
X4wN_KO_R4 <- read.delim("output/featurecounts_hg38/4wN_KO_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X4wN_KO_R4$Geneid <- gsub("\\..*", "", X4wN_KO_R4$Geneid)
X4wN_KO_R4_geneSymbol <- merge(X4wN_KO_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.4wN_KO_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X4wN_KO_R4_count_summary <- X4wN_KO_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.4wN_KO_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X4wN_KO_R4_count_matrix <- as.matrix(X4wN_KO_R4_count_summary$median_count)
rownames(X4wN_KO_R4_count_matrix) <- X4wN_KO_R4_count_summary$external_gene_name
colnames(X4wN_KO_R4_count_matrix) <- "4wN_KO_R4"






# NPC, 2dN and 4wN expressionSet
all_counts <- cbind(XNPC_WT_R1_count_matrix, XNPC_WT_R2_count_matrix, XNPC_WT_R3_count_matrix, 
XNPC_HET_R1_count_matrix, XNPC_HET_R2_count_matrix, XNPC_HET_R3_count_matrix,
XNPC_KO_R1_count_matrix, XNPC_KO_R2_count_matrix, XNPC_KO_R3_count_matrix,
X2dN_WT_R1_count_matrix, X2dN_WT_R2_count_matrix, X2dN_WT_R3_count_matrix, 
X2dN_HET_R1_count_matrix, X2dN_HET_R2_count_matrix, X2dN_HET_R3_count_matrix,
X2dN_KO_R1_count_matrix, X2dN_KO_R2_count_matrix, X2dN_KO_R3_count_matrix,
X4wN_WT_R1_count_matrix, X4wN_WT_R2_count_matrix,
X4wN_HET_R1_count_matrix, X4wN_HET_R2_count_matrix, X4wN_HET_R3_count_matrix, X4wN_HET_R4_count_matrix,
X4wN_KO_R1_count_matrix, X4wN_KO_R2_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
all_NPC_2dN_4wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)
### Then convert to expression matrix
all_NPC_2dN_4wN.bulk.mtx = exprs(all_NPC_2dN_4wN.bulkeset)



# MuSiC1 deconvolution

## Estimate cell type proportions
# TO AVOID SUBSCRIPTY OUT OF BOUND ERROR:
### subset bulk features that are only present in scRNAseq assay
###### Assuming 'featureNames' can be used to retrieve the features from both ExpressionSets
common_features <- intersect(featureNames(all_NPC_2dN_4wN.bulkeset), featureNames(CHOOSE_full_dataset_srt.sceset))
# Subset the WT_HET_8wN.bulkeset to only include common features
all_NPC_2dN_4wN.bulkeset_subset <- all_NPC_2dN_4wN.bulkeset[common_features, ]
all_NPC_2dN_4wN.bulk.mtx = exprs(all_NPC_2dN_4wN.bulkeset_subset)
Est.prop.all_NPC_2dN_4wN = music_prop(bulk.mtx = all_NPC_2dN_4wN.bulk.mtx, CHOOSE_full_dataset_srt.SingleCellExperiment, clusters = 'celltype_ctrl_transfer', samples = 'gRNA', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "INP", "IP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "OPC", "oRG", "RG","vRG"), verbose = TRUE)


## write output
console_output <- capture.output(print(Est.prop.all_NPC_2dN_4wN$Est.prop.weighted))
writeLines(console_output, "output/deconv/MuSiC-Est.prop.weighted-all_NPC_2dN_4wN.txt")





```


# Bisque RNAseq bulk deconvolution

[Paper](https://www.nature.com/articles/s41467-020-15816-6) and [github](https://github.com/cozygene/bisque) and [docs](https://cran.r-project.org/web/packages/BisqueRNA/BisqueRNA.pdf) and [tuto](https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html) and [tuto2](https://cozygene.r-universe.dev/BisqueRNA/doc/manual.html)


Install it in R within `scRNAseqV3`


```R
# library
library("BisqueRNA")
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") 
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")
library("biomaRt")








# Import raw RNAseq read counts
# code for ensembl to genesymbol conv
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <-  useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "uswest") # uswest useast asia
ensembl <-  useEnsembl(biomart = 'ensembl', 
                       dataset = 'hsapiens_gene_ensembl',
                       version = 110)
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
genes <- X8wN_WT_R1$Geneid
#### Get the mapping from Ensembl ID to gene symbol
genes_mapped  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = genes,
                      mart = ensembl)
####

### import featureCounts output
#### 8wN WT R1
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
X8wN_WT_R1_geneSymbol <- merge(X8wN_WT_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam')
## Calculate median as gene dupplicated with the conversino!
X8wN_WT_R1_count_summary <- X8wN_WT_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R1_count_matrix <- as.matrix(X8wN_WT_R1_count_summary$median_count)
rownames(X8wN_WT_R1_count_matrix) <- X8wN_WT_R1_count_summary$external_gene_name
colnames(X8wN_WT_R1_count_matrix) <- "8wN_WT_R1"



#### 8wN WT R2
X8wN_WT_R2 <- read.delim("output/featurecounts_hg38/8wN_WT_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R2$Geneid <- gsub("\\..*", "", X8wN_WT_R2$Geneid)
X8wN_WT_R2_geneSymbol <- merge(X8wN_WT_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R2_count_summary <- X8wN_WT_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R2_count_matrix <- as.matrix(X8wN_WT_R2_count_summary$median_count)
rownames(X8wN_WT_R2_count_matrix) <- X8wN_WT_R2_count_summary$external_gene_name
colnames(X8wN_WT_R2_count_matrix) <- "8wN_WT_R2"


#### 8wN WT R3
X8wN_WT_R3 <- read.delim("output/featurecounts_hg38/8wN_WT_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R3$Geneid <- gsub("\\..*", "", X8wN_WT_R3$Geneid)
X8wN_WT_R3_geneSymbol <- merge(X8wN_WT_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R3_count_summary <- X8wN_WT_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R3_count_matrix <- as.matrix(X8wN_WT_R3_count_summary$median_count)
rownames(X8wN_WT_R3_count_matrix) <- X8wN_WT_R3_count_summary$external_gene_name
colnames(X8wN_WT_R3_count_matrix) <- "8wN_WT_R3"


#### 8wN WT R4
X8wN_WT_R4 <- read.delim("output/featurecounts_hg38/8wN_WT_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R4$Geneid <- gsub("\\..*", "", X8wN_WT_R4$Geneid)
X8wN_WT_R4_geneSymbol <- merge(X8wN_WT_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R4_count_summary <- X8wN_WT_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R4_count_matrix <- as.matrix(X8wN_WT_R4_count_summary$median_count)
rownames(X8wN_WT_R4_count_matrix) <- X8wN_WT_R4_count_summary$external_gene_name
colnames(X8wN_WT_R4_count_matrix) <- "8wN_WT_R4"


#### 8wN HET R1
X8wN_HET_R1 <- read.delim("output/featurecounts_hg38/8wN_HET_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R1$Geneid <- gsub("\\..*", "", X8wN_HET_R1$Geneid)
X8wN_HET_R1_geneSymbol <- merge(X8wN_HET_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R1_count_summary <- X8wN_HET_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R1_count_matrix <- as.matrix(X8wN_HET_R1_count_summary$median_count)
rownames(X8wN_HET_R1_count_matrix) <- X8wN_HET_R1_count_summary$external_gene_name
colnames(X8wN_HET_R1_count_matrix) <- "8wN_HET_R1"



#### 8wN HET R2
X8wN_HET_R2 <- read.delim("output/featurecounts_hg38/8wN_HET_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R2$Geneid <- gsub("\\..*", "", X8wN_HET_R2$Geneid)
X8wN_HET_R2_geneSymbol <- merge(X8wN_HET_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R2_count_summary <- X8wN_HET_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R2_count_matrix <- as.matrix(X8wN_HET_R2_count_summary$median_count)
rownames(X8wN_HET_R2_count_matrix) <- X8wN_HET_R2_count_summary$external_gene_name
colnames(X8wN_HET_R2_count_matrix) <- "8wN_HET_R2"


#### 8wN HET R3
X8wN_HET_R3 <- read.delim("output/featurecounts_hg38/8wN_HET_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R3$Geneid <- gsub("\\..*", "", X8wN_HET_R3$Geneid)
X8wN_HET_R3_geneSymbol <- merge(X8wN_HET_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R3_count_summary <- X8wN_HET_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R3_count_matrix <- as.matrix(X8wN_HET_R3_count_summary$median_count)
rownames(X8wN_HET_R3_count_matrix) <- X8wN_HET_R3_count_summary$external_gene_name
colnames(X8wN_HET_R3_count_matrix) <- "8wN_HET_R3"


#### 8wN HET R4
X8wN_HET_R4 <- read.delim("output/featurecounts_hg38/8wN_HET_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R4$Geneid <- gsub("\\..*", "", X8wN_HET_R4$Geneid)
X8wN_HET_R4_geneSymbol <- merge(X8wN_HET_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R4_count_summary <- X8wN_HET_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R4_count_matrix <- as.matrix(X8wN_HET_R4_count_summary$median_count)
rownames(X8wN_HET_R4_count_matrix) <- X8wN_HET_R4_count_summary$external_gene_name
colnames(X8wN_HET_R4_count_matrix) <- "8wN_HET_R4"


#### 8wN KO R1
X8wN_KO_R1 <- read.delim("output/featurecounts_hg38/8wN_KO_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R1$Geneid <- gsub("\\..*", "", X8wN_KO_R1$Geneid)
X8wN_KO_R1_geneSymbol <- merge(X8wN_KO_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R1_count_summary <- X8wN_KO_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R1_count_matrix <- as.matrix(X8wN_KO_R1_count_summary$median_count)
rownames(X8wN_KO_R1_count_matrix) <- X8wN_KO_R1_count_summary$external_gene_name
colnames(X8wN_KO_R1_count_matrix) <- "8wN_KO_R1"



#### 8wN KO R2
X8wN_KO_R2 <- read.delim("output/featurecounts_hg38/8wN_KO_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R2$Geneid <- gsub("\\..*", "", X8wN_KO_R2$Geneid)
X8wN_KO_R2_geneSymbol <- merge(X8wN_KO_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R2_count_summary <- X8wN_KO_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R2_count_matrix <- as.matrix(X8wN_KO_R2_count_summary$median_count)
rownames(X8wN_KO_R2_count_matrix) <- X8wN_KO_R2_count_summary$external_gene_name
colnames(X8wN_KO_R2_count_matrix) <- "8wN_KO_R2"



#### 8wN KO R3
X8wN_KO_R3 <- read.delim("output/featurecounts_hg38/8wN_KO_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R3$Geneid <- gsub("\\..*", "", X8wN_KO_R3$Geneid)
X8wN_KO_R3_geneSymbol <- merge(X8wN_KO_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R3_count_summary <- X8wN_KO_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R3_count_matrix <- as.matrix(X8wN_KO_R3_count_summary$median_count)
rownames(X8wN_KO_R3_count_matrix) <- X8wN_KO_R3_count_summary$external_gene_name
colnames(X8wN_KO_R3_count_matrix) <- "8wN_KO_R3"



#### 8wN KO R4
X8wN_KO_R4 <- read.delim("output/featurecounts_hg38/8wN_KO_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R4$Geneid <- gsub("\\..*", "", X8wN_KO_R4$Geneid)
X8wN_KO_R4_geneSymbol <- merge(X8wN_KO_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R4_count_summary <- X8wN_KO_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R4_count_matrix <- as.matrix(X8wN_KO_R4_count_summary$median_count)
rownames(X8wN_KO_R4_count_matrix) <- X8wN_KO_R4_count_summary$external_gene_name
colnames(X8wN_KO_R4_count_matrix) <- "8wN_KO_R4"



# WT expressionSet
all_counts <- cbind(X8wN_WT_R1_count_matrix, X8wN_WT_R2_count_matrix, X8wN_WT_R3_count_matrix, X8wN_WT_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("WT", "WT","WT", "WT"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
WT_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)

# HET expressionSet
all_counts <- cbind(X8wN_HET_R1_count_matrix, X8wN_HET_R2_count_matrix, X8wN_HET_R3_count_matrix, X8wN_HET_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("HET", "HET","HET", "HET"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
HET_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)  

# KO expressionSet
all_counts <- cbind(X8wN_KO_R1_count_matrix, X8wN_KO_R2_count_matrix, X8wN_KO_R3_count_matrix, X8wN_KO_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("KO", "KO","KO", "KO"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
KO_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)  

# ALL expressionSet
all_counts <- cbind(X8wN_WT_R1_count_matrix, X8wN_WT_R2_count_matrix, X8wN_WT_R3_count_matrix, X8wN_WT_R4_count_matrix, X8wN_HET_R1_count_matrix, X8wN_HET_R2_count_matrix, X8wN_HET_R3_count_matrix, X8wN_HET_R4_count_matrix, X8wN_KO_R1_count_matrix, X8wN_KO_R2_count_matrix, X8wN_KO_R3_count_matrix, X8wN_KO_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("WT", "WT","WT", "WT", "HET", "HET","HET", "HET", "KO", "KO","KO", "KO"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
X8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)



# V1 
# scRNASeq import; filter control
## import seurat object
CHOOSE_CTRL_annot_srt <- readRDS(file = "output/deconv/CHOOSE_CTRL_annot_srt.rds")
DefaultAssay(CHOOSE_CTRL_annot_srt) <- "RNA" # 
pdf("output/deconv/UMAP_CHOOSE_CTRL_annot_srt.celltype_cl_coarse2.pdf", width=10, height=6)
DimPlot(CHOOSE_CTRL_annot_srt, reduction = "umap", group.by = "celltype_cl_coarse2", label=TRUE)
dev.off()
## isolate the cells of interest; control one
cells_to_keep <- WhichCells(CHOOSE_CTRL_annot_srt, expression = celltype_cl_coarse2 != "NA")
CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2 <- subset(CHOOSE_CTRL_annot_srt, cells = cells_to_keep)
pdf("output/deconv/UMAP_CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.pdf", width=10, height=6)
DimPlot(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2, reduction = "umap", group.by = "celltype_cl_coarse2", label=TRUE)
dev.off()
pdf("output/deconv/UMAP_CHOOSE_CTRL_annot_srt.subset_pseudotime_ranks.pdf", width=10, height=6)
FeaturePlot(
  object = CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2,
  features = 'pseudotime_ranks', 
  reduction = 'umap', 
  pt.size = 1, # Adjust point size if needed
  cols = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")) # Use reverse Spectral palette for color gradient
)
dev.off()
## Transform scRNAseq data in ExpressionSet Class
CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.sceset = ExpressionSet(assayData = as.matrix(GetAssayData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2)), phenoData =  new("AnnotatedDataFrame",CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2@meta.data))
### Add WT to all entries of our scRNAseq
pData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.sceset)$group <- rep("WT", nrow(pData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.sceset)))

#--> This method V1 fail as sneed several scRNAseq individuals


# V2 scRNASeq import all
CHOOSE_full_dataset_srt <- readRDS(file = "output/deconv/CHOOSE_full_dataset_srt.rds")
DefaultAssay(CHOOSE_full_dataset_srt) <- "RNA" # 
pdf("output/deconv/UMAP_CHOOSE_full_dataset_srt.celltype_cl_coarse2.pdf", width=10, height=6)
DimPlot(CHOOSE_full_dataset_srt, reduction = "umap", group.by = "celltype_cl_coarse2", label=TRUE)
dev.off()
pdf("output/deconv/UMAP_CHOOSE_full_dataset_srt.gRNA.pdf", width=10, height=6)
DimPlot(CHOOSE_full_dataset_srt, reduction = "umap", group.by = "gRNA", label=TRUE)
dev.off()
pdf("output/deconv/UMAP_CHOOSE_full_dataset_srt.celltype_cl_refined.pdf", width=10, height=6)
DimPlot(CHOOSE_full_dataset_srt, reduction = "umap", group.by = "celltype_cl_refined", label=TRUE)
dev.off()
pdf("output/deconv/UMAP_CHOOSE_full_dataset_srt.celltype_jf_refined.pdf", width=10, height=6)
DimPlot(CHOOSE_full_dataset_srt, reduction = "umap", group.by = "celltype_jf_refined", label=TRUE)
dev.off()
pdf("output/deconv/UMAP_CHOOSE_full_dataset_srt.celltype_ctrl_transfer.pdf", width=10, height=6)
DimPlot(CHOOSE_full_dataset_srt, reduction = "umap", group.by = "celltype_ctrl_transfer", label=TRUE)
dev.off()
pdf("output/deconv/UMAP_CHOOSE_full_dataset_srt.celltype_ctrl_transfer_coarse.pdf", width=10, height=6)
DimPlot(CHOOSE_full_dataset_srt, reduction = "umap", group.by = "celltype_ctrl_transfer_coarse", label=TRUE)
dev.off()
pdf("output/deconv/UMAP_CHOOSE_full_dataset_srt.celltype_jf_ctrl.pdf", width=10, height=6)
DimPlot(CHOOSE_full_dataset_srt, reduction = "umap", group.by = "celltype_jf_ctrl", label=TRUE)
dev.off()
pdf("output/deconv/UMAP_CHOOSE_full_dataset_srt.lineage.pdf", width=10, height=6)
DimPlot(CHOOSE_full_dataset_srt, reduction = "umap", group.by = "lineage", label=TRUE)
dev.off()

## celltype_ctrl_transfer is GOOD
## Transform scRNAseq data in ExpressionSet Class
CHOOSE_full_dataset_srt.sceset = ExpressionSet(assayData = as.matrix(GetAssayData(CHOOSE_full_dataset_srt)), phenoData =  new("AnnotatedDataFrame",CHOOSE_full_dataset_srt@meta.data))

# run reference based deconvolution _ celltype_ctrl_transfer
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = WT_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "celltype_ctrl_transfer", subject.names = "gRNA")
ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)

res.HET <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = HET_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "celltype_ctrl_transfer", subject.names = "gRNA")
ref.based.estimates.HET <- res.HET$bulk.props
knitr::kable(ref.based.estimates.HET, digits=2)

res.KO <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = KO_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "celltype_ctrl_transfer", subject.names = "gRNA")
ref.based.estimates.KO <- res.KO$bulk.props
knitr::kable(ref.based.estimates.KO, digits=2)

res.all <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = X8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "celltype_ctrl_transfer", subject.names = "gRNA")
ref.based.estimates.all <- res.all$bulk.props
knitr::kable(ref.based.estimates.all, digits=2)


## write output
console_output <- capture.output(print(knitr::kable(ref.based.estimates, digits=2)))
writeLines(console_output, "output/deconv/bisque-8wN_WT-ReferenceBasedDecomposition.txt")
console_output <- capture.output(print(knitr::kable(ref.based.estimates.HET, digits=2)))
writeLines(console_output, "output/deconv/bisque-8wN_HET-ReferenceBasedDecomposition.txt")
console_output <- capture.output(print(knitr::kable(ref.based.estimates.KO, digits=2)))
writeLines(console_output, "output/deconv/bisque-8wN_KO-ReferenceBasedDecomposition.txt")
console_output <- capture.output(print(knitr::kable(ref.based.estimates.all, digits=2)))
writeLines(console_output, "output/deconv/bisque-8wN_all-ReferenceBasedDecomposition.txt")




# run reference based deconvolution _ lineage
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = WT_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "lineage", subject.names = "gRNA")
ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)

res.HET <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = HET_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "lineage", subject.names = "gRNA")
ref.based.estimates.HET <- res.HET$bulk.props
knitr::kable(ref.based.estimates.HET, digits=2)

res.KO <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = KO_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "lineage", subject.names = "gRNA")
ref.based.estimates.KO <- res.KO$bulk.props
knitr::kable(ref.based.estimates.KO, digits=2)

res.all <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = X8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "lineage", subject.names = "gRNA")
ref.based.estimates.all <- res.all$bulk.props
knitr::kable(ref.based.estimates.all, digits=2)


## write output
console_output <- capture.output(print(knitr::kable(ref.based.estimates, digits=2)))
writeLines(console_output, "output/deconv/bisque-8wN_WT-ReferenceBasedDecomposition-lineage.txt")
console_output <- capture.output(print(knitr::kable(ref.based.estimates.HET, digits=2)))
writeLines(console_output, "output/deconv/bisque-8wN_HET-ReferenceBasedDecomposition-lineage.txt")
console_output <- capture.output(print(knitr::kable(ref.based.estimates.KO, digits=2)))
writeLines(console_output, "output/deconv/bisque-8wN_KO-ReferenceBasedDecomposition-lineage.txt")
console_output <- capture.output(print(knitr::kable(ref.based.estimates.all, digits=2)))
writeLines(console_output, "output/deconv/bisque-8wN_all-ReferenceBasedDecomposition-lineage.txt")


```

--> Fail as need several individuals scRNAseq data... Weird Give another try using the full dataset (with the 36 deleted genes)

--> celltype_ctrl_transfer is the best (all cells are included and annotated!)
--> It work but no significant difference detected... Still same tendency for Astrocyte tho!



# Granulator bullk deconvolution

[tuto](https://bioconductor.org/packages/release/bioc/vignettes/granulator/inst/doc/granulator.html)

Seems granulator needs a reference profile (Average gene expression within each cell types); seems to work with already processed data. So let's use the BisqueRNA tool to create **normalized bulk** and **signature profile**
- normalized bulk will be extracted from BisqueRNA output `res$transformed.bulk`
*--> The norm bulk are not good as contained negative value!! Instead let's use TPM*
- signature will be genrate using BisqueRNA with `GenerateSCReference`(https://cran.r-project.org/web/packages/BisqueRNA/BisqueRNA.pdf)

--> input files generated within `scRNAseqV3`; then granulator need a separate conda env... (as R install failed)



## Granulator input generation

```R
# library
library("BisqueRNA")
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") 
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")
library("biomaRt")


# Import raw RNAseq read counts
# code for ensembl to genesymbol conv
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <-  useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "uswest") # uswest useast asia
ensembl <-  useEnsembl(biomart = 'ensembl', 
                       dataset = 'hsapiens_gene_ensembl',
                       version = 110)
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
genes <- X8wN_WT_R1$Geneid
#### Get the mapping from Ensembl ID to gene symbol
genes_mapped  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = genes,
                      mart = ensembl)
####

### import featureCounts output
#### 8wN WT R1
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
X8wN_WT_R1_geneSymbol <- merge(X8wN_WT_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam')
## Calculate median as gene dupplicated with the conversino!
X8wN_WT_R1_count_summary <- X8wN_WT_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R1_count_matrix <- as.matrix(X8wN_WT_R1_count_summary$median_count)
rownames(X8wN_WT_R1_count_matrix) <- X8wN_WT_R1_count_summary$external_gene_name
colnames(X8wN_WT_R1_count_matrix) <- "8wN_WT_R1"



#### 8wN WT R2
X8wN_WT_R2 <- read.delim("output/featurecounts_hg38/8wN_WT_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R2$Geneid <- gsub("\\..*", "", X8wN_WT_R2$Geneid)
X8wN_WT_R2_geneSymbol <- merge(X8wN_WT_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R2_count_summary <- X8wN_WT_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R2_count_matrix <- as.matrix(X8wN_WT_R2_count_summary$median_count)
rownames(X8wN_WT_R2_count_matrix) <- X8wN_WT_R2_count_summary$external_gene_name
colnames(X8wN_WT_R2_count_matrix) <- "8wN_WT_R2"


#### 8wN WT R3
X8wN_WT_R3 <- read.delim("output/featurecounts_hg38/8wN_WT_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R3$Geneid <- gsub("\\..*", "", X8wN_WT_R3$Geneid)
X8wN_WT_R3_geneSymbol <- merge(X8wN_WT_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R3_count_summary <- X8wN_WT_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R3_count_matrix <- as.matrix(X8wN_WT_R3_count_summary$median_count)
rownames(X8wN_WT_R3_count_matrix) <- X8wN_WT_R3_count_summary$external_gene_name
colnames(X8wN_WT_R3_count_matrix) <- "8wN_WT_R3"


#### 8wN WT R4
X8wN_WT_R4 <- read.delim("output/featurecounts_hg38/8wN_WT_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R4$Geneid <- gsub("\\..*", "", X8wN_WT_R4$Geneid)
X8wN_WT_R4_geneSymbol <- merge(X8wN_WT_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R4_count_summary <- X8wN_WT_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R4_count_matrix <- as.matrix(X8wN_WT_R4_count_summary$median_count)
rownames(X8wN_WT_R4_count_matrix) <- X8wN_WT_R4_count_summary$external_gene_name
colnames(X8wN_WT_R4_count_matrix) <- "8wN_WT_R4"


#### 8wN HET R1
X8wN_HET_R1 <- read.delim("output/featurecounts_hg38/8wN_HET_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R1$Geneid <- gsub("\\..*", "", X8wN_HET_R1$Geneid)
X8wN_HET_R1_geneSymbol <- merge(X8wN_HET_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R1_count_summary <- X8wN_HET_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R1_count_matrix <- as.matrix(X8wN_HET_R1_count_summary$median_count)
rownames(X8wN_HET_R1_count_matrix) <- X8wN_HET_R1_count_summary$external_gene_name
colnames(X8wN_HET_R1_count_matrix) <- "8wN_HET_R1"



#### 8wN HET R2
X8wN_HET_R2 <- read.delim("output/featurecounts_hg38/8wN_HET_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R2$Geneid <- gsub("\\..*", "", X8wN_HET_R2$Geneid)
X8wN_HET_R2_geneSymbol <- merge(X8wN_HET_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R2_count_summary <- X8wN_HET_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R2_count_matrix <- as.matrix(X8wN_HET_R2_count_summary$median_count)
rownames(X8wN_HET_R2_count_matrix) <- X8wN_HET_R2_count_summary$external_gene_name
colnames(X8wN_HET_R2_count_matrix) <- "8wN_HET_R2"


#### 8wN HET R3
X8wN_HET_R3 <- read.delim("output/featurecounts_hg38/8wN_HET_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R3$Geneid <- gsub("\\..*", "", X8wN_HET_R3$Geneid)
X8wN_HET_R3_geneSymbol <- merge(X8wN_HET_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R3_count_summary <- X8wN_HET_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R3_count_matrix <- as.matrix(X8wN_HET_R3_count_summary$median_count)
rownames(X8wN_HET_R3_count_matrix) <- X8wN_HET_R3_count_summary$external_gene_name
colnames(X8wN_HET_R3_count_matrix) <- "8wN_HET_R3"


#### 8wN HET R4
X8wN_HET_R4 <- read.delim("output/featurecounts_hg38/8wN_HET_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R4$Geneid <- gsub("\\..*", "", X8wN_HET_R4$Geneid)
X8wN_HET_R4_geneSymbol <- merge(X8wN_HET_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R4_count_summary <- X8wN_HET_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R4_count_matrix <- as.matrix(X8wN_HET_R4_count_summary$median_count)
rownames(X8wN_HET_R4_count_matrix) <- X8wN_HET_R4_count_summary$external_gene_name
colnames(X8wN_HET_R4_count_matrix) <- "8wN_HET_R4"


#### 8wN KO R1
X8wN_KO_R1 <- read.delim("output/featurecounts_hg38/8wN_KO_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R1$Geneid <- gsub("\\..*", "", X8wN_KO_R1$Geneid)
X8wN_KO_R1_geneSymbol <- merge(X8wN_KO_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R1_count_summary <- X8wN_KO_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R1_count_matrix <- as.matrix(X8wN_KO_R1_count_summary$median_count)
rownames(X8wN_KO_R1_count_matrix) <- X8wN_KO_R1_count_summary$external_gene_name
colnames(X8wN_KO_R1_count_matrix) <- "8wN_KO_R1"



#### 8wN KO R2
X8wN_KO_R2 <- read.delim("output/featurecounts_hg38/8wN_KO_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R2$Geneid <- gsub("\\..*", "", X8wN_KO_R2$Geneid)
X8wN_KO_R2_geneSymbol <- merge(X8wN_KO_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R2_count_summary <- X8wN_KO_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R2_count_matrix <- as.matrix(X8wN_KO_R2_count_summary$median_count)
rownames(X8wN_KO_R2_count_matrix) <- X8wN_KO_R2_count_summary$external_gene_name
colnames(X8wN_KO_R2_count_matrix) <- "8wN_KO_R2"



#### 8wN KO R3
X8wN_KO_R3 <- read.delim("output/featurecounts_hg38/8wN_KO_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R3$Geneid <- gsub("\\..*", "", X8wN_KO_R3$Geneid)
X8wN_KO_R3_geneSymbol <- merge(X8wN_KO_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R3_count_summary <- X8wN_KO_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R3_count_matrix <- as.matrix(X8wN_KO_R3_count_summary$median_count)
rownames(X8wN_KO_R3_count_matrix) <- X8wN_KO_R3_count_summary$external_gene_name
colnames(X8wN_KO_R3_count_matrix) <- "8wN_KO_R3"



#### 8wN KO R4
X8wN_KO_R4 <- read.delim("output/featurecounts_hg38/8wN_KO_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R4$Geneid <- gsub("\\..*", "", X8wN_KO_R4$Geneid)
X8wN_KO_R4_geneSymbol <- merge(X8wN_KO_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R4_count_summary <- X8wN_KO_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R4_count_matrix <- as.matrix(X8wN_KO_R4_count_summary$median_count)
rownames(X8wN_KO_R4_count_matrix) <- X8wN_KO_R4_count_summary$external_gene_name
colnames(X8wN_KO_R4_count_matrix) <- "8wN_KO_R4"



# WT expressionSet
all_counts <- cbind(X8wN_WT_R1_count_matrix, X8wN_WT_R2_count_matrix, X8wN_WT_R3_count_matrix, X8wN_WT_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("WT", "WT","WT", "WT"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
WT_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)

# HET expressionSet
all_counts <- cbind(X8wN_HET_R1_count_matrix, X8wN_HET_R2_count_matrix, X8wN_HET_R3_count_matrix, X8wN_HET_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("HET", "HET","HET", "HET"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
HET_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)  

# KO expressionSet
all_counts <- cbind(X8wN_KO_R1_count_matrix, X8wN_KO_R2_count_matrix, X8wN_KO_R3_count_matrix, X8wN_KO_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("KO", "KO","KO", "KO"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
KO_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)  

# ALL expressionSet
all_counts <- cbind(X8wN_WT_R1_count_matrix, X8wN_WT_R2_count_matrix, X8wN_WT_R3_count_matrix, X8wN_WT_R4_count_matrix, X8wN_HET_R1_count_matrix, X8wN_HET_R2_count_matrix, X8wN_HET_R3_count_matrix, X8wN_HET_R4_count_matrix, X8wN_KO_R1_count_matrix, X8wN_KO_R2_count_matrix, X8wN_KO_R3_count_matrix, X8wN_KO_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("WT", "WT","WT", "WT", "HET", "HET","HET", "HET", "KO", "KO","KO", "KO"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
X8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)



# V2 scRNASeq import all
CHOOSE_full_dataset_srt <- readRDS(file = "output/deconv/CHOOSE_full_dataset_srt.rds")
DefaultAssay(CHOOSE_full_dataset_srt) <- "RNA" #


## celltype_ctrl_transfer is GOOD
## Transform scRNAseq data in ExpressionSet Class
CHOOSE_full_dataset_srt.sceset = ExpressionSet(assayData = as.matrix(GetAssayData(CHOOSE_full_dataset_srt)), phenoData =  new("AnnotatedDataFrame",CHOOSE_full_dataset_srt@meta.data))



# run reference based deconvolution _ celltype_ctrl_transfer
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = WT_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "celltype_ctrl_transfer", subject.names = "gRNA")


res.HET <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = HET_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "celltype_ctrl_transfer", subject.names = "gRNA")


res.KO <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = KO_8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "celltype_ctrl_transfer", subject.names = "gRNA")

res.all <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = X8wN.bulkeset, sc.eset = CHOOSE_full_dataset_srt.sceset, markers=NULL, use.overlap=FALSE, cell.types = "celltype_ctrl_transfer", subject.names = "gRNA")

# collect cpm (normalized bulk RNA seq)
WT_bulk_CPM = res$transformed.bulk
HET_bulk_CPM = res.HET$transformed.bulk
KO_bulk_CPM = res.KO$transformed.bulk
all_bulk_CPM = res.all$transformed.bulk

# generate scRNAseq signature

CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature = GenerateSCReference(sc.eset = CHOOSE_full_dataset_srt.sceset, cell.types = "celltype_ctrl_transfer")



## write output files

write.table(CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature, file = "output/deconv/CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(WT_bulk_CPM, file = "output/deconv/WT_bulk_CPM.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(HET_bulk_CPM, file = "output/deconv/HET_bulk_CPM.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(KO_bulk_CPM, file = "output/deconv/KO_bulk_CPM.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(all_bulk_CPM, file = "output/deconv/all_bulk_CPM.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)






# Generate tpm files with geneSymbol


tpm_all_sample <- read_csv("output/tpm_hg38/tpm_all_sample.txt") %>%
  dplyr::select(-'...1')

tpm_all_sample$Geneid <- gsub("\\..*", "", tpm_all_sample$Geneid)
tpm_all_sample_geneSymbol <- merge(tpm_all_sample, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na()


### Output table
write.table(tpm_all_sample_geneSymbol, file = "output/tpm_hg38/tpm_all_sample_geneSymbol.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


```

## Granulator deconvolution


Install granulator on R\4.3.1 like in the [tuto](https://bioconductor.org/packages/release/bioc/vignettes/granulator/inst/doc/granulator.html#installation)


--> Lot of fail installing it with Bioconda, even when using the R4.3.1 version from the tutorial, I ended up installing it through conda with `conda create -n granulator r-base=4.3.1 bioconductor-granulator -c conda-forge -c bioconda`




```bash
conda activate granulator
```



```R
# packages
library("granulator")
library("tidyverse")



# import files _ V1

CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature <- read.delim("output/deconv/CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature.txt",   header = TRUE,  row.names = 1)
WT_bulk_CPM <- read.delim("output/deconv/WT_bulk_CPM.txt",   header = TRUE,  row.names = 1)
HET_bulk_CPM <- read.delim("output/deconv/HET_bulk_CPM.txt",   header = TRUE,  row.names = 1)
KO_bulk_CPM <- read.delim("output/deconv/KO_bulk_CPM.txt",   header = TRUE,  row.names = 1)
all_bulk_CPM <- read.delim("output/deconv/all_bulk_CPM.txt",   header = TRUE,  row.names = 1)

# import tpm _ V2
all_TPM <- read.delim("output/tpm_hg38/tpm_all_sample_geneSymbol.txt",   header = TRUE,  row.names = 1)
## WT
WT_bulk_TPM <- all_TPM %>%
  dplyr::select(external_gene_name, X8wN_WT_R1, X8wN_WT_R2, X8wN_WT_R3, X8wN_WT_R4) %>%
  as.tibble()
WT_bulk_matrix <- as.matrix(WT_bulk_TPM[,-1])  # Exclude the first column before converting
rownames(WT_bulk_matrix) <- WT_bulk_TPM$external_gene_name
## HET
HET_bulk_TPM <- all_TPM %>%
  dplyr::select(external_gene_name, X8wN_HET_R1, X8wN_HET_R2, X8wN_HET_R3, X8wN_HET_R4) %>%
  as.tibble()
HET_bulk_matrix <- as.matrix(HET_bulk_TPM[,-1])  # Exclude the first column before converting
rownames(HET_bulk_matrix) <- HET_bulk_TPM$external_gene_name
## KO
KO_bulk_TPM <- all_TPM %>%
  dplyr::select(external_gene_name, X8wN_KO_R1, X8wN_KO_R2, X8wN_KO_R3, X8wN_KO_R4) %>%
  as.tibble()
KO_bulk_matrix <- as.matrix(KO_bulk_TPM[,-1])  # Exclude the first column before converting
rownames(KO_bulk_matrix) <- KO_bulk_TPM$external_gene_name



# run granulator
## WT
decon <- deconvolute(m = as.matrix(WT_bulk_matrix), sigMatrix = as.matrix(CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature) )
decon <- deconvolute(m = as.matrix(WT_bulk_matrix), sigMatrix = as.matrix(CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature) , methods = "dtangle")
## HET
decon <- deconvolute(m = as.matrix(HET_bulk_matrix), sigMatrix = as.matrix(CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature), methods = "dtangle")
## KO
decon <- deconvolute(m = as.matrix(KO_bulk_matrix), sigMatrix = as.matrix(CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature), methods = "dtangle")


## Save output
console_output <- capture.output(print(decon))
writeLines(console_output, "output/deconv/granulator_WT_dtangle.txt")
console_output <- capture.output(print(decon))
writeLines(console_output, "output/deconv/granulator_HET_dtangle.txt")
console_output <- capture.output(print(decon))
writeLines(console_output, "output/deconv/granulator_KO_dtangle.txt")

## granulator plot
pdf("output/deconv/granulator-plot_similarity_CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature.pdf", width=14, height=20)  
plot_similarity(sigMatrix=as.matrix(CHOOSE_full_dataset_srt.sceset.celltype_ctrl_transfer_signature))
dev.off()
pdf("output/deconv/granulator_WT-plot_deconvolute.pdf", width=14, height=20)  
plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = FALSE)
dev.off()


## save environemnt

save.image(file="output/deconv/granulator_WT.RData")
load("output/deconv/granulator_WT.RData")




```


--> only **dTangle** work; the other mthod fail; gave weird results




## decvonvolution data vizualization


--> MuSiC1 (`output/deconv/MuSiC-Est.prop.weighted-*.txt`) and Bisque outputs transferred in local to generate tidy txt file via Xcell = and transferred back to the cluster as `output/deconv/deconv_Lietal2023.txt` (proporiton for all time points/all genotpes)




```bash
conda activate deseq2
```



```R
library("tidyverse")

deconv_Lietal2023 = read.table("output/deconv/deconv_Lietal2023.txt", 
                                           header = TRUE) %>%
                               as_tibble()

deconv_Lietal2023_MuSiC = deconv_Lietal2023 %>% filter(method == "MuSiC") %>%
  dplyr::select(-method) %>%
  pivot_longer(cols = -sample, names_to = "tissue", values_to = "proportion") %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep ="_") %>%
  add_column(method = "MuSiC")
deconv_Lietal2023_bisque = deconv_Lietal2023 %>% filter(method == "bisque") %>%
  dplyr::select(-method) %>%
  pivot_longer(cols = -sample, names_to = "tissue", values_to = "proportion") %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep ="_") %>%
  add_column(method = "bisque")

deconv_Lietal2023_tidy = deconv_Lietal2023_MuSiC %>%
  bind_rows(deconv_Lietal2023_bisque)



# calculate the mean and standard error for each tissue-genotype combination
deconv_Lietal2023_stat <- deconv_Lietal2023_tidy %>%
  group_by(tissue, genotype,method) %>%
  summarize(
    mean_proportion = mean(proportion),
    se_proportion = sd(proportion) / sqrt(n()),
    .groups = 'drop'
  )


deconv_Lietal2023_stat$genotype <-
  factor(deconv_Lietal2023_stat$genotype,
         c("WT", "HET", "KO"))


# Create the plot
pdf("output/deconv/deconv_Lietal2023_barplot.pdf", width=18, height=10)  
ggplot(deconv_Lietal2023_stat, aes(x = tissue, y = mean_proportion, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(ymin = mean_proportion - se_proportion, ymax = mean_proportion + se_proportion),
    position = position_dodge(0.9),
    width = 0.25
  ) +
  scale_fill_manual(values = c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  facet_wrap(~method) +
  theme_bw() +
  labs(x = "Tissue", y = "Proportion", fill = "Genotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


```





