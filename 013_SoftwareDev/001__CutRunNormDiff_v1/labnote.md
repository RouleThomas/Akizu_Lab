# Project

Generate a tool that will take output from alignment (BAM files) and generate:
- Normalized bigwig files
- Coordinates of diff. bound regions using MACS2/DESEQ2
    - volcano plot of region gain and lost the mark (1 dot = 1 region)
    

The tool will be divided in two main commands:

***1) Generation of normalized bigwig files***
- *Input*: BAM files (uniquely aligned reads)
- Convert BAM → bigWig (bigWigToBedGraph)
- Convert bigWig → bedGraph (bigWigToBedGraph) 
- Remove blacklisted regions (bedtools intersect)
- Identify local signal maxima (Python)
- Calculate the 99th percentile of signal (Python)
- Apply scaling factor to blacklist-free bedGraph (Python)
- Generate normalized bigWig files (bedGraphToBigWig)

***2) Diff. binding***

**Peak-based mode**
- *Input*: Normalized bedGraph files + peaks in bed format
- Define experimental design (conditions / replicates)
- Identify consensus peaks across samples
- Quantify signal in peaks (deepTools computeMatrix)
- Perform differential binding analysis (DESeq2)



# Workflow

***1) Manage package installation dependencies***
- Make a new conda env and install EVERYTHING I will call
- SAVE environment.yml so that anyone can install it
- Ask someone to install it as test

***2) Build the software***


--> Potential tool name: *NORMDiff* (NORMalization & DIFFerential binding) or *NORMDB* (NORMalization & Differential Binding)

--> Potential magazine to submit: Bioinformatics (Oxford Academic), GigaSciences, PLOS Computational Biology, NAR Genomics and Bioinformatics (or less great: F1000Research)


# Package installation dependencies

Let's try to create a clean conda env with everything needed installed; this include the following:
- bigWigToBedGraph, BedGraphTobigWig
- bedtools
- MACS2
- deepTools
- DESEQ2



```bash
conda create -n normdb_v1 -c conda-forge -c bioconda -y \
  python=3.11 \
  samtools bedtools \
  deeptools \
  macs2 \
  ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph \
  r-base=4.3 \
  r-tidyverse \
  r-optparse r-data.table \
  bioconductor-deseq2 \
  bioconductor-edger \
  bioconductor-enhancedvolcano
```

FAIL, lets try using mamba

```bash
conda create -n mamba -c conda-forge -y mamba

# conda run = conda activate / run command / conda deactivate
conda run -n mamba mamba create -n normdb_v2 -c conda-forge -c bioconda -y \
  python=3.11 \
  samtools bedtools \
  deeptools \
  macs2 \
  ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph \
  r-base=4.3 \
  r-tidyverse r-optparse r-data.table \
  bioconductor-deseq2 bioconductor-edger bioconductor-enhancedvolcano


conda run -n mamba mamba install -n normdb_v2 -c conda-forge -y pandas numpy
conda run -n mamba mamba install -n normdb_v2 -c bioconda -y bioconductor-apeglm


conda activate normdb_v2

# Quick test
command -v bigWigToBedGraph bedGraphToBigWig bedtools diffreps macs2 bamCoverage
Rscript -e 'library(tidyverse); library(DESeq2); library(edgeR); library(EnhancedVolcano); cat("OK\n")'
```

--> All good!

**Export system/package version information for reproducibility:**
- Export an exact “lock” of what is installed with `conda list --explicit > meta/normdb_v2-conda-explicit.txt`
Then people can re-install it with `conda create -n normdb --file normdb_v2-conda-explicit.txt`
- Export a GitHub-friendly with `conda env export --no-builds > meta/normdb_v2-environment.yml`



# Generation of norm bigwig files

## Test1

Let's try to make a command that take uniquely aligned reads (BAM) as input and that output normalized bigwigs, following the local maxima method; here is how the command should look; two modes, paired and single end mode:

**Suggested CLI**: 

```bash

# Paired end mode
normdb normalize \
  --meta samples.tsv \
  --outdir output/normdb_norm \
  --blacklist hg38-blacklist.v2.bed \
  --chrom-sizes GRCh38_chrom_sizes.tab \
  --threads 7 \
  --mode PE \
  --reference auto


# Single end mode
normdb normalize \
  --meta samples.tsv \
  --outdir output/normdb_norm \
  --blacklist hg38-blacklist.v2.bed \
  --chrom-sizes GRCh38_chrom_sizes.tab \
  --threads 7 \
  --mode SE \
  --se-fragment-length 200 \
  --reference auto


# samples.tsv
sample_id	bam	condition	target
WT_R1	output/bowtie2/WT_R1.bam	WT	H3K27me3
WT_R2	output/bowtie2/WT_R2.bam	WT	H3K27me3
KO_R1	output/bowtie2/KO_R1.bam	KO	H3K27me3
KO_R2	output/bowtie2/KO_R2.bam	KO	H3K27me3

```

- meta: meta sample information, tell where to find the bam files and also which comparison to perform (ie. genotype). Need to have `.bam` and index `.bam.bai`
- outdir: indicate where to save the output files, folders will be automatically created if they do not exist



```bash
# Create script
nano scripts/normdb_norm_v1.py

# Make it executable
chmod +x scripts/normdb_norm_v1.py

# Test run
scripts/normdb_norm_v1.py normalize \
  --meta samples.tsv \
  --outdir output/normdb_norm \
  --blacklist hg38-blacklist.v2.bed \
  --chrom-sizes GRCh38_chrom_sizes.tab \
  --threads 7 \
  --mode PE \
  --reference auto

#--> Seems to run
```

Let's do a real test with data from `001*/018*` CutRun for H3K27me3/EZH2 WT/KO; copy the bam files into `output/bowtie2`



```bash
# meta file
sample_id	bam	condition	target
ESC_WT_H3K27me3_R1	output/bowtie2/ESC_WT_H3K27me3_R1.unique.dupmark.sorted.bam	WT	H3K27me3
ESC_WT_H3K27me3_R2	output/bowtie2/ESC_WT_H3K27me3_R2.unique.dupmark.sorted.bam	WT	H3K27me3
ESC_WT_H3K27me3_R3	output/bowtie2/ESC_WT_H3K27me3_R3.unique.dupmark.sorted.bam	WT	H3K27me3
ESC_KO_H3K27me3_R1	output/bowtie2/ESC_KO_H3K27me3_R1.unique.dupmark.sorted.bam	KO	H3K27me3
ESC_KO_H3K27me3_R2	output/bowtie2/ESC_KO_H3K27me3_R2.unique.dupmark.sorted.bam	KO	H3K27me3
ESC_KO_H3K27me3_R3	output/bowtie2/ESC_KO_H3K27me3_R3.unique.dupmark.sorted.bam	KO	H3K27me3
ESC_WT_EZH2_R1	output/bowtie2/ESC_WT_EZH2_R1.unique.dupmark.sorted.bam	WT	EZH2
ESC_WT_EZH2_R2	output/bowtie2/ESC_WT_EZH2_R2.unique.dupmark.sorted.bam	WT	EZH2
ESC_WT_EZH2_R3	output/bowtie2/ESC_WT_EZH2_R3.unique.dupmark.sorted.bam	WT	EZH2
ESC_KO_EZH2_R1	output/bowtie2/ESC_KO_EZH2_R1.unique.dupmark.sorted.bam	KO	EZH2
ESC_KO_EZH2_R2	output/bowtie2/ESC_KO_EZH2_R2.unique.dupmark.sorted.bam	KO	EZH2
ESC_KO_EZH2_R3	output/bowtie2/ESC_KO_EZH2_R3.unique.dupmark.sorted.bam	KO	EZH2




# Run code

scripts/normdb_norm_v1.py normalize \
  --meta meta/samples-001018__WTKO_H3K27me3EZH2.tsv \
  --outdir output/normdb_norm-001018__WTKO_H3K27me3EZH2 \
  --blacklist meta/hg38-blacklist.v2.bed \
  --chrom-sizes meta/GRCh38_chrom_sizes.tab \
  --threads 7 \
  --mode PE \
  --reference auto

#--> Too long, batch

conda activate normdb_v2

sbatch scripts/normdb_norm_v1-001018__WTKO_H3K27me3EZH2.sh # 65410635 xxx



```

--> XXXY HERE !! check if work






# Diff. binding


Let's automate the diff. binding part; input needed:
- Region in which perform the diff. binding count (ie. peak, gene promoter region?) in BED format
- Normalized bedGraph files





Output will include:
- Results table (TSV): log2FC, pvalue, padj, ...
- MA plot (log2FC vs mean signal): Can spot signal/noise issue
- Volcano plot: Quick look whether more gain/lost 
- heatmap of sample correlation: QC plot


**Suggested CLI**: 


```bash
normdb diffbind \
  --meta samples.tsv \
  --regions regions.bed \
  --bigwig-dir output/normdb_norm/06_normalized_bigwig \
  --contrast condition:KO:WT \
  --outdir output/normdb_diffbind \
  --alpha 0.05 \
  --lfc 0
```


- sample: Same as the one used to generate the bigwig
- contrast: <column_name>:<test_level>:<reference_level> so: condition:KO:WT; so test is KO vs WT; positive FC = more in KO




## Test1



```bash
# Create script
nano scripts/normdb_diffbind_v1.py

# Make it executable
chmod +x scripts/normdb_diffbind_v1.py

# Test run
scripts/normdb_diffbind_v1.py \
  --meta samples.tsv \
  --regions regions.bed \
  --bigwig-dir output/normdb_norm/06_normalized_bigwig \
  --contrast condition:KO:WT \
  --outdir output/normdb_diffbind \
  --alpha 0.05 \
  --lfc 0.58 \
  --target H3K27me3

#--> Seems to run
```




Let's do a real test with data from `001*/018*` CutRun for H3K27me3/EZH2 WT/KO; copy the macs2 peak files into `output/macs2/` and count on this.

XXXY HE RE!!!



```bash
# meta file




# Run code

scripts/normdb_diffbind_v1.py \
  --meta samples.tsv \
  --regions regions.bed \
  --bigwig-dir output/normdb_norm/06_normalized_bigwig \
  --contrast condition:KO:WT \
  --outdir output/normdb_diffbind \
  --alpha 0.05 \
  --lfc 0.58 \
  --target H3K27me3


```










xxxxxxx



So to count signal in peak, I did this, some random exmaple maybe no:


### H3K27me3
## qvalue 2.3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.sh # 51089151 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.sh # 51089153 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.sh # 51089157 ok
#### KO
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.sh # 51089158 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.sh # 51089161 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.sh # 51089163 ok
#### OEKO





is this kind of, maybe not the exact command  for sample name but that is the idea


computeMatrix scale-regions -S output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bw \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 100 \
    --regionBodyLength 100




Then DESEQ2 comparions is like this:


### H3K27me3 - merge100bp qval2.3 - R DESEQ2


```bash
conda activate deseq2
```


```R
library("tidyverse")
library("DESeq2")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot <- read.delim("output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)



# import SCORE 
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
  
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, score per row, coordinate and row

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R1")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R2")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R3")




# Tidy into a single tibble
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp = SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3)




######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsKO = SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score)) %>%
  filter(!startsWith(peakID, "chrX")) # 112,254 to 97200


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsKO %>%
  mutate(replicate = paste0(genotype, "_", replicate)) %>%  # Create unique column names
  select(-genotype) %>%  # Remove genotype column (since it's now part of replicate)
  pivot_wider(names_from = replicate, values_from = median_score, values_fill = 0)  
  


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
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKO, -peakID), pull(countData_WTvsKO, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKO_raw <- SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsKO %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKO_raw, -sample), pull(colData_WTvsKO_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 100 # I TESTED 100, 250, 500: VERY COMPARABLE! So pick 100
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.58 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot)
# Export result 
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()




# Identify gain lost peak and genes
## 1) Keep only significant peaks, then keep only strong effects (|log2FC| >= threshold)
res_sig <- res_tibble %>%
  filter(padj < 0.05) %>%
  mutate(effect_dir = case_when(
    log2FoldChange >=  0.58 ~ "pos",
    log2FoldChange <= -0.58 ~ "neg",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(effect_dir))   # drop weak effects
## 2) Determine gene-level direction (Gain / Lost / Mix) from the remaining peaks
gene_direction <- res_sig %>%
  group_by(geneSymbol) %>%
  summarise(
    n_pos = sum(effect_dir == "pos"),
    n_neg = sum(effect_dir == "neg"),
    direction = case_when(
      n_pos > 0 & n_neg > 0 ~ "Mix",
      n_pos > 0             ~ "Gain",
      n_neg > 0             ~ "Lost",
      TRUE                  ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(direction))
### Here remove the 'Mix' genes
res_kept_peaks <- res_sig %>%
  inner_join(gene_direction %>% select(geneSymbol, direction), by = "geneSymbol")

# HERE WITHOUT MIX GENES:
upregulated_genes <- res_kept_peaks %>%
  dplyr::filter(direction == "Gain") %>%
  dplyr::select(geneSymbol) %>% unique()
downregulated_genes <- res_kept_peaks %>%
  dplyr::filter(direction == "Lost") %>%
  dplyr::select(geneSymbol) %>% unique()
write.table(upregulated_genes, file = "output/edgeR/upregulatedNoMix_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated_genes, file = "output/edgeR/downregulatedNoMix_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



# HERE KEEPING MIX GENES:
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.58)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.58)


# Save as coordinates for bed deeptools - including Mix peaks

write_tsv( upregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.bed", col_names = FALSE)
write_tsv( downregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.bed", col_names = FALSE)






Please generate a complete code that will automatically do this using this CLI or something comaprable



normdb diffbind \
  --meta samples.tsv \
  --regions regions.bed \
  --bigwig-dir output/normdb_norm/06_normalized_bigwig \
  --contrast condition:KO:WT \
  --outdir output/normdb_diffbind \
  --alpha 0.05 \
  --lfc 0