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

## Test1 - 001018 CutRun

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

sbatch scripts/normdb_norm_v1-001018__WTKO_H3K27me3EZH2.sh # 65410635 ok (125min)


sbatch scripts/normdb_norm_v2-001018__WTKO_H3K27me3EZH2.sh # 65633087 ok



```

--> `scripts/normdb_norm_v1.py normalize` works GREAT!!! Lets make a tiny modification v2 so that it output the SF we apply to each sample
  --> `scripts/normdb_norm_v2.py normalize` ALL GOOD!



## Test2 - 001002 ChIPseq

Let's test ESC WT vs KO for H3K27me3 from `001/002` 3 bio rep WT and 2 bio rep KO



```bash
# meta file
sample_id	bam	condition	target
ESC_WT_H3K27me3_R1	output/bowtie2/ESC_WT_H3K27me3_R1.CHIP.unique.dupmark.sorted.bam	WT	H3K27me3
ESC_WT_H3K27me3_R2	output/bowtie2/ESC_WT_H3K27me3_R2.CHIP.unique.dupmark.sorted.bam	WT	H3K27me3
ESC_WT_H3K27me3_R3	output/bowtie2/ESC_WT_H3K27me3_R3.CHIP.unique.dupmark.sorted.bam	WT	H3K27me3
ESC_KO_H3K27me3_R1	output/bowtie2/ESC_KO_H3K27me3_R1.CHIP.unique.dupmark.sorted.bam	KO	H3K27me3
ESC_KO_H3K27me3_R2	output/bowtie2/ESC_KO_H3K27me3_R2.CHIP.unique.dupmark.sorted.bam	KO	H3K27me3


conda activate normdb_v2

sbatch scripts/normdb_norm_v2-001002__WTKO_H3K27me3.sh # 65800833 xxx



```













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




## Test1 - 001018 CutRun



```bash
# Create script
nano scripts/normdb_diffbind_v1.py

# Make it executable
chmod +x scripts/normdb_diffbind_v1.py
chmod +x scripts/normdb_diffbind_v2.py

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




Let's do a real test with data from `001*/018*` CutRun for H3K27me3/EZH2 WT/KO; copy the macs2 peak files (`output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed`) into `output/macs2/` and count on this.




```bash
# meta file




# Run code
conda activate normdb_v2

## Test run for H3K27me3
scripts/normdb_diffbind_v1.py \
  --meta meta/samples-001018__WTKO_H3K27me3EZH2.tsv \
  --regions output/macs2/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed \
  --bigwig-dir output/normdb_norm-001018__WTKO_H3K27me3EZH2/06_normalized_bigwig \
  --contrast condition:KO:WT \
  --outdir output/normdb_diff-001018__WTKO_H3K27me3 \
  --alpha 0.05 \
  --lfc 0.58 \
  --target H3K27me3

scripts/normdb_diffbind_v2.py \
  --meta meta/samples-001018__WTKO_H3K27me3EZH2.tsv \
  --regions output/macs2/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed \
  --bigwig-dir output/normdb_norm-001018__WTKO_H3K27me3EZH2/06_normalized_bigwig \
  --contrast condition:KO:WT \
  --outdir output/normdb_diff-001018__WTKO_H3K27me3 \
  --alpha 0.05 \
  --lfc 0.58 \
  --target H3K27me3 \
  --min-counts 100

## Test run for EZH2
scripts/normdb_diffbind_v2.py \
  --meta meta/samples-001018__WTKO_H3K27me3EZH2.tsv \
  --regions output/macs2/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge100bp.bed \
  --bigwig-dir output/normdb_norm-001018__WTKO_H3K27me3EZH2/06_normalized_bigwig \
  --contrast condition:KO:WT \
  --outdir output/normdb_diff-001018__WTKO_EZH2 \
  --alpha 0.05 \
  --lfc 0.58 \
  --target EZH2 \
  --min-counts 100





```


--> `scripts/normdb_diffbind_v1.py` works GREAT!!! Some improvement: export gain/lost signif list, correct bug for coordinate; improve readability of PDF, + add option to remove low reads.
  --> `scripts/normdb_diffbind_v2.py` PERFECT!












