# Pipeline from the paper
- Library prep: Chromium Single Cell Gene Expression system, Single Cell 3' Reagent v2/v3 kits (10X Genomics)
- Sequencing: Illumina HiSeq 4000 for library sequencing
- Barcode processing, mapping, UMI counting, dimension reduction: Cell Ranger v3.0.2, aligned to human GRCh38 reference genome, gene annotations/counting with Ensembl version 93
- Further analysis: Seurat 3.0.3, filtered feature-barcode matrices, remove low-expressed genes/cells, calculate top 2000 variable genes (vst method), calculate mitochondrial transcripts percentage
- Data integration: pre-computed anchorsets, regress out cell cycle effects, perform PCA
Dimensional reduction: UMAP for visualization
- Clustering: SNN modularity optimization-based clustering to identify cell groups
- Visualization: ggplot2, rgl for data visualizations
- Demultiplexing: BD Single-cell Multiplexing Kit, classify sample origin based on highest count per cell
- Cell hashing: oligo-tagged antibodies for cell demultiplexing
- Run velocyto.py annotator: for each mapped bam file, use default parameters for 10X Genomics technology, same gtf file for intron-exon annotation
- Process loom objects: velocyto.R v.0.17, use UMAP embeddings for cell-cell distance calculation, create final velocity plots
RNA velocity estimation: performed with default parameters




# Some helpful tutorial 
- [Here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) ressource for the **10X Chromium Cell Ranger** 
- Data for tutorial is [here](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8734990)




# Import and fastqc 
**Data importation from SRA NCBI**
```bash
module load SRA-Toolkit/2.10.5-centos_linux64

# Use custom script to import, compress and perform fastqc 
python ../../Master/scripts/Import_Compress_QC_V5.py -i SRR8734990 -t P -r 5dOrg
```
--> My script failed, need troubleshoot (I updated a new version V6 but need to be tested...). Maybe because that is a scRNAseq data...

So let's do the old-fashion way:
```bash
fasterq-dump SRR8734990 -S
```
--> It fail for disk-space issue

lets try sbatch with the command inside...
```bash
sbatch scripts/download_SRR8734990.sh # 11839151
```
Fail with: 
```
2023-04-04T16:02:05 fasterq-dump.2.10.5 err: cmn_iter.c cmn_read_uint8_array( #160612353 ).VCursorCellDataDirect() -> RC(rcPS,rcCondition,rcWaiting,rcTimeout,rcExhausted) 
2023-04-04T16:02:05 fasterq-dump.2.10.5 err: row #160612353 : READ.len(134) != QUALITY.len(0) (D) 
2023-04-04T16:02:05 fasterq-dump.2.10.5 fatal: SIGNAL - Segmentation fault 
fasterq-dump (PID 1026649) quit with error code 1
```
Try increase memory (200G instead of 50G) and use --split-files instead of -S. 
It seems that even though it is written paired end, I only have 1 file... !

**--> Tried with  `--include-technical -S` and it works!!!! I have 3 files!!!**

- SRR8734990_S1_L001_R1_001.fastq file contains the Cell Barcode and UMI. In 10x Genomics 3' single-cell RNA-seq protocols, the Read 1 (R1) is typically used for this purpose, and its length is 26bp, which matches the length of sequences in your file.

- SRR8734990_S1_L001_R2_001.fastq file is the Read 2 (R2), which is typically used for gene expression information in 10x Genomics protocols. The length of sequences in this file is 100bp, which is a common read length for R2 in these protocols.

- SRR8734990_S1_L001_I1_001.fastq file likely contains the index reads (I1), which are typically used to distinguish different samples in a multiplexed sequencing run. The length of sequences in this file is 8bp, which is a common length for index reads.



# Cell Ranger pipeline (for 10x)
## install cell ranger

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_in


```bash
cd /scr1/users/roulet/Akizu_Lab/Master/software

curl -o cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1685510849&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODU1MTA4NDl9fX1dfQ__&Signature=gZ2~OED0f8HmpQ29tBzD9Lqqnl4OXcCuSNFV1mpYorVsqGNadHRY-tZ-BYFpXIF93H~jxwphUdEVtrmXUmgYpgjFspGrVhqMb8z6tiidEPCFWcOcTlLKN6qhnsAT6MqeTdznLXkiXsTP5888o7XLCPxOdGO31eSlqmLUhNFUPbFddsziDfgDuSm3WNt57aA7BHQu129PCBnEpb5r9OMRKyKrry2Uh~E5lyllqKNoHV1MXa8sZll0U9Du1uspJ3Vda19LGVycwObiGsIOxdtlS2gyhUXi-9IZJO3TlkOZqwLq7n-3afpNoseCdYSEBU2f0C3VwfuO6OPWArK7cvCq1g__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar -zxvf cellranger-7.1.0.tar.gz

# add cellranger to our PATH
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software/cellranger-7.1.0
# Restart terminal
which cellranger
```
## Download 10x human reference genome
From [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
```bash
cd meta
curl -O meta/https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz
```


## Generate sc feature counts for a single library
[Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) and [cell ranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count)


generate single cell feature counts for a single library:

```bash
# fasq have to be .fastq.gz; so lets zip them
gzip input/SRR8734990_* # it took >2hrs

# rename so that they respect naming convention
mv input/SRR8734990_1.fastq.gz input/SRR8734990_S1_L001_R1_001.fastq.gz
mv input/SRR8734990_2.fastq.gz input/SRR8734990_S1_L001_R2_001.fastq.gz
mv input/SRR8734990_3.fastq.gz input/SRR8734990_S1_L001_I1_001.fastq.gz

# ex of command
cellranger count --id=count \
                   --transcriptome=meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/ \
                   --sample=SRR8734990

# run into sbatch
sbatch scripts/cellranger_count.sh # 890789


```
- *NOTE: it is important to rename our files to respect cell ranger nomenclature*

--> XXX







## Then play with our data  (clustering, gene expr)
[cellranger reanalyze](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/reanalyze)

cellranger reanalyze takes feature-barcode matrices produced by cellranger count, cellranger multi, or cellranger aggr and reruns the dimensionality reduction, clustering, and gene expression algorithms using tunable parameter settings.


Use .cloupe file and load it in lpupe [browser](https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser)



## Seurat for clustering

[Tutorial](https://satijalab.org/seurat/articles/get_started.html)

[Official seurat website](https://satijalab.org/seurat/); [here tuto mouse brain 10X Genomics](https://satijalab.org/seurat/articles/seurat5_sketch_analysis.html)


XXX Cellrang to seurat objecct : https://satijalab.org/seurat/reference/read10x
































# Big workshop
Great tutorial with plenty of ressources and courses and wrkshop [here](https://hbctraining.github.io/scRNA-seq/).


XXX : https://hbctraining.github.io/scRNA-seq/lessons/02_SC_generation_of_count_matrix.html

## Install prerequired packages

Many fail upon installation, notably tidyverse in R, even though I follow what I did for deseq2 lol... So let's, copy our deseq2 environment that works great and have plenty of R stuff already installed and working

```bash
conda create --name scRNAseq --clone deseq2
conda activate scRNAseq
```

In R, install addititonal packages packages one by one:
```R
# package
install.packages("devtools")
install.packages("Seurat")

# bioconductor package
BiocManager::install("SingleCellExperiment")
BiocManager::install("ensembldb")

# load them
library("tidyverse")
library("Matrix")
library("RCurl")
library("scales")
library("cowplot")
library("devtools")
library("Seurat")
library("AnnotationHub")
library("SingleCellExperiment")
library("ensembldb")
```

--> The conda env is working smoothly! All packages can be loaded!

## Getting started

Protocol for 3’ end sequencing --> droplet-based methods (eg. inDrops, 10X Genomics, and Drop-seq)

### Generating the count matrix from the raw sequencing data
[Tuto](https://hbctraining.github.io/scRNA-seq/lessons/02_SC_generation_of_count_matrix.html) 

*NOTE: sometime files are in BCL; to transform into fastq use: `bcl2fastq`* 


If using 10X Genomics library preparation method, then the [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline would be used for all of the above steps:
- Formatting reads and filtering noisy cellular barcodes
- Demultiplexing the samples
- Mapping/pseudo-mapping to transcriptome
- Collapsing UMIs and quantification of reads

Otherwise, done with `umis` (if 1 sample sequenced) or `zUMIs` (if several samples sequenced)

### Auqlity control step
https://hbctraining.github.io/scRNA-seq/lessons/03_SC_quality_control-setup.html



XXX