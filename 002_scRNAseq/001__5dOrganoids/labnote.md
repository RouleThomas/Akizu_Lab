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
- 




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





**fastqc**


XXX

# Cell Ranger pipeline 
## Generate sc feature counts for a single library






