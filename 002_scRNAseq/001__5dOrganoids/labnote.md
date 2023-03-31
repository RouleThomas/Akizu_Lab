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
--> My script failed, need troubleshoot (I updated a new version V5 but need to be tested...)

XXX I launch it in interactive and leave, let see!







OR DO ::
```bash
fasterq-dump SRR8734990
```

**fastqc**


XXX

# Cell Ranger pipeline 
## Generate sc feature counts for a single library







