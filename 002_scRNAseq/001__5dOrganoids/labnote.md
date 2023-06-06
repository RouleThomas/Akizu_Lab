Data from [Faustino Martins et al 2020; Self-Organizing 3D Human Trunk Neuromuscular Organoids](https://www.sciencedirect.com/science/article/pii/S1934590919305259)


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
sbatch scripts/cellranger_count.sh # 890789 # ok


```
- *NOTE: it is important to rename our files to respect cell ranger nomenclature*

--> run succesfully!



## Then play with our data  (clustering, gene expr)
[cellranger reanalyze](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/reanalyze)

cellranger reanalyze takes feature-barcode matrices produced by cellranger count, cellranger multi, or cellranger aggr and reruns the dimensionality reduction, clustering, and gene expression algorithms using tunable parameter settings.


Use .cloupe file and load it in lpupe [browser](https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser)


--> .cloupe is usefull for quick first look; but boring to change parameters...

Let's use Seurat comonly used after processing instead


## Seurat for clustering

[Tutorial](https://satijalab.org/seurat/articles/get_started.html)

### Setup the Seurat object

```bash
conda activate scRNAseq`
```


```R
library(dplyr)
library(Seurat)
library(patchwork)

# Load the 10x dataset
orga5d.data <- Read10X(data.dir = "count/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
orga5d <- CreateSeuratObject(counts = orga5d.data, project = "pbmc3k", min.cells = 3, min.features = 200) # PAPER used  5 / 500
orga5d

# Check percent reads mapped to Mitchondrial M genome
orga5d[["percent.mt"]] <- PercentageFeatureSet(orga5d, pattern = "^MT-")

# Vizualize QC metrics
pdf("output/seurat/QC_VlnPlot.pdf", width=10, height=6)
VlnPlot(orga5d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(orga5d, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(orga5d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("output/seurat/QC_FeatureScatter.pdf", width=10, height=4)
plot1 + plot2
dev.off()

# Set additional treshold after vizualization
orga5d_sub <- subset(orga5d, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 20)

# Normalization
orga5d_sub_nom <- NormalizeData(orga5d_sub, normalization.method = "LogNormalize", scale.factor = 10000)
```
- *NOTE: Read10X data.dir is the `filtered_feature_bc_matrix` folder in `counts/outs`*
- *NOTE: `CreateSeuratObject` we can filter the min nb of cells and features*
- *NOTE: 10% mitochondrial genome per cell is pretty standard (>20% not great); even though in Seurat tutorial remove cell that have >5%!!*
- *NOTE: `LogNormalize` that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.*

filter cells that have: 
- unique feature counts over 2,500 or less than 200
- >5% mitochondrial counts


### Identification of highly variable features (feature selection)


calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 


```R
# select the top 2000 variable genes with vst method
orga5d_sub_nom_variable <- FindVariableFeatures(orga5d_sub_nom, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(orga5d_sub_nom_variable), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(orga5d_sub_nom_variable)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("output/seurat/VariableFeaturePlot.pdf", width=10, height=4)
plot1 + plot2
dev.off()
```

### Scaling the data

Apply a **linear transformation (‘scaling’)** that is a standard pre-processing step prior to dimensional reduction techniques like PCA:
- Shifts the expression of each gene, so that the **mean expression across cells is 0**
- Scales the expression of each gene, so that the **variance across cells is 1** (gives equal weight in downstream analyses, so that highly-expressed genes do not dominate)


```R
all.genes <- rownames(orga5d_sub_nom_variable)
orga5d_sub_nom_variable_allScaled <- ScaleData(orga5d_sub_nom_variable, features = all.genes) 
```
*NOTE: scaling is here perform on ALL features not the 2000 first; thks to the `features = all.genes`*


### Perform linear dimensional reduction (PCA for QC)

The idea here is to understand what explain heterogeneity in our data; so we need to pick, select relevant PC (like PC1 or 2...). # methods:
- by eye
- JackStraw
- Elbow

#### By eye

PCA on the scaled data. By default, **only the previously determined variable features are used as input**, but can be defined using features argument if you wish to choose a different subset.


```R
# Calculate the PCA
orga5d_sub_nom_variable_allScaled_PCA <- RunPCA(orga5d_sub_nom_variable_allScaled, features = VariableFeatures(object = orga5d_sub_nom_variable_allScaled))

# Vizualize the 1st 15 PC
pdf("output/seurat/DimHeatmap.pdf", width=10, height=20)
DimHeatmap(orga5d_sub_nom_variable_allScaled_PCA, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

```
*NOTE: `DimHeatmap()` allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features are ordered according to their PCA scores. Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets.* 

--> After 17 it is weird

#### JackStraw and Elbow

With this method, generate the plot, and identify the drop-off pvalue = this is the treshold for relevant PC

```R
# Compute JackStraw score
orga5d_sub_nom_variable_allScaled_PCA_Jack <- JackStraw(orga5d_sub_nom_variable_allScaled_PCA, num.replicate = 100)
orga5d_sub_nom_variable_allScaled_PCA_JackScore <- ScoreJackStraw(orga5d_sub_nom_variable_allScaled_PCA_Jack, dims = 1:20)

# Generate plot
pdf("output/seurat/JackStraw.pdf", width=10, height=10)
JackStrawPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore, dims = 1:20)
dev.off()

pdf("output/seurat/Elbow.pdf", width=10, height=10)
ElbowPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore)
dev.off()
```


--> JAckstraw dimensionality is 15-16; Elbow 16-18

--> all together let's pick 17 as treshold





### Cluster the cells

partitioning the cellular distance matrix into clusters based on the dimensionality of the dataset (for us; 17)

```R
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17 <- FindNeighbors(orga5d_sub_nom_variable_allScaled_PCA_JackScore, dims = 1:17)
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster <- FindClusters(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17, resolution = 0.2)
```


*NOTE: For `FindClusters(pbmc, resolution = 0.5)`; parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents() function.*


**--> It gave 12 clusters with resolution 0.5 and with 0.2; 6 clusters**



### Run non-linear dimensional reduction (UMAP/tSNE)


non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. --> to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. 

As **input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.**


```R
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP <- RunUMAP(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster, dims = 1:17)

pdf("output/seurat/UMAP.pdf", width=10, height=10)
DimPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, reduction = "umap", label = TRUE)
dev.off()

# save
saveRDS(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, file = "output/seurat/orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.rds")
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP <- readRDS(file = "output/seurat/orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.rds")
```

--> 6 clusters looks good

### Finding differentially expressed features (cluster biomarkers)


find markers that define clusters via differential expression

```R
# find all markers of cluster 1
cluster1.markers <- FindMarkers(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 2 from clusters 1 and 3
cluster2.markers <- FindMarkers(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, ident.1 = 2, ident.2 = c(1, 3), min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.markers <- FindAllMarkers(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

# Tools for vizualization
## Violin
pdf("output/seurat/VlnPlot_TopMarker.pdf", width=10, height=10)
VlnPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, features = c("CDH6", "NEFM", "HIST1H4C", "CNTN5", "PLCG2", "CRABP1"))
dev.off()

## Umap expr
pdf("output/seurat/FeaturePlot_TopMarker.pdf", width=10, height=10)
FeaturePlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, features = c("CDH6", "NEFM", "HIST1H4C", "CNTN5", "PLCG2", "CRABP1"))
dev.off()

pdf("output/seurat/FeaturePlot_PaperMarker.pdf", width=10, height=20)
FeaturePlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, features = c("FOXC1", "SIX1", "TWIST1", "SOX9", "PAX6", "CDH6", "NEUROG2", "ELAVL4"), ncol = 2)
dev.off()

## heatmap top 20 markers expr
top20 = orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) 

pdf("output/seurat/DoHeatmap_TopMarker.pdf", width=10, height=10)
DoHeatmap(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, features = top20$gene) + NoLegend()
dev.off()

# Rename cluster in agreement
new.cluster.ids <- c("Neuroectodermal progenitors_0", "Mesodermal progenitors_1", "Mesodermal progenitors_2", "Mesodermal progenitors_3", "Mesodermal progenitors_4", "Neurons")
names(new.cluster.ids) <- levels(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP)
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP_rename <- RenameIdents(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, new.cluster.ids)

pdf("output/seurat/UMAP_label.pdf", width=10, height=10)
DimPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP_rename, reduction = "umap", label = TRUE, pt.size = 0.7, label.size = 6) + NoLegend()
dev.off()
```

- *NOTE: `min.pct` requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups.*
- *NOTE: we can set `min.pct` and `thresh.test` to 0 to maybe hgave more DEGs but potentially false positive!*
- *NOTE: **Alternative statistical test** can be used ! see [here](https://satijalab.org/seurat/articles/de_vignette.html)*


## Let's repeat the analysis now; try obtain 4 clusters as in the paper...


```R
library(dplyr)
library(Seurat)
library(patchwork)

# Load the 10x dataset
orga5d.data <- Read10X(data.dir = "count/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
orga5d <- CreateSeuratObject(counts = orga5d.data, project = "pbmc3k", min.cells = 5, min.features = 500) # PAPER used  5 / 500
orga5d

# Check percent reads mapped to Mitchondrial M genome
orga5d[["percent.mt"]] <- PercentageFeatureSet(orga5d, pattern = "^MT-")

# Vizualize QC metrics
pdf("output/seurat/QC_VlnPlot_paper.pdf", width=10, height=6)
VlnPlot(orga5d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(orga5d, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(orga5d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("output/seurat/QC_FeatureScatter_paper.pdf", width=10, height=4)
plot1 + plot2
dev.off()

# Set additional treshold after vizualization
orga5d_sub <- subset(orga5d, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

# Normalization
orga5d_sub_nom <- NormalizeData(orga5d_sub, normalization.method = "LogNormalize", scale.factor = 10000)
```
- *NOTE: Read10X data.dir is the `filtered_feature_bc_matrix` folder in `counts/outs`*
- *NOTE: `CreateSeuratObject` we can filter the min nb of cells and features*
- *NOTE: 10% mitochondrial genome per cell is pretty standard (>20% not great); even though in Seurat tutorial remove cell that have >5%!!*
- *NOTE: `LogNormalize` that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.*

--> I tried to not filter cells; no mor than `min.cells = 5, min.features = 500` as in the paper but the resulting clustering is bad;

--> let's try this `nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25`



### Identification of highly variable features (feature selection)


calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 


```R
# select the top 2000 variable genes with vst method
orga5d_sub_nom_variable <- FindVariableFeatures(orga5d_sub_nom, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(orga5d_sub_nom_variable), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(orga5d_sub_nom_variable)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("output/seurat/VariableFeaturePlot_paper.pdf", width=10, height=4)
plot1 + plot2
dev.off()
```

### Scaling the data

Apply a **linear transformation (‘scaling’)** that is a standard pre-processing step prior to dimensional reduction techniques like PCA:
- Shifts the expression of each gene, so that the **mean expression across cells is 0**
- Scales the expression of each gene, so that the **variance across cells is 1** (gives equal weight in downstream analyses, so that highly-expressed genes do not dominate)


```R
all.genes <- rownames(orga5d_sub_nom_variable)
orga5d_sub_nom_variable_allScaled <- ScaleData(orga5d_sub_nom_variable, features = all.genes) 
```
*NOTE: scaling is here perform on ALL features not the 2000 first; thks to the `features = all.genes`*


### Perform linear dimensional reduction (PCA for QC)

The idea here is to understand what explain heterogeneity in our data; so we need to pick, select relevant PC (like PC1 or 2...). # methods:
- by eye
- JackStraw
- Elbow

#### By eye

PCA on the scaled data. By default, **only the previously determined variable features are used as input**, but can be defined using features argument if you wish to choose a different subset.


```R
# Calculate the PCA
orga5d_sub_nom_variable_allScaled_PCA <- RunPCA(orga5d_sub_nom_variable_allScaled, features = VariableFeatures(object = orga5d_sub_nom_variable_allScaled))

# Vizualize the 1st 15 PC
pdf("output/seurat/DimHeatmap_paper.pdf", width=10, height=20)
DimHeatmap(orga5d_sub_nom_variable_allScaled_PCA, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

```
*NOTE: `DimHeatmap()` allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features are ordered according to their PCA scores. Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets.* 





#### JackStraw and Elbow

With this method, generate the plot, and identify the drop-off pvalue = this is the treshold for relevant PC

```R
# Compute JackStraw score
orga5d_sub_nom_variable_allScaled_PCA_Jack <- JackStraw(orga5d_sub_nom_variable_allScaled_PCA, num.replicate = 100)
orga5d_sub_nom_variable_allScaled_PCA_JackScore <- ScoreJackStraw(orga5d_sub_nom_variable_allScaled_PCA_Jack, dims = 1:20)

# Generate plot
pdf("output/seurat/JackStraw_paper.pdf", width=10, height=10)
JackStrawPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore, dims = 1:20)
dev.off()

pdf("output/seurat/Elbow_paper.pdf", width=10, height=10)
ElbowPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore)
dev.off()
```
--> Jackstraw: 13 ; ELbow 7;  so let's pick 10 (in agreement with 'by eyte'); without filtering

--> Jackstraw:  16 ; ELbow 7;  so let's pick 17 (in agreement with 'by eyte'); with the more stringenat parameter 


### Cluster the cells

partitioning the cellular distance matrix into clusters based on the dimensionality of the dataset (for us; 17)

```R
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17 <- FindNeighbors(orga5d_sub_nom_variable_allScaled_PCA_JackScore, dims = 1:17)
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster <- FindClusters(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17, resolution = 0.2)
```


*NOTE: For `FindClusters(pbmc, resolution = 0.5)`; parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents() function.*


**--> It gave 12 clusters with resolution 0.5 and with 0.2; 6 clusters**
### Run non-linear dimensional reduction (UMAP/tSNE)


non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. --> to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. 

As **input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.**


```R
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP <- RunUMAP(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster, dims = 1:17)

pdf("output/seurat/UMAP_paper.pdf", width=10, height=10)
DimPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, reduction = "umap", label = TRUE)
dev.off()

# save
saveRDS(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, file = "output/seurat/orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.rds")
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP <- readRDS(file = "output/seurat/orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.rds")
```

--> 6 clusters looks good

### Finding differentially expressed features (cluster biomarkers)


find markers that define clusters via differential expression

```R
# find all markers of cluster 1
cluster1.markers <- FindMarkers(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 2 from clusters 1 and 3
cluster2.markers <- FindMarkers(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, ident.1 = 2, ident.2 = c(1, 3), min.pct = 0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.markers <- FindAllMarkers(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

# Tools for vizualization
## Violin
pdf("output/seurat/VlnPlot_TopMarker_paper.pdf", width=10, height=10)
VlnPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, features = c("RBFOX1", "CDH6", "NEFM", "VEGFD", "HIST1H4C", "UBE2C", "CNTN5", "GPC6", "PLCG2", "TSPYL2", "CRABP1", "DLL3"))
dev.off()

## Umap expr
pdf("output/seurat/FeaturePlot_TopMarker_paper.pdf", width=10, height=10)
FeaturePlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, features = c("RBFOX1", "CDH6", "NEFM", "VEGFD", "HIST1H4C", "UBE2C", "CNTN5", "GPC6", "PLCG2", "TSPYL2", "CRABP1", "DLL3"))
dev.off()


## heatmap top 20 markers expr
top20 = orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) 

pdf("output/seurat/DoHeatmap_TopMarker_paper.pdf", width=10, height=10)
DoHeatmap(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, features = top20$gene) + NoLegend()
dev.off()

# look marker from paper

pdf("output/seurat/FeaturePlot_PaperMarker_paper.pdf", width=10, height=10)
FeaturePlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, features = c("FOXC1", "SIX1", "TWIST1", "SOX9", "PAX6", "CDH6", "NEUROG2", "ELAVL4"), ncol = 2)
dev.off()



# Rename cluster in agreement
new.cluster.ids <- c("Neuroectodermal progenitors_0", "Mesodermal progenitors_1", "Mesodermal progenitors_2", "Mesodermal progenitors_3", "Mesodermal progenitors_4", "Neurons")
names(new.cluster.ids) <- levels(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP)
orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP_rename <- RenameIdents(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP, new.cluster.ids)

pdf("output/seurat/UMAP_label.pdf", width=10, height=10)
DimPlot(orga5d_sub_nom_variable_allScaled_PCA_JackScore_17_cluster_UMAP_rename, reduction = "umap", label = TRUE, pt.size = 0.7, label.size = 6) + NoLegend()
dev.off()
```

- *NOTE: `min.pct` requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups.*
- *NOTE: we can set `min.pct` and `thresh.test` to 0 to maybe hgave more DEGs but potentially false positive!*
- *NOTE: **Alternative statistical test** can be used ! see [here](https://satijalab.org/seurat/articles/de_vignette.html)*




--> Not sure what is the best... Maybe the 1st analysis was better. I am not able to have 4 clusters as in the paper.


# New tuto with Sanger institute
Let's follow this [tuto](https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html) that seems to provide more detail on Seurat

Preliminary steps, starting with the `filtered_feature_bc_matrix` from Cellranger 10x:

- ambient RNA correction using `soupX`
- doublet detection using `scrublet`

First step is to apply [soupX](https://github.com/constantAmateur/SoupX) for ambiant RNA contamination:


```R
# install.packages('SoupX')
library("SoupX")

# Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('count/outs')

# Assess % of conta
pdf("output/soupX/autoEstCont_5dOrga.pdf", width=10, height=10)
sc = autoEstCont(sc)
dev.off()

# Generate the corrected matrix
out = adjustCounts(sc)

# Save the matrix
save(out, file = "output/soupX/out.RData")
```
**NOTE: to make it work provide the path to the all 10x analyses! Not the filtered**



Second step is to detect doublet using [scrublet](https://github.com/swolock/scrublet)
```bash
# isntallation
conda activate scRNAseq

pip install scrublet
```
Let's try to create a custom python script (`scrublet.py`) to use scrublet: 
```bash
python3 scrublet.py [input_path] [output_path]

python3 scripts/scrublet_doublets.py count/outs/filtered_feature_bc_matrix output/scrublet
```
--> The script will output if each cell is a doublet or not!!!
**NOTE: do not name any file `scrublet.smthg` or it will fail!**


Then pursue the analyses! Now with our corrected matrix: 
https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html



```R
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")


# Load the soupX corrected counts
load("output/soupX/out.RData")

orga5d <- CreateSeuratObject(counts = out, project = "orga5d") # 

# QUALITY CONTROL add mitochondrial and Ribosomal conta
orga5d[["percent.mt"]] <- PercentageFeatureSet(orga5d, pattern = "^MT-")
orga5d[["percent.rb"]] <- PercentageFeatureSet(orga5d, pattern = "^RP[SL]")

# ADD doublet information (scrublet)
doublets <- read.table("output/doublets/scrublet.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
orga5d <- AddMetaData(orga5d,doublets)
orga5d$Doublet_score <- as.numeric(orga5d$Doublet_score) # make score as numeric
head(orga5d[[]])

# Quality check
pdf("output/seurat_sanger/VlnPlot_QC_5dOrga.pdf", width=10, height=6)
VlnPlot(orga5d, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat_sanger/FeatureScatter_QC_RNAcount_mt_5dOrga.pdf", width=5, height=5)
FeatureScatter(orga5d, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf("output/seurat_sanger/FeatureScatter_QC_RNAcount_RNAfeature_5dOrga.pdf", width=5, height=5)
FeatureScatter(orga5d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
pdf("output/seurat_sanger/FeatureScatter_QC_RNAcount_rb_5dOrga.pdf", width=5, height=5)
FeatureScatter(orga5d, feature1 = "nCount_RNA", feature2 = "percent.rb")
dev.off()
pdf("output/seurat_sanger/FeatureScatter_QC_rb_mt_5dOrga.pdf", width=5, height=5)
FeatureScatter(orga5d, feature1 = "percent.rb", feature2 = "percent.mt")
dev.off()
pdf("output/seurat_sanger/FeatureScatter_QC_RNAfeature_doublet_5dOrga.pdf", width=5, height=5)
FeatureScatter(orga5d, feature1 = "nFeature_RNA", feature2 = "Doublet_score")
dev.off()

# After seeing the plot; add QC information in our seurat object
orga5d[['QC']] <- ifelse(orga5d@meta.data$Is_doublet == 'True','Doublet','Pass')
orga5d[['QC']] <- ifelse(orga5d@meta.data$nFeature_RNA < 500 & orga5d@meta.data$QC == 'Pass','Low_nFeature',orga5d@meta.data$QC)
orga5d[['QC']] <- ifelse(orga5d@meta.data$nFeature_RNA < 500 & orga5d@meta.data$QC != 'Pass' & orga5d@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',orga5d@meta.data$QC,sep = ','),orga5d@meta.data$QC)
orga5d[['QC']] <- ifelse(orga5d@meta.data$percent.mt > 15 & orga5d@meta.data$QC == 'Pass','High_MT',orga5d@meta.data$QC)
orga5d[['QC']] <- ifelse(orga5d@meta.data$nFeature_RNA < 500 & orga5d@meta.data$QC != 'Pass' & orga5d@meta.data$QC != 'High_MT',paste('High_MT',orga5d@meta.data$QC,sep = ','),orga5d@meta.data$QC)
table(orga5d[['QC']])

# Quality check after QC filtering
pdf("output/seurat_sanger/VlnPlot_QCPass_5dOrga.pdf", width=10, height=6)
VlnPlot(subset(orga5d, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

# normalize data
orga5d <- NormalizeData(orga5d, normalization.method = "LogNormalize", scale.factor = 10000)

# Discover the 2000 first more variable genes
orga5d <- FindVariableFeatures(orga5d, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(orga5d), 10)
top10 

## plot the 10 first variable
pdf("output/seurat_sanger/VariableFeaturePlot_top10_5dOrga.pdf", width=10, height=10)
plot1 <- VariableFeaturePlot(orga5d)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
dev.off()

# scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(orga5d)
orga5d <- ScaleData(orga5d, features = all.genes)

# PCA
orga5d <- RunPCA(orga5d, features = VariableFeatures(object = orga5d))

# investigate to find optimal nb of dimension

## Vizualize the 1st 20 PC
pdf("output/seurat_sanger/DimHeatmap.pdf", width=10, height=20)
DimHeatmap(orga5d, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

## Compute JackStraw score
orga5d_Jack <- JackStraw(orga5d, num.replicate = 100)
orga5d_JackScore <- ScoreJackStraw(orga5d_Jack, dims = 1:20)

## Generate plot
pdf("output/seurat_sanger/JackStraw.pdf", width=10, height=10)
JackStrawPlot(orga5d_JackScore, dims = 1:20)
dev.off()

pdf("output/seurat_sanger/Elbow.pdf", width=10, height=10)
ElbowPlot(orga5d)
dev.off()

# clustering
orga5d <- FindNeighbors(orga5d, dims = 1:7)
orga5d <- FindClusters(orga5d, resolution = 0.2)

orga5d <- RunUMAP(orga5d, dims = 1:7, verbose = F)
table(orga5d@meta.data$seurat_clusters) # to check cluster size

pdf("output/seurat_sanger/Umap.pdf", width=10, height=10)
DimPlot(orga5d,label.size = 4,repel = T,label = T)
dev.off()

# Check minor cell population to make sure clustering is good
pdf("output/seurat_sanger/FeaturePlot_Papercluster4.pdf", width=10, height=5)
FeaturePlot(orga5d, features = c("NEUROG2", "ELAVL4"))
dev.off()

# Check parameter to see if one clsuter to do technical stuff
pdf("output/seurat_sanger/FeaturePlot_Doublet_score.pdf", width=10, height=5)
FeaturePlot(orga5d, features = "Doublet_score") & theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat_sanger/FeaturePlot_mt.pdf", width=10, height=5)
FeaturePlot(orga5d, features = "percent.mt") & theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat_sanger/FeaturePlot_RNAfeature.pdf", width=10, height=5)
FeaturePlot(orga5d, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))
dev.off()

# Remove cells that did not path QC
orga5d_QCPass <- subset(orga5d, subset = QC == 'Pass')

pdf("output/seurat_sanger/Umap_QCPass.pdf", width=10, height=10)
DimPlot(orga5d_QCPass,label.size = 4,repel = T,label = T)
dev.off()

# Check cell cycles
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

orga5d_QCPass <- CellCycleScoring(orga5d_QCPass, s.features = s.genes, g2m.features = g2m.genes)
table(orga5d_QCPass[[]]$Phase)

# Let's see if cluster define by technical difference (even after QC filtering)

## percent mt
pdf("output/seurat_sanger/FeaturePlot_QCPass_mt.pdf", width=10, height=10)
FeaturePlot(orga5d_QCPass,features = "percent.mt",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat_sanger/VlnPlot_QCPass_mt.pdf", width=10, height=10)
VlnPlot(orga5d_QCPass,features = "percent.mt") & theme(plot.title = element_text(size=10))
dev.off()
## percent rb
pdf("output/seurat_sanger/FeaturePlot_QCPass_rb.pdf", width=10, height=10)
FeaturePlot(orga5d_QCPass,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat_sanger/VlnPlot_QCPass_rb.pdf", width=10, height=10)
VlnPlot(orga5d_QCPass,features = "percent.rb") & theme(plot.title = element_text(size=10))
dev.off()
## RNA
pdf("output/seurat_sanger/VlnPlot_QCPass_RNA.pdf", width=10, height=10)
VlnPlot(orga5d_QCPass,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## cell cycle
pdf("output/seurat_sanger/FeaturePlot_QCPass_cellCycle.pdf", width=10, height=10)
FeaturePlot(orga5d_QCPass,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()  
pdf("output/seurat_sanger/VlnPlot_QCPass_cellCycle.pdf", width=10, height=10)
VlnPlot(orga5d_QCPass,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()

# SCTransform normalization and clustering
orga5d_QCPass <- SCTransform(orga5d_QCPass, method = "glmGamPoi", ncells = 6941,     # HERE PASTE NB OF CELLS AFTER QC!!! NB OF CEL IN THE orga5d_QCPass
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
orga5d_QCPass

# SCTransform PCA 
orga5d_QCPass <- RunPCA(orga5d_QCPass, verbose = F)
orga5d_QCPass <- RunUMAP(orga5d_QCPass, dims = 1:30, verbose = F)
orga5d_QCPass <- FindNeighbors(orga5d_QCPass, k.param = 15, dims = 1:30, verbose = F) # k.param can be changed
orga5d_QCPass <- FindClusters(orga5d_QCPass, verbose = F, resolution = 0.2, algorithm = 4) # algorithm can be changed
table(orga5d_QCPass[[]]$seurat_clusters)

## Plot
pdf("output/seurat_sanger/Umap_QCPass_SCT.pdf", width=10, height=10)
DimPlot(orga5d_QCPass, label = T)
dev.off()

## Check marker genes
pdf("output/seurat_sanger/FeaturePlot_Papercluster1_SCT.pdf", width=10, height=5)
FeaturePlot(orga5d_QCPass, features = c("FOXC1", "SIX1"))
dev.off()
pdf("output/seurat_sanger/FeaturePlot_Papercluster2_SCT.pdf", width=10, height=5)
FeaturePlot(orga5d_QCPass, features = c("TWIST1", "SOX9"))
dev.off()
pdf("output/seurat_sanger/FeaturePlot_Papercluster3_SCT.pdf", width=10, height=5)
FeaturePlot(orga5d_QCPass, features = c("PAX6", "CDH6"))
dev.off()
pdf("output/seurat_sanger/FeaturePlot_Papercluster4_SCT.pdf", width=10, height=5)
FeaturePlot(orga5d_QCPass, features = c("NEUROG2", "ELAVL4"))
dev.off()

# Save the seurat
save(orga5d_QCPass, file = "output/seurat_sanger/orga5d_QCPass.RData")
load("output/seurat_sanger/orga5d_QCPass.RData")

# Differential expression and marker selection
## set back  active assay to “RNA” and re-do norm and scaling (as we filter out many cells)
DefaultAssay(orga5d_QCPass) <- "RNA"
orga5d_QCPass <- NormalizeData(orga5d_QCPass)
orga5d_QCPass <- FindVariableFeatures(orga5d_QCPass, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(orga5d_QCPass)
orga5d_QCPass <- ScaleData(orga5d_QCPass, features = all.genes)

## DEGs
all.markers <- FindAllMarkers(orga5d_QCPass, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers, file = "output/seurat_sanger/all_markers.tsv", sep = "\t", quote = F, row.names = T)


## check markers
dim(all.markers)
table(all.markers$cluster)
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers


write.table( as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)), file = "output/seurat_sanger/top20_markers.tsv", sep = "\t", quote = F, row.names = T)



# Cell type annotation using SingleR
## Get reference datasets
monaco.ref <- celldex::MonacoImmuneData()
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()

## Convert Seurat object into SCE
sce <- as.SingleCellExperiment(DietSeurat(orga5d_QCPass))
sce

## Analyse
monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)

## check results
table(monaco.main$pruned.labels)
table(hpca.main$pruned.labels)
table(dice.main$pruned.labels)
table(monaco.fine$pruned.labels)
table(hpca.fine$pruned.labels)
table(dice.fine$pruned.labels)

# Add annotation to seurat object
orga5d_QCPass@meta.data$hpca.main   <- hpca.main$pruned.labels
orga5d_QCPass@meta.data$hpca.fine   <- hpca.fine$pruned.labels

orga5d_QCPass@meta.data$monaco.main <- monaco.main$pruned.labels
orga5d_QCPass@meta.data$monaco.fine <- monaco.fine$pruned.labels

orga5d_QCPass@meta.data$dice.main   <- dice.main$pruned.labels
orga5d_QCPass@meta.data$dice.fine   <- dice.fine$pruned.labels

# cluster with name
pdf("output/seurat_sanger/UMAP_SCT_annotated_hpca.main.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "hpca.main")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_hpca.fine.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "hpca.fine")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_monaco.main.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "monaco.main")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_monaco.fine.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "monaco.fine")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_dice.main.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "dice.main")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_dice.fine.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "dice.fine")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()
# manual annotation
# Rename cluster in agreement
## Re run this to re-generate SCTransform

orga5d_QCPass <- SCTransform(orga5d_QCPass, method = "glmGamPoi", ncells = 6941,     # HERE PASTE NB OF CELLS AFTER QC!!! NB OF CEL IN THE orga5d_QCPass
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
orga5d_QCPass

# SCTransform PCA 
orga5d_QCPass <- RunPCA(orga5d_QCPass, verbose = F)
orga5d_QCPass <- RunUMAP(orga5d_QCPass, dims = 1:30, verbose = F)
orga5d_QCPass <- FindNeighbors(orga5d_QCPass, k.param = 15, dims = 1:30, verbose = F) # k.param can be changed
orga5d_QCPass <- FindClusters(orga5d_QCPass, verbose = F, resolution = 0.2, algorithm = 4) # algorithm can be changed
table(orga5d_QCPass[[]]$seurat_clusters)


new.cluster.ids <- c("Mesodermal progenitors", "Neuroectodermal progenitors", "Mesodermal progenitors- Sclerotome", "Unkowkn_1", "Unkowkn_2", "Neurons")
names(new.cluster.ids) <- levels(orga5d_QCPass)
orga5d_QCPass <- RenameIdents(orga5d_QCPass, new.cluster.ids)

pdf("output/seurat_sanger/UMAP_label_manual_fromSCTransform.pdf", width=10, height=10)
DimPlot(orga5d_QCPass, reduction = "umap", label = TRUE, pt.size = 0.7, label.size = 6) + NoLegend()
dev.off()
##


# Find what the cluster 4 and 5 could be

## percent mt
pdf("output/seurat_sanger/FeaturePlot_QCPass_mt_SCT.pdf", width=10, height=10)
FeaturePlot(orga5d_QCPass,features = "percent.mt",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat_sanger/VlnPlot_QCPass_mt_SCT.pdf", width=10, height=10)
VlnPlot(orga5d_QCPass,features = "percent.mt") & theme(plot.title = element_text(size=10))
dev.off()
## percent rb
pdf("output/seurat_sanger/FeaturePlot_QCPass_rb_SCT.pdf", width=10, height=10)
FeaturePlot(orga5d_QCPass,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat_sanger/VlnPlot_QCPass_rb_SCT.pdf", width=10, height=10)
VlnPlot(orga5d_QCPass,features = "percent.rb") & theme(plot.title = element_text(size=10))
dev.off()
## RNA
pdf("output/seurat_sanger/VlnPlot_QCPass_RNA_SCT.pdf", width=10, height=10)
VlnPlot(orga5d_QCPass,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## cell cycle
pdf("output/seurat_sanger/FeaturePlot_QCPass_cellCycle_SCT.pdf", width=10, height=10)
FeaturePlot(orga5d_QCPass,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()  
pdf("output/seurat_sanger/VlnPlot_QCPass_cellCycle_SCT.pdf", width=10, height=10)
VlnPlot(orga5d_QCPass,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()

# Some genes proposed by ChatGPT
pdf("output/seurat_sanger/FeaturePlot_SCT_ChatGPT.pdf", width=20, height=20)
FeaturePlot(orga5d_QCPass, features = c("NES", "PAX6", "PAX7", "MYOD1", "NEUROD1", "MAP2", "SYN1", "MKI67", "PCNA", "BIRC5", "CD3D", "CD19", "PTPRC"))
dev.off()

# Save the seurat if need do more plots
save(orga5d_QCPass, file = "output/seurat_sanger/orga5d_QCPass_final.RData")
load("output/seurat_sanger/orga5d_QCPass_final.RData")


```
- **NOTE: the treshold for QC can be fine-tune depending on our data (this was value from tutorial look likes it fits well with our data)**
- **NOTE: Conventional normalization is to scale it to 10,000 (as if all cells have 10k UMIs overall), and log2-transform the obtained values**
- **NOTE: Additional steps here is the SCTransform normalization and clustering; improve signal noise/ratio:  Single SCTransform command replaces NormalizeData, ScaleData, and FindVariableFeatures --> Take as input the QCpass!!!**
- **NOTE: when using SCTransform we can often use more PC !**
- **NOTE: If need arises, we can separate some clusters manualy. There are also clustering methods geared towards indentification of rare cell populations. Let’s try using fewer neighbors in the KNN graph, combined with Leiden algorithm (now default in scanpy) and slightly increased resolution: `srat <- FindNeighbors(srat, dims = 1:30, k.param = 15, verbose = F)` and  `srat <- FindClusters(srat, verbose = F, algorithm = 4, resolution = 0.9)`**
- **NOTE: at FindClusters can test different alorithm! see [here](https://satijalab.org/seurat/reference/findclusters); to use algorithm 4 I installed it with ` conda install -c conda-forge leidenalg` and `conda install -c anaconda numpy` and `conda install -c anaconda pandas`**
- **NOTE: Reccommended to do the differential expression on the RNA assay, and not the SCTransform**
- **NOTE: singleR automatically assign cell types, infor [here](http://bioconductor.org/books/release/SingleRBook/)** 
- NOTE: at the end; looking at genes it does NOT change something if we look at the seurat (orga5d_QCPass) generated from RNA or VSCTransform


--> The finer cell types annotations are you after, the harder they are to get reliably. Good to compare the different databases** 


To make sure R is using the correct Python version; do:
```bash
which python
```
Check which python is using; and also if I install stuff, where it will be installed. Then force R to use the correct version

```R
library(reticulate)
use_python('~/anaconda3/envs/scRNAseq/bin/python')
```

--> Elbow drop at 7/12; Jack at 13 --> let's pick 7

--> I tested different resolution parameter until I can individualize the neurons cluster; seems 7 dim and 0.2 res is good (6 clusters including neurons)

--> checking technical difference after QC filtering shows that cluster 3 may be composed of dead unhealthy cells

--> If we observe crazy changes of cell cycle score or MT, we could correct this with SCTransform. Not done here as not the case

--> For cleaner code, name `srat` the object and not `orga5d_QCPass`

--> The automatic annotation fail in my case, not super informative; but could be , more info [here](https://bioconductor.org/books/3.12/SingleRBook/)!!

--> Hard to decipher what are cluster 4 and 5 ! No clear pattern when looking at QC metrics... 



# 50 days organoids analysis
Here are the 2 files:
- SRR8734991 organoids 1-3 
- SRR10914868	organoids 4

## Download
```bash
sbatch scripts/download_SRR8734991.sh # 894908
sbatch scripts/download_SRR10914868.sh # 894909
```

## Count matrix with cellranger 10x

There is 2 files but cellrange can handle it, different indexes has been used. Let's generate the count matrix with 10x:
(sequencing machine generated several outputs, which is pretty common in high-throughput sequencing methods like the Illumina platforms)


```bash
# Rename files to respect cell ranger nomenclature
mv input_50dOrga/SRR8734991_1.fastq input_50dOrga/SRR8734991_S1_L001_R1_001.fastq
mv input_50dOrga/SRR8734991_2.fastq input_50dOrga/SRR8734991_S1_L001_R2_001.fastq
mv input_50dOrga/SRR8734991_3.fastq input_50dOrga/SRR8734991_S1_L001_I1_001.fastq

mv input_50dOrga/SRR10914868_1.fastq input_50dOrga/SRR10914868_S1_L001_R1_001.fastq
mv input_50dOrga/SRR10914868_2.fastq input_50dOrga/SRR10914868_S1_L001_R2_001.fastq
mv input_50dOrga/SRR10914868_3.fastq input_50dOrga/SRR10914868_S1_L001_I1_001.fastq

# Compress them
gzip input_50dOrga/SRR*

# start counting
conda activate scRNAseq
sbatch scripts/cellranger_count_50dOrga.sh # 
```

















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