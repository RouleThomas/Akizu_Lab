# Project

For Conchi' project, 3D gastruloid paper (24, and 72hrs); response to reviewer, task from email 20250626 ("we are working on the Stem Cell Reports paper-> back to gastruloids"):

- Compare our untreated 72hr gastruloid with human gastruloid from this [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7615353/)
    - Data:  ArrayExpress under accession code: E-MTAB-9388

--> Re-analyze their data; and perform scRNAseq projection: project our data onto their, and see if we have cell types in common (part of the method describe [here](https://www.science.org/doi/10.1126/sciadv.ado1350#sec-4))



# Data download

From [here](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9388), I click on ENA ID and arrived [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB40781)



There is many FASTQ files because it used MART-seq2 dataset, not 10x Genomics. That means:
- Each cell is processed individually in a well.
- One pair of FASTQ files = one single cell (R1 and R2 for paired-end reads).
- Since they sequenced ~1195 cells, you're seeing ~1195 pairs of FASTQ files.


Instead, let's download already processed file from [here](https://github.com/ScialdoneLab/human-gastrula-shiny), there is the raw reads and cell information! --> Data transfer to `input/`





# Analyzis with the processed data


Let's try to directly used the processed data from the authors


```bash
conda activate scRNAseqV2
```

```R
# packages
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
library("ggpubr")

set.seed(42)


# import humanGastrula .rds files

raw_reads <- readRDS("input/raw_reads.rds")
umap_info <- readRDS("input/umap.rds")



# Set cell names as rownames
rownames(raw_reads) <- umap_info$cell_name

# Convert to seurat
humanGastrula <- CreateSeuratObject(counts = t(raw_reads))  # transpose since genes are in columns

# Build UMAP matrix
umap_mat <- as.matrix(umap_info[, c("X", "X1")])
rownames(umap_mat) <- umap_info$cell_name
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")

# Attach UMAP to Seurat
humanGastrula[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAP_", assay = "RNA")

# Make sure rownames are cell names
rownames(umap_info) <- umap_info$cell_name

# Add desired metadata columns to the Seurat object
humanGastrula$cluster_id <- umap_info$cluster_id
humanGastrula$sub_cluster <- umap_info$sub_cluster


head(humanGastrula@meta.data)


# Re-perform clustering

DefaultAssay(humanGastrula) <- "RNA"

## NORMALIZE AND SCALE DATA 
humanGastrula <- NormalizeData(humanGastrula, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(humanGastrula)
humanGastrula <- ScaleData(humanGastrula, features = all.genes) # zero-centres and scales it


## Find variable features
humanGastrula = FindVariableFeatures(humanGastrula, nfeatures = 3000)

humanGastrula <- RunPCA(humanGastrula, verbose = FALSE, npcs = 25)
humanGastrula <- RunUMAP(humanGastrula, reduction = "pca", dims = 1:25, verbose = FALSE)
humanGastrula <- FindNeighbors(humanGastrula, reduction = "pca", k.param = 15, dims = 1:25)
humanGastrula <- FindClusters(humanGastrula, resolution = 0.3, verbose = FALSE, algorithm = 4)



# Plot to confirm it work

pdf("output/seurat/UMAP_humanGastrula.pdf", width=10, height=6)
DimPlot(humanGastrula, reduction = "umap", label=TRUE, group.by = "cluster_id")
dev.off()


# Check some genes

ParaxialMesoderm= c("TBX6", "MESP2", "HES7", "MSGN1", "NKX1.2")
pdf("output/seurat/FeaturePlot_SCT_humanGastrula-dim25-ParaxialMesoderm.pdf", width=10, height=10)
FeaturePlot(humanGastrula, features = ParaxialMesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()






# Import our data

humangastruloid_dim25kparam15res02 <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid.combined.sct_V2-dim25kparam15res02.rds")
humangastruloid_dim25kparam50res07 <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid.combined.sct_V2-dim25kparam50res07.rds")




# Do scRNAseq projection (reference humanGasutrla and query humangastruloid_dim25kparam15res02) - V1 
# Set assay
DefaultAssay(humanGastrula) <- "RNA"
DefaultAssay(humangastruloid_dim25kparam15res02) <- "RNA"


# Find anchors using PCA projection
anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = humangastruloid_dim25kparam15res02,
  normalization.method = "LogNormalize", # use "SCT" if both objects are SCTransformed; use "LogNormalize" otherwise
  reduction = "pcaproject",
  dims = 1:25,
  k.anchor = 100,
  k.filter = 500
)

# Transfer labels (e.g. cluster_id from reference)
predictions <- TransferData(
  anchorset = anchors,
  refdata = humanGastrula$cluster_id,  # or sub_cluster if you prefer
  dims = 1:25,
  k.weight = 100
)

# Add predicted labels to query object
humangastruloid_dim25kparam15res02 <- AddMetaData(humangastruloid_dim25kparam15res02, metadata = predictions)

pdf("output/seurat/UMAP_humanGastrula-Projected_humangastruloid_dim25kparam15res02.pdf", width=10, height=6)
DimPlot(humangastruloid_dim25kparam15res02, group.by = "predicted.id", label = TRUE)
dev.off()

#--> V1 is the opposite as done in the paper: we annotate our cells using the other cell type annotation


# DO OVERLAY with v1

humanGastrula <- RunUMAP(
  humangastruloid_dim25kparam15res02,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(humanGastrula),
  rownames(humangastruloid_dim25kparam15res02)
)

anchors <- FindTransferAnchors(
  reference = humangastruloid_dim25kparam15res02,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

XXXY HERE TO RUN!!!

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = humangastruloid_dim25kparam15res02$cluster.annot,
  dims = 1:25,
  k.weight = 100
)
humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)
# Step 3: Project query onto reference UMAP
humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = humangastruloid_dim25kparam15res02,
  query = humanGastrula,
  refdata = list(cluster_id = "cluster.annot"),
  reference.reduction = "pca",
  reduction.model = "umap"
)




# Save to PDF
pdf("output/seurat/UMAP_humanGastrula-overlay_v1.pdf", width=10, height=6)
print(g_overlay)
dev.off()







# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = humangastruloid_dim25kparam15res02,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = humanGastrula$cluster_id,
  dims = 1:25,
  k.weight = 100
)

humangastruloid_dim25kparam15res02 <- AddMetaData(humangastruloid_dim25kparam15res02, metadata = predictions)


# Step 3: Project query onto reference UMAP
humangastruloid_dim25kparam15res02 <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = humangastruloid_dim25kparam15res02,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
# Panel 1: Human gastrula
p1 <- DimPlot(humanGastrula, reduction = "umap", group.by = "cluster_id", label = TRUE, pt.size = 1) +
  ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(humangastruloid_dim25kparam15res02, reduction = "ref.umap", group.by = "cluster.annot", label = TRUE, pt.size = 1) +
  ggtitle("Human gastruloid 72hr")

# Panel 3: True overlay

# Plot reference in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(humanGastrula, reduction = "umap", group.by = "dummy_group", cols = "lightgray", pt.size = 1) +
  NoLegend()

# Plot query separately to extract point positions + colors
p_query <- DimPlot(humangastruloid_dim25kparam15res02, reduction = "ref.umap", group.by = "predicted.id", pt.size = 1) +
  NoAxes() + NoLegend()

# Extract layers and combine manually
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(data = query_layer, aes(x = x, y = y), color = query_layer$colour, size = 1) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humangastruloid_dim25kparam15res02$prediction.score.max,
  group = "72hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humangastruloid_dim25kparam15res02$prediction.score.max,
  predicted_id = humangastruloid_dim25kparam15res02$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


```

--> Pretty well recapitulated!!! I generated a plot that say the prediction score for each cluster, it say how well each of the cluster is recapitulated in our data (>0.8 can be considered very high).












# Data analyzis

XXX MAYBE NOT NEEDED AS I TRY THE PROCESSED DATA


Follow the same method as we used in Conchi paper, BUT use their QC filetering parameters; we need clustering as close to what they publish.








