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





# Analyzis with the processed data - version1 messy


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
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)
humangastruloid_dim25kparam15res02 <- RunUMAP(
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
  normalization.method = "SCT",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)



# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = humangastruloid_dim25kparam15res02$cluster.annot,
  dims = 1:25,
  k.weight = 100
)
humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)
# Step 3: Project query onto reference UMAP
## Need to re-run this, not sure why...

humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = humangastruloid_dim25kparam15res02,
  query = humanGastrula,
  refdata = list(cluster_id = "cluster.annot"),
  reference.reduction = "pca",
  reduction.model = "umap"
)


# Plot reference in gray
humangastruloid_dim25kparam15res02$dummy_group <- "Reference"
p_ref <- DimPlot(humangastruloid_dim25kparam15res02, reduction = "umap", group.by = "dummy_group", cols = "lightgray", pt.size = 1) +
  NoLegend()

# Plot query separately to extract point positions + colors
p_query <- DimPlot(humanGastrula, reduction = "ref.umap", group.by = "predicted.id", pt.size = 1) +
  NoAxes() + NoLegend()

# Extract layers and combine manually
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(data = query_layer, aes(x = x, y = y), color = query_layer$colour, size = 1) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))


# Save to PDF
pdf("output/seurat/UMAP_humanGastrula-overlay_v1.pdf", width=10, height=6)
print(g_overlay)
dev.off()




# Step 4: Generate UMAP plots
# Panel 1: Human gastrula
p1 <- DimPlot(humangastruloid_dim25kparam15res02, reduction = "umap", group.by = "cluster.annot", label = TRUE, pt.size = 1) +
  ggtitle("Human gastruloid 72hr")

# Panel 2: Projected gastruloid
p2 <- DimPlot(humanGastrula, reduction = "ref.umap", group.by = "cluster_id", label = TRUE, pt.size = 1) +
  ggtitle("Human gastrula")

# Panel 3: True overlay

# Plot reference in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(humanGastrula, reduction = "umap", group.by = "dummy_group", cols = "lightgray", pt.size = 1) +
  NoLegend()

# Plot query separately to extract point positions + colors
p_query <- DimPlot(humanGastrula, reduction = "ref.umap", group.by = "predicted.id", pt.size = 1) +
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
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay_v1.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
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




# Analyzis with the processed data - version2 clean - human gastruloid 72hrs UNTREATED projection


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


DefaultAssay(humangastruloid_dim25kparam15res02) <- "RNA"


# Do scRNAseq projection (reference humanGasutrla and query humangastruloid_dim25kparam15res02) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
humangastruloid_dim25kparam15res02 <- RunUMAP(
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
  reference = humangastruloid_dim25kparam15res02 ,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

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

# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(humangastruloid_dim25kparam15res02$cluster.annot),
  unique(humanGastrula$predicted.id)
)

# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)

# Panel 1: Human gastruloid
p1 <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "umap",
  group.by = "cluster.annot",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 2: Projected gastrula
p2 <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 3: Overlay
# Reference in gray
humangastruloid_dim25kparam15res02$dummy_group <- "Reference"
p_ref <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Query plotted separately with correct colors
p_query <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version2.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  group = "Human gastrula"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  predicted_id = humanGastrula$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(humangastruloid_dim25kparam15res02) <- "RNA"


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
  features = shared_features,
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

# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(humangastruloid_dim25kparam15res02$predicted.id)
)

# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 3: Overlay

# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Plot projected query separately with colors
p_query <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-version2.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()




# Step 4: Generate UMAP plots - USING OUR GASTRULOID 72hr ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(humangastruloid_dim25kparam15res02$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "ref.umap",
  group.by = "cluster.annot",
  label = FALSE,
  pt.size = 1,
  cols = c(
    "CardiacMesoderm" = "#F8766D", # red
    "CardiacProgenitors" = "#AEA200",
    "NascentMesoderm1" = "#00A6FF",
    "Endoderm" = "#00C1A7",
    "Ectoderm" = "#00BD5C",
    "NascentMesoderm2" = "#EF67EB"
  )
) + ggtitle("Human gastruloid 72hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "ref.umap",
  group.by = "cluster.annot",
  pt.size = 1,
  cols = c(
    "CardiacMesoderm" = "#F8766D", # red
    "CardiacProgenitors" = "#AEA200",
    "NascentMesoderm1" = "#00A6FF",
    "Endoderm" = "#00C1A7",
    "Ectoderm" = "#00BD5C",
    "NascentMesoderm2" = "#EF67EB"
  )
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotation-version2.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()


### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(humangastruloid_dim25kparam15res02, file = "output/seurat/humangastruloid_dim25kparam15res02_scRNAseqProjectionversion2.rds")
humangastruloid_dim25kparam15res02 <- readRDS("output/seurat/humangastruloid_dim25kparam15res02_scRNAseqProjectionversion2.rds")
# 
pdf("output/seurat/UMAP_humanGastrula-query-TBXT-version2.pdf", width = 7, height = 7)
FeaturePlot(
  humangastruloid_dim25kparam15res02,
  features = "TBXT",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 1, cols = c("grey", "red")
) + ggtitle("TBXT expression in human gastruloid 72hr")
dev.off()
pdf("output/seurat/UMAP_humanGastrula-query-EOMES-version2.pdf", width = 7, height = 7)
FeaturePlot(
  humangastruloid_dim25kparam15res02,
  features = "EOMES",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 0.5, cols = c("grey", "red")
) + ggtitle("EOMES expression in human gastruloid 72hr")
dev.off()

pdf("output/seurat/UMAP_humanGastrula-query-TFAP2A-version2.pdf", width = 7, height = 7)
FeaturePlot(
  humangastruloid_dim25kparam15res02,
  features = "TFAP2A",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 3, cols = c("grey", "red")
) + ggtitle("TFAP2A expression in human gastruloid 72hr")
dev.off()
pdf("output/seurat/UMAP_humanGastrula-query-KDR-version2.pdf", width = 7, height = 7)
FeaturePlot(
  humangastruloid_dim25kparam15res02,
  features = "KDR",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 5, cols = c("grey", "red")
) + ggtitle("KDR expression in human gastruloid 72hr")
dev.off()



############################################################





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

--> GOOD for both order













# Analyzis with the processed data - version2 clean - human gastruloid 24hrs UNTREATED projection


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
humangastruloid24hr <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid24hr.combined.sct_25dim.rds")



DefaultAssay(humangastruloid24hr) <- "RNA"


# Do scRNAseq projection (reference humanGasutrla and query humangastruloid24hr) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
humangastruloid24hr <- RunUMAP(
  humangastruloid24hr,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(humanGastrula),
  rownames(humangastruloid24hr)
)



anchors <- FindTransferAnchors(
  reference = humangastruloid24hr ,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = humangastruloid24hr$cluster.annot,
  dims = 1:25,
  k.weight = 100
)

humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)


# Step 3: Project query onto reference UMAP
humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = humangastruloid24hr,
  query = humanGastrula,
  refdata = list(cluster_id = "cluster.annot"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(humangastruloid24hr$cluster.annot),
  unique(humanGastrula$predicted.id)
)

# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)

# Panel 1: Human gastruloid
p1 <- DimPlot(
  humangastruloid24hr,
  reduction = "umap",
  group.by = "cluster.annot",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 2: Projected gastrula
p2 <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 3: Overlay
# Reference in gray
humangastruloid24hr$dummy_group <- "Reference"
p_ref <- DimPlot(
  humangastruloid24hr,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Query plotted separately with correct colors
p_query <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version2-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  group = "Human gastrula"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-prediction_score-24hr.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  predicted_id = humanGastrula$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-prediction_score_predicted_id-24hr.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(humangastruloid24hr) <- "RNA"


# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = humangastruloid24hr,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
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

humangastruloid24hr <- AddMetaData(humangastruloid24hr, metadata = predictions)


# Step 3: Project query onto reference UMAP
humangastruloid24hr <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = humangastruloid24hr,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(humangastruloid24hr$predicted.id)
)

# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  humangastruloid24hr,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 3: Overlay

# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Plot projected query separately with colors
p_query <- DimPlot(
  humangastruloid24hr,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-version2-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()


# Step 4: Generate UMAP plots - USING OUR GASTRULOID 72hr ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(humangastruloid24hr$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  humangastruloid24hr,
  reduction = "ref.umap",
  group.by = "cluster.annot",
  label = FALSE,
  pt.size = 1,
  cols = c(
    "Nascent_Mesoderm" = "#F8766D", # red
    "Primitive_Streak" = "#AEA200",
    "Epiblast" = "#00A6FF",
    "Endoderm" = "#00C1A7",
    "Unkown" = "#00BD5C"
  )
) + ggtitle("Human gastruloid 72hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  humangastruloid24hr,
  reduction = "ref.umap",
  group.by = "cluster.annot",
  pt.size = 1,
  cols = c(
    "Nascent_Mesoderm" = "#F8766D", # red
    "Primitive_Streak" = "#AEA200",
    "Epiblast" = "#00A6FF",
    "Endoderm" = "#00C1A7",
    "Unkown" = "#00BD5C"
  )
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotation-version2-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()

xxx

### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(humangastruloid24hr, file = "output/seurat/humangastruloid24hr_scRNAseqProjectionversion2.rds")
humangastruloid24hr <- readRDS("output/seurat/humangastruloid24hr_scRNAseqProjectionversion2.rds")
# 
pdf("output/seurat/UMAP_humanGastrula-query-TBXT-version2-24hr.pdf", width = 7, height = 7)
FeaturePlot(
  humangastruloid24hr,
  features = "TBXT",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 1, cols = c("grey", "red")
) + ggtitle("TBXT expression in human gastruloid 24hr")
dev.off()


############################################################





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humangastruloid24hr$prediction.score.max,
  group = "72hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score-24hr.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humangastruloid24hr$prediction.score.max,
  predicted_id = humangastruloid24hr$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-24hr.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()

```

--> Looks great! 24h sample with less advanced tissue than 72hr gastruloid




# Analyzis with the processed data - version3 MERGE UNTREATED, DASATINIB, XMU - human gastruloid 24hrs projection


Let's try to directly used the processed data from the authors


```bash
conda activate scRNAseqV3
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






# Import our data - human gastruloid 24h merge 
GASTRU_24h_merge <- readRDS(file = "../003__YAP1/output/seurat/GASTRU_24h_merge-dim20kparam30res03.rds")



DefaultAssay(GASTRU_24h_merge) <- "RNA"


# Do scRNAseq projection (reference humanGasutrla and query GASTRU_24h_merge) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
GASTRU_24h_merge <- RunUMAP(
  GASTRU_24h_merge,
  dims = 1:20,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(humanGastrula),
  rownames(GASTRU_24h_merge)
)



anchors <- FindTransferAnchors(
  reference = GASTRU_24h_merge ,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:20,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = GASTRU_24h_merge$seurat_clusters,
  dims = 1:20,
  k.weight = 100
)

humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)


# Step 3: Project query onto reference UMAP
humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = GASTRU_24h_merge,
  query = humanGastrula,
  refdata = list(cluster_id = "seurat_clusters"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(GASTRU_24h_merge$seurat_clusters),
  unique(humanGastrula$predicted.id)
)

# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)

# Panel 1: Human gastruloid
p1 <- DimPlot(
  GASTRU_24h_merge,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 24hr merge")

# Panel 2: Projected gastrula
p2 <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 3: Overlay
# Reference in gray
GASTRU_24h_merge$dummy_group <- "Reference"
p_ref <- DimPlot(
  GASTRU_24h_merge,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Query plotted separately with correct colors
p_query <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version3-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  group = "Human gastrula"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version3-prediction_score-24hr.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  predicted_id = humanGastrula$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version3-prediction_score_predicted_id-24hr.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(GASTRU_24h_merge) <- "RNA"


# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = GASTRU_24h_merge,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
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

GASTRU_24h_merge <- AddMetaData(GASTRU_24h_merge, metadata = predictions)


# Step 3: Project query onto reference UMAP
GASTRU_24h_merge <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = GASTRU_24h_merge,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_24h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 24hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
# Plot projected query separately with colors
p_query <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()
# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-version3-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()


pdf("output/seurat/UMAP_humanGastrula-reference_noLabel-version3-24hr.pdf", width = 10, height = 7)
DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = FALSE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
dev.off()




# Generate UMAP plots - USING OUR GASTRULOID 24hr CLUSTER ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_24h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 24hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotation-version3-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()
## Condition split
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationSplitCondition-version3-24hr.pdf", width = 30, height = 7)
DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  split.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 24hr")
dev.off()



# Generate UMAP plots - USING OUR GASTRULOID 24hr CONDITION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_24h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 24hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "condition",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationCondition-version3-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()
## Condition split
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationConditionSplit-version3-24hr.pdf", width = 30, height = 7)
DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  split.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 24hr")
dev.off()












xxx

### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(GASTRU_24h_merge, file = "output/seurat/GASTRU_24h_merge_scRNAseqProjectionversion2.rds")
GASTRU_24h_merge <- readRDS("output/seurat/GASTRU_24h_merge_scRNAseqProjectionversion2.rds")
# 
pdf("output/seurat/UMAP_humanGastrula-query-TBXT-version2-24hr.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_24h_merge,
  features = "TBXT",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 1, cols = c("grey", "red")
) + ggtitle("TBXT expression in human gastruloid 24hr")
dev.off()


############################################################





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  group = "24hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score-version3-24hr.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  predicted_id = GASTRU_24h_merge$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-version3-24hr.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df3 <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  predicted_id = GASTRU_24h_merge$predicted.id,
  condition = GASTRU_24h_merge$condition,
  condition = factor(GASTRU_24h_merge$condition, levels = c("UNTREATED", "DASATINIB", "XMU"))
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-splitCondition-version4-24hr.pdf", width = 6, height = 5)
ggplot(df3, aes(x = predicted_id, y = prediction_score, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = "Score") +
  scale_fill_manual(values = c("UNTREATED" = "#999999", "DASATINIB" = "#E69F00", "XMU" = "#56B4E9"))  # Custom colors
dev.off()

```


--> Work great, UNTREATED (more epiblast), DASATINIB, XMU (no Emergent_Mesoderm!) projected differently.




# Analyzis with the processed data - version3 MERGE UNTREATED, DASATINIB, XMU - human gastruloid 72hrs projection


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
GASTRU_72h_merge <- readRDS(file = "../003__YAP1/output/seurat/GASTRU_72h_merge-dim20kparam30res03.rds")


DefaultAssay(GASTRU_72h_merge) <- "RNA"


# Do scRNAseq projection (reference humanGasutrla and query GASTRU_72h_merge) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
GASTRU_72h_merge <- RunUMAP(
  GASTRU_72h_merge,
  dims = 1:20,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(humanGastrula),
  rownames(GASTRU_72h_merge)
)



anchors <- FindTransferAnchors(
  reference = GASTRU_72h_merge ,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:20,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = GASTRU_72h_merge$seurat_clusters,
  dims = 1:20,
  k.weight = 100
)

humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)


# Step 3: Project query onto reference UMAP
humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = GASTRU_72h_merge,
  query = humanGastrula,
  refdata = list(cluster_id = "seurat_clusters"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(GASTRU_72h_merge$seurat_clusters),
  unique(humanGastrula$predicted.id)
)
# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)
# Panel 1: Human gastruloid
p1 <- DimPlot(
  GASTRU_72h_merge,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")
# Panel 2: Projected gastrula
p2 <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
# Panel 3: Overlay
# Reference in gray
GASTRU_72h_merge$dummy_group <- "Reference"
p_ref <- DimPlot(
  GASTRU_72h_merge,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
# Query plotted separately with correct colors
p_query <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()
# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version3.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  group = "Human gastrula"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version3-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  predicted_id = humanGastrula$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version3-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(GASTRU_72h_merge) <- "RNA"


# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = GASTRU_72h_merge,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
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

GASTRU_72h_merge <- AddMetaData(GASTRU_72h_merge, metadata = predictions)


# Step 3: Project query onto reference UMAP
GASTRU_72h_merge <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = GASTRU_72h_merge,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge$predicted.id)
)

# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 3: Overlay

# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Plot projected query separately with colors
p_query <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-version3.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()




# Generate UMAP plots - USING OUR GASTRULOID 72hr ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotation-version3.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Generate UMAP plots - USING OUR GASTRULOID 24hr CONDITION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "condition",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationCondition-version3.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()
## Condition split
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationConditionSplit-version3.pdf", width = 30, height = 7)
DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  split.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
dev.off()











### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(GASTRU_72h_merge, file = "output/seurat/GASTRU_72h_merge_scRNAseqProjectionversion2.rds")
GASTRU_72h_merge <- readRDS("output/seurat/GASTRU_72h_merge_scRNAseqProjectionversion2.rds")
# 
pdf("output/seurat/UMAP_humanGastrula-query-TBXT-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "TBXT",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 1, cols = c("grey", "red")
) + ggtitle("TBXT expression in human gastruloid 72hr")
dev.off()
pdf("output/seurat/UMAP_humanGastrula-query-EOMES-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "EOMES",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 0.5, cols = c("grey", "red")
) + ggtitle("EOMES expression in human gastruloid 72hr")
dev.off()

pdf("output/seurat/UMAP_humanGastrula-query-TFAP2A-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "TFAP2A",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 3, cols = c("grey", "red")
) + ggtitle("TFAP2A expression in human gastruloid 72hr")
dev.off()
pdf("output/seurat/UMAP_humanGastrula-query-KDR-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "KDR",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 5, cols = c("grey", "red")
) + ggtitle("KDR expression in human gastruloid 72hr")
dev.off()



############################################################





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = GASTRU_72h_merge$prediction.score.max,
  group = "72hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score-version3.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = GASTRU_72h_merge$prediction.score.max,
  predicted_id = GASTRU_72h_merge$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-version3.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()



df3 <- data.frame(
  prediction_score = GASTRU_72h_merge$prediction.score.max,
  predicted_id = GASTRU_72h_merge$predicted.id,
  condition = GASTRU_72h_merge$condition,
  condition = factor(GASTRU_72h_merge$condition, levels = c("UNTREATED", "DASATINIB", "XMU"))
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-splitCondition-version3.pdf", width = 6, height = 5)
ggplot(df3, aes(x = predicted_id, y = prediction_score, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = "Score") +
  scale_fill_manual(values = c("UNTREATED" = "#999999", "DASATINIB" = "#E69F00", "XMU" = "#56B4E9"))  # Custom colors
dev.off()

```



--> XMU very different, more Advanced_Mesoderm
--> UNTREATED got more Endoderm







# Analyzis with the processed data - version4 MERGE UNTREATED/DASATINIB (pastQC!) and XMU - human gastruloid 24hrs projection


Let's try to directly used the processed data from the authors


```bash
conda activate scRNAseqV3
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






# Import our data - human gastruloid 24h merge 
GASTRU_24h_merge <- readRDS(file = "../../002_scRNAseq/003__YAP1/output/seurat/GASTRU_24h_merge-dim25kparam30res03_QCkeptUNDASA.rds")



DefaultAssay(GASTRU_24h_merge) <- "RNA"


# Do scRNAseq projection (reference humanGasutrla and query GASTRU_24h_merge) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
GASTRU_24h_merge <- RunUMAP(
  GASTRU_24h_merge,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(humanGastrula),
  rownames(GASTRU_24h_merge)
)



anchors <- FindTransferAnchors(
  reference = GASTRU_24h_merge ,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = GASTRU_24h_merge$seurat_clusters,
  dims = 1:25,
  k.weight = 100
)

humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)


# Step 3: Project query onto reference UMAP
humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = GASTRU_24h_merge,
  query = humanGastrula,
  refdata = list(cluster_id = "seurat_clusters"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(GASTRU_24h_merge$seurat_clusters),
  unique(humanGastrula$predicted.id)
)

# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)

# Panel 1: Human gastruloid
p1 <- DimPlot(
  GASTRU_24h_merge,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 24hr merge")

# Panel 2: Projected gastrula
p2 <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 3: Overlay
# Reference in gray
GASTRU_24h_merge$dummy_group <- "Reference"
p_ref <- DimPlot(
  GASTRU_24h_merge,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Query plotted separately with correct colors
p_query <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  group = "Human gastrula"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4-prediction_score-24hr.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  predicted_id = humanGastrula$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4-prediction_score_predicted_id-24hr.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(GASTRU_24h_merge) <- "RNA"


# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = GASTRU_24h_merge,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
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

GASTRU_24h_merge <- AddMetaData(GASTRU_24h_merge, metadata = predictions)


# Step 3: Project query onto reference UMAP
GASTRU_24h_merge <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = GASTRU_24h_merge,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_24h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 24hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
# Plot projected query separately with colors
p_query <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()
# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-version4-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()


pdf("output/seurat/UMAP_humanGastrula-reference_noLabel-version4-24hr.pdf", width = 10, height = 7)
DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = FALSE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
dev.off()






# Generate UMAP plots - USING OUR GASTRULOID 24hr CONDITION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_24h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 24hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  group.by = "condition",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationCondition-version4-24hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()
## Condition split
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationConditionSplit-version4-24hr.pdf", width = 30, height = 7)
DimPlot(
  GASTRU_24h_merge,
  reduction = "ref.umap",
  split.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 24hr")
dev.off()












xxxy below not mod

### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(GASTRU_24h_merge, file = "output/seurat/GASTRU_24h_merge_scRNAseqProjectionversion2.rds")
GASTRU_24h_merge <- readRDS("output/seurat/GASTRU_24h_merge_scRNAseqProjectionversion2.rds")
# 
pdf("output/seurat/UMAP_humanGastrula-query-TBXT-version2-24hr.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_24h_merge,
  features = "TBXT",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 1, cols = c("grey", "red")
) + ggtitle("TBXT expression in human gastruloid 24hr")
dev.off()


############################################################





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  group = "24hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score-version3-24hr.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  predicted_id = GASTRU_24h_merge$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-version3-24hr.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df3 <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  predicted_id = GASTRU_24h_merge$predicted.id,
  condition = GASTRU_24h_merge$condition,
  condition = factor(GASTRU_24h_merge$condition, levels = c("UNTREATED", "DASATINIB", "XMU"))
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-splitCondition-version4-24hr.pdf", width = 6, height = 5)
ggplot(df3, aes(x = predicted_id, y = prediction_score, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = "Score") +
  scale_fill_manual(values = c("UNTREATED" = "#999999", "DASATINIB" = "#E69F00", "XMU" = "#56B4E9"))  # Custom colors
dev.off()

```








# Analyzis with the processed data - version4 MERGE UNTREATED (pastQC!) ONLY - human gastruloid 24hrs projection


UNTREATED condition only!

```bash
conda activate scRNAseqV3
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






# Import our data - human gastruloid 24h merge 
GASTRU_24h_merge <- readRDS(file = "../../002_scRNAseq/003__YAP1/output/seurat/GASTRU_24h_merge-dim25kparam30res03_QCkeptUNDASA.rds")

GASTRU_24h_merge_UNTREATED <- subset(GASTRU_24h_merge, condition == "UNTREATED")


DefaultAssay(GASTRU_24h_merge_UNTREATED) <- "RNA"


# Do scRNAseq projection (reference humanGasutrla and query GASTRU_24h_merge) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
GASTRU_24h_merge_UNTREATED <- RunUMAP(
  GASTRU_24h_merge_UNTREATED,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(humanGastrula),
  rownames(GASTRU_24h_merge_UNTREATED)
)



anchors <- FindTransferAnchors(
  reference = GASTRU_24h_merge_UNTREATED ,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = GASTRU_24h_merge_UNTREATED$seurat_clusters,
  dims = 1:25,
  k.weight = 100
)

humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)


# Step 3: Project query onto reference UMAP
humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = GASTRU_24h_merge_UNTREATED,
  query = humanGastrula,
  refdata = list(cluster_id = "seurat_clusters"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(GASTRU_24h_merge_UNTREATED$seurat_clusters),
  unique(humanGastrula$predicted.id)
)

# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)

# Panel 1: Human gastruloid
p1 <- DimPlot(
  GASTRU_24h_merge_UNTREATED,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 24hr merge")

# Panel 2: Projected gastrula
p2 <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 3: Overlay
# Reference in gray
GASTRU_24h_merge_UNTREATED$dummy_group <- "Reference"
p_ref <- DimPlot(
  GASTRU_24h_merge_UNTREATED,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Query plotted separately with correct colors
p_query <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_24hr_UNTREATED.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  group = "Human gastrula"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_24hr_UNTREATED-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  predicted_id = humanGastrula$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_24hr_UNTREATED-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(GASTRU_24h_merge_UNTREATED) <- "RNA"


# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = GASTRU_24h_merge_UNTREATED,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
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

GASTRU_24h_merge_UNTREATED <- AddMetaData(GASTRU_24h_merge_UNTREATED, metadata = predictions)


# Step 3: Project query onto reference UMAP
GASTRU_24h_merge_UNTREATED <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = GASTRU_24h_merge_UNTREATED,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_24h_merge_UNTREATED$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_24h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 24hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
# Plot projected query separately with colors
p_query <- DimPlot(
  GASTRU_24h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()
# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-version4_24hr_UNTREATED.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()


pdf("output/seurat/UMAP_humanGastrula-reference_noLabel-version4_24hr_UNTREATED.pdf", width = 10, height = 7)
DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = FALSE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
dev.off()






# Generate UMAP plots - USING OUR GASTRULOID 24hr CONDITION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_24h_merge_UNTREATED$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_24h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 24hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_24h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "condition",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationCondition-version4_24hr_UNTREATED.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()
## Condition split
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationConditionSplit-version4_24hr_UNTREATED.pdf", width = 30, height = 7)
DimPlot(
  GASTRU_24h_merge_UNTREATED,
  reduction = "ref.umap",
  split.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 24hr")
dev.off()












xxxy below not mod

### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(GASTRU_24h_merge, file = "output/seurat/GASTRU_24h_merge_scRNAseqProjectionversion2.rds")
GASTRU_24h_merge <- readRDS("output/seurat/GASTRU_24h_merge_scRNAseqProjectionversion2.rds")
# 
pdf("output/seurat/UMAP_humanGastrula-query-TBXT-version2-24hr.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_24h_merge,
  features = "TBXT",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 1, cols = c("grey", "red")
) + ggtitle("TBXT expression in human gastruloid 24hr")
dev.off()


############################################################





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  group = "24hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score-version3-24hr.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  predicted_id = GASTRU_24h_merge$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-version3-24hr.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df3 <- data.frame(
  prediction_score = GASTRU_24h_merge$prediction.score.max,
  predicted_id = GASTRU_24h_merge$predicted.id,
  condition = GASTRU_24h_merge$condition,
  condition = factor(GASTRU_24h_merge$condition, levels = c("UNTREATED", "DASATINIB", "XMU"))
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-splitCondition-version4-24hr.pdf", width = 6, height = 5)
ggplot(df3, aes(x = predicted_id, y = prediction_score, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = "Score") +
  scale_fill_manual(values = c("UNTREATED" = "#999999", "DASATINIB" = "#E69F00", "XMU" = "#56B4E9"))  # Custom colors
dev.off()

```












# Analyzis with the processed data - version4 MERGE UNTREATED/DASATINIB (pastQC!) XMU - human gastruloid 72hrs projection


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
GASTRU_72h_merge <- readRDS(file = "../../002_scRNAseq/003__YAP1/output/seurat/GASTRU_72h_merge-dim25kparam30res03_QCkeptUNDASA.rds")




DefaultAssay(GASTRU_72h_merge) <- "RNA"


# Do scRNAseq projection (reference humanGasutrla and query GASTRU_72h_merge) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
GASTRU_72h_merge <- RunUMAP(
  GASTRU_72h_merge,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(humanGastrula),
  rownames(GASTRU_72h_merge)
)



anchors <- FindTransferAnchors(
  reference = GASTRU_72h_merge ,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = GASTRU_72h_merge$seurat_clusters,
  dims = 1:25,
  k.weight = 100
)

humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)


# Step 3: Project query onto reference UMAP
humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = GASTRU_72h_merge,
  query = humanGastrula,
  refdata = list(cluster_id = "seurat_clusters"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(GASTRU_72h_merge$seurat_clusters),
  unique(humanGastrula$predicted.id)
)
# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)
# Panel 1: Human gastruloid
p1 <- DimPlot(
  GASTRU_72h_merge,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")
# Panel 2: Projected gastrula
p2 <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
# Panel 3: Overlay
# Reference in gray
GASTRU_72h_merge$dummy_group <- "Reference"
p_ref <- DimPlot(
  GASTRU_72h_merge,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
# Query plotted separately with correct colors
p_query <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()
# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_72hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  group = "Human gastrula"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_72hr-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  predicted_id = humanGastrula$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_72hr-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(GASTRU_72h_merge) <- "RNA"


# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = GASTRU_72h_merge,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
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

GASTRU_72h_merge <- AddMetaData(GASTRU_72h_merge, metadata = predictions)


# Step 3: Project query onto reference UMAP
GASTRU_72h_merge <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = GASTRU_72h_merge,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge$predicted.id)
)

# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 3: Overlay

# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Plot projected query separately with colors
p_query <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-version4_72hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()




# Generate UMAP plots - USING OUR GASTRULOID 72hr ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotation-version4_72hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Generate UMAP plots - USING OUR GASTRULOID 24hr CONDITION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  group.by = "condition",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationCondition-version4_72hr.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()
## Condition split
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationConditionSplit-version4_72hr.pdf", width = 30, height = 7)
DimPlot(
  GASTRU_72h_merge,
  reduction = "ref.umap",
  split.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
dev.off()











### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(GASTRU_72h_merge, file = "output/seurat/GASTRU_72h_merge_scRNAseqProjectionversion2.rds")
GASTRU_72h_merge <- readRDS("output/seurat/GASTRU_72h_merge_scRNAseqProjectionversion2.rds")
# 
pdf("output/seurat/UMAP_humanGastrula-query-TBXT-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "TBXT",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 1, cols = c("grey", "red")
) + ggtitle("TBXT expression in human gastruloid 72hr")
dev.off()
pdf("output/seurat/UMAP_humanGastrula-query-EOMES-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "EOMES",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 0.5, cols = c("grey", "red")
) + ggtitle("EOMES expression in human gastruloid 72hr")
dev.off()

pdf("output/seurat/UMAP_humanGastrula-query-TFAP2A-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "TFAP2A",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 3, cols = c("grey", "red")
) + ggtitle("TFAP2A expression in human gastruloid 72hr")
dev.off()
pdf("output/seurat/UMAP_humanGastrula-query-KDR-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "KDR",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 5, cols = c("grey", "red")
) + ggtitle("KDR expression in human gastruloid 72hr")
dev.off()



############################################################





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = GASTRU_72h_merge$prediction.score.max,
  group = "72hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score-version4_72hr.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = GASTRU_72h_merge$prediction.score.max,
  predicted_id = GASTRU_72h_merge$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-version4_72hr.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()



df3 <- data.frame(
  prediction_score = GASTRU_72h_merge$prediction.score.max,
  predicted_id = GASTRU_72h_merge$predicted.id,
  condition = GASTRU_72h_merge$condition,
  condition = factor(GASTRU_72h_merge$condition, levels = c("UNTREATED", "DASATINIB", "XMU"))
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-splitCondition-version4_72hr.pdf", width = 6, height = 5)
ggplot(df3, aes(x = predicted_id, y = prediction_score, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = "Score") +
  scale_fill_manual(values = c("UNTREATED" = "#999999", "DASATINIB" = "#E69F00", "XMU" = "#56B4E9"))  # Custom colors
dev.off()

```













# Analyzis with the processed data - version4 MERGE UNTREATED (pastQC!) ONLY - human gastruloid 72hrs projection


UNTREATED condition only!


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
GASTRU_72h_merge <- readRDS(file = "../../002_scRNAseq/003__YAP1/output/seurat/GASTRU_72h_merge-dim25kparam30res03_QCkeptUNDASA.rds")

GASTRU_72h_merge_UNTREATED <- subset(GASTRU_72h_merge, condition == "UNTREATED")



DefaultAssay(GASTRU_72h_merge_UNTREATED) <- "RNA"


# Do scRNAseq projection (reference humanGasutrla and query GASTRU_72h_merge) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
GASTRU_72h_merge_UNTREATED <- RunUMAP(
  GASTRU_72h_merge_UNTREATED,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(humanGastrula),
  rownames(GASTRU_72h_merge_UNTREATED)
)



anchors <- FindTransferAnchors(
  reference = GASTRU_72h_merge_UNTREATED ,
  query = humanGastrula,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = GASTRU_72h_merge_UNTREATED$seurat_clusters,
  dims = 1:25,
  k.weight = 100
)

humanGastrula <- AddMetaData(humanGastrula, metadata = predictions)


# Step 3: Project query onto reference UMAP
humanGastrula <- MapQuery(
  anchorset = anchors,
  reference = GASTRU_72h_merge_UNTREATED,
  query = humanGastrula,
  refdata = list(cluster_id = "seurat_clusters"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(GASTRU_72h_merge_UNTREATED$seurat_clusters),
  unique(humanGastrula$predicted.id)
)
# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)
# Panel 1: Human gastruloid
p1 <- DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")
# Panel 2: Projected gastrula
p2 <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")
# Panel 3: Overlay
# Reference in gray
GASTRU_72h_merge_UNTREATED$dummy_group <- "Reference"
p_ref <- DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
# Query plotted separately with correct colors
p_query <- DimPlot(
  humanGastrula,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()
# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_72hr_UNTREATED.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  group = "Human gastrula"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_72hr_UNTREATED-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humanGastrula$prediction.score.max,
  predicted_id = humanGastrula$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-order1-version4_72hr_UNTREATED-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(GASTRU_72h_merge_UNTREATED) <- "RNA"


# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = GASTRU_72h_merge_UNTREATED,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
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

GASTRU_72h_merge_UNTREATED <- AddMetaData(GASTRU_72h_merge_UNTREATED, metadata = predictions)


# Step 3: Project query onto reference UMAP
GASTRU_72h_merge_UNTREATED <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = GASTRU_72h_merge_UNTREATED,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge_UNTREATED$predicted.id)
)

# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 3: Overlay

# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Plot projected query separately with colors
p_query <- DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-version4_72hr_UNTREATED.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()




# Generate UMAP plots - USING OUR GASTRULOID 72hr ANNOTATION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge_UNTREATED$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "seurat_clusters",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotation-version4_72hr_UNTREATED.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Generate UMAP plots - USING OUR GASTRULOID 24hr CONDITION
all_clusters <- union(
  unique(humanGastrula$cluster_id),
  unique(GASTRU_72h_merge_UNTREATED$predicted.id)
)
# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "cluster_id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
# Panel 3: Overlay
# Plot reference (humanGastrula) in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(
  humanGastrula,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()
p_query <- DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "ref.umap",
  group.by = "condition",
  pt.size = 1
) + NoAxes() + NoLegend()
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]
g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
# Step 4: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationCondition-version4_72hr_UNTREATED.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()
## Condition split
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-annotationConditionSplit-version4_72hr_UNTREATED.pdf", width = 30, height = 7)
DimPlot(
  GASTRU_72h_merge_UNTREATED,
  reduction = "ref.umap",
  split.by = "condition",
  label = FALSE,
  pt.size = 1
) + ggtitle("Human gastruloid 72hr")
dev.off()











### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(GASTRU_72h_merge, file = "output/seurat/GASTRU_72h_merge_scRNAseqProjectionversion2.rds")
GASTRU_72h_merge <- readRDS("output/seurat/GASTRU_72h_merge_scRNAseqProjectionversion2.rds")
# 
pdf("output/seurat/UMAP_humanGastrula-query-TBXT-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "TBXT",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 1, cols = c("grey", "red")
) + ggtitle("TBXT expression in human gastruloid 72hr")
dev.off()
pdf("output/seurat/UMAP_humanGastrula-query-EOMES-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "EOMES",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 0.5, cols = c("grey", "red")
) + ggtitle("EOMES expression in human gastruloid 72hr")
dev.off()

pdf("output/seurat/UMAP_humanGastrula-query-TFAP2A-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "TFAP2A",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 3, cols = c("grey", "red")
) + ggtitle("TFAP2A expression in human gastruloid 72hr")
dev.off()
pdf("output/seurat/UMAP_humanGastrula-query-KDR-version2.pdf", width = 7, height = 7)
FeaturePlot(
  GASTRU_72h_merge,
  features = "KDR",
  reduction = "ref.umap",
  pt.size = 1, max.cutoff = 5, cols = c("grey", "red")
) + ggtitle("KDR expression in human gastruloid 72hr")
dev.off()



############################################################





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = GASTRU_72h_merge_UNTREATED$prediction.score.max,
  group = "72hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score-version4_72hr_UNTREATED.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = GASTRU_72h_merge_UNTREATED$prediction.score.max,
  predicted_id = GASTRU_72h_merge_UNTREATED$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-version4_72hr_UNTREATED.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()



df3 <- data.frame(
  prediction_score = GASTRU_72h_merge_UNTREATED$prediction.score.max,
  predicted_id = GASTRU_72h_merge_UNTREATED$predicted.id,
  condition = GASTRU_72h_merge_UNTREATED$condition,
  condition = factor(GASTRU_72h_merge_UNTREATED$condition, levels = c("UNTREATED", "DASATINIB", "XMU"))
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id-splitCondition-version4_72hr_UNTREATED.pdf", width = 6, height = 5)
ggplot(df3, aes(x = predicted_id, y = prediction_score, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = "Score") +
  scale_fill_manual(values = c("UNTREATED" = "#999999", "DASATINIB" = "#E69F00", "XMU" = "#56B4E9"))  # Custom colors
dev.off()

```





