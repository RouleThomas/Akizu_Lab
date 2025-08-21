#! /usr/bin/env Rscript



library("Seurat")
library("dplyr")
library("Matrix")
library("ggplot2")
set.seed(42)


# import 

 
Yao_Cortex <- readRDS(file = "output/seurat/Yao_Cortex-10X_nuclei_v3_AIBS-30dim.rds") # 
DefaultAssay(Yao_Cortex) <- "RNA"

WT_Kcnc1_p35_CX_1step <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CX_1step-version2dim35kparam15res065.sct_V1_numeric.rds") # 
DefaultAssay(WT_Kcnc1_p35_CX_1step) <- "RNA"





# Find anchors using PCA projection
shared_features <- intersect(
  rownames(Yao_Cortex),
  rownames(WT_Kcnc1_p35_CX_1step)
)





##############################################################
# The other way ##############################################
##############################################################


# Step 1: Find anchors (you may have already done this)
Yao_Cortex <- RunUMAP(
  Yao_Cortex,
  dims = 1:30,
  reduction = "pca",
  return.model = TRUE  
)

shared_features <- intersect(
  rownames(Yao_Cortex),
  rownames(WT_Kcnc1_p35_CX_1step)
)

anchors <- FindTransferAnchors(
  reference = Yao_Cortex,
  query = WT_Kcnc1_p35_CX_1step,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",
  reduction = "pcaproject",
  dims = 1:30,
  features = shared_features,
  k.anchor = 100, # 100 50 250
  k.filter = 500 # 500
)

#--> TOO LONG!!! Let's run in slurm job without downsamnpling!




# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = Yao_Cortex$cell_type,
  dims = 1:30,
  k.weight = 100
)




WT_Kcnc1_p35_CX_1step <- AddMetaData(WT_Kcnc1_p35_CX_1step, metadata = predictions)


# Step 3: Project query onto reference UMAP
WT_Kcnc1_p35_CX_1step <- MapQuery(
  anchorset = anchors,
  reference = Yao_Cortex,
  query = WT_Kcnc1_p35_CX_1step,
  refdata = list(cluster_id = "cell_type"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots - USING HUMAN GASTRULA ANNOTATION
all_clusters <- union(
  unique(Yao_Cortex$cell_type),
  unique(WT_Kcnc1_p35_CX_1step$predicted.id)
)

# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Yao_Cortex
p1 <- DimPlot(
  Yao_Cortex,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Yao_Cortex")

# Panel 2: Projected cortex p35
p2 <- DimPlot(
  WT_Kcnc1_p35_CX_1step,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("WT_Kcnc1_p35_CX")

# Panel 3: Overlay

# Plot reference (Yao_Cortex) in gray
Yao_Cortex$dummy_group <- "Reference"
p_ref <- DimPlot(
  Yao_Cortex,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Plot projected query separately with colors
p_query <- DimPlot(
  WT_Kcnc1_p35_CX_1step,
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
pdf("output/seurat/UMAP_Yao_Cortex-reference_query_overlay-order2-p35_CX.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()


# SAVE .rds file



saveRDS(Yao_Cortex, file = "output/seurat/Yao_Cortex-10X_nuclei_v3_AIBS-30dim-order2-p35_CX.rds") 
saveRDS(WT_Kcnc1_p35_CX_1step, file = "output/seurat/WT_Kcnc1_p35_CX_1step-10X_nuclei_v3_AIBS-30dim-order2-p35_CX.rds") 



