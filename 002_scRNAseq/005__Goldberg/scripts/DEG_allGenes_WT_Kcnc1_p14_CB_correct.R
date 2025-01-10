#! /usr/bin/env Rscript



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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

# import Seurat object
WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V5_numeric.rds") # regMtRbCount with QC_V3; after Naiara review Goldberg_V2.pptx; QCV3dim30kparam50res035 with name V2 and intenreuron corrected
set.seed(42)


# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "RNA"

WT_Kcnc1_p14_CB_1step.sct <- NormalizeData(WT_Kcnc1_p14_CB_1step.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(WT_Kcnc1_p14_CB_1step.sct)
WT_Kcnc1_p14_CB_1step.sct <- ScaleData(WT_Kcnc1_p14_CB_1step.sct, features = all.genes) # zero-centres and scales it


WT_Kcnc1_p14_CB_1step.sct$celltype.stim <- paste(WT_Kcnc1_p14_CB_1step.sct$cluster.annot, WT_Kcnc1_p14_CB_1step.sct$condition,
    sep = "-")
Idents(WT_Kcnc1_p14_CB_1step.sct) <- "celltype.stim"

# Define the list of clusters for comparison
clusters <- c(
  "PLI23",
  "Granular_1",
  "MLI1",
  "Granular_2",
  "Granular_3",
  "Granular_4",
  "MLI2",
  "Endothelial",
  "Granular_5",
  "Astrocyte",
  "OPC",
  "Bergmann_Glia",
  "PLI12",
  "Oligodendrocyte",
  "Mix_Microglia_Meningeal",
  "Endothelial_Mural" ,
  "Purkinje",
  "Golgi",
  "Unipolar_Brush",
  "Choroid_Plexus"
)
# Initialize an empty list to store the results for each cluster
cluster_markers <- list()
# Loop through each cluster and perform FindMarkers
for (cluster in clusters) {
  cat("Processing cluster:", cluster, "\n")
  # Run FindMarkers for each cluster comparing Kcnc1 vs WT condition
  markers <- FindMarkers(WT_Kcnc1_p14_CB_1step.sct, 
                         ident.1 = paste(cluster, "Kcnc1", sep = "-"), 
                         ident.2 = paste(cluster, "WT", sep = "-"), 
                         verbose = TRUE, 
                         test.use = "wilcox",
                         logfc.threshold = -Inf,
                         min.pct = -Inf,
                         min.diff.pct = -Inf,
                         assay = "RNA")
  # Store the result in the list
  cluster_markers[[cluster]] <- markers
  # Save the result as a text file
  output_filename <- paste0("output/seurat/", cluster, "-Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes_correct.txt")
  write.table(markers, file = output_filename, sep = "\t", quote = FALSE, row.names = TRUE)
}



# Light interactive testing


Granular_4 <- FindMarkers(WT_Kcnc1_p14_CB_1step.sct, 
                         ident.1 = "Granular_4-Kcnc1", 
                         ident.2 = "Granular_4-WT", 
                         verbose = TRUE, 
                         test.use = "wilcox",
                         logfc.threshold = 0,
                         min.pct = 0,
                         assay = "RNA")

Purkinje <- FindMarkers(WT_Kcnc1_p14_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT", 
                         verbose = TRUE, 
                         test.use = "poisson",
                         logfc.threshold = 0,
                         min.pct = 0,
                         assay = "RNA")  


Purkinje <- FindMarkers(WT_Kcnc1_p14_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT", 
                         verbose = TRUE, 
                         test.use = "wilcox",
                         logfc.threshold = 0,
                         min.pct = 0,
                         assay = "RNA",
                         slot= "scale.data")                         



Purkinje_poisson_p14 <- FindMarkers(WT_Kcnc1_p14_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT", 
                         verbose = TRUE, 
                         test.use = "poisson",
                         logfc.threshold = 0,
                         min.pct = 0,
                         assay = "RNA", # Specify the RNA assay (default for raw counts)
                         slot = "counts") # Use raw UMI counts


MLI2_poisson_p14 <- FindMarkers(WT_Kcnc1_p14_CB_1step.sct, 
                         ident.1 = "MLI2-Kcnc1", 
                         ident.2 = "MLI2-WT", 
                         verbose = TRUE, 
                         test.use = "poisson",
                         logfc.threshold = 0,
                         min.pct = 0,
                         assay = "RNA", # Specify the RNA assay (default for raw counts)
                         slot = "counts") # Use raw UMI counts


pdf("output/seurat/test_MLI2_poisson_p14.pdf", width = 8, height = 6)
MLI2_poisson_p14$p_val_adj[MLI2_poisson_p14$p_val_adj == 0] <- 1e-300
MLI2_poisson_p14$log10_pval_adj <- -log10(MLI2_poisson_p14$p_val_adj)
ggplot(MLI2_poisson_p14, aes(x = avg_log2FC, y = log10_pval_adj)) +
  geom_point(color = "black", size = 1) +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
  labs(title = "Volcano Plot", x = "log2 Fold Change (avg_diff)", y = "-log10 Adjusted p-value") +
  theme_minimal()
dev.off()     



