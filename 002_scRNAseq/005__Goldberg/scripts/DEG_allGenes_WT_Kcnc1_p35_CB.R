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
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_label.rds")




# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "RNA"

WT_Kcnc1_p35_CB_1step.sct$celltype.stim <- paste(WT_Kcnc1_p35_CB_1step.sct$cluster.annot, WT_Kcnc1_p35_CB_1step.sct$condition,
    sep = "-")
Idents(WT_Kcnc1_p35_CB_1step.sct) <- "celltype.stim"

# Define the list of clusters for comparison
clusters <- c(
  "Granular",
  "Interneuron",
  "MLI1",
  "MLI2_1",
  "MLI2_2",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergman_Glia",
  "OPC",
  "Meningeal",
  "Endothelial",
  "Choroid_Plexus",
  "Endothelial_Stalk"
)
# Initialize an empty list to store the results for each cluster
cluster_markers <- list()
# Loop through each cluster and perform FindMarkers
for (cluster in clusters) {
  cat("Processing cluster:", cluster, "\n")
  # Run FindMarkers for each cluster comparing Kcnc1 vs WT condition
  markers <- FindMarkers(WT_Kcnc1_p35_CB_1step.sct, 
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
  output_filename <- paste0("output/seurat/", cluster, "-Kcnc1_response_p35_CB_QCV3dim50kparam50res03_allGenes.txt")
  write.table(markers, file = output_filename, sep = "\t", quote = FALSE, row.names = TRUE)
}


