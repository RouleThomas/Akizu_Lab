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
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V2_label_ReProcess.rds")



# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "RNA"

WT_Kcnc1_p35_CB_1step.sct <- NormalizeData(WT_Kcnc1_p35_CB_1step.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(WT_Kcnc1_p35_CB_1step.sct)
WT_Kcnc1_p35_CB_1step.sct <- ScaleData(WT_Kcnc1_p35_CB_1step.sct, features = all.genes) # zero-centres and scales it




WT_Kcnc1_p35_CB_1step.sct$celltype.stim <- paste(WT_Kcnc1_p35_CB_1step.sct$cluster.annot, WT_Kcnc1_p35_CB_1step.sct$condition,
    sep = "-")
Idents(WT_Kcnc1_p35_CB_1step.sct) <- "celltype.stim"

# Define the list of clusters for comparison
clusters <- c(
  "Granule",
  "MLI1",
  "PLI23" ,
  "MLI2" ,
  "Endothelial_Stalk",
  "Astrocyte" ,
  "PLI12" ,
  "Golgi" ,
  "Unipolar_Brush" ,
  "Bergman_Glia",
  "Endothelial",
  "Choroid_Plexus",
  "Purkinje",
  "Meningeal",
  "Unknown_Neuron_Subpop",
  "OPC" 
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
  output_filename <- paste0("output/seurat/", cluster, "-Kcnc1_response_p35_CB_QCV2dim50kparam20res03_allGenes_correct.txt")
  write.table(markers, file = output_filename, sep = "\t", quote = FALSE, row.names = TRUE)
}


