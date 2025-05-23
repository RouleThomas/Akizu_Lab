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
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V1_label.rds") # QC_V4 with PLI12 PLI23



# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "RNA"

WT_Kcnc1_p180_CB_1step.sct <- NormalizeData(WT_Kcnc1_p180_CB_1step.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(WT_Kcnc1_p180_CB_1step.sct)
WT_Kcnc1_p180_CB_1step.sct <- ScaleData(WT_Kcnc1_p180_CB_1step.sct, features = all.genes) # zero-centres and scales it




WT_Kcnc1_p180_CB_1step.sct$celltype.stim <- paste(WT_Kcnc1_p180_CB_1step.sct$cluster.annot, WT_Kcnc1_p180_CB_1step.sct$condition,
    sep = "-")
Idents(WT_Kcnc1_p180_CB_1step.sct) <- "celltype.stim"

# Define the list of clusters for comparison
clusters <- c(
  "Granular_1",
  "Granular_2",
  "MLI1",
  "Granular_3",
  "MLI2",
  "PLI23",
  "Astrocyte",
  "PLI12",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)
# Initialize an empty list to store the results for each cluster
cluster_markers <- list()
# Loop through each cluster and perform FindMarkers
for (cluster in clusters) {
  cat("Processing cluster:", cluster, "\n")
  # Run FindMarkers for each cluster comparing Kcnc1 vs WT condition
  markers <- FindMarkers(WT_Kcnc1_p180_CB_1step.sct, 
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
  output_filename <- paste0("output/seurat/", cluster, "-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct.txt")
  write.table(markers, file = output_filename, sep = "\t", quote = FALSE, row.names = TRUE)
}





# Light interactive testing
## Count up and down
cluster_markers <- list()

# Loop through each cluster to load the saved files
for (cluster in clusters) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct.txt")
  
  if (file.exists(file_path)) {
    cat("Loading file for cluster:", cluster, "\n")
    # Read the file and store it in the list
    cluster_markers[[cluster]] <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)
  } else {
    cat("File not found for cluster:", cluster, "\n")
  }
}
regulation_summary <- tibble(cluster = character(), up = numeric(), down = numeric())

# Loop through the clusters and calculate up and downregulated genes
for (cluster in names(cluster_markers)) {
  markers <- cluster_markers[[cluster]]
  
  # Count upregulated genes
  upregulated <- markers %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
    nrow()
  
  # Count downregulated genes
  downregulated <- markers %>%
    filter(p_val_adj < 0.05, avg_log2FC < -0.25) %>%
    nrow()
  
  # Add the results to the tibble
  regulation_summary <- regulation_summary %>%
    add_row(cluster = cluster, up = upregulated, down = downregulated)
}

# Print the summary
print(regulation_summary)



