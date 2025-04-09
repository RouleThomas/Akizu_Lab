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
WT_Kcnc1_p14_CX_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CX_1step-version2dim30kparam50res07.sct_V1_label.rds") # 
set.seed(42)


# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(WT_Kcnc1_p14_CX_1step.sct) <- "RNA"

WT_Kcnc1_p14_CX_1step.sct <- NormalizeData(WT_Kcnc1_p14_CX_1step.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(WT_Kcnc1_p14_CX_1step.sct)
WT_Kcnc1_p14_CX_1step.sct <- ScaleData(WT_Kcnc1_p14_CX_1step.sct, features = all.genes) # zero-centres and scales it


WT_Kcnc1_p14_CX_1step.sct$celltype.stim <- paste(WT_Kcnc1_p14_CX_1step.sct$cluster.annot, WT_Kcnc1_p14_CX_1step.sct$condition,
    sep = "-")
Idents(WT_Kcnc1_p14_CX_1step.sct) <- "celltype.stim"


# Define the list of clusters for comparison
clusters <- c(
  "L2L3_1",
  "L2L3_2",
  "L6__LT",
  "L4L5",
  "cluster5",
  "L2_L3__inh",
  "GABA_Vipr2",
  "L6__IT",
  "L5_1",
  "Microglia",
  "L5_2",
  "GABA_Pvalb",
  "GABA_Vip",
  "GABA_Sst",
  "L6B",
  "Bergman",
  "OPC",
  "cluster18",
  "GABA_SstCalb2",
  "L5__PT",
  "Endothelial",
  "L2__IT",
  "GABA_Lamp5",
  "Oligodendrocyte",
  "cluster25",
  "cluster26",
  "L5_L6__NP",
  "L6__inh",
  "L6_Car3",
  "Meningeal",
  "cluster31",
  "L5__IT",
  "cluster33",
  "cluster34"
)
# Initialize an empty list to store the results for each cluster
cluster_markers <- list()
# Loop through each cluster and perform FindMarkers
for (cluster in clusters) {
  cat("Processing cluster:", cluster, "\n")
  # Run FindMarkers for each cluster comparing Kcnc1 vs WT condition
  markers <- FindMarkers(WT_Kcnc1_p14_CX_1step.sct, 
                         ident.1 = paste(cluster, "Kcnc1", sep = "-"), 
                         ident.2 = paste(cluster, "WT", sep = "-"), 
                         verbose = TRUE, 
                         test.use = "wilcox",
                         logfc.threshold = -Inf,
                         min.pct = -Inf,
                         min.diff.pct = -Inf,
                         assay = "RNA") # 
  # Store the result in the list
  cluster_markers[[cluster]] <- markers
  # Save the result as a text file
  output_filename <- paste0("output/seurat/", cluster, "-Kcnc1_response_p14_CX_version2dim30kparam50res07_allGenes_default.txt")
  write.table(markers, file = output_filename, sep = "\t", quote = FALSE, row.names = TRUE)
}



# Volcano Plot Generation

library(ggplot2)

# List of clusters (used to locate files)
clusters <- c(
  "L2L3_1",
  "L2L3_2",
  "L6__LT",
  "L4L5",
  "cluster5",
  "L2_L3__inh",
  "GABA_Vipr2",
  "L6__IT",
  "L5_1",
  "Microglia",
  "L5_2",
  "GABA_Pvalb",
  "GABA_Vip",
  "GABA_Sst",
  "L6B",
  "Bergman",
  "OPC",
  "cluster18",
  "GABA_SstCalb2",
  "L5__PT",
  "Endothelial",
  "L2__IT",
  "GABA_Lamp5",
  "Oligodendrocyte",
  "cluster25",
  "cluster26",
  "L5_L6__NP",
  "L6__inh",
  "L6_Car3",
  "Meningeal",
  "cluster31",
  "L5__IT",
  "cluster33",
  "cluster34"
)

# Generate volcano plots for each saved file
for (cluster in clusters) {
  cat("Generating volcano plot for cluster:", cluster, "\n")
  
  # Define the input file path
  input_filename <- paste0("output/seurat/", cluster, "-Kcnc1_response_p14_CX_version2dim30kparam50res07_allGenes_default.txt")
  
  # Check if the file exists
  if (file.exists(input_filename)) {
    # Load the DE results
    markers <- read.table(input_filename, header = TRUE, sep = "\t", row.names = 1)
    
    # Adjust p-values to avoid Inf
    markers$p_val_adj[markers$p_val_adj == 0] <- 1e-300
    markers$log10_pval_adj <- -log10(markers$p_val_adj)
    
    # Create the volcano plot
    volcano_plot <- ggplot(markers, aes(x = avg_log2FC, y = log10_pval_adj)) +
      geom_point(color = "black", size = 1) +
    geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(-0.25, 0.25), color = "red", linetype = "dashed", size = 0.5) +
      labs(
        title = paste("Volcano Plot -", cluster),
        x = "log2 Fold Change (avg_log2FC)",
        y = "-log10 Adjusted p-value"
      ) +
      theme_minimal()
    
    # Save the volcano plot as a PDF
    pdf_filename <- paste0("output/seurat/VolcanoPlot_", cluster, "_p14_CX_version2dim30kparam50res07_allGenes_default.pdf")
    pdf(pdf_filename, width = 7, height = 6)
    print(volcano_plot)
    dev.off()
  } else {
    cat("File not found for cluster:", cluster, "\n")
  }
}






# Light test interactive

# List of clusters (used to locate files)
clusters <- c(
 "L2__IT",
  "L2L3_1",
  "L2L3_2",
  "L2_L3__inh",
  "L4L5",
  "L5_1",
  "L5_2",
  "L5__PT",
  "L5__IT",
  "L5_L6__NP",
  "L6__inh",
  "L6_Car3",
  "L6__LT",
  "L6__IT",
  "L6B",

  "GABA_Vipr2",
  "GABA_Pvalb",
  "GABA_Vip",
  "GABA_Sst",
  "GABA_SstCalb2",
  "GABA_Lamp5",

  "Microglia",
  "Bergman",
  "OPC",
  "Oligodendrocyte",
  "Endothelial",
  "Meningeal",

  "cluster5",
  "cluster18",
  "cluster25",
  "cluster26",
  "cluster31",
  "cluster33",
  "cluster34"
)

## Count up and down
cluster_markers <- list()

# Loop through each cluster to load the saved files
for (cluster in clusters) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p14_CX_version2dim30kparam50res07_allGenes_default.txt")
  
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


