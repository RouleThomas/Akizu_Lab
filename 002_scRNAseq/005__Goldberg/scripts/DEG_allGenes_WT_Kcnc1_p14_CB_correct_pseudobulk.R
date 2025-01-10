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



### Transform seurat to pseudobulk seurat
pseudo_WT_Kcnc1_p14_CB_1step.sct <- AggregateExpression(WT_Kcnc1_p14_CB_1step.sct, assays = "RNA", return.seurat = T, group.by = c("condition", "orig.ident", "cluster.annot"))
### each 'cell' is a genotype-condition-celltype pseudobulk profile
tail(Cells(pseudo_WT_Kcnc1_p14_CB_1step.sct))
head(Cells(pseudo_WT_Kcnc1_p14_CB_1step.sct))
### Extract cell names and update the pseudo seurat with $informations
#### Extract cell names
cell_names <- Cells(pseudo_WT_Kcnc1_p14_CB_1step.sct)
#### Split the cell names into parts
#### The split pattern below splits into 3 parts: condition, orig.ident, and cluster.annot
cell_parts <- strsplit(cell_names, "_", fixed = TRUE)
#### Extract condition, orig.ident, and cluster.annot from each cell name
condition <- sapply(cell_parts, function(x) x[1]) # The first part is the condition
orig_ident <- sapply(cell_parts, function(x) paste(x[2:length(x) - 1], collapse = "_")) # Middle parts are orig.ident, which may contain underscores
cluster_annot <- sapply(cell_parts, function(x) x[length(x)]) # The last part is the cluster annotation
#### Add extracted data as metadata to the Seurat object
pseudo_WT_Kcnc1_p14_CB_1step.sct$condition <- condition
pseudo_WT_Kcnc1_p14_CB_1step.sct$orig.ident <- orig_ident
pseudo_WT_Kcnc1_p14_CB_1step.sct$cluster.annot <- cluster_annot

### Prepare and perform pseudobulk DEG
pseudo_WT_Kcnc1_p14_CB_1step.sct$celltype.stim <- paste(pseudo_WT_Kcnc1_p14_CB_1step.sct$cluster.annot, pseudo_WT_Kcnc1_p14_CB_1step.sct$condition, sep = "_")

Idents(pseudo_WT_Kcnc1_p14_CB_1step.sct) <- "celltype.stim"
DefaultAssay(pseudo_WT_Kcnc1_p14_CB_1step.sct) <- "RNA"

### round the counts to be used by deseq2 (require integer)
#### Access the raw counts matrix in the RNA assay
counts_matrix <- GetAssayData(pseudo_WT_Kcnc1_p14_CB_1step.sct, slot = "counts", assay = "RNA")
#### Round the counts to the nearest integer
rounded_counts_matrix <- round(counts_matrix)
#### Replace the counts in the Seurat object with the rounded counts
pseudo_WT_Kcnc1_p14_CB_1step.sct <- SetAssayData(pseudo_WT_Kcnc1_p14_CB_1step.sct, slot = "counts", new.data = rounded_counts_matrix, assay = "RNA")





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


Granular_1_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p14_CB_1step.sct, 
                         ident.1 = "Granular_1_Kcnc1", 
                         ident.2 = "Granular_1_WT",
                         test.use = "DESeq2")


PLI23_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p14_CB_1step.sct, 
                         ident.1 = "PLI23_Kcnc1", 
                         ident.2 = "PLI23_WT",
                         test.use = "DESeq2")


Granule_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Granule_Kcnc1", 
                         ident.2 = "Granule_WT",
                         test.use = "DESeq2")
Astrocyte_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Astrocyte_Kcnc1", 
                         ident.2 = "Astrocyte_WT",
                         test.use = "DESeq2")
Meningeal_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Meningeal_Kcnc1", 
                         ident.2 = "Meningeal_WT",
                         test.use = "DESeq2")
cat("Upregulated:", sum(Purkinje_pseudobulk$p_val_adj < 0.05 & Purkinje_pseudobulk$avg_log2FC > 0), 
    "\nDownregulated:", sum(Purkinje_pseudobulk$p_val_adj < 0.05 & Purkinje_pseudobulk$avg_log2FC < 0), "\n")
#-=-> DESEQ2 to be run on the raw counts of RNA assay!                         
#--> Few to no DEGs





























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



