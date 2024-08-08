#! /usr/bin/env Rscript


library("Signac")
library("Seurat")
library("tidyverse")
library("EnsDb.Hsapiens.v86") # hg38
library("EnsDb.Mmusculus.v79") # mm10


set.seed(42)

# import files
RNA_WT_Bap1KO.sct <- readRDS(file = "output/seurat/RNA_WT_Bap1KO.sct_V1_numeric.rds")
ATAC_WT_Bap1KO_integrated.method1 <- readRDS(file = "output/Signac/ATAC_WT_Bap1KO_integrated.method1.rds")


DefaultAssay(RNA_WT_Bap1KO.sct) <- "RNA"
DefaultAssay(ATAC_WT_Bap1KO_integrated.method1) <- "peaks"

# Run FindMultiModalNeighbors
RNA_ATAC_WT_Bap1KO <- FindMultiModalNeighbors(reduction.list = list(RNA_WT_Bap1KO.sct, ATAC_WT_Bap1KO_integrated.method1), dims.list = list(1:20, 2:30)) # Indicate nb dims to use (for RNA we used 20; 30 for ATAC)

# Save 
saveRDS(RNA_ATAC_WT_Bap1KO, file = "output/Signac/FindMultiModalNeighbors_RNA_ATAC_WT_Bap1KO_V1.rds")


