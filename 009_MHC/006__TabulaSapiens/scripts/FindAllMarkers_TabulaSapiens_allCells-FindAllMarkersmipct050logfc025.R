#! /usr/bin/env Rscript



# library

library("Seurat")
library("BPCells")
library("Azimuth")
library("biomaRt")
library("Matrix")
library("ggplot2")
set.seed(42)

# load rds object

TabulaSapiens_allCells <- readRDS(file = "output/seurat/TabulaSapiens_allCells-dim30.rds") # 



Idents(TabulaSapiens_allCells) <- "cell_type"
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(TabulaSapiens_allCells) <- "RNA"
TabulaSapiens_allCells <- NormalizeData(TabulaSapiens_allCells, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(TabulaSapiens_allCells)
TabulaSapiens_allCells <- ScaleData(TabulaSapiens_allCells, features = all.genes) # zero-centres and scales it

all_markers <- FindAllMarkers(TabulaSapiens_allCells, assay = "RNA", only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_TabulaSapiens_allCells-dim30-FindAllMarkersmipct050logfc025-cell_type.txt", sep = "\t", quote = FALSE, row.names = TRUE)




