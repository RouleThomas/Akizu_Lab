#! /usr/bin/env Rscript


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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes



set.seed(42)



# UNTREATED and DASATINIB together

humangastruloid.combined.sct <- readRDS(file = "output/seurat/humangastruloid2472hr.dim30kparam15res04.rds")




# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(humangastruloid.combined.sct) <- "RNA"

humangastruloid.combined.sct <- NormalizeData(humangastruloid.combined.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(humangastruloid.combined.sct)
humangastruloid.combined.sct <- ScaleData(humangastruloid.combined.sct, features = all.genes) # zero-centres and scales it


Idents(humangastruloid.combined.sct) <- "cluster.annot"

# cell type marker all genes


all_markers_allGenes <- FindAllMarkers(humangastruloid.combined.sct, assay = "RNA",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf)
write.table(all_markers_allGenes, file = "output/seurat/srat_humangastruloid.combined.sct-dim30kparam15res04-all_markers_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)


