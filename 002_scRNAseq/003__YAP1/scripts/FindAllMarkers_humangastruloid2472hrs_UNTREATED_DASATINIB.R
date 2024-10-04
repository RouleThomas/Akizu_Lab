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


humangastruloid.combined.sct <- readRDS(file = "output/seurat/humangastruloid2472hr.dim30kparam15res03.rds")

new.cluster.ids <- c(
  "Neurogenic_Progenitors",
  "Muscle_Progenitors",
  "Epiblast_ESC",
  "Endoderm",
  "Progenitor_Cell_Undefined",
  "Nascent_Mesoderm",
  "Progenitor_Cell_Mitotic",
  "Cardiac_Progenitors"
)

names(new.cluster.ids) <- levels(humangastruloid.combined.sct)
humangastruloid.combined.sct <- RenameIdents(humangastruloid.combined.sct, new.cluster.ids)
humangastruloid.combined.sct$cluster.annot <- Idents(humangastruloid.combined.sct) # create a new slot in my seurat object


### Find all markers separated UNTREATED and DASATINIB
humangastruloid_UNTREATED.combined.sct <- subset(humangastruloid.combined.sct, treatment == "UNTREATED")
humangastruloid_DASATINIB.combined.sct <- subset(humangastruloid.combined.sct, treatment == "DASATINIB")


all_markers_UNTREATED <- FindAllMarkers(humangastruloid_UNTREATED.combined.sct,
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")
write.table(all_markers_UNTREATED, file = "output/seurat/srat_humangastruloid.combined.sct-dim30kparam30res03_QCV2-all_markers_UNTREATED_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE) #  

all_markers_DASATINIB <- FindAllMarkers(humangastruloid_DASATINIB.combined.sct,
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")
write.table(all_markers_DASATINIB, file = "output/seurat/srat_humangastruloid.combined.sct-dim30kparam30res03_QCV2-all_markers_DASATINIB_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE) # 



