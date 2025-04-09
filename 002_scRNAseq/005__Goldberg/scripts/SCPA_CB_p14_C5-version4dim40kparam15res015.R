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
library("gprofiler2") 
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")
set.seed(42)



# Data import 


WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step-version4dim40kparam15res015.sct_V1_label.rds") # 
set.seed(42)


DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "RNA" # Recommended 



# Test different Pathway collections and generate enrichment plot for each cell types (C2 = Pathway, C5 = ontology )
## import Pathways
pathways <- msigdbr("Mus musculus", "C5") %>%          # !!!!!! CHANGE HERE PATHWAYS !!!!!!
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
"Granule",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Purkinje",
  "Golgi",
  "UnipolarBrush",
  "Oligodendrocyte",
  "Astrocyte",
  "Bergman",
  "Endothelial",
  "Meningeal",
  "ChoroidPlexus",
  "Unknown"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(WT_Kcnc1_p14_CB_1step.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  Kcnc1 <- seurat_extract(WT_Kcnc1_p14_CB_1step.sct,
                           meta1 = "condition", value_meta1 = "Kcnc1",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, Kcnc1),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_CB_p14_C5_version4dim40kparam15res015_", cluster, ".txt")       # !!!!!! CHANGE HERE PATHWAYS !!!!!!
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}
# 



