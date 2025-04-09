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

WT_Kcnc1_p14_CX_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CX_1step-version2dim30kparam50res07.sct_V1_label.rds") # 
set.seed(42)


DefaultAssay(WT_Kcnc1_p14_CX_1step.sct) <- "RNA" # Recommended 



# Test different Pathway collections and generate enrichment plot for each cell types (C2 = Pathway, C5 = ontology )
## import Pathways
pathways <- msigdbr("Mus musculus", "C2") %>%          # !!!!!! CHANGE HERE PATHWAYS !!!!!!
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
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
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(WT_Kcnc1_p14_CX_1step.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  Kcnc1 <- seurat_extract(WT_Kcnc1_p14_CX_1step.sct,
                           meta1 = "condition", value_meta1 = "Kcnc1",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, Kcnc1),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_CX_p14_C2_version2dim30kparam50res07_", cluster, ".txt")       # !!!!!! CHANGE HERE PATHWAYS !!!!!!
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}
#--> Long ~2hrs




