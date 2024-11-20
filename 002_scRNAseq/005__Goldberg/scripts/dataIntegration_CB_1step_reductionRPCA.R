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



set.seed(42)

# P14
## LOAD QCV3 samples
load_seurat_objects <- function(file_paths) {
  seurat_objects <- list()
  for (file_path in file_paths) {
    sample_name <- gsub("_V3_numeric.rds", "", basename(file_path))
    seurat_objects[[sample_name]] <- readRDS(file_path)
  }
  return(seurat_objects)
}
output_dir <- "output/seurat/"
file_paths <- list.files(output_dir, pattern = "_V3_numeric.rds$", full.names = TRUE)
# Call the function to load the Seurat objects
seurat_objects <- load_seurat_objects(file_paths)
# Loop through the list and assign each Seurat object to a variable with the same name
for (sample_name in names(seurat_objects)) {
  assign(sample_name, seurat_objects[[sample_name]])
}

# SCT Transform
WT_p14_CB_Rep1$replicate <- "Rep1"
WT_p14_CB_Rep2$replicate <- "Rep2"
WT_p14_CB_Rep3$replicate <- "Rep3"

WT_p14_CB_Rep1$condition <- "WT"
WT_p14_CB_Rep2$condition <- "WT"
WT_p14_CB_Rep3$condition <- "WT"

Kcnc1_p14_CB_Rep1$replicate <- "Rep1"
Kcnc1_p14_CB_Rep2$replicate <- "Rep2"
Kcnc1_p14_CB_Rep3$replicate <- "Rep3"

Kcnc1_p14_CB_Rep1$condition <- "Kcnc1"
Kcnc1_p14_CB_Rep2$condition <- "Kcnc1"
Kcnc1_p14_CB_Rep3$condition <- "Kcnc1"



WT_p14_CB_Rep1_SCT <- SCTransform(WT_p14_CB_Rep1, method = "glmGamPoi", ncells = 12369, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p14_CB_Rep2_SCT <- SCTransform(WT_p14_CB_Rep2, method = "glmGamPoi", ncells = 13414, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p14_CB_Rep3_SCT <- SCTransform(WT_p14_CB_Rep3, method = "glmGamPoi", ncells = 13181, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p14_CB_Rep1_SCT <- SCTransform(Kcnc1_p14_CB_Rep1, method = "glmGamPoi", ncells = 10382, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p14_CB_Rep2_SCT <- SCTransform(Kcnc1_p14_CB_Rep2, method = "glmGamPoi", ncells = 10934, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p14_CB_Rep3_SCT <- SCTransform(Kcnc1_p14_CB_Rep3, method = "glmGamPoi", ncells = 15577, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 






# P35

load_seurat_objects <- function(file_paths) {
  seurat_objects <- list()
  for (file_path in file_paths) {
    sample_name <- gsub("_V2_ReProcess_numeric.rds", "", basename(file_path))
    seurat_objects[[sample_name]] <- readRDS(file_path)
  }
  return(seurat_objects)
}
output_dir <- "output/seurat/"
file_paths <- list.files(output_dir, pattern = "_V2_ReProcess_numeric.rds$", full.names = TRUE)
# Call the function to load the Seurat objects
seurat_objects <- load_seurat_objects(file_paths)
# Loop through the list and assign each Seurat object to a variable with the same name
for (sample_name in names(seurat_objects)) {
  assign(sample_name, seurat_objects[[sample_name]])
}




WT_p35_CB_Rep1$replicate <- "Rep1"
WT_p35_CB_Rep2$replicate <- "Rep2"
WT_p35_CB_Rep3$replicate <- "Rep3"

WT_p35_CB_Rep1$condition <- "WT"
WT_p35_CB_Rep2$condition <- "WT"
WT_p35_CB_Rep3$condition <- "WT"

Kcnc1_p35_CB_Rep1$replicate <- "Rep1"
Kcnc1_p35_CB_Rep2$replicate <- "Rep2"  
Kcnc1_p35_CB_Rep3$replicate <- "Rep3"

Kcnc1_p35_CB_Rep1$condition <- "Kcnc1"
Kcnc1_p35_CB_Rep2$condition <- "Kcnc1" 
Kcnc1_p35_CB_Rep3$condition <- "Kcnc1"


WT_p35_CB_Rep1_SCT <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 7299, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2_SCT <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 10683, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3_SCT <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 13664, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1_SCT <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 10264, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2_SCT <- SCTransform(Kcnc1_p35_CB_Rep2, method = "glmGamPoi", ncells = 33623, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3_SCT <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 16447, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 




# P180

load_seurat_objects <- function(file_paths) {
  seurat_objects <- list()
  for (file_path in file_paths) {
    sample_name <- gsub("_V4_numeric.rds", "", basename(file_path))
    seurat_objects[[sample_name]] <- readRDS(file_path)
  }
  return(seurat_objects)
}
output_dir <- "output/seurat/"
file_paths <- list.files(output_dir, pattern = "_V4_numeric.rds$", full.names = TRUE)
# Call the function to load the Seurat objects
seurat_objects <- load_seurat_objects(file_paths)
# Loop through the list and assign each Seurat object to a variable with the same name
for (sample_name in names(seurat_objects)) {
  assign(sample_name, seurat_objects[[sample_name]])
}




WT_p180_CB_Rep1$replicate <- "Rep1"
WT_p180_CB_Rep2$replicate <- "Rep2"
WT_p180_CB_Rep3$replicate <- "Rep3"

WT_p180_CB_Rep1$condition <- "WT"
WT_p180_CB_Rep2$condition <- "WT"
WT_p180_CB_Rep3$condition <- "WT"

Kcnc1_p180_CB_Rep1$replicate <- "Rep1"
Kcnc1_p180_CB_Rep2$replicate <- "Rep2"
Kcnc1_p180_CB_Rep3$replicate <- "Rep3"

Kcnc1_p180_CB_Rep1$condition <- "Kcnc1"
Kcnc1_p180_CB_Rep2$condition <- "Kcnc1"
Kcnc1_p180_CB_Rep3$condition <- "Kcnc1"


WT_p180_CB_Rep1_SCT <- SCTransform(WT_p180_CB_Rep1, method = "glmGamPoi", ncells = 10375, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p180_CB_Rep2_SCT <- SCTransform(WT_p180_CB_Rep2, method = "glmGamPoi", ncells = 10468, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p180_CB_Rep3_SCT <- SCTransform(WT_p180_CB_Rep3, method = "glmGamPoi", ncells = 11641, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep1_SCT <- SCTransform(Kcnc1_p180_CB_Rep1, method = "glmGamPoi", ncells = 10864, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep2_SCT <- SCTransform(Kcnc1_p180_CB_Rep2, method = "glmGamPoi", ncells = 11931, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep3_SCT <- SCTransform(Kcnc1_p180_CB_Rep3, method = "glmGamPoi", ncells = 12882, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 






# Integrate all SCT transform samples in 1-step


srat.list <- list(WT_p14_CB_Rep1_SCT = WT_p14_CB_Rep1_SCT, WT_p14_CB_Rep2_SCT = WT_p14_CB_Rep2_SCT, WT_p14_CB_Rep3_SCT = WT_p14_CB_Rep3_SCT, Kcnc1_p14_CB_Rep1_SCT = Kcnc1_p14_CB_Rep1_SCT, Kcnc1_p14_CB_Rep2_SCT = Kcnc1_p14_CB_Rep2_SCT, Kcnc1_p14_CB_Rep3_SCT = Kcnc1_p14_CB_Rep3_SCT, WT_p35_CB_Rep1_SCT = WT_p35_CB_Rep1_SCT, WT_p35_CB_Rep2_SCT = WT_p35_CB_Rep2_SCT, WT_p35_CB_Rep3_SCT = WT_p35_CB_Rep3_SCT, Kcnc1_p35_CB_Rep1_SCT = Kcnc1_p35_CB_Rep1_SCT, Kcnc1_p35_CB_Rep2_SCT = Kcnc1_p35_CB_Rep2_SCT, Kcnc1_p35_CB_Rep3_SCT = Kcnc1_p35_CB_Rep3_SCT, WT_p180_CB_Rep1_SCT = WT_p180_CB_Rep1_SCT, WT_p180_CB_Rep2_SCT = WT_p180_CB_Rep2_SCT, WT_p180_CB_Rep3_SCT = WT_p180_CB_Rep3_SCT, Kcnc1_p180_CB_Rep1_SCT = Kcnc1_p180_CB_Rep1_SCT, Kcnc1_p180_CB_Rep2_SCT = Kcnc1_p180_CB_Rep2_SCT, Kcnc1_p180_CB_Rep3_SCT = Kcnc1_p180_CB_Rep3_SCT)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_Kcnc1_CB_1step.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features, reduction="rpca") # HERE ADDED  `reduction="rpca"`
WT_Kcnc1_CB_1step.sct <- IntegrateData(anchorset = WT_Kcnc1_CB_1step.anchors, normalization.method = "SCT")


#### UMAP
DefaultAssay(WT_Kcnc1_CB_1step.sct) <- "integrated"

WT_Kcnc1_CB_1step.sct <- RunPCA(WT_Kcnc1_CB_1step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_CB_1step.sct <- RunUMAP(WT_Kcnc1_CB_1step.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_Kcnc1_CB_1step.sct <- FindNeighbors(WT_Kcnc1_CB_1step.sct, reduction = "pca", k.param = 20, dims = 1:50)
WT_Kcnc1_CB_1step.sct <- FindClusters(WT_Kcnc1_CB_1step.sct, resolution = 0.2, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_CB_1step.sct$condition <- factor(WT_Kcnc1_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-dim50kparam20res02.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()


# SAVE OUTPUT
saveRDS(WT_Kcnc1_CB_1step.sct, file = "output/seurat/WT_Kcnc1_CB_1step-dim50kparam20res02.rds")





