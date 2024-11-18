#! /usr/bin/env Rscript



# install.packages('SoupX')
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



### Data integration
WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V5_numeric.rds")
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V2_label_ReProcess.rds") # THIS IS THE LAST ONE!
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V1_label.rds") # QC_V4 with PLI12 PLI23


# 2 STEP DATA INTEGRATION #################################
WT_Kcnc1_p14_CB_1step.sct <- SCTransform(WT_Kcnc1_p14_CB_1step.sct, method = "glmGamPoi", ncells = 75857, verbose = TRUE, variable.features.n = 3000) 
WT_Kcnc1_p35_CB_1step.sct <- SCTransform(WT_Kcnc1_p35_CB_1step.sct, method = "glmGamPoi", ncells = 91980, verbose = TRUE, variable.features.n = 3000) 
WT_Kcnc1_p180_CB_1step.sct <- SCTransform(WT_Kcnc1_p180_CB_1step.sct, method = "glmGamPoi", ncells = 68161, verbose = TRUE, variable.features.n = 3000) 

srat.list <- list(WT_Kcnc1_p14_CB_1step.sct = WT_Kcnc1_p14_CB_1step.sct, WT_Kcnc1_p35_CB_1step.sct = WT_Kcnc1_p35_CB_1step.sct, WT_Kcnc1_p180_CB_1step.sct = WT_Kcnc1_p180_CB_1step.sct)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

WT_Kcnc1_CB.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_CB_2step.sct <- IntegrateData(anchorset = WT_Kcnc1_CB.anchors, normalization.method = "SCT")

#### UMAP
DefaultAssay(WT_Kcnc1_CB_2step.sct) <- "integrated"

WT_Kcnc1_CB_2step.sct <- RunPCA(WT_Kcnc1_CB_2step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_CB_2step.sct <- RunUMAP(WT_Kcnc1_CB_2step.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_Kcnc1_CB_2step.sct <- FindNeighbors(WT_Kcnc1_CB_2step.sct, reduction = "pca", k.param = 20, dims = 1:50)
WT_Kcnc1_CB_2step.sct <- FindClusters(WT_Kcnc1_CB_2step.sct, resolution = 0.2, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_CB_2step.sct$condition <- factor(WT_Kcnc1_CB_2step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1_CB-2stepIntegrationRegressNotRepeatedregMtRbCou-dim50kparam20res02.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_CB_2step.sct, reduction = "umap", label=TRUE)
dev.off()




# SAVE OUTPUT
saveRDS(WT_Kcnc1_CB_2step.sct, file = "output/seurat/WT_Kcnc1_CB_2step-dim50kparam20res02.rds")






