# Project

Test the BD Rhapsody pipeline for scRNAseq. Compare result with already processed data.

Design:
- Normoxia vs hypoxia at fetal (neonatal) stage
    - 2 Bio rep per condition (1 male, 1 female)


--> Previous marker genes used by Ana can be found at `011*/docs/Annotation.xlsx`

# Data access

Access Cristancho lab folder from CHOP computer: `/mnt/isilon/cristancho_data`

Access data from [Seven Bridge - BD Rhapsody](https://igor.sbgenomics.com/home)
--> Create [Seven Bridge](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/setup/sbg/top_sbg_setup.html) account and ask data file access to Paulo 


**Download seurat object**
```bash
wget -O output/seurat/AC_WTA_SMK_index_library_Seurat.rds \
'https://[LINK BD]'
```

Available rds files (in *Seven Bridge*:  `Analysis/Tasks`):
- **Jul 16, 2025 15:40**: `ATACMultiomewithST-CC-additionalSMK_Seurat.rds`
- **May 12, 2025 13:54**: `ATACMultiomewithST-CC_Seurat.rds` (not good)
- **Apr 15, 2025 19:16**: `AC_WTA_SMK_index_library_Seurat.rds` (not good)

--> Use the last one (it got the extra SMK reads ("barcode reads"))

**Download fragment files**
```bash
wget -O output/seurat/ATACMultiomewithST-CC-additionalSMK_ATAC_Fragments.bed.gz \
'https://[LINK BD]'

wget -O output/seurat/ATACMultiomewithST-CC-additionalSMK_ATAC_Fragments.bed.gz.tbi \
'https://[LINK BD]'
```
- **Jul 16, 2025 15:40**: `ATACMultiomewithST-CC-additionalSMK_ATAC_Fragments.bed.gz`



# Data analysis -V1



```bash
conda activate SignacV5
module load hdf5
```

```R
set.seed(42)

# library
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()

library("Signac")
library("Seurat")
library("tidyverse")

# load seurat object
ATACMultiomewithST_SMK <- readRDS(file = "output/seurat/ATACMultiomewithST-CC-additionalSMK_Seurat.rds")

# Add mit and rib genes information
ATACMultiomewithST_SMK[["percent.mt"]] <- PercentageFeatureSet(ATACMultiomewithST_SMK, pattern = "^mt-")
ATACMultiomewithST_SMK[["percent.rb"]] <- PercentageFeatureSet(ATACMultiomewithST_SMK, pattern = "^Rp[sl]")

##################################
# QC RNA assay ###################
##################################


# QC plots
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nFeature_RNA.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nFeature_RNA_02000.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,2000)
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nFeature_RNA-groupSample_Name.pdf", width = 40, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nFeature_RNA"), ncol = 4, pt.size = 0.1, group.by= "Sample_Name") 
dev.off()

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nCount_RNA.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nCount_RNA_010000.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,10000)
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nCount_RNA_05000.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1) + ylim(0,5000)
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-nCount_RNA-groupSample_Name.pdf", width = 40, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1, group.by= "Sample_Name") 
dev.off()

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-percent.mt.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("percent.mt"), ncol = 4, pt.size = 0.1) 
dev.off()

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-percent.rb.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()



# Filter seurat object

ATACMultiomewithST_SMK_QCv1 <- subset(
  ATACMultiomewithST_SMK,
  subset =
    nFeature_RNA > 250 &
    nFeature_RNA < 3000 &
    nCount_RNA > 250 &
    nCount_RNA < 7500 &
    percent.mt < 10 &
    percent.rb < 5 &
    Sample_Name != "Multiplet" &
    Sample_Name != "Undetermined"
)
#--> 27,173 to 23,322 cells


pdf("output/seurat/UMAP-ATACMultiomewithST_SMK_QCv1-default.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_QCv1, reduction = "umap", group.by= "Sample_Name" , label=FALSE)
dev.off()
#--> Sample integration is bad; strong experimental effect




##################################
# QC peaks (ATAC) assay ###################
##################################

# Path to your fragments file
fragpath <- "output/seurat/ATACMultiomewithST-CC-additionalSMK_ATAC_Fragments.bed.gz"
# Update ATAC assay with correct file path
ATACMultiomewithST_SMK_QCv1[['peaks']]@fragments <- list(
  CreateFragmentObject(
    path = fragpath,
    cellnames = colnames(ATACMultiomewithST_SMK_QCv1),
    validate.fragments = TRUE
  )
)

DefaultAssay(ATACMultiomewithST_SMK_QCv1) <- "peaks" # 

# compute nucleosome signal score per cell
ATACMultiomewithST_SMK_QCv1 <- NucleosomeSignal(object = ATACMultiomewithST_SMK_QCv1)
# compute TSS enrichment score per cell
ATACMultiomewithST_SMK_QCv1 <- TSSEnrichment(object = ATACMultiomewithST_SMK_QCv1, fast = FALSE)


# QC plots
pdf("output/seurat/DensityScatter-ATACMultiomewithST_SMK_QCv1.pdf", width=5, height=5)
DensityScatter(ATACMultiomewithST_SMK_QCv1, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
#--> Bottom red line around TSS.enrichment =2 (test this treshold for high low)

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_QCv1-nCount_peaks.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nCount_peaks'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_QCv1-nCount_peaks_030000.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nCount_peaks'),
  pt.size = 0.1,
  ncol = 5 ) + ylim(0,30000)
dev.off()

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_QCv1-TSSenrichment.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_QCv1-nucleosome_signal.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()




ATACMultiomewithST_SMK_QCv1$high.tss <- ifelse(ATACMultiomewithST_SMK_QCv1$TSS.enrichment > 3, 'High', 'Low') # 3 is comonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot-ATACMultiomewithST_SMK_QCv1-TSSenrichment3.pdf", width=12, height=6)
TSSPlot(ATACMultiomewithST_SMK_QCv1, group.by = 'high.tss') + NoLegend()
dev.off()

ATACMultiomewithST_SMK_QCv1$high.tss <- ifelse(ATACMultiomewithST_SMK_QCv1$TSS.enrichment > 2, 'High', 'Low') # 3 is comonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot-ATACMultiomewithST_SMK_QCv1-TSSenrichment2.pdf", width=12, height=6)
TSSPlot(ATACMultiomewithST_SMK_QCv1, group.by = 'high.tss') + NoLegend()
dev.off()

#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment; let's pick TSS.enrichment > 2



## subset cells that pass QC
ATACMultiomewithST_SMK_QCv2 <- subset(
  x = ATACMultiomewithST_SMK_QCv1,
  subset = nCount_peaks > 100 &
    nCount_peaks < 25000 &
    TSS.enrichment > 2 &
    TSS.enrichment < 11 &
    nucleosome_signal > 0.11 &
    nucleosome_signal < 0.28
)
#--> 23,322 to 22,921 cells
ATACMultiomewithST_SMK_QCv2





###############################################################
# Sample integration with SCT #################################
###############################################################
# Individualize sample

Normoxia_F <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Normoxia_F") # 5098
Normoxia_M <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Normoxia_M") # 5284

Hypoxia_F <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Hypoxia_F") # 7412
Hypoxia_M <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Hypoxia_M") # 5127


# Add sex, condition metrics
Normoxia_F$sex <- "Female"
Normoxia_M$sex <- "Male"
Hypoxia_F$sex <- "Female"
Hypoxia_M$sex <- "Male"
Normoxia_F$condition <- "Normoxia"
Normoxia_M$condition <- "Normoxia"
Hypoxia_F$condition <- "Hypoxia"
Hypoxia_M$condition <- "Hypoxia"

DefaultAssay(Normoxia_F) <- "RNA"
DefaultAssay(Normoxia_M) <- "RNA"
DefaultAssay(Hypoxia_F) <- "RNA"
DefaultAssay(Hypoxia_M) <- "RNA"


# SCTransform
Normoxia_F <- SCTransform(Normoxia_F, method = "glmGamPoi", ncells = 5098, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Normoxia_M <- SCTransform(Normoxia_M, method = "glmGamPoi", ncells = 5284, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Hypoxia_F <- SCTransform(Hypoxia_F, method = "glmGamPoi", ncells = 7412, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Hypoxia_M <- SCTransform(Hypoxia_M, method = "glmGamPoi", ncells = 5127, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
#--> glmGamPoi not found but affect only speed; so it s ok!



# SCT integration
srat.list <- list(Normoxia_F = Normoxia_F, Normoxia_M = Normoxia_M, Hypoxia_F = Hypoxia_F, Hypoxia_M = Hypoxia_M)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
ATACMultiomewithST_SMK_QCv2.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
ATACMultiomewithST_SMK_QCv2.sct <- IntegrateData(anchorset = ATACMultiomewithST_SMK_QCv2.anchors, normalization.method = "SCT")


#### UMAP
DefaultAssay(ATACMultiomewithST_SMK_QCv2.sct) <- "integrated"

ATACMultiomewithST_SMK_QCv2.sct <- RunPCA(ATACMultiomewithST_SMK_QCv2.sct, verbose = FALSE, npcs = 50)
ATACMultiomewithST_SMK_QCv2.sct <- RunUMAP(ATACMultiomewithST_SMK_QCv2.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
ATACMultiomewithST_SMK_QCv2.sct <- FindNeighbors(ATACMultiomewithST_SMK_QCv2.sct, reduction = "pca", k.param = 30, dims = 1:50)
ATACMultiomewithST_SMK_QCv2.sct <- FindClusters(ATACMultiomewithST_SMK_QCv2.sct, resolution = 0.3, verbose = FALSE, algorithm = 4) # method = "igraph" needed for large nb of cells


ATACMultiomewithST_SMK_QCv2.sct$condition <- factor(ATACMultiomewithST_SMK_QCv2.sct$condition, levels = c("Normoxia", "Hypoxia")) # 

pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_QCv2-dim50kparam30res03.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=TRUE)
dev.off()

pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_QCv2-dim50kparam30res03-groupSex.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=TRUE, group.by= "sex" )
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_QCv2-dim50kparam30res03-groupSex.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=FALSE, group.by= "sex" )
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_QCv2-dim50kparam30res03-groupCondition.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=FALSE, group.by= "condition" )
dev.off()

###########################################################################
# saveRDS(ATACMultiomewithST_SMK_QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_QCv2-dim50kparam30res03.rds") 
ATACMultiomewithST_SMK_QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_QCv2-dim50kparam30res03.rds")
###########################################################################



# From the annotation.xlsx file
NSC= c("Mki67", "Top2a", "Nrxn3", "Dlx6os1")
NPC= c("Eomes", "Slc1a3", "Sema3a", "Slc17a6")
In1= c("Gad2", "Rbfox1", "Nrxn3", "Lhx6", "Satb2", "Sema3a")
In2= c("Gad2", "Sst", "Nrxn3", "Lhx6", "Sema3a", "Reln")
In3= c("Gad2", "Adarb2", "Nrxn3", "Dlx6os1", "Sema3a", "Slc17a6", "Reln")
In4= c("Gad2", "Isl1", "Nrxn3", "Sema3a", "Ldb2")
In5= c("Gad2", "Drd2", "Tle4")
Glut1= c("Neurod2", "Neurod6", "Rbfox1", "Tle4", "Zfpm2", "Ldb2", "Rbfox3", "Adarb2", "Nrxn3")
Glut2= c("Neurod2", "Neurod6", "Rbfox1", "Satb2", "Dok5", "Ldb2", "Nrxn3", "Tle4", "Rbfox3")
Glut3= c("Neurod2", "Rbfox1", "Unc5d", "Sema3c","Satb2", "Slc17a6", "Neurod6", "Sema3a")
Glut4= c("Rbfox1", "Rbfox3", "Zfpm2", "Sema3a", "Slc17a6", "Neurod2", "Neurod6")
Glut5= c("Rbfox1", "Ldb2", "Dok5", "Satb2", "Reln")
Glut6= c("Rbfox3", "Slc1a3")
Neuron= c("Reln", "Slc17a6")
OPC= c("Olig1", "Olig2", "Pdgfra", "Sox10", "Zfpm2", "Slc1a3")
RG= c("Mki67", "Top2a", "Zfpm2", "Fabp7", "Slc1a3", "Tnc")
Mg= c("C1qa", "C1qb")
Astro= c("Slc1a3", "Fabp7", "Zfpm2", "Nrxn3")
Endo= c("Col4a1", "Col4a2", "Cldn5", "Zfpm2")


# Summary main genes
NSC= c("Mki67", "Top2a", "Nrxn3", "Dlx6os1")
NPC= c("Eomes", "Slc1a3", "Sema3a", "Slc17a6")
In= c("Gad2", "Rbfox1", "Nrxn3", "Lhx6", "Satb2", "Sema3a", "Sst",  "Reln", "Adarb2",  "Dlx6os1", "Sema3a", "Slc17a6", "Isl1", "Ldb2", "Drd2", "Tle4")
Glut= c("Neurod2", "Neurod6", "Rbfox1", "Tle4", "Zfpm2", "Ldb2", "Rbfox3", "Adarb2", "Nrxn3", "Unc5d", "Sema3c","Satb2", "Slc17a6", "Dok5", "Sema3a", "Reln", "Slc1a3")
Neuron= c("Reln", "Slc17a6")
OPC= c("Olig1", "Olig2", "Pdgfra", "Sox10", "Zfpm2", "Slc1a3")
RG= c("Mki67", "Top2a", "Zfpm2", "Fabp7", "Slc1a3", "Tnc")
Mg= c("C1qa", "C1qb")
Astro= c("Slc1a3", "Fabp7", "Zfpm2", "Nrxn3")
Endo= c("Col4a1", "Col4a2", "Cldn5", "Zfpm2")


XXXY HERE !!!



DefaultAssay(ATACMultiomewithST_SMK_QCv2.sct) <- "SCT"

pdf("output/seurat/FeaturePlot-dim50-ListdotPLot.pdf", width=30, height=60)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = c(   "Gabra6","Pax6", # Granular_1
  "Sorcs3", "Ptprk", # MLI1
  "Nxph1", "Cdh22", # MLI2
  "Klhl1", "Gfra2", "Aldh1a3", # PLI12
  "Galntl6", "Kcnc2", # PLI23
  "Pax2", # Golgi
  "Eomes", # Unipolar_Brush
  "Calb1", "Slc1a6", "Car8", # Purkinje
  "Zeb2", # Astrocyte
  "Aqp4", "Slc39a12", # Bergmann_Glia
  "Mbp", "Mag", "Plp1", # Oligodendrocyte
  "Aldoc", "Cnp", # OPC
  "Itgam", "Cx3cr1", # Mix_Microglia_Meningeal
  "Ptgds", "Dcn", # Endothelial
  "Lef1", "Notum", "Apcdd1", # Endothelial_Mural
  "Dlc1", "Pdgfrb", # Choroid_Plexus
  "Kl",  "Ttr",
  "Eomes", "Rgs6", "Tafa2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()



















```



- NOTE: BD seurat object is partially filtered (ie. Min number of `nCount_RNA` and `nFeature_RNA`) + doublet already identified in `Sample_Name`
    --> Confirmed by XXX






