# Project

Test the BD Rhapsody pipeline for scRNAseq. Compare result with already processed data.

Design:
- Normoxia vs hypoxia at fetal (neonatal) stage
    - 2 Bio rep per condition (1 male, 1 female)


--> Previous marker genes used by Ana can be found at `011*/docs/Annotation.xlsx`

--> BD Rhapsody work with mm39 (not mm10!)

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



# Data analysis - V1 - not stringeant QC



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
#--> 27,173 to 23,322 cells; First test: not enough clean; some cluster form becasue of high nCount_RNA; and nFeature_RNA and percent.mt


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
ATACMultiomewithST_SMK_QCv2.sct <- FindClusters(ATACMultiomewithST_SMK_QCv2.sct, resolution = 0.2, verbose = FALSE, algorithm = 4) # method = "igraph" needed for large nb of cells


ATACMultiomewithST_SMK_QCv2.sct$condition <- factor(ATACMultiomewithST_SMK_QCv2.sct$condition, levels = c("Normoxia", "Hypoxia")) # 

pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_QCv2-dim50kparam30res01.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=TRUE)
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_QCv2-dim50-groupSex.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=TRUE, group.by= "sex" )
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_QCv2-dim50-groupCondition.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=FALSE, group.by= "condition" )
dev.off()

pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_QCv2-dim50-nFeature_RNA.pdf", width=5, height=5)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_QCv2-dim50-nCount_RNA.pdf", width=5, height=5)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=FALSE, features = "nCount_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_QCv2-dim50-percentmt.pdf", width=5, height=5)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_QCv2-dim50-percentrb.pdf", width=5, height=5)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, reduction = "umap", label=FALSE, features = "percent.rb")
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
Astro= c("Slc1a3", "Fabp7", "Zfpm2", "Nrxn3", "Tnc")
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
Astro= c("Slc1a3", "Fabp7", "Zfpm2", "Nrxn3", "Tnc")
Endo= c("Col4a1", "Col4a2", "Cldn5", "Zfpm2")

DefaultAssay(ATACMultiomewithST_SMK_QCv2.sct) <- "SCT"

pdf("output/seurat/FeaturePlot-dim50-NSC.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = NSC, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-NPC.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = NPC, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-In.pdf", width=25, height=10)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = In, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-Glut.pdf", width=30, height=20)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = Glut, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-Neuron.pdf", width=15, height=5)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = Neuron, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-OPC.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = OPC, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-RG.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = RG, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-Mg.pdf", width=15, height=5)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = Mg, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-Astro.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = Astro, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-dim50-Endo.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = Endo, max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot-dim50-AstroGoldberg.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_QCv2.sct, features = c("Aqp4", "Slc39a12"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

#--> More cleaning needed; some cluster forming because of high RNA and mit content

```









# Data analysis - V2 - more stringeant QC



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
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-percent.mt-groupSample_Name.pdf", width = 40, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("percent.mt"), ncol = 4, pt.size = 0.1, group.by= "Sample_Name") 
dev.off()

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-percent.rb.pdf", width = 10, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("percent.rb"), ncol = 4, pt.size = 0.1) 
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK-percent.rb-groupSample_Name.pdf", width = 40, height = 6)
VlnPlot(ATACMultiomewithST_SMK, features = c("percent.rb"), ncol = 4, pt.size = 0.1, group.by= "Sample_Name") 
dev.off()


# Filter seurat object

ATACMultiomewithST_SMK_QCv1 <- subset(
  ATACMultiomewithST_SMK,
  subset =
    nFeature_RNA > 250 &
    nFeature_RNA < 2000 &
    nCount_RNA > 250 &
    nCount_RNA < 5000 &
    percent.mt < 7.5 &
    percent.rb < 5 &
    Sample_Name != "Multiplet" &
    Sample_Name != "Undetermined"
)
#--> 27,173 to 21,685 cells



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
pdf("output/seurat/DensityScatter-ATACMultiomewithST_SMK_V2QCv1.pdf", width=5, height=5)
DensityScatter(ATACMultiomewithST_SMK_QCv1, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
#--> Bottom red line around TSS.enrichment =2 (test this treshold for high low)

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-nCount_peaks.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nCount_peaks'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-nCount_peaks-groupSample_Name.pdf", width = 40, height = 6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nCount_peaks'),
  pt.size = 0.1,
  ncol = 5, group.by= "Sample_Name" )
dev.off()


pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-nCount_peaks_030000.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nCount_peaks'),
  pt.size = 0.1,
  ncol = 5 ) + ylim(0,30000)
dev.off()

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-TSSenrichment.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-TSSenrichment-groupSample_Name.pdf", width = 40, height = 6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('TSS.enrichment'),
  pt.size = 0.1,
  ncol = 5, group.by= "Sample_Name" )
dev.off()

pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-nucleosome_signal.pdf", width=12, height=6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()
pdf("output/seurat/VlnPlot-QC-ATACMultiomewithST_SMK_V2QCv1-nucleosome_signal-groupSample_Name.pdf", width = 40, height = 6)
VlnPlot(
  object = ATACMultiomewithST_SMK_QCv1,
  features = c('nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5, group.by= "Sample_Name" )
dev.off()



ATACMultiomewithST_SMK_QCv1$high.tss <- ifelse(ATACMultiomewithST_SMK_QCv1$TSS.enrichment > 3, 'High', 'Low') # 3 is comonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot-ATACMultiomewithST_SMK_QCv1-TSSenrichment3.pdf", width=12, height=6)
TSSPlot(ATACMultiomewithST_SMK_QCv1, group.by = 'high.tss') + NoLegend()
dev.off()

ATACMultiomewithST_SMK_QCv1$high.tss <- ifelse(ATACMultiomewithST_SMK_QCv1$TSS.enrichment > 2, 'High', 'Low') # 3 is comonly used but could be adjusted; based on violin and below plot
pdf("output/seurat/TSSPlot-ATACMultiomewithST_SMK_QCv1-TSSenrichment2.pdf", width=12, height=6)
TSSPlot(ATACMultiomewithST_SMK_QCv1, group.by = 'high.tss') + NoLegend()
dev.off()

#--> Let's keep 2



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
#--> 21,685 to 21,333 cells
ATACMultiomewithST_SMK_QCv2




###############################################################
# Sample integration with SCT #################################
###############################################################
# Individualize sample

Normoxia_F <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Normoxia_F") # 4488
Normoxia_M <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Normoxia_M") # 5178

Hypoxia_F <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Hypoxia_F") # 6881
Hypoxia_M <- subset(ATACMultiomewithST_SMK_QCv2, subset = Sample_Name == "Hypoxia_M") # 4786


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
Normoxia_F <- SCTransform(Normoxia_F, method = "glmGamPoi", ncells = 4488, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Normoxia_M <- SCTransform(Normoxia_M, method = "glmGamPoi", ncells = 5178, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Hypoxia_F <- SCTransform(Hypoxia_F, method = "glmGamPoi", ncells = 6881, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
Hypoxia_M <- SCTransform(Hypoxia_M, method = "glmGamPoi", ncells = 4786, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "nFeature_RNA")) 
#--> glmGamPoi not found but affect only speed; so it s ok!



# SCT integration
srat.list <- list(Normoxia_F = Normoxia_F, Normoxia_M = Normoxia_M, Hypoxia_F = Hypoxia_F, Hypoxia_M = Hypoxia_M)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
ATACMultiomewithST_SMK_V2QCv2.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)

ATACMultiomewithST_SMK_V2QCv2.sct <- IntegrateData(anchorset = ATACMultiomewithST_SMK_V2QCv2.anchors, normalization.method = "SCT")


#### UMAP
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "integrated"

ATACMultiomewithST_SMK_V2QCv2.sct <- RunPCA(ATACMultiomewithST_SMK_V2QCv2.sct, verbose = FALSE, npcs = 30)
ATACMultiomewithST_SMK_V2QCv2.sct <- RunUMAP(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
ATACMultiomewithST_SMK_V2QCv2.sct <- FindNeighbors(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "pca", k.param = 30, dims = 1:30)
ATACMultiomewithST_SMK_V2QCv2.sct <- FindClusters(ATACMultiomewithST_SMK_V2QCv2.sct, resolution = 0.4, verbose = FALSE, algorithm = 4) # method = "igraph" needed for large nb of cells


ATACMultiomewithST_SMK_V2QCv2.sct$condition <- factor(ATACMultiomewithST_SMK_V2QCv2.sct$condition, levels = c("Normoxia", "Hypoxia")) # 

pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=TRUE)
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30-groupSex.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=TRUE, group.by= "sex" )
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30-groupCondition.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=FALSE, group.by= "condition" )
dev.off()

pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-nFeature_RNA.pdf", width=5, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-nCount_RNA.pdf", width=5, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=FALSE, features = "nCount_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-percentmt.pdf", width=5, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-percentrb.pdf", width=5, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", label=FALSE, features = "percent.rb")
dev.off()  


# mitochondria content histogram


pdf("output/seurat/geom_histogram-ATACMultiomewithST_SMK_V2QCv2-dim30-BINpercentrb.pdf", width=5, height=5)
ggplot(ATACMultiomewithST_SMK_V2QCv2.sct@meta.data, aes(x = percent.mt)) +
  geom_histogram(binwidth = 0.5) +   # 1% mt bins; change to 0.5, 2, etc.
  labs(x = "Percent mitochondrial reads (percent.mt)",
       y = "Cell count",
       title = "Cell count by mitochondrial content") +
  theme_classic()
dev.off()  



###########################################################################
# saveRDS(ATACMultiomewithST_SMK_V2QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.rds") 
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.rds")
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
Astro= c("Slc1a3", "Fabp7", "Zfpm2", "Nrxn3", "Tnc")
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
Astro= c("Slc1a3", "Fabp7", "Zfpm2", "Nrxn3", "Tnc")
Endo= c("Col4a1", "Col4a2", "Cldn5", "Zfpm2")

DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "SCT"

pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-NSC.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = NSC, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-NPC.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = NPC, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-In.pdf", width=25, height=10)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = In, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Glut.pdf", width=30, height=20)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = Glut, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Neuron.pdf", width=15, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = Neuron, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-OPC.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = OPC, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-RG.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = RG, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Mg.pdf", width=15, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = Mg, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Astro.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = Astro, max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Endo.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = Endo, max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-AstroGoldberg.pdf", width=10, height=10)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = c("Aqp4", "Slc39a12"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Neurod2.pdf", width=6, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = c("Neurod2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Neurod6.pdf", width=6, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = c("Neurod6"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Gad2.pdf", width=6, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = c("Gad2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Nrxn3.pdf", width=6, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = c("Nrxn3"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-C1qb.pdf", width=6, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = c("C1qb"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-dim30-Mki67.pdf", width=6, height=5)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = c("Mki67"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()





#-->  dim30kparam30res04 seems the best option!


################################################
# cluster naming _ dim30kparam30res04 ##########
################################################

Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "seurat_clusters"


Cluster1= Glut1= Neurod2, Neurod6, Rbfox1
Cluster2= Glut2= Neurod2, Neurod6, Rbfox1
Cluster3= Glut3= Neurod2, Neurod6, Rbfox1
Cluster4= Glut4= Neurod2, Neurod6, Rbfox1
Cluster5= NSC= Mki67, Top2a, Dlx6os1
Cluster6= RG= Fabp7, Slc1a3, Tnc
Cluster7= In1= Gad2, Nrxn3, Lhx6
Cluster8= In2= Gad2, Nrxn3, Lhx6
Cluster9= In3= Gad2, Nrxn3, Lhx6
Cluster10= Glut5= Neurod2, Neurod6, Rbfox1
Cluster11= Glut6= Neurod2, Neurod6, Rbfox1
Cluster12= In4= Gad2, Nrxn3, Lhx6
Cluster13= Astro= Aqp4
Cluster14= NPC= Eomes
Cluster15= OPC= Olig1, Olig2, Pdgfra
Cluster16= Neuron= Reln, Slc17a6
Cluster17= Endo= Col4a1, Col4a2
Cluster18= Mg= C1qa, C1qb






new.cluster.ids <- c(
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "NSC",
  "RG",
  "In1",
  "In2",
  "In3",
  "Glut5",
  "Glut6",
  "In4",
  "Astro",
  "NPC",
  "OPC",
  "Neuron",
  "Endo",
  "Mg"
)

names(new.cluster.ids) <- levels(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- RenameIdents(ATACMultiomewithST_SMK_V2QCv2.sct, new.cluster.ids)
ATACMultiomewithST_SMK_V2QCv2.sct$cluster.annot <- Idents(ATACMultiomewithST_SMK_V2QCv2.sct) # create a new slot in my seurat object



pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-label.pdf", width=15, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 6)
dev.off()




pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelNoSplit.pdf", width=7, height=6)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5)
dev.off()


################################################
# marker genes expression _ dim30kparam30res04 #
################################################

# dotplot of main marker genes


all_markers <- c(
  "Aqp4",
  "Col4a1", "Col4a2",
  "Neurod2", "Neurod6", "Rbfox1",
  "Gad2", "Nrxn3", "Lhx6",
  "C1qa", "C1qb",
  "Reln", "Slc17a6",
  "Eomes",
  "Mki67", "Top2a", "Dlx6os1",
  "Olig1", "Olig2", "Pdgfra",
  "Fabp7", "Slc1a3", "Tnc"
)


markers_Cristancho <- c(
  "Neurod6", "Rbfox1", "Tle4", "Zfpm2", "Satb2", "Dok5", "Unc5d", "Sema3c", "Sema3a", "Nrxn3", "Reln", "Ldb2", "Gad2", "Lhx6", "Sst", "Adarb2", "Isl1", "Drd2", "Eomes", "Mki67", "Top2a", "Tnc", "Olig2", "Slc1a3", "Fabp7", "Col4a1", "Cldn5", "C1qb", "C1qa"
)



levels(ATACMultiomewithST_SMK_V2QCv2.sct) <- c(
  "Astro",
  "Endo",
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "Glut5",
  "Glut6",
  "In1",
  "In2",
  "In3",
  "In4",
  "Mg",
  "Neuron",
  "NPC",
  "NSC",
  "OPC",
  "RG"
)



pdf("output/seurat/DotPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-all_markers.pdf", width=8, height=4.5)
DotPlot(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "SCT", features = all_markers, cols = c("#1A00FF", "#FF0000")) + RotatedAxis() + scale_y_discrete(limits = rev)
dev.off()

pdf("output/seurat/DotPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-markersCristancho.pdf", width=9, height=4.5)
DotPlot(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "SCT", features = markers_Cristancho, cols = c("#1A00FF", "#FF0000")) + RotatedAxis() + scale_y_discrete(limits = rev)
dev.off()




################################################
# DEG _ dim30kparam30res04 #####################
################################################

# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "RNA"

ATACMultiomewithST_SMK_V2QCv2.sct <- NormalizeData(ATACMultiomewithST_SMK_V2QCv2.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- ScaleData(ATACMultiomewithST_SMK_V2QCv2.sct, features = all.genes) # zero-centres and scales it


## what genes change in different conditions for cells of the same type

ATACMultiomewithST_SMK_V2QCv2.sct$celltype.stim <- paste(ATACMultiomewithST_SMK_V2QCv2.sct$cluster.annot, ATACMultiomewithST_SMK_V2QCv2.sct$condition,
    sep = "-")
Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "celltype.stim"


## Automation::
cell_types <- c(  "Astro",
  "Endo",
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "Glut5",
  "Glut6",
  "In1",
  "In2",
  "In3",
  "In4",
  "Mg",
  "Neuron",
  "NPC",
  "NSC",
  "OPC",
  "RG")


for (cell_type in cell_types) {
  response_name <- paste(cell_type, "Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.response", sep = ".")
  ident_1 <- paste0(cell_type, "-Hypoxia")
  ident_2 <- paste0(cell_type, "-Normoxia")

  file_name <- paste0("output/seurat/", cell_type,
                      "-Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt")

  message("Running DE for cluster ", cell_type, " ...")

  tryCatch(
    {
      response <- FindMarkers(
        ATACMultiomewithST_SMK_V2QCv2.sct,
        assay = "RNA",
        ident.1 = ident_1,
        ident.2 = ident_2,
        verbose = FALSE
      )
      
      print(head(response, 15))
      write.table(response, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE)
    },
    error = function(e) {
      message("⚠️  Skipping cluster ", cell_type, ": ", conditionMessage(e))
    }
  )
}


#--> DEGs save as output/seurat/[CLUSTERNAME]-Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt





# DEGs number in dotplot
DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())
## List of cell types
cell_types <- c(   "Astro",
  "Endo",
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "Glut5",
  "Glut6",
  "In1",
  "In2",
  "In3",
  "In4",
  "Mg",
  "Neuron",
  "NPC",
  "NSC",
  "OPC",
  "RG")
## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/seurat/", cell_type, "-Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt", sep = "")
  deg_data <- read.table(file_name, header = TRUE, sep = "\t") ## Read the DEGs data
  num_degs <- sum(deg_data$p_val_adj < 0.05 & (deg_data$avg_log2FC > 0.25 | deg_data$avg_log2FC < -0.25)) ## Count the number of significant DEGs
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))  ## Append to the summary table
}
## Dotplot
DEG_count= DEG_count %>%
  mutate(Cell_Type = factor(Cell_Type, levels = c( 
  "Astro",
  "Endo",
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "Glut5",
  "Glut6",
  "In1",
  "In2",
  "In3",
  "In4",
  "Mg",
  "Neuron",
  "NPC",
  "NSC",
  "OPC",
  "RG") ) ) 
DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = rev(levels(DEG_count$Cell_Type)))

# Generate the dot plot
pdf("output/seurat/Dotplot-DEG_count-Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.pdf", width=9, height=4)
ggplot(DEG_count, aes(x = 1, y = Cell_Type, color = Cell_Type)) +
  geom_point(aes(size = Num_DEGs), alpha = 0.8) +
  scale_size(range = c(2, 10)) +
  theme_void() + # Removes all gridlines, axes, and background elements
  theme(
    axis.text.x = element_blank(), # Removes x-axis labels
    axis.ticks.x = element_blank(), # Removes x-axis ticks
    axis.title.x = element_blank(), # Removes x-axis title
    axis.text.y = element_text(size = 10, hjust = 0),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "left", # Keeps the legend visible
    legend.title = element_text(size = 12, face = "bold"), # Adjusts legend title style
    legend.text = element_text(size = 10) # Adjusts legend text style
  ) +
  labs(
    y = "Cell Type",
    color = "Cell Type",
    size = "Number of DEGs"
  )
dev.off()


# Output number up and down ############################################


DEG_count <- data.frame(
  Cell_Type = character(),
  Num_Up_DEGs = integer(),
  Num_Down_DEGs = integer(),
  DEG_total = integer()
)
cell_types <- c(
  "Astro","Endo","Glut1","Glut2","Glut3","Glut4","Glut5","Glut6",
  "In1","In2","In3","In4","Mg","Neuron","NPC","NSC","OPC","RG"
)


for (cell_type in cell_types) {
  file_name <- paste0(
    "output/seurat/",
    cell_type,
    "-Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt"
  )
  # Skip missing files
  if (!file.exists(file_name)) {
    message("⚠️ Missing: ", file_name, " — skipping")
    DEG_count <- rbind(
      DEG_count,
      data.frame(Cell_Type = cell_type, Num_Up_DEGs = 0,
                 Num_Down_DEGs = 0, DEG_total = 0)
    )
    next
  }
  # Read file
  deg_data <- read.table(file_name, header = TRUE, sep = "\t")
  # Count DEGs
  num_up   <- sum(deg_data$p_val_adj < 0.05 & deg_data$avg_log2FC >  0.25)
  num_down <- sum(deg_data$p_val_adj < 0.05 & deg_data$avg_log2FC < -0.25)
  num_total <- num_up + num_down
  # Store
  DEG_count <- rbind(
    DEG_count,
    data.frame(Cell_Type = cell_type,
               Num_Up_DEGs = num_up,
               Num_Down_DEGs = num_down,
               DEG_total = num_total)
  )
}

# Set factor order
DEG_count <- DEG_count %>%
  mutate(Cell_Type = factor(Cell_Type, levels = rev(cell_types)))
cat("Total Upregulated DEGs:", sum(DEG_count$Num_Up_DEGs), "\n")
cat("Total Downregulated DEGs:", sum(DEG_count$Num_Down_DEGs), "\n")
cat("Total DEGs (Up + Down):", sum(DEG_count$DEG_total), "\n")
DEG_count
########################################################



# DEGs number colored in a UMAP
Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "cluster.annot"

DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())
## List of cell types
cell_types <- c(     "Astro","Endo","Glut1","Glut2","Glut3","Glut4","Glut5","Glut6",
  "In1","In2","In3","In4","Mg","Neuron","NPC","NSC","OPC","RG")
## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/seurat/", cell_type, "-Hypoxia-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt", sep = "") # CHANGE FILE HERE
  deg_data <- read.table(file_name, header = TRUE, sep = "\t") ## Read the DEGs data
  num_degs <- sum(deg_data$p_val_adj < 0.05 & (deg_data$avg_log2FC > 0.25 | deg_data$avg_log2FC < -0.25)) ## Count the number of significant DEGs
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))  ## Append to the summary table
}

DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = c(    "ImmatureGranule",
  "Astro","Endo","Glut1","Glut2","Glut3","Glut4","Glut5","Glut6",
  "In1","In2","In3","In4","Mg","Neuron","NPC","NSC","OPC","RG")) 
  
  
# Add DEG information to my seurat object - DEG_count
cell_clusters <- ATACMultiomewithST_SMK_V2QCv2.sct@meta.data$cluster.annot
names(cell_clusters) <- rownames(ATACMultiomewithST_SMK_V2QCv2.sct@meta.data)
DEG_named_vector <- DEG_count$Num_DEGs[match(cell_clusters, DEG_count$Cell_Type)]
names(DEG_named_vector) <- names(cell_clusters)
# Integrate DEG values into the Seurat object
ATACMultiomewithST_SMK_V2QCv2.sct <- AddMetaData(ATACMultiomewithST_SMK_V2QCv2.sct, metadata = DEG_named_vector, col.name = "DEG")

# Add values on the heatmap
## Extract UMAP coordinates
umap_coordinates <- as.data.frame(ATACMultiomewithST_SMK_V2QCv2.sct@reductions$umap@cell.embeddings)
umap_coordinates$cluster <- ATACMultiomewithST_SMK_V2QCv2.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DEG_count %>% dplyr::rename( "cluster"="Cell_Type"))
## Create a UMAP plot colored by DEG values, with cluster DEG counts as text annotations
pdf("output/seurat/FeaturePlot-ATACMultiomewithST_SMK_V2QCv2-padj05fc025_numeric_dim30kparam30res04.pdf", width=6, height=6)
FeaturePlot(ATACMultiomewithST_SMK_V2QCv2.sct, features = "DEG", pt.size = 0.5, reduction = "umap") +
  scale_colour_viridis(option="magma") + # option="magma"
  geom_text(data = cluster_centers, aes(x = umap_1, y = umap_2, label = Num_DEGs), 
            size = 4, color = "red", fontface = "bold") 
dev.off()






#######################################################
# Unbiased cell type markers _ dim30kparam30res04 #####
#######################################################



# Unbiased cell type marker genes
Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "cluster.annot"
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "RNA"
ATACMultiomewithST_SMK_V2QCv2.sct <- NormalizeData(ATACMultiomewithST_SMK_V2QCv2.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- ScaleData(ATACMultiomewithST_SMK_V2QCv2.sct, features = all.genes) # zero-centres and scales it

all_markers <- FindAllMarkers(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/FindAllMarkers-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)




###########################################################################
# saveRDS(ATACMultiomewithST_SMK_V2QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1.rds") 
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1.rds")
###########################################################################






######################################################################################################
## Pre-processing ATAC _ dim30kparam30res04  ####################################################################
######################################################################################################



###  pre-processing and dimensional reductio
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "peaks"
ATACMultiomewithST_SMK_V2QCv2.sct <- RunTFIDF(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- FindTopFeatures(ATACMultiomewithST_SMK_V2QCv2.sct, min.cutoff = 'q0')
ATACMultiomewithST_SMK_V2QCv2.sct <- RunSVD(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- RunUMAP(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = 'lsi', dims = 2:40, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # We exclude the first dimension as this is typically correlated with sequencing depth


### WNN graph, representing a weighted combination of RNA and ATAC-seq modalities

ATACMultiomewithST_SMK_V2QCv2.sct <- FindMultiModalNeighbors(ATACMultiomewithST_SMK_V2QCv2.sct, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:40))
ATACMultiomewithST_SMK_V2QCv2.sct <- RunUMAP(ATACMultiomewithST_SMK_V2QCv2.sct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ATACMultiomewithST_SMK_V2QCv2.sct <- FindClusters(ATACMultiomewithST_SMK_V2QCv2.sct, graph.name = "wsnn", algorithm = 3, verbose = TRUE)

# plot gene expression, ATAC-seq, or WNN analysis
p1 <- DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", group.by = "cluster.annot", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap.atac", group.by = "cluster.annot", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "wnn.umap", group.by = "cluster.annot", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")

pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-RNAATACWNN.pdf", width=14, height=5)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()



ATACMultiomewithST_SMK_V2QCv2.sct$condition <- factor(ATACMultiomewithST_SMK_V2QCv2.sct$condition, levels = c("Normoxia", "Hypoxia")) # Reorder untreated 1st

pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-WNN.pdf", width=5, height=5)
p = DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "wnn.umap", group.by = "cluster.annot", label = FALSE, label.size = 4, repel = TRUE) + ggtitle("Cell type") + NoLegend()
LabelClusters(p, id = "cluster.annot", fontface = "bold", color = "black", size = 3)
dev.off()

pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-RNAcondition.pdf", width=5, height=5)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap", group.by = "condition", label = FALSE, cols = c("blue", "red")) + ggtitle("RNA")  + NoLegend()
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-ATACcondition.pdf", width=5, height=5)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "umap.atac", group.by = "condition", label = FALSE, cols = c("blue", "red")) + ggtitle("ATAC")  + NoLegend()
dev.off()
pdf("output/seurat/DimPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04_ATACdim240-WNNcondition.pdf", width=5, height=5)
DimPlot(ATACMultiomewithST_SMK_V2QCv2.sct, reduction = "wnn.umap", group.by = "condition", label = FALSE, cols = c("blue", "red")) + ggtitle("WNN")  + NoLegend()
dev.off()



#--> Different dim for ATAC tested; 2:40 best (notably for separating Endo and Mg)


# Calculate gene activity (count ATAC peak within gene and promoter)
GeneActivity = GeneActivity(
  ATACMultiomewithST_SMK_V2QCv2.sct,
  assay = "peaks",
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = NULL,
  max.width = NULL,
  process_n = 2000,
  gene.id = FALSE,
  verbose = TRUE
)
## Add gene activity as new assay

ATACMultiomewithST_SMK_V2QCv2.sct[["GeneActivity"]] <- CreateAssayObject(counts = GeneActivity)


###########################################################################
# saveRDS(ATACMultiomewithST_SMK_V2QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivity.rds") 
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivity.rds")
###########################################################################




# Marker genes GeneActivity


all_markers <- c(
  "Aqp4",
  "Col4a1", "Col4a2",
  "Neurod2", "Neurod6", "Rbfox1",
  "Gad2", "Nrxn3", "Lhx6",
  "C1qa", "C1qb",
  "Reln", "Slc17a6",
  "Eomes",
  "Mki67", "Top2a", "Dlx6os1",
  "Olig1", "Olig2", "Pdgfra",
  "Fabp7", "Slc1a3", "Tnc"
)


markers_ATAC_Cristancho <- c(
  "Neurod6", "Tle4", "Zfpm2", "Satb2", "Dok5", "Sema3c", "Sema3a",  "Reln", "Ldb2", "Gad2", "Lhx6", "Sst", "Isl1", "Drd2", "Eomes", "Mki67", "Top2a", "Tnc", "Olig2", "Slc1a3", "Fabp7", "Col4a1", "Cldn5", "C1qb", "C1qa"
)
# "Rbfox1", "Unc5d", "Nrxn3"  removed

Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "cluster.annot"


levels(ATACMultiomewithST_SMK_V2QCv2.sct) <- c(
  "Astro",
  "Endo",
  "Glut1",
  "Glut2",
  "Glut3",
  "Glut4",
  "Glut5",
  "Glut6",
  "In1",
  "In2",
  "In3",
  "In4",
  "Mg",
  "Neuron",
  "NPC",
  "NSC",
  "OPC",
  "RG"
)



pdf("output/seurat/DotPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-all_markers-GeneActivity.pdf", width=8, height=4.5)
DotPlot(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "GeneActivity", features = all_markers, cols = c("#0A7B83", "#E78F0B")) + RotatedAxis() + scale_y_discrete(limits = rev)
dev.off()

pdf("output/seurat/DotPlot-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-markers_ATAC_Cristancho-GeneActivity.pdf", width=9, height=4.5)
DotPlot(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "GeneActivity", features = markers_ATAC_Cristancho, cols = c("#0A7B83", "#E78F0B")) + RotatedAxis() + scale_y_discrete(limits = rev)
dev.off()






# Identify diff access regions (DAR)
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- 'peaks'

Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "celltype.stim"

### Here is to automatize the process clsuter per cluster:
#### Define the cluster pairs for comparison
cluster_pairs <- list(
  c("Astro-Normoxia",   "Astro-Hypoxia",   "Astro"),
  c("Endo-Normoxia",   "Endo-Hypoxia",    "Endo"),
  c("Glut1-Normoxia",  "Glut1-Hypoxia",   "Glut1"),
  c("Glut2-Normoxia",  "Glut2-Hypoxia",   "Glut2"),
  c("Glut3-Normoxia",  "Glut3-Hypoxia",   "Glut3"),
  c("Glut4-Normoxia",  "Glut4-Hypoxia",   "Glut4"),
  c("Glut5-Normoxia",  "Glut5-Hypoxia",   "Glut5"),
  c("Glut6-Normoxia",  "Glut6-Hypoxia",   "Glut6"),
  c("In1-Normoxia",    "In1-Hypoxia",     "In1"),
  c("In2-Normoxia",    "In2-Hypoxia",     "In2"),
  c("In3-Normoxia",    "In3-Hypoxia",     "In3"),
  c("In4-Normoxia",    "In4-Hypoxia",     "In4"),
  c("Mg-Normoxia",     "Mg-Hypoxia",      "Mg"),
  c("Neuron-Normoxia", "Neuron-Hypoxia",  "Neuron"),
  c("NPC-Normoxia",    "NPC-Hypoxia",     "NPC"),
  c("NSC-Normoxia",    "NSC-Hypoxia",     "NSC"),
  c("OPC-Normoxia",    "OPC-Hypoxia",     "OPC"),
  c("RG-Normoxia",     "RG-Hypoxia",      "RG")
)


## Function to run FindMarkers and return a tibble with cluster and gene information
run_find_markers <- function(pair) {
  ident_1 <- pair[1]
  ident_2 <- pair[2]
  cluster_name <- pair[3]
  
  # Run FindMarkers for each pair of clusters
  markers <- FindMarkers(
    object = ATACMultiomewithST_SMK_V2QCv2.sct,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = 'wilcox',
    min.pct = 0.1
  ) %>%
  rownames_to_column(var = "query_region") %>%  # Add gene names as a column
  add_column(cluster = cluster_name)    # Add cluster name as a new column
  
  return(markers)
}

## Combine all FindMarkers results into one tibble
all_markers_tibble <- map_dfr(cluster_pairs, run_find_markers) %>%
  as_tibble()
## save output
write.table(all_markers_tibble, file = "output/seurat/DAR_peaks-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


## Find closest gene to DAR
# PREREQUISETE #################
# Download genome file mm39 and gene version v106 - work with BD Rhapsody
### Link peaks to genes - Find peaks that are correlated with the expression of nearby genes
#### Install the BS genome (GC content, region lengths, and dinucleotide base frequencies for regions in the assay and add to the feature metadata)
# BiocManager::install(c("AnnotationHub"))

library("AnnotationHub")
library("ensembldb")

ah <- AnnotationHub()

# Query EnsDb mouse + Ensembl release 106 (GRCm39/mm39-era)
hits <- query(ah, c("EnsDb", "Mus musculus", "Ensembl", "106"))
hits
EnsDb.Mmusculus.v106 <- ah[["AH100674"]]   # Replace with correct names()
EnsDb.Mmusculus.v106
library("BSgenome.Mmusculus.UCSC.mm39") # BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")

## Add genomic information to Seurat object
### Convert rownames of ATAC counts to GRanges
grange.counts <- StringToGRanges(rownames(ATACMultiomewithST_SMK_V2QCv2.sct[["peaks"]]), sep = c(":", "-"))
### Filter for standard chromosomes
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- ATACMultiomewithST_SMK_V2QCv2.sct[as.vector(grange.use), ]
### Get annotations for the mouse genome
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v106) # EnsDb.Mmusculus.v79 if using mm10
### Adjust chromosome naming style
seqlevelsStyle(annotations) <- 'UCSC'
###################################################

# Identify closest gene to DAR
DAR_closestGene <- ClosestFeature(ATACMultiomewithST_SMK_V2QCv2.sct, regions = all_markers_tibble$query_region, annotation = annotations)
DAR_genes = all_markers_tibble %>% 
  left_join(DAR_closestGene) %>%
  as_tibble() %>%
  unique()

## save output
write.table(DAR_genes, file = "output/seurat/DAR_genes-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########### COUNT UP DOWN PEAK PER CLUSTER ##################################################################
all_markers_tibble = read_tsv("output/seurat/DAR_peaks-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt")

all_markers_tibble %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarize(
    nb_upregulated = sum(avg_log2FC > 0.58),
    nb_downregulated = sum(avg_log2FC < -0.58)
  )

DAR_genes = read_tsv("output/seurat/DAR_genes-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt")

dar_gene_summary <- DAR_genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select(cluster, gene_name, avg_log2FC) %>%
  dplyr::filter(!is.na(gene_name), gene_name != "NA") %>%          # optional cleanup
  group_by(cluster, gene_name) %>%
  summarise(
    any_open  = any(avg_log2FC >  0.58),
    any_close = any(avg_log2FC < -0.58),
    .groups = "drop"
  ) %>%
  mutate(
    class = case_when(
      any_open & any_close ~ "Mix",
      any_open             ~ "Opening",
      any_close            ~ "Closing",
      TRUE                 ~ "NoStrongChange"
    )
  )

## Count genes once per cluster, by class
dar_counts <- dar_gene_summary %>%
  dplyr::filter(class != "NoStrongChange") %>%                     # keep only strong changes, optional
  count(cluster, class) %>%
  tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0) %>%
  mutate(total = Opening + Closing + Mix)

dar_counts



##############################################################################################################
# Link peaks to genes - Find peaks that are correlated with the expression of nearby genes

#### Troubleshooting section ####################################
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "peaks"
ATACMultiomewithST_SMK_V2QCv2.sct <- RegionStats(
  object = ATACMultiomewithST_SMK_V2QCv2.sct,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  assay = "peaks",
  verbose = TRUE
)
gr <- granges(ATACMultiomewithST_SMK_V2QCv2.sct[["peaks"]])
c(max_chr4 = max(end(gr[seqnames(gr)=="chr4"])),
  mm10_chr4 = seqlengths(BSgenome.Mmusculus.UCSC.mm10)["chr4"])
gr <- granges(ATACMultiomewithST_SMK_V2QCv2.sct[["peaks"]])
c(max_chr4 = max(end(gr[seqnames(gr)=="chr4"])),
  mm10_chr4 = seqlengths(BSgenome.Mmusculus.UCSC.mm10)["chr4"],
  mm39_chr4 = seqlengths(BSgenome.Mmusculus.UCSC.mm39)["chr4"])

#--> Some peaks, notably at chr4 are invalid: occur after the end of chr4! Probably issue with BD Rhapsody calling peak at extremities... 
#--> Using mm39 solve the issue!!!
################################################################


ATACMultiomewithST_SMK_V2QCv2.sct <- RegionStats(
  object = ATACMultiomewithST_SMK_V2QCv2.sct,
  genome = BSgenome.Mmusculus.UCSC.mm39,
  assay = "peaks",
  verbose = TRUE
)


## Run  LinkPeaks
ATACMultiomewithST_SMK_V2QCv2.sct = LinkPeaks(
  ATACMultiomewithST_SMK_V2QCv2.sct,
  peak.assay = "peaks",
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  genes.use = NULL,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  gene.id = FALSE,
  verbose = TRUE
)




###########################################################################
# saveRDS(ATACMultiomewithST_SMK_V2QCv2.sct, file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks.rds") 
ATACMultiomewithST_SMK_V2QCv2.sct <- readRDS(file = "output/seurat/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks.rds")
###########################################################################



Links = as_tibble(Links(ATACMultiomewithST_SMK_V2QCv2.sct))



### Adjust the pvalue and select positive corr as in the https://www.nature.com/articles/s41467-024-45199-x#Fig2 paper
Links$adjusted_pvalue <- p.adjust(Links$pvalue, method = "BH")
#write.table(Links, file = c("output/seurat/LinkPeaks-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04.txt"),sep="\t", quote=FALSE, row.names=FALSE)

### Isolate signif genes
Links_signif = Links %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%  
  dplyr::select(gene) %>% # 1,942 Link Signif
  unique()
### Reorder the gene based on their cluster max expression

Links %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%  
  dplyr::select(peak) %>% # 3,740 Link peak Signif
  unique()

## log norm and Scale GeneActivity assay
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "GeneActivity"

ATACMultiomewithST_SMK_V2QCv2.sct <- NormalizeData(ATACMultiomewithST_SMK_V2QCv2.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(ATACMultiomewithST_SMK_V2QCv2.sct)
ATACMultiomewithST_SMK_V2QCv2.sct <- ScaleData(ATACMultiomewithST_SMK_V2QCv2.sct, features = all.genes) # zero-centres and scales it



###### Find all markers 
Idents(ATACMultiomewithST_SMK_V2QCv2.sct) <- "cluster.annot"
DefaultAssay(ATACMultiomewithST_SMK_V2QCv2.sct) <- "RNA"

Links_markers <- FindAllMarkers(ATACMultiomewithST_SMK_V2QCv2.sct, features = Links_signif$gene, assay = "RNA", only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
###### Identify in which cluster the Links_markers gene is highly express
Links_markers_pval= as_tibble(Links_markers) %>%
  group_by(gene) %>%
  dplyr::filter(p_val == min(p_val)) %>%
  dplyr::select(gene, cluster)
#write.table(Links_markers, file = "output/Signac/ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks-Links_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# plot heatmap
pdf("output/seurat/DoHeatmap-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks-LinksPadj05Score0-SCTscaledata.pdf", width=8, height=4)
DoHeatmap(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "SCT", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE)
dev.off()
pdf("output/seurat/DoHeatmap-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks-LinksPadj05Score0-GeneActivityscaledata.pdf", width=8, height=4)
DoHeatmap(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE)
dev.off()

pdf("output/seurat/DoHeatmap-ATACMultiomewithST_SMK_V2QCv2-dim30kparam30res04-labelV1GeneActivityLinkPeaks-LinksPadj05Score0-GeneActivityscaledata1.pdf", width=8, height=4)
DoHeatmap(ATACMultiomewithST_SMK_V2QCv2.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE, disp.max = 1.25, disp.min = -2)
dev.off()

# , disp.max = 2, disp.min = -2)








```





- NOTE: BD seurat object is partially filtered (ie. Min number of `nCount_RNA` and `nFeature_RNA`) + doublet already identified in `Sample_Name`

--> RNA and ATAC marker genes show very comparable expression as found by Ana, BD Rhapsody seems promising to me!













