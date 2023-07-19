# Overview

Collaboration with Estaras lab for scRNAseq analysis. 

Two scRNAseq datasets:
- **Human gastruloids** using RUES hESCs treated with DASATINIB (YAP Inhibitor) (two samples: Untreated + DASATINIB) at 72hr timepoint 
- E7.75 Control vs cYAPKO **Embryos** (two samples: Control + cYAPKO) 

Objectives:
- Overall QC of the data (nb of RNA/genes detected, doublet? RNA conta? % mitochondrial genes)
- Cleaning and quick/raw clustering of each samples inficivually; make sure we are able to detect the 3 main expected clusters (endoderm, epiderm, mesoderm)
--> Send update to Estaras lab
- Put data together and identify common marker genes in control and treated sample in each scRNAseq datasets
- Perform clean clustering to have the 3 main expected clusters (endoderm, epiderm, mesoderm); in each datasets and conditions
- Perform DEGs analysis to identify, upregulated and downregulated in each cluster 
- Confirm expected YAP1 target are deregulated
--> Send update to Estaras lab


# Quick 1st analysis; data per data

Let's process the data individually as if each sample were individual, unrelated, scRNAseq experiment; in order to make sure we identify the 3 main cluster in each samples.

## Counting with cellranger count

```bash 
conda activate scRNAseq
which cellranger

# Run counting sample per sample
sbatch scripts/cellranger_count_humangastruloid_UNTREATED72hr.sh # 2367676 ok (~12hrs with 500g)
sbatch scripts/cellranger_count_humangastruloid_DASATINIB72hr.sh # 2367686 ok (~12hrs with 500g)
sbatch scripts/cellranger_count_embryo_control.sh # 2368763 ok (~24hrs with 300g)
sbatch scripts/cellranger_count_embryo_cYAPKO.sh # 2368769 ok (~24hrs with 300g)
```
--> Run succesfullly

## RNA contamination and doublet detection
- doublet detection using [scrublet](https://github.com/swolock/scrublet) **on the filtered matrix**
- ambient RNA correction using `soupX` in R before generating the Seurat object

```bash
srun --mem=500g --pty bash -l
conda deactivate # base environment needed
python3 scrublet.py [input_path] [output_path]
# Run doublet detection/scrublet sample per sample
python3 scripts/scrublet_doublets.py humangastruloid_UNTREATED72hr/outs/filtered_feature_bc_matrix output/doublets/humangastruloid_UNTREATED72hr.tsv
python3 scripts/scrublet_doublets.py humangastruloid_DASATINIB72hr/outs/filtered_feature_bc_matrix output/doublets/humangastruloid_DASATINIB72hr.tsv
```
Doublet detection score:
- humangastruloid_UNTREATED72hr: 0% doublet
- humangastruloid_DASATINIB72hr: 34.2% doublet
- XXX
--> Successfully assigned doublet

# humangastruloid analysis in Seurat

Then `conda activate scRNAseq` and go into R and filtered out **RNA contamination and start with SEURAT**.

We will apply SCT transform V2 for normalization; as [here](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html)

```R
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


# soupX decontamination
## Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('humangastruloid_UNTREATED72hr/outs') 
sc = load10X('humangastruloid_DASATINIB72hr/outs') 

## Assess % of conta
pdf("output/soupX/autoEstCont_humangastruloid_UNTREATED72hr.pdf", width=10, height=10)
pdf("output/soupX/autoEstCont_humangastruloid_DASATINIB72hr.pdf", width=10, height=10)
sc = autoEstCont(sc)
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/out_humangastruloid_UNTREATED72hr.RData")
save(out, file = "output/soupX/out_humangastruloid_DASATINIB72hr.RData")
## Load the matrix and Create SEURAT object
load("output/soupX/out_humangastruloid_UNTREATED72hr.RData")
srat_UNTREATED72hr <- CreateSeuratObject(counts = out, project = "UNTREATED72hr") # 36,601 features across 6,080 samples

load("output/soupX/out_humangastruloid_DASATINIB72hr.RData")
srat_DASATINIB72hr <- CreateSeuratObject(counts = out, project = "DASATINIB72hr") # 36,601 features across 9,938 samples

# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
srat_UNTREATED72hr[["percent.mt"]] <- PercentageFeatureSet(srat_UNTREATED72hr, pattern = "^MT-")
srat_UNTREATED72hr[["percent.rb"]] <- PercentageFeatureSet(srat_UNTREATED72hr, pattern = "^RP[SL]")

srat_DASATINIB72hr[["percent.mt"]] <- PercentageFeatureSet(srat_DASATINIB72hr, pattern = "^MT-")
srat_DASATINIB72hr[["percent.rb"]] <- PercentageFeatureSet(srat_DASATINIB72hr, pattern = "^RP[SL]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/humangastruloid_UNTREATED72hr.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_UNTREATED72hr <- AddMetaData(srat_UNTREATED72hr,doublets)
srat_UNTREATED72hr$Doublet_score <- as.numeric(srat_UNTREATED72hr$Doublet_score) # make score as numeric
head(srat_UNTREATED72hr[[]])

doublets <- read.table("output/doublets/humangastruloid_DASATINIB72hr.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_DASATINIB72hr <- AddMetaData(srat_DASATINIB72hr,doublets)
srat_DASATINIB72hr$Doublet_score <- as.numeric(srat_DASATINIB72hr$Doublet_score) # make score as numeric
head(srat_DASATINIB72hr[[]])

## Plot
pdf("output/seurat/VlnPlot_QC_UNTREATED72hr.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QC_DASATINIB72hr.pdf", width=10, height=6)
VlnPlot(srat_DASATINIB72hr, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_UNTREATED72hr.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_DASATINIB72hr.pdf", width=5, height=5)
FeatureScatter(srat_DASATINIB72hr, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_UNTREATED72hr.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_DASATINIB72hr.pdf", width=5, height=5)
FeatureScatter(srat_DASATINIB72hr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
pdf("output/seurat/FeatureScatter_QC_rb_mt_UNTREATED72hr.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_rb_mt_DASATINIB72hr.pdf", width=5, height=5)
FeatureScatter(srat_DASATINIB72hr, feature1 = "percent.rb", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_mt_doublet_UNTREATED72hr.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_mt_doublet_DASATINIB72hr.pdf", width=5, height=5)
FeatureScatter(srat_DASATINIB72hr, feature1 = "percent.mt", feature2 = "Doublet_score")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_UNTREATED72hr.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_DASATINIB72hr.pdf", width=5, height=5)
FeatureScatter(srat_DASATINIB72hr, feature1 = "nFeature_RNA", feature2 = "Doublet_score")
dev.off()
# After seeing the plot; add QC information in our seurat object
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 500 & srat_UNTREATED72hr@meta.data$QC == 'Pass','Low_nFeature',srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 500 & srat_UNTREATED72hr@meta.data$QC != 'Pass' & srat_UNTREATED72hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_UNTREATED72hr@meta.data$QC,sep = ','),srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.mt > 15 & srat_UNTREATED72hr@meta.data$QC == 'Pass','High_MT',srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 500 & srat_UNTREATED72hr@meta.data$QC != 'Pass' & srat_UNTREATED72hr@meta.data$QC != 'High_MT',paste('High_MT',srat_UNTREATED72hr@meta.data$QC,sep = ','),srat_UNTREATED72hr@meta.data$QC)
table(srat_UNTREATED72hr[['QC']])

## OPTIONAL FILTERING OF LOW RIBO CONTENT CELL
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.rb < 5 & srat_UNTREATED72hr@meta.data$QC == 'Pass', 'Low_ribo', srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.rb < 5 & srat_UNTREATED72hr@meta.data$QC != 'Pass' & srat_UNTREATED72hr@meta.data$QC != 'Low_ribo', paste('Low_ribo', srat_UNTREATED72hr@meta.data$QC, sep = ','), srat_UNTREATED72hr@meta.data$QC)
table(srat_UNTREATED72hr[['QC']])
##


srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$nFeature_RNA < 500 & srat_DASATINIB72hr@meta.data$QC == 'Pass','Low_nFeature',srat_DASATINIB72hr@meta.data$QC)
srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$nFeature_RNA < 500 & srat_DASATINIB72hr@meta.data$QC != 'Pass' & srat_DASATINIB72hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_DASATINIB72hr@meta.data$QC,sep = ','),srat_DASATINIB72hr@meta.data$QC)
srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$percent.mt > 15 & srat_DASATINIB72hr@meta.data$QC == 'Pass','High_MT',srat_DASATINIB72hr@meta.data$QC)
srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$nFeature_RNA < 500 & srat_DASATINIB72hr@meta.data$QC != 'Pass' & srat_DASATINIB72hr@meta.data$QC != 'High_MT',paste('High_MT',srat_DASATINIB72hr@meta.data$QC,sep = ','),srat_DASATINIB72hr@meta.data$QC)
table(srat_DASATINIB72hr[['QC']])
# Quality check after QC filtering
pdf("output/seurat/VlnPlot_QCPass_UNTREATED72hr.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QCPass_DASATINIB72hr.pdf", width=10, height=6)
VlnPlot(subset(srat_DASATINIB72hr, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()



## subset my seurat object to only analyze the cells that pass the QC
srat_UNTREATED72hr <- subset(srat_UNTREATED72hr, subset = QC == 'Pass')
srat_DASATINIB72hr <- subset(srat_DASATINIB72hr, subset = QC == 'Pass')
srat_UNTREATED72hr$condition <- "UNTREATED72hr"
srat_DASATINIB72hr$condition <- "DASATINIB72hr"

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srat_UNTREATED72hr <- CellCycleScoring(srat_UNTREATED72hr, s.features = s.genes, g2m.features = g2m.genes)
table(srat_UNTREATED72hr[[]]$Phase)


# Test different normalization method until we found one that cluster well endoderm, mesoderm and ectoderm in control

## normalize and run dimensionality reduction on control dataset with SCTransform V2
### Default parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE)
### Fine-tune parameters (from v1)
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, vst.flavor = "v2", method = "glmGamPoi", ncells = 5487,     # HERE PASTE NB OF CELLS AFTER QC!!! NB OF CEL IN THE srat
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 15, dims = 1:30, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) without regression
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, vst.flavor = "v2", method = "glmGamPoi", ncells = 5487, verbose = F) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 15, dims = 1:30, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression only
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, vst.flavor = "v2", method = "glmGamPoi", ncells = 5487, vars.to.regress = c("percent.mt"), verbose = F) %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 15, dims = 1:50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.7, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5487, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 15, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.7, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5487, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F) %>%
    RunPCA(npcs = 25, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:25, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:25, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.2, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with n features chagnes
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5487, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 15, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:15, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with UMAP parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5487, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 15, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:15, n.neighbors = 50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:15, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with UMAP parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5487, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 15, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:15, spread = 2, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:15, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with UMAP parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5487, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 15, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 50, dims = 1:15, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with UMAP parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5487, vars.to.regress = c("percent.mt"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 40, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with UMAP parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5847, vars.to.regress = c("percent.mt","S.Score","G2M.Score","percent.rb"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 10, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.1, verbose = FALSE, algorithm = 4)

### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with UMAP parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5847, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.2, verbose = FALSE, algorithm = 4)


pdf("output/seurat/UMAP_SCT_UNTREATED72hr.pdf", width=10, height=6)
DimPlot(srat_UNTREATED72hr, label = T, repel = T) + ggtitle("UNTREATED72hr_Unsupervised clustering")
dev.off()

## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "NES", "SOX2", "SOX9", "BRN2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "GSC", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP", "SOX17")

DefaultAssay(srat_UNTREATED72hr) <- "SCT"
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_endoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_ectoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
# --> NOT AMAZING with Default; DIM 10,30,50 AND RES 0.3 0.5 TESTED
# --> Tested with parameters used in the non v2 VST transformation and perform far best!


# other parameters

## percent mt
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_mt.pdf", width=10, height=10)
FeaturePlot(srat_UNTREATED72hr,features = "percent.mt",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()
## percent rb
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_rb.pdf", width=10, height=10)
FeaturePlot(srat_UNTREATED72hr,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_rb.pdf", width=10, height=10)
VlnPlot(srat_UNTREATED72hr,features = c("percent.rb")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## RNA
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_count_feature.pdf", width=10, height=10)
VlnPlot(srat_UNTREATED72hr,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## cell cycle
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_cellCycle.pdf", width=10, height=10)
FeaturePlot(srat_UNTREATED72hr,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()  










# Let's try now the non SCTransform method

## LogNormalize ##
srat_UNTREATED72hr <- NormalizeData(srat_UNTREATED72hr, normalization.method = "LogNormalize", scale.factor = 10000)
## Discover the 2000 first more variable genes
srat_UNTREATED72hr <- FindVariableFeatures(srat_UNTREATED72hr, selection.method = "vst", nfeatures = 3000)
## scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(srat_UNTREATED72hr)
srat_UNTREATED72hr <- ScaleData(srat_UNTREATED72hr, features = all.genes)
## PCA
srat_UNTREATED72hr <- RunPCA(srat_UNTREATED72hr, features = VariableFeatures(object = srat_UNTREATED72hr))
## investigate to find optimal nb of dimension
### Vizualize the 1st 20 PC
pdf("output/seurat/DimHeatmap_LogNormalize.pdf", width=10, height=40)
DimHeatmap(srat_UNTREATED72hr, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

### Compute JackStraw score
srat_UNTREATED72hr <- JackStraw(srat_UNTREATED72hr, num.replicate = 100)
srat_UNTREATED72hr <- ScoreJackStraw(srat_UNTREATED72hr, dims = 1:30)

## Generate plot
pdf("output/seurat/JackStraw_LogNormalize.pdf", width=10, height=10)
JackStrawPlot(srat_UNTREATED72hr, dims = 1:20)  # 20
dev.off()
## --> Jack drop at 17
pdf("output/seurat/Elbow_LogNormalize.pdf", width=10, height=10)
ElbowPlot(srat_UNTREATED72hr) # 8 
dev.off()
## --> Elbow drop at 4-7 or 10


# clustering
srat_UNTREATED72hr <- FindNeighbors(srat_UNTREATED72hr, dims = 1:30, k.param = 50)
srat_UNTREATED72hr <- FindClusters(srat_UNTREATED72hr, resolution = 0.2)
srat_UNTREATED72hr <- RunUMAP(srat_UNTREATED72hr, dims = 1:30, verbose = F)
table(srat_UNTREATED72hr@meta.data$seurat_clusters) # to check cluster size

pdf("output/seurat/Umap_LogNormalize.pdf", width=10, height=10)
DimPlot(srat_UNTREATED72hr,label.size = 4,repel = T,label = T)
dev.off()


## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "NES", "SOX2", "SOX9", "BRN2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "GSC", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP", "SOX17")

pdf("output/seurat/FeaturePlot_LogNormalize_UNTREATED72hr_endoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_UNTREATED72hr_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_UNTREATED72hr_ectoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()






## CLR ##
srat_UNTREATED72hr <- NormalizeData(srat_UNTREATED72hr, normalization.method = "CLR", scale.factor = 10000)
## Discover the 2000 first more variable genes
srat_UNTREATED72hr <- FindVariableFeatures(srat_UNTREATED72hr, selection.method = "vst", nfeatures = 2000)
## scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(srat_UNTREATED72hr)
srat_UNTREATED72hr <- ScaleData(srat_UNTREATED72hr, features = all.genes)
## PCA
srat_UNTREATED72hr <- RunPCA(srat_UNTREATED72hr, features = VariableFeatures(object = srat_UNTREATED72hr))

# clustering
srat_UNTREATED72hr <- FindNeighbors(srat_UNTREATED72hr, dims = 1:30, k.param = 30)
srat_UNTREATED72hr <- FindClusters(srat_UNTREATED72hr, resolution = 0.2)
srat_UNTREATED72hr <- RunUMAP(srat_UNTREATED72hr, dims = 1:30, verbose = F)
table(srat_UNTREATED72hr@meta.data$seurat_clusters) # to check cluster size

pdf("output/seurat/Umap_CLR.pdf", width=10, height=10)
DimPlot(srat_UNTREATED72hr,label.size = 4,repel = T,label = T)
dev.off()


## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "NES", "SOX2", "SOX9", "BRN2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "GSC", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP", "SOX17")

pdf("output/seurat/FeaturePlot_CLR_UNTREATED72hr_endoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_CLR_UNTREATED72hr_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_CLR_UNTREATED72hr_ectoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()



## RC ##
srat_UNTREATED72hr <- NormalizeData(srat_UNTREATED72hr, normalization.method = "RC", scale.factor = 10000)
## Discover the 2000 first more variable genes
srat_UNTREATED72hr <- FindVariableFeatures(srat_UNTREATED72hr, selection.method = "vst", nfeatures = 2000)
## scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(srat_UNTREATED72hr)
srat_UNTREATED72hr <- ScaleData(srat_UNTREATED72hr, features = all.genes)
## PCA
srat_UNTREATED72hr <- RunPCA(srat_UNTREATED72hr, features = VariableFeatures(object = srat_UNTREATED72hr))

# clustering
srat_UNTREATED72hr <- FindNeighbors(srat_UNTREATED72hr, dims = 1:30, k.param = 30)
srat_UNTREATED72hr <- FindClusters(srat_UNTREATED72hr, resolution = 0.2)
srat_UNTREATED72hr <- RunUMAP(srat_UNTREATED72hr, dims = 1:30, verbose = F)
table(srat_UNTREATED72hr@meta.data$seurat_clusters) # to check cluster size

pdf("output/seurat/Umap_RC.pdf", width=10, height=10)
DimPlot(srat_UNTREATED72hr,label.size = 4,repel = T,label = T)
dev.off()


## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "NES", "SOX2", "SOX9", "BRN2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "GSC", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP", "SOX17")

pdf("output/seurat/FeaturePlot_RC_UNTREATED72hr_endoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_RC_UNTREATED72hr_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_RC_UNTREATED72hr_ectoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()




## LogNormalize disp ##
srat_UNTREATED72hr <- NormalizeData(srat_UNTREATED72hr, normalization.method = "LogNormalize", scale.factor = 10000)
## Discover the 2000 first more variable genes
srat_UNTREATED72hr <- FindVariableFeatures(srat_UNTREATED72hr, selection.method = "disp", nfeatures = 2000)
## scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(srat_UNTREATED72hr)
srat_UNTREATED72hr <- ScaleData(srat_UNTREATED72hr, features = all.genes)
## PCA
srat_UNTREATED72hr <- RunPCA(srat_UNTREATED72hr, features = VariableFeatures(object = srat_UNTREATED72hr))



# clustering
srat_UNTREATED72hr <- FindNeighbors(srat_UNTREATED72hr, dims = 1:30, k.param = 30)
srat_UNTREATED72hr <- FindClusters(srat_UNTREATED72hr, resolution = 0.2)
srat_UNTREATED72hr <- RunUMAP(srat_UNTREATED72hr, dims = 1:30, verbose = F)
table(srat_UNTREATED72hr@meta.data$seurat_clusters) # to check cluster size

pdf("output/seurat/Umap_LogNormalize_disp.pdf", width=10, height=10)
DimPlot(srat_UNTREATED72hr,label.size = 4,repel = T,label = T)
dev.off()


## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "NES", "SOX2", "SOX9", "BRN2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "GSC", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP", "SOX17")

pdf("output/seurat/FeaturePlot_LogNormalize_disp_UNTREATED72hr_endoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_disp_UNTREATED72hr_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_disp_UNTREATED72hr_ectoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()





## LogNormalize fine tuning ##
srat_UNTREATED72hr <- NormalizeData(srat_UNTREATED72hr, normalization.method = "LogNormalize", scale.factor = 10000)
## Discover the 2000 first more variable genes
srat_UNTREATED72hr <- FindVariableFeatures(srat_UNTREATED72hr, selection.method = "vst", nfeatures = 2000)
## scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(srat_UNTREATED72hr)
srat_UNTREATED72hr <- ScaleData(srat_UNTREATED72hr, features = all.genes)
## PCA
srat_UNTREATED72hr <- RunPCA(srat_UNTREATED72hr, features = VariableFeatures(object = srat_UNTREATED72hr))



# clustering
srat_UNTREATED72hr <- FindNeighbors(srat_UNTREATED72hr, dims = 1:5, k.param = 30)
srat_UNTREATED72hr <- FindClusters(srat_UNTREATED72hr, resolution = 0.2)
srat_UNTREATED72hr <- RunUMAP(srat_UNTREATED72hr, dims = 1:5, verbose = F)
table(srat_UNTREATED72hr@meta.data$seurat_clusters) # to check cluster size

pdf("output/seurat/Umap_LogNormalize.pdf", width=10, height=10)
DimPlot(srat_UNTREATED72hr,label.size = 4,repel = T,label = T)
dev.off()


## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "NES", "SOX2", "SOX9", "BRN2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "GSC", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP", "SOX17")

pdf("output/seurat/FeaturePlot_LogNormalize_UNTREATED72hr_endoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_UNTREATED72hr_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_UNTREATED72hr_ectoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()




















## integrate the two datasets using pearson residuals
srat_DASATINIB72hr <- SCTransform(srat_DASATINIB72hr, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)
srat.list <- list(srat_UNTREATED72hr = srat_UNTREATED72hr, srat_DASATINIB72hr = srat_DASATINIB72hr)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

humangastruloid.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
humangastruloid.combined.sct <- IntegrateData(anchorset = humangastruloid.anchors, normalization.method = "SCT")



# Perform integrated analysis
humangastruloid.combined.sct <- RunPCA(humangastruloid.combined.sct, verbose = FALSE)
humangastruloid.combined.sct <- RunUMAP(humangastruloid.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
humangastruloid.combined.sct <- FindNeighbors(humangastruloid.combined.sct, reduction = "pca", dims = 1:30)
humangastruloid.combined.sct <- FindClusters(humangastruloid.combined.sct, resolution = 0.3)



pdf("output/seurat/UMAP_raw_UNTREATED72hr_DASATINIB72hr.pdf", width=10, height=6)
DimPlot(humangastruloid.combined.sct, reduction = "umap", split.by = "condition")
dev.off()



### BACKUP ### TOO LONG , fuck it
# Save workspace to .RData file
save.image("output/seurat/R_session_humangastruloid.RData")
# Load workspace from .RData file
load("output/seurat/R_session_humangastruloid.RData")
### BACKUP ###




# differential expressed genes across conditions


XXX




















# normalize data
srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)

# Discover the 2000 first more variable genes
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(srat), 10)
top10 

## plot the 10 first variable
pdf("output/seurat/VariableFeaturePlot_top10_50dOrga.pdf", width=10, height=10)
plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
dev.off()

# scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)

# PCA
srat <- RunPCA(srat, features = VariableFeatures(object = srat))

# investigate to find optimal nb of dimension

## Vizualize the 1st 20 PC
pdf("output/seurat/DimHeatmap.pdf", width=10, height=20)
DimHeatmap(srat, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

## Compute JackStraw score
srat <- JackStraw(srat, num.replicate = 100)
srat <- ScoreJackStraw(srat, dims = 1:20)

## Generate plot
pdf("output/seurat/JackStraw.pdf", width=10, height=10)
JackStrawPlot(srat, dims = 1:20)  # 20
dev.off()

pdf("output/seurat/Elbow.pdf", width=10, height=10)
ElbowPlot(srat) # 8 
dev.off()

# clustering
srat <- FindNeighbors(srat, dims = 1:10)
srat <- FindClusters(srat, resolution = 0.1)

srat <- RunUMAP(srat, dims = 1:10, verbose = F)
table(srat@meta.data$seurat_clusters) # to check cluster size

pdf("output/seurat/Umap.pdf", width=10, height=10)
DimPlot(srat,label.size = 4,repel = T,label = T)
dev.off()

# Check minor cell population to make sure clustering is good
pdf("output/seurat/FeaturePlot_Papercluster4.pdf", width=10, height=5)
FeaturePlot(srat, features = c("XXX"))
dev.off()

# Check parameter to see if one clsuter to do technical stuff
pdf("output/seurat/FeaturePlot_Doublet_score.pdf", width=10, height=10)
FeaturePlot(srat, features = "Doublet_score") & theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat/FeaturePlot_mt.pdf", width=10, height=10)
FeaturePlot(srat, features = "percent.mt") & theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat/FeaturePlot_RNAfeature.pdf", width=10, height=10)
FeaturePlot(srat, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))
dev.off()

# Remove cells that did not path QC
srat <- subset(srat, subset = QC == 'Pass')

pdf("output/seurat/Umap_QCPass.pdf", width=10, height=10)
DimPlot(srat,label.size = 4,repel = T,label = T)
dev.off()

# Check cell cycles
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)
table(srat[[]]$Phase)

# Let's see if cluster define by technical difference (even after QC filtering)

## percent mt
pdf("output/seurat/FeaturePlot_QCPass_mt.pdf", width=10, height=10)
FeaturePlot(srat,features = "percent.mt",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat/VlnPlot_QCPass_mt.pdf", width=10, height=10)
VlnPlot(srat,features = "percent.mt") & theme(plot.title = element_text(size=10))
dev.off()
## percent rb
pdf("output/seurat/FeaturePlot_QCPass_rb.pdf", width=10, height=10)
FeaturePlot(srat,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
dev.off()
pdf("output/seurat/VlnPlot_QCPass_rb.pdf", width=10, height=10)
VlnPlot(srat,features = "percent.rb") & theme(plot.title = element_text(size=10))
dev.off()
## RNA
pdf("output/seurat/VlnPlot_QCPass_RNA.pdf", width=10, height=10)
VlnPlot(srat,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## cell cycle
pdf("output/seurat/FeaturePlot_QCPass_cellCycle.pdf", width=10, height=10)
FeaturePlot(srat,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()  
pdf("output/seurat/VlnPlot_QCPass_cellCycle.pdf", width=10, height=10)
VlnPlot(srat,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()



# SCTransform normalization and clustering
dim(srat)[2] # to know nb of cells
srat <- SCTransform(srat, method = "glmGamPoi", ncells = 20862,     # HERE PASTE NB OF CELLS AFTER QC!!! NB OF CEL IN THE srat
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F) 



# SCTransform PCA 
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, k.param = 15, dims = 1:30, verbose = F) # k.param can be changed: A larger value will create more connections between cells and may result in less distinct clustering
srat <- FindClusters(srat, verbose = F, resolution = 0.1, algorithm = 4) # algorithm can be changed
table(srat[[]]$seurat_clusters)


## Plot
pdf("output/seurat/Umap_QCPass_SCT.pdf", width=10, height=10)
DimPlot(srat, label = T)
dev.off()

## Check marker genes
pdf("output/seurat/FeaturePlot_Papercluster_muscle_myogenicProgenitors_SCT.pdf", width=10, height=10)
FeaturePlot(srat, features = c("PAX7", "PITX2", "MSC","MKI67","MYF5","SYTL2"))
dev.off()
pdf("output/seurat/FeaturePlot_Papercluster_muscle_myocytes_SCT.pdf", width=10, height=10)
FeaturePlot(srat, features = c("MYMX", "SGCA", "MYOD1", "MYOG", "CDH15", "CHRNB1"))
dev.off()
pdf("output/seurat/FeaturePlot_Papercluster_muscle_skeletalMuscleFibers_SCT.pdf", width=10, height=10)
FeaturePlot(srat, features = c("CHRND", "TTN", "MYBPH", "ACTA1", "ACTN2", "UNC45B"))
dev.off()
pdf("output/seurat/FeaturePlot_Papercluster_neural_meuralProgenitors_SCT.pdf", width=10, height=10)
FeaturePlot(srat, features = c("PAX6", "HES5", "HOPX", "FABP7"))
dev.off()
pdf("output/seurat/FeaturePlot_Papercluster_neural_immatureNeurons_SCT.pdf", width=10, height=10)
FeaturePlot(srat, features = c("HES6", "DLL3", "ASCL1"))
dev.off()
pdf("output/seurat/FeaturePlot_Papercluster_neural_matureSpinalCordNeurons_SCT.pdf", width=10, height=10)
FeaturePlot(srat, features = c("NHLH1", "ELAVL4", "ELAVL3", "STMN2", "NSG2"))
dev.off()
pdf("output/seurat/FeaturePlot_Papercluster_neural_neuralCrestDerivatives_SCT.pdf", width=10, height=10)
FeaturePlot(srat, features = c("TWIST1", "COL5A1", "SERPINF1", "PDGFRA", "PDGFRB"))
dev.off()
pdf("output/seurat/FeaturePlot_Papercluster_neural_neural_gliaCells_SCT.pdf", width=10, height=10)
FeaturePlot(srat, features = c("GFAP", "DKK1", "EFNB3", "SCG2"))
dev.off()




# Save the seurat
save(srat, file = "output/seurat/srat.RData")


# Differential expression and marker selection
## set back  active assay to “RNA” and re-do norm and scaling (as we filter out many cells)
DefaultAssay(srat) <- "RNA"
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)

## DEGs
all.markers <- FindAllMarkers(srat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers, file = "output/seurat/all_markers.tsv", sep = "\t", quote = F, row.names = T)


## check markers
dim(all.markers)
table(all.markers$cluster)
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers


write.table(as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)), file = "output/seurat/top20_markers.tsv", sep = "\t", quote = F, row.names = T)



# Cell type annotation using SingleR
## Get reference datasets
monaco.ref <- celldex::MonacoImmuneData()
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()

## Convert Seurat object into SCE
sce <- as.SingleCellExperiment(DietSeurat(orga5d_QCPass))
sce

## Analyse
monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)

## check results
table(monaco.main$pruned.labels)
table(hpca.main$pruned.labels)
table(dice.main$pruned.labels)
table(monaco.fine$pruned.labels)
table(hpca.fine$pruned.labels)
table(dice.fine$pruned.labels)

# Add annotation to seurat object
orga5d_QCPass@meta.data$hpca.main   <- hpca.main$pruned.labels
orga5d_QCPass@meta.data$hpca.fine   <- hpca.fine$pruned.labels

orga5d_QCPass@meta.data$monaco.main <- monaco.main$pruned.labels
orga5d_QCPass@meta.data$monaco.fine <- monaco.fine$pruned.labels

orga5d_QCPass@meta.data$dice.main   <- dice.main$pruned.labels
orga5d_QCPass@meta.data$dice.fine   <- dice.fine$pruned.labels

# cluster with name
pdf("output/seurat_sanger/UMAP_SCT_annotated_hpca.main.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "hpca.main")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_hpca.fine.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "hpca.fine")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_monaco.main.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "monaco.main")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_monaco.fine.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "monaco.fine")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_dice.main.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "dice.main")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat_sanger/UMAP_SCT_annotated_dice.fine.pdf", width=10, height=10)
orga5d_QCPass <- SetIdent(orga5d_QCPass, value = "dice.fine")
DimPlot(orga5d_QCPass, label = T , repel = T, label.size = 3)
dev.off()


# manual annotation
# Rename cluster in agreement
## Re run this to re-generate SCTransform

srat <- SCTransform(srat, method = "glmGamPoi", ncells = 20862,     # HERE PASTE NB OF CELLS AFTER QC!!! NB OF CEL IN THE srat
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F) 

# SCTransform PCA 
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, k.param = 15, dims = 1:30, verbose = F) # k.param can be changed: A larger value will create more connections between cells and may result in less distinct clustering
srat <- FindClusters(srat, verbose = F, resolution = 0.1, algorithm = 4) # algorithm can be changed
table(srat[[]]$seurat_clusters)




new.cluster.ids <- c("Neural_NeuralCrestDerivatives_1", "Neural_ImmatureNeurons_MatureSpinalCords", "Neural_NeuralCrestDerivatives_2", "Neural_GliaCells", "Muscle_MyogenicProgenitors", "Muscle_Myocytes", "Unknown_1", "Neural_NeuralProgenitors", "Muscle_SkeletalMuscleFibers", "Unknown_2")
names(new.cluster.ids) <- levels(srat)
srat <- RenameIdents(srat, new.cluster.ids)

pdf("output/seurat/UMAP_label_manual_fromSCTransform.pdf", width=10, height=10)
DimPlot(srat, reduction = "umap", label = TRUE, pt.size = 0.7, label.size = 6) + NoLegend()
dev.off()
##
```



For human gastruloid, the optimal parameter seems to be with the following:

XXX check from ppt XXX
