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

*NOTE: data downloaded in local from basespace and then imported into the cluster*

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
# Re-run count for embryo using mice genome
cd ../meta
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
sbatch scripts/cellranger_count_embryo_control_mice.sh # 2982903 ok
sbatch scripts/cellranger_count_embryo_cYAPKO_mice.sh # 2982920 ok
```
--> Run succesfullly; re-run embryo with mice genome... (same folder but label as `*_mice`)

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
python3 scripts/scrublet_doublets.py embryo_control_e775_mice/outs/filtered_feature_bc_matrix output/doublets/embryo_control.tsv
python3 scripts/scrublet_doublets.py embryo_cYAPKO_e775_mice/outs/filtered_feature_bc_matrix output/doublets/embryo_cYAPKO.tsv
```
Doublet detection score:
- humangastruloid_UNTREATED72hr: 0% doublet
- humangastruloid_DASATINIB72hr: 34.2% doublet
- embryo_control: 0.1% doublet
- embryo_cYAPKO: 3.6%
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
## After seeing the plot; add QC information in our seurat object
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 500 & srat_UNTREATED72hr@meta.data$QC == 'Pass','Low_nFeature',srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 500 & srat_UNTREATED72hr@meta.data$QC != 'Pass' & srat_UNTREATED72hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_UNTREATED72hr@meta.data$QC,sep = ','),srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.mt > 15 & srat_UNTREATED72hr@meta.data$QC == 'Pass','High_MT',srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 500 & srat_UNTREATED72hr@meta.data$QC != 'Pass' & srat_UNTREATED72hr@meta.data$QC != 'High_MT',paste('High_MT',srat_UNTREATED72hr@meta.data$QC,sep = ','),srat_UNTREATED72hr@meta.data$QC)
table(srat_UNTREATED72hr[['QC']])
## 



## OPTIONAL FILTERING OF LOW RIBO CONTENT CELL
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.rb < 5 & srat_UNTREATED72hr@meta.data$QC == 'Pass', 'Low_ribo', srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.rb < 5 & srat_UNTREATED72hr@meta.data$QC != 'Pass' & srat_UNTREATED72hr@meta.data$QC != 'Low_ribo', paste('Low_ribo', srat_UNTREATED72hr@meta.data$QC, sep = ','), srat_UNTREATED72hr@meta.data$QC)
table(srat_UNTREATED72hr[['QC']])
##
## V2 QC with low RNA count filter out
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 1000 & srat_UNTREATED72hr@meta.data$QC == 'Pass', 'Low_nFeature', srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.mt > 15 & srat_UNTREATED72hr@meta.data$QC == 'Pass', 'High_MT', srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nCount_RNA < 5000 & srat_UNTREATED72hr@meta.data$QC == 'Pass', 'Low_nCount', srat_UNTREATED72hr@meta.data$QC)
# Check other QC tags if current QC is not Pass
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 1000 & srat_UNTREATED72hr@meta.data$QC != 'Pass', paste(srat_UNTREATED72hr@meta.data$QC, 'Low_nFeature', sep = ','), srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.mt > 15 & srat_UNTREATED72hr@meta.data$QC != 'Pass', paste(srat_UNTREATED72hr@meta.data$QC, 'High_MT', sep = ','), srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nCount_RNA < 5000 & srat_UNTREATED72hr@meta.data$QC != 'Pass', paste(srat_UNTREATED72hr@meta.data$QC, 'Low_nCount', sep = ','), srat_UNTREATED72hr@meta.data$QC)

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
VlnPlot(subset(srat_UNTREATED72hr, subset = QC == 'Pass'), 
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
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5847, vars.to.regress = c("percent.mt","percent.rb"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with UMAP parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5847, verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.3, verbose = FALSE, algorithm = 4)
### Fine-tune parameters (from v1) with percent mt regression and not vst.flavor v2 with UMAP parameter
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, vst.flavor = "v2", method = "glmGamPoi", ncells = 5738, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.2, verbose = FALSE, algorithm = 4)

### Optimal parameter right now:
set.seed(42)
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5847, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = F, variable.features.n = 3000) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 50, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.25, verbose = FALSE, algorithm = 4)


pdf("output/seurat/UMAP_SCT_UNTREATED72hr.pdf", width=10, height=6)
pdf("output/seurat/UMAP_SCT_UNTREATED72hr_raw.pdf", width=10, height=6)

DimPlot(srat_UNTREATED72hr, label = T, repel = T) + ggtitle("UNTREATED72hr_Unsupervised clustering")
dev.off()
table(srat_UNTREATED72hr@meta.data$seurat_clusters)
## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "SOX2", "SOX9", "POU3F2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP")

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


DefaultAssay(srat_UNTREATED72hr) <- "SCT"
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_endoderm_raw.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_mesoderm_raw.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_ectoderm_raw.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_NODAL.pdf", width=5, height=10)
FeaturePlot(srat_UNTREATED72hr, features = "NODAL", max.cutoff = 3, cols = c("grey", "red"))
dev.off()
# All in dotplot
ecto_meso_endo <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "SOX2", "SOX9", "POU3F2","IRX2", "SOX21","NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1","KIT", "PRDM1", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP")

levels(srat_UNTREATED72hr) <- c("Ectoderm", "Mesoderm_1","Mesoderm_2","Mesoderm_3","Endoderm")# reorder clusters
levels(srat_UNTREATED72hr) <- c("Endoderm", "Mesoderm_3", "Mesoderm_2","Mesoderm_1","Ectoderm")
pdf("output/seurat/DotPlot_SCT_UNTREATED72hr.pdf", width=20, height=3)
DotPlot(srat_UNTREATED72hr, features = ecto_meso_endo, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# All in heatmap
pdf("output/seurat/DoHeatmap_SCT_UNTREATED72hr.pdf", width=20, height=3)
DoHeatmap(srat_UNTREATED72hr, assay = "SCT", slot = "scale.data", features = ecto_meso_endo) + NoLegend()
dev.off()

## QC plot
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_count_feature.pdf", width=10, height=7)
VlnPlot(srat_UNTREATED72hr,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## percent mt and rb
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_mt.pdf", width=10, height=7)
VlnPlot(srat_UNTREATED72hr,features = c("percent.mt", "percent.rb")) & 
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









# Aditional list of genes for troubleshooting:
### sub-mesoderm marker from ChatGPT

paraxial_mesoderm <- c("TBX6", "PAX3", "MESP2")
intermediate_mesoderm <- c("LHX1", "PAX2", "WT1")
lateral_plate_mesoderm <- c("MEOX1", "NKX2.5", "TBX5")
chordamesoderm <- c("BRACHYURY", "FOXA2", "SHH")
hemangioblasts <- c("KDR", "TAL1", "LMO2")

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_paraxial_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = paraxial_mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_intermediate_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = intermediate_mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_lateral_plate_mesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = lateral_plate_mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_chordamesoderm.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = chordamesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_hemangioblasts.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = hemangioblasts, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

### other marker from this paper (Minn2020): https://elifesciences.org/articles/59445
ectoderm <- c("SOX2", "NES", "VIM", "ID3")
epiblast <- c("NANOG", "POU5F1", "DPPA4", "NODAL", "GDF3", "TDGF1")
mesoderm <- c("T", "MIXL1", "EOMES", "KDR", "DLL3", "LHX1", "APLNR", "TBX6", "MESP1", "HAS2", "PDGFRA")
endoderm <- c("SOX17", "PRDM1", "FOXA2", "GATA6")
PGC <- c("CXCR4", "NANOS3", "TFAP2C")
TE <- c("CDX2", "GATA3", "KRT7", "GATA2")
amnion <- c("TBX3", "TFAP2A", "HAND1", "WNT6")

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_ectoderm_Minn2020.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_epiblast_Minn2020.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = epiblast, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_mesoderm_Minn2020.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_endoderm_Minn2020.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_PGC_Minn2020.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = PGC, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_TE_Minn2020.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = TE, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_amnion_Minn2020.pdf", width=20, height=20)
FeaturePlot(srat_UNTREATED72hr, features = amnion, max.cutoff = 3, cols = c("grey", "red"))
dev.off()



# Rename cluster

# RAW: new.cluster.ids <- c("Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Endoderm", "Ectoderm", "Mesoderm_4")

new.cluster.ids <- c("Mesoderm_1", "Mesoderm_2", "Endoderm", "Ectoderm", "Mesoderm_3")
names(new.cluster.ids) <- levels(srat_UNTREATED72hr)
srat_UNTREATED72hr <- RenameIdents(srat_UNTREATED72hr, new.cluster.ids)

pdf("output/seurat/UMAP_SCT_UNTREATED72hr_raw.pdf", width=10, height=6)
pdf("output/seurat/UMAP_SCT_UNTREATED72hr.pdf", width=10, height=6)
DimPlot(srat_UNTREATED72hr, label = T, repel = T, pt.size = 0.7, label.size = 6) + ggtitle("UNTREATED72hr_Unsupervised clustering")
dev.off()



## Display the top 10 marker genes of each cluster; unbiased
### Find all markers 
all_markers <- FindAllMarkers(srat_UNTREATED72hr, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_UNTREATED72hr_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
##
top10 <- all_markers %>% 
  mutate(cluster = factor(cluster, levels = c("Ectoderm", "Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Endoderm"))) %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Ectoderm", "Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Endoderm")))
# View the top 10 markers for each cluster
print(top10)
write.table(top10, file = "output/seurat/srat_UNTREATED72hr_top10_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# Visualize the top 10 markers for each cluster
marker_genes <- unique(top10$gene)
levels(srat_UNTREATED72hr) <- c("Endoderm", "Mesoderm_3", "Mesoderm_2","Mesoderm_1","Ectoderm")
pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_top10.pdf", width=20, height=3)
DotPlot(srat_UNTREATED72hr, features = marker_genes, cols = c("grey", "red")) + RotatedAxis()
dev.off()







## integrate the two datasets using pearson residuals
srat_DASATINIB72hr <- SCTransform(srat_DASATINIB72hr, verbose = FALSE) %>%  # ,vst.flavor = "v2" 
    RunPCA(npcs = 40, verbose = FALSE)
srat.list <- list(srat_UNTREATED72hr = srat_UNTREATED72hr, srat_DASATINIB72hr = srat_DASATINIB72hr)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

humangastruloid.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
humangastruloid.combined.sct <- IntegrateData(anchorset = humangastruloid.anchors, normalization.method = "SCT")



# Perform integrated analysis
humangastruloid.combined.sct <- RunPCA(humangastruloid.combined.sct, verbose = FALSE, npcs = 40)
humangastruloid.combined.sct <- RunUMAP(humangastruloid.combined.sct, reduction = "pca", dims = 1:40, verbose = FALSE)
humangastruloid.combined.sct <- FindNeighbors(humangastruloid.combined.sct, reduction = "pca", k.param = 50, dims = 1:40)
humangastruloid.combined.sct <- FindClusters(humangastruloid.combined.sct, resolution = 0.25, verbose = FALSE, algorithm = 4)
#RAW: humangastruloid.combined.sct <- FindClusters(humangastruloid.combined.sct, resolution = 0.3, verbose = FALSE, algorithm = 4)

humangastruloid.combined.sct$condition <- factor(humangastruloid.combined.sct$condition, levels = c("UNTREATED72hr", "DASATINIB72hr")) # Reorder untreated 1st

pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr.pdf", width=10, height=6)
pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr_raw.pdf", width=10, height=6)

DimPlot(humangastruloid.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()



# Rename cluster
new.cluster.ids <- c("Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Mesoderm_4","Endoderm", "Ectoderm","Mesoderm_5")
# RAW: new.cluster.ids <- c("Mesoderm_1","Mesoderm_2", "Mesoderm_3","Mesoderm_4", "Endoderm", "Ectoderm", "Mesoderm_5", "Mesoderm_6", "Mesoderm_7")
names(new.cluster.ids) <- levels(humangastruloid.combined.sct)
humangastruloid.combined.sct <- RenameIdents(humangastruloid.combined.sct, new.cluster.ids)

humangastruloid.combined.sct$cluster.annot <- Idents(humangastruloid.combined.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr_label.pdf", width=10, height=6)
pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr_label_raw.pdf", width=10, height=6)

DimPlot(humangastruloid.combined.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 5)
dev.off()

#overlapping condition
pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr_label_overlap.pdf", width=6, height=5)
DimPlot(humangastruloid.combined.sct, reduction = "umap", group.by = "condition", pt.size = 0.000001, cols = c("red","blue"))
dev.off()

# count nb of cells in each cluster
untreated_clusters <- table(Idents(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "UNTREATED72hr"])
print(untreated_clusters)
dasatinib_clusters <- table(Idents(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "DASATINIB72hr"])
print(dasatinib_clusters)

# statisctics to assess if distribution of cells among the clusters is significantly different between the two conditions
## Get the counts of cells in each cluster for each condition
untreated_clusters <- table(Idents(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "UNTREATED72hr"])
dasatinib_clusters <- table(Idents(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "DASATINIB72hr"])
## Combine the counts into a matrix
cluster_counts <- rbind(untreated_clusters, dasatinib_clusters)
rownames(cluster_counts) <- c("UNTREATED72hr", "DASATINIB72hr")
## Perform the chi-square test
chisq_test_result <- chisq.test(cluster_counts)
## Print the result for overall difference
print(chisq_test_result)

# Statiscitcs for in-cluster difference (which cluster has significant different nb of cells):
total_untreated <- sum(untreated_clusters)
total_dasatinib <- sum(dasatinib_clusters)

## Initialize a vector to store p-values
p_values <- numeric(length(untreated_clusters))
## Loop over each cluster and perform a chi-square test
for (i in seq_along(untreated_clusters)) {
  # Build contingency table for the i-th cluster
  contingency_table <- matrix(c(untreated_clusters[i], total_untreated - untreated_clusters[i], 
                                dasatinib_clusters[i], total_dasatinib - dasatinib_clusters[i]), 
                              nrow = 2)
  # Perform chi-square test
  chi_square_result <- chisq.test(contingency_table)
  # Store the p-value
  p_values[i] <- chi_square_result$p.value
}
## Perform Bonferroni correction
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
## Print adjusted p-values
names(adjusted_p_values) <- names(untreated_clusters)
print(adjusted_p_values)


## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "NES", "SOX2", "SOX9", "BRN2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "GSC", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP", "SOX17")

DefaultAssay(humangastruloid.combined.sct) <- "SCT"
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_endoderm.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = endoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_mesoderm.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_ectoderm.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

DefaultAssay(humangastruloid.combined.sct) <- "SCT"
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_endoderm_raw.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = endoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_mesoderm_raw.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_ectoderm_raw.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


# All in dotplot
ecto_meso_endo <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "SOX2", "SOX9", "POU3F2","IRX2", "SOX21","NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1","KIT", "PRDM1", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP")

levels(humangastruloid.combined.sct) <- c("Ectoderm", "Mesoderm_1","Mesoderm_2","Mesoderm_3","Endoderm")# reorder clusters
levels(humangastruloid.combined.sct) <- c("Endoderm", "Mesoderm_5","Mesoderm_4","Mesoderm_3", "Mesoderm_2","Mesoderm_1","Ectoderm")
pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr.pdf", width=20, height=5)
DotPlot(humangastruloid.combined.sct, features = ecto_meso_endo, cols = c("blue","red"), split.by = 'condition') + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_ident.pdf", width=14, height=3)
DotPlot(humangastruloid.combined.sct, assay = "SCT", features = ecto_meso_endo, cols = c("grey", "red")) + RotatedAxis()
dev.off()

# All in heatmap
pdf("output/seurat/DoHeatmap_SCT_UNTREATED72hr_DASATINIB72hr.pdf", width=10, height=5)
DoHeatmap(humangastruloid.combined.sct, assay = "integrated", slot = "scale.data", features = ecto_meso_endo, group.by = "ident", group.bar = TRUE)
dev.off()


pluripotency <- c("NANOG", "POU5F1", "DNMT3B", "CDH1", "SOX2") # from this https://www.sciencedirect.com/science/article/pii/S2213671121006512#fig1

pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_pluripotency.pdf", width=10, height=5)
DotPlot(humangastruloid.combined.sct, features = pluripotency, cols = c("blue","red"), split.by = 'condition') + RotatedAxis()
dev.off()


### Check YAP1 bulk-regulated genes:
YAP1_NODAL <- c("CER1", "TDGF1", "FOXH1", "NODAL", "DACT1", "TDGF1P3", "CFC1", "ACVR1C", "CITED2", "DACT2", "ACVR1B", "DMRT1", "CFC1B", "DAND5", "TGIF2", "SMAD3") 

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_NODAL.pdf", width=10, height=80)
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_NODAL_raw.pdf", width=10, height=80)

FeaturePlot(humangastruloid.combined.sct, features = YAP1_NODAL, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_NODAL.pdf", width=10, height=5)
DotPlot(humangastruloid.combined.sct, features = YAP1_NODAL, cols = c("blue","red"), split.by = 'condition') + RotatedAxis()
dev.off()


pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_CTGF.pdf", width=10, height=5)
FeaturePlot(humangastruloid.combined.sct, features = "CCN2", max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

#### From ppt


pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_NODAL.pdf", width=10, height=5)
DotPlot(humangastruloid.combined.sct, features = c("SHH", "DMRT1", "CITED2", 'DACT2', "TGIF2", "ACVR1B", 'SMAD3', 'DAND5', 'SMAD2','CER1', 'NODAL', 'FOXH1', 'DACT1', 'TDGF1', 'TDGF1P3', 'ACVR1C', 'CFC1', 'CFC1B'), cols = c("blue","red"), split.by = 'condition') + RotatedAxis()
dev.off()


pdf("output/seurat/DoHeatmap_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_NODAL.pdf", width=10, height=5)
DoHeatmap(humangastruloid.combined.sct, assay = "integrated", slot = "scale.data", features = c("SHH", "DMRT1", "CITED2", 'DACT2', "TGIF2", "ACVR1B", 'SMAD3', 'DAND5', 'SMAD2','CER1', 'NODAL', 'FOXH1', 'DACT1', 'TDGF1', 'TDGF1P3', 'ACVR1C', 'CFC1', 'CFC1B'), group.by = "condition", group.bar = TRUE)
dev.off()


pdf("output/seurat/RidgePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_NODAL.pdf", width=10, height=25)
RidgePlot(humangastruloid.combined.sct, features = c('CITED2','TGIF2', 'ACVR1B', 'SMAD3', 'SMAD2','FOXH1', 'DACT1', 'CFC1'), cols = c("red","blue"), group.by = 'condition', ncol =1) + RotatedAxis()
dev.off()

## From Estaras 2017  10.1101/gad.307512.117 
### I filter at least 100 count for WT or YAP and pick the most qvalue then FC

pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_Top10RepressInduced.pdf", width=10, height=5)
RidgePlot(humangastruloid.combined.sct, features = c("LINC00458", "AMOTL2", "RAB17", "TAGLN", "TPM1", "ACTA1", "MYL9", "ARHGAP23", "CHCHD2", "YAP1", "TXNIP", "DDIT4", "CHST8", "NODAL", "GRM4", "IGFBP2", "PLBD1", "B2M", "CLIP2", "ERVH48-1"), cols = c("blue","red"), split.by = "condition") + RotatedAxis()
dev.off()






"PAX8-AS1","EBF2"

### other representation to troubleshoot:
humangastruloid.combined.sct.UNTREATED72hr = humangastruloid.combined.sct$condition == "UNTREATED72hr"
pdf("output/seurat/RidgePlot_SCT_UNTREATED72hr_YAP1_NODAL.pdf", width=50, height=100)
RidgePlot(humangastruloid.combined.sct.UNTREATED72hr, ncol = 1, features = YAP1_NODAL)
dev.off()
###

# differential expressed genes across conditions
## what genes change in different conditions for cells of the same type

humangastruloid.combined.sct$celltype.stim <- paste(humangastruloid.combined.sct$cluster.annot, humangastruloid.combined.sct$condition,
    sep = "-")
Idents(humangastruloid.combined.sct) <- "celltype.stim"

# use SCT corrected count for DEGs
humangastruloid.combined.sct <- PrepSCTFindMarkers(humangastruloid.combined.sct)
## DEGs for each cell type: "Endoderm", "Mesoderm_5","Mesoderm_4","Mesoderm_3", "Mesoderm_2","Mesoderm_1","Ectoderm"
Endoderm.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Endoderm-DASATINIB72hr", ident.2 = "Endoderm-UNTREATED72hr",
    verbose = FALSE)
head(Endoderm.DASATINIB.response, n = 15)
write.table(Endoderm.DASATINIB.response, file = "output/seurat/Endoderm-DASATINIB_response.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_1.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_1-DASATINIB72hr", ident.2 = "Mesoderm_1-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_1.DASATINIB.response, n = 15)
write.table(Mesoderm_1.DASATINIB.response, file = "output/seurat/Mesoderm_1-DASATINIB_response.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_2.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_2-DASATINIB72hr", ident.2 = "Mesoderm_2-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_2.DASATINIB.response, n = 15)
write.table(Mesoderm_2.DASATINIB.response, file = "output/seurat/Mesoderm_2-DASATINIB_response.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_3.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_3-DASATINIB72hr", ident.2 = "Mesoderm_3-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_3.DASATINIB.response, n = 15)
write.table(Mesoderm_3.DASATINIB.response, file = "output/seurat/Mesoderm_3-DASATINIB_response.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_4.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_4-DASATINIB72hr", ident.2 = "Mesoderm_4-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_4.DASATINIB.response, n = 15)
write.table(Mesoderm_4.DASATINIB.response, file = "output/seurat/Mesoderm_4-DASATINIB_response.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_5.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_5-DASATINIB72hr", ident.2 = "Mesoderm_5-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_5.DASATINIB.response, n = 15)
write.table(Mesoderm_5.DASATINIB.response, file = "output/seurat/Mesoderm_5-DASATINIB_response.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Ectoderm.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Ectoderm-DASATINIB72hr", ident.2 = "Ectoderm-UNTREATED72hr",
    verbose = FALSE)
head(Ectoderm.DASATINIB.response, n = 15)
write.table(Ectoderm.DASATINIB.response, file = "output/seurat/Ectoderm-DASATINIB_response.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# PLOT norm for DEGs count
## humangastruloid.combined.sct$condition <- factor(humangastruloid.combined.sct$condition, levels = c("UNTREATED72hr", "DASATINIB72hr")) # Reorder untreated 1st

Idents(humangastruloid.combined.sct) <- "cluster.annot"
DefaultAssay(humangastruloid.combined.sct) <- "SCT"
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_norm_YAP1_NODAL.pdf", width=7, height=50)
FeaturePlot(humangastruloid.combined.sct, features = c("SHH", "DMRT1", "CITED2", 'DACT2', "TGIF2", "ACVR1B", 'SMAD3', 'DAND5', 'SMAD2','CER1', 'NODAL', 'FOXH1', 'DACT1', 'TDGF1', 'TDGF1P3', 'ACVR1C', 'CFC1', 'CFC1B'), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_norm_hippoYAP.pdf", width=7, height=25)
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_norm_hippoYAP_raw.pdf", width=7, height=25)

FeaturePlot(humangastruloid.combined.sct, features = c("TEAD1", "TEAD4", "CCN2", "CCN1","AREG","MYC","GLI2","VIM","AXL","BIRC5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/DoHeatmap_SCT_UNTREATED72hr_DASATINIB72hr_norm_YAP1_NODAL.pdf", width=10, height=5)
DoHeatmap(humangastruloid.combined.sct, features = c("SHH", "DMRT1", "CITED2", 'DACT2', "TGIF2", "ACVR1B", 'SMAD3', 'DAND5', 'SMAD2','CER1', 'NODAL', 'FOXH1', 'DACT1', 'TDGF1', 'TDGF1P3', 'ACVR1C', 'CFC1', 'CFC1B'), group.by = "condition", group.bar = TRUE)
dev.off()





## Display the top 10 marker genes of each cluster; this version display 10 to 20 genes; either the top10 are shared between condition
### Find all markers 
all_markers <- FindAllMarkers(humangastruloid.combined.sct, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_UNTREATED72hr_DASATINIB72hr_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
##
top10 <- all_markers %>% 
  mutate(cluster = factor(cluster, levels = c("Ectoderm-UNTREATED72hr","Ectoderm-DASATINIB72hr", "Mesoderm_1-UNTREATED72hr", "Mesoderm_1-DASATINIB72hr", "Mesoderm_2-UNTREATED72hr", "Mesoderm_2-DASATINIB72hr", "Mesoderm_3-UNTREATED72hr", "Mesoderm_3-DASATINIB72hr", "Mesoderm_4-UNTREATED72hr", "Mesoderm_4-DASATINIB72hr", "Mesoderm_5-UNTREATED72hr", "Mesoderm_5-DASATINIB72hr", "Endoderm-UNTREATED72hr", "Endoderm-DASATINIB72hr" ))) %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Ectoderm-UNTREATED72hr","Ectoderm-DASATINIB72hr", "Mesoderm_1-UNTREATED72hr", "Mesoderm_1-DASATINIB72hr", "Mesoderm_2-UNTREATED72hr", "Mesoderm_2-DASATINIB72hr", "Mesoderm_3-UNTREATED72hr", "Mesoderm_3-DASATINIB72hr", "Mesoderm_4-UNTREATED72hr", "Mesoderm_4-DASATINIB72hr", "Mesoderm_5-UNTREATED72hr", "Mesoderm_5-DASATINIB72hr", "Endoderm-UNTREATED72hr", "Endoderm-DASATINIB72hr")))
# View the top 10 markers for each cluster
print(top10)
write.table(top10, file = "output/seurat/srat_UNTREATED72hr_DASATINIB72hr_top10_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# Visualize the top 10 markers for each cluster
marker_genes <- unique(top10$gene)
levels(humangastruloid.combined.sct) <- c("Endoderm-DASATINIB72hr", "Endoderm-UNTREATED72hr", "Mesoderm_5-DASATINIB72hr", "Mesoderm_5-UNTREATED72hr", "Mesoderm_4-DASATINIB72hr", "Mesoderm_4-UNTREATED72hr", "Mesoderm_3-DASATINIB72hr", "Mesoderm_3-UNTREATED72hr", "Mesoderm_2-DASATINIB72hr", "Mesoderm_2-UNTREATED72hr", "Mesoderm_1-DASATINIB72hr", "Mesoderm_1-UNTREATED72hr", "Ectoderm-DASATINIB72hr", "Ectoderm-UNTREATED72hr")
pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_top10.pdf", width=26, height=5)
DotPlot(humangastruloid.combined.sct, features = marker_genes, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# Display the top 10 CONSERVED marker genes of each cluster
## DEGs cluster versus all other
Ectoderm.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Ectoderm", grouping.var = "condition", verbose = FALSE) %>% mutate(cluster = "Ectoderm")
Mesoderm_1.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_1", grouping.var = "condition", verbose = FALSE) %>% mutate(cluster = "Mesoderm_1")
Mesoderm_2.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_2", grouping.var = "condition", verbose = FALSE) %>% mutate(cluster = "Mesoderm_2")
Mesoderm_3.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_3", grouping.var = "condition", verbose = FALSE) %>% mutate(cluster = "Mesoderm_3")
Mesoderm_4.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_4", grouping.var = "condition", verbose = FALSE) %>% mutate(cluster = "Mesoderm_4")
Mesoderm_5.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Mesoderm_5", grouping.var = "condition", verbose = FALSE) %>% mutate(cluster = "Mesoderm_5")
Endoderm.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "SCT", ident.1 = "Endoderm", grouping.var = "condition", verbose = FALSE) %>% mutate(cluster = "Endoderm")
## Combine all conserved markers into one data frame
all_conserved <- bind_rows(Ectoderm.conserved, Mesoderm_1.conserved, Mesoderm_2.conserved, Mesoderm_3.conserved, Mesoderm_4.conserved, Mesoderm_5.conserved, Endoderm.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Ectoderm", "Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Mesoderm_5", "Endoderm"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Ectoderm", "Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Mesoderm_5", "Endoderm")))


## Write the top 10 conserved markers for each cluster to a file
write.table(top10_conserved, file = "output/seurat/srat_top10_conserved_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Visualize the top 10 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(humangastruloid.combined.sct) <- c("Endoderm", "Mesoderm_5", "Mesoderm_4", "Mesoderm_3", "Mesoderm_2", "Mesoderm_1", "Ectoderm")
pdf("output/seurat/DotPlot_SCT_top10_conserved.pdf", width=26, height=4)
DotPlot(humangastruloid.combined.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()





# Past troubleshoot for control data to check which normalization method was the best (in the end SCT V1 is the best)
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

```


### Troubleshooting for **metap** package installation (required for `FindConservedMarkers()`)
```R
# Install metap package
install.packages('metap')
fail with dependencies for 'qqconf' and 'mutoss'
install.packages('qqconf') # FAIL AGAIN
```
--> Need to install through conda the fft3 C libraries

```bash
conda install -c eumetsat fftw3
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/home/roulet/anaconda3/envs/scRNAseq/lib/pkgconfig/ # help conda to find the fftw3
conda install -c bioconda r-mutoss
```
--> Now metap can be loaded; but is messed up all the other package because of Matrix1.4.0 is load and need >1.5; so I remove Matrix with `conda uninstall r-matrix=1.4.1` and re-install lot of R packages in the following order: `install.packages('Matrix')` and `install.packages('MASS')` and `install.packages('codetools')` and `install.packages('survival')` and `BiocManager::install("BiocGenerics")` and `BiocManager::install("Biobase")` and `install.packages('mvtnorm')` and `BiocManager::install("multtest")` and `install.packages('mutoss')`

--> **NOTE: the export may need to be re-run at each session I want to use metap**



#### Main conclusion for human gastruloid_analysis V1
- the optimal parameter seems to be SCT Transform version1 (not vst.flavor V2) with 40 dim and 0.25 resolutions. The
- More cells belong to mesoderm; identification of 5 sub-clusters of mesoderm
- Endoderm is the more clearly sperated cluster
- Dasatanib treatment efficacy is not clear; target genes from Hippo-YAP and NODAL are not or very lowly affected

Some issues:
- Different parameters (nb of dimensions, and regression) used for SCTransform in Untreated and Dasatinib and integrated together --> Better to use same parameters for both
- Downsampling is shit (lose information randomly!); but we may give a try (see [here](https://www.biostars.org/p/423349/))
- For DEG (FindMarkers) use the RNA assay and not the SCT assay!

In the V2 version I will do the integration analysis more accurately:

# human gastruloid_analysis V2_integration only

Workflow (from [here](https://github.com/satijalab/seurat/issues/1836); and [here](https://github.com/satijalab/seurat/discussions/4032) ):
- QC filtering
- default assay to "RNA" and normalize and scale data with `NormalizeData` and `ScaleData` (normalised counts to calculate the cell cycle scores with CellCycleScoring that later I regress out in SCTransform and secondly, I want to have my normalised raw counts for DE expression later on when I have integrated the dataset)
- Run SCT normalization on each dataset (with same parameter)
- dataset integration
- PCA/UMAP/FindClusters/FindNeighbors (on the default "integrated" assay)
- Change again default assay to "RNA" and run again `NormalizeData` and `ScaleData` so that integrated data are normalized; and generate FeaturePlots and perform DEGs (**for DEG use assay RNA and slot data**)
(to make **featurePlot/ViolinPlot use RNA assay NOT SURE SLOT?**; and use scale.data from SCT assay for DoHeatmap )


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




## After seeing the plot; add QC information in our seurat object
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 2000 & srat_UNTREATED72hr@meta.data$QC == 'Pass','Low_nFeature',srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 2000 & srat_UNTREATED72hr@meta.data$QC != 'Pass' & srat_UNTREATED72hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_UNTREATED72hr@meta.data$QC,sep = ','),srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$percent.mt > 15 & srat_UNTREATED72hr@meta.data$QC == 'Pass','High_MT',srat_UNTREATED72hr@meta.data$QC)
srat_UNTREATED72hr[['QC']] <- ifelse(srat_UNTREATED72hr@meta.data$nFeature_RNA < 2000 & srat_UNTREATED72hr@meta.data$QC != 'Pass' & srat_UNTREATED72hr@meta.data$QC != 'High_MT',paste('High_MT',srat_UNTREATED72hr@meta.data$QC,sep = ','),srat_UNTREATED72hr@meta.data$QC)
table(srat_UNTREATED72hr[['QC']])
## 



srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$nFeature_RNA < 2000 & srat_DASATINIB72hr@meta.data$QC == 'Pass','Low_nFeature',srat_DASATINIB72hr@meta.data$QC)
srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$nFeature_RNA < 2000 & srat_DASATINIB72hr@meta.data$QC != 'Pass' & srat_DASATINIB72hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_DASATINIB72hr@meta.data$QC,sep = ','),srat_DASATINIB72hr@meta.data$QC)
srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$percent.mt > 15 & srat_DASATINIB72hr@meta.data$QC == 'Pass','High_MT',srat_DASATINIB72hr@meta.data$QC)
srat_DASATINIB72hr[['QC']] <- ifelse(srat_DASATINIB72hr@meta.data$nFeature_RNA < 2000 & srat_DASATINIB72hr@meta.data$QC != 'Pass' & srat_DASATINIB72hr@meta.data$QC != 'High_MT',paste('High_MT',srat_DASATINIB72hr@meta.data$QC,sep = ','),srat_DASATINIB72hr@meta.data$QC)
table(srat_DASATINIB72hr[['QC']])




## subset my seurat object to only analyze the cells that pass the QC
srat_UNTREATED72hr <- subset(srat_UNTREATED72hr, subset = QC == 'Pass')
srat_DASATINIB72hr <- subset(srat_DASATINIB72hr, subset = QC == 'Pass')
srat_UNTREATED72hr$condition <- "UNTREATED72hr"
srat_DASATINIB72hr$condition <- "DASATINIB72hr"

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
srat_UNTREATED72hr <- NormalizeData(srat_UNTREATED72hr, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_UNTREATED72hr)
srat_UNTREATED72hr <- ScaleData(srat_UNTREATED72hr, features = all.genes) # zero-centres and scales it

srat_DASATINIB72hr <- NormalizeData(srat_DASATINIB72hr, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_DASATINIB72hr)
srat_DASATINIB72hr <- ScaleData(srat_DASATINIB72hr, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
srat_UNTREATED72hr <- CellCycleScoring(srat_UNTREATED72hr, s.features = s.genes, g2m.features = g2m.genes)
table(srat_UNTREATED72hr[[]]$Phase)
srat_DASATINIB72hr <- CellCycleScoring(srat_DASATINIB72hr, s.features = s.genes, g2m.features = g2m.genes)
table(srat_DASATINIB72hr[[]]$Phase)

set.seed(42)

# Run SCTransform
## Version OK with 2000 treshold RNA
srat_UNTREATED72hr <- SCTransform(srat_UNTREATED72hr, method = "glmGamPoi", ncells = 5713, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 25, verbose = FALSE)

srat_DASATINIB72hr <- SCTransform(srat_DASATINIB72hr, method = "glmGamPoi", ncells = 7148, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 25, verbose = FALSE)


# Data integration (check active assay is 'SCT')
srat.list <- list(srat_UNTREATED72hr = srat_UNTREATED72hr, srat_DASATINIB72hr = srat_DASATINIB72hr)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

humangastruloid.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
humangastruloid.combined.sct <- IntegrateData(anchorset = humangastruloid.anchors, normalization.method = "SCT")



# Perform integrated analysis (check active assay is 'integrated')
## XXX RE FINE A BIT BUT ALMOST PERFECT !!! tested 0.6 not good; play with k.param!! was 40 dim at sctransform XXX

DefaultAssay(humangastruloid.combined.sct) <- "integrated"

humangastruloid.combined.sct <- RunPCA(humangastruloid.combined.sct, verbose = FALSE, npcs = 25)
humangastruloid.combined.sct <- RunUMAP(humangastruloid.combined.sct, reduction = "pca", dims = 1:25, verbose = FALSE)
humangastruloid.combined.sct <- FindNeighbors(humangastruloid.combined.sct, reduction = "pca", k.param = 15, dims = 1:25)
humangastruloid.combined.sct <- FindClusters(humangastruloid.combined.sct, resolution = 0.2, verbose = FALSE, algorithm = 4)
#RAW: humangastruloid.combined.sct <- FindClusters(humangastruloid.combined.sct, resolution = 0.3, verbose = FALSE, algorithm = 4)

humangastruloid.combined.sct$condition <- factor(humangastruloid.combined.sct$condition, levels = c("UNTREATED72hr", "DASATINIB72hr")) # Reorder untreated 1st

pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr_V2.pdf", width=10, height=6)
pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr_V2_test.pdf", width=10, height=6)
DimPlot(humangastruloid.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()



## Check some marker genes 
### Marker gene list
ectoderm <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "NES", "SOX2", "SOX9", "BRN2","IRX2", "SOX21")
mesoderm <- c("NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1")
endoderm <- c("KIT", "PRDM1", "GSC", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP", "SOX17")

DefaultAssay(humangastruloid.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_endoderm_V2.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = endoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_mesoderm_V2.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_ectoderm_V2.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


DefaultAssay(humangastruloid.combined.sct) <- "RNA" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_endoderm_V2_RNA.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = endoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_mesoderm_V2_RNA.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_ectoderm_V2_RNA.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


## Unbiased markers:
DefaultAssay(humangastruloid.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

unbiased = c('GABRB3', 'WNT8A', 'LIX1', 'PLAT', 'POLQ', 'CDC20', 'APOA2')
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_unbiased_V2.pdf", width=10, height=40)
FeaturePlot(humangastruloid.combined.sct, features = unbiased, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


## QC plot
## percent mt and rb
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_DASATINIB72hr_mt_rb_V2.pdf", width=10, height=7)
VlnPlot(humangastruloid.combined.sct,features = c("percent.mt", "percent.rb"), split.by = "condition") & 
  theme(plot.title = element_text(size=10))
dev.off()
## RNA
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_DASATINIB72hr_count_feature_V2.pdf", width=10, height=10)
VlnPlot(humangastruloid.combined.sct,features = c("nCount_RNA","nFeature_RNA"), split.by = "condition") & 
  theme(plot.title = element_text(size=10))
dev.off()
## cell cycle
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_DASATINIB72hr_cellCycle_V2.pdf", width=15, height=7)
VlnPlot(humangastruloid.combined.sct,features = c("S.Score","G2M.Score"), split.by = "condition") & 
  theme(plot.title = element_text(size=10))
dev.off()  
pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_cellCycle_V2.pdf", width=10, height=10)
FeaturePlot(humangastruloid.combined.sct,features = c("S.Score","G2M.Score"), split.by = "condition") & 
  theme(plot.title = element_text(size=10))
dev.off()  

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_Phase_V2.pdf", width=10, height=10)
DimPlot(humangastruloid.combined.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  

# Rename cluster
new.cluster.ids <- c("Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Endoderm","Ectoderm","Mesoderm_4")
names(new.cluster.ids) <- levels(humangastruloid.combined.sct)
humangastruloid.combined.sct <- RenameIdents(humangastruloid.combined.sct, new.cluster.ids)

humangastruloid.combined.sct$cluster.annot <- Idents(humangastruloid.combined.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr_label_V2.pdf", width=10, height=6)

DimPlot(humangastruloid.combined.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 5)
dev.off()

#overlapping condition
pdf("output/seurat/UMAP_UNTREATED72hr_DASATINIB72hr_label_overlap_V2.pdf", width=6, height=5)
DimPlot(humangastruloid.combined.sct, reduction = "umap", group.by = "condition", pt.size = 0.000001, cols = c("blue","red"))
dev.off()


# All in dotplot
DefaultAssay(humangastruloid.combined.sct) <- "SCT"

ecto_meso_endo <- c("PAX6", "OTX2", "NES", "TFAP2A", "DLK1", "LHX2", "RAX", "HES5", "SOX1", "SOX2", "SOX9", "POU3F2","IRX2", "SOX21","NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR", "MIXL1", "CXCR4", "ETV3", "TBX6", "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "NODAL", "TDGF1", "GDF3", "PAX2", "DLL1", "MEOX1", "MEIS2", "Vg1","KIT", "PRDM1", "TTR", "SOX17", "GSC", "CER1", "IRX1", "EOMES", "SHH", "FOXA2", "Cdx1", "GATA4", "Cdx4", "GATA6", "HHEX", "LEFTY1", "AFP")

levels(humangastruloid.combined.sct) <- c("Endoderm", "Mesoderm_4","Mesoderm_3", "Mesoderm_2","Mesoderm_1","Ectoderm")
pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_ident_V2.pdf", width=14, height=3)
DotPlot(humangastruloid.combined.sct, assay = "SCT", features = ecto_meso_endo, cols = c("grey", "red")) + RotatedAxis()
dev.off()


ecto_meso_endo_subset <- c("NES", "TFAP2A", "SOX2", "SOX9", "IRX2", "NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR",  "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "MEIS2", "Vg1","KIT", "PRDM1", "TTR", "SOX17", "IRX1", "FOXA2", "Cdx1", "GATA4", "GATA6", "HHEX")

ecto_meso_endo_subset <- c("NES", "TFAP2A", "SOX2", "SOX9", "IRX2", "NCAM1", "HAND1", "PDGFRA", "CDX2", "APLNR",  "LHX1", "KDR", "ID2","APLINR", "T", "MESP1", "HAS2", "SNAJ1", "SNAJ2", "MEIS2", "Vg1","KIT", "PRDM1", "TTR", "SOX17", "IRX1", "FOXA2", "Cdx1", "GATA4", "GATA6", "HHEX")

levels(humangastruloid.combined.sct) <- c("Endoderm", "Mesoderm_4","Mesoderm_3", "Mesoderm_2","Mesoderm_1","Ectoderm")
pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_ident_subset_V2.pdf", width=14, height=3)
DotPlot(humangastruloid.combined.sct, assay = "SCT", features = ecto_meso_endo_subset, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# count nb of cells in each cluster
untreated_clusters <- table(Idents(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "UNTREATED72hr"])
print(untreated_clusters)
dasatinib_clusters <- table(Idents(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "DASATINIB72hr"])
print(dasatinib_clusters)

# statisctics to assess if distribution of cells among the clusters is significantly different between the two conditions
## Get the counts of cells in each cluster for each condition
untreated_clusters <- table(Idents(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "UNTREATED72hr"])
dasatinib_clusters <- table(Idents(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "DASATINIB72hr"])
## Combine the counts into a matrix
cluster_counts <- rbind(untreated_clusters, dasatinib_clusters)
rownames(cluster_counts) <- c("UNTREATED72hr", "DASATINIB72hr")
## Perform the chi-square test
chisq_test_result <- chisq.test(cluster_counts)
## Print the result for overall difference
print(chisq_test_result)

# Statiscitcs for in-cluster difference (which cluster has significant different nb of cells):
total_untreated <- sum(untreated_clusters)
total_dasatinib <- sum(dasatinib_clusters)

## Initialize a vector to store p-values
p_values <- numeric(length(untreated_clusters))
## Loop over each cluster and perform a chi-square test
for (i in seq_along(untreated_clusters)) {
  # Build contingency table for the i-th cluster
  contingency_table <- matrix(c(untreated_clusters[i], total_untreated - untreated_clusters[i], 
                                dasatinib_clusters[i], total_dasatinib - dasatinib_clusters[i]), 
                              nrow = 2)
  # Perform chi-square test
  chi_square_result <- chisq.test(contingency_table)
  # Store the p-value
  p_values[i] <- chi_square_result$p.value
}
## Perform Bonferroni correction
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
## Print adjusted p-values
names(adjusted_p_values) <- names(untreated_clusters)
print(adjusted_p_values)


## Downsampling with bootstrap to compare the nb of cell per cell types_V1
# --> Version works great! But not beautiful
### Identify the unique clusters
unique_clusters <- unique(Idents(humangastruloid.combined.sct))

### Create empty matrices to store cell counts
untreated_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
dasatinib_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(untreated_clusters_counts) <- unique_clusters
colnames(dasatinib_clusters_counts) <- unique_clusters


######## THIS CODE IS NOT GOOD, check the one after  ########
### Loop through 100 iterations
for (i in 1:100) {
  # Downsampling
  humangastruloid.combined.sct_DASATINIB72hr_downsample <- sample(humangastruloid.combined.sct_DASATINIB72hr, 5713)
  humangastruloid.combined.sct_integrated_downsample <- humangastruloid.combined.sct[,c(humangastruloid.combined.sct_UNTREATED72hr, humangastruloid.combined.sct_DASATINIB72hr_downsample)]

  # Count nb of cells in each cluster
  untreated_clusters <- table(Idents(humangastruloid.combined.sct_integrated_downsample)[humangastruloid.combined.sct_integrated_downsample$condition == "UNTREATED72hr"])
  dasatinib_clusters <- table(Idents(humangastruloid.combined.sct_integrated_downsample)[humangastruloid.combined.sct_integrated_downsample$condition == "DASATINIB72hr"])

  # Align the counts with the unique clusters
  untreated_clusters_counts[i, names(untreated_clusters)] <- as.numeric(untreated_clusters)
  dasatinib_clusters_counts[i, names(dasatinib_clusters)] <- as.numeric(dasatinib_clusters)
}
### Calculate mean counts
mean_untreated_clusters <- colMeans(untreated_clusters_counts)
mean_dasatinib_clusters <- colMeans(dasatinib_clusters_counts)
### Combine the counts into a matrix
cluster_counts <- rbind(mean_untreated_clusters, mean_dasatinib_clusters)
rownames(cluster_counts) <- c("UNTREATED72hr", "DASATINIB72hr")
#### You can then proceed with your chi-square tests and other statistics as before
#### To visualize the mean number of cells per cell type
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling10_V2.pdf", width=10, height=6)
barplot(as.matrix(cluster_counts), beside=TRUE, legend.text=rownames(cluster_counts))
dev.off()
################

### V2
# I think correct code to subset:
humangastruloid.combined.sct_UNTREATED72hr <- which(humangastruloid.combined.sct$orig.ident == 'UNTREATED72hr')
humangastruloid.combined.sct_DASATINIB72hr <- which(humangastruloid.combined.sct$orig.ident == 'DASATINIB72hr')
#
library(tidyverse)

### Identify the unique clusters
unique_clusters <- unique(Idents(humangastruloid.combined.sct))

### Create empty matrices to store cell counts
untreated_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
dasatinib_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(untreated_clusters_counts) <- unique_clusters
colnames(dasatinib_clusters_counts) <- unique_clusters

### Loop through 100 iterations
for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  humangastruloid.combined.sct_DASATINIB72hr_downsample <- sample(humangastruloid.combined.sct_DASATINIB72hr, 5713)
  humangastruloid.combined.sct_integrated_downsample <- humangastruloid.combined.sct[,c(humangastruloid.combined.sct_UNTREATED72hr, humangastruloid.combined.sct_DASATINIB72hr_downsample)]

  # Count nb of cells in each cluster
  untreated_clusters <- table(Idents(humangastruloid.combined.sct_integrated_downsample)[humangastruloid.combined.sct_integrated_downsample$condition == "UNTREATED72hr"])
  dasatinib_clusters <- table(Idents(humangastruloid.combined.sct_integrated_downsample)[humangastruloid.combined.sct_integrated_downsample$condition == "DASATINIB72hr"])

  # Align the counts with the unique clusters
  untreated_clusters_counts[i, names(untreated_clusters)] <- as.numeric(untreated_clusters)
  dasatinib_clusters_counts[i, names(dasatinib_clusters)] <- as.numeric(dasatinib_clusters)
}

### Calculate mean and standard error
mean_untreated_clusters <- colMeans(untreated_clusters_counts)
mean_dasatinib_clusters <- colMeans(dasatinib_clusters_counts)
std_error_dasatinib_clusters <- apply(dasatinib_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("Untreated", "Dasatinib")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- untreated_clusters_counts[j,i]
    contingency_table[1,2] <- dasatinib_clusters_counts[j,i]
    contingency_table[2,1] <- sum(untreated_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(dasatinib_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_untreated_clusters),
  untreated = mean_untreated_clusters,
  dasatinib = mean_dasatinib_clusters,
  std_error_dasatinib = std_error_dasatinib_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_dasatinib, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "UNTREATED72hr", "DASATINIB72hr"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.01, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("UNTREATED72hr", "DASATINIB72hr")) # Reorder untreated 1st

# Plotting using ggplot2
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling10_clean_V2.pdf", width=6, height=4)
ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = filter(plot_data, condition == "DASATINIB72hr"),
    aes(label = significance, y = value + std_error_dasatinib),
    vjust = -0.8,
    position = position_dodge(0.9), size = 8
  ) +
  scale_fill_manual(values = c("UNTREATED72hr" = "blue", "DASATINIB72hr" = "red")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()









### Check YAP1 NODAL bulk-regulated genes:
DefaultAssay(humangastruloid.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

YAP1_NODAL <- c("SHH", "DMRT1", "CITED2", 'DACT2', "TGIF2", "ACVR1B", 'SMAD3', 'DAND5', 'SMAD2','CER1', 'NODAL', 'FOXH1', 'DACT1', 'TDGF1', 'TDGF3', 'ACVR1C', 'CFC1', 'CFC1B', 'WNT3A' , 'YAP1') 

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_NODAL_V2.pdf", width=10, height=100)
FeaturePlot(humangastruloid.combined.sct, features = YAP1_NODAL, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_test.pdf", width=10, height=60)
FeaturePlot(humangastruloid.combined.sct, features = c("MYL4", "S100A10", "S100A11", "SDK1", "NPC2", "KRT19", "MT-CO1"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()




### Check YAP1 hippo bulk-regulated genes:
DefaultAssay(humangastruloid.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_hippo_V2.pdf", width=10, height=50)

FeaturePlot(humangastruloid.combined.sct, features = c("TEAD1", "TEAD4", "CCN2", "CCN1","AREG","MYC","GLI2","VIM","AXL","BIRC5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()



### Check ANXA family from bulk-regulated genes (from after 0801 meeting)



pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_ANXA_V2.pdf", width=10, height=35)
FeaturePlot(humangastruloid.combined.sct, features = c("ANXA1","AMOTL2","CCN1","CCN2","TAGLN","ANXA3","CHCHD2"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


### Check Gata6 and CDX2 (from after 0801 meeting)
DefaultAssay(humangastruloid.combined.sct) <- "SCT"

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_GATA6_CDX2_V2.pdf", width=10, height=10)
FeaturePlot(humangastruloid.combined.sct, features = c("GATA6","CDX2"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

DefaultAssay(humangastruloid.combined.sct) <- "RNA"
humangastruloid.combined.sct$condition <- factor(humangastruloid.combined.sct$condition, levels = c("UNTREATED72hr", "DASATINIB72hr"))
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_DASATINIB72hr_GATA6_CDX2_V2.pdf", width=10, height=7)
VlnPlot(humangastruloid.combined.sct, features = c("GATA6","CDX2"), split.by = "condition", log = TRUE) + theme(legend.position = 'top')
dev.off()

##### overlay both genes
DefaultAssay(humangastruloid.combined.sct) <- "SCT"
humangastruloid.untreated <- subset(humangastruloid.combined.sct, subset = condition == "UNTREATED72hr")

# Define a color scale for each gene
transparent_red <- scales::alpha("red", 0.15)
transparent_blue <- scales::alpha("blue", 0.15)
# Generate plots without displaying
p1 <- FeaturePlot(humangastruloid.untreated , features = "GATA6", max.cutoff = 3, cols = c("low" = "transparent", "high" = transparent_red), pt.size = 1.5, raster = FALSE)
p2 <- FeaturePlot(humangastruloid.untreated , features = "CDX2", max.cutoff = 3, cols = c("low" = "transparent", "high" = transparent_blue), pt.size = 1.5, raster = FALSE)
combined_plot <- CombinePlots(plots = list(p1, p2))


pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_GATA6_CDX2_overlay_V2.pdf", width=10, height=5)
print(combined_plot)
dev.off()


humangastruloid.dasatinib <- subset(humangastruloid.combined.sct, subset = condition == "DASATINIB72hr")

# Define a color scale for each gene
transparent_red <- scales::alpha("red", 0.15)
transparent_blue <- scales::alpha("blue", 0.15)
# Generate plots without displaying
p1 <- FeaturePlot(humangastruloid.dasatinib , features = "GATA6", max.cutoff = 3, cols = c("low" = "transparent", "high" = transparent_red), pt.size = 1.5, raster = FALSE)
p2 <- FeaturePlot(humangastruloid.dasatinib , features = "CDX2", max.cutoff = 3, cols = c("low" = "transparent", "high" = transparent_blue), pt.size = 1.5, raster = FALSE)
combined_plot <- CombinePlots(plots = list(p1, p2))


pdf("output/seurat/FeaturePlot_SCT_DASATINIB72hr_GATA6_CDX2_overlay_V2.pdf", width=10, height=5)
print(combined_plot)
dev.off()




### Check EZH2 and JMJD3/KDM6B (from after 0801 meeting)
DefaultAssay(humangastruloid.combined.sct) <- "SCT"

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_EZH2_KDM6B_V2.pdf", width=10, height=10)
FeaturePlot(humangastruloid.combined.sct, features = c("KDM6B","EZH2"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


##### overlay both genes
DefaultAssay(humangastruloid.combined.sct) <- "SCT"
humangastruloid.untreated <- subset(humangastruloid.combined.sct, subset = condition == "UNTREATED72hr")

# Define a color scale for each gene
transparent_red <- scales::alpha("red", 0.15)
transparent_blue <- scales::alpha("blue", 0.15)
# Generate plots without displaying
p1 <- FeaturePlot(humangastruloid.untreated , features = "KDM6B", max.cutoff = 3, cols = c("low" = "transparent", "high" = transparent_red), pt.size = 1.5, raster = FALSE)
p2 <- FeaturePlot(humangastruloid.untreated , features = "EZH2", max.cutoff = 3, cols = c("low" = "transparent", "high" = transparent_blue), pt.size = 1.5, raster = FALSE)
combined_plot <- CombinePlots(plots = list(p1, p2))


pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_EZH2_KDM6B_overlay_V2.pdf", width=10, height=5)
print(combined_plot)
dev.off()


humangastruloid.dasatinib <- subset(humangastruloid.combined.sct, subset = condition == "DASATINIB72hr")

# Define a color scale for each gene
transparent_red <- scales::alpha("red", 0.15)
transparent_blue <- scales::alpha("blue", 0.15)
# Generate plots without displaying
p1 <- FeaturePlot(humangastruloid.dasatinib , features = "KDM6B", max.cutoff = 3, cols = c("low" = "transparent", "high" = transparent_red), pt.size = 1.5, raster = FALSE)
p2 <- FeaturePlot(humangastruloid.dasatinib , features = "EZH2", max.cutoff = 3, cols = c("low" = "transparent", "high" = transparent_blue), pt.size = 1.5, raster = FALSE)
combined_plot <- CombinePlots(plots = list(p1, p2))


pdf("output/seurat/FeaturePlot_SCT_DASATINIB72hr_EZH2_KDM6B_overlay_V2.pdf", width=10, height=5)
print(combined_plot)
dev.off()











### Check top 10 downredulated genes in YAP1KO bulk-regulated genes (table shared after 0801 meeting):
## intitial list top 1-10: c('ZNF558','LINC02693','NXPH2','C6orf141','ZNF736','GCNT4','HHLA1','GABRA1','TEK','ZNF502')
## top10-20: c('LHX8','SCG2','SLFN12','LHFPL6','TNNT2','LINC00458','TXNRD2','CXCL14','AMOTL2','TMC1')
## top20-30: c('RGS5','MYLPF','VEPH1','MYL7','SLFN13','UPK3B','COL12A1','RIPOR2','IRX2','NR4A2')

## only the 10 top detected in our assay: c('LINC02693','NXPH2','TEK','LHX8','LHFPL6','TNNT2','LINC00458','RGS5','VEPH1','MYL7')


pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_top10downBulk_V2.pdf", width=10, height=45)
FeaturePlot(humangastruloid.combined.sct, features = c('LINC02693','NXPH2','TEK','LHX8','LHFPL6','TNNT2','LINC00458','RGS5','VEPH1','MYL7'), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()




### Check top 10 upredulated genes in YAP1KO bulk-regulated genes (table shared after 0801 meeting):
## intitial list top 1-10: c('LOC100134868','EIF1AY','DPYS','ZNF662','FAM228A','PAX8-AS1','RFPL2','LOXHD1','RBM46','EBF2')

## list top 10-20: c('ZNF572','ZNF385D-AS1','SCNN1B','C22orf42','LINC00654','SLC6A20','BMP3','SHISA3','EDNRB','STPG3-AS1')
## list top 20-30: c('STPG3','CTSF','POU3F4','ZXDA','MIR924HG','CCL28','ZNF726','NEURL1','CD1D','ISLR2')
## list top 30-40: c('VMO1','SKOR2','MAEL','CMTM5','FOXC1','PRAM1','GAD2','MADCAM1','HHEX','NPNT')
## list top 40-50: c('SPI1','OLIG1','LRRC37A11P','TDRD12','KCNJ2','LTB','SLC9A3','ZNF541','PLA2G2A','GSTM5')
## list top 50-60: c('LOC100131347','CR1L','HCG4B','PRR15L','NDUFA4L2','EPHB3','NKAPL','DUSP15','UTF1','RELN')
## list top 60-70: c('C1QL2','ARHGEF15','PLXDC1','SP5','FLRT1','GUCA1A','GPT','ADAMTSL1','CRLF1','TFF3')
## list top 70-80: c('MGARP','LGR5','BDKRB2','IL23A','FAAH','NTNG2','SLC7A4','TRIM4','TXNIP','TMIGD2')



## only the 10 top detected in our assay: c('PAX8-AS1','EBF2','EDNRB','MIR924HG','HHEX','NPNT','EPHB3',RELN,SP5,LGR5)


pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_top10upBulk_V2.pdf", width=10, height=45)
FeaturePlot(humangastruloid.combined.sct, features = c('PAX8-AS1','EBF2','EDNRB','MIR924HG','HHEX','NPNT','EPHB3','RELN','SP5','LGR5'), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()







# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

DefaultAssay(humangastruloid.combined.sct) <- "RNA"

humangastruloid.combined.sct <- NormalizeData(humangastruloid.combined.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(humangastruloid.combined.sct)
humangastruloid.combined.sct <- ScaleData(humangastruloid.combined.sct, features = all.genes) # zero-centres and scales it


## what genes change in different conditions for cells of the same type

humangastruloid.combined.sct$celltype.stim <- paste(humangastruloid.combined.sct$cluster.annot, humangastruloid.combined.sct$condition,
    sep = "-")
Idents(humangastruloid.combined.sct) <- "celltype.stim"

# use RNA corrected count for DEGs
## humangastruloid.combined.sct <- PrepSCTFindMarkers(humangastruloid.combined.sct)

## DEGs for each cell type: "Endoderm", "Mesoderm_5","Mesoderm_4","Mesoderm_3", "Mesoderm_2","Mesoderm_1","Ectoderm"
Endoderm.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Endoderm-DASATINIB72hr", ident.2 = "Endoderm-UNTREATED72hr",
    verbose = FALSE)
head(Endoderm.DASATINIB.response, n = 15)
write.table(Endoderm.DASATINIB.response, file = "output/seurat/Endoderm-DASATINIB_response_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_1.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Mesoderm_1-DASATINIB72hr", ident.2 = "Mesoderm_1-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_1.DASATINIB.response, n = 15)
write.table(Mesoderm_1.DASATINIB.response, file = "output/seurat/Mesoderm_1-DASATINIB_response_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_2.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Mesoderm_2-DASATINIB72hr", ident.2 = "Mesoderm_2-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_2.DASATINIB.response, n = 15)
write.table(Mesoderm_2.DASATINIB.response, file = "output/seurat/Mesoderm_2-DASATINIB_response_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_3.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Mesoderm_3-DASATINIB72hr", ident.2 = "Mesoderm_3-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_3.DASATINIB.response, n = 15)
write.table(Mesoderm_3.DASATINIB.response, file = "output/seurat/Mesoderm_3-DASATINIB_response_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Mesoderm_4.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Mesoderm_4-DASATINIB72hr", ident.2 = "Mesoderm_4-UNTREATED72hr",
    verbose = FALSE)
head(Mesoderm_4.DASATINIB.response, n = 15)
write.table(Mesoderm_4.DASATINIB.response, file = "output/seurat/Mesoderm_4-DASATINIB_response_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)

Ectoderm.DASATINIB.response <- FindMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Ectoderm-DASATINIB72hr", ident.2 = "Ectoderm-UNTREATED72hr",
    verbose = FALSE)
head(Ectoderm.DASATINIB.response, n = 15)
write.table(Ectoderm.DASATINIB.response, file = "output/seurat/Ectoderm-DASATINIB_response_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)






## Display the top 10 marker genes of each cluster; this version display 10 to 20 genes; either the top10 are shared between condition
### Find all markers 
all_markers <- FindAllMarkers(humangastruloid.combined.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_UNTREATED72hr_DASATINIB72hr_all_markers_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)
##
top10 <- all_markers %>% 
  mutate(cluster = factor(cluster, levels = c("Ectoderm-UNTREATED72hr","Ectoderm-DASATINIB72hr", "Mesoderm_1-UNTREATED72hr", "Mesoderm_1-DASATINIB72hr", "Mesoderm_2-UNTREATED72hr", "Mesoderm_2-DASATINIB72hr", "Mesoderm_3-UNTREATED72hr", "Mesoderm_3-DASATINIB72hr", "Mesoderm_4-UNTREATED72hr", "Mesoderm_4-DASATINIB72hr", "Mesoderm_5-UNTREATED72hr", "Mesoderm_5-DASATINIB72hr", "Endoderm-UNTREATED72hr", "Endoderm-DASATINIB72hr" ))) %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Ectoderm-UNTREATED72hr","Ectoderm-DASATINIB72hr", "Mesoderm_1-UNTREATED72hr", "Mesoderm_1-DASATINIB72hr", "Mesoderm_2-UNTREATED72hr", "Mesoderm_2-DASATINIB72hr", "Mesoderm_3-UNTREATED72hr", "Mesoderm_3-DASATINIB72hr", "Mesoderm_4-UNTREATED72hr", "Mesoderm_4-DASATINIB72hr", "Mesoderm_5-UNTREATED72hr", "Mesoderm_5-DASATINIB72hr", "Endoderm-UNTREATED72hr", "Endoderm-DASATINIB72hr")))
# View the top 10 markers for each cluster
print(top10)
write.table(top10, file = "output/seurat/srat_UNTREATED72hr_DASATINIB72hr_top10_markers_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# Visualize the top 10 markers for each cluster
marker_genes <- unique(top10$gene)
levels(humangastruloid.combined.sct) <- c("Endoderm-DASATINIB72hr", "Endoderm-UNTREATED72hr", "Mesoderm_5-DASATINIB72hr", "Mesoderm_5-UNTREATED72hr", "Mesoderm_4-DASATINIB72hr", "Mesoderm_4-UNTREATED72hr", "Mesoderm_3-DASATINIB72hr", "Mesoderm_3-UNTREATED72hr", "Mesoderm_2-DASATINIB72hr", "Mesoderm_2-UNTREATED72hr", "Mesoderm_1-DASATINIB72hr", "Mesoderm_1-UNTREATED72hr", "Ectoderm-DASATINIB72hr", "Ectoderm-UNTREATED72hr")
pdf("output/seurat/DotPlot_SCT_UNTREATED72hr_DASATINIB72hr_top10_V2.pdf", width=26, height=5)
DotPlot(humangastruloid.combined.sct, features = marker_genes, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# Display the top 10 CONSERVED marker genes of each cluster
Idents(humangastruloid.combined.sct) <- "cluster.annot"

## DEGs cluster versus all other
Ectoderm.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Ectoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Ectoderm")
Mesoderm_1.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Mesoderm_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mesoderm_1")
Mesoderm_2.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Mesoderm_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mesoderm_2")
Mesoderm_3.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Mesoderm_3", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mesoderm_3")
Mesoderm_4.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Mesoderm_4", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mesoderm_4")
Endoderm.conserved <- FindConservedMarkers(humangastruloid.combined.sct, assay = "RNA", ident.1 = "Endoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Endoderm")
## Combine all conserved markers into one data frame
all_conserved <- bind_rows(Ectoderm.conserved, Mesoderm_1.conserved, Mesoderm_2.conserved, Mesoderm_3.conserved, Mesoderm_4.conserved, Endoderm.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Ectoderm", "Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Endoderm"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Ectoderm", "Mesoderm_1", "Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Endoderm")))


## Write the top 10 conserved markers for each cluster to a file
write.table(top10_conserved, file = "output/seurat/srat_top10_conserved_markers_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Visualize the top 10 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(humangastruloid.combined.sct) <- c("Endoderm", "Mesoderm_4", "Mesoderm_3", "Mesoderm_2", "Mesoderm_1", "Ectoderm")
pdf("output/seurat/DotPlot_SCT_top10_conserved_V2.pdf", width=26, height=4)
DotPlot(humangastruloid.combined.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# save
saveRDS(humangastruloid.combined.sct, file = "output/seurat/humangastruloid.combined.sct_V2.rds")
humangastruloid.combined.sct <- readRDS(file = "output/seurat/humangastruloid.combined.sct_V2.rds")







# test EasyCellType
# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Hs.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_UNTREATED72hr_DASATINIB72hr_all_markers_V2.txt", header = TRUE, row.names = 1)
### Filter either WT or cYAPKO
all_markers <- all_markers[grepl("UNTREATED72hr$", all_markers$cluster), ]

## Convert geneSymbol to EntrezID
all_markers$entrezid <- mapIds(org.Hs.eg.db,
                           keys=all_markers$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
all_markers <- na.omit(all_markers)

## Sort the datafram (data frame containing Entrez IDs, clusters and expression scores)

all_markers_sort <- data.frame(gene=all_markers$entrezid, cluster=all_markers$cluster, 
                      score=all_markers$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(all_markers_sort[, 1:3])

## Run the enrihcment analysis
annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or Panglaodb or Clustermole
                    species="Human", #  Human or Mouse
                    tissue=c("Embryo", "Ee"), p_cut=0.3,   # many other tissue available
                    test="GSEA")    # GSEA or Fisher?

annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or Panglaodb or Clustermole
                    species="Human", #  Human or Mouse
                    tissue=c("Embryo", "Embryoid body", "Embryonic brain", "Embryonic stem cell", "Embryonic prefrontal cortex"), p_cut=0.3,   # many other tissue available
                    test="GSEA")    # GSEA or Fisher?
### save output
annot.GSEA_df <- data.frame(annot.GSEA$Mesoderm_4)

write.table(annot.GSEA_df, "output/seurat/EasyCellType_GSEA_Mesoderm_4_UNTREATED72hr_allEmbryoTissue_V2.txt", sep="\t", row.names=FALSE, quote=FALSE)
###


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or Panglaodb or Clustermole
                    species="Human", #  Human or Mouse
                    p_cut=0.3,   # many other tissue available
                    test="GSEA")  

## plots
pdf("output/seurat/EasyCellType_dotplot_UNTREATED72hr_V2.pdf", width=10, height=5)
pdf("output/seurat/EasyCellType_dotplot_UNTREATED72hr_allEmbryoTissue_V2.pdf", width=5, height=5)
pdf("output/seurat/EasyCellType_dotplot_UNTREATED72hr_noTissue_V2.pdf", width=5, height=5)

plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or Panglaodb or Clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Embryo"),   # many other tissue available
                    test="fisher")    # GSEA or fisher

pdf("output/seurat/EasyCellType_dotplot_SCT_control_fisher.pdf", width=6, height=8)
plot_dot(test="fisher", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



pdf("output/seurat/EasyCellType_barplot_SCT_control_cYAPKO.pdf", width=10, height=5)
plot_bar(test="GSEA", annot.GSEA)
dev.off()







# Functional analysis GO pathway
humangastruloid.combined.sct <- readRDS(file = "output/seurat/humangastruloid.combined.sct_V2.rds")

library("ReactomeGSA")
library("ggrepel")
library("RColorBrewer")

DefaultAssay(humangastruloid.combined.sct) <- "RNA" # For 

# separate condition from my seurat object
## Subset Seurat object based on condition
humangastruloid.combined.sct.WT <- subset(humangastruloid.combined.sct, subset = condition == "WT")
humangastruloid.combined.sct.cYAPKO <- subset(humangastruloid.combined.sct, subset = condition == "cYAPKO")

gsva_result <- analyse_sc_clusters(humangastruloid.combined.sct.cYAPKO, verbose = TRUE)

gsva_result <- analyse_sc_clusters(humangastruloid.combined.sct, verbose = TRUE)
pathway_expression <- pathways(gsva_result)
## maximum difference in expression for every pathway
### find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min
### sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

## Plot
### Expression for a single pathway
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])
### Heatmap pathway
pdf("output/seurat/ReactomeGSA_heatmap_UNTREATED72hr_DASATINIB72hr.pdf", width=15, height=10)
plot_gsva_heatmap(gsva_result, max_pathways = 25, margins = c(12,40), truncate_names = FALSE, col = colorRampPalette(c("blue", "white", "red"))(100)) # ,   scale = "row"
dev.off()

### Pathway-level PCA
pdf("output/seurat/ReactomeGSA_PCA_UNTREATED72hr_DASATINIB72hr.pdf", width=15, height=10)
plot_gsva_pca(gsva_result) +
  geom_text_repel(aes(label = sample), 
                   box.padding = 0.35, 
                   point.padding = 0.5, 
                   segment.color = 'grey50')
dev.off()



## SCPA
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
humangastruloid.combined.sct <- readRDS(file = "output/seurat/humangastruloid.combined.sct_V2.rds")

DefaultAssay(humangastruloid.combined.sct) <- "RNA" # Recommended 


# Isolate untreated and dasaatinib
pathways <- msigdbr("Homo sapiens", "C2") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise
pathways$REACTOME_SIGNALING_BY_RETINOIC_ACID$Genes

# Compare pathway activity between condition within  
cell_types <- unique(humangastruloid.combined.sct$cluster.annot)
humangastruloid.combined.sct_split <- SplitObject(humangastruloid.combined.sct, split.by = "condition")

# SCPA comparison

# Code to save output for each cell type comparison
clusters = c(
"Mesoderm_1",
"Mesoderm_2",
"Mesoderm_3",
"Endoderm",
"Ectoderm",
"Mesoderm_4"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  UNTREATED72hr <- seurat_extract(humangastruloid.combined.sct,
                       meta1 = "condition", value_meta1 = "UNTREATED72hr",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  DASATINIB72hr <- seurat_extract(humangastruloid.combined.sct,
                           meta1 = "condition", value_meta1 = "DASATINIB72hr",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  UNTREATED72hr_DASATINIB72hr <- compare_pathways(samples = list(UNTREATED72hr, DASATINIB72hr),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_humangastruloid_", cluster, ".txt")
  write.table(UNTREATED72hr_DASATINIB72hr, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}


# BALOW CLEAN CODE; dotplot>enrichplot>heatmap>violinplot
# DOT PLOT pepresetnation
## load all the comparison for each cell type (FC qval information)
clusters = c(
"Mesoderm_1",
"Mesoderm_2",
"Mesoderm_3",
"Endoderm",
"Ectoderm",
"Mesoderm_4"
)
## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_humangastruloid_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster))
# --> import the gene pathways (used table from Conchi `Pathwyas of interest*.xlsx`)
pathways = c(
"BIOCARTA_WNT_PATHWAY",
"KEGG_WNT_SIGNALING_PATHWAY",
"WP_WNT_SIGNALING",
"WP_WNT_SIGNALING_PATHWAY",
"WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY",
"PID_WNT_CANONICAL_PATHWAY",
"PID_WNT_NONCANONICAL_PATHWAY",
"PID_WNT_SIGNALING_PATHWAY",
"WILLERT_WNT_SIGNALING",
"WNT_SIGNALING",
"REACTOME_SIGNALING_BY_WNT",
"REACTOME_WNT5A_DEPENDENT_INTERNALIZATION_OF_FZD4",
"REACTOME_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING",
"REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING",
"REACTOME_NEGATIVE_REGULATION_OF_TCF_DEPENDENT_SIGNALING_BY_WNT_LIGAND_ANTAGONISTS",
"REACTOME_DEACTIVATION_OF_THE_BETA_CATENIN_TRANSACTIVATING_COMPLEX",
"REACTOME_FORMATION_OF_THE_BETA_CATENIN_TCF_TRANSACTIVATING_COMPLEX",
"REACTOME_DEGRADATION_OF_BETA_CATENIN_BY_THE_DESTRUCTION_COMPLEX",
"REACTOME_DEGRADATION_OF_AXIN",
"REACTOME_BETA_CATENIN_PHOSPHORYLATION_CASCADE",
"REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT",
"REACTOME_NEGATIVE_REGULATION_OF_TCF_DEPENDENT_SIGNALING_BY_WNT_LIGAND_ANTAGONISTS",
"REACTOME_WNT5A_DEPENDENT_INTERNALIZATION_OF_FZD4"
)
pathways = c(
"WP_TGFBETA_RECEPTOR_SIGNALING",
"WP_TGFBETA_SIGNALING_PATHWAY",
"KEGG_TGF_BETA_SIGNALING_PATHWAY",
"BIOCARTA_TGFB_PATHWAY",
"KARAKAS_TGFB1_SIGNALING",
"PID_TGFBR_PATHWAY",
"REACTOME_SIGNALING_BY_NODAL",
"REACTOME_DOWNREGULATION_OF_SMAD2_3_SMAD4_TRANSCRIPTIONAL_ACTIVITY",
"REACTOME_TRANSCRIPTIONAL_ACTIVITY_OF_SMAD2_SMAD3_SMAD4_HETEROTRIMER",
"REACTOME_SMAD2_SMAD3_SMAD4_HETEROTRIMER_REGULATES_TRANSCRIPTION",
"REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
"REACTOME_SIGNALING_BY_ACTIVIN",
"REACTOME_DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING",
"REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS",
"REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",
"REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
"WP_CANONICAL_AND_NONCANONICAL_TGFB_SIGNALING"
)
pathways = c(
"REACTOME_SIGNALING_BY_HIPPO",
"REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION",
"WP_HIPPO_SIGNALING_REGULATION_PATHWAYS",
"WP_MECHANOREGULATION_AND_PATHOLOGY_OF_YAPTAZ_VIA_HIPPO_AND_NONHIPPO_MECHANISMS",
"WP_HIPPOYAP_SIGNALING_PATHWAY",
"REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION",
"WP_MECHANOREGULATION_AND_PATHOLOGY_OF_YAPTAZ_VIA_HIPPO_AND_NONHIPPO_MECHANISMS"
)
pathways = c(
"REACTOME_SIGNALING_BY_RETINOIC_ACID",
"REACTOME_THE_CANONICAL_RETINOID_CYCLE_IN_RODS_TWILIGHT_VISION",
"KEGG_RETINOL_METABOLISM",
"PID_RETINOIC_ACID_PATHWAY",
"WP_4HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION"
)
pathways = c(
"REACTOME_SIGNALING_BY_NOTCH",
"REACTOME_SIGNALING_BY_NOTCH1",
"REACTOME_SIGNALING_BY_NOTCH2",
"REACTOME_SIGNALING_BY_NOTCH3",
"REACTOME_SIGNALING_BY_NOTCH4",
"REACTOME_NOTCH_HLH_TRANSCRIPTION_PATHWAY" ,
"REACTOME_PRE_NOTCH_PROCESSING_IN_GOLGI",
"REACTOME_PRE_NOTCH_EXPRESSION_AND_PROCESSING",
"REACTOME_NOTCH2_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",
"REACTOME_NOTCH3_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",
"REACTOME_ACTIVATED_NOTCH1_TRANSMITS_SIGNAL_TO_THE_NUCLEUS",
"REACTOME_NOTCH3_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION",
"REACTOME_NOTCH4_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION",
"KEGG_NOTCH_SIGNALING_PATHWAY" ,
"PID_NOTCH_PATHWAY" ,
"WP_NOTCH_SIGNALING" ,
"WP_NOTCH_SIGNALING_PATHWAY" ,
"WP_CANONICAL_AND_NONCANONICAL_NOTCH_SIGNALING"
)
pathways = c(
"REACTOME_SIGNALING_BY_BMP",
"WP_BMP_SIGNALING_IN_EYELID_DEVELOPMENT",
"PID_BMP_PATHWAY"
)
pathways = c(
"REACTOME_HDMS_DEMETHYLATE_HISTONES",
"REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA",
"REACTOME_HDACS_DEACETYLATE_HISTONES",
"REACTOME_RMTS_METHYLATE_HISTONE_ARGININES",
"REACTOME_HATS_ACETYLATE_HISTONES",
"REACTOME_PKMTS_METHYLATE_HISTONE_LYSINES",
"WP_INTERACTOME_OF_POLYCOMB_REPRESSIVE_COMPLEX_2_PRC2",
"BIOCARTA_HDAC_PATHWAY"
)
pathways = c(
"WP_ENDODERM_DIFFERENTIATION",
"WP_MESODERMAL_COMMITMENT_PATHWAY",
"WP_ECTODERM_DIFFERENTIATION"
)
pathways = c(
"REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_DOPAMINE_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_SEROTONIN_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
"WP_COMPLEMENT_SYSTEM_IN_NEURONAL_DEVELOPMENT_AND_PLASTICITY",
"REACTOME_MECP2_REGULATES_NEURONAL_RECEPTORS_AND_CHANNELS",
"REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS",
"WP_MAJOR_RECEPTORS_TARGETED_BY_EPINEPHRINE_AND_NOREPINEPHRINE",
"REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
"REACTOME_NEUREXINS_AND_NEUROLIGINS",
"REACTOME_ASSEMBLY_AND_CELL_SURFACE_PRESENTATION_OF_NMDA_RECEPTORS",
"REACTOME_NEGATIVE_REGULATION_OF_NMDA_RECEPTOR_MEDIATED_NEURONAL_TRANSMISSION",
"REACTOME_ACTIVATION_OF_NMDA_RECEPTORS_AND_POSTSYNAPTIC_EVENTS",
"KEGG_NEUROTROPHIN_SIGNALING_PATHWAY",
"REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES"
)
pathways = c(
"REACTOME_BLOOD_GROUP_SYSTEMS_BIOSYNTHESIS",
"REACTOME_HEME_SIGNALING",
"REACTOME_SCAVENGING_OF_HEME_FROM_PLASMA",
"REACTOME_RESPONSE_OF_EIF2AK1_HRI_TO_HEME_DEFICIENCY",
"KEGG_VEGF_SIGNALING_PATHWAY",
"REACTOME_SIGNALING_BY_VEGF",
"REACTOME_VEGFR2_MEDIATED_VASCULAR_PERMEABILITY",
"REACTOME_VEGFR2_MEDIATED_CELL_PROLIFERATION",
"BIOCARTA_VEGF_PATHWAY",
"WP_VEGFAVEGFR2_SIGNALING_PATHWAY",
"PID_VEGFR1_PATHWAY",
"PID_VEGFR1_2_PATHWAY"
)
pathways = c(
"REACTOME_SIGNALING_BY_FGFR",
"REACTOME_SIGNALING_BY_FGFR1",
"REACTOME_FGFR1_LIGAND_BINDING_AND_ACTIVATION",
"REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR1",
"REACTOME_SHC_MEDIATED_CASCADE_FGFR1",
"REACTOME_FRS_MEDIATED_FGFR1_SIGNALING",
"REACTOME_PI_3K_CASCADE_FGFR1",
"REACTOME_NEGATIVE_REGULATION_OF_FGFR1_SIGNALING",
"REACTOME_SIGNALING_BY_FGFR2",
"REACTOME_FGFR2_LIGAND_BINDING_AND_ACTIVATION",
"REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR2",
"REACTOME_SHC_MEDIATED_CASCADE_FGFR2",
"REACTOME_FRS_MEDIATED_FGFR2_SIGNALING",
"REACTOME_PI_3K_CASCADE_FGFR2",
"REACTOME_PHOSPHOLIPASE_C_MEDIATED_CASCADE_FGFR2",
"REACTOME_NEGATIVE_REGULATION_OF_FGFR2_SIGNALING",
"REACTOME_FGFR2_ALTERNATIVE_SPLICING",
"REACTOME_SIGNALING_BY_FGFR3",
"REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR3",
"REACTOME_SHC_MEDIATED_CASCADE_FGFR3",
"REACTOME_FRS_MEDIATED_FGFR3_SIGNALING",
"REACTOME_PI_3K_CASCADE_FGFR3",
"REACTOME_NEGATIVE_REGULATION_OF_FGFR3_SIGNALING",
"REACTOME_SIGNALING_BY_FGFR4",
"REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR4",
"REACTOME_SHC_MEDIATED_CASCADE_FGFR4",
"REACTOME_FRS_MEDIATED_FGFR4_SIGNALING",
"REACTOME_PI_3K_CASCADE_FGFR4",
"REACTOME_PHOSPHOLIPASE_C_MEDIATED_CASCADE_FGFR4",
"REACTOME_NEGATIVE_REGULATION_OF_FGFR4_SIGNALING"
)
all_data_pathways = as_tibble(all_data) %>% 
  filter(Pathway %in% pathways)

# tidy the df
all_data_pathways_tidy <- all_data_pathways %>%
  filter(qval >= 1.4) 

# Dot plot
pdf("output/Pathway/dotplot_WNT_signaling_humangastruloid.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_NODAL_TGFB_signaling_humangastruloid.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_HIPPO_signaling_humangastruloid.pdf", width=12, height=3)
pdf("output/Pathway/dotplot_RA_signaling_humangastruloid.pdf", width=12, height=2)
pdf("output/Pathway/dotplot_Notch_signaling_humangastruloid.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_BMP_signaling_humangastruloid.pdf", width=8, height=2)
pdf("output/Pathway/dotplot_Epigenetics_humangastruloid.pdf", width=10, height=3)
pdf("output/Pathway/dotplot_EndoMesoEcto_humangastruloid.pdf", width=10, height=3)
pdf("output/Pathway/dotplot_Neuro_humangastruloid.pdf", width=12, height=4)
pdf("output/Pathway/dotplot_Blood_humangastruloid.pdf", width=12, height=4)
pdf("output/Pathway/dotplot_FGF_humangastruloid.pdf", width=10, height=6)

ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, color = -FC), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8), guide = "legend") +  # Adjust the range for size if needed
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, guide = "colourbar") +
  theme_bw() +
  labs(size = "q-value", color = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# REFINE COLOR
custom_color <- function(fc_value){
  ifelse(fc_value >= 5, "strong_blue", 
         ifelse(fc_value > 2, "light_blue",
                ifelse(fc_value <= -5, "strong_red", 
                       ifelse(fc_value < -2, "light_red", "grey"))))
}

# Add a column for this custom color
all_data_pathways_tidy <- all_data_pathways_tidy %>%
  mutate(custom_col = sapply(FC, custom_color))

pdf("output/Pathway/dotplot_WNT_signaling_FCtresh_humangastruloid.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_NODAL_TGFB_signaling_FCtresh_humangastruloid.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_HIPPO_signaling_FCtresh_humangastruloid.pdf", width=12, height=3)
pdf("output/Pathway/dotplot_RA_signaling_FCtresh_humangastruloid.pdf", width=12, height=2)
pdf("output/Pathway/dotplot_Notch_signaling_FCtresh_humangastruloid.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_BMP_signaling_FCtresh_humangastruloid.pdf", width=12, height=2)
pdf("output/Pathway/dotplot_Epigenetics_FCtresh_humangastruloid.pdf", width=12, height=3)
pdf("output/Pathway/dotplot_EndoMesoEcto_FCtresh_humangastruloid.pdf", width=12, height=2)
pdf("output/Pathway/dotplot_Neuro_FCtresh_humangastruloid.pdf", width=12, height=4)
pdf("output/Pathway/dotplot_Blood_FCtresh_humangastruloid.pdf", width=12, height=4)
pdf("output/Pathway/dotplot_FGF_FCtresh_humangastruloid.pdf", width=12, height=6)


ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, color = custom_col), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = c("strong_blue" = "dodgerblue4",
                                "light_blue" = "lightblue2",
                                "grey" = "grey",
                                "light_red" = "indianred1",
                                "strong_red" = "red3")) +
  theme_bw() +
  labs(size = "q-value", color = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




# heatmap
pathways_heatmap <- msigdbr("Homo sapiens", "C2") %>%
format_pathways()
names(pathways_heatmap) <- sapply(pathways_heatmap, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

##pick pathway to represent
pathways_heatmap$REACTOME_SIGNALING_BY_FGFR$Genes
genes_of_interest <- pathways_heatmap$REACTOME_SIGNALING_BY_FGFR$Genes
##

genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(humangastruloid.combined.sct@assays$RNA@data)]

### extract UNTREATED72hr and DASATINIB72hr gene expression values from RNA assay
UNTREATED72hr_expression <- humangastruloid.combined.sct@assays$RNA@data[genes_of_interest, colnames(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "UNTREATED72hr" & humangastruloid.combined.sct$cluster.annot == "Endoderm"]]
DASATINIB72hr_expression <- humangastruloid.combined.sct@assays$RNA@data[genes_of_interest, colnames(humangastruloid.combined.sct)[humangastruloid.combined.sct$condition == "DASATINIB72hr" & humangastruloid.combined.sct$cluster.annot == "Endoderm"]]

### mean expression values for each gene
UNTREATED72hr_mean <- rowMeans(UNTREATED72hr_expression)
DASATINIB72hr_mean <- rowMeans(DASATINIB72hr_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  UNTREATED72hr = UNTREATED72hr_mean,
  DASATINIB72hr = DASATINIB72hr_mean
) %>% 
pivot_longer(cols = c(UNTREATED72hr, DASATINIB72hr), names_to = "Condition", values_to = "Expression")

## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(UNTREATED72hr_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)

# if few genes:
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_WNT_Mesoderm_2_humangastruloid.pdf", width=10, height=3)

pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS_Mesoderm_1_humangastruloid.pdf", width=10, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_NODAL_Mesoderm_3_humangastruloid.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_HIPPO_SIGNALING_REGULATION_PATHWAYS_Endoderm_humangastruloid.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_RETINOIC_ACID_Mesoderm_3_humangastruloid.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_PID_BMP_PATHWAY_Endoderm_humangastruloid.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_BIOCARTA_HDAC_PATHWAY_Mesoderm_2_humangastruloid.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_MESODERMAL_COMMITMENT_PATHWAY_Mixed_Mesoderm.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_FGFR_Blood_Progenitor_1.pdf", width=6, height=3)

library("viridis")
ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()

# if many genesL
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_WNT_Mesoderm_2_humangastruloid.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS_Mesoderm_1_humangastruloid.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_HIPPO_SIGNALING_REGULATION_PATHWAYS_Endoderm_humangastruloid.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_NOTCH_Endoderm_humangastruloid.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_NOTCH_Mesoderm_2_humangastruloid.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_HATS_ACETYLATE_HISTONES_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_ENDODERM_DIFFERENTIATION_Mesoderm_2_humangastruloid.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION_Mesoderm_2_humangastruloid.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_VEGFAVEGFR2_SIGNALING_PATHWAY_Mesoderm_2_humangastruloid.pdf", width=5, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_FGFR_Endoderm_humangastruloid.pdf", width=3, height=2)

library("viridis")
ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") 
dev.off()




# GSEA plot 
library("fgsea")
## Re-calculate DEGs keeping ALL genes
humangastruloid.combined.sct$celltype.stim <- paste(humangastruloid.combined.sct$cluster.annot, humangastruloid.combined.sct$condition,
    sep = "-")
Idents(humangastruloid.combined.sct) <- "celltype.stim"

Mesoderm_3 <- FindMarkers(humangastruloid.combined.sct, ident.1 = "Mesoderm_3-DASATINIB72hr", ident.2 = "Mesoderm_3-UNTREATED72hr",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 

### save output
write.table(Mesoderm_3, file = "output/seurat/Mesoderm_3-DASATINIB72hr_response_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)



Endoderm <- read.delim("output/seurat/Endoderm-DASATINIB72hr_response_allGenes.txt", header = TRUE, row.names = 1)
Mesoderm_2 <- read.delim("output/seurat/Mesoderm_2-DASATINIB72hr_response_allGenes.txt", header = TRUE, row.names = 1)
Mesoderm_1 <- read.delim("output/seurat/Mesoderm_1-DASATINIB72hr_response_allGenes.txt", header = TRUE, row.names = 1)
Mesoderm_3 <- read.delim("output/seurat/Mesoderm_3-DASATINIB72hr_response_allGenes.txt", header = TRUE, row.names = 1)


###

## load list of genes to test
pathways <- msigdbr("Homo sapiens", "C2") 
fgsea_sets <- pathways %>% split(x = .$gene_symbol, f = .$gs_name)

## Rank genes based on FC
genes <- Endoderm %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

## plot GSEA
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_WNT_Mesoderm_2_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_KEGG_WNT_SIGNALING_PATHWAY_Endoderm_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS_Mesoderm_1_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_NODAL_Mesoderm_3_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_WP_HIPPO_SIGNALING_REGULATION_PATHWAYS_Endoderm_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_RETINOIC_ACID_Mesoderm_3_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_NOTCH_Endoderm_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_NOTCH_Mesoderm_2_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_PID_BMP_PATHWAY_Endoderm_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_BIOCARTA_HDAC_PATHWAY_Mesoderm_2_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_WP_ENDODERM_DIFFERENTIATION_Mesoderm_2_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION_Mesoderm_2_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_WP_VEGFAVEGFR2_SIGNALING_PATHWAY_Mesoderm_2_humangastruloid.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_FGFR_Endoderm_humangastruloid.pdf", width=5, height=3)


plotEnrichment(fgsea_sets[["REACTOME_SIGNALING_BY_FGFR"]],
               ranks) + labs(title="REACTOME_SIGNALING_BY_FGFR-Endoderm") +
               theme_bw()
dev.off()






```

I now use the same number of dimensions for SCTransform and data integration steps (UMAP, neighbor, cluster...) for Untreated and Dasatinib + used RNA assay for all DEGs analysis (eg. FindMarker) 

--> For vizualization both SCT count and RNA (normalize and scaled) lead to similar results. I used and show SCT count per default

--> Clustering is optimal, 4 Mesoderm sub-clusters (could be identify more carefully looking at our unbiased marker genes; notably Mesoderm_2) and ectoderm and mesoderm
----> Similarly, Ectoderm was hard to individualize but optimization of nb of dimensions, k.param and resolution worked
----> Endoderm, very easy to individualize

--> Downsampling and bootstraps have been used to identify the overall changes in cell type distribution
----> Revealed Decrease Endoderm; and Increase Mesoderm_3_4 with the treatment

--> Could test the Seurat integration v5 that allow user to compare different integration method in a one line code; see here XXX

**To do next:**
- identify the Mesoderm sub-types (check unbiased marker gene list + Conchi meeting)
- Test different integration method (if Conchi not fully satisfied)
- Provide clean gene list to Conchi (see what she want after meeting 20230801)


--> Share to Conchi the Conserved Marker list (`srat_all_conserved_markers_V2.xlsx`). To avoid confusion, I did some filtering: For each cell type; I only keep log2FC positive (= correspond to gene more highly express in this cell types) and I told her to filter per pvalue which is the max_pvalue. Like this, she will only see the highly express genes in each cluster





# Embryo analysis in Seurat

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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

# soupX decontamination
## Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('embryo_control_e775_mice/outs') 
sc = load10X('embryo_cYAPKO_e775_mice/outs') 

## Assess % of conta
pdf("output/soupX/autoEstCont_embryo_control_e775_mice.pdf", width=10, height=10)
pdf("output/soupX/autoEstCont_embryo_cYAPKO_e775_mice.pdf", width=10, height=10)
sc = autoEstCont(sc)
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/out_embryo_control_e775_mice.RData")
save(out, file = "output/soupX/out_embryo_cYAPKO_e775_mice.RData")
## Load the matrix and Create SEURAT object
load("output/soupX/out_embryo_control_e775_mice.RData")
srat_WT <- CreateSeuratObject(counts = out, project = "WT") # 32,285 features across 5,568 samples

load("output/soupX/out_embryo_cYAPKO_e775_mice.RData")
srat_cYAPKO <- CreateSeuratObject(counts = out, project = "cYAPKO") # 32,285 features across 4,504 samples

# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
srat_WT[["percent.mt"]] <- PercentageFeatureSet(srat_WT, pattern = "^mt-")
srat_WT[["percent.rb"]] <- PercentageFeatureSet(srat_WT, pattern = "^Rp[sl]")

srat_cYAPKO[["percent.mt"]] <- PercentageFeatureSet(srat_cYAPKO, pattern = "^mt-")
srat_cYAPKO[["percent.rb"]] <- PercentageFeatureSet(srat_cYAPKO, pattern = "^Rp[sl]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/embryo_control.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_WT <- AddMetaData(srat_WT,doublets)
srat_WT$Doublet_score <- as.numeric(srat_WT$Doublet_score) # make score as numeric
head(srat_WT[[]])

doublets <- read.table("output/doublets/embryo_cYAPKO.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_cYAPKO <- AddMetaData(srat_cYAPKO,doublets)
srat_cYAPKO$Doublet_score <- as.numeric(srat_cYAPKO$Doublet_score) # make score as numeric
head(srat_cYAPKO[[]])



## Plot
pdf("output/seurat/VlnPlot_QC_embryo_control.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QC_embryo_cYAPKO.pdf", width=10, height=6)
VlnPlot(srat_cYAPKO, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_cYAPKO.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_control.pdf", width=5, height=5)
FeatureScatter(srat_WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_cYAPKO.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_control.pdf", width=5, height=5)
FeatureScatter(srat_WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
pdf("output/seurat/FeatureScatter_QC_rb_mt_cYAPKO.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_rb_mt_control.pdf", width=5, height=5)
FeatureScatter(srat_WT, feature1 = "percent.rb", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_mt_doublet_cYAPKO.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_mt_doublet_control.pdf", width=5, height=5)
FeatureScatter(srat_WT, feature1 = "percent.mt", feature2 = "Doublet_score")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_cYAPKO.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_control.pdf", width=5, height=5)
FeatureScatter(srat_WT, feature1 = "nFeature_RNA", feature2 = "Doublet_score")
dev.off()



## After seeing the plot; add QC information in our seurat object
## V1 QC; optimal with vst V1; dim 19 k param 70 es 0.9 --> THE WINNER pro winner
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC == 'Pass','Low_nFeature',srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_WT@meta.data$QC,sep = ','),srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$percent.mt > 25 & srat_WT@meta.data$QC == 'Pass','High_MT',srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'High_MT',paste('High_MT',srat_WT@meta.data$QC,sep = ','),srat_WT@meta.data$QC)
table(srat_WT[['QC']])
## 
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC == 'Pass','Low_nFeature',srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_cYAPKO@meta.data$QC,sep = ','),srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$percent.mt > 25 & srat_cYAPKO@meta.data$QC == 'Pass','High_MT',srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'High_MT',paste('High_MT',srat_cYAPKO@meta.data$QC,sep = ','),srat_cYAPKO@meta.data$QC)
table(srat_cYAPKO[['QC']])


## V2 more stringeant; define after looking at integration
### Mark doublets
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
### Mark cells with nFeature_RNA < 2000 as Low_nFeature
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC == 'Pass', 'Low_nFeature', srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', srat_WT@meta.data$QC, sep = ','), srat_WT@meta.data$QC)
### Mark cells with percent.mt > 15 as High_MT
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$percent.mt > 15 & srat_WT@meta.data$QC == 'Pass', 'High_MT', srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$percent.mt > 15 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'High_MT', paste('High_MT', srat_WT@meta.data$QC, sep = ','), srat_WT@meta.data$QC)
### NEW: Mark cells with nFeature_RNA > 10000 as High_nFeature
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA > 10000 & srat_WT@meta.data$QC == 'Pass', 'High_nFeature', srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA > 10000 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'High_nFeature', paste('High_nFeature', srat_WT@meta.data$QC, sep = ','), srat_WT@meta.data$QC)
table(srat_WT[['QC']])
### Mark doublets
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
### Mark cells with nFeature_RNA < 2000 as Low_nFeature
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC == 'Pass', 'Low_nFeature', srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', srat_cYAPKO@meta.data$QC, sep = ','), srat_cYAPKO@meta.data$QC)
### Mark cells with percent.mt > 15 as High_MT
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$percent.mt > 15 & srat_cYAPKO@meta.data$QC == 'Pass', 'High_MT', srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$percent.mt > 15 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'High_MT', paste('High_MT', srat_cYAPKO@meta.data$QC, sep = ','), srat_cYAPKO@meta.data$QC)
### NEW: Mark cells with nFeature_RNA > 10000 as High_nFeature
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA > 10000 & srat_cYAPKO@meta.data$QC == 'Pass', 'High_nFeature', srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA > 10000 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'High_nFeature', paste('High_nFeature', srat_cYAPKO@meta.data$QC, sep = ','), srat_cYAPKO@meta.data$QC)
table(srat_cYAPKO[['QC']])


## V2=3  even more stringeant; define after looking at integration
### Mark doublets
# Mark doublets
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
# Mark cells with nFeature_RNA < 2000 as Low_nFeature
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC == 'Pass', 'Low_nFeature', srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', srat_WT@meta.data$QC, sep = ','), srat_WT@meta.data$QC)
# Mark cells with percent.mt > 15 as High_MT
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$percent.mt > 15 & srat_WT@meta.data$QC == 'Pass', 'High_MT', srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$percent.mt > 15 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'High_MT', paste('High_MT', srat_WT@meta.data$QC, sep = ','), srat_WT@meta.data$QC)
# Mark cells with nFeature_RNA > 9000 as High_nFeature
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA > 8500 & srat_WT@meta.data$QC == 'Pass', 'High_nFeature', srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA > 8500 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'High_nFeature', paste('High_nFeature', srat_WT@meta.data$QC, sep = ','), srat_WT@meta.data$QC)
# NEW: Mark cells with percent.rb > 35 as High_RB
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$percent.rb > 35 & srat_WT@meta.data$QC == 'Pass', 'High_RB', srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$percent.rb > 35 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'High_RB', paste('High_RB', srat_WT@meta.data$QC, sep = ','), srat_WT@meta.data$QC)
# NEW: Mark cells with nCount_RNA > 125000 as High_nCount
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nCount_RNA > 125000 & srat_WT@meta.data$QC == 'Pass', 'High_nCount', srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nCount_RNA > 125000 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'High_nCount', paste('High_nCount', srat_WT@meta.data$QC, sep = ','), srat_WT@meta.data$QC)
# View the table
table(srat_WT[['QC']])

### Mark doublets
# Mark doublets
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
# Mark cells with nFeature_RNA < 2000 as Low_nFeature
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC == 'Pass', 'Low_nFeature', srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', srat_cYAPKO@meta.data$QC, sep = ','), srat_cYAPKO@meta.data$QC)
# Mark cells with percent.mt > 15 as High_MT
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$percent.mt > 15 & srat_cYAPKO@meta.data$QC == 'Pass', 'High_MT', srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$percent.mt > 15 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'High_MT', paste('High_MT', srat_cYAPKO@meta.data$QC, sep = ','), srat_cYAPKO@meta.data$QC)
# Mark cells with nFeature_RNA > 9000 as High_nFeature
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA > 8500 & srat_cYAPKO@meta.data$QC == 'Pass', 'High_nFeature', srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA > 8500 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'High_nFeature', paste('High_nFeature', srat_cYAPKO@meta.data$QC, sep = ','), srat_cYAPKO@meta.data$QC)
# NEW: Mark cells with percent.rb > 35 as High_RB
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$percent.rb > 35 & srat_cYAPKO@meta.data$QC == 'Pass', 'High_RB', srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$percent.rb > 35 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'High_RB', paste('High_RB', srat_cYAPKO@meta.data$QC, sep = ','), srat_cYAPKO@meta.data$QC)
# NEW: Mark cells with nCount_RNA > 125000 as High_nCount
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nCount_RNA > 125000 & srat_cYAPKO@meta.data$QC == 'Pass', 'High_nCount', srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nCount_RNA > 125000 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'High_nCount', paste('High_nCount', srat_cYAPKO@meta.data$QC, sep = ','), srat_cYAPKO@meta.data$QC)
# View the table
table(srat_cYAPKO[['QC']])




## subset my seurat object to only analyze the cells that pass the QC
srat_WT <- subset(srat_WT, subset = QC == 'Pass')
srat_cYAPKO <- subset(srat_cYAPKO, subset = QC == 'Pass')
srat_WT$condition <- "WT"
srat_cYAPKO$condition <- "cYAPKO"



mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
srat_WT <- NormalizeData(srat_WT, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_WT)
srat_WT <- ScaleData(srat_WT, features = all.genes) # zero-centres and scales it

srat_cYAPKO <- NormalizeData(srat_cYAPKO, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_cYAPKO)
srat_cYAPKO <- ScaleData(srat_cYAPKO, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
srat_WT <- CellCycleScoring(srat_WT, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_WT[[]]$Phase)
srat_cYAPKO <- CellCycleScoring(srat_cYAPKO, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_cYAPKO[[]]$Phase)

set.seed(42)





# Run SCTransform
## Version OK with 2000 treshold RNA
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 50, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.25, verbose = FALSE, algorithm = 4)
## perform OK; better than with 50 or 25 dim




## investigate to find optimal nb of dimension
### Vizualize the 1st 20 PC --> around 20 is good or 12
pdf("output/seurat/DimHeatmap_SCTV1_control.pdf", width=10, height=40)
DimHeatmap(srat_WT, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()
pdf("output/seurat/Elbow_SCTV1_control.pdf", width=10, height=10)
ElbowPlot(srat_WT) # 13, 15 
dev.off()
## so between 12-20 ( around 12-13 perfect)



srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 12, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:12, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:12, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.25, verbose = FALSE, algorithm = 4)
## looks worst without regtression
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4612, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000, vst.flavor = "v2") %>% 
    RunPCA(npcs = 10, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:10, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:10, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.25, verbose = FALSE, algorithm = 4)
## looks better with vst.flavor = v2
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4612, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000, vst.flavor = "v2") %>% 
    RunPCA(npcs = 50, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.2, verbose = FALSE, algorithm = 4)
## seems perform OK. As the other
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4612, vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000, vst.flavor = "v2") %>% 
    RunPCA(npcs = 50, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:50, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.2, verbose = FALSE, algorithm = 4)
## seems perform OK. 
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4612, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000, vst.flavor = "v2") %>% 
    RunPCA(npcs = 20, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:20, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.2, verbose = FALSE, algorithm = 4)
## seems perform good. 
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4612, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","nCount_RNA","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 2000) %>% 
    RunPCA(npcs = 15, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:15, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.2, verbose = FALSE, algorithm = 4)
## seems perform good.



## This is the best for now 
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score")) %>% 
    RunPCA(npcs = 19, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:19, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:19, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.5, verbose = FALSE, algorithm = 4)
##


srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score")) %>% 
    RunPCA(npcs = 30, verbose = FALSE) %>%  # CHANGE HERE NB OF DIM
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindNeighbors(reduction = "pca", k.param = 30, dims = 1:30, verbose = FALSE) %>% # CHANGE HERE NB OF DIM
    FindClusters(resolution = 0.5, verbose = FALSE, algorithm = 4)

# Run on srat_WT condition only and find optimal parameters; expend this to both then??


pdf("output/seurat/UMAP_control.pdf", width=10, height=6)
DimPlot(srat_WT, reduction = "umap", label=TRUE)
dev.off()


## Check some marker genes 
### Marker gene list from Conchi
epiblast = c("Nanog", "Pou5f1", "Dppa3") # watch out DPPA => Dppa3
endoderm = c("Gata2", "Gata1", "Hesx1", "Gata6", "Eomes", "Foxa2", "Sox17", "Prdm1", "Cxcr4", "Foxa1", "Stat3", "Hnf4a", "Tfcp2l1", "Foxm1", "Foxa2", "Foxa3") # watch out Crcr4 => Cxcr4; Tfcp21 => Tfcp2l1 
ectoderm = c("Rarb", "Sox9", "Sox1", "Sirt6", "Runx2", "Foxg1", "Sox2", "Nes", "Vim", "Id3", "Sox3", "Pou3f3", "Pou5f1", "Pou2f1") # Brn-1 => Pou3f3; 4-Oct =>  Pou5f1; Oct11 =>  Pou2f1
mesoderm = c("Lef1", "Cdx2", "Hoxa1", "Mixl1", "Sp5", "T", "Hmgb3", "Gata5", "Hmga2", "Smad1", "Tal1", "Foxf1", "Gata2", "Nodal", "Kdr", "Dll1", "Lhx1", "Aplnr", "Tbx6", "Mesp1", "Has2", "Pdgfra", "Hand1", "Lef1", "Gata6") # Left1 => Lef1 
### Updated marker gene list from Alex
Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # OK
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # OK
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # OK
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # OK
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # OK
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # OK
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # OK
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # OK
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # OK
Caudal_Neurectoderm = c("Olig3","Hes3","Lrrc7","Cntfr") # OK
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # OK
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # OK
Intermediate_Mesoderm = c("Bik","Arg1","Meox1","Ism1") # OK (very similar to Somitic_Mesoderm)
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # OK 
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # OK  (very similar to Somitic_Mesoderm and Intermediate_Mesoderm)
Notocord = c("Shh","Noto","Vtn","Fgg") # OK 
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # OK (very few cells)
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3") # OK could be better

## Watch out to Mixed_Mesoderm, Somitic_Mesoderm, Intermediate_Mesoderm; hardest to separate


DefaultAssay(srat_WT) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_control_Epiblast_PrimStreak.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Epiblast_PrimStreak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_ExE_Ectoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Nascent_Mesoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Somitic_Mesoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Somitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Caudal_Mesoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Caudal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Haematodenothelial_progenitors.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Haematodenothelial_progenitors, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Paraxial_Mesoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Paraxial_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Mesenchyme.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Mesenchyme, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Blood_Progenitor_1.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Blood_Progenitor_1, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Caudal_Neurectoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Caudal_Neurectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Pharyngeal_Mesoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Pharyngeal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Blood_Progenitor_2.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Blood_Progenitor_2, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Intermediate_Mesoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Intermediate_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Surface_Ectoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Surface_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Mixed_Mesoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Mixed_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Notocord.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Notocord, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Gut.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Gut, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Primordial_Germ_Cells.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Primordial_Germ_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()



### Marker gene list from Paper (https://omi-top.shinyapps.io/rna3/)
Epithelial_cells = c("Agr2", "Cldn10", "Cldn7", "Cyr61", "Epcam", "Krt15", "Krt18", "Krt5", "Krt8", "Ly6a", "Spint", "Spint2", "Sftpc", "Spp1")
Ectoderm = c("Cdh1", "Crabp2", "Dnmt3b", "Eno1", "Epcam", "Foxa2", "Hmga1", "Krt18", "Ldha", "Malat1", "Meg3", "Phlda2", "Pim2", "Pou5f1", "Trap1a", "Utf1")
Endoderm= c("Apoa1", "Car4", "Foxa2", "Phgdh", "Rbp4", "Rrm2", "Ttr")
Mesoderm= c("Cer1", "Lefty2", "Lin28a", "Mesp1", "Msgn1", "Phlda2")
Muscle= c("Actc1", "Atp2a1", "Cdh15", "Chrna1", "Msc")
Neuron_related = c("Bcam", "Cited1", "Col9a1", "Crabp2", "Eya1", "Mest", "Pax2", "Pclaf", "Rtn1", "Top2a", "Tubb3", "Uncx")
Immune_and_Blood_cells = c("Afp", "Alox12", "Alas2", "Arhgdib", "Asf1b", "C1qc", "C1qb", "Cpa3", "Coro1a", "Cd200r3", "Col1a1", "Dcn", "Emcn", "F12", "Fam214b", "Hba-a1", "Hba-a2", "Hba-bt", "Hba-x", "Hbb-y", "Hbb-bt", "Hp", "Ikzf3", "Il13", "Igll1", "Kif18a", "Krt18", "Ms4a6c", "Mfap2", "Nrep", "Phospho1", "Pdcd1", "Pou5f1", "Plac8", "Plvap", "Rad51ap1", "Rbp4", "Rac2", "Rpl18", "Srgn", "Slfn14", "Steap3", "Tnnt2", "Tmsb4x", "Tubb3", "Tyrobp", "Tnn2", "Unc93b1") 
Endothelial_cells = c("Adgrg6", "Ccl21a", "Cd34", "Cdh5", "Col1a1", "Cxcl1", "Cryab", "Egfl7", "Fabp4", "Gja4", "Icam2", "Plvap", "Ramp2", "Rrm2", "Rspo3", "Snca", "Sparc", "Upk3b")



DefaultAssay(srat_WT) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_control_Epithelial_cells.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Epithelial_cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Ectoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Endoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Mesoderm.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Muscle.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Muscle, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Neuron_related.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Neuron_related, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Immune_and_Blood_cells.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Immune_and_Blood_cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_Endothelial_cells.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = Endothelial_cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()





# Automatic cell-type annotation

https://biostatistics.mdanderson.org/shinyapps/EasyCellType/




## QC plot
## percent mt and rb
pdf("output/seurat/VlnPlot_SCT_control_mt_rb.pdf", width=10, height=7)
VlnPlot(srat_WT,features = c("percent.mt", "percent.rb")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## RNA
pdf("output/seurat/VlnPlot_SCT_control_count_feature.pdf", width=10, height=10)
VlnPlot(srat_WT,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## cell cycle
pdf("output/seurat/VlnPlot_SCT_control_cellCycle.pdf", width=10, height=10)
VlnPlot(srat_WT,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()  


## LogNormalize ##
srat_WT <- NormalizeData(srat_WT, normalization.method = "LogNormalize", scale.factor = 10000)
## Discover the 2000 first more variable genes
srat_WT <- FindVariableFeatures(srat_WT, selection.method = "vst", nfeatures = 2000)
## scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(srat_WT)
srat_WT <- ScaleData(srat_WT, features = all.genes,vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score","S.Score","G2M.Score"))
## PCA
srat_WT <- RunPCA(srat_WT, features = VariableFeatures(object = srat_WT))


# clustering
srat_WT <- FindNeighbors(srat_WT, dims = 1:15, k.param = 30)
srat_WT <- FindClusters(srat_WT, resolution = 0.2)
srat_WT <- RunUMAP(srat_WT, dims = 1:15, verbose = F)
table(srat_WT@meta.data$seurat_clusters) # to check cluster size

pdf("output/seurat/UMAP_LogNormalize_control.pdf", width=10, height=10)
DimPlot(srat_WT,label.size = 4,repel = T,label = T)
dev.off()



## Check some marker genes 
### Marker gene list
epiblast = c("Nanog", "Pou5f1", "Dppa3") # watch out DPPA => Dppa3
endoderm = c("Gata2", "Gata1", "Hesx1", "Gata6", "Eomes", "Foxa2", "Sox17", "Prdm1", "Cxcr4", "Foxa1", "Stat3", "Hnf4a", "Tfcp2l1", "Foxm1", "Foxa2", "Foxa3") # watch out Crcr4 => Cxcr4; Tfcp21 => Tfcp2l1 
ectoderm = c("Rarb", "Sox9", "Sox1", "Sirt6", "Runx2", "Foxg1", "Sox2", "Nes", "Vim", "Id3", "Sox3", "Pou3f3", "Pou5f1", "Pou2f1") # Brn-1 => Pou3f3; 4-Oct =>  Pou5f1; Oct11 =>  Pou2f1
mesoderm = c("Lef1", "Cdx2", "Hoxa1", "Mixl1", "Sp5", "T", "Hmgb3", "Gata5", "Hmga2", "Smad1", "Tal1", "Foxf1", "Gata2", "Nodal", "Kdr", "Dll1", "Lhx1", "Aplnr", "Tbx6", "Mesp1", "Has2", "Pdgfra", "Hand1", "Lef1", "Gata6") # Left1 => Lef1 



DefaultAssay(srat_WT) <- "RNA" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_LogNormalize_control_epiblast.pdf", width=10, height=10)
FeaturePlot(srat_WT, features = epiblast, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_control_endoderm.pdf", width=40, height=50)
FeaturePlot(srat_WT, features = endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_control_mesoderm.pdf", width=40, height=50)
FeaturePlot(srat_WT, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_LogNormalize_control_ectoderm.pdf", width=40, height=50)
FeaturePlot(srat_WT, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()



## QC plot
## percent mt and rb
pdf("output/seurat/VlnPlot_SCT_control_mt_rb.pdf", width=10, height=7)
VlnPlot(srat_WT,features = c("percent.mt", "percent.rb")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## RNA
pdf("output/seurat/VlnPlot_SCT_control_count_feature.pdf", width=10, height=10)
VlnPlot(srat_WT,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))
dev.off()
## cell cycle
pdf("output/seurat/VlnPlot_SCT_control_cellCycle.pdf", width=10, height=10)
VlnPlot(srat_WT,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()  







# Run SCTransform now with integration
## SCT V2
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000, vst.flavor = "v2") %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO <- SCTransform(srat_cYAPKO, method = "glmGamPoi", ncells = 3621, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000, vst.flavor = "v2") %>%
    RunPCA(npcs = 19, verbose = FALSE)

## SCT V1 (better) with stringeant QC (V2) 
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4745, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO <- SCTransform(srat_cYAPKO, method = "glmGamPoi", ncells = 3433, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)


## SCT V1 (better)  --> THE WINNER
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO <- SCTransform(srat_cYAPKO, method = "glmGamPoi", ncells = 3621, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)

## TEST feature regression; FAIL
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, vars.to.regress = c("nFeature_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO <- SCTransform(srat_cYAPKO, method = "glmGamPoi", ncells = 3621, vars.to.regress = c("nFeature_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)


    
## TEST RNA regression  --> THE WINNER PRO WINNER !!! RNA regression 
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO <- SCTransform(srat_cYAPKO, method = "glmGamPoi", ncells = 3621, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(srat_WT = srat_WT, srat_cYAPKO = srat_cYAPKO)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

embryo.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
embryo.combined.sct <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT")



# Perform integrated analysis (check active assay is 'integrated')
# pretty good
DefaultAssay(embryo.combined.sct) <- "integrated"
embryo.combined.sct <- RunPCA(embryo.combined.sct, verbose = FALSE, npcs = 19)
embryo.combined.sct <- RunUMAP(embryo.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryo.combined.sct <- FindNeighbors(embryo.combined.sct, reduction = "pca", k.param = 40, dims = 1:19)
embryo.combined.sct <- FindClusters(embryo.combined.sct, resolution = 0.7, verbose = FALSE, algorithm = 4)
embryo.combined.sct$condition <- factor(embryo.combined.sct$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st
pdf("output/seurat/UMAP_control_cYAPKO.pdf", width=10, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()
#
# Pretty good
DefaultAssay(embryo.combined.sct) <- "integrated"
embryo.combined.sct <- RunPCA(embryo.combined.sct, verbose = FALSE, npcs = 19)
embryo.combined.sct <- RunUMAP(embryo.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryo.combined.sct <- FindNeighbors(embryo.combined.sct, reduction = "pca", k.param = 20, dims = 1:19)
embryo.combined.sct <- FindClusters(embryo.combined.sct, resolution = 0.35, verbose = FALSE, algorithm = 4)
embryo.combined.sct$condition <- factor(embryo.combined.sct$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st
pdf("output/seurat/UMAP_control_cYAPKO.pdf", width=10, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()
##
DefaultAssay(embryo.combined.sct) <- "integrated"
embryo.combined.sct <- RunPCA(embryo.combined.sct, verbose = FALSE, npcs = 19)
embryo.combined.sct <- RunUMAP(embryo.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryo.combined.sct <- FindNeighbors(embryo.combined.sct, reduction = "pca", k.param = 20.5, dims = 1:19)
embryo.combined.sct <- FindClusters(embryo.combined.sct, resolution = 0.38, verbose = FALSE, algorithm = 4)
embryo.combined.sct$condition <- factor(embryo.combined.sct$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st
pdf("output/seurat/UMAP_control_cYAPKO.pdf", width=10, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()


# BEST FOR NOW WITHOUT RNA REGRESS
## individualize germ cells; but fail at separating   --> THE WINNER
DefaultAssay(embryo.combined.sct) <- "integrated"

embryo.combined.sct <- RunPCA(embryo.combined.sct, verbose = FALSE, npcs = 19)
embryo.combined.sct <- RunUMAP(embryo.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryo.combined.sct <- FindNeighbors(embryo.combined.sct, reduction = "pca", k.param = 20, dims = 1:19)
embryo.combined.sct <- FindClusters(embryo.combined.sct, resolution = 0.4, verbose = FALSE, algorithm = 4)

embryo.combined.sct$condition <- factor(embryo.combined.sct$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st

pdf("output/seurat/UMAP_control_cYAPKO.pdf", width=10, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()





## PLAY WITH RESOLUTION AFTER CONCHI MEETING 20231005
## individualize germ cells; but fail at separating  
DefaultAssay(embryo.combined.sct) <- "integrated"

embryo.combined.sct <- RunPCA(embryo.combined.sct, verbose = FALSE, npcs = 19)
embryo.combined.sct <- RunUMAP(embryo.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryo.combined.sct <- FindNeighbors(embryo.combined.sct, reduction = "pca", k.param = 5, dims = 1:19)
embryo.combined.sct <- FindClusters(embryo.combined.sct, resolution = 0.3, verbose = FALSE, algorithm = 4)

embryo.combined.sct$condition <- factor(embryo.combined.sct$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st

pdf("output/seurat/UMAP_control_cYAPKO_V2clust.pdf", width=10, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()



## PRO WINNER ! With RNA regression
## individualize germ cells; but fail at separating  
DefaultAssay(embryo.combined.sct) <- "integrated"

embryo.combined.sct <- RunPCA(embryo.combined.sct, verbose = FALSE, npcs = 19)
embryo.combined.sct <- RunUMAP(embryo.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryo.combined.sct <- FindNeighbors(embryo.combined.sct, reduction = "pca", k.param = 20, dims = 1:19)
embryo.combined.sct <- FindClusters(embryo.combined.sct, resolution = 0.65, verbose = FALSE, algorithm = 4)

embryo.combined.sct$condition <- factor(embryo.combined.sct$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st

pdf("output/seurat/UMAP_control_cYAPKO.pdf", width=10, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()




## Check some marker genes 
### Marker gene list
epiblast = c("Nanog", "Pou5f1", "DPPA")
endoderm = c("Gata2", "Gata1", "Hesx1", "Gata6", "Eomes", "Foxa2", "Sox17", "Prdm1", "Crcr4", "Foxa1", "Stat3", "Hnf4a", "Tfcp21", "Foxm1", "Foxa2", "Foxa3")
ectoderm = c("Rarb", "Sox9", "Sox1", "Sirt6", "Runx2", "Foxg1", "Sox2", "Nes", "VIM", "ID3", "Sox3", "Brn1", "4-Oct", "Oct11")
mesoderm = c("Lef1", "Cdx2", "Hoxa1", "Mixl1", "Sp5", "T", "Hmgb3", "Gata5", "Hmga2", "Smad1", "Tal1", "Foxf1", "Gata2", "NODAL", "KDR", "DLL1", "LXH1", "APLNR", "TBX6", "MESP1", "HAS2", "PDGFRA", "Hand1", "Left1", "Gata6")


DefaultAssay(embryo.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_epiblast.pdf", width=10, height=30)
FeaturePlot(embryo.combined.sct, features = epiblast, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_endoderm.pdf", width=10, height=80)
FeaturePlot(embryo.combined.sct, features = endoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_mesoderm.pdf", width=10, height=80)
FeaturePlot(embryo.combined.sct, features = mesoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_ectoderm.pdf", width=10, height=80)
FeaturePlot(embryo.combined.sct, features = ectoderm, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

### Updated marker gene list from Alex # 19 dim res 0.5 \ 19 dim  res 0.4 $ 19 dim res 0.7 kparam 70 * 19 dim res 0.9 kparam 70 + WINNER % pro-winner
Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # 1 \  1 $ 1 + 1 %  1
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # 3 \ 4 $ 3  + 4 %  2
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # 2 \ 3 $ 8  + 3 %  7
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # 4 \ 2 $ 2 + 2 %  3
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # 5 \ 5 $ 5 +  5 %  4
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # 2 \ 3 $ 4  +  3 %  8
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # 7 \ 7 $ 6  +    %  5
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # 8 \ 8 $ 9  +    %  9
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # 11  \  $ 11 +   %  11
Caudal_Neurectoderm = c("Olig3","Hes3","Lrrc7","Cntfr") # 6  \   $ 5  +  %  4
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # 5  \   $ 7 + %  6
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # 9  \  $  10 + %  12
Intermediate_Mesoderm = c("Bik","Arg1","Meox1","Ism1") # 4 \  2  $  2 + %  mix of 3 5 6 (let's not name it as more general than 3 5 6)
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # 12  \  $  12 + %  13
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # 4 and 10 \  2  $  2 + %  10 (and 3)
Notocord = c("Shh","Noto","Vtn","Fgg") #  13 \   $  13 ( wierd separated) + %  14
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # 14  \   $  14 + %  15
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3") # between 12 and 13    + %  18

## would put intermediate and somatic Mesoderm and Mixed together
### NEED TRY RES 0.45 so that we habe mesoderm together but separate the  Haematodenothelial_progenitors and Nascent_Mesoderm

DefaultAssay(embryo.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Epiblast_PrimStreak.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Epiblast_PrimStreak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_ExE_Ectoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Nascent_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Somitic_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Somitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Caudal_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Caudal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Haematodenothelial_progenitors.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Haematodenothelial_progenitors, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Paraxial_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Paraxial_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Mesenchyme.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Mesenchyme, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Blood_Progenitor_1.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Blood_Progenitor_1, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Caudal_Neurectoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Caudal_Neurectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Pharyngeal_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Pharyngeal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Blood_Progenitor_2.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Blood_Progenitor_2, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Intermediate_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Intermediate_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Surface_Ectoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Surface_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Mixed_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Mixed_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Notocord.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Notocord, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Gut.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Gut, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Primordial_Germ_Cells.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Primordial_Germ_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()



pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_YAP1.pdf", width=10, height=5)
FeaturePlot(embryo.combined.sct, features = "Yap1", max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()




## QC plot
## percent mt and rb
pdf("output/seurat/VlnPlot_SCT_control_cYAPKO_mt_rb.pdf", width=15, height=7)
VlnPlot(embryo.combined.sct,features = c("percent.mt", "percent.rb"), split.by = "condition") & 
  theme(plot.title = element_text(size=10))
dev.off()
## RNA
pdf("output/seurat/VlnPlot_SCT_control_cYAPKO_count_feature.pdf", width=15, height=10)
VlnPlot(embryo.combined.sct,features = c("nCount_RNA","nFeature_RNA"), split.by = "condition") & 
  theme(plot.title = element_text(size=10))
dev.off()
## cell cycle
pdf("output/seurat/VlnPlot_SCT_control_cYAPKO_cellCycle.pdf", width=15, height=10)
VlnPlot(embryo.combined.sct,features = c("S.Score","G2M.Score"), split.by = "condition") & 
  theme(plot.title = element_text(size=10))
dev.off()  





# Rename cluster
Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # cluster_1
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # cluster_2
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # cluster_3
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # cluster_4
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # cluster_5
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # cluster_6
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # cluster_7
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # cluster_8
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # cluster_9
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # cluster_10 (and 3)
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # cluster_11
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # cluster_12
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # cluster_13
Notocord = c("Shh","Noto","Vtn","Fgg") #  cluster_14
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # cluster_15
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3") # cluster_ 18
# Intermediate_Mesoderm = c("Bik","Arg1","Meox1","Ism1") # mix of 3 5 6 (let's not name it as more general than 3 5 6)



new.cluster.ids <- c(
  "Epiblast_PrimStreak", 
  "ExE_Ectoderm", 
  "Somitic_Mesoderm", 
  "Caudal_Mesoderm", 
  "Paraxial_Mesoderm", 
  "Pharyngeal_Mesoderm", 
  "Nascent_Mesoderm", 
  "Haematodenothelial_progenitors", 
  "Mesenchyme", 
  "Mixed_Mesoderm", 
  "Blood_Progenitor_1", 
  "Blood_Progenitor_2", 
  "Surface_Ectoderm", 
  "Notocord", 
  "Gut", 
  "Unknow_1",
  "Unknow_2", 
  "Primordial_Germ_Cells"
)

names(new.cluster.ids) <- levels(embryo.combined.sct)
embryo.combined.sct <- RenameIdents(embryo.combined.sct, new.cluster.ids)

embryo.combined.sct$cluster.annot <- Idents(embryo.combined.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_control_cYAPKO_label.pdf", width=12, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

#overlapping condition
pdf("output/seurat/UMAP_control_cYAPKO_label_overlap.pdf", width=6, height=5)
DimPlot(embryo.combined.sct, reduction = "umap", group.by = "condition", pt.size = 0.000001, cols = c("blue","red"))
dev.off()


# All in dotplot
DefaultAssay(embryo.combined.sct) <- "SCT"

Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # cluster_1
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # cluster_2
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # cluster_3
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # cluster_4
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # cluster_5
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # cluster_6
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # cluster_7
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # cluster_8
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # cluster_9
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # cluster_10 (and 3)
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # cluster_11
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # cluster_12
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # cluster_13
Notocord = c("Shh","Noto","Vtn","Fgg") #  cluster_14
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # cluster_15
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3")  # cluster_18


all_markers <- c(
  "Hesx1", "Otx2", "Pax2", "Otx2os1",
  "Tex19.1", "Elf5", "Wnt7b", "Dnmt3l",
  "Meox1", "Aldh1a2", "Synm", "Pcdh8",
  "Wnt3a", "Nkx1-2", "Apln", "Rxrg",
  "Ebf2", "Ptgfr", "Ebf3", "Col23a1",
  "Nkx2-5", "Tbx5", "Mef2c", "Myocd",
  "Col9a1", "Mesp1", "Osr1",
  "Hoxa9", "Hoxa10", "Tbx4", "Pitx1",
  "Tdo2", "Adamts2", "Colec11", "Snai2",
  "Tbx6", "Dll1", "Hes7", "Fgf17",
  "Sox18", "Esam", "Rhoj", "Flt4",
  "Gypa", "Gata1", "Cited4", "Gfi1b",
  "Tfap2a", "Npnt", "Wnt4", "Id4",
  "Shh", "Noto", "Vtn", "Fgg",
  "Pga5", "Hs3cst1", "Wfdc1", "Islr2",
  "Dppa3", "Sprr2a3", "Irf1", "Ifitm3"
)


levels(embryo.combined.sct) <- c("Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak")
pdf("output/seurat/DotPlot_SCT_control_cYAPKO_ident.pdf", width=16.5, height=4.5)
DotPlot(embryo.combined.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()



# count nb of cells in each cluster
control_clusters <- table(Idents(embryo.combined.sct)[embryo.combined.sct$condition == "WT"])
print(control_clusters)
cYAPKO_clusters <- table(Idents(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO"])
print(cYAPKO_clusters)

# statisctics to assess if distribution of cells among the clusters is significantly different between the two conditions
## Get the counts of cells in each cluster for each condition
control_clusters <- table(Idents(embryo.combined.sct)[embryo.combined.sct$condition == "WT"])
cYAPKO_clusters <- table(Idents(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO"])
## Combine the counts into a matrix
cluster_counts <- rbind(control_clusters, cYAPKO_clusters)
rownames(cluster_counts) <- c("WT", "cYAPKO")
## Perform the chi-square test
chisq_test_result <- chisq.test(cluster_counts)
## Print the result for overall difference
print(chisq_test_result)

# Statiscitcs for in-cluster difference (which cluster has significant different nb of cells):
total_untreated <- sum(control_clusters)
total_dasatinib <- sum(cYAPKO_clusters)

## Initialize a vector to store p-values
p_values <- numeric(length(control_clusters))
## Loop over each cluster and perform a chi-square test
for (i in seq_along(control_clusters)) {
  # Build contingency table for the i-th cluster
  contingency_table <- matrix(c(control_clusters[i], total_untreated - control_clusters[i], 
                                cYAPKO_clusters[i], total_dasatinib - cYAPKO_clusters[i]), 
                              nrow = 2)
  # Perform chi-square test
  chi_square_result <- chisq.test(contingency_table)
  # Store the p-value
  p_values[i] <- chi_square_result$p.value
}
## Perform Bonferroni correction
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
## Print adjusted p-values
names(adjusted_p_values) <- names(control_clusters)
print(adjusted_p_values)


## Downsampling with bootstrap to compare the nb of cell per cell types

library(tidyverse)

### Identify the unique clusters
unique_clusters <- unique(Idents(embryo.combined.sct))

### Create empty matrices to store cell counts
control_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
cYAPKO_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(control_clusters_counts) <- unique_clusters
colnames(cYAPKO_clusters_counts) <- unique_clusters

### Loop through 100 iterations
embryo.combined.sct_WT <- which(embryo.combined.sct$orig.ident == 'WT')
embryo.combined.sct_cYAPKO <- which(embryo.combined.sct$orig.ident == 'cYAPKO')

for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  embryo.combined.sct_WT_downsample <- sample(embryo.combined.sct_WT, 3621)
  embryo.combined.sct_integrated_downsample <- embryo.combined.sct[,c(embryo.combined.sct_cYAPKO, embryo.combined.sct_WT_downsample)]

  # Count nb of cells in each cluster
  control_clusters <- table(Idents(embryo.combined.sct_integrated_downsample)[embryo.combined.sct_integrated_downsample$condition == "WT"])
  cYAPKO_clusters <- table(Idents(embryo.combined.sct_integrated_downsample)[embryo.combined.sct_integrated_downsample$condition == "cYAPKO"])

  # Align the counts with the unique clusters
  control_clusters_counts[i, names(control_clusters)] <- as.numeric(control_clusters)
  cYAPKO_clusters_counts[i, names(cYAPKO_clusters)] <- as.numeric(cYAPKO_clusters)
}



### Calculate mean and standard error
mean_control_clusters <- colMeans(control_clusters_counts)
mean_cYAPKO_clusters <- colMeans(cYAPKO_clusters_counts)
std_error_WT_clusters <- apply(control_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("WT", "cYAPKO")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- control_clusters_counts[j,i]
    contingency_table[1,2] <- cYAPKO_clusters_counts[j,i]
    contingency_table[2,1] <- sum(control_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(cYAPKO_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_control_clusters),
  untreated = mean_control_clusters,
  dasatinib = mean_cYAPKO_clusters,
  std_error_WT = std_error_WT_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_WT, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "WT", "cYAPKO"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.01, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st
plot_data$cluster <- factor(plot_data$cluster, levels = c("Epiblast_PrimStreak",  # Early developmental stage
  "ExE_Ectoderm",         # Extraembryonic ectoderm; related to the epiblast
  "Surface_Ectoderm",     # Forms from ectoderm
  "Notocord",             # Important mesodermal structure
  "Nascent_Mesoderm",     # Early mesoderm
  "Paraxial_Mesoderm",    # Forms alongside notochord
  "Somitic_Mesoderm",     # Gives rise to somites
  "Pharyngeal_Mesoderm",  # Part of head mesoderm
  "Caudal_Mesoderm",      # Posterior mesoderm
  "Mixed_Mesoderm",       # Intermediate mesodermal stage
  "Mesenchyme",           # Loose, undifferentiated cells often derived from mesoderm
  "Haematodenothelial_progenitors", # Blood and endothelial progenitors
  "Blood_Progenitor_1",   # Blood precursors
  "Blood_Progenitor_2",   # Blood precursors
  "Gut",                  # Derived from endoderm
  "Primordial_Germ_Cells",# Early germ cell precursors
  "Unknow_1",             # Unknown categories placed at the end
  "Unknow_2")) # Reorder untreated 1st


# Plotting using ggplot2
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling10_clean_embryo.pdf", width=9, height=4)
ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = filter(plot_data, condition == "cYAPKO"),
    aes(label = significance, y = value + std_error_WT_clusters),
    vjust = -0.8,
    position = position_dodge(0.9), size = 5
  ) +
  scale_fill_manual(values = c("WT" = "blue", "cYAPKO" = "red")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,900)
dev.off()





# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

DefaultAssay(embryo.combined.sct) <- "RNA"

embryo.combined.sct <- NormalizeData(embryo.combined.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(embryo.combined.sct)
embryo.combined.sct <- ScaleData(embryo.combined.sct, features = all.genes) # zero-centres and scales it


## what genes change in different conditions for cells of the same type

embryo.combined.sct$celltype.stim <- paste(embryo.combined.sct$cluster.annot, embryo.combined.sct$condition,
    sep = "-")
Idents(embryo.combined.sct) <- "celltype.stim"

# use RNA corrected count for DEGs
## embryo.combined.sct <- PrepSCTFindMarkers(embryo.combined.sct)

## DEGs for each cell type:  
c("Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak")

Primordial_Germ_Cells.cYAPKO.response <- FindMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Primordial_Germ_Cells-cYAPKO", ident.2 = "Primordial_Germ_Cells-WT",
    verbose = FALSE)
head(Primordial_Germ_Cells.cYAPKO.response, n = 15)
write.table(Primordial_Germ_Cells.cYAPKO.response, file = "output/seurat/Primordial_Germ_Cells-cYAPKO_response_test.txt", sep = "\t", quote = FALSE, row.names = TRUE)


## Automation::
cell_types <- c("Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak")

for (cell_type in cell_types) {
  response_name <- paste(cell_type, "cYAPKO.response", sep = ".")
  ident_1 <- paste(cell_type, "-cYAPKO", sep = "")
  ident_2 <- paste(cell_type, "-WT", sep = "")

  response <- FindMarkers(embryo.combined.sct, assay = "RNA", ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
  
  print(head(response, n = 15))
  
  file_name <- paste("output/seurat/", cell_type, "-cYAPKO_response.txt", sep = "")
  write.table(response, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE)
}




## Display the top 10 marker genes of each cluster; this version display 10 to 20 genes; either the top10 are shared between condition
### Find all markers 
all_markers <- FindAllMarkers(embryo.combined.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_WT_cYAPKO_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# all_markers <- read.delim("output/seurat/srat_WT_cYAPKO_all_markers.txt", header = TRUE, row.names = 1)
##
top10 <- all_markers %>% 
  mutate(cluster = factor(cluster, levels = c(
  "Epiblast_PrimStreak-WT", "Epiblast_PrimStreak-cYAPKO",
  "ExE_Ectoderm-WT", "ExE_Ectoderm-cYAPKO",
  "Surface_Ectoderm-WT", "Surface_Ectoderm-cYAPKO",
  "Notocord-WT", "Notocord-cYAPKO",
  "Nascent_Mesoderm-WT", "Nascent_Mesoderm-cYAPKO",
  "Paraxial_Mesoderm-WT", "Paraxial_Mesoderm-cYAPKO",
  "Somitic_Mesoderm-WT", "Somitic_Mesoderm-cYAPKO",
  "Pharyngeal_Mesoderm-WT", "Pharyngeal_Mesoderm-cYAPKO",
  "Caudal_Mesoderm-WT", "Caudal_Mesoderm-cYAPKO",
  "Mixed_Mesoderm-WT", "Mixed_Mesoderm-cYAPKO",
  "Mesenchyme-WT", "Mesenchyme-cYAPKO",
  "Haematodenothelial_progenitors-WT", "Haematodenothelial_progenitors-cYAPKO",
  "Blood_Progenitor_1-WT", "Blood_Progenitor_1-cYAPKO",
  "Blood_Progenitor_2-WT", "Blood_Progenitor_2-cYAPKO",
  "Gut-WT", "Gut-cYAPKO",
  "Primordial_Germ_Cells-WT", "Primordial_Germ_Cells-cYAPKO",
  "Unknow_1-WT", "Unknow_1-cYAPKO",
  "Unknow_2-WT", "Unknow_2-cYAPKO"
))) %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>% 
  ungroup() %>% 
  arrange(match(cluster, c(
  "Epiblast_PrimStreak-WT", "Epiblast_PrimStreak-cYAPKO",
  "ExE_Ectoderm-WT", "ExE_Ectoderm-cYAPKO",
  "Surface_Ectoderm-WT", "Surface_Ectoderm-cYAPKO",
  "Notocord-WT", "Notocord-cYAPKO",
  "Nascent_Mesoderm-WT", "Nascent_Mesoderm-cYAPKO",
  "Paraxial_Mesoderm-WT", "Paraxial_Mesoderm-cYAPKO",
  "Somitic_Mesoderm-WT", "Somitic_Mesoderm-cYAPKO",
  "Pharyngeal_Mesoderm-WT", "Pharyngeal_Mesoderm-cYAPKO",
  "Caudal_Mesoderm-WT", "Caudal_Mesoderm-cYAPKO",
  "Mixed_Mesoderm-WT", "Mixed_Mesoderm-cYAPKO",
  "Mesenchyme-WT", "Mesenchyme-cYAPKO",
  "Haematodenothelial_progenitors-WT", "Haematodenothelial_progenitors-cYAPKO",
  "Blood_Progenitor_1-WT", "Blood_Progenitor_1-cYAPKO",
  "Blood_Progenitor_2-WT", "Blood_Progenitor_2-cYAPKO",
  "Gut-WT", "Gut-cYAPKO",
  "Primordial_Germ_Cells-WT", "Primordial_Germ_Cells-cYAPKO",
  "Unknow_1-WT", "Unknow_1-cYAPKO",
  "Unknow_2-WT", "Unknow_2-cYAPKO"
)))
# View the top 10 markers for each cluster
print(top10)
write.table(top10, file = "output/seurat/srat_WT_cYAPKO_top10_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# Visualize the top 10 markers for each cluster
marker_genes <- unique(top10$gene)
levels(embryo.combined.sct) <- c(
  "Unknow_2-cYAPKO", "Unknow_2-WT",
  "Unknow_1-cYAPKO", "Unknow_1-WT",
  "Primordial_Germ_Cells-cYAPKO", "Primordial_Germ_Cells-WT",
  "Gut-cYAPKO", "Gut-WT",
  "Blood_Progenitor_2-cYAPKO", "Blood_Progenitor_2-WT",
  "Blood_Progenitor_1-cYAPKO", "Blood_Progenitor_1-WT",
  "Haematodenothelial_progenitors-cYAPKO", "Haematodenothelial_progenitors-WT",
  "Mesenchyme-cYAPKO", "Mesenchyme-WT",
  "Mixed_Mesoderm-cYAPKO", "Mixed_Mesoderm-WT",
  "Caudal_Mesoderm-cYAPKO", "Caudal_Mesoderm-WT",
  "Pharyngeal_Mesoderm-cYAPKO", "Pharyngeal_Mesoderm-WT",
  "Somitic_Mesoderm-cYAPKO", "Somitic_Mesoderm-WT",
  "Paraxial_Mesoderm-cYAPKO", "Paraxial_Mesoderm-WT",
  "Nascent_Mesoderm-cYAPKO", "Nascent_Mesoderm-WT",
  "Notocord-cYAPKO", "Notocord-WT",
  "Surface_Ectoderm-cYAPKO", "Surface_Ectoderm-WT",
  "ExE_Ectoderm-cYAPKO", "ExE_Ectoderm-WT",
  "Epiblast_PrimStreak-cYAPKO", "Epiblast_PrimStreak-WT"
)
pdf("output/seurat/DotPlot_SCT_WT_cYAPKO_top10.pdf", width=50, height=10)
DotPlot(embryo.combined.sct, features = marker_genes, cols = c("grey", "red")) + RotatedAxis()
dev.off()


## display gene list in dotplot for both WT and cYPAKO
## !! may need to re-run FindAllMarkers and levels !!
gene_list= c("Tal1", "Gata1", "Vamp5", "Sox18", "Kdr")

pdf("output/seurat/DotPlot_SCT_WT_cYAPKO_geneList20230831.pdf", width=8, height=10)
DotPlot(embryo.combined.sct, features = gene_list, cols = c("grey", "red")) + RotatedAxis()
dev.off()

# Display the top 10 CONSERVED marker genes of each cluster
Idents(embryo.combined.sct) <- "cluster.annot"

## DEGs cluster versus all other
Primordial_Germ_Cells.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Primordial_Germ_Cells", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Primordial_Germ_Cells")
Unknow_2.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Unknow_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknow_2")
Unknow_1.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Unknow_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknow_1")
Gut.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Gut", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Gut")
Notocord.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Notocord", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Notocord")
Surface_Ectoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Surface_Ectoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Surface_Ectoderm")
Blood_Progenitor_2.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Blood_Progenitor_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_2")
Blood_Progenitor_1.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Blood_Progenitor_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_1")
Mixed_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Mixed_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mixed_Mesoderm")
Mesenchyme.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Mesenchyme", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mesenchyme")
Haematodenothelial_progenitors.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Haematodenothelial_progenitors", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Haematodenothelial_progenitors")
Nascent_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Nascent_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Nascent_Mesoderm")
Pharyngeal_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Pharyngeal_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Pharyngeal_Mesoderm")
Paraxial_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Paraxial_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Paraxial_Mesoderm")
Caudal_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Caudal_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Caudal_Mesoderm")
Somitic_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Somitic_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Somitic_Mesoderm")
ExE_Ectoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "ExE_Ectoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "ExE_Ectoderm")
Epiblast_PrimStreak.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Epiblast_PrimStreak", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Epiblast_PrimStreak")





## Combine all conserved markers into one data frame
all_conserved <- bind_rows(Primordial_Germ_Cells.conserved, Unknow_2.conserved, Unknow_1.conserved, Gut.conserved, Notocord.conserved, Surface_Ectoderm.conserved, Blood_Progenitor_2.conserved, Blood_Progenitor_1.conserved, Mixed_Mesoderm.conserved, Mesenchyme.conserved, Haematodenothelial_progenitors.conserved, Nascent_Mesoderm.conserved, Pharyngeal_Mesoderm.conserved, Paraxial_Mesoderm.conserved, Caudal_Mesoderm.conserved, Somitic_Mesoderm.conserved, ExE_Ectoderm.conserved, Epiblast_PrimStreak.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers_embryo.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak")))

## Find the top 3 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 3) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak")))



## Write the top 10 conserved markers for each cluster to a file
write.table(top10_conserved, file = "output/seurat/srat_top10_conserved_markers_embryo.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Visualize the top 10/3 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(embryo.combined.sct) <- c("Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak")
pdf("output/seurat/DotPlot_SCT_top10_conserved_embryo.pdf", width=35, height=5)
pdf("output/seurat/DotPlot_SCT_top3_conserved_embryo.pdf", width=18, height=5)

DotPlot(embryo.combined.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# save
## saveRDS(embryo.combined.sct, file = "output/seurat/embryo.combined.sct.rds")
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")


# Check some genes
DefaultAssay(embryo.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

## post 20231005 Conchi meeting
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Aldh1a2_Cyp26a1.pdf", width=10, height=13)
FeaturePlot(embryo.combined.sct, features = c("Aldh1a2", "Cyp26a1"), max.cutoff = 10, cols = c("grey", "red"), split.by = "condition")
dev.off()

## change cutoff res for highly express genes
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Birc5_Chchd2.pdf", width=10, height=15)
FeaturePlot(embryo.combined.sct, features = c("Birc5","Chchd2"), min.cutoff = 1, max.cutoff = 100, cols = c("grey", "red"), split.by = "condition")
dev.off()


## Loupe Browser increase in KO
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_increaseInLoupe.pdf", width=10, height=22)
FeaturePlot(embryo.combined.sct, features = c("Tal1", "Gata2", "Cited1", "Sox18"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()






## PRC2-related genes

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_PRC2.pdf", width=10, height=22)
FeaturePlot(embryo.combined.sct, features = c("Ezh2", "Kdm6b", "Gata6", "Cdx2","T"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

## Hand1
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_HAND1.pdf", width=10, height=5)
FeaturePlot(embryo.combined.sct, features = c("Hand1"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


### Check YAP1 hippo bulk-regulated genes:
DefaultAssay(humangastruloid.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_YAP1_hippo.pdf", width=7, height=30)

# CTGF, CRY61, AREG, MYC, GLI2, VIM, AXL and BIRC5

FeaturePlot(embryo.combined.sct, features = c("Tead1", "Tead4", "Ccn2", "Ccn1","Areg","Myc","Gli2","Vim","Axl","Birc5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()



### Check YAP1 NODAL bulk-regulated genes:

YAP1_NODAL <- c("Shh", "Dmrt1", "Cited2", 'Dact2', "Tgif2", "Acvr1b", 'Smad3', 'Dand5', 'Smad2','Cer1', 'Nodal', 'Foxh1', 'Dact1', 'Tdgf1', 'Cripto3', 'Acvr1c', 'Cfc1', 'Cfc1b') 

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_YAP1_NODAL.pdf", width=10, height=80)
FeaturePlot(embryo.combined.sct, features = YAP1_NODAL, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

### Check ANXA NODAL bulk-regulated genes:


pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_ANXA.pdf", width=10, height=35)
FeaturePlot(embryo.combined.sct, features = c("Anxa1","Amotl2","Ccn1","Ccn2","Tagln","Anxa3","Chchd2"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

## cell cycle


pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_cellCycle.pdf", width=10, height=5)
FeaturePlot(embryo.combined.sct,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()  



# Check marker genes
DefaultAssay(embryo.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA


pdf("output/seurat/FeaturePlot_SCT_control_test.pdf", width=10, height=30)
FeaturePlot(embryo.combined.sct, features = c("Eef1a1", "Rps29", "Tuba1b", "Tmsb10"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()  


pdf("output/seurat/FeaturePlot_SCT_control_test.pdf", width=10, height=30)
FeaturePlot(embryo.combined.sct, features = c("Hmgn2", "Rps29", "Hsp90aa1", "Tmsb10"), max.cutoff = 10, cols = c("grey", "red"), split.by = "condition")
dev.off()  

# test EasyCellType
# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_WT_cYAPKO_all_markers.txt", header = TRUE, row.names = 1)
### Filter either WT or cYAPKO
all_markers <- all_markers[grepl("WT$", all_markers$cluster), ]
all_markers <- all_markers[grepl("cYAPKO$", all_markers$cluster), ]

## Convert geneSymbol to EntrezID
all_markers$entrezid <- mapIds(org.Mm.eg.db,
                           keys=all_markers$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
all_markers <- na.omit(all_markers)

## Sort the datafram (data frame containing Entrez IDs, clusters and expression scores)

all_markers_sort <- data.frame(gene=all_markers$entrezid, cluster=all_markers$cluster, 
                      score=all_markers$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(all_markers_sort[, 1:3])

## Run the enrihcment analysis
annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or Panglaodb or Clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Embryo"), p_cut=0.3,   # many other tissue available
                    test="GSEA")    # GSEA or Fisher?


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or Panglaodb or Clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Embryo", "Embryoid body", "Embryonic stem cell", "Embryonic heart"), p_cut=0.3,   # many other tissue available
                    test="GSEA")    # GSEA or Fisher?


## plots
pdf("output/seurat/EasyCellType_dotplot_SCT_cYAPKO.pdf", width=10, height=5)
pdf("output/seurat/EasyCellType_dotplot_SCT_control.pdf", width=6, height=5)
pdf("output/seurat/EasyCellType_dotplot_SCT_control_noTissue.pdf", width=6, height=8)
pdf("output/seurat/EasyCellType_dotplot_SCT_control_clustermole.pdf", width=6, height=8)
pdf("output/seurat/EasyCellType_dotplot_SCT_control_allEmbryoTissue.pdf", width=6, height=5)
plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or Panglaodb or Clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Embryo"),   # many other tissue available
                    test="fisher")    # GSEA or fisher

pdf("output/seurat/EasyCellType_dotplot_SCT_control_fisher.pdf", width=6, height=8)
plot_dot(test="fisher", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



pdf("output/seurat/EasyCellType_barplot_SCT_control_cYAPKO.pdf", width=10, height=5)
plot_bar(test="GSEA", annot.GSEA)
dev.off()


# test singleR 

# Cell type annotation using SingleR
## Get reference datasets
hpca.ref <- celldex::HumanPrimaryCellAtlasData()


## Convert Seurat object into SCE
sce <- as.SingleCellExperiment(DietSeurat(embryo.combined.sct))
sce

## Analyse
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)


## check results

table(hpca.main$pruned.labels)


table(hpca.fine$pruned.labels)


# Add annotation to seurat object
embryo.combined.sct@meta.data$hpca.main   <- hpca.main$pruned.labels
embryo.combined.sct@meta.data$hpca.fine   <- hpca.fine$pruned.labels

# cluster with name
pdf("output/seurat/UMAP_SCT_control_cYAPKO_hpca.main.pdf", width=15, height=10)
embryo.combined.sct <- SetIdent(embryo.combined.sct, value = "hpca.main")
DimPlot(embryo.combined.sct, label = T , repel = T, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_SCT_control_cYAPKO_hpca.fine.pdf", width=35, height=30)
embryo.combined.sct <- SetIdent(embryo.combined.sct, value = "hpca.fine")
DimPlot(embryo.combined.sct, label = T , repel = T, label.size = 3)
dev.off()










# Functional analysis GO pathway
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")

library("ReactomeGSA")
library("ggrepel")
library("RColorBrewer")

DefaultAssay(embryo.combined.sct) <- "RNA" # For 

# separate condition from my seurat object
## Subset Seurat object based on condition
embryo.combined.sct.WT <- subset(embryo.combined.sct, subset = condition == "WT")
embryo.combined.sct.cYAPKO <- subset(embryo.combined.sct, subset = condition == "cYAPKO")

gsva_result <- analyse_sc_clusters(embryo.combined.sct.cYAPKO, verbose = TRUE)

gsva_result <- analyse_sc_clusters(embryo.combined.sct, verbose = TRUE)
pathway_expression <- pathways(gsva_result)
## maximum difference in expression for every pathway
### find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min
### sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]


## Save outputs
write.table(pathway_expression, "output/seurat/pathway_expression_control_cYAPKO.txt", sep="\t", row.names=TRUE, quote=FALSE)
write.table(max_difference, "output/seurat/max_difference_pathway_expression_control_cYAPKO.txt", sep="\t", row.names=TRUE, quote=FALSE)



## Plot
### Expression for a single pathway
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])

pdf("output/seurat/ReactomeGSA_SignalingRA_control_cYAPKO.pdf", width=15, height=10)
plot_gsva_pathway(gsva_result, pathway_id = "R-HSA-5362517")
dev.off()


pdf("output/seurat/ReactomeGSA_AlanineMet_control_cYAPKO.pdf", width=15, height=10)
plot_gsva_pathway(gsva_result, pathway_id = "R-HSA-8964540")
dev.off()
### Heatmap pathway


pdf("output/seurat/ReactomeGSA_heatmap_control.pdf", width=15, height=10)
pdf("output/seurat/ReactomeGSA_heatmap_cYAPKO.pdf", width=15, height=10)
pdf("output/seurat/ReactomeGSA_heatmap_control_cYAPKO.pdf", width=15, height=10)
plot_gsva_heatmap(gsva_result, max_pathways = 20, margins = c(12,40), truncate_names = FALSE, col = colorRampPalette(c("blue", "white", "red"))(100)) # ,   scale = "row"
dev.off()

pdf("output/seurat/ReactomeGSA_DevelopmentalBiology_heatmap_control_cYAPKO.pdf", width=15, height=10)
plot_gsva_heatmap(gsva_result, pathway_ids = c("R-HSA-9754189", "R-HSA-9758919", "R-HSA-9758920", "R-HSA-9761174", "R-HSA-9793380", "R-HSA-9796292", "R-HSA-9823730", "R-HSA-9733709", "R-HSA-525793", "R-HSA-9675108", "R-HSA-168256", "R-HSA-109582", "R-HSA-983231", "R-HSA-397014", "R-HSA-112316", "R-HSA-556833", "R-HSA-71387", "R-HSA-8963743", "R-HSA-453274", "R-HSA-453279", "R-HSA-6805567", "R-HSA-1181150", "R-HSA-5576891","R-HSA-8964540", "R-HSA-432030", "R-HSA-192456", "R-HSA-1433617", "R-HSA-916853", "R-HSA-380612", "R-HSA-5662702", "R-HSA-622323", "R-HSA-390651" ), margins = c(12,40), truncate_names = FALSE, col = colorRampPalette(c("blue", "white", "red"))(100),   scale = "row") # 
dev.off()


R-HSA-9754189.4 Germ layer formation at gastrulation
R-HSA-9758919.3 Epithelial-Mesenchymal Transition (EMT) during gastrulation
R-HSA-9758920.3 Formation of lateral plate mesoderm
R-HSA-9761174.1 Formation of intermediate mesoderm
R-HSA-9793380.3 Formation of paraxial mesoderm
R-HSA-9796292.3 Formation of axial mesoderm
R-HSA-9823730.2 Formation of definitive endoderm
R-HSA-9733709.1 Cardiogenesis
R-HSA-525793.3 Myogenesis
R-HSA-9675108.3 Nervous system development
R-HSA-168256  Immune System
R-HSA-109582 Hemostasis
R-HSA-983231 Factors involved in megakaryocyte development and platelet production
R-HSA-397014 Muscle contraction
R-HSA-112316 Neuronal System
R-HSA-556833 Metabolism of lipids
R-HSA-71387 Metabolism of carbohydrates
R-HSA-8963743 Digestion and absorption
R-HSA-453274 Mitotic G2-G2/M phases (proliferative)
R-HSA-453279 Mitotic G1 phase and G1/S transition
R-HSA-6805567 Keratinization
R-HSA-1181150 Signaling by NODAL
R-HSA-5576891 Cardiac conduction
	


R-HSA-8964540	Alanine metabolism
R-HSA-9758920	Formation of lateral plate mesoderm
R-HSA-432030	Transport of glycerol from adipocytes to the liver by Aquaporins
R-HSA-192456	Digestion of dietary lipid
R-HSA-1433617	Regulation of signaling by NODAL
R-HSA-916853	Degradation of GABA
R-HSA-380612	Metabolism of serotonin
R-HSA-5662702	Melanin biosynthesis
R-HSA-622323	Presynaptic nicotinic acetylcholine receptors
R-HSA-390651	Dopamine receptors

### Pathway-level PCA
pdf("output/seurat/ReactomeGSA_PCA_control.pdf", width=15, height=10)
pdf("output/seurat/ReactomeGSA_PCA_cYAPKO.pdf", width=15, height=10)
pdf("output/seurat/ReactomeGSA_PCA_control_cYAPKO.pdf", width=15, height=10)
plot_gsva_pca(gsva_result) +
  geom_text_repel(aes(label = sample), 
                   box.padding = 0.35, 
                   point.padding = 0.5, 
                   segment.color = 'grey50')
dev.off()


## Compare WT and cYAPKO using SCPA
### installation (in scRNAseqV1)
# devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
# devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
# devtools::install_github("jackbibby1/SCPA")
# --> Last fail so installed dependencies
# BiocManager::install(c("clustermole", "ComplexHeatmap", "singscore"))
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
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")

DefaultAssay(embryo.combined.sct) <- "RNA" # Recommended 


# Isolate WT and cYAPKO
## These do not work as for human
pathways <- "output/Pathway/h_k_r_go_pid_reg_wik.csv" ## combination of canonical pathways, gene ontology pathways, and regulatory pathways from MSigDB
pathways <- "output/Pathway/combined_metabolic_pathways.csv" ## 

## These for mice dl from msigdb
#### if single file
pathways <- get_paths("output/Pathway/REACTOME_SIGNALING_BY_RETINOIC_ACID.v2023.1.Mm.gmt")
names(pathways) <- sapply(pathways, function(x) x$Pathway[1])
pathways$REACTOME_SIGNALING_BY_RETINOIC_ACID
#### if multiple files
pathways <- msigdbr("Mus musculus", "C2") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise
pathways$REACTOME_SIGNALING_BY_RETINOIC_ACID$Genes

# Compare pathway activity between condition within  
cell_types <- unique(embryo.combined.sct$cluster.annot)
embryo.combined.sct_split <- SplitObject(embryo.combined.sct, split.by = "condition")

# SCPA comparison
## keep all columns
scpa_out_all <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryo.combined.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryo.combined.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}
## keep only PAthway and qval column (Best for heatmap qval representation)
scpa_out <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryo.combined.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryo.combined.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}
saveRDS(scpa_out, file = "output/Pathway/scpa_out_eachCellTypes.rds")

console_output <- capture.output(print(scpa_out))
writeLines(console_output, "output/Pathway/scpa_out_console.txt")

## Save output as R object
# saveRDS(scpa_out, file = "output/Pathway/scpa_out_eachCellTypes.rds")
# scpa_out_all <- readRDS("output/Pathway/scpa_out_all_eachCellTypes.rds")

# plot output
## Filter pathways with a qval of > 2 in any comparison
scpa_out_filter <- scpa_out %>% 
  purrr::reduce(full_join, by = "Pathway") %>% 
  set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
  select(c("Pathway", grep("_qval", colnames(.)))) %>%
  filter_all(any_vars(. > 2)) %>%
  column_to_rownames("Pathway")
## pathway to highlight
highlight_paths <- c("REACTOME_SIGNALING_BY_RETINOIC_ACID")
highlight_paths <- c("REACTOME_SIGNALING_BY_RETINOIC_ACID","KEGG_RETINOL_METABOLISM","PID_RETINOIC_ACID_PATHWAY","WP_4HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION")

position <- which(rownames(scpa_out_filter) %in% highlight_paths)
row_an <- rowAnnotation(Genes = anno_mark(at = which(rownames(scpa_out_filter) %in% highlight_paths),
                                          labels = rownames(scpa_out_filter)[position],
                                          labels_gp = gpar(fontsize = 7),
                                          link_width = unit(2.5, "mm"),
                                          padding = unit(1, "mm"),
                                          link_gp = gpar(lwd = 0.5)))
col_hm <- colorRamp2(colors = c("blue", "white", "red"), breaks = c(0, 3, 6))

pdf("output/Pathway/SCPA_embryo_msigdb_C2_cl5.pdf", width=10, height=30)
pdf("output/Pathway/SCPA_embryo_retino_msigdb_C2_cl5.pdf", width=10, height=30)
Heatmap(scpa_out_filter,
        col = col_hm,
        name = "Qval",
        show_row_names = F,
        right_annotation = row_an,
        column_names_gp = gpar(fontsize = 8),
        border = T,
        column_km = 5,
        row_km = 5)
dev.off()

# Heatmap representation
scpa_out_tibble <- scpa_out_filter %>% 
  rownames_to_column(var = "Pathway") %>%
  as_tibble() %>%
  pivot_longer(cols = -Pathway, 
               names_to = "cluster", 
               values_to = "qval")


scpa_out_tibble_REACTOME_SIGNALING_BY_RETINOIC_ACID = scpa_out_tibble %>%
  filter(Pathway == "REACTOME_SIGNALING_BY_RETINOIC_ACID")
scpa_out_tibble_REACTOME_SIGNALING_BY_RETINOIC_ACID = scpa_out_tibble %>%
  filter(Pathway %in% c("REACTOME_SIGNALING_BY_RETINOIC_ACID","KEGG_RETINOL_METABOLISM","PID_RETINOIC_ACID_PATHWAY","WP_4HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION"))

scpa_out_all_tibble_REACTOME_SIGNALING_BY_RETINOIC_ACID = scpa_out_all %>% 
  purrr::reduce(full_join, by = "Pathway") %>% 
  set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
  select(c("Pathway", grep("_qval", colnames(.)))) %>%
  filter_all(any_vars(. > 2)) %>%
  column_to_rownames("Pathway") %>% 
  rownames_to_column(var = "Pathway") %>%
  as_tibble() %>%
  pivot_longer(cols = -Pathway, 
               names_to = "cluster", 
               values_to = "qval") %>%
  filter(Pathway == "REACTOME_SIGNALING_BY_RETINOIC_ACID")
# --> Negative FC = MORE ACTIVE in cYAPKO

pdf("output/Pathway/heatmap_embryo_msigdb_REACTOME_SIGNALING_BY_RETINOIC_ACID.pdf", width=5, height=5)
pdf("output/Pathway/heatmap_embryo_msigdb_retino.pdf", width=10, height=10)
ggplot(scpa_out_tibble_REACTOME_SIGNALING_BY_RETINOIC_ACID, aes(x=cluster, y=Pathway, fill=qval)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=1.4, name="qval") + # 1.4 = qval 0.05
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()


# Per Cell-type/Cluster comparison; on the signficiant
## c("Epiblast_PrimStreak", "Paraxial_Mesoderm", "Somitic_Mesoderm","Caudal_Mesoderm", "Nascent_Mesoderm", "Pharyngeal_Mesoderm", Mesenchyme")
WT <- seurat_extract(embryo.combined.sct,
                          meta1 = "condition", value_meta1 = "WT",
                          meta2 = "cluster.annot", value_meta2 = "Paraxial_Mesoderm")
cYAPKO <- seurat_extract(embryo.combined.sct,
                            meta1 = "condition", value_meta1 = "cYAPKO",
                            meta2 = "cluster.annot", value_meta2 = "Paraxial_Mesoderm")


WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),   # list(population1,population2) FC = population 2 - population 1
                             pathways = pathways,           # so FC >1 = less active in cYAPKO
                             parallel = TRUE, cores = 8)


WT_cYAPKO_filter <- WT_cYAPKO %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))
# plot single apthway
ra_path <- WT_cYAPKO_filter %>% 
  filter(grepl(pattern = "REACTOME_SIGNALING_BY_RETINOIC_ACID", ignore.case = T, x = Pathway))
ra_path <- WT_cYAPKO_filter %>% 
  filter(grepl(pattern = "PID_RETINOIC_ACID_PATHWAY", ignore.case = T, x = Pathway))
ra_path <- WT_cYAPKO_filter %>% 
  filter(grepl(pattern = "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE", ignore.case = T, x = Pathway))

pdf("output/Pathway/plot_embryo_msigdb_REACTOME_SIGNALING_BY_RETINOIC_ACID_Epiblast_PrimStreak.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_REACTOME_SIGNALING_BY_RETINOIC_ACID_Paraxial_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_REACTOME_SIGNALING_BY_RETINOIC_ACID_Somitic_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_PID_RETINOIC_ACID_PATHWAY_Caudal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_PID_RETINOIC_ACID_PATHWAY_Nascent_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_PID_RETINOIC_ACID_PATHWAY_Pharyngeal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_Mesenchyme.pdf", width=5, height=5)

ggplot(WT_cYAPKO_filter, aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = WT_cYAPKO_filter$color, stroke = 0.3) +
  geom_point(data = ra_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
dev.off()


pdf("output/Pathway/plotrank_embryo_msigdb_REACTOME_SIGNALING_BY_RETINOIC_ACID_Epiblast_PrimStreak.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_REACTOME_SIGNALING_BY_RETINOIC_ACID_Paraxial_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_REACTOME_SIGNALING_BY_RETINOIC_ACID_Somitic_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_PID_RETINOIC_ACID_PATHWAY_Caudal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_PID_RETINOIC_ACID_PATHWAY_Nascent_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_PID_RETINOIC_ACID_PATHWAY_Pharyngeal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_Mesenchyme.pdf", width=5, height=5)

plot_rank(WT_cYAPKO_filter, "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
                highlight_point_size = 3.5, highlight_point_color = "orangered2", label_pathway = FALSE)
dev.off()
# plot multiple apthway
ra_path <- WT_cYAPKO_filter %>% 
  filter(Pathway %in% c("REACTOME_SIGNALING_BY_RETINOIC_ACID","KEGG_RETINOL_METABOLISM","PID_RETINOIC_ACID_PATHWAY","WP_4HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION"))

pdf("output/Pathway/plotrank_embryo_msigdb_retino_Epiblast_PrimStreak.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_retino_Paraxial_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_retino_Somitic_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_retino_Caudal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_retino_Nascent_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plotrank_embryo_msigdb_retino_Pharyngeal_Mesoderm.pdf", width=5, height=5)

plot_rank(WT_cYAPKO_filter, c("REACTOME_SIGNALING_BY_RETINOIC_ACID","KEGG_RETINOL_METABOLISM","PID_RETINOIC_ACID_PATHWAY","WP_4HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION"),
                highlight_point_size = 3.5, highlight_point_color = "orangered2", label_pathway = TRUE)
dev.off()

pdf("output/Pathway/plot_embryo_msigdb_retino_Epiblast_PrimStreak.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_retino_Paraxial_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_retino_Somitic_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_retino_Caudal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_retino_Nascent_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_msigdb_retino_Pharyngeal_Mesoderm.pdf", width=5, height=5)

ggplot(WT_cYAPKO_filter, aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = WT_cYAPKO_filter$color, stroke = 0.3) +
  geom_point(data = ra_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
dev.off()



# Code to save output for each cell type comparison
clusters = c(
"Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(embryo.combined.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  cYAPKO <- seurat_extract(embryo.combined.sct,
                           meta1 = "condition", value_meta1 = "cYAPKO",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_", cluster, ".txt")
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Display non enriched but qvalue significant pathways
## heatmap
### extract Pathway gene name
pathways$REACTOME_SIGNALING_BY_RETINOIC_ACID$Genes
pathways$REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE$Genes
pathways$REACTOME_SIGNALING_BY_WNT$Genes


genes_of_interest <- pathways$REACTOME_SIGNALING_BY_WNT$Genes
### extract WT and cYAPKO gene expression values from RNA assay
WT_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "WT" & embryo.combined.sct$cluster.annot == "Paraxial_Mesoderm"]]
cYAPKO_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO" & embryo.combined.sct$cluster.annot == "Paraxial_Mesoderm"]]

### mean expression values for each gene
WT_mean <- rowMeans(WT_expression)
cYAPKO_mean <- rowMeans(cYAPKO_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  WT = WT_mean,
  cYAPKO = cYAPKO_mean
) %>% 
pivot_longer(cols = c(WT, cYAPKO), names_to = "Condition", values_to = "Expression")

## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(WT_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)


pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_RETINOIC_ACID_Paraxial_Mesoderm.pdf", width=10, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_Mesenchyme.pdf", width=20, height=3)

ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()


## CDF plot representation
### Additional filtering removing the gene not express in both genotypes
WT_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "WT" & embryo.combined.sct$cluster.annot == "Mesenchyme"]]
cYAPKO_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO" & embryo.combined.sct$cluster.annot == "Mesenchyme"]]


### Identify and remove the genes not expressing in WT and cYAPKO; NOT MANDATORY !!!
WT_tibble <- WT_expression %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  mutate(Sum = rowSums(select(., -Gene)),
         Expressed_WT = ifelse(Sum > 0, "yes", "no"))

cYAPKO_tibble <- cYAPKO_expression %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  mutate(Sum = rowSums(select(., -Gene)),
         Expressed_cYAPKO = ifelse(Sum > 0, "yes", "no"))
# Combine the tibbles
combined_tibble <- WT_tibble %>%
  select(Gene, Expressed_WT) %>%
  left_join(cYAPKO_tibble %>% select(Gene, Expressed_cYAPKO), by = "Gene") %>%
  replace_na(list(Expressed_WT = "no", Expressed_cYAPKO = "no")) %>%
  mutate(Expressed = ifelse(Expressed_WT == "yes" | Expressed_cYAPKO == "yes", "yes", "no"))
genes_to_keep_tibble <- combined_tibble %>% filter(Expressed == "yes")
# Extract the list of genes to keep
genes_to_keep <- genes_to_keep_tibble$Gene

WT_expression_filtered <- WT_expression[genes_to_keep, ]
cYAPKO_expression_filtered <- cYAPKO_expression[genes_to_keep, ]

#### Pick either keep zero or not
WT_values <- as.vector(WT_expression)
cYAPKO_values <- as.vector(cYAPKO_expression)

WT_values <- as.vector(WT_expression_filtered)
cYAPKO_values <- as.vector(cYAPKO_expression_filtered)
####


data_for_cdf <- data.frame(
  Expression = c(WT_values, cYAPKO_values),
  Condition = c(rep("WT", length(WT_values)), rep("cYAPKO", length(cYAPKO_values)))
)

# Filtering the data
# data_for_cdf <- data_for_cdf %>% filter(Expression > 0)
# Conduct the KS test
ks_result <- ks.test(data_for_cdf$Expression[data_for_cdf$Condition == "WT"],
                    data_for_cdf$Expression[data_for_cdf$Condition == "cYAPKO"])
# Extract p-value and statistic from KS test
p_value <- ks_result$p.value
statistic <- ks_result$statistic
# Plotting
pdf("output/Pathway/cdf_REACTOME_SIGNALING_BY_RETINOIC_ACID_Paraxial_Mesoderm.pdf", width=3, height=3)
pdf("output/Pathway/cdf_REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE_Mesenchyme.pdf", width=3, height=3)
ggplot(data_for_cdf, aes(x=Expression, color=Condition)) +
  stat_ecdf(geom = "step") +
  theme_bw() +
  labs(x="norm counts", y="Cumulative Proportion") +
  theme(legend.position = "bottom") +
  annotate("text", x = max(data_for_cdf$Expression) * 0.5, y = 0.2, 
           label = paste0("D = ", round(statistic, 5), "\np-value = ", format(p_value, digits=3, scientific=TRUE)), hjust = 0)
dev.off()




# BALOW CLEAN CODE; dotplot>enrichplot>heatmap>violinplot
# DOT PLOT pepresetnation
## load all the comparison for each cell type (FC qval information)
clusters = c(
"Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak"
)
## One by one
SCPA_Primordial_Germ_Cells <- read.delim("output/Pathway/SCPA_Primordial_Germ_Cells.txt", header = TRUE) %>%
  add_column(cluster = "Primordial_Germ_Cells")
Unknow_2 <- read.delim("output/Pathway/SCPA_Unknow_2.txt", header = TRUE) %>%
  add_column(cluster = "Unknow_2")
## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster))
# --> import the gene pathways (used table from Conchi `Pathwyas of interest*.xlsx`)
pathways = c(
"BIOCARTA_WNT_PATHWAY",
"KEGG_WNT_SIGNALING_PATHWAY",
"WP_WNT_SIGNALING",
"WP_WNT_SIGNALING_PATHWAY",
"WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY",
"PID_WNT_CANONICAL_PATHWAY",
"PID_WNT_NONCANONICAL_PATHWAY",
"PID_WNT_SIGNALING_PATHWAY",
"WILLERT_WNT_SIGNALING",
"WNT_SIGNALING",
"REACTOME_SIGNALING_BY_WNT",
"REACTOME_WNT5A_DEPENDENT_INTERNALIZATION_OF_FZD4",
"REACTOME_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING",
"REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING",
"REACTOME_NEGATIVE_REGULATION_OF_TCF_DEPENDENT_SIGNALING_BY_WNT_LIGAND_ANTAGONISTS",
"REACTOME_DEACTIVATION_OF_THE_BETA_CATENIN_TRANSACTIVATING_COMPLEX",
"REACTOME_FORMATION_OF_THE_BETA_CATENIN_TCF_TRANSACTIVATING_COMPLEX",
"REACTOME_DEGRADATION_OF_BETA_CATENIN_BY_THE_DESTRUCTION_COMPLEX",
"REACTOME_DEGRADATION_OF_AXIN",
"REACTOME_BETA_CATENIN_PHOSPHORYLATION_CASCADE",
"REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT",
"REACTOME_NEGATIVE_REGULATION_OF_TCF_DEPENDENT_SIGNALING_BY_WNT_LIGAND_ANTAGONISTS",
"REACTOME_WNT5A_DEPENDENT_INTERNALIZATION_OF_FZD4"
)
pathways = c(
"WP_TGFBETA_RECEPTOR_SIGNALING",
"WP_TGFBETA_SIGNALING_PATHWAY",
"KEGG_TGF_BETA_SIGNALING_PATHWAY",
"BIOCARTA_TGFB_PATHWAY",
"KARAKAS_TGFB1_SIGNALING",
"PID_TGFBR_PATHWAY",
"REACTOME_SIGNALING_BY_NODAL",
"REACTOME_DOWNREGULATION_OF_SMAD2_3_SMAD4_TRANSCRIPTIONAL_ACTIVITY",
"REACTOME_TRANSCRIPTIONAL_ACTIVITY_OF_SMAD2_SMAD3_SMAD4_HETEROTRIMER",
"REACTOME_SMAD2_SMAD3_SMAD4_HETEROTRIMER_REGULATES_TRANSCRIPTION",
"REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
"REACTOME_SIGNALING_BY_ACTIVIN",
"REACTOME_DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING",
"REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS",
"REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",
"REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
"WP_CANONICAL_AND_NONCANONICAL_TGFB_SIGNALING"
)
pathways = c(
"REACTOME_SIGNALING_BY_HIPPO",
"REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION",
"WP_HIPPO_SIGNALING_REGULATION_PATHWAYS",
"WP_MECHANOREGULATION_AND_PATHOLOGY_OF_YAPTAZ_VIA_HIPPO_AND_NONHIPPO_MECHANISMS",
"WP_HIPPOYAP_SIGNALING_PATHWAY",
"REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION",
"WP_MECHANOREGULATION_AND_PATHOLOGY_OF_YAPTAZ_VIA_HIPPO_AND_NONHIPPO_MECHANISMS"
)
pathways = c(
"REACTOME_SIGNALING_BY_RETINOIC_ACID",
"REACTOME_THE_CANONICAL_RETINOID_CYCLE_IN_RODS_TWILIGHT_VISION",
"KEGG_RETINOL_METABOLISM",
"PID_RETINOIC_ACID_PATHWAY",
"WP_4HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION"
)
pathways = c(
"REACTOME_SIGNALING_BY_NOTCH",
"REACTOME_SIGNALING_BY_NOTCH1",
"REACTOME_SIGNALING_BY_NOTCH2",
"REACTOME_SIGNALING_BY_NOTCH3",
"REACTOME_SIGNALING_BY_NOTCH4",
"REACTOME_NOTCH_HLH_TRANSCRIPTION_PATHWAY" ,
"REACTOME_PRE_NOTCH_PROCESSING_IN_GOLGI",
"REACTOME_PRE_NOTCH_EXPRESSION_AND_PROCESSING",
"REACTOME_NOTCH2_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",
"REACTOME_NOTCH3_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",
"REACTOME_ACTIVATED_NOTCH1_TRANSMITS_SIGNAL_TO_THE_NUCLEUS",
"REACTOME_NOTCH3_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION",
"REACTOME_NOTCH4_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION",
"KEGG_NOTCH_SIGNALING_PATHWAY" ,
"PID_NOTCH_PATHWAY" ,
"WP_NOTCH_SIGNALING" ,
"WP_NOTCH_SIGNALING_PATHWAY" ,
"WP_CANONICAL_AND_NONCANONICAL_NOTCH_SIGNALING"
)
pathways = c(
"REACTOME_SIGNALING_BY_BMP",
"WP_BMP_SIGNALING_IN_EYELID_DEVELOPMENT",
"PID_BMP_PATHWAY"
)
pathways = c(
"REACTOME_HDMS_DEMETHYLATE_HISTONES",
"REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA",
"REACTOME_HDACS_DEACETYLATE_HISTONES",
"REACTOME_RMTS_METHYLATE_HISTONE_ARGININES",
"REACTOME_HATS_ACETYLATE_HISTONES",
"REACTOME_PKMTS_METHYLATE_HISTONE_LYSINES",
"WP_INTERACTOME_OF_POLYCOMB_REPRESSIVE_COMPLEX_2_PRC2",
"BIOCARTA_HDAC_PATHWAY"
)
pathways = c(
"WP_ENDODERM_DIFFERENTIATION",
"WP_MESODERMAL_COMMITMENT_PATHWAY",
"WP_ECTODERM_DIFFERENTIATION"
)
pathways = c(
"REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_DOPAMINE_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_SEROTONIN_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE",
"REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
"WP_COMPLEMENT_SYSTEM_IN_NEURONAL_DEVELOPMENT_AND_PLASTICITY",
"REACTOME_MECP2_REGULATES_NEURONAL_RECEPTORS_AND_CHANNELS",
"REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS",
"WP_MAJOR_RECEPTORS_TARGETED_BY_EPINEPHRINE_AND_NOREPINEPHRINE",
"REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
"REACTOME_NEUREXINS_AND_NEUROLIGINS",
"REACTOME_ASSEMBLY_AND_CELL_SURFACE_PRESENTATION_OF_NMDA_RECEPTORS",
"REACTOME_NEGATIVE_REGULATION_OF_NMDA_RECEPTOR_MEDIATED_NEURONAL_TRANSMISSION",
"REACTOME_ACTIVATION_OF_NMDA_RECEPTORS_AND_POSTSYNAPTIC_EVENTS",
"KEGG_NEUROTROPHIN_SIGNALING_PATHWAY",
"REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES"
)
pathways = c(
"REACTOME_BLOOD_GROUP_SYSTEMS_BIOSYNTHESIS",
"REACTOME_HEME_SIGNALING",
"REACTOME_SCAVENGING_OF_HEME_FROM_PLASMA",
"REACTOME_RESPONSE_OF_EIF2AK1_HRI_TO_HEME_DEFICIENCY",
"KEGG_VEGF_SIGNALING_PATHWAY",
"REACTOME_SIGNALING_BY_VEGF",
"REACTOME_VEGFR2_MEDIATED_VASCULAR_PERMEABILITY",
"REACTOME_VEGFR2_MEDIATED_CELL_PROLIFERATION",
"BIOCARTA_VEGF_PATHWAY",
"WP_VEGFAVEGFR2_SIGNALING_PATHWAY",
"PID_VEGFR1_PATHWAY",
"PID_VEGFR1_2_PATHWAY"
)
pathways = c(
"REACTOME_SIGNALING_BY_FGFR",
"REACTOME_SIGNALING_BY_FGFR1",
"REACTOME_FGFR1_LIGAND_BINDING_AND_ACTIVATION",
"REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR1",
"REACTOME_SHC_MEDIATED_CASCADE_FGFR1",
"REACTOME_FRS_MEDIATED_FGFR1_SIGNALING",
"REACTOME_PI_3K_CASCADE_FGFR1",
"REACTOME_NEGATIVE_REGULATION_OF_FGFR1_SIGNALING",
"REACTOME_SIGNALING_BY_FGFR2",
"REACTOME_FGFR2_LIGAND_BINDING_AND_ACTIVATION",
"REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR2",
"REACTOME_SHC_MEDIATED_CASCADE_FGFR2",
"REACTOME_FRS_MEDIATED_FGFR2_SIGNALING",
"REACTOME_PI_3K_CASCADE_FGFR2",
"REACTOME_PHOSPHOLIPASE_C_MEDIATED_CASCADE_FGFR2",
"REACTOME_NEGATIVE_REGULATION_OF_FGFR2_SIGNALING",
"REACTOME_FGFR2_ALTERNATIVE_SPLICING",
"REACTOME_SIGNALING_BY_FGFR3",
"REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR3",
"REACTOME_SHC_MEDIATED_CASCADE_FGFR3",
"REACTOME_FRS_MEDIATED_FGFR3_SIGNALING",
"REACTOME_PI_3K_CASCADE_FGFR3",
"REACTOME_NEGATIVE_REGULATION_OF_FGFR3_SIGNALING",
"REACTOME_SIGNALING_BY_FGFR4",
"REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR4",
"REACTOME_SHC_MEDIATED_CASCADE_FGFR4",
"REACTOME_FRS_MEDIATED_FGFR4_SIGNALING",
"REACTOME_PI_3K_CASCADE_FGFR4",
"REACTOME_PHOSPHOLIPASE_C_MEDIATED_CASCADE_FGFR4",
"REACTOME_NEGATIVE_REGULATION_OF_FGFR4_SIGNALING"
)
all_data_pathways = as_tibble(all_data) %>% 
  filter(Pathway %in% pathways)

# tidy the df
all_data_pathways_tidy <- all_data_pathways %>%
  filter(qval >= 1.4) 

# Dot plot
pdf("output/Pathway/dotplot_WNT_signaling.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_NODAL_TGFB_signaling.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_HIPPO_signaling.pdf", width=12, height=3)
pdf("output/Pathway/dotplot_RA_signaling.pdf", width=12, height=2)
pdf("output/Pathway/dotplot_Notch_signaling.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_BMP_signaling.pdf", width=8, height=2)
pdf("output/Pathway/dotplot_Epigenetics.pdf", width=10, height=3)
pdf("output/Pathway/dotplot_EndoMesoEcto.pdf", width=10, height=3)
pdf("output/Pathway/dotplot_Neuro.pdf", width=12, height=4)
pdf("output/Pathway/dotplot_Blood.pdf", width=12, height=4)
pdf("output/Pathway/dotplot_FGF.pdf", width=10, height=6)

ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, color = -FC), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8), guide = "legend") +  # Adjust the range for size if needed
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, guide = "colourbar") +
  theme_bw() +
  labs(size = "q-value", color = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# REFINE COLOR
custom_color <- function(fc_value){
  ifelse(fc_value >= 5, "strong_blue", 
         ifelse(fc_value > 2, "light_blue",
                ifelse(fc_value <= -5, "strong_red", 
                       ifelse(fc_value < -2, "light_red", "grey"))))
}

# Add a column for this custom color
all_data_pathways_tidy <- all_data_pathways_tidy %>%
  mutate(custom_col = sapply(FC, custom_color))

pdf("output/Pathway/dotplot_WNT_signaling_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_NODAL_TGFB_signaling_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_HIPPO_signaling_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_RA_signaling_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_Notch_signaling_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_BMP_signaling_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_Epigenetics_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_EndoMesoEcto_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_Neuro_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_Blood_FCtresh.pdf", width=12, height=6)
pdf("output/Pathway/dotplot_FGF_FCtresh.pdf", width=12, height=6)


ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, color = custom_col), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = c("strong_blue" = "dodgerblue4",
                                "light_blue" = "lightblue2",
                                "grey" = "grey",
                                "light_red" = "indianred1",
                                "strong_red" = "red3")) +
  theme_bw() +
  labs(size = "q-value", color = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# ENRCIHMENT PLOT

all_data_pathways_tidy_pa <- all_data_pathways_tidy %>% 
  filter(Pathway %in% c("REACTOME_SIGNALING_BY_FGFR"),
         cluster == "Blood_Progenitor_1") %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))


all_data_all_pa = as_tibble(all_data) %>% 
  filter(cluster == "Blood_Progenitor_1") %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))


pdf("output/Pathway/plot_embryo_WNT_signaling_Paraxial_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_TGFB_signaling_Nascent_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_REACTOME_SIGNALING_BY_NODAL_Pharyngeal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_REACTOME_SIGNALING_BY_HIPPO_Mesenchyme.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_PID_RETINOIC_ACID_PATHWAY_Pharyngeal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_REACTOME_SIGNALING_BY_NOTCH_Paraxial_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_WP_BMP_SIGNALING_IN_EYELID_DEVELOPMENT_Paraxial_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_REACTOME_HATS_ACETYLATE_HISTONES_Pharyngeal_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_WP_MESODERMAL_COMMITMENT_PATHWAY_Mixed_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION_Nascent_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_WP_VEGFAVEGFR2_SIGNALING_PATHWAY_Nascent_Mesoderm.pdf", width=5, height=5)
pdf("output/Pathway/plot_embryo_REACTOME_SIGNALING_BY_FGFR_Blood_Progenitor_1.pdf", width=5, height=5)


ggplot(all_data_all_pa, aes(x = -FC, y = qval, fill = -FC)) +
  geom_point(aes(color = -FC), cex = 2, shape = 21, stroke = 0.3) +  # Color based on -FC
  geom_point(data = all_data_pathways_tidy_pa, aes(x = -FC, y = qval), shape = 21, cex = 2.8, fill = "forestgreen", color = "black", stroke = 0.3) +
 # geom_label_repel(data = all_data_pathways_tidy_pa, aes(label = Pathway, x = -FC, y = qval), box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50') +
  xlim(-25, 25) +
  ylim(0, 7.5) +
  xlab("Enrichment") +
  ylab("qval") +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, guide = "colourbar") +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, guide = "colourbar") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1) +
  geom_vline(xintercept = 0, linetype = "solid", col = 'black', lwd = 0.3) +
  geom_hline(yintercept = 1.4, linetype = "dashed", col = 'black', lwd = 0.3)
dev.off()


# heatmap
pathways_heatmap <- msigdbr("Mus musculus", "C2") %>%
format_pathways()
names(pathways_heatmap) <- sapply(pathways_heatmap, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

##pick pathway to representL
pathways_heatmap$REACTOME_SIGNALING_BY_FGFR$Genes
genes_of_interest <- pathways_heatmap$REACTOME_SIGNALING_BY_FGFR$Genes
##

genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(embryo.combined.sct@assays$RNA@data)]

### extract WT and cYAPKO gene expression values from RNA assay
WT_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "WT" & embryo.combined.sct$cluster.annot == "Blood_Progenitor_1"]]
cYAPKO_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO" & embryo.combined.sct$cluster.annot == "Blood_Progenitor_1"]]

### mean expression values for each gene
WT_mean <- rowMeans(WT_expression)
cYAPKO_mean <- rowMeans(cYAPKO_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  WT = WT_mean,
  cYAPKO = cYAPKO_mean
) %>% 
pivot_longer(cols = c(WT, cYAPKO), names_to = "Condition", values_to = "Expression")

## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(WT_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)

# if few genes:
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_WNT_Paraxial_Mesoderm.pdf", width=10, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS_Nascent_Mesoderm.pdf", width=10, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_NODAL_Pharyngeal_Mesoderm.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_HIPPO_Mesenchyme.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_PID_RETINOIC_ACID_PATHWAY_Pharyngeal_Mesoderm.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_BMP_SIGNALING_IN_EYELID_DEVELOPMENT_Paraxial_Mesoderm.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_HATS_ACETYLATE_HISTONES_Pharyngeal_Mesoderm.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_MESODERMAL_COMMITMENT_PATHWAY_Mixed_Mesoderm.pdf", width=6, height=3)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_FGFR_Blood_Progenitor_1.pdf", width=6, height=3)

library("viridis")
ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()

# if many genesL
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_WNT_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_NOTCH_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_HATS_ACETYLATE_HISTONES_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_MESODERMAL_COMMITMENT_PATHWAY_Mixed_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_WP_VEGFAVEGFR2_SIGNALING_PATHWAY_Nascent_Mesoderm.pdf", width=5, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_REACTOME_SIGNALING_BY_FGFR_Blood_Progenitor_1.pdf", width=3, height=2)

ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") 
dev.off()

# VIOLIN if enrichment
data_for_plot$Condition <-
  factor(data_for_plot$Condition,
         c("WT", "cYAPKO"))

test_result <- wilcox.test(Expression ~ Condition, data = data_for_plot)
# Extract p-value
p_val <- test_result$p.value

pdf("output/Pathway/violin_embryo_RNAexpression_REACTOME_SIGNALING_BY_WNT_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_REACTOME_SIGNALING_BY_NODAL_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_REACTOME_SIGNALING_BY_HIPPO_Mesenchyme.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_PID_RETINOIC_ACID_PATHWAY_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_REACTOME_SIGNALING_BY_NOTCH_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_WP_BMP_SIGNALING_IN_EYELID_DEVELOPMENT_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_REACTOME_HATS_ACETYLATE_HISTONES_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_WP_MESODERMAL_COMMITMENT_PATHWAY_Mixed_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_WP_VEGFAVEGFR2_SIGNALING_PATHWAY_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/violin_embryo_RNAexpression_REACTOME_SIGNALING_BY_FGFR_Blood_Progenitor_1.pdf", width=3, height=2)

ggplot(data_for_plot, aes(x=Condition, y=Expression, fill=Condition)) +
  geom_violin(scale="width", trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  labs(y="Expression") +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("WT" = "blue", "cYAPKO" = "red"))+
  theme_bw() +
  geom_signif(comparisons = list(c("WT", "cYAPKO")), y_position = max(data_for_plot$Expression) + 0.5, test = "wilcox.test", map_signif_level = TRUE)
dev.off()


# GSEA plot
library("fgsea")
## Re-calculate DEGs keeping ALL genes
embryo.combined.sct$celltype.stim <- paste(embryo.combined.sct$cluster.annot, embryo.combined.sct$condition,
    sep = "-")
Idents(embryo.combined.sct) <- "celltype.stim"

Paraxial_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Paraxial_Mesoderm-cYAPKO", ident.2 = "Paraxial_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Nascent_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Nascent_Mesoderm-cYAPKO", ident.2 = "Nascent_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Pharyngeal_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Pharyngeal_Mesoderm-cYAPKO", ident.2 = "Pharyngeal_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Mesenchyme <- FindMarkers(embryo.combined.sct, ident.1 = "Mesenchyme-cYAPKO", ident.2 = "Mesenchyme-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
Mixed_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Mixed_Mesoderm-cYAPKO", ident.2 = "Mixed_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
Blood_Progenitor_1 <- FindMarkers(embryo.combined.sct, ident.1 = "Blood_Progenitor_1-cYAPKO", ident.2 = "Blood_Progenitor_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
### save output
write.table(Paraxial_Mesoderm, file = "output/seurat/Paraxial_Mesoderm-cYAPKO_response_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
Paraxial_Mesoderm <- read.delim("output/seurat/Paraxial_Mesoderm-cYAPKO_response_allGenes.txt", header = TRUE, row.names = 1)

write.table(Nascent_Mesoderm, file = "output/seurat/Nascent_Mesoderm-cYAPKO_response_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
Nascent_Mesoderm <- read.delim("output/seurat/Nascent_Mesoderm-cYAPKO_response_allGenes.txt", header = TRUE, row.names = 1)

write.table(Pharyngeal_Mesoderm, file = "output/seurat/Pharyngeal_Mesoderm-cYAPKO_response_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
Pharyngeal_Mesoderm <- read.delim("output/seurat/Pharyngeal_Mesoderm-cYAPKO_response_allGenes.txt", header = TRUE, row.names = 1)
write.table(Mesenchyme, file = "output/seurat/Mesenchyme-cYAPKO_response_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
Mesenchyme <- read.delim("output/seurat/Mesenchyme-cYAPKO_response_allGenes.txt", header = TRUE, row.names = 1)

write.table(Mixed_Mesoderm, file = "output/seurat/Mixed_Mesoderm-cYAPKO_response_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Blood_Progenitor_1, file = "output/seurat/Blood_Progenitor_1-cYAPKO_response_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
###

## load list of genes to test
pathways <- msigdbr("Mus musculus", "C2") 
fgsea_sets <- pathways %>% split(x = .$gene_symbol, f = .$gs_name)

## Rank genes based on FC
genes <- Blood_Progenitor_1 %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

## plot GSEA
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_WNT_Paraxial_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS_Nascent_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_NODAL_Pharyngeal_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_HIPPO_Mesenchyme.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_PID_RETINOIC_ACID_PATHWAY_Pharyngeal_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_NOTCH_Paraxial_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_WP_BMP_SIGNALING_IN_EYELID_DEVELOPMENT_Paraxial_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_HATS_ACETYLATE_HISTONES_Pharyngeal_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_WP_MESODERMAL_COMMITMENT_PATHWAY_Mixed_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION_Nascent_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_WP_VEGFAVEGFR2_SIGNALING_PATHWAY_Nascent_Mesoderm.pdf", width=5, height=3)
pdf("output/Pathway/GSEA_REACTOME_SIGNALING_BY_FGFR_Blood_Progenitor_1.pdf", width=5, height=3)

plotEnrichment(fgsea_sets[["REACTOME_SIGNALING_BY_FGFR"]],
               ranks) + labs(title="REACTOME_SIGNALING_BY_FGFR-Blood_Progenitor_1") +
               theme_bw()
dev.off()






# Test different Pathway collections and generate enrichment plot for each cell types (C8, C5, H)
## import Pathways
pathways <- msigdbr("Mus musculus", "H") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
"Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(embryo.combined.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  cYAPKO <- seurat_extract(embryo.combined.sct,
                           meta1 = "condition", value_meta1 = "cYAPKO",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_H_", cluster, ".txt")       # CHANGE HERE PATHWAYS !!!!!!
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

### Plot to count how many pos/neg FC
SCPA_H <- list()
for (cluster in clusters) {
  # Create the filename based on the cluster name
  input_filename <- paste0("output/Pathway/SCPA_H_", cluster, ".txt")
  
  # Import the data
  current_data <- read.delim(input_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Store the data in the list with the cluster name as the list name
  SCPA_H[[cluster]] <- current_data
}
SCPA_C5 <- list()
for (cluster in clusters) {
  # Create the filename based on the cluster name
  input_filename <- paste0("output/Pathway/SCPA_C5_", cluster, ".txt")
  
  # Import the data
  current_data <- read.delim(input_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Store the data in the list with the cluster name as the list name
  SCPA_C5[[cluster]] <- current_data
}
SCPA_C2 <- list()
for (cluster in clusters) {
  # Create the filename based on the cluster name
  input_filename <- paste0("output/Pathway/SCPA_", cluster, ".txt")
  
  # Import the data
  current_data <- read.delim(input_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Store the data in the list with the cluster name as the list name
  SCPA_C2[[cluster]] <- current_data
}


### Tidy the data; convert list into tibble

tidy_C2 <- map2_df(SCPA_C2, names(SCPA_C2), ~ {
    data_frame <- .x
    data_frame$cluster <- .y
    return(data_frame)
  }) %>%
  add_column(db = "C2")
tidy_C5 <- map2_df(SCPA_C5, names(SCPA_C5), ~ {
    data_frame <- .x
    data_frame$cluster <- .y
    return(data_frame)
  }) %>%
  add_column(db = "C5")
tidy_H <- map2_df(SCPA_H, names(SCPA_H), ~ {
    data_frame <- .x
    data_frame$cluster <- .y
    return(data_frame)
  }) %>%
  add_column(db = "H")

tidy_all = tidy_C2 %>%
  bind_rows(tidy_C5) %>%
  bind_rows(tidy_H) %>% 
  as_tibble() %>%
  mutate(FC_direction = ifelse(FC < 0, "positive", "negative"))

# Plot
pdf("output/Pathway/Count_Pathway_FC_direction.pdf", width=10, height=3)
tidy_all %>%
  filter(qval > 1.4) %>%
  group_by(cluster, db, FC_direction) %>%
  summarise(n=n()) %>%
  ggplot(aes(x = FC_direction, y = n, fill = cluster)) +  # Moved fill inside aes
  geom_bar(stat = "identity", position = "dodge") +  # Add dodging
  facet_wrap(~db, scale = "free") +  # Facet by db
  labs(title = "Distribution of FC_direction per cluster",
       y = "Count", x = "FC Direction") +
  theme_bw()
dev.off()



tidy_all = tidy_C2 %>%
  bind_rows(tidy_C5) %>%
  bind_rows(tidy_H) %>% 
  as_tibble() %>%
  mutate(FC_direction = case_when(
    FC < -4 ~ "positive",
    FC > 4 ~ "negative",
    TRUE ~ "neutral"  # for values between -1 and 1
  ))


# Plot
pdf("output/Pathway/Count_Pathway_FC_direction_4treshold.pdf", width=10, height=3)
tidy_all %>%
  filter(qval > 1.4) %>%
  group_by(cluster, db, FC_direction) %>%
  summarise(n=n()) %>%
  ggplot(aes(x = FC_direction, y = n, fill = cluster)) +  # Moved fill inside aes
  geom_bar(stat = "identity", position = "dodge") +  # Add dodging
  facet_wrap(~db, scale = "free") +  # Facet by db
  labs(title = "Distribution of FC_direction per cluster",
       y = "Count", x = "FC Direction") +
  theme_bw()
dev.off()


# separate if more cells in WT or KO
max_cells_df <- data.frame (cluster  = c("Primordial_Germ_Cells","Unknow_2","Unknow_1", "Gut", "Notocord","Surface_Ectoderm", "Blood_Progenitor_2", "Blood_Progenitor_1","Mixed_Mesoderm","Mesenchyme","Haematodenothelial_progenitors","Nascent_Mesoderm","Pharyngeal_Mesoderm","Paraxial_Mesoderm","Caudal_Mesoderm","Somitic_Mesoderm","ExE_Ectoderm","Epiblast_PrimStreak"),
                  max_cells = c("WT", "cYAPKO","WT", "cYAPKO", "WT","WT","cYAPKO", "WT", "WT","WT","WT","WT","WT","WT","WT","WT", "cYAPKO","WT")
                  )

tidy_all = tidy_C2 %>%
  bind_rows(tidy_C5) %>%
  bind_rows(tidy_H) %>% 
  as_tibble() %>%
  mutate(FC_direction = ifelse(FC < 0, "positive", "negative"))
                 
# Plot
pdf("output/Pathway/Count_Pathway_FC_direction_maxCells.pdf", width=10, height=10)
tidy_all %>%
  filter(qval > 1.4) %>%
  group_by(cluster, db, FC_direction) %>%
  summarise(n=n()) %>% 
    left_join(max_cells_df) %>%
  ggplot(aes(x = FC_direction, y = n, fill = cluster)) +  # Moved fill inside aes
  geom_bar(stat = "identity", position = "dodge") +  # Add dodging
  facet_grid(max_cells~db, scale = "free") +  # Facet by db
  labs(title = "Distribution of FC_direction per cluster",
       y = "Count", x = "FC Direction") +
  theme_bw()
dev.off()


# SAME PLOT WITH SCT NORM, lets just do C2
DefaultAssay(embryo.combined.sct) <- "SCT"


# Test different Pathway collections and generate enrichment plot for each cell types (C8, C5, H)
## import Pathways
pathways <- msigdbr("Mus musculus", "C2") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
"Primordial_Germ_Cells",
  "Unknow_2",
  "Unknow_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm",
  "Epiblast_PrimStreak"
)


### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(embryo.combined.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  cYAPKO <- seurat_extract(embryo.combined.sct,
                           meta1 = "condition", value_meta1 = "cYAPKO",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCT_SCPA_C2_", cluster, ".txt")       # CHANGE HERE PATHWAYS !!!!!!
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

### Plot to count how many pos/neg FC
SCT_SCPA_C2 <- list()
for (cluster in clusters) {
  # Create the filename based on the cluster name
  input_filename <- paste0("output/Pathway/SCT_SCPA_C2_", cluster, ".txt")
  
  # Import the data
  current_data <- read.delim(input_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Store the data in the list with the cluster name as the list name
  SCT_SCPA_C2[[cluster]] <- current_data
}



### Tidy the data; convert list into tibble

tidy_SCT_C2 <- map2_df(SCT_SCPA_C2, names(SCT_SCPA_C2), ~ {
    data_frame <- .x
    data_frame$cluster <- .y
    return(data_frame)
  }) %>%
  add_column(db = "C2") %>% 
  as_tibble() %>%
  mutate(FC_direction = case_when(
    FC < -0 ~ "positive",
    FC > 0 ~ "negative",
    TRUE ~ "neutral"  # for values between -1 and 1
  ))

# Plot
pdf("output/Pathway/Count_Pathway_FC_direction_SCT_C2.pdf", width=6, height=3)
tidy_SCT_C2 %>%
  filter(qval > 1.4) %>%
  group_by(cluster, db, FC_direction) %>%
  summarise(n=n()) %>%
  ggplot(aes(x = FC_direction, y = n, fill = cluster)) +  # Moved fill inside aes
  geom_bar(stat = "identity", position = "dodge") +  # Add dodging
  facet_wrap(~db, scale = "free") +  # Facet by db
  labs(title = "Distribution of FC_direction per cluster",
       y = "Count", x = "FC Direction") +
  theme_bw()
dev.off()


# inverse pop1 and pop2:

#### Extract data for WT and cYAPKO based on current value
WT <- seurat_extract(embryo.combined.sct,
                      meta1 = "condition", value_meta1 = "WT",
                      meta2 = "cluster.annot", value_meta2 = "Blood_Progenitor_1")

cYAPKO <- seurat_extract(embryo.combined.sct,
                          meta1 = "condition", value_meta1 = "cYAPKO",
                          meta2 = "cluster.annot", value_meta2 = "Blood_Progenitor_1")

##### Compare pathways
cYAPKO_WT <- compare_pathways(samples = list(cYAPKO, WT),    # initially WT, cYPAKO
                              pathways = pathways,
                              parallel = TRUE, cores = 8)

write.table(cYAPKO_WT, file = "output/Pathway/SCPA_Blood_Progenitor_1_cYAPKO_WT", sep = "\t", quote = FALSE, row.names = FALSE)



cYAPKO_WT_filter <- cYAPKO_WT %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

pdf("output/Pathway/plot_embryo_msigdb_Blood_Progenitor_1_cYAPKO_WT.pdf", width=5, height=5)
ggplot(cYAPKO_WT_filter, aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = cYAPKO_WT_filter$color, stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
dev.off()



# change downsample on Blood_Progenitor_2 to 30 cells

#### Extract data for WT and cYAPKO based on current value
WT <- seurat_extract(embryo.combined.sct,
                      meta1 = "condition", value_meta1 = "WT",
                      meta2 = "cluster.annot", value_meta2 = "Blood_Progenitor_2")

cYAPKO <- seurat_extract(embryo.combined.sct,
                          meta1 = "condition", value_meta1 = "cYAPKO",
                          meta2 = "cluster.annot", value_meta2 = "Blood_Progenitor_2")

##### Compare pathways
WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),    # initially WT, cYPAKO
                              pathways = pathways,
                              downsample = 30,
                              parallel = TRUE, cores = 8)

write.table(WT_cYAPKO, file = "output/Pathway/SCPA_Blood_Progenitor_2_downsample30", sep = "\t", quote = FALSE, row.names = FALSE)



WT_cYAPKO_filter <- WT_cYAPKO %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

pdf("output/Pathway/plot_embryo_msigdb_Blood_Progenitor_2_downsample30.pdf", width=5, height=5)
ggplot(WT_cYAPKO_filter, aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = WT_cYAPKO_filter$color, stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
dev.off()


# Rank and color FC with treshold 5
### import all C2 comparisons

SCPA_C2 <- list()
for (cluster in clusters) {
  # Create the filename based on the cluster name
  input_filename <- paste0("output/Pathway/SCPA_", cluster, ".txt")
  
  # Import the data
  current_data <- read.delim(input_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Store the data in the list with the cluster name as the list name
  SCPA_C2[[cluster]] <- current_data
}

### convert list into tibble
tidy_SCT_C2 <- map2_df(SCPA_C2, names(SCPA_C2), ~ {
    data_frame <- .x
    data_frame$cluster <- .y
    return(data_frame)
  }) %>%
  as_tibble() %>%
  mutate(FC_direction = case_when(
    FC < -2.5 ~ "positive",
    FC > 2.5 ~ "negative",
    TRUE ~ "neutral"  # for values between -1 and 1
  ))

tidy_SCT_C2_Blood_Progenitor_1 = 
  tidy_SCT_C2 %>%
  filter(cluster == "Epiblast_PrimStreak", qval > 1.4) %>%
  arrange(qval) %>%
  mutate(rank = (row_number() - 1) / (n() - 1) * 100)

### plot rank


#### color scale
pdf("output/Pathway/plot_rank-Epiblast_PrimStreak.pdf", width=5, height=5)
tidy_SCT_C2_Blood_Progenitor_1 %>%
  ggplot(., aes(x = qval, y = rank)) +
  geom_jitter(aes(color = -FC, alpha = FC_direction), size=0.1, width = 0.2, height = 2) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, guide = "colourbar") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 1, "neutral" = 0.5)) +
  labs(x = "Qval", y = "Pathway rank", size = "q-value", color = "Fold Change") +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept =95)
dev.off()



#### color
color_palette <- c(
  "positive" = "red",
  "negative" = "blue",
  "neutral"  = "grey50"
)
alpha_values <- c(
  "positive" = 1,
  "negative" = 1,
  "neutral"  = 0.5
)

pdf("output/Pathway/plot_rank_Blood_Progenitor_1.pdf", width=5, height=5)
tidy_SCT_C2_Blood_Progenitor_1 %>%
  ggplot(., aes(qval, rank)) +
  geom_jitter(aes(color = FC_direction, alpha = FC_direction), size = 0.1,width = 0.2, height = 2) +
  scale_color_manual(values = color_palette) + 
  scale_alpha_manual(values = alpha_values) +
  labs(x = "Qval", y = "Pathway rank") +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept =95)
dev.off()


# Comand to ouptut all the genes including in the C2 db

pathways <- msigdbr("Homo sapiens", "C2") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

long_df <- do.call(rbind, pathways) # This will bind all list of tibble into a single tible!
long_df <- long_df %>%
  group_by(Pathway) %>%
  mutate(id = row_number()) %>%
  ungroup()
wide_df <- tidyr::spread(long_df, key = Pathway, value = Genes, fill = NA)

filtered_pathway_1 = c("BIOCARTA_WNT_PATHWAY",
  "KEGG_WNT_SIGNALING_PATHWAY",
  "WP_WNT_SIGNALING",
  "WP_WNT_SIGNALING_PATHWAY",
  "WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY",
  "PID_WNT_CANONICAL_PATHWAY",
  "PID_WNT_NONCANONICAL_PATHWAY",
  "PID_WNT_SIGNALING_PATHWAY",
  "WILLERT_WNT_SIGNALING",
  "WNT_SIGNALING",
  "REACTOME_SIGNALING_BY_WNT",
  "REACTOME_WNT5A_DEPENDENT_INTERNALIZATION_OF_FZD4",
  "REACTOME_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING",
  "REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING",
  "REACTOME_NEGATIVE_REGULATION_OF_TCF_DEPENDENT_SIGNALING_BY_WNT_LIGAND_ANTAGONISTS",
  "REACTOME_DEACTIVATION_OF_THE_BETA_CATENIN_TRANSACTIVATING_COMPLEX",
  "REACTOME_FORMATION_OF_THE_BETA_CATENIN_TCF_TRANSACTIVATING_COMPLEX",
  "REACTOME_DEGRADATION_OF_BETA_CATENIN_BY_THE_DESTRUCTION_COMPLEX",
  "REACTOME_DEGRADATION_OF_AXIN",
  "REACTOME_BETA_CATENIN_PHOSPHORYLATION_CASCADE",
  "REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT",
  "REACTOME_NEGATIVE_REGULATION_OF_TCF_DEPENDENT_SIGNALING_BY_WNT_LIGAND_ANTAGONISTS",
  "REACTOME_WNT5A_DEPENDENT_INTERNALIZATION_OF_FZD4",
  "WP_TGFBETA_RECEPTOR_SIGNALING",
  "WP_TGFBETA_SIGNALING_PATHWAY",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "BIOCARTA_TGFB_PATHWAY",
  "KARAKAS_TGFB1_SIGNALING",
  "PID_TGFBR_PATHWAY",
  "REACTOME_SIGNALING_BY_NODAL",
  "REACTOME_DOWNREGULATION_OF_SMAD2_3_SMAD4_TRANSCRIPTIONAL_ACTIVITY",
  "REACTOME_TRANSCRIPTIONAL_ACTIVITY_OF_SMAD2_SMAD3_SMAD4_HETEROTRIMER",
  "REACTOME_SMAD2_SMAD3_SMAD4_HETEROTRIMER_REGULATES_TRANSCRIPTION",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  "REACTOME_SIGNALING_BY_ACTIVIN")

filtered_pathway_2 = c(  "REACTOME_DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING",
  "REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS",
  "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",
  "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
  "WP_CANONICAL_AND_NONCANONICAL_TGFB_SIGNALING",
  "REACTOME_SIGNALING_BY_HIPPO",
  "REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION",
  "WP_HIPPO_SIGNALING_REGULATION_PATHWAYS",
  "WP_MECHANOREGULATION_AND_PATHOLOGY_OF_YAPTAZ_VIA_HIPPO_AND_NONHIPPO_MECHANISMS",
  "WP_HIPPOYAP_SIGNALING_PATHWAY",
  "REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION",
  "WP_MECHANOREGULATION_AND_PATHOLOGY_OF_YAPTAZ_VIA_HIPPO_AND_NONHIPPO_MECHANISMS",
  "REACTOME_SIGNALING_BY_RETINOIC_ACID",
  "REACTOME_THE_CANONICAL_RETINOID_CYCLE_IN_RODS_TWILIGHT_VISION",
  "KEGG_RETINOL_METABOLISM",
  "PID_RETINOIC_ACID_PATHWAY",
  "WP_4HYDROXYTAMOXIFEN_DEXAMETHASONE_AND_RETINOIC_ACIDS_REGULATION_OF_P27_EXPRESSION",
  "REACTOME_SIGNALING_BY_NOTCH",
  "REACTOME_SIGNALING_BY_NOTCH1",
  "REACTOME_SIGNALING_BY_NOTCH2",
  "REACTOME_SIGNALING_BY_NOTCH3",
  "REACTOME_SIGNALING_BY_NOTCH4",
  "REACTOME_NOTCH_HLH_TRANSCRIPTION_PATHWAY" ,
  "REACTOME_PRE_NOTCH_PROCESSING_IN_GOLGI",
  "REACTOME_PRE_NOTCH_EXPRESSION_AND_PROCESSING",
  "REACTOME_NOTCH2_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",
  "REACTOME_NOTCH3_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",
  "REACTOME_ACTIVATED_NOTCH1_TRANSMITS_SIGNAL_TO_THE_NUCLEUS",
  "REACTOME_NOTCH3_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION",
  "REACTOME_NOTCH4_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION",
  "KEGG_NOTCH_SIGNALING_PATHWAY" ,
  "PID_NOTCH_PATHWAY" ,
  "WP_NOTCH_SIGNALING" ,
  "WP_NOTCH_SIGNALING_PATHWAY" ,
  "WP_CANONICAL_AND_NONCANONICAL_NOTCH_SIGNALING",
  "REACTOME_SIGNALING_BY_BMP",
  "WP_BMP_SIGNALING_IN_EYELID_DEVELOPMENT",
  "PID_BMP_PATHWAY")

filtered_pathway_3 = c( "REACTOME_HDMS_DEMETHYLATE_HISTONES",
  "REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA",
  "REACTOME_HDACS_DEACETYLATE_HISTONES",
  "REACTOME_RMTS_METHYLATE_HISTONE_ARGININES",
  "REACTOME_HATS_ACETYLATE_HISTONES",
  "REACTOME_PKMTS_METHYLATE_HISTONE_LYSINES",
  "WP_INTERACTOME_OF_POLYCOMB_REPRESSIVE_COMPLEX_2_PRC2",
  "BIOCARTA_HDAC_PATHWAY",
  "WP_ENDODERM_DIFFERENTIATION",
  "WP_MESODERMAL_COMMITMENT_PATHWAY",
  "WP_ECTODERM_DIFFERENTIATION",
  "REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_DOPAMINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_SEROTONIN_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
  "WP_COMPLEMENT_SYSTEM_IN_NEURONAL_DEVELOPMENT_AND_PLASTICITY",
  "REACTOME_MECP2_REGULATES_NEURONAL_RECEPTORS_AND_CHANNELS",
  "REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS",
  "WP_MAJOR_RECEPTORS_TARGETED_BY_EPINEPHRINE_AND_NOREPINEPHRINE",
  "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
  "REACTOME_NEUREXINS_AND_NEUROLIGINS",
  "REACTOME_ASSEMBLY_AND_CELL_SURFACE_PRESENTATION_OF_NMDA_RECEPTORS",
  "REACTOME_NEGATIVE_REGULATION_OF_NMDA_RECEPTOR_MEDIATED_NEURONAL_TRANSMISSION",
  "REACTOME_ACTIVATION_OF_NMDA_RECEPTORS_AND_POSTSYNAPTIC_EVENTS",
  "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY",
  "REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES",
  "REACTOME_BLOOD_GROUP_SYSTEMS_BIOSYNTHESIS",
  "REACTOME_HEME_SIGNALING",
  "REACTOME_SCAVENGING_OF_HEME_FROM_PLASMA",
  "REACTOME_RESPONSE_OF_EIF2AK1_HRI_TO_HEME_DEFICIENCY",
  "KEGG_VEGF_SIGNALING_PATHWAY")

filtered_pathway_4 = c( 
  "REACTOME_SIGNALING_BY_VEGF",
  "REACTOME_VEGFR2_MEDIATED_VASCULAR_PERMEABILITY",
  "REACTOME_VEGFR2_MEDIATED_CELL_PROLIFERATION",
  "BIOCARTA_VEGF_PATHWAY",
  "WP_VEGFAVEGFR2_SIGNALING_PATHWAY",
  "PID_VEGFR1_PATHWAY",
  "PID_VEGFR1_2_PATHWAY",
  "REACTOME_SIGNALING_BY_FGFR",
  "REACTOME_SIGNALING_BY_FGFR1",
  "REACTOME_FGFR1_LIGAND_BINDING_AND_ACTIVATION",
  "REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR1",
  "REACTOME_SHC_MEDIATED_CASCADE_FGFR1",
  "REACTOME_FRS_MEDIATED_FGFR1_SIGNALING",
  "REACTOME_PI_3K_CASCADE_FGFR1",
  "REACTOME_NEGATIVE_REGULATION_OF_FGFR1_SIGNALING",
  "REACTOME_SIGNALING_BY_FGFR2",
  "REACTOME_FGFR2_LIGAND_BINDING_AND_ACTIVATION",
  "REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR2",
  "REACTOME_FRS_MEDIATED_FGFR2_SIGNALING",
  "REACTOME_PI_3K_CASCADE_FGFR2",
  "REACTOME_PHOSPHOLIPASE_C_MEDIATED_CASCADE_FGFR2",
  "REACTOME_NEGATIVE_REGULATION_OF_FGFR2_SIGNALING",
  "REACTOME_FGFR2_ALTERNATIVE_SPLICING",
  "REACTOME_SIGNALING_BY_FGFR3",
  "REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR3",
  "REACTOME_SHC_MEDIATED_CASCADE_FGFR3",
  "REACTOME_FRS_MEDIATED_FGFR3_SIGNALING",
  "REACTOME_PI_3K_CASCADE_FGFR3",
  "REACTOME_NEGATIVE_REGULATION_OF_FGFR3_SIGNALING",
  "REACTOME_SIGNALING_BY_FGFR4",
  "REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR4",
  "REACTOME_SHC_MEDIATED_CASCADE_FGFR4",
  "REACTOME_FRS_MEDIATED_FGFR4_SIGNALING",
  "REACTOME_PI_3K_CASCADE_FGFR4",
  "REACTOME_PHOSPHOLIPASE_C_MEDIATED_CASCADE_FGFR4",
  "REACTOME_NEGATIVE_REGULATION_OF_FGFR4_SIGNALING")

filtered_pathway_all = c(filtered_pathway_1, filtered_pathway_2, filtered_pathway_3, filtered_pathway_4)


wide_df_filter = wide_df %>%
  dplyr::select( filtered_pathway_all )
write.table(wide_df_filter, file = "output/Pathway/msigdbr_Mus_Musculus_C2_allCategories.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(wide_df_filter, file = "output/Pathway/msigdbr_Homo_Sapiens_C2_allCategories.txt", sep = "\t", quote = FALSE, row.names = TRUE)


```

--> **To reproduce code; follow "PRO WINNER" or "WINNER PRO"...**

--> For automatic cell type annotation; the [EasyCellType] [shiny app](https://biostatistics.mdanderson.org/shinyapps/EasyCellType/) has been tested. 
----> Works great! Can play with the pval cutoff; between 0.3-0.5 is ok

--> For **GSEA**; seems only use log2FC ranking is ok (some people recomend "p-value multiplied with the sign of log-fold change for the ranking" but obverall majority recomend only log2FC see [here](https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_05_dge.html#Gene_Set_Enrichment_Analysis_(GSEA))).
----> I followed [this](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html) and [this](https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_05_dge.html#Gene_Set_Enrichment_Analysis_(GSEA)), tutorials for data vizualization and prep, respectively
----> We need to re-calculate DEGs FC because previous method do NOT give information for all the genes (it is to avoid multi-testing correction); recommended parameter found [here](https://github.com/satijalab/seurat/issues/532)


--> Share to Conchi the Conserved Marker list (`srat_all_conserved_markers_embryo.txt`). To avoid confusion, I did some filtering: For each cell type; I only keep log2FC positive (= correspond to gene more highly express in this cell types) and I told her to filter per pvalue which is the max_pvalue. Like this, she will only see the highly express genes in each cluster
----> Conchi provided an updated list (check Dropbox)

--> Fucntional analysis
----> ReactomeGSEA to help cell type annotation (not for WT cYAPKO)
----> [SCPA](https://jackbibby1.github.io/SCPA/) for comparison
------> Gene list from Retinoic Acid downloaded from [this](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_SIGNALING_BY_RETINOIC_ACID.html); and then converted into mice ortholog [online](https://biit.cs.ut.ee/gprofiler/orth); file `Pathway/gProfiler_hsapiens_mmusculus_9-6-2023_3-52-54.csv`
------> Alternative method: [this]https://jackbibby1.github.io/SCPA/articles/seurat_comparison.html and [this](https://jackbibby1.github.io/SCPA/articles/comparing_two_populations.html). Here I compare all the cell population together but compare my condition
------> [This](https://jackbibby1.github.io/SCPA/articles/disease_comparison.html) is the tutorial to follow; disease comparison
------> following this discussion for mice:
- download gmt file for [Retinoic Acid](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/REACTOME_SIGNALING_BY_RETINOIC_ACID.html)

--> Seems that most pathways are ACTIVATED! Weird, might be a biase? That is occuring in each cell type. To verify:
(from [this](https://jackbibby1.github.io/SCPA/reference/compare_seurat.html):*The FC output is generated from a running sum of mean changes in gene expression from all genes of the pathway. It's calculated from average pathway expression in population1 - population2, so a negative FC means the pathway is higher in population2.*)
----> Test other pathway db
----> Test using SCT normalization

- *NOTE: it last forever, so I used parralell core processing `srun --mem=500g --cpus-per-task=10 --pty bash -l`*
- *NOTE: I work on `conda activate scRNAseqV1`*
- *NOTE: Pretty cool paper that uses pathway enrichment to guide clustering! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9985338/*


Comand to print how many DEGs in each cell type:
```bash
awk -F'\t' '$6 < 0.05' output/seurat/Primordial_Germ_Cells-cYAPKO_response.txt | wc -l # 11
awk -F'\t' '$6 < 0.05' output/seurat/Unknow_2-cYAPKO_response.txt | wc -l # 1
awk -F'\t' '$6 < 0.05' output/seurat/Unknow_1-cYAPKO_response.txt | wc -l # 49
awk -F'\t' '$6 < 0.05' output/seurat/Gut-cYAPKO_response.txt | wc -l # 5
awk -F'\t' '$6 < 0.05' output/seurat/Notocord-cYAPKO_response.txt | wc -l # 116
awk -F'\t' '$6 < 0.05' output/seurat/Surface_Ectoderm-cYAPKO_response.txt | wc -l # 293
awk -F'\t' '$6 < 0.05' output/seurat/Blood_Progenitor_2-cYAPKO_response.txt | wc -l # 129
awk -F'\t' '$6 < 0.05' output/seurat/Blood_Progenitor_1-cYAPKO_response.txt | wc -l # 516
awk -F'\t' '$6 < 0.05' output/seurat/Mixed_Mesoderm-cYAPKO_response.txt | wc -l # 410
awk -F'\t' '$6 < 0.05' output/seurat/Mesenchyme-cYAPKO_response.txt | wc -l # 895
awk -F'\t' '$6 < 0.05' output/seurat/Haematodenothelial_progenitors-cYAPKO_response.txt | wc -l # 547
awk -F'\t' '$6 < 0.05' output/seurat/Nascent_Mesoderm-cYAPKO_response.txt | wc -l # 527
awk -F'\t' '$6 < 0.05' output/seurat/Pharyngeal_Mesoderm-cYAPKO_response.txt | wc -l # 621
awk -F'\t' '$6 < 0.05' output/seurat/Paraxial_Mesoderm-cYAPKO_response.txt | wc -l # 585
awk -F'\t' '$6 < 0.05' output/seurat/Caudal_Mesoderm-cYAPKO_response.txt | wc -l # 434
awk -F'\t' '$6 < 0.05' output/seurat/Somitic_Mesoderm-cYAPKO_response.txt | wc -l # 515
awk -F'\t' '$6 < 0.05' output/seurat/ExE_Ectoderm-cYAPKO_response.txt | wc -l # 429
awk -F'\t' '$6 < 0.05' output/seurat/Epiblast_PrimStreak-cYAPKO_response.txt | wc -l # 465
```


--> Let's now start all over again with the new clustering version V2 (post Conchi meeting 20231005)

# embryo; start over again with new clustering V2 (post Conchi meeting 20231005)

Here is below the complete code (copy paste and modified from V1):

```bash
conda activate scRNAseqV2
```

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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

## Load the matrix and Create SEURAT object
load("output/soupX/out_embryo_control_e775_mice.RData")
srat_WT <- CreateSeuratObject(counts = out, project = "WT") # 32,285 features across 5,568 samples

load("output/soupX/out_embryo_cYAPKO_e775_mice.RData")
srat_cYAPKO <- CreateSeuratObject(counts = out, project = "cYAPKO") # 32,285 features across 4,504 samples


# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
srat_WT[["percent.mt"]] <- PercentageFeatureSet(srat_WT, pattern = "^mt-")
srat_WT[["percent.rb"]] <- PercentageFeatureSet(srat_WT, pattern = "^Rp[sl]")

srat_cYAPKO[["percent.mt"]] <- PercentageFeatureSet(srat_cYAPKO, pattern = "^mt-")
srat_cYAPKO[["percent.rb"]] <- PercentageFeatureSet(srat_cYAPKO, pattern = "^Rp[sl]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/embryo_control.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_WT <- AddMetaData(srat_WT,doublets)
srat_WT$Doublet_score <- as.numeric(srat_WT$Doublet_score) # make score as numeric
head(srat_WT[[]])

doublets <- read.table("output/doublets/embryo_cYAPKO.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_cYAPKO <- AddMetaData(srat_cYAPKO,doublets)
srat_cYAPKO$Doublet_score <- as.numeric(srat_cYAPKO$Doublet_score) # make score as numeric
head(srat_cYAPKO[[]])


## After seeing the plot; add QC information in our seurat object
## V1 QC; optimal with vst V1; dim 19 k param 70 es 0.9 --> THE WINNER pro winner
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC == 'Pass','Low_nFeature',srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_WT@meta.data$QC,sep = ','),srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$percent.mt > 25 & srat_WT@meta.data$QC == 'Pass','High_MT',srat_WT@meta.data$QC)
srat_WT[['QC']] <- ifelse(srat_WT@meta.data$nFeature_RNA < 2000 & srat_WT@meta.data$QC != 'Pass' & srat_WT@meta.data$QC != 'High_MT',paste('High_MT',srat_WT@meta.data$QC,sep = ','),srat_WT@meta.data$QC)
table(srat_WT[['QC']])
## 
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC == 'Pass','Low_nFeature',srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_cYAPKO@meta.data$QC,sep = ','),srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$percent.mt > 25 & srat_cYAPKO@meta.data$QC == 'Pass','High_MT',srat_cYAPKO@meta.data$QC)
srat_cYAPKO[['QC']] <- ifelse(srat_cYAPKO@meta.data$nFeature_RNA < 2000 & srat_cYAPKO@meta.data$QC != 'Pass' & srat_cYAPKO@meta.data$QC != 'High_MT',paste('High_MT',srat_cYAPKO@meta.data$QC,sep = ','),srat_cYAPKO@meta.data$QC)
table(srat_cYAPKO[['QC']])




## subset my seurat object to only analyze the cells that pass the QC
srat_WT <- subset(srat_WT, subset = QC == 'Pass')
srat_cYAPKO <- subset(srat_cYAPKO, subset = QC == 'Pass')
srat_WT$condition <- "WT"
srat_cYAPKO$condition <- "cYAPKO"



mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
srat_WT <- NormalizeData(srat_WT, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_WT)
srat_WT <- ScaleData(srat_WT, features = all.genes) # zero-centres and scales it

srat_cYAPKO <- NormalizeData(srat_cYAPKO, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_cYAPKO)
srat_cYAPKO <- ScaleData(srat_cYAPKO, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
srat_WT <- CellCycleScoring(srat_WT, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_WT[[]]$Phase)
srat_cYAPKO <- CellCycleScoring(srat_cYAPKO, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_cYAPKO[[]]$Phase)

set.seed(42)



set.seed(42)

## TEST RNA regression  --> THE WINNER PRO WINNER !!! RNA regression 
srat_WT <- SCTransform(srat_WT, method = "glmGamPoi", ncells = 4812, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO <- SCTransform(srat_cYAPKO, method = "glmGamPoi", ncells = 3621, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(srat_WT = srat_WT, srat_cYAPKO = srat_cYAPKO)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

embryo.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
embryo.combined.sct <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT")

set.seed(42)

## PLAY WITH RESOLUTION AFTER CONCHI MEETING 20231005
## individualize germ cells; but fail at separating  
DefaultAssay(embryo.combined.sct) <- "integrated"

embryo.combined.sct <- RunPCA(embryo.combined.sct, verbose = FALSE, npcs = 19)
embryo.combined.sct <- RunUMAP(embryo.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryo.combined.sct <- FindNeighbors(embryo.combined.sct, reduction = "pca", k.param = 5, dims = 1:19)
embryo.combined.sct <- FindClusters(embryo.combined.sct, resolution = 0.3, verbose = FALSE, algorithm = 4)

embryo.combined.sct$condition <- factor(embryo.combined.sct$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st

pdf("output/seurat/UMAP_control_cYAPKO_V2clust__.pdf", width=10, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()





### Updated marker gene list from Alex + WINNER pro-winner
Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # 1 \  1 $ 1 + 1 %  1
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # 3 \ 4 $ 3  + 4 %  2
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # 2 \ 3 $ 8  + 3 %  7
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # 4 \ 2 $ 2 + 2 %  3
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # 5 \ 5 $ 5 +  5 %  4
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # 2 \ 3 $ 4  +  3 %  8
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # 7 \ 7 $ 6  +    %  5
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # 8 \ 8 $ 9  +    %  9
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # 11  \  $ 11 +   %  11
Caudal_Neurectoderm = c("Olig3","Hes3","Lrrc7","Cntfr") # 6  \   $ 5  +  %  4
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # 5  \   $ 7 + %  6
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # 9  \  $  10 + %  12
Intermediate_Mesoderm = c("Bik","Arg1","Meox1","Ism1") # 4 \  2  $  2 + %  mix of 3 5 6 (let's not name it as more general than 3 5 6)
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # 12  \  $  12 + %  13
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # 4 and 10 \  2  $  2 + %  10 (and 3)
Notocord = c("Shh","Noto","Vtn","Fgg") #  13 \   $  13 ( wierd separated) + %  14
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # 14  \   $  14 + %  15
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3") # between 12 and 13    + %  18



DefaultAssay(embryo.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Epiblast_PrimStreak_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Epiblast_PrimStreak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_ExE_Ectoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Nascent_Mesoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Somitic_Mesoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Somitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Caudal_Mesoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Caudal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Haematodenothelial_progenitors_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Haematodenothelial_progenitors, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Paraxial_Mesoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Paraxial_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Mesenchyme_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Mesenchyme, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Blood_Progenitor_1_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Blood_Progenitor_1, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Caudal_Neurectoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Caudal_Neurectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Pharyngeal_Mesoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Pharyngeal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Blood_Progenitor_2_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Blood_Progenitor_2, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Intermediate_Mesoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Intermediate_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Surface_Ectoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Surface_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Mixed_Mesoderm_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Mixed_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Notocord_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Notocord, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Gut_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Gut, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Primordial_Germ_Cells_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Primordial_Germ_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Noto_V2clust.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = "Noto", max.cutoff = 3, cols = c("grey", "red"))
dev.off()






# Rename cluster
Epiblast_PrimStreak  # cluster_1
ExE_Ectoderm_1 # cluster_2
Nascent_Mesoderm  # cluster_3
Somitic_Mesoderm  # cluster_4
Caudal_Mesoderm # cluster_5
Paraxial_Mesoderm  # cluster_6
Pharyngeal_Mesoderm # cluster_7
Haematodenothelial_progenitors  # cluster_8
Mesenchyme # cluster_9
Blood_Progenitor_1 # cluster_10 (and 3)
Mixed_Mesoderm # cluster_11
Blood_Progenitor_2 # cluster_12
Surface_Ectoderm  # cluster_13
Gut #  cluster_14
Unknown_1 # cluster_15
Notocord # cluster_16
Unknown_2 # cluster_17
Primordial_Germ_Cells # cluster_18
ExE_Ectoderm_2 # cluster_19



new.cluster.ids <- c(
  "Epiblast_PrimStreak", 
  "ExE_Ectoderm_1", 
  "Nascent_Mesoderm", 
  "Somitic_Mesoderm", 
  "Caudal_Mesoderm", 
  "Paraxial_Mesoderm", 
  "Pharyngeal_Mesoderm", 
  "Haematodenothelial_progenitors", 
  "Mesenchyme", 
  "Blood_Progenitor_1", 
  "Mixed_Mesoderm", 
  "Blood_Progenitor_2", 
  "Surface_Ectoderm", 
  "Gut", 
  "Unknown_1", 
  "Notocord", 
  "Unknown_2", 
  "Primordial_Germ_Cells", 
  "ExE_Ectoderm_2"
)

names(new.cluster.ids) <- levels(embryo.combined.sct)
embryo.combined.sct <- RenameIdents(embryo.combined.sct, new.cluster.ids)

embryo.combined.sct$cluster.annot <- Idents(embryo.combined.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_control_cYAPKO_label_V2clust.pdf", width=12, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

#overlapping condition
pdf("output/seurat/UMAP_control_cYAPKO_label_overlap_V2clust.pdf", width=6, height=5)
DimPlot(embryo.combined.sct, reduction = "umap", group.by = "condition", pt.size = 0.000001, cols = c("blue","red"))
dev.off()


# All in dotplot
DefaultAssay(embryo.combined.sct) <- "SCT"

Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # cluster_1
ExE_Ectoderm_1 = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # cluster_2
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # cluster_3
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # cluster_4
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # cluster_5
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # cluster_6
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # cluster_7
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # cluster_8
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # cluster_9
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # cluster_10 (and 3)
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # cluster_11
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # cluster_12
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # cluster_13
Notocord = c("Shh","Noto","Vtn","Fgg") #  cluster_14
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # cluster_15
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3")  # cluster_18


all_markers <- c(
  "Hesx1", "Otx2", "Pax2", "Otx2os1",
  "Tex19.1", "Elf5", "Wnt7b", "Dnmt3l",
  "Meox1", "Aldh1a2", "Synm", "Pcdh8",
  "Wnt3a", "Nkx1-2", "Apln", "Rxrg",
  "Ebf2", "Ptgfr", "Ebf3", "Col23a1",
  "Nkx2-5", "Tbx5", "Mef2c", "Myocd",
  "Col9a1", "Mesp1", "Osr1",
  "Hoxa9", "Hoxa10", "Tbx4", "Pitx1",
  "Tdo2", "Adamts2", "Colec11", "Snai2",
  "Tbx6", "Dll1", "Hes7", "Fgf17",
  "Sox18", "Esam", "Rhoj", "Flt4",
  "Gypa", "Gata1", "Cited4", "Gfi1b",
  "Tfap2a", "Npnt", "Wnt4", "Id4",
  "Shh", "Noto", "Vtn", "Fgg",
  "Pga5", "Hs3cst1", "Wfdc1", "Islr2",
  "Dppa3", "Sprr2a3", "Irf1", "Ifitm3"
)


levels(embryo.combined.sct) <- c(
  "ExE_Ectoderm_2",
  "Primordial_Germ_Cells",
  "Unknown_2",
  "Notocord",
  "Unknown_1",
  "Gut",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Mixed_Mesoderm",
  "Blood_Progenitor_1",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "Nascent_Mesoderm",
  "ExE_Ectoderm_1",
  "Epiblast_PrimStreak"
)
pdf("output/seurat/DotPlot_SCT_control_cYAPKO_ident_V2clust.pdf", width=16.5, height=4.5)
DotPlot(embryo.combined.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()





## Downsampling with bootstrap to compare the nb of cell per cell types

library("tidyverse")

### Identify the unique clusters
unique_clusters <- unique(Idents(embryo.combined.sct))

### Create empty matrices to store cell counts
control_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
cYAPKO_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(control_clusters_counts) <- unique_clusters
colnames(cYAPKO_clusters_counts) <- unique_clusters

### Loop through 100 iterations
embryo.combined.sct_WT <- which(embryo.combined.sct$orig.ident == 'WT')
embryo.combined.sct_cYAPKO <- which(embryo.combined.sct$orig.ident == 'cYAPKO')

for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  embryo.combined.sct_WT_downsample <- sample(embryo.combined.sct_WT, 3621)
  embryo.combined.sct_integrated_downsample <- embryo.combined.sct[,c(embryo.combined.sct_cYAPKO, embryo.combined.sct_WT_downsample)]

  # Count nb of cells in each cluster
  control_clusters <- table(Idents(embryo.combined.sct_integrated_downsample)[embryo.combined.sct_integrated_downsample$condition == "WT"])
  cYAPKO_clusters <- table(Idents(embryo.combined.sct_integrated_downsample)[embryo.combined.sct_integrated_downsample$condition == "cYAPKO"])

  # Align the counts with the unique clusters
  control_clusters_counts[i, names(control_clusters)] <- as.numeric(control_clusters)
  cYAPKO_clusters_counts[i, names(cYAPKO_clusters)] <- as.numeric(cYAPKO_clusters)
}



### Calculate mean and standard error
mean_control_clusters <- colMeans(control_clusters_counts)
mean_cYAPKO_clusters <- colMeans(cYAPKO_clusters_counts)
std_error_WT_clusters <- apply(control_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("WT", "cYAPKO")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- control_clusters_counts[j,i]
    contingency_table[1,2] <- cYAPKO_clusters_counts[j,i]
    contingency_table[2,1] <- sum(control_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(cYAPKO_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_control_clusters),
  untreated = mean_control_clusters,
  dasatinib = mean_cYAPKO_clusters,
  std_error_WT = std_error_WT_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_WT, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "WT", "cYAPKO"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.05, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st
plot_data$cluster <- factor(plot_data$cluster, levels = c("Epiblast_PrimStreak",  # Early developmental stage
  "ExE_Ectoderm_1",         # Extraembryonic ectoderm; related to the epiblast
  "ExE_Ectoderm_2",         # Extraembryonic ectoderm; related to the epiblast
  "Surface_Ectoderm",     # Forms from ectoderm
  "Notocord",             # Important mesodermal structure
  "Nascent_Mesoderm",     # Early mesoderm
  "Paraxial_Mesoderm",    # Forms alongside notochord
  "Somitic_Mesoderm",     # Gives rise to somites
  "Pharyngeal_Mesoderm",  # Part of head mesoderm
  "Caudal_Mesoderm",      # Posterior mesoderm
  "Mixed_Mesoderm",       # Intermediate mesodermal stage
  "Mesenchyme",           # Loose, undifferentiated cells often derived from mesoderm
  "Haematodenothelial_progenitors", # Blood and endothelial progenitors
  "Blood_Progenitor_1",   # Blood precursors
  "Blood_Progenitor_2",   # Blood precursors
  "Gut",                  # Derived from endoderm
  "Primordial_Germ_Cells",# Early germ cell precursors
  "Unknown_1",             # Unknown categories placed at the end
  "Unknown_2")) # Reorder untreated 1st


# Plotting using ggplot2
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling10_clean_embryo_V2clust.pdf", width=9, height=4)
ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = filter(plot_data, condition == "cYAPKO"),
    aes(label = significance, y = value + std_error_WT_clusters),
    vjust = -0.8,
    position = position_dodge(0.9), size = 5
  ) +
  scale_fill_manual(values = c("WT" = "blue", "cYAPKO" = "red")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,900)
dev.off()




# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

DefaultAssay(embryo.combined.sct) <- "RNA"

embryo.combined.sct <- NormalizeData(embryo.combined.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(embryo.combined.sct)
embryo.combined.sct <- ScaleData(embryo.combined.sct, features = all.genes) # zero-centres and scales it


## what genes change in different conditions for cells of the same type

embryo.combined.sct$celltype.stim <- paste(embryo.combined.sct$cluster.annot, embryo.combined.sct$condition,
    sep = "-")
Idents(embryo.combined.sct) <- "celltype.stim"

# use RNA corrected count for DEGs
embryo.combined.sct <- PrepSCTFindMarkers(embryo.combined.sct)


## Automation::
cell_types <- c(
  "Primordial_Germ_Cells",
  "Unknown_2",
  "Notocord",
  "Unknown_1",
  "Gut",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Mixed_Mesoderm",
  "Blood_Progenitor_1",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "Nascent_Mesoderm",
  "ExE_Ectoderm_1",
  "Epiblast_PrimStreak")

for (cell_type in cell_types) {
  response_name <- paste(cell_type, "cYAPKO.response", sep = ".")
  ident_1 <- paste(cell_type, "-cYAPKO", sep = "")
  ident_2 <- paste(cell_type, "-WT", sep = "")

  response <- FindMarkers(embryo.combined.sct, assay = "RNA", ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
  
  print(head(response, n = 15))
  
  file_name <- paste("output/seurat/", cell_type, "-cYAPKO_response_V2clust.txt", sep = "")
  write.table(response, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE)
}



### Find all markers 
all_markers <- FindAllMarkers(embryo.combined.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_WT_cYAPKO_all_markers_V2clust.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# Display the top 10 CONSERVED marker genes of each cluster
Idents(embryo.combined.sct) <- "cluster.annot"

## DEGs cluster versus all other
Primordial_Germ_Cells.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Primordial_Germ_Cells", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Primordial_Germ_Cells")
Unknow_2.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Unknown_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknown_2")
Unknow_1.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Unknown_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknown_1")
Gut.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Gut", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Gut")
Notocord.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Notocord", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Notocord")
Surface_Ectoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Surface_Ectoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Surface_Ectoderm")
Blood_Progenitor_2.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Blood_Progenitor_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_2")
Blood_Progenitor_1.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Blood_Progenitor_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_1")
Mixed_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Mixed_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mixed_Mesoderm")
Mesenchyme.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Mesenchyme", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mesenchyme")
Haematodenothelial_progenitors.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Haematodenothelial_progenitors", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Haematodenothelial_progenitors")
Nascent_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Nascent_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Nascent_Mesoderm")
Pharyngeal_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Pharyngeal_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Pharyngeal_Mesoderm")
Paraxial_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Paraxial_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Paraxial_Mesoderm")
Caudal_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Caudal_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Caudal_Mesoderm")
Somitic_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Somitic_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Somitic_Mesoderm")
ExE_Ectoderm_1.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "ExE_Ectoderm_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "ExE_Ectoderm_1")
ExE_Ectoderm_2.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "ExE_Ectoderm_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "ExE_Ectoderm_2")
Epiblast_PrimStreak.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Epiblast_PrimStreak", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Epiblast_PrimStreak")


Unknown_2.conserved = Unknow_2.conserved
Unknown_1.conserved = Unknow_1.conserved


## Combine all conserved markers into one data frame
all_conserved <- bind_rows(Primordial_Germ_Cells.conserved, Unknown_2.conserved, Unknown_1.conserved, Gut.conserved, Notocord.conserved, Surface_Ectoderm.conserved, Blood_Progenitor_2.conserved, Blood_Progenitor_1.conserved, Mixed_Mesoderm.conserved, Mesenchyme.conserved, Haematodenothelial_progenitors.conserved, Nascent_Mesoderm.conserved, Pharyngeal_Mesoderm.conserved, Paraxial_Mesoderm.conserved, Caudal_Mesoderm.conserved, Somitic_Mesoderm.conserved, ExE_Ectoderm_1.conserved, ExE_Ectoderm_2.conserved, Epiblast_PrimStreak.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers_embryo_V2clust.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak")))

## Find the top 3 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 3) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak")))



## Write the top 10 conserved markers for each cluster to a file
write.table(top10_conserved, file = "output/seurat/srat_top10_conserved_markers_embryo_V2clust.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Visualize the top 10/3 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(embryo.combined.sct) <- c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak")

pdf("output/seurat/DotPlot_SCT_top3_conserved_embryo_V2clust.pdf", width=18, height=5)
DotPlot(embryo.combined.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# save
## saveRDS(embryo.combined.sct, file = "output/seurat/embryo.combined.sct_V2clust.rds")
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct_V2clust.rds")




# Check some genes
DefaultAssay(embryo.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

## post 20231005 Conchi meeting
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Aldh1a2_Cyp26a1_V2clust.pdf", width=10, height=13)
FeaturePlot(embryo.combined.sct, features = c("Aldh1a2", "Cyp26a1"), max.cutoff = 10, cols = c("grey", "red"), split.by = "condition")
dev.off()

## post 20240110 Conchi email
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_20240110geneList_V2clust.pdf", width=10, height=75)
FeaturePlot(embryo.combined.sct, features = c("Crabp2", "Hand1", "Rbp1", "Aldh1a2", "Cyp26a1", "Hoxa1", "Fgf8", "Nr2f2", "Tnnt2", "Actc1", "Myl6", "Myl7", "Myl3", "Tbx1", "Tbx5", "Mab21l2", "Lbh"), max.cutoff = 10, cols = c("grey", "red"), split.by = "condition")
dev.off()



# Compare WT and cYAPKO using SCPA ##########################################

### installation (in scRNAseqV1)
# devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
# devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
# devtools::install_github("jackbibby1/SCPA")
# --> Last fail so installed dependencies
# BiocManager::install(c("clustermole", "ComplexHeatmap", "singscore"))
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
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct_V2clust.rds")
DefaultAssay(embryo.combined.sct) <- "RNA" # Recommended 



# Import manually curated gene list as a list
#### Need to folllow this format
pathways <- msigdbr("Mus musculus", "C2") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise
pathways$REACTOME_SIGNALING_BY_RETINOIC_ACID$Genes
### Convert table to list of genes
Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

pathways <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt"),
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt"),
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt"),
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")
)

pathways$Hippo$Genes


# Compare pathway activity between condition within  
cell_types <- unique(embryo.combined.sct$cluster.annot)
embryo.combined.sct_split <- SplitObject(embryo.combined.sct, split.by = "condition")

# SCPA comparison
## keep all columns
scpa_out_all <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryo.combined.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryo.combined.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_all[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}
## keep only PAthway and qval column (Best for heatmap qval representation)
scpa_out <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryo.combined.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryo.combined.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}


saveRDS(scpa_out_all, file = "output/Pathway/scpa_out_eachCellTypes_ManualGeneLists.rds")

console_output <- capture.output(print(scpa_out_all))
writeLines(console_output, "output/Pathway/scpa_out_console_ManualGeneLists.txt")

## Save output as R object
# saveRDS(scpa_out, file = "output/Pathway/scpa_out_eachCellTypes.rds")
# scpa_out_all <- readRDS("output/Pathway/scpa_out_all_eachCellTypes.rds")


# plot output

# Code to save output for each cell type comparison
clusters = c(
  "Epiblast_PrimStreak", 
  "ExE_Ectoderm_1", 
  "Nascent_Mesoderm", 
  "Somitic_Mesoderm", 
  "Caudal_Mesoderm", 
  "Paraxial_Mesoderm", 
  "Pharyngeal_Mesoderm", 
  "Haematodenothelial_progenitors", 
  "Mesenchyme", 
  "Blood_Progenitor_1", 
  "Mixed_Mesoderm", 
  "Blood_Progenitor_2", 
  "Surface_Ectoderm", 
  "Gut", 
  "Unknown_1", 
  "Notocord", 
  "Unknown_2", 
  "Primordial_Germ_Cells", 
  "ExE_Ectoderm_2"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(embryo.combined.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  cYAPKO <- seurat_extract(embryo.combined.sct,
                           meta1 = "condition", value_meta1 = "cYAPKO",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_V2clust_ManualGeneLists_", cluster, ".txt")
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}



# BALOW CLEAN CODE; dotplot>enrichplot>heatmap>violinplot
# DOT PLOT pepresetnation
## load all the comparison for each cell type (FC qval information)
clusters = c(
  "Epiblast_PrimStreak", 
  "ExE_Ectoderm_1", 
  "Nascent_Mesoderm", 
  "Somitic_Mesoderm", 
  "Caudal_Mesoderm", 
  "Paraxial_Mesoderm", 
  "Pharyngeal_Mesoderm", 
  "Haematodenothelial_progenitors", 
  "Mesenchyme", 
  "Blood_Progenitor_1", 
  "Mixed_Mesoderm", 
  "Blood_Progenitor_2", 
  "Surface_Ectoderm", 
  "Gut", 
  "Unknown_1", 
  "Notocord", 
  "Unknown_2", 
  "Primordial_Germ_Cells"
)
## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_V2clust_ManualGeneLists_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster))
# --> import the gene pathways (used table from Conchi `Pathwyas of interest*.xlsx`)
pathways = c("Hippo", "Nodal_TGF", "WNT", "Vitamin_A")
all_data_pathways = as_tibble(all_data) %>% 
  filter(Pathway %in% pathways)

# tidy the df
all_data_pathways_tidy <- all_data_pathways %>%
  filter(qval >= 1.4) 



# REFINE COLOR
custom_color <- function(fc_value){
  ifelse(fc_value >= 5, "strong_blue", 
         ifelse(fc_value > 2, "light_blue",
                ifelse(fc_value <= -5, "strong_red", 
                       ifelse(fc_value < -2, "light_red", "grey"))))
}

# Add a column for this custom color
all_data_pathways_tidy <- all_data_pathways_tidy %>%
  mutate(custom_col = sapply(FC, custom_color))

pdf("output/Pathway/dotplot_V2clust_ManualGeneLists_Hippo_FCtresh.pdf", width=10, height=4)
ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, color = custom_col), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = c("strong_blue" = "dodgerblue4",
                                "light_blue" = "lightblue2",
                                "grey" = "grey",
                                "light_red" = "indianred1",
                                "strong_red" = "red3")) +
  theme_bw() +
  labs(size = "q-value", color = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



pdf("output/Pathway/dotplot_V2clust_ManualGeneLists_grey.pdf", width=10, height=4)
ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, fill = "black"), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  labs(size = "q-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# heatmap

Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

pathways <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt"),
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt"),
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt"),
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")
)

pathways$Hippo$Genes



genes_of_interest <- pathways$Hippo$Genes
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(embryo.combined.sct@assays$RNA@data)]
### extract WT and cYAPKO gene expression values from RNA assay
WT_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "WT" & embryo.combined.sct$cluster.annot == "Unknown_1"]]
cYAPKO_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO" & embryo.combined.sct$cluster.annot == "Unknown_1"]]
### mean expression values for each gene
WT_mean <- rowMeans(WT_expression)
cYAPKO_mean <- rowMeans(cYAPKO_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  WT = WT_mean,
  cYAPKO = cYAPKO_mean
) %>% 
pivot_longer(cols = c(WT, cYAPKO), names_to = "Condition", values_to = "Expression")
## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(WT_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)




# if few genes:
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Blood_Progenitor_1.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Caudal_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Epiblast_PrimStreak.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Haematodenothelial_progenitors.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Mesenchyme.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Mixed_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Nascent_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Paraxial_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Pharyngeal_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Somitic_Mesoderm.pdf", width=6, height=2)


library("viridis")
ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()

# if many genesL
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Blood_Progenitor_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Blood_Progenitor_2.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Caudal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Epiblast_PrimStreak.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_ExE_Ectoderm_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Haematodenothelial_progenitors.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Mesenchyme.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Mixed_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Somitic_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Surface_Ectoderm.pdf", width=3, height=2)

pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Blood_Progenitor_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Caudal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Epiblast_PrimStreak.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_ExE_Ectoderm_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Haematodenothelial_progenitors.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Mesenchyme.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Mixed_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Somitic_Mesoderm.pdf", width=3, height=2)

pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Blood_Progenitor_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Blood_Progenitor_2.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Caudal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Epiblast_PrimStreak.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_ExE_Ectoderm_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Haematodenothelial_progenitors.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Mesenchyme.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Mixed_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Somitic_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Surface_Ectoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Unknown_1.pdf", width=3, height=2)


ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") 
dev.off()








# GSEA plot
library("fgsea")
## Re-calculate DEGs keeping ALL genes
embryo.combined.sct$celltype.stim <- paste(embryo.combined.sct$cluster.annot, embryo.combined.sct$condition,
    sep = "-")
Idents(embryo.combined.sct) <- "celltype.stim"

Paraxial_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Paraxial_Mesoderm-cYAPKO", ident.2 = "Paraxial_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Nascent_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Nascent_Mesoderm-cYAPKO", ident.2 = "Nascent_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Pharyngeal_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Pharyngeal_Mesoderm-cYAPKO", ident.2 = "Pharyngeal_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Mesenchyme <- FindMarkers(embryo.combined.sct, ident.1 = "Mesenchyme-cYAPKO", ident.2 = "Mesenchyme-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
Mixed_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Mixed_Mesoderm-cYAPKO", ident.2 = "Mixed_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
Blood_Progenitor_1 <- FindMarkers(embryo.combined.sct, ident.1 = "Blood_Progenitor_1-cYAPKO", ident.2 = "Blood_Progenitor_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Epiblast_PrimStreak <- FindMarkers(embryo.combined.sct, ident.1 = "Epiblast_PrimStreak-cYAPKO", ident.2 = "Epiblast_PrimStreak-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
ExE_Ectoderm_1 <- FindMarkers(embryo.combined.sct, ident.1 = "ExE_Ectoderm_1-cYAPKO", ident.2 = "ExE_Ectoderm_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Somitic_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Somitic_Mesoderm-cYAPKO", ident.2 = "Somitic_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Caudal_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Caudal_Mesoderm-cYAPKO", ident.2 = "Caudal_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Haematodenothelial_progenitors <- FindMarkers(embryo.combined.sct, ident.1 = "Haematodenothelial_progenitors-cYAPKO", ident.2 = "Haematodenothelial_progenitors-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Blood_Progenitor_2 <- FindMarkers(embryo.combined.sct, ident.1 = "Blood_Progenitor_2-cYAPKO", ident.2 = "Blood_Progenitor_2-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Surface_Ectoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Surface_Ectoderm-cYAPKO", ident.2 = "Surface_Ectoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Gut <- FindMarkers(embryo.combined.sct, ident.1 = "Gut-cYAPKO", ident.2 = "Gut-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Unknown_1 <- FindMarkers(embryo.combined.sct, ident.1 = "Unknown_1-cYAPKO", ident.2 = "Unknown_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Unknown_2 <- FindMarkers(embryo.combined.sct, ident.1 = "Unknown_2-cYAPKO", ident.2 = "Unknown_2-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Notocord <- FindMarkers(embryo.combined.sct, ident.1 = "Notocord-cYAPKO", ident.2 = "Notocord-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Primordial_Germ_Cells <- FindMarkers(embryo.combined.sct, ident.1 = "Primordial_Germ_Cells-cYAPKO", ident.2 = "Primordial_Germ_Cells-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  





### save output
#### write.table(Primordial_Germ_Cells, file = "output/seurat/Primordial_Germ_Cells-cYAPKO_response_V2clust_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#### import all clsuter DEGs output :
cluster_types <- c("Epiblast_PrimStreak", "ExE_Ectoderm_1", "Nascent_Mesoderm", 
                   "Somitic_Mesoderm", "Caudal_Mesoderm", "Paraxial_Mesoderm", 
                   "Pharyngeal_Mesoderm", "Haematodenothelial_progenitors", 
                   "Mesenchyme", "Blood_Progenitor_1", "Mixed_Mesoderm", 
                   "Blood_Progenitor_2", "Surface_Ectoderm", "Gut", "Unknown_1", 
                   "Notocord", "Unknown_2", "Primordial_Germ_Cells")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-cYAPKO_response_V2clust_allGenes.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}


## load list of genes to test
Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

fgsea_sets <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt")$Genes,
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt")$Genes,
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt")$Genes,
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")$Genes
)





## Rank genes based on FC
genes <- Blood_Progenitor_1 %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(ES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -NES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()






## plot GSEA
pdf("output/Pathway/GSEA_V2clust_ManualGeneLists_Hippo-Primordial_Germ_Cells.pdf", width=5, height=3)

plotEnrichment(fgsea_sets[["Hippo"]],
               ranks) + labs(title="Hippo-Primordial_Germ_Cells") +
               theme_bw()
dev.off()





# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c("Epiblast_PrimStreak", "ExE_Ectoderm_1", "Nascent_Mesoderm", 
                   "Somitic_Mesoderm", "Caudal_Mesoderm", "Paraxial_Mesoderm", 
                   "Pharyngeal_Mesoderm", "Haematodenothelial_progenitors", 
                   "Mesenchyme", "Blood_Progenitor_1", "Mixed_Mesoderm", 
                   "Blood_Progenitor_2", "Surface_Ectoderm", "Gut", "Unknown_1", 
                   "Notocord", "Unknown_2", "Primordial_Germ_Cells")

## Initialize an empty list to store the results for each cluster type
all_results <- list()
## Loop over each cluster type
for (cluster in cluster_types) {
  
  # Extract genes for the current cluster
  genes <- get(cluster) %>% 
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC)
  
  ranks <- deframe(genes)
  
  # Run GSEA for the current cluster
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(ES))
  
  # Extract summary table and add cluster column
  fgseaResTidy_summary = fgseaResTidy %>% 
    dplyr::select(pathway, pval, padj, ES, size, NES) %>%
    mutate(cluster = cluster) %>%
    arrange(padj) %>% 
    head()
  
  # Store results in the list
  all_results[[cluster]] <- fgseaResTidy_summary
}
## Combine results from all cluster types into one table
final_results <- bind_rows(all_results, .id = "cluster")


write.table(final_results, file = c("output/Pathway/gsea_output_V2clust_ManualGeneLists.txt"), sep = "\t", quote = FALSE, row.names = FALSE)



# Heatmap all GSEA


pdf("output/Pathway/heatmap_gsea_padj.pdf", width=5, height=5)
ggplot(final_results, aes(x=cluster, y=pathway, fill=ES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Enrichment\nScore") +
  geom_text(aes(label=sprintf("%.2f", ES)), 
            color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()



pdf("output/Pathway/heatmap_gsea_pval_greyTile.pdf", width=5, height=5)
ggplot(final_results, aes(x=cluster, y=pathway)) + 
  geom_tile(aes(fill = ifelse(pval <= 0.05, ES, NA)), color = "black") +  # Conditional fill based on significance
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Enrichment\nScore", na.value="grey") +
  # geom_text(aes(label=sprintf("%.2f", ES)), 
  #           color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change between pvalue, qvalue,p.adjust
  #           size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1

dev.off()




# SCPA and GSEA plot
gsea_output = final_results %>% 
  rename("pathway" = "Pathway") %>%
  dplyr::select(Pathway, pval, padj, ES, cluster)
scpa_output = all_data_pathways_tidy  %>%
  dplyr::select(Pathway, qval, cluster)

scpa_gsea = gsea_output %>%
  left_join(scpa_output)

## Filter data for GSEA padj < 0.05
filtered_data <- scpa_gsea
## Assign color based on your criteria
filtered_data$is_significant <- filtered_data$padj < 0.05

pdf("output/Pathway/dotplot_scpa_gsea_V2clust.pdf", width=10, height=4)

ggplot(filtered_data, aes(x=cluster, y=Pathway, size=qval)) +
  geom_point(aes(color=ifelse(is_significant, ES, NA), 
                 shape=is_significant), alpha=0.7) +
  scale_color_gradient2(low="blue", high="red", midpoint=0, na.value="black", name = "GSEA Enrichment score") +
  scale_size_continuous(name="SCPA qval") +  # Title for size scale
  scale_shape_manual(values=c(`TRUE`=16, `FALSE`=16)) + # Using circle shape for both
  theme_light() +
  labs(title="SCPA and GSEA", x="Cluster", y="Pathway") +
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

```








# embryo; start over again with final cluster name (pre- publication of cardiac paper) _V3

--> I copy below the entire code from `# embryo; start over again with new clustering V2 (post Conchi meeting 20231005)`; need to re-run everything as cluster name changed...


Here is below the complete code (copy paste and modified from V1):

```bash
conda activate scRNAseqV2
```

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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

# LOAD previous backup to keep same UMAP shape
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct_V2clust.rds")

pdf("output/seurat/UMAP_control_cYAPKO_V3.pdf", width=10, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()





### Updated marker gene list FINAL (in pink in file docs/annotated clusters and markers.pptx)
### original
Neural_Ectoderm = Pax2, Six3, Hesx1, Otx2, Otx1, Otx2os1
ExE_Ectoderm = Elf5, Wnt7b, Dnmt3l, Tex19.1
Presomitic_Mesoderm = TBX6, Meox1, Dll1, aldh1a2
Neuromesodermal_Progenitors = Nkx1-2, T, Apln, Wnt3a, Rxrg # NMP
Second_Heart_Field  = Tbx1, Nr2f2, cyp26a1, Ebf2, Ptgfr, Ebf3, Col23a1
First_Heart_Field = Tbx5, Nkx2-5, Myocd
Cardiac_Mesoderm = Col9a1, Mesp1, Osr1, Cited1
Posterior_Lateral_Plate_Mesoderm = Hoxa9, Hoxa10, Tbx4, Pitx1
Juxta_Cardiac_Field = Tagln2, Postn, Hand1, Acta2, Snai2, Mab21l2
Somitic_Mesoderm = Hes7, Tbx6, aldh1a2, Dll1, Fgf17
Hemangioblasts = Tal1, sox18, etv2, Esam, Rhoj, Flt4, Sox18
Hematopoietic_Stem_Cells = Runx1, Gypa, Gata1, Cited4, Gfi1b # HSC
Surface_Ectoderm = Tfap2a, Npnt, Wnt4, Id4
Notochord = Shh, Noto, Vtn, Fgg
Endoderm = Foxa2, Sox17, Grem2, sox7, Pga5, Hs3cst1, Wfdc1, Islr2
Visceral_Endoderm = Slc39a8, Amot, Apoa4, Sox17, Sox7, Foxa2, Gata4
Ectoderm_undefined = Six3, Otx2, tfap2a
Primordial_Germ_Cells = Dppa3

### removing duplicate genes
Neural_Ectoderm = Pax2, Six3, Hesx1, Otx2, Otx1, Otx2os1
ExE_Ectoderm = Elf5, Wnt7b, Dnmt3l, Tex19.1

Presomitic_Mesoderm = Meox1, Dll1, aldh1a2 # Removed Tbx6
Somitic_Mesoderm = Hes7, Tbx6, Fgf17 # Removed aldh1a2 and Dll1

Neuromesodermal_Progenitors = Nkx1-2, T, Apln, Wnt3a, Rxrg # NMP
Second_Heart_Field  = Tbx1, Nr2f2, cyp26a1, Ebf2, Ptgfr, Ebf3, Col23a1
First_Heart_Field = Tbx5, Nkx2-5, Myocd
Cardiac_Mesoderm = Col9a1, Mesp1, Osr1, Cited1
Posterior_Lateral_Plate_Mesoderm = Hoxa9, Hoxa10, Tbx4, Pitx1
Juxta_Cardiac_Field = Tagln2, Postn, Hand1, Acta2, Snai2, Mab21l2
Hemangioblasts = Tal1, sox18, etv2, Esam, Rhoj, Flt4
Hematopoietic_Stem_Cells = Runx1, Gypa, Gata1, Cited4, Gfi1b # HSC
Notochord = Shh, Noto, Vtn, Fgg

Endoderm = Sox17, Grem2, Pga5, Hs3cst1, Wfdc1, Islr2 # Removed Sox7 and Foxa2
Visceral_Endoderm = Slc39a8, Amot, Apoa4, Sox7, Foxa2, Gata4 # Removed Sox17

Surface_Ectoderm = Tfap2a, Npnt, Wnt4, Id4
Ectoderm_undefined = Six3 # Otx2 and tfap2a removed

Primordial_Germ_Cells = Dppa3


#--> For theg enes in double, do plot and decide where to put them; that's the one in 'removed'



Neural_Ectoderm = c("Pax2", "Six3", "Hesx1", "Otx2", "Otx1", "Otx2os1")
ExE_Ectoderm = c("Elf5", "Wnt7b", "Dnmt3l", "Tex19.1")
Presomitic_Mesoderm = c("Meox1", "Dll1", "Aldh1a2", "Tbx6")  
Somitic_Mesoderm = c("Hes7", "Fgf17") # Removed aldh1a2, Tbx6 and Dll1
Neuromesodermal_Progenitors = c("Nkx1-2", "T", "Apln", "Wnt3a", "Rxrg") # NMP
Second_Heart_Field  = c("Tbx1", "Nr2f2", "Cyp26a1", "Ebf2", "Ptgfr", "Ebf3", "Col23a1")
First_Heart_Field = c("Tbx5", "Nkx2-5", "Myocd")
Cardiac_Mesoderm = c("Col9a1", "Mesp1", "Osr1", "Cited1")
Posterior_Lateral_Plate_Mesoderm = c("Hoxa9", "Hoxa10", "Tbx4", "Pitx1")
Juxta_Cardiac_Field = c("Tagln2", "Postn", "Hand1", "Acta2", "Snai2", "Mab21l2")
Hemangioblasts = c("Tal1", "Sox18", "Etv2", "Esam", "Rhoj", "Flt4")
Hematopoietic_Stem_Cells = c("Runx1", "Gypa", "Gata1", "Cited4", "Gfi1b") # HSC
Notochord = c("Shh", "Noto", "Vtn", "Fgg")
Endoderm = c("Sox17", "Grem2", "Pga5", "Hs3st1", "Wfdc1", "Islr2", "Sox7", "Foxa2", "Sox17") 
Visceral_Endoderm = c("Slc39a8", "Amot", "Apoa4", "Gata4") # Removed , Sox7, Foxa2 and Sox17
Surface_Ectoderm = c("Tfap2a", "Npnt", "Wnt4", "Id4")
Ectoderm_undefined = c("Six3") # Otx2 and tfap2a removed
Primordial_Germ_Cells = c("Dppa3")







DefaultAssay(embryo.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Neural_Ectoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Neural_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_ExE_Ectoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Presomitic_Mesoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Presomitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Somitic_Mesoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Somitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Neuromesodermal_Progenitors_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Neuromesodermal_Progenitors, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Second_Heart_Field_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Second_Heart_Field, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_First_Heart_Field_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = First_Heart_Field, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Cardiac_Mesoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Cardiac_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Posterior_Lateral_Plate_Mesoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Posterior_Lateral_Plate_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Juxta_Cardiac_Field_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Juxta_Cardiac_Field, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Hemangioblasts_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Hemangioblasts, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Hematopoietic_Stem_Cells_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Hematopoietic_Stem_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Notochord_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Notochord, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Endoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Visceral_Endoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Visceral_Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Surface_Ectoderm_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Surface_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Ectoderm_undefined_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Ectoderm_undefined, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Primordial_Germ_Cells_V3.pdf", width=10, height=10)
FeaturePlot(embryo.combined.sct, features = Primordial_Germ_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()









# Rename cluster
Epiblast_PrimStreak  # cluster_1
ExE_Ectoderm_1 # cluster_2
Nascent_Mesoderm  # cluster_3
Somitic_Mesoderm  # cluster_4
Caudal_Mesoderm # cluster_5
Paraxial_Mesoderm  # cluster_6
Pharyngeal_Mesoderm # cluster_7
Haematodenothelial_progenitors  # cluster_8
Mesenchyme # cluster_9
Blood_Progenitor_1 # cluster_10 (and 3)
Mixed_Mesoderm # cluster_11
Blood_Progenitor_2 # cluster_12
Surface_Ectoderm  # cluster_13
Gut #  cluster_14
Unknown_1 # cluster_15
Notocord # cluster_16
Unknown_2 # cluster_17
Primordial_Germ_Cells # cluster_18
ExE_Ectoderm_2 # cluster_19














embryo.combined.sct <- SetIdent(embryo.combined.sct, value = "seurat_clusters") # reverse back identiy to cluster number


pdf("output/seurat/UMAP_control_cYAPKO_noLabel_V3noSplit.pdf", width=6, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 4)
dev.off()



new.cluster.ids <- c(
  "Neural_Ectoderm", 
  "ExE_Ectoderm_1", 
  "Cardiac_Mesoderm", 
  "Presomitic_Mesoderm", 
  "Neuromesodermal_Progenitors", 
  "Second_Heart_Field", 
  "First_Heart_Field", 
  "Posterior_Lateral_Plate_Mesoderm", 
  "Juxta_Cardiac_Field", 
  "Hemangioblasts", 
  "Somitic_Mesoderm", 
  "Hematopoietic_Stem_Cells", 
  "Surface_Ectoderm", 
  "Endoderm", 
  "Ectoderm_undefined", 
  "Notochord", 
  "Visceral_Endoderm", 
  "Primordial_Germ_Cells", 
  "ExE_Ectoderm_2"
)

names(new.cluster.ids) <- levels(embryo.combined.sct)
embryo.combined.sct <- RenameIdents(embryo.combined.sct, new.cluster.ids)

embryo.combined.sct$cluster.annot <- Idents(embryo.combined.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_control_cYAPKO_label_V3.pdf", width=12, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_control_cYAPKO_label_V3noSplit.pdf", width=8, height=6)
DimPlot(embryo.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3.5)
dev.off()

#overlapping condition
pdf("output/seurat/UMAP_control_cYAPKO_label_overlap_V3.pdf", width=6, height=5)
DimPlot(embryo.combined.sct, reduction = "umap", group.by = "condition", pt.size = 0.000001, cols = c("blue","red"))
dev.off()

embryo.combined.sct <- SetIdent(embryo.combined.sct, value = "seurat_clusters") # reverse back identiy to cluster number


# Heatmap with marker genes
DefaultAssay(embryo.combined.sct) <- "SCT"


as.character(1) = c("Pax2", "Six3", "Hesx1", "Otx2", "Otx1", "Otx2os1")
as.character(2) = c("Elf5", "Wnt7b", "Dnmt3l", "Tex19.1")
as.character(3) = c("Col9a1", "Mesp1", "Osr1", "Cited1")
as.character(4) = c("Meox1", "Dll1", "Aldh1a2", "Tbx6")  
as.character(5) = c("Nkx1-2", "T", "Apln", "Wnt3a", "Rxrg") # NMP
as.character(6)  = c("Tbx1", "Nr2f2", "Cyp26a1", "Ebf2", "Ptgfr", "Ebf3", "Col23a1")
as.character(7) = c("Tbx5", "Nkx2-5", "Myocd")
as.character(8) = c("Hoxa9", "Hoxa10", "Tbx4", "Pitx1")
as.character(9) = c("Tagln2", "Postn", "Hand1", "Acta2", "Snai2", "Mab21l2")
as.character(10) = c("Tal1", "Sox18", "Etv2", "Esam", "Rhoj", "Flt4")
as.character(11) = c("Hes7", "Fgf17") # Removed aldh1a2, Tbx6 and Dll1
as.character(12) = c("Runx1", "Gypa", "Gata1", "Cited4", "Gfi1b") # HSC
as.character(13) = c("Tfap2a", "Npnt", "Wnt4", "Id4")
as.character(14) = c("Sox17", "Grem2", "Pga5", "Hs3st1", "Wfdc1", "Islr2", "Sox7", "Foxa2", "Sox17") 
as.character(15) = c("Six3") # Otx2 and tfap2a removed
as.character(16) = c("Shh", "Noto", "Vtn", "Fgg")
as.character(17) = c("Slc39a8", "Amot", "Apoa4", "Gata4") # Removed , Sox7, Foxa2 and Sox17
as.character(18) = c("Dppa3")


all_markers <- c(
  "Pax2", "Six3", "Hesx1", "Otx2", "Otx1", "Otx2os1",
  "Elf5", "Wnt7b", "Dnmt3l", "Tex19.1",
  "Col9a1", "Mesp1", "Osr1", "Cited1",
  "Meox1", "Dll1", "Aldh1a2", "Tbx6",
  "Nkx1-2", "T", "Apln", "Wnt3a", "Rxrg",
  "Tbx1", "Nr2f2", "Cyp26a1", "Ebf2", "Ptgfr", "Ebf3", "Col23a1",
  "Tbx5", "Nkx2-5", "Myocd",
  "Hoxa9", "Hoxa10", "Tbx4", "Pitx1",
  "Tagln2", "Postn", "Hand1", "Acta2", "Snai2", "Mab21l2",
  "Tal1", "Sox18", "Etv2", "Esam", "Rhoj", "Flt4",
  "Hes7", "Fgf17",
  "Runx1", "Gypa", "Gata1", "Cited4", "Gfi1b",
  "Tfap2a", "Npnt", "Wnt4", "Id4",
  "Sox17", "Grem2", "Pga5", "Hs3st1", "Wfdc1", "Islr2", "Sox7", "Foxa2", "Sox17", 
  "Six3",
  "Shh", "Noto", "Vtn", "Fgg",
  "Slc39a8", "Amot", "Apoa4", "Gata4",
  "Dppa3"
)

pdf("output/seurat/DoHeatmap_SCT_control_cYAPKO_ident_V3sct.pdf", width=8, height=9)
DoHeatmap(embryo.combined.sct, assay = "SCT", features = all_markers,angle = 0, hjust = 0.5)
dev.off()

pdf("output/seurat/DoHeatmap_SCT_control_cYAPKO_ident_V3.pdf", width=8, height=9)
DoHeatmap(embryo.combined.sct, assay = "RNA", features = all_markers,angle = 0, hjust = 0.5)
dev.off()
pdf("output/seurat/DoHeatmap_SCT_control_cYAPKO_ident_V3large.pdf", width=12, height=9)
DoHeatmap(embryo.combined.sct, assay = "RNA", features = all_markers,angle = 0, hjust = 0.5)
dev.off()

# All in dotplot



all_markers <- c(
  "Pax2", "Six3", "Hesx1", "Otx2", "Otx1", "Otx2os1",
  "Elf5", "Wnt7b", "Dnmt3l", "Tex19.1",
  "Col9a1", "Mesp1", "Osr1", "Cited1",
  "Meox1", "Dll1", "Aldh1a2", "Tbx6",
  "Nkx1-2", "T", "Apln", "Wnt3a", "Rxrg",
  "Tbx1", "Nr2f2", "Cyp26a1", "Ebf2", "Ptgfr", "Ebf3", "Col23a1",
  "Tbx5", "Nkx2-5", "Myocd",
  "Hoxa9", "Hoxa10", "Tbx4", "Pitx1",
  "Tagln2", "Postn", "Hand1", "Acta2", "Snai2", "Mab21l2",
  "Tal1", "Sox18", "Etv2", "Esam", "Rhoj", "Flt4",
  "Hes7", "Fgf17",
  "Runx1", "Gypa", "Gata1", "Cited4", "Gfi1b",
  "Tfap2a", "Npnt", "Wnt4", "Id4",
  "Sox17", "Grem2", "Pga5", "Hs3st1", "Wfdc1", "Islr2", "Sox7", "Foxa2", 
  "Shh", "Noto", "Vtn", "Fgg",
  "Slc39a8", "Amot", "Apoa4", "Gata4",
  "Dppa3"
)


levels(embryo.combined.sct) <- c("19",
  "18",
  "17",
  "16",
  "15",
  "14",
  "13",
  "12",
  "11",
  "10",
  "9",
  "8",
  "7",
  "6",
  "5",
  "4",
  "3",
  "2",
  "1")

pdf("output/seurat/DotPlot_SCT_control_cYAPKO_ident_V3.pdf", width=16.5, height=4.5)
DotPlot(embryo.combined.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_control_cYAPKO_ident_V3green.pdf", width=16.5, height=4.5)
DotPlot(embryo.combined.sct, assay = "SCT", features = all_markers, cols = c("grey95", "green")) + RotatedAxis()
dev.off()

## dotplot vetical orientation
levels(embryo.combined.sct) <- c("1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7",
  "8",
  "9",
  "10",
  "11",
  "12",
  "13",
  "14",
  "15",
  "16",
  "17",
  "18",
  "19")

pdf("output/seurat/DotPlot_SCT_control_cYAPKO_ident_V3vertical.pdf", width=16.5, height=5)
DotPlot(embryo.combined.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()
pdf("output/seurat/DotPlot_SCT_control_cYAPKO_ident_V3verticalgreen.pdf", width=16.5, height=5)
DotPlot(embryo.combined.sct, assay = "SCT", features = all_markers, cols = c("grey95", "green")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()

## Downsampling with bootstrap to compare the nb of cell per cell types

library("tidyverse")

### Identify the unique clusters
unique_clusters <- unique(Idents(embryo.combined.sct))

### Create empty matrices to store cell counts
control_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
cYAPKO_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(control_clusters_counts) <- unique_clusters
colnames(cYAPKO_clusters_counts) <- unique_clusters

### Loop through 100 iterations
embryo.combined.sct_WT <- which(embryo.combined.sct$orig.ident == 'WT')
embryo.combined.sct_cYAPKO <- which(embryo.combined.sct$orig.ident == 'cYAPKO')

for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  embryo.combined.sct_WT_downsample <- sample(embryo.combined.sct_WT, 3621)
  embryo.combined.sct_integrated_downsample <- embryo.combined.sct[,c(embryo.combined.sct_cYAPKO, embryo.combined.sct_WT_downsample)]

  # Count nb of cells in each cluster
  control_clusters <- table(Idents(embryo.combined.sct_integrated_downsample)[embryo.combined.sct_integrated_downsample$condition == "WT"])
  cYAPKO_clusters <- table(Idents(embryo.combined.sct_integrated_downsample)[embryo.combined.sct_integrated_downsample$condition == "cYAPKO"])

  # Align the counts with the unique clusters
  control_clusters_counts[i, names(control_clusters)] <- as.numeric(control_clusters)
  cYAPKO_clusters_counts[i, names(cYAPKO_clusters)] <- as.numeric(cYAPKO_clusters)
}

XXXXX here 

### Calculate mean and standard error
mean_control_clusters <- colMeans(control_clusters_counts)
mean_cYAPKO_clusters <- colMeans(cYAPKO_clusters_counts)
std_error_WT_clusters <- apply(control_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("WT", "cYAPKO")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- control_clusters_counts[j,i]
    contingency_table[1,2] <- cYAPKO_clusters_counts[j,i]
    contingency_table[2,1] <- sum(control_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(cYAPKO_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_control_clusters),
  untreated = mean_control_clusters,
  dasatinib = mean_cYAPKO_clusters,
  std_error_WT = std_error_WT_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_WT, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "WT", "cYAPKO"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.05, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st
plot_data$cluster <- factor(plot_data$cluster, levels = c("Epiblast_PrimStreak",  # Early developmental stage
  "ExE_Ectoderm_1",         # Extraembryonic ectoderm; related to the epiblast
  "ExE_Ectoderm_2",         # Extraembryonic ectoderm; related to the epiblast
  "Surface_Ectoderm",     # Forms from ectoderm
  "Notocord",             # Important mesodermal structure
  "Nascent_Mesoderm",     # Early mesoderm
  "Paraxial_Mesoderm",    # Forms alongside notochord
  "Somitic_Mesoderm",     # Gives rise to somites
  "Pharyngeal_Mesoderm",  # Part of head mesoderm
  "Caudal_Mesoderm",      # Posterior mesoderm
  "Mixed_Mesoderm",       # Intermediate mesodermal stage
  "Mesenchyme",           # Loose, undifferentiated cells often derived from mesoderm
  "Haematodenothelial_progenitors", # Blood and endothelial progenitors
  "Blood_Progenitor_1",   # Blood precursors
  "Blood_Progenitor_2",   # Blood precursors
  "Gut",                  # Derived from endoderm
  "Primordial_Germ_Cells",# Early germ cell precursors
  "Unknown_1",             # Unknown categories placed at the end
  "Unknown_2")) # Reorder untreated 1st


# Plotting using ggplot2
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling10_clean_embryo_V2clust.pdf", width=9, height=4)
ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = filter(plot_data, condition == "cYAPKO"),
    aes(label = significance, y = value + std_error_WT_clusters),
    vjust = -0.8,
    position = position_dodge(0.9), size = 5
  ) +
  scale_fill_manual(values = c("WT" = "blue", "cYAPKO" = "red")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,900)
dev.off()




# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

DefaultAssay(embryo.combined.sct) <- "RNA"

embryo.combined.sct <- NormalizeData(embryo.combined.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(embryo.combined.sct)
embryo.combined.sct <- ScaleData(embryo.combined.sct, features = all.genes) # zero-centres and scales it


## what genes change in different conditions for cells of the same type

embryo.combined.sct$celltype.stim <- paste(embryo.combined.sct$cluster.annot, embryo.combined.sct$condition,
    sep = "-")
Idents(embryo.combined.sct) <- "celltype.stim"

# use RNA corrected count for DEGs
embryo.combined.sct <- PrepSCTFindMarkers(embryo.combined.sct)


## Automation::
cell_types <- c(
  "Primordial_Germ_Cells",
  "Unknown_2",
  "Notocord",
  "Unknown_1",
  "Gut",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Mixed_Mesoderm",
  "Blood_Progenitor_1",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "Nascent_Mesoderm",
  "ExE_Ectoderm_1",
  "Epiblast_PrimStreak")

for (cell_type in cell_types) {
  response_name <- paste(cell_type, "cYAPKO.response", sep = ".")
  ident_1 <- paste(cell_type, "-cYAPKO", sep = "")
  ident_2 <- paste(cell_type, "-WT", sep = "")

  response <- FindMarkers(embryo.combined.sct, assay = "RNA", ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
  
  print(head(response, n = 15))
  
  file_name <- paste("output/seurat/", cell_type, "-cYAPKO_response_V2clust.txt", sep = "")
  write.table(response, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE)
}



### Find all markers 
all_markers <- FindAllMarkers(embryo.combined.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_WT_cYAPKO_all_markers_V2clust.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# Display the top 10 CONSERVED marker genes of each cluster
Idents(embryo.combined.sct) <- "cluster.annot"

## DEGs cluster versus all other
Primordial_Germ_Cells.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Primordial_Germ_Cells", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Primordial_Germ_Cells")
Unknow_2.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Unknown_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknown_2")
Unknow_1.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Unknown_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknown_1")
Gut.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Gut", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Gut")
Notocord.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Notocord", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Notocord")
Surface_Ectoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Surface_Ectoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Surface_Ectoderm")
Blood_Progenitor_2.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Blood_Progenitor_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_2")
Blood_Progenitor_1.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Blood_Progenitor_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_1")
Mixed_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Mixed_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mixed_Mesoderm")
Mesenchyme.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Mesenchyme", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mesenchyme")
Haematodenothelial_progenitors.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Haematodenothelial_progenitors", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Haematodenothelial_progenitors")
Nascent_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Nascent_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Nascent_Mesoderm")
Pharyngeal_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Pharyngeal_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Pharyngeal_Mesoderm")
Paraxial_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Paraxial_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Paraxial_Mesoderm")
Caudal_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Caudal_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Caudal_Mesoderm")
Somitic_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Somitic_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Somitic_Mesoderm")
ExE_Ectoderm_1.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "ExE_Ectoderm_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "ExE_Ectoderm_1")
ExE_Ectoderm_2.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "ExE_Ectoderm_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "ExE_Ectoderm_2")
Epiblast_PrimStreak.conserved <- FindConservedMarkers(embryo.combined.sct, assay = "RNA", ident.1 = "Epiblast_PrimStreak", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Epiblast_PrimStreak")


Unknown_2.conserved = Unknow_2.conserved
Unknown_1.conserved = Unknow_1.conserved


## Combine all conserved markers into one data frame
all_conserved <- bind_rows(Primordial_Germ_Cells.conserved, Unknown_2.conserved, Unknown_1.conserved, Gut.conserved, Notocord.conserved, Surface_Ectoderm.conserved, Blood_Progenitor_2.conserved, Blood_Progenitor_1.conserved, Mixed_Mesoderm.conserved, Mesenchyme.conserved, Haematodenothelial_progenitors.conserved, Nascent_Mesoderm.conserved, Pharyngeal_Mesoderm.conserved, Paraxial_Mesoderm.conserved, Caudal_Mesoderm.conserved, Somitic_Mesoderm.conserved, ExE_Ectoderm_1.conserved, ExE_Ectoderm_2.conserved, Epiblast_PrimStreak.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers_embryo_V2clust.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak")))

## Find the top 3 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 3) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak")))



## Write the top 10 conserved markers for each cluster to a file
write.table(top10_conserved, file = "output/seurat/srat_top10_conserved_markers_embryo_V2clust.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Visualize the top 10/3 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(embryo.combined.sct) <- c("Primordial_Germ_Cells",
  "Unknown_2",
  "Unknown_1",
  "Gut",
  "Notocord",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "Blood_Progenitor_1",
  "Mixed_Mesoderm",
  "Mesenchyme",
  "Haematodenothelial_progenitors",
  "Nascent_Mesoderm",
  "Pharyngeal_Mesoderm",
  "Paraxial_Mesoderm",
  "Caudal_Mesoderm",
  "Somitic_Mesoderm",
  "ExE_Ectoderm_1",
  "ExE_Ectoderm_2",
  "Epiblast_PrimStreak")

pdf("output/seurat/DotPlot_SCT_top3_conserved_embryo_V2clust.pdf", width=18, height=5)
DotPlot(embryo.combined.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# save
## saveRDS(embryo.combined.sct, file = "output/seurat/embryo.combined.sct_V2clust.rds")
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct_V2clust.rds")




# Check some genes
DefaultAssay(embryo.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

## post 20231005 Conchi meeting
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_Aldh1a2_Cyp26a1_V2clust.pdf", width=10, height=13)
FeaturePlot(embryo.combined.sct, features = c("Aldh1a2", "Cyp26a1"), max.cutoff = 10, cols = c("grey", "red"), split.by = "condition")
dev.off()

## post 20240110 Conchi email
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_20240110geneList_V2clust.pdf", width=10, height=75)
FeaturePlot(embryo.combined.sct, features = c("Crabp2", "Hand1", "Rbp1", "Aldh1a2", "Cyp26a1", "Hoxa1", "Fgf8", "Nr2f2", "Tnnt2", "Actc1", "Myl6", "Myl7", "Myl3", "Tbx1", "Tbx5", "Mab21l2", "Lbh"), max.cutoff = 10, cols = c("grey", "red"), split.by = "condition")
dev.off()



# Compare WT and cYAPKO using SCPA ##########################################

### installation (in scRNAseqV1)
# devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
# devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
# devtools::install_github("jackbibby1/SCPA")
# --> Last fail so installed dependencies
# BiocManager::install(c("clustermole", "ComplexHeatmap", "singscore"))
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
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct_V2clust.rds")
DefaultAssay(embryo.combined.sct) <- "RNA" # Recommended 



# Import manually curated gene list as a list
#### Need to folllow this format
pathways <- msigdbr("Mus musculus", "C2") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise
pathways$REACTOME_SIGNALING_BY_RETINOIC_ACID$Genes
### Convert table to list of genes
Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

pathways <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt"),
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt"),
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt"),
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")
)

pathways$Hippo$Genes


# Compare pathway activity between condition within  
cell_types <- unique(embryo.combined.sct$cluster.annot)
embryo.combined.sct_split <- SplitObject(embryo.combined.sct, split.by = "condition")

# SCPA comparison
## keep all columns
scpa_out_all <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryo.combined.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryo.combined.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_all[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}
## keep only PAthway and qval column (Best for heatmap qval representation)
scpa_out <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryo.combined.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryo.combined.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}


saveRDS(scpa_out_all, file = "output/Pathway/scpa_out_eachCellTypes_ManualGeneLists.rds")

console_output <- capture.output(print(scpa_out_all))
writeLines(console_output, "output/Pathway/scpa_out_console_ManualGeneLists.txt")

## Save output as R object
# saveRDS(scpa_out, file = "output/Pathway/scpa_out_eachCellTypes.rds")
# scpa_out_all <- readRDS("output/Pathway/scpa_out_all_eachCellTypes.rds")


# plot output

# Code to save output for each cell type comparison
clusters = c(
  "Epiblast_PrimStreak", 
  "ExE_Ectoderm_1", 
  "Nascent_Mesoderm", 
  "Somitic_Mesoderm", 
  "Caudal_Mesoderm", 
  "Paraxial_Mesoderm", 
  "Pharyngeal_Mesoderm", 
  "Haematodenothelial_progenitors", 
  "Mesenchyme", 
  "Blood_Progenitor_1", 
  "Mixed_Mesoderm", 
  "Blood_Progenitor_2", 
  "Surface_Ectoderm", 
  "Gut", 
  "Unknown_1", 
  "Notocord", 
  "Unknown_2", 
  "Primordial_Germ_Cells", 
  "ExE_Ectoderm_2"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(embryo.combined.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  cYAPKO <- seurat_extract(embryo.combined.sct,
                           meta1 = "condition", value_meta1 = "cYAPKO",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_V2clust_ManualGeneLists_", cluster, ".txt")
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}



# BALOW CLEAN CODE; dotplot>enrichplot>heatmap>violinplot
# DOT PLOT pepresetnation
## load all the comparison for each cell type (FC qval information)
clusters = c(
  "Epiblast_PrimStreak", 
  "ExE_Ectoderm_1", 
  "Nascent_Mesoderm", 
  "Somitic_Mesoderm", 
  "Caudal_Mesoderm", 
  "Paraxial_Mesoderm", 
  "Pharyngeal_Mesoderm", 
  "Haematodenothelial_progenitors", 
  "Mesenchyme", 
  "Blood_Progenitor_1", 
  "Mixed_Mesoderm", 
  "Blood_Progenitor_2", 
  "Surface_Ectoderm", 
  "Gut", 
  "Unknown_1", 
  "Notocord", 
  "Unknown_2", 
  "Primordial_Germ_Cells"
)
## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_V2clust_ManualGeneLists_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster))
# --> import the gene pathways (used table from Conchi `Pathwyas of interest*.xlsx`)
pathways = c("Hippo", "Nodal_TGF", "WNT", "Vitamin_A")
all_data_pathways = as_tibble(all_data) %>% 
  filter(Pathway %in% pathways)

# tidy the df
all_data_pathways_tidy <- all_data_pathways %>%
  filter(qval >= 1.4) 



# REFINE COLOR
custom_color <- function(fc_value){
  ifelse(fc_value >= 5, "strong_blue", 
         ifelse(fc_value > 2, "light_blue",
                ifelse(fc_value <= -5, "strong_red", 
                       ifelse(fc_value < -2, "light_red", "grey"))))
}

# Add a column for this custom color
all_data_pathways_tidy <- all_data_pathways_tidy %>%
  mutate(custom_col = sapply(FC, custom_color))

pdf("output/Pathway/dotplot_V2clust_ManualGeneLists_Hippo_FCtresh.pdf", width=10, height=4)
ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, color = custom_col), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = c("strong_blue" = "dodgerblue4",
                                "light_blue" = "lightblue2",
                                "grey" = "grey",
                                "light_red" = "indianred1",
                                "strong_red" = "red3")) +
  theme_bw() +
  labs(size = "q-value", color = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



pdf("output/Pathway/dotplot_V2clust_ManualGeneLists_grey.pdf", width=10, height=4)
ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, fill = "black"), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  labs(size = "q-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# heatmap

Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

pathways <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt"),
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt"),
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt"),
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")
)

pathways$Hippo$Genes



genes_of_interest <- pathways$Hippo$Genes
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(embryo.combined.sct@assays$RNA@data)]
### extract WT and cYAPKO gene expression values from RNA assay
WT_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "WT" & embryo.combined.sct$cluster.annot == "Unknown_1"]]
cYAPKO_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO" & embryo.combined.sct$cluster.annot == "Unknown_1"]]
### mean expression values for each gene
WT_mean <- rowMeans(WT_expression)
cYAPKO_mean <- rowMeans(cYAPKO_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  WT = WT_mean,
  cYAPKO = cYAPKO_mean
) %>% 
pivot_longer(cols = c(WT, cYAPKO), names_to = "Condition", values_to = "Expression")
## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(WT_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)




# if few genes:
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Blood_Progenitor_1.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Caudal_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Epiblast_PrimStreak.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Haematodenothelial_progenitors.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Mesenchyme.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Mixed_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Nascent_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Paraxial_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Pharyngeal_Mesoderm.pdf", width=6, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Nodal_TGF_Somitic_Mesoderm.pdf", width=6, height=2)


library("viridis")
ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()

# if many genesL
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Blood_Progenitor_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Blood_Progenitor_2.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Caudal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Epiblast_PrimStreak.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_ExE_Ectoderm_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Haematodenothelial_progenitors.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Mesenchyme.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Mixed_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Somitic_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_WNT_Surface_Ectoderm.pdf", width=3, height=2)

pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Blood_Progenitor_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Caudal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Epiblast_PrimStreak.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_ExE_Ectoderm_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Haematodenothelial_progenitors.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Mesenchyme.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Mixed_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_Vitamin_A_Somitic_Mesoderm.pdf", width=3, height=2)

pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Blood_Progenitor_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Blood_Progenitor_2.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Caudal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Epiblast_PrimStreak.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_ExE_Ectoderm_1.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Haematodenothelial_progenitors.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Mesenchyme.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Mixed_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Nascent_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Paraxial_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Pharyngeal_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Somitic_Mesoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Surface_Ectoderm.pdf", width=3, height=2)
pdf("output/Pathway/heatmap_embryo_RNAexpression_V2clust_ManualGeneLists_hippo_Unknown_1.pdf", width=3, height=2)


ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") 
dev.off()








# GSEA plot
library("fgsea")
## Re-calculate DEGs keeping ALL genes
embryo.combined.sct$celltype.stim <- paste(embryo.combined.sct$cluster.annot, embryo.combined.sct$condition,
    sep = "-")
Idents(embryo.combined.sct) <- "celltype.stim"

Paraxial_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Paraxial_Mesoderm-cYAPKO", ident.2 = "Paraxial_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Nascent_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Nascent_Mesoderm-cYAPKO", ident.2 = "Nascent_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Pharyngeal_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Pharyngeal_Mesoderm-cYAPKO", ident.2 = "Pharyngeal_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Mesenchyme <- FindMarkers(embryo.combined.sct, ident.1 = "Mesenchyme-cYAPKO", ident.2 = "Mesenchyme-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
Mixed_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Mixed_Mesoderm-cYAPKO", ident.2 = "Mixed_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
Blood_Progenitor_1 <- FindMarkers(embryo.combined.sct, ident.1 = "Blood_Progenitor_1-cYAPKO", ident.2 = "Blood_Progenitor_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Epiblast_PrimStreak <- FindMarkers(embryo.combined.sct, ident.1 = "Epiblast_PrimStreak-cYAPKO", ident.2 = "Epiblast_PrimStreak-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
ExE_Ectoderm_1 <- FindMarkers(embryo.combined.sct, ident.1 = "ExE_Ectoderm_1-cYAPKO", ident.2 = "ExE_Ectoderm_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Somitic_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Somitic_Mesoderm-cYAPKO", ident.2 = "Somitic_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Caudal_Mesoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Caudal_Mesoderm-cYAPKO", ident.2 = "Caudal_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Haematodenothelial_progenitors <- FindMarkers(embryo.combined.sct, ident.1 = "Haematodenothelial_progenitors-cYAPKO", ident.2 = "Haematodenothelial_progenitors-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Blood_Progenitor_2 <- FindMarkers(embryo.combined.sct, ident.1 = "Blood_Progenitor_2-cYAPKO", ident.2 = "Blood_Progenitor_2-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Surface_Ectoderm <- FindMarkers(embryo.combined.sct, ident.1 = "Surface_Ectoderm-cYAPKO", ident.2 = "Surface_Ectoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Gut <- FindMarkers(embryo.combined.sct, ident.1 = "Gut-cYAPKO", ident.2 = "Gut-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Unknown_1 <- FindMarkers(embryo.combined.sct, ident.1 = "Unknown_1-cYAPKO", ident.2 = "Unknown_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Unknown_2 <- FindMarkers(embryo.combined.sct, ident.1 = "Unknown_2-cYAPKO", ident.2 = "Unknown_2-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Notocord <- FindMarkers(embryo.combined.sct, ident.1 = "Notocord-cYAPKO", ident.2 = "Notocord-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Primordial_Germ_Cells <- FindMarkers(embryo.combined.sct, ident.1 = "Primordial_Germ_Cells-cYAPKO", ident.2 = "Primordial_Germ_Cells-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  





### save output
#### write.table(Primordial_Germ_Cells, file = "output/seurat/Primordial_Germ_Cells-cYAPKO_response_V2clust_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#### import all clsuter DEGs output :
cluster_types <- c("Epiblast_PrimStreak", "ExE_Ectoderm_1", "Nascent_Mesoderm", 
                   "Somitic_Mesoderm", "Caudal_Mesoderm", "Paraxial_Mesoderm", 
                   "Pharyngeal_Mesoderm", "Haematodenothelial_progenitors", 
                   "Mesenchyme", "Blood_Progenitor_1", "Mixed_Mesoderm", 
                   "Blood_Progenitor_2", "Surface_Ectoderm", "Gut", "Unknown_1", 
                   "Notocord", "Unknown_2", "Primordial_Germ_Cells")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-cYAPKO_response_V2clust_allGenes.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}


## load list of genes to test
Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

fgsea_sets <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt")$Genes,
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt")$Genes,
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt")$Genes,
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")$Genes
)





## Rank genes based on FC
genes <- Blood_Progenitor_1 %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(ES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -NES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()






## plot GSEA
pdf("output/Pathway/GSEA_V2clust_ManualGeneLists_Hippo-Primordial_Germ_Cells.pdf", width=5, height=3)

plotEnrichment(fgsea_sets[["Hippo"]],
               ranks) + labs(title="Hippo-Primordial_Germ_Cells") +
               theme_bw()
dev.off()





# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c("Epiblast_PrimStreak", "ExE_Ectoderm_1", "Nascent_Mesoderm", 
                   "Somitic_Mesoderm", "Caudal_Mesoderm", "Paraxial_Mesoderm", 
                   "Pharyngeal_Mesoderm", "Haematodenothelial_progenitors", 
                   "Mesenchyme", "Blood_Progenitor_1", "Mixed_Mesoderm", 
                   "Blood_Progenitor_2", "Surface_Ectoderm", "Gut", "Unknown_1", 
                   "Notocord", "Unknown_2", "Primordial_Germ_Cells")

## Initialize an empty list to store the results for each cluster type
all_results <- list()
## Loop over each cluster type
for (cluster in cluster_types) {
  
  # Extract genes for the current cluster
  genes <- get(cluster) %>% 
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC)
  
  ranks <- deframe(genes)
  
  # Run GSEA for the current cluster
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(ES))
  
  # Extract summary table and add cluster column
  fgseaResTidy_summary = fgseaResTidy %>% 
    dplyr::select(pathway, pval, padj, ES, size, NES) %>%
    mutate(cluster = cluster) %>%
    arrange(padj) %>% 
    head()
  
  # Store results in the list
  all_results[[cluster]] <- fgseaResTidy_summary
}
## Combine results from all cluster types into one table
final_results <- bind_rows(all_results, .id = "cluster")


write.table(final_results, file = c("output/Pathway/gsea_output_V2clust_ManualGeneLists.txt"), sep = "\t", quote = FALSE, row.names = FALSE)



# Heatmap all GSEA


pdf("output/Pathway/heatmap_gsea_padj.pdf", width=5, height=5)
ggplot(final_results, aes(x=cluster, y=pathway, fill=ES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Enrichment\nScore") +
  geom_text(aes(label=sprintf("%.2f", ES)), 
            color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()



pdf("output/Pathway/heatmap_gsea_pval_greyTile.pdf", width=5, height=5)
ggplot(final_results, aes(x=cluster, y=pathway)) + 
  geom_tile(aes(fill = ifelse(pval <= 0.05, ES, NA)), color = "black") +  # Conditional fill based on significance
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Enrichment\nScore", na.value="grey") +
  # geom_text(aes(label=sprintf("%.2f", ES)), 
  #           color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change between pvalue, qvalue,p.adjust
  #           size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1

dev.off()




# SCPA and GSEA plot
gsea_output = final_results %>% 
  rename("pathway" = "Pathway") %>%
  dplyr::select(Pathway, pval, padj, ES, cluster)
scpa_output = all_data_pathways_tidy  %>%
  dplyr::select(Pathway, qval, cluster)

scpa_gsea = gsea_output %>%
  left_join(scpa_output)

## Filter data for GSEA padj < 0.05
filtered_data <- scpa_gsea
## Assign color based on your criteria
filtered_data$is_significant <- filtered_data$padj < 0.05

pdf("output/Pathway/dotplot_scpa_gsea_V2clust.pdf", width=10, height=4)

ggplot(filtered_data, aes(x=cluster, y=Pathway, size=qval)) +
  geom_point(aes(color=ifelse(is_significant, ES, NA), 
                 shape=is_significant), alpha=0.7) +
  scale_color_gradient2(low="blue", high="red", midpoint=0, na.value="black", name = "GSEA Enrichment score") +
  scale_size_continuous(name="SCPA qval") +  # Title for size scale
  scale_shape_manual(values=c(`TRUE`=16, `FALSE`=16)) + # Using circle shape for both
  theme_light() +
  labs(title="SCPA and GSEA", x="Cluster", y="Pathway") +
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

```

























# Shinny app


Created with [ShinyCell](https://github.com/SGDDNB/ShinyCell); and follow [shinyapps](https://www.shinyapps.io/) to put it online 

```bash
conda activate scRNAseqV2
# conda install -c anaconda hdf5
```

```R
# installation
## devtools::install_github("SGDDNB/ShinyCell")
## install.packages('rsconnect')


# Packages
library("Seurat")
library("ShinyCell")
library("rsconnect")

# Data import EMBRYO
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")
DefaultAssay(embryo.combined.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(embryo.combined.sct)

makeShinyApp(embryo.combined.sct, scConf, gene.mapping = TRUE,
             shiny.title = "Embryo_V1",
             shiny.dir = "shinyApp_embryo_V1/") 

rsconnect::deployApp('shinyApp_embryo_V1')


# Generate Shiny app V2 new clustering Conchi meeting 20231005
scConf = createConfig(embryo.combined.sct)

makeShinyApp(embryo.combined.sct, scConf, gene.mapping = TRUE,
             shiny.title = "Embryo_V2",
             shiny.dir = "shinyApp_embryo_V2/") 

rsconnect::deployApp('shinyApp_embryo_V2')



# Data import HUMAN
humangastruloid.combined.sct <- readRDS(file = "output/seurat/humangastruloid.combined.sct_V2.rds")
DefaultAssay(humangastruloid.combined.sct) <- "RNA" # Recommended 


# Generate Shiny app
scConf = createConfig(humangastruloid.combined.sct)

makeShinyApp(humangastruloid.combined.sct, scConf, gene.mapping = TRUE,
             shiny.title = "Humangastruloid_V1",
             shiny.dir = "shinyApp_humangastruloid_V1/") 

rsconnect::deployApp('shinyApp_humangastruloid_V1')




```

--> Shiny apps:
- https://roulethomas.shinyapps.io/shinyapp_embryo_v1/ 
- https://roulethomas.shinyapps.io/shinyapp_humangastruloid_V1/ 


# Pseudotime

Let's use Monocle3 for pseudotime analsyis; some tutorial here:
- https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
- http://cole-trapnell-lab.github.io/monocle-release/docs/#recommended-analysis-protocol


And I'm gonna follow the [Bioinformagician stuff](https://www.youtube.com/watch?v=iq4T_uzMFcY)

--> Installation was painfull starting with scRNAseqV2 so let's create a new separated conda environment for Monocle as `monocle`

----> Only things that work was installation through mamba (help found [here](https://bioconda.github.io/recipes/r-monocle3/README.html))

```bash
conda create -n monocle3

conda install -c conda-forge mamba

conda config --add channels bioconda
conda config --add channels conda-forge

mamba create -n monocle3 r-base=4.1 r-monocle3

   mamba install r-monocle3
and update with::
   mamba update r-monocle3

## THEN JUST
conda activate monocle3

## new one to install SeuratWrappers
conda activate monocle3_V1

conda install -c bu_cnio r-seuratwrappers 
```

I needed to install in additino in R Seurat (install.packages) and R.utils and install.packages("remotes") and  seurat-wrappers and and remotes::install_github('satijalab/seurat-wrappers') and 

packageurl <- "https://cran.r-project.org/src/contrib/Archive/igraph/igraph_1.4.3.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

```bash
conda create --name monocle3_V1 --clone monocle3
conda activate monocle3_V1
```




Here is a 1st version using the integrated assay; not sure that is the best as the DEGs then failed...

```R
# Installation
library("igraph") # important to load it first so that v1.4.3 and not v1.5 is used (v1.5 is loaded via Seurat per default!)
library("monocle3")
library("Seurat")
library("SeuratWrappers")

# Data import EMBRYO
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct_V2clust.rds")

DefaultAssay(embryo.combined.sct) <- "integrated" # According to ucDAvis tutorial I think should be in integrative mode

# convert data to seurat object to cell_data_set
cds <- as.cell_data_set(embryo.combined.sct)
cds <- cluster_cells(cds, resolution=1e-3) # also tries  1e-5 and 1e-1 but WEIRD


pdf("output/monocle3/plot_cells_cluster_embryo.pdf", width=5, height=5)
pdf("output/monocle3/plot_cells_cluster_res1e5_embryo.pdf", width=5, height=5)
pdf("output/monocle3/plot_cells_cluster_res1e1_embryo.pdf", width=5, height=5)
pdf("output/monocle3/plot_cells_cluster_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
dev.off()

pdf("output/monocle3/plot_cells_partition_embryo.pdf", width=5, height=5)
pdf("output/monocle3/plot_cells_partition_res1e5_embryo.pdf", width=5, height=5)
pdf("output/monocle3/plot_cells_partition_res1e1_embryo.pdf", width=5, height=5)
pdf("output/monocle3/plot_cells_partition_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
dev.off()


# subsetting partition
integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

# Trajectory analysis
## RAW
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

pdf("output/monocle3/plot_cells_trajectory_partition1_embryo.pdf", width=5, height=5)
pdf("output/monocle3/plot_cells_trajectory_partition1_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()

## COLORED BY PSEUDOTIME

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 1])) # for V1
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 3])) # for V2cvlust
pdf("output/monocle3/plot_cells_trajectory_partition1_Root1_embryo.pdf", width=5, height=5)
pdf("output/monocle3/plot_cells_trajectory_partition1_Root1_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")
dev.off()


## PLOT as seurat object
integrated.sub <- as.Seurat(cds, assay = NULL)
pdf("output/monocle3/FeaturePlot_trajectory_partition1_Root1_embryo.pdf", width=5, height=5)
pdf("output/monocle3/FeaturePlot_trajectory_partition1_Root1_embryo_V2clust.pdf", width=5, height=5)
FeaturePlot(integrated.sub, "monocle3_pseudotime")
dev.off()


# Pseudotime differential genes
cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 1)


## Save the output table
write.table(cds_graph_test_results, file = "output/monocle3/cds_graph_test_results.txt", sep="\t", row.names=TRUE, quote=FALSE)





```

--> All information from WT and cYAPKO is used; I think that is OK as we want to see pseudotime no comparison between gentoypes?

--> I think I should NOT use the ibntegrative assay; otherwise it result VERY FEW genes only 3000 from the integration at the DEGs step... With RNA I will have all genes



2nd try using RNA and separate WT and cYAPKO
--> Then we could do Venn diagram to check the gene specific pseudotime-diff in WT and cYPAKO and plot their expression along pseudotime in x to stregnthen differences!!


```R
# Installation
library("igraph") # important to load it first so that v1.4.3 and not v1.5 is used (v1.5 is loaded via Seurat per default!)
library("monocle3")
library("Seurat")
library("SeuratWrappers")

# Data import EMBRYO
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")
DefaultAssay(embryo.combined.sct) <- "RNA" # According to ucDAvis tutorial I think should be in integrative mode

# subset WT and cYPAKO

embryo.combined.sct.WT <- subset(embryo.combined.sct, subset = condition == "WT")
embryo.combined.sct.cYAPKO <- subset(embryo.combined.sct, subset = condition == "cYAPKO")


# convert data to seurat object to cell_data_set
cds.WT <- as.cell_data_set(embryo.combined.sct.WT)
cds.WT <- cluster_cells(cds.WT, resolution=1e-3) # also tries  1e-5 and 1e-1 but WEIRD

pdf("output/monocle3/plot_cells_RNA_cluster_embryo_WT.pdf", width=5, height=5)
plot_cells(cds.WT, color_cells_by = "cluster", show_trajectory_graph = FALSE)
dev.off()

pdf("output/monocle3/plot_cells_RNA_partition_embryo_WT.pdf", width=5, height=5)
plot_cells(cds.WT, color_cells_by = "partition", show_trajectory_graph = FALSE)
dev.off()



cds.cYAPKO <- as.cell_data_set(embryo.combined.sct.cYAPKO)
cds.cYAPKO <- cluster_cells(cds.cYAPKO, resolution=1e-3) # also tries  1e-5 and 1e-1 but WEIRD

pdf("output/monocle3/plot_cells_RNA_cluster_embryo_cYAPKO.pdf", width=5, height=5)
plot_cells(cds.cYAPKO, color_cells_by = "cluster", show_trajectory_graph = FALSE)
dev.off()

pdf("output/monocle3/plot_cells_RNA_partition_embryo_cYAPKO.pdf", width=5, height=5)
plot_cells(cds.cYAPKO, color_cells_by = "partition", show_trajectory_graph = FALSE)
dev.off()




```

--> Not clear how to handle WT and cYAPKO as with the same resolution the partition is different; thus the trajectory will not be study on the same group of cells....


## monocle3 with RNA assay - E7.75
Here is a 3rd version using RNA assay; `conda activate monocle3_V1`



```R
# Installation
library("igraph") # important to load it first so that v1.4.3 and not v1.5 is used (v1.5 is loaded via Seurat per default!)
library("monocle3")
library("Seurat")
library("SeuratWrappers")

# Data import EMBRYO
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")
DefaultAssay(embryo.combined.sct) <- "RNA" # According to ucDAvis tutorial I think should be in integrative mode

# convert data to seurat object to cell_data_set
cds <- as.cell_data_set(embryo.combined.sct)
cds <- cluster_cells(cds, resolution=1e-3) # also tries  1e-5 and 1e-1 but WEIRD


pdf("output/monocle3/plot_cells_RNA_cluster_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
dev.off()

pdf("output/monocle3/plot_cells_RNA_partition_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
dev.off()


# subsetting partition
integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

# Trajectory analysis
## RAW
cds <- learn_graph(cds, use_partition = TRUE, verbose = TRUE)

pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()

## COLORED BY PSEUDOTIME

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 3]))
pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_Root1_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")
dev.off()


## PLOT as seurat object
integrated.sub <- as.Seurat(cds, assay = NULL)
pdf("output/monocle3/FeaturePlot_RNA_trajectory_partition1_Root1_embryo_V2clust.pdf", width=5, height=5)
FeaturePlot(integrated.sub, "monocle3_pseudotime")
dev.off()


# Pseudotime differential genes (1 hour!!!)
cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)


## Save the output table
# write.table(cds_graph_test_results, file = "output/monocle3/cds_graph_test_results_RNA_V2clust.txt", sep="\t", row.names=TRUE, quote=FALSE)
# cds_graph_test_results <- read.table("output/monocle3/cds_graph_test_results_RNA.txt", header = TRUE, sep = "\t", row.names = 1, quote = "")

# Check expression of some genes
rowData(cds)$gene_short_name <- row.names(rowData(cds))
head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_Root1_TopDEGsGenes_embryo.pdf", width=5, height=5)
plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)
dev.off()

pdf("output/monocle3/FeaturePlot_RNA_trajectory_partition1_Root1_TopDEGsGenes_embryo.pdf", width=10, height=12)
FeaturePlot(integrated.sub, features = head(deg_ids), max.cutoff = 10, cols = c("grey", "red"))
dev.off()  


## Plot gene as a function of pseudotime
#### Import DEGs
cds_graph_test_results <- read.table("output/monocle3/cds_graph_test_results_RNA.txt", header = TRUE, sep = "\t", row.names = 1, quote = "")
rowData(cds)$gene_short_name <- row.names(rowData(cds))
head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))
### Identify top 10 morans_I genes
ordered_results <- cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE), ]
top_morans_I <- rownames(ordered_results)[1:10]

partition_1_cds <- as.cell_data_set(integrated.sub, assay = "RNA")
partition_1_cds <- estimate_size_factors(partition_1_cds)


partition_1_cds <- cds[rowData(cds)$gene_short_name %in% top_morans_I,
                       colData(cds)$monocle3_partitions %in% c("1")]
partition_1_cds <- order_cells(partition_1_cds)

partition_1_cds <- partition_1_cds[,Matrix::colSums(exprs(partition_1_cds)) != 0]
partition_1_cds <- estimate_size_factors(partition_1_cds)


pdf("output/monocle3/plot_genes_in_pseudotime_trajectory_partition1_top_morans_I_embryo.pdf", width=10, height=12)
plot_genes_in_pseudotime(partition_1_cds,
                         color_cells_by="cluster.annot",
                         min_expr=0.5)
dev.off()

### Identify top 10 qvalue genes
ordered_results <- cds_graph_test_results[order(cds_graph_test_results$q_value, decreasing = FALSE), ]
top_qvalue <- rownames(ordered_results)[1:10]

partition_1_cds <- as.cell_data_set(integrated.sub, assay = "RNA")
partition_1_cds <- estimate_size_factors(partition_1_cds)


partition_1_cds <- cds[rowData(cds)$gene_short_name %in% top_qvalue,
                       colData(cds)$monocle3_partitions %in% c("1")]
partition_1_cds <- order_cells(partition_1_cds)

partition_1_cds <- partition_1_cds[,Matrix::colSums(exprs(partition_1_cds)) != 0]
partition_1_cds <- estimate_size_factors(partition_1_cds)


pdf("output/monocle3/plot_genes_in_pseudotime_trajectory_partition1_top_qvalue_embryo.pdf", width=10, height=12)
plot_genes_in_pseudotime(partition_1_cds,
                         color_cells_by="cluster.annot",
                         min_expr=0.5)
dev.off()


# GEnerate heatmap



```

--> Using RNA assay the clustering is slightly better; more in agreement with what we have in our previous Seurat UMAP

--> We got this time DEGs genes!

--> There is in the [UCDavis tutorial](https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed) a part where it generate heatmap with module of co-expressed genes but not sure that is necessary for us + it FAIL!

--> cds_graph_test_results output check morans_I column which is between 0 and 1 and 1 = positive correlation with pseudotime (help [here](https://cole-trapnell-lab.github.io/monocle3/docs/differential/) )
----> The one with the highest qvalue does not seems to be the most interested; the high morans_I look more responsive to the pseudotime



## monocle3 with RNA assay - integration E7 and E7.75
Here is a 3rd version using RNA assay; `conda activate monocle3_V1`



```R
# Installation
library("igraph") # important to load it first so that v1.4.3 and not v1.5 is used (v1.5 is loaded via Seurat per default!)
library("monocle3")
library("Seurat")
library("SeuratWrappers")

# Data import EMBRYO
embryo.combined.E7E775.sct <- readRDS(file = "output/seurat/embryo.combined.E7E775.sct.rds")
DefaultAssay(embryo.combined.E7E775.sct) <- "RNA" # According to ucDAvis tutorial I think should be in integrative mode

# convert data to seurat object to cell_data_set
cds <- as.cell_data_set(embryo.combined.E7E775.sct)
cds <- cluster_cells(cds, resolution=1e-3) # also tries  1e-5 show one big cluster with all cells


pdf("output/monocle3/plot_cells_RNA_cluster_embryo_E7E775.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
dev.off()

pdf("output/monocle3/plot_cells_RNA_partition_embryo_E7E775.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
dev.off()


# subsetting partition
integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1) # change partition/trajectory
cds <- as.cell_data_set(integrated.sub)

# Trajectory analysis
## RAW
cds <- learn_graph(cds, use_partition = TRUE, verbose = TRUE)

pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_embryo_E7E775.pdf", width=5, height=5)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()

## COLORED BY PSEUDOTIME

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 1])) # ROOT cluster
pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_Root1_embryo_E7E775.pdf", width=5, height=5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")
dev.off()


X
#################### BELOW NOT MODIFIED ####################
X


## PLOT as seurat object
integrated.sub <- as.Seurat(cds, assay = NULL)
pdf("output/monocle3/FeaturePlot_RNA_trajectory_partition1_Root1_embryo_V2clust.pdf", width=5, height=5)
FeaturePlot(integrated.sub, "monocle3_pseudotime")
dev.off()


# Pseudotime differential genes (1 hour!!!)
cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)


## Save the output table
# write.table(cds_graph_test_results, file = "output/monocle3/cds_graph_test_results_RNA_V2clust.txt", sep="\t", row.names=TRUE, quote=FALSE)
# cds_graph_test_results <- read.table("output/monocle3/cds_graph_test_results_RNA.txt", header = TRUE, sep = "\t", row.names = 1, quote = "")

# Check expression of some genes
rowData(cds)$gene_short_name <- row.names(rowData(cds))
head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_Root1_TopDEGsGenes_embryo.pdf", width=5, height=5)
plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)
dev.off()

pdf("output/monocle3/FeaturePlot_RNA_trajectory_partition1_Root1_TopDEGsGenes_embryo.pdf", width=10, height=12)
FeaturePlot(integrated.sub, features = head(deg_ids), max.cutoff = 10, cols = c("grey", "red"))
dev.off()  


## Plot gene as a function of pseudotime
#### Import DEGs
cds_graph_test_results <- read.table("output/monocle3/cds_graph_test_results_RNA.txt", header = TRUE, sep = "\t", row.names = 1, quote = "")
rowData(cds)$gene_short_name <- row.names(rowData(cds))
head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))
### Identify top 10 morans_I genes
ordered_results <- cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE), ]
top_morans_I <- rownames(ordered_results)[1:10]

partition_1_cds <- as.cell_data_set(integrated.sub, assay = "RNA")
partition_1_cds <- estimate_size_factors(partition_1_cds)


partition_1_cds <- cds[rowData(cds)$gene_short_name %in% top_morans_I,
                       colData(cds)$monocle3_partitions %in% c("1")]
partition_1_cds <- order_cells(partition_1_cds)

partition_1_cds <- partition_1_cds[,Matrix::colSums(exprs(partition_1_cds)) != 0]
partition_1_cds <- estimate_size_factors(partition_1_cds)


pdf("output/monocle3/plot_genes_in_pseudotime_trajectory_partition1_top_morans_I_embryo.pdf", width=10, height=12)
plot_genes_in_pseudotime(partition_1_cds,
                         color_cells_by="cluster.annot",
                         min_expr=0.5)
dev.off()

### Identify top 10 qvalue genes
ordered_results <- cds_graph_test_results[order(cds_graph_test_results$q_value, decreasing = FALSE), ]
top_qvalue <- rownames(ordered_results)[1:10]

partition_1_cds <- as.cell_data_set(integrated.sub, assay = "RNA")
partition_1_cds <- estimate_size_factors(partition_1_cds)


partition_1_cds <- cds[rowData(cds)$gene_short_name %in% top_qvalue,
                       colData(cds)$monocle3_partitions %in% c("1")]
partition_1_cds <- order_cells(partition_1_cds)

partition_1_cds <- partition_1_cds[,Matrix::colSums(exprs(partition_1_cds)) != 0]
partition_1_cds <- estimate_size_factors(partition_1_cds)


pdf("output/monocle3/plot_genes_in_pseudotime_trajectory_partition1_top_qvalue_embryo.pdf", width=10, height=12)
plot_genes_in_pseudotime(partition_1_cds,
                         color_cells_by="cluster.annot",
                         min_expr=0.5)
dev.off()


# GEnerate heatmap



```





## Condiments workflow to compare condition - E775

Let's try the Condiments workflow, specifically designed for this purpose:
- https://www.biorxiv.org/content/10.1101/2021.03.09.433671v1.full
- https://hectorrdb.github.io/condimentsPaper/ 



The installation of both Seurat and condiments was EXTREMELY PAINFUL! What work was to install 1st condiments then install `Seurat_v4.0.3` in `R\4.1.0` and then some packages individually from `tidyverse` ; as follow:

```bash
conda create --name condiments_V5 --clone condiments_V1 # clone condiments_V1 where only condiments is installed
conda activate condiments_V5
conda activate condiments_V6
```

```R
# package installation 
## install.packages("remotes")
## remotes::install_github("cran/spatstat.core")
## remotes::install_version("Seurat", "4.0.3")
## install.packages("magrittr")
## install.packages("magrittr")
## install.packages("dplyr")
## BiocManager::install("DelayedMatrixStats")
## BiocManager::install("tradeSeq")


# packages
library("condiments")
library("Seurat")
library("magrittr") # to use pipe
library("dplyr") # to use bind_cols and sample_frac
library("SingleCellExperiment") # for reducedDims
library("ggplot2")
library("slingshot")
library("DelayedMatrixStats")
library("tidyr")
library("tradeSeq")
library("cowplot")
library("scales")
library("pheatmap")

# Data import EMBRYO
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")
embryo.combined.sct <- readRDS(file = "output/seurat/embryo.combined.sct_V2clust.rds") 

DefaultAssay(embryo.combined.sct) <- "RNA" # According to condiments workflow

# convert to SingleCellExperiment
embryo <- as.SingleCellExperiment(embryo.combined.sct, assay = "RNA")


# tidy

df <- bind_cols(
  as.data.frame(reducedDims(embryo)$UMAP),
  as.data.frame(colData(embryo)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
## genotype overlap
pdf("output/condiments/UMAP_condition_embryo.pdf", width=5, height=5)
pdf("output/condiments/UMAP_condition_embryo_V2clust.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("blue", "red")) + # Specify colors here
  labs(col = "Genotypes") +
  theme_classic()
dev.off()

## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$condition,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

pdf("output/condiments/UMAP_imbalance_score_embryo.pdf", width=5, height=5)
pdf("output/condiments/UMAP_imbalance_score_embryo_V2clust.pdf", width=5, height=5)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores") +
  theme_classic()
dev.off()


#  Trajectory Inference and Differential Topology
set.seed(42)

## PLOT with Separate trajectories
embryo <- slingshot(embryo, reducedDim = 'UMAP',
                 clusterLabels = colData(embryo)$cluster.annot,
                 start.clus = 'Epiblast_PrimStreak', approx_points = 100)

set.seed(42)
topologyTest(SlingshotDataSet(embryo), embryo$condition) # KS_mean / 0.01 / 0.006288382 / 0.9789538


sdss <- slingshot_conditions(SlingshotDataSet(embryo), embryo$condition)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "condition")

pdf("output/condiments/UMAP_trajectory_separated_embryo.pdf", width=5, height=5)
pdf("output/condiments/UMAP_trajectory_separated_embryo_V2clust.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  theme_classic()
dev.off()

##### NEED TO MODIFY THE CODE BELOW TO ANNOTATE THE DIFFERENT TRAJECTORIES
pdf("output/condiments/UMAP_trajectory_separated_trajAnnotated_embryo.pdf", width=5, height=5)
pdf("output/condiments/UMAP_trajectory_separated_trajAnnotated_embryo_V2clust.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  annotate("text", x = -10, y = 6, label = "Lineage1", size = 5) +
  annotate("text", x = -7, y = -2.7, label = "Lineage2", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) +
  NULL
dev.off()
####### 






## PLOT with common trajectories
df_2 <- bind_cols(
  as.data.frame(reducedDim(embryo, "UMAP")),
  slingPseudotime(embryo) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(embryo) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst),
        # pst = max(pst) - pst)
)
curves <- slingCurves(embryo, as.df = TRUE)

pdf("output/condiments/UMAP_trajectory_common_embryo.pdf", width=5, height=5)
pdf("output/condiments/UMAP_trajectory_common_embryo_V2clust.pdf", width=5, height=5)

ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1.5) +
  theme_classic()
dev.off()


### With label
#### Calculate midpoint for each trajectory to place the label
curves_midpoints <- curves %>%
  group_by(Lineage) %>%
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))

pdf("output/condiments/UMAP_trajectory_common_label_embryo.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_midpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Add labels
  theme_classic()
dev.off()


curves_endpoints <- curves %>%
  group_by(Lineage) %>%
  arrange(Order) %>%
  top_n(1, Order) # Get the top/last ordered point for each group
pdf("output/condiments/UMAP_trajectory_common_label_embryo_V2clust.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_endpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
  theme_classic()
dev.off()


# Differential Progression
progressionTest(embryo, conditions = embryo$condition, lineages = TRUE)

prog_res <- progressionTest(embryo, conditions = embryo$condition, lineages = TRUE)


df_3 <-  slingPseudotime(embryo) %>% as.data.frame() 

df_3$condition <- embryo$condition
df_3 <- df_3 %>% 
  pivot_longer(-condition, names_to = "Lineage",
               values_to = "pst") %>%
  filter(!is.na(pst))

pdf("output/condiments/densityPlot_trajectory_lineages_embryo.pdf", width=10, height=5)
pdf("output/condiments/densityPlot_trajectory_lineages_embryo_V2clust.pdf", width=10, height=5)

ggplot(df_3, aes(x = pst)) +
  geom_density(alpha = .8, aes(fill = condition), col = "transparent") +
  geom_density(aes(col = condition), fill = "transparent", size = 1.5) +
  labs(x = "Pseudotime", fill = "condition") +
  facet_wrap(~Lineage, scales = "free", nrow=2) +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("blue", "red"))
  )) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()
dev.off()


#### ->  save.image(file="output/condiments/condiments_embryo_V2clust.RData")
### load("output/condiments/condiments_embryo_V2clust.RData")
set.seed(42)

#  Differential expression

## Identify the needed number of knots 

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 8
icMat <- evaluateK(counts = embryo, sds = SlingshotDataSet(embryo), 
                   conditions = factor(embryo$condition),
                   nGenes = 300, parallel = FALSE, BPPARAM = BPPARAM, k = 3:7) # set parallel = FALSE otherwise the code never end!
icMat
### !!! --> EXAMINE icMat to determine the optimal nb of knots for the GAM  !!! TOO LONG FUCK IT, use 6
## Fit GAM with the indicated nb of knots (4 to 7 works for most data according to developers) https://github.com/statOmics/tradeSeq/issues/54


# THE fitGAM() CODE IS TOO LONG TO RUN, SO HAS BEEN INTRODUCED INTO A RSCRIPT INSTEAD AND IMAGE AS BEEN SAVED AS condiments_embryo_V2.RData. But it fail too...


# Try running the DEGs trajectory per trajectory
### FROM THIS https://github.com/statOmics/tradeSeq/issues/64 :

#### Let's try to run the DEGs trajectory per trajectory
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,1]
cellweights <- slingCurveWeights(embryo) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj1 <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, 
     nknots = 6,
     sce = TRUE
   )
### IT WORK!!! Estimated run 24hours which is OK !!! Let's run this in slurm job traj per traj


###### --------------
# Try running the DEGs on the top 10k features (the 10k most variable genes) and trajectory per trajectory; This does not reduce a lot the time... With 10 featuers = 5hours; versus 24hours for all...

## Isolate the top 10k
embryo.combined.sct.10k <- FindVariableFeatures(embryo.combined.sct, selection.method = "vst", nfeatures = 10)

#### Let's try to run the DEGs trajectory per trajectory
counts <- embryo.combined.sct.10k[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct.10k$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,1]
cellweights <- slingCurveWeights(embryo) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj1 <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, 
     nknots = 6,
     sce = TRUE
   )

###### --------------


### Once slurm jobs finished for each trajectory; load traj1 <- readRDS("output/condiments/traj1.rds") OK

traj1 <- readRDS("output/condiments/traj1.rds")
traj2 <- readRDS("output/condiments/traj2.rds")
traj3 <- readRDS("output/condiments/traj3.rds")
traj4 <- readRDS("output/condiments/traj4.rds")
traj5 <- readRDS("output/condiments/traj5.rds")
traj6 <- readRDS("output/condiments/traj6.rds")
traj7 <- readRDS("output/condiments/traj7.rds")
traj8 <- readRDS("output/condiments/traj8.rds")

## DEGs between condition
condRes_traj5 <- conditionTest(traj5)
condRes_traj8_l2fc2 <- conditionTest(traj8, l2fc = log2(2)) # 

condRes_traj1_l2fc4 <- conditionTest(traj1, l2fc = log2(4)) # let s prefer to use this one
condRes_traj2_l2fc4 <- conditionTest(traj2, l2fc = log2(4)) # let s prefer to use this one




# Correct the pvalue with fdr

condRes_traj2_l2fc4$padj <- p.adjust(condRes_traj2_l2fc4$pvalue, "fdr")



### Save output tables
condRes_traj2_l2fc4$gene <- rownames(condRes_traj2_l2fc4) # create new column l;abel gene; as matrix before
condRes_traj2_l2fc4 <- condRes_traj2_l2fc4[, c(ncol(condRes_traj2_l2fc4), 1:(ncol(condRes_traj2_l2fc4)-1))] # just to put gene column 1st
write.table(condRes_traj2_l2fc4, file = c("output/condiments/condRes_traj2_l2fc4.txt"),sep="\t", quote=FALSE, row.names=FALSE)





# Heatmap clutering DEGs per traj _ REVISED METHOD
## import DEGs
condRes_traj1 <- read.table("output/condiments/condRes_traj1_l2fc4.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
condRes_traj1 <- read.table("output/condiments/condRes_traj1_l2fc2.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

## Isolate significant DEGs and transform into a vector
conditionGenes_traj1_vector <- condRes_traj1 %>% 
  filter(padj <= 0.05) %>%
  pull(gene)
  
# Predict smoothed values
yhatSmooth <- 
  predictSmooth(traj1, gene = conditionGenes_traj1_vector, nPoints = 50, tidy = FALSE) %>%
  log1p()

yhatSmoothScaled <- t(apply(yhatSmooth, 1, scales::rescale))
combinedData <- yhatSmoothScaled[, c(51:100, 1:50)]
# Generate heatmap with clustering
# Perform hierarchical clustering
hc <- hclust(dist(combinedData))
clusters <- cutree(hc, k=10)
# Create an annotation data frame for the rows based on cluster assignments
annotation_row <- data.frame(Cluster = factor(clusters))
# Define colors for each cluster
# 20
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen", "grey", "darkred", "darkblue", "gold", "darkgray", "lightblue", "lightgreen", "lightcoral", "lightpink", "lightcyan"))(20),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 10
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen" ))(10),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 8
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown" ))(8),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 7
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink" ))(7),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 6
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange"))(6),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 5
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple"))(5),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 4
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow" ))(4),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Generate the heatmap
pdf("output/condiments/clustered_heatmap_traj1.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_l2fc4_cl10.pdf", width=8, height=10)

pdf("output/condiments/clustered_heatmap_traj4_l2fc2_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj5_l2fc2_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj3_l2fc2_cl7.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_l2fc2_cl10.pdf", width=8, height=10)

pheatmap(combinedData,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Trajectory 1 - Hierarchical Clustering",
  legend = TRUE,
  cutree_rows = 10,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
)
dev.off()


# Line plots
library("reshape2")
library("stringr")
# Assuming yhatSmoothScaled contains your smoothed gene expression data
# Convert the yhatSmoothScaled data to a dataframe
df <- as.data.frame(yhatSmoothScaled)
df$Gene <- rownames(df)
# Transform the data into a long format
df_long <- melt(df, id.vars = "Gene", variable.name = "Pseudotime", value.name = "Expression")
# Attach the cluster information to the data frame
df$Cluster <- factor(clusters[df$Gene])
df_long$Cluster <- df$Cluster[match(df_long$Gene, df$Gene)]

# Extract condition column
df_long$Condition <- str_extract(df_long$Pseudotime, "condition[[:alnum:]]+")
df_long$Condition <- ifelse(str_detect(df_long$Condition, "WT"), "WT", "KO")

# Extract the point value and convert it to numeric
df_long$Updated_Pseudotime <- as.numeric(str_extract(df_long$Pseudotime, "(?<=point)\\d+"))

# Define colors for the conditions
color_map <- c("WT" = "blue", "KO" = "red")

# Plot using ggplot
pdf("output/condiments/clustered_linePlot_traj1_l2fc4_cl10.pdf", width=10, height=5)

pdf("output/condiments/clustered_linePlot_traj4_l2fc2_cl4.pdf", width=10, height=5)

ggplot(df_long, aes(x = as.numeric(Updated_Pseudotime), y = Expression, group = Gene)) + 
  geom_line(data = subset(df_long, Condition == "WT"), aes(color = Condition), alpha = 0.5) +
  geom_line(data = subset(df_long, Condition == "KO"), aes(color = Condition), alpha = 0.5) +
  scale_color_manual(values = color_map) + 
  facet_wrap(~Cluster, scales = "free_y", nrow=2) +
  theme_bw() +
  labs(title = "Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")

dev.off()

# Plot using ggplot
pdf("output/condiments/smoothed_linePlot_traj1_l2fc4_cl10_smooth.pdf", width=10, height=5)
pdf("output/condiments/smoothed_linePlot_traj2_l2fc2_cl10_smooth.pdf", width=10, height=5)

pdf("output/condiments/smoothed_linePlot_traj4_l2fc2_cl4_smooth.pdf", width=10, height=5)

ggplot(df_long, aes(x = Updated_Pseudotime, y = Expression, color = Condition)) + 
  geom_smooth(method = "loess", se = TRUE, span = 0.5) + 
  scale_color_manual(values = color_map) + 
  facet_wrap(~Cluster, scales = "free_y", nrow=2) +
  theme_bw() +
  labs(title = "Smoothed Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")
dev.off()




### Export gene list from each cluster
## Create a data frame with gene names and their respective cluster assignments
output_df <- data.frame(
  gene = rownames(combinedData),
  cluster = clusters
)

# Write the data frame to a .txt file
write.table(output_df, 
            file = "output/condiments/gene_clusters_traj1_l2fc2_cl10.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

# Check some genes individually
## FOR LINEAGE 1
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,1]
cellweights <- slingCurveWeights(embryo) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj1_Aff3.pdf", width=8, height=4)
plotSmoothers(traj1, sub_counts, gene = "Aff3", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()


## FOR LINEAGE 2
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,2]
cellweights <- slingCurveWeights(embryo) [,2]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj2_Gm28376.pdf", width=8, height=4)
plotSmoothers(traj2, sub_counts, gene = "Gm28376", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()

## FOR LINEAGE 3
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,3]
cellweights <- slingCurveWeights(embryo) [,3]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj3_Fcrla.pdf", width=8, height=4)
plotSmoothers(traj3, sub_counts, gene = "Fcrla", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()

## FOR LINEAGE 4
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,4]
cellweights <- slingCurveWeights(embryo) [,4]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj4_Hist2h2ac.pdf", width=8, height=4)
plotSmoothers(traj4, sub_counts, gene = "Hist2h2ac", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()



```

--> I did all the tutorial from [here](https://hectorrdb.github.io/condimentsPaper/articles/Fibrosis.html) 
----> This part Differential fate selection FAIL; but not sure it is usefull

--> prog_res? Lineage all pvalue signif! (paste in the ppt Conchi 20231005)

--> The icMAt is too long! I gave up running it. Author recommend using 6 knots (https://github.com/statOmics/tradeSeq/issues/121). Other discussion about [time](https://github.com/statOmics/tradeSeq/issues/41)

--> Same for the fitGAM, too long, so I run it into a Rscript as follow:
----> In the end; **what worked is to run trajectory per trajectory; it should last around 24-72hours per trajectory. Using the top 10k features do not reduce time significantly.**
----> Parrelization can be stuck sometime so I ll run it without...

--> The pseudotime diff looks weirk so let s use log2fc files anmd not the raw not log2fc filtered


```bash
conda activate condiments_V5
# test
sbatch scripts/fitGAM_6knots.sh # 5756977 fail
sbatch scripts/fitGAM_6knots_parralell_subset.sh # 6104946 fail

# trajectory per trajectory (all features, no parralelization)
sbatch scripts/fitGAM_6knots_traj1.sh # 6110340 ok
sbatch scripts/fitGAM_6knots_traj2.sh # 6110351 ok
sbatch scripts/fitGAM_6knots_traj3.sh # 6110359 ok
sbatch scripts/fitGAM_6knots_traj4.sh # 6110363 ok
sbatch scripts/fitGAM_6knots_traj5.sh # 6110370 ok
sbatch scripts/fitGAM_6knots_traj6.sh # 6110374 ok
sbatch scripts/fitGAM_6knots_traj7.sh # 6110376 ok
sbatch scripts/fitGAM_6knots_traj8.sh # 6110378 ok

```

--> `scripts/fitGAM_6knots.sh` 24hrs still running... FAIL time limit after 6 days LOL!!!! same for `fitGAM_6knots_parralell_subset.sh`

--> Without parralell processing trajectory per traj works great!! 24-72hrs to run

to check:
- https://www.bioconductor.org/packages/devel/bioc/vignettes/condiments/inst/doc/condiments.html
- For [plots](https://bioconductor.org/packages/devel/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html)




## Condiments workflow to compare condition - E7 and E775 integrated

```bash
conda activate condiments_V6
```

```R
# package installation 
## install.packages("remotes")
## remotes::install_github("cran/spatstat.core")
## remotes::install_version("Seurat", "4.0.3")
## install.packages("magrittr")
## install.packages("magrittr")
## install.packages("dplyr")
## BiocManager::install("DelayedMatrixStats")
## BiocManager::install("tradeSeq")


# packages
library("condiments")
library("Seurat")
library("magrittr") # to use pipe
library("dplyr") # to use bind_cols and sample_frac
library("SingleCellExperiment") # for reducedDims
library("ggplot2")
library("slingshot")
library("DelayedMatrixStats")
library("tidyr")
library("tradeSeq")
library("cowplot")
library("scales")
library("pheatmap")

# Data import EMBRYO
embryo.combined.E7E775.sct <- readRDS(file = "output/seurat/embryo.combined.E7E775.sct.rds") 

DefaultAssay(embryo.combined.E7E775.sct) <- "RNA" # According to condiments workflow

# convert to SingleCellExperiment
embryo <- as.SingleCellExperiment(embryo.combined.E7E775.sct, assay = "RNA")


# tidy

df <- bind_cols(
  as.data.frame(reducedDims(embryo)$UMAP),
  as.data.frame(colData(embryo)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
## genotype overlap
pdf("output/condiments/UMAP_condition_embryo_E7E775.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("blue", "red")) + # Specify colors here
  labs(col = "Genotypes") +
  theme_classic()
dev.off()

## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$condition,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

pdf("output/condiments/UMAP_imbalance_score_embryo_E7E775.pdf", width=5, height=5)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores") +
  theme_classic()
dev.off()


#  Trajectory Inference and Differential Topology
set.seed(42)

## PLOT with Separate trajectories
embryo <- slingshot(embryo, reducedDim = 'UMAP',
                 clusterLabels = colData(embryo)$cluster.annot,
                 start.clus = 'Epiblast_1', approx_points = 100)

### fine tune nb of trajectories
#### end.clusV1
embryo <- slingshot(embryo, reducedDim = 'UMAP',
                 clusterLabels = colData(embryo)$cluster.annot,
                 start.clus = 'Epiblast_1',
                 end.clus = c("Blood_Progenitor_2", "ExE_Ectoderm_2", "Endoderm", "Notocord", "Paraxial_Mesoderm", "Pharyngeal_Mesoderm", "Mesenchyme" ),
                 approx_points = 100)
###
set.seed(42)
topologyTest(SlingshotDataSet(embryo), embryo$condition) # KS_mean / 0.01 / 0.007581268 / 0.7581611


sdss <- slingshot_conditions(SlingshotDataSet(embryo), embryo$condition)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "condition")

# pdf("output/condiments/UMAP_trajectory_separated_embryo_E7E775.pdf", width=5, height=5)
pdf("output/condiments/UMAP_trajectory_separated_embryo_E7E775_end.clusV1.pdf", width=5, height=5)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  theme_classic()
dev.off()

##### NEED TO MODIFY THE CODE BELOW TO ANNOTATE THE DIFFERENT TRAJECTORIES
pdf("output/condiments/UMAP_trajectory_separated_trajAnnotated_embryo.pdf", width=5, height=5)
pdf("output/condiments/UMAP_trajectory_separated_trajAnnotated_embryo_V2clust.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  annotate("text", x = -10, y = 6, label = "Lineage1", size = 5) +
  annotate("text", x = -7, y = -2.7, label = "Lineage2", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) +
  NULL
dev.off()
####### 






## PLOT with common trajectories
df_2 <- bind_cols(
  as.data.frame(reducedDim(embryo, "UMAP")),
  slingPseudotime(embryo) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(embryo) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst),
        # pst = max(pst) - pst)
)
curves <- slingCurves(embryo, as.df = TRUE)

pdf("output/condiments/UMAP_trajectory_common_embryo_E7E775.pdf", width=5, height=5)

ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1.5) +
  theme_classic()
dev.off()


### With label
#### Calculate midpoint for each trajectory to place the label
curves_midpoints <- curves %>%
  group_by(Lineage) %>%
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))

pdf("output/condiments/UMAP_trajectory_common_label_embryo.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_midpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Add labels
  theme_classic()
dev.off()


curves_endpoints <- curves %>%
  group_by(Lineage) %>%
  arrange(Order) %>%
  top_n(1, Order) # Get the top/last ordered point for each group
pdf("output/condiments/UMAP_trajectory_common_label_embryo_V2clust.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_endpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
  theme_classic()
dev.off()


# Differential Progression
progressionTest(embryo, conditions = embryo$condition, lineages = TRUE)

prog_res <- progressionTest(embryo, conditions = embryo$condition, lineages = TRUE)


df_3 <-  slingPseudotime(embryo) %>% as.data.frame() 

df_3$condition <- embryo$condition
df_3 <- df_3 %>% 
  pivot_longer(-condition, names_to = "Lineage",
               values_to = "pst") %>%
  filter(!is.na(pst))

pdf("output/condiments/densityPlot_trajectory_lineages_embryo.pdf", width=10, height=5)
pdf("output/condiments/densityPlot_trajectory_lineages_embryo_V2clust.pdf", width=10, height=5)

ggplot(df_3, aes(x = pst)) +
  geom_density(alpha = .8, aes(fill = condition), col = "transparent") +
  geom_density(aes(col = condition), fill = "transparent", size = 1.5) +
  labs(x = "Pseudotime", fill = "condition") +
  facet_wrap(~Lineage, scales = "free", nrow=2) +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("blue", "red"))
  )) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()
dev.off()


#### ->  save.image(file="output/condiments/condiments_embryo_V2clust.RData")
### load("output/condiments/condiments_embryo_V2clust.RData")
set.seed(42)

#  Differential expression

## Identify the needed number of knots 

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 8
icMat <- evaluateK(counts = embryo, sds = SlingshotDataSet(embryo), 
                   conditions = factor(embryo$condition),
                   nGenes = 300, parallel = FALSE, BPPARAM = BPPARAM, k = 3:7) # set parallel = FALSE otherwise the code never end!
icMat
### !!! --> EXAMINE icMat to determine the optimal nb of knots for the GAM  !!! TOO LONG FUCK IT, use 6
## Fit GAM with the indicated nb of knots (4 to 7 works for most data according to developers) https://github.com/statOmics/tradeSeq/issues/54


# THE fitGAM() CODE IS TOO LONG TO RUN, SO HAS BEEN INTRODUCED INTO A RSCRIPT INSTEAD AND IMAGE AS BEEN SAVED AS condiments_embryo_V2.RData. But it fail too...


# Try running the DEGs trajectory per trajectory
### FROM THIS https://github.com/statOmics/tradeSeq/issues/64 :

#### Let's try to run the DEGs trajectory per trajectory
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,1]
cellweights <- slingCurveWeights(embryo) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj1 <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, 
     nknots = 6,
     sce = TRUE
   )
### IT WORK!!! Estimated run 24hours which is OK !!! Let's run this in slurm job traj per traj


###### --------------
# Try running the DEGs on the top 10k features (the 10k most variable genes) and trajectory per trajectory; This does not reduce a lot the time... With 10 featuers = 5hours; versus 24hours for all...

## Isolate the top 10k
embryo.combined.sct.10k <- FindVariableFeatures(embryo.combined.sct, selection.method = "vst", nfeatures = 10)

#### Let's try to run the DEGs trajectory per trajectory
counts <- embryo.combined.sct.10k[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct.10k$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,1]
cellweights <- slingCurveWeights(embryo) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj1 <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, 
     nknots = 6,
     sce = TRUE
   )

###### --------------


### Once slurm jobs finished for each trajectory; load traj1 <- readRDS("output/condiments/traj1.rds") OK

traj1 <- readRDS("output/condiments/traj1.rds")
traj2 <- readRDS("output/condiments/traj2.rds")
traj3 <- readRDS("output/condiments/traj3.rds")
traj4 <- readRDS("output/condiments/traj4.rds")
traj5 <- readRDS("output/condiments/traj5.rds")
traj6 <- readRDS("output/condiments/traj6.rds")
traj7 <- readRDS("output/condiments/traj7.rds")
traj8 <- readRDS("output/condiments/traj8.rds")

## DEGs between condition
condRes_traj5 <- conditionTest(traj5)
condRes_traj8_l2fc2 <- conditionTest(traj8, l2fc = log2(2)) # 

condRes_traj1_l2fc4 <- conditionTest(traj1, l2fc = log2(4)) # let s prefer to use this one
condRes_traj2_l2fc4 <- conditionTest(traj2, l2fc = log2(4)) # let s prefer to use this one




# Correct the pvalue with fdr

condRes_traj2_l2fc4$padj <- p.adjust(condRes_traj2_l2fc4$pvalue, "fdr")



### Save output tables
condRes_traj2_l2fc4$gene <- rownames(condRes_traj2_l2fc4) # create new column l;abel gene; as matrix before
condRes_traj2_l2fc4 <- condRes_traj2_l2fc4[, c(ncol(condRes_traj2_l2fc4), 1:(ncol(condRes_traj2_l2fc4)-1))] # just to put gene column 1st
write.table(condRes_traj2_l2fc4, file = c("output/condiments/condRes_traj2_l2fc4.txt"),sep="\t", quote=FALSE, row.names=FALSE)





# Heatmap clutering DEGs per traj _ REVISED METHOD
## import DEGs
condRes_traj1 <- read.table("output/condiments/condRes_traj1_l2fc4.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
condRes_traj1 <- read.table("output/condiments/condRes_traj1_l2fc2.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

## Isolate significant DEGs and transform into a vector
conditionGenes_traj1_vector <- condRes_traj1 %>% 
  filter(padj <= 0.05) %>%
  pull(gene)
  
# Predict smoothed values
yhatSmooth <- 
  predictSmooth(traj1, gene = conditionGenes_traj1_vector, nPoints = 50, tidy = FALSE) %>%
  log1p()

yhatSmoothScaled <- t(apply(yhatSmooth, 1, scales::rescale))
combinedData <- yhatSmoothScaled[, c(51:100, 1:50)]
# Generate heatmap with clustering
# Perform hierarchical clustering
hc <- hclust(dist(combinedData))
clusters <- cutree(hc, k=10)
# Create an annotation data frame for the rows based on cluster assignments
annotation_row <- data.frame(Cluster = factor(clusters))
# Define colors for each cluster
# 20
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen", "grey", "darkred", "darkblue", "gold", "darkgray", "lightblue", "lightgreen", "lightcoral", "lightpink", "lightcyan"))(20),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 10
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen" ))(10),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 8
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown" ))(8),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 7
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink" ))(7),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 6
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange"))(6),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 5
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple"))(5),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 4
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow" ))(4),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Generate the heatmap
pdf("output/condiments/clustered_heatmap_traj1.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_l2fc4_cl10.pdf", width=8, height=10)

pdf("output/condiments/clustered_heatmap_traj4_l2fc2_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj5_l2fc2_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj3_l2fc2_cl7.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_l2fc2_cl10.pdf", width=8, height=10)

pheatmap(combinedData,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Trajectory 1 - Hierarchical Clustering",
  legend = TRUE,
  cutree_rows = 10,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
)
dev.off()


# Line plots
library("reshape2")
library("stringr")
# Assuming yhatSmoothScaled contains your smoothed gene expression data
# Convert the yhatSmoothScaled data to a dataframe
df <- as.data.frame(yhatSmoothScaled)
df$Gene <- rownames(df)
# Transform the data into a long format
df_long <- melt(df, id.vars = "Gene", variable.name = "Pseudotime", value.name = "Expression")
# Attach the cluster information to the data frame
df$Cluster <- factor(clusters[df$Gene])
df_long$Cluster <- df$Cluster[match(df_long$Gene, df$Gene)]

# Extract condition column
df_long$Condition <- str_extract(df_long$Pseudotime, "condition[[:alnum:]]+")
df_long$Condition <- ifelse(str_detect(df_long$Condition, "WT"), "WT", "KO")

# Extract the point value and convert it to numeric
df_long$Updated_Pseudotime <- as.numeric(str_extract(df_long$Pseudotime, "(?<=point)\\d+"))

# Define colors for the conditions
color_map <- c("WT" = "blue", "KO" = "red")

# Plot using ggplot
pdf("output/condiments/clustered_linePlot_traj1_l2fc4_cl10.pdf", width=10, height=5)

pdf("output/condiments/clustered_linePlot_traj4_l2fc2_cl4.pdf", width=10, height=5)

ggplot(df_long, aes(x = as.numeric(Updated_Pseudotime), y = Expression, group = Gene)) + 
  geom_line(data = subset(df_long, Condition == "WT"), aes(color = Condition), alpha = 0.5) +
  geom_line(data = subset(df_long, Condition == "KO"), aes(color = Condition), alpha = 0.5) +
  scale_color_manual(values = color_map) + 
  facet_wrap(~Cluster, scales = "free_y", nrow=2) +
  theme_bw() +
  labs(title = "Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")

dev.off()

# Plot using ggplot
pdf("output/condiments/smoothed_linePlot_traj1_l2fc4_cl10_smooth.pdf", width=10, height=5)
pdf("output/condiments/smoothed_linePlot_traj2_l2fc2_cl10_smooth.pdf", width=10, height=5)

pdf("output/condiments/smoothed_linePlot_traj4_l2fc2_cl4_smooth.pdf", width=10, height=5)

ggplot(df_long, aes(x = Updated_Pseudotime, y = Expression, color = Condition)) + 
  geom_smooth(method = "loess", se = TRUE, span = 0.5) + 
  scale_color_manual(values = color_map) + 
  facet_wrap(~Cluster, scales = "free_y", nrow=2) +
  theme_bw() +
  labs(title = "Smoothed Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")
dev.off()




### Export gene list from each cluster
## Create a data frame with gene names and their respective cluster assignments
output_df <- data.frame(
  gene = rownames(combinedData),
  cluster = clusters
)

# Write the data frame to a .txt file
write.table(output_df, 
            file = "output/condiments/gene_clusters_traj1_l2fc2_cl10.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

# Check some genes individually
## FOR LINEAGE 1
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,1]
cellweights <- slingCurveWeights(embryo) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj1_Aff3.pdf", width=8, height=4)
plotSmoothers(traj1, sub_counts, gene = "Aff3", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()


## FOR LINEAGE 2
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,2]
cellweights <- slingCurveWeights(embryo) [,2]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj2_Gm28376.pdf", width=8, height=4)
plotSmoothers(traj2, sub_counts, gene = "Gm28376", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()

## FOR LINEAGE 3
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,3]
cellweights <- slingCurveWeights(embryo) [,3]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj3_Fcrla.pdf", width=8, height=4)
plotSmoothers(traj3, sub_counts, gene = "Fcrla", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()

## FOR LINEAGE 4
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,4]
cellweights <- slingCurveWeights(embryo) [,4]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj4_Hist2h2ac.pdf", width=8, height=4)
plotSmoothers(traj4, sub_counts, gene = "Hist2h2ac", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()



```

--> I did all the tutorial from [here](https://hectorrdb.github.io/condimentsPaper/articles/Fibrosis.html) 
----> This part Differential fate selection FAIL; but not sure it is usefull

--> prog_res? Lineage all pvalue signif! (paste in the ppt Conchi 20231005)

--> The icMAt is too long! I gave up running it. Author recommend using 6 knots (https://github.com/statOmics/tradeSeq/issues/121). Other discussion about [time](https://github.com/statOmics/tradeSeq/issues/41)

--> Same for the fitGAM, too long, so I run it into a Rscript as follow:
----> In the end; **what worked is to run trajectory per trajectory; it should last around 24-72hours per trajectory. Using the top 10k features do not reduce time significantly.**
----> Parrelization can be stuck sometime so I ll run it without...

--> The pseudotime diff looks weirk so let s use log2fc files anmd not the raw not log2fc filtered


```bash
conda activate condiments_V5
# test
sbatch scripts/fitGAM_6knots.sh # 5756977 fail
sbatch scripts/fitGAM_6knots_parralell_subset.sh # 6104946 fail

# trajectory per trajectory (all features, no parralelization)
sbatch scripts/fitGAM_6knots_traj1.sh # 6110340 ok
sbatch scripts/fitGAM_6knots_traj2.sh # 6110351 ok
sbatch scripts/fitGAM_6knots_traj3.sh # 6110359 ok
sbatch scripts/fitGAM_6knots_traj4.sh # 6110363 ok
sbatch scripts/fitGAM_6knots_traj5.sh # 6110370 ok
sbatch scripts/fitGAM_6knots_traj6.sh # 6110374 ok
sbatch scripts/fitGAM_6knots_traj7.sh # 6110376 ok
sbatch scripts/fitGAM_6knots_traj8.sh # 6110378 ok

```

--> `scripts/fitGAM_6knots.sh` 24hrs still running... FAIL time limit after 6 days LOL!!!! same for `fitGAM_6knots_parralell_subset.sh`

--> Without parralell processing trajectory per traj works great!! 24-72hrs to run

to check:
- https://www.bioconductor.org/packages/devel/bioc/vignettes/condiments/inst/doc/condiments.html
- For [plots](https://bioconductor.org/packages/devel/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html)






## Fuicntional analysis on pseudotime clustered genes

This analysis is done in conda activate deseq2 to use fgsea packages

--> Functional analysis following [this](https://hectorrdb.github.io/condimentsPaper/articles/TGFB.html); for BP (C5 BP) and pathway KEGG



```R
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Mm.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")



# import gene set files and convert as entrezid
gene_clusters_traj_cl = read.table(file = c("output/condiments/gene_clusters_traj3_l2fc2_cl7.txt"),   # CHAGNE TITLE  !!!!
                                 header=TRUE, sep="\t", stringsAsFactors=FALSE) %>% 
  as_tibble() %>%
  mutate(entrez_id = mapIds(org.Mm.eg.db, keys = gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first") ) %>%
  dplyr::filter(cluster == 7) %>%       # !!!!!!!!  CHAGNE TITLE   !!!!!!!!!!!!!!!!
  dplyr::select(entrez_id) %>%
  na.omit()

# transform into vector
gene_clusters_traj_cl_vector <- pull(gene_clusters_traj_cl, entrez_id)

# Traj1
## cluster per cluster
enrichKEGG <- enrichKEGG(gene   = gene_clusters_traj_cl$entrez_id ,
                         organism = "mmu",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
pdf("output/condiments/dotplot_KEGG_gene_clusters_traj3_cl7.pdf", width=7, height=7)  # !!!!!!!!  CHAGNE TITLE !!!!!!!!!!!!
dotplot(enrichKEGG, showCategory=20)
dev.off()


ego <- enrichGO(gene = gene_clusters_traj_cl_vector, 
                keyType = "ENTREZID",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Mm.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)           
pdf("output/condiments/dotplot_BP_gene_clusters_traj3_cl7.pdf", width=7, height=7) # !!!!!!!!  CHAGNE TITLE !!!!!!!!!!!!
dotplot(ego, showCategory=20)
dev.off()























## all cluster together THERE IS THE CODE, BUT THAT IS VERY CONFUSIGN SO FUCK IT
### Initialize a list to store the gene clusters for each iteration
gene_clusters_list <- list()

### Loop through clusters 1 to 10
for (i in 1:10) {
  gene_clusters_list[[paste0("gene_clusters_traj1_cl", i)]] <- read.table(file = "output/condiments/gene_clusters_traj1_l2fc2_cl10.txt",
                                                         header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
    as_tibble() %>%
    mutate(entrez_id = mapIds(org.Mm.eg.db, keys = gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>%
    dplyr::filter(cluster == i) %>%
    dplyr::select(entrez_id) %>%
    na.omit()
}

### Extract vectors for each cluster
vectors_list <- lapply(gene_clusters_list, function(x) pull(x, entrez_id))
names(vectors_list) <- paste0("gene_clusters_traj1_cl", 1:10, "_vector")

### Combine the vectors into a list
entrez_list <- list(cluster_3 = cluster_3_vector, 
                    cluster_4 = cluster_4_vector, 
                    cluster_6 = cluster_6_vector, 
                    cluster_7 = cluster_7_vector, 
                    cluster_8 = cluster_8_vector)


## GO
compGO <- compareCluster(geneCluster   = vectors_list,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Mm.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
pdf("output/condiments/dotplot_BP_gene_clusters_traj1.pdf", width=10, height=25)
dotplot(compGO, showCategory = 15, title = "GO_Biological Process Enrichment Analysis")
dev.off()

```


# New time point 20231206
Fodler label as `12.5.23 E7 cYAP KO and 24h human gastruloids treated with DASATINIB` in basespace

***NOTE: the 24hr UNTREATED are .gz.temp! Weird; I just remove the .temp...***

## Counting with cellranger count

```bash 
conda activate scRNAseq
which cellranger

# Run counting sample per sample for humangastruloid
sbatch scripts/cellranger_count_humangastruloid_UNTREATED24hr.sh # 9010720 FAIL; with QCNOTfail 9047403 ok
sbatch scripts/cellranger_count_humangastruloid_DASATINIB24hr.sh # 9010738 ok

# Run count for embryo using mice genome
sbatch scripts/cellranger_count_embryo_control_E7_mice.sh # 9010616 FAIL; with QCNOTfail 9016750 ok
sbatch scripts/cellranger_count_embryo_cYAPKO_E7_mice.sh # 9010675 ok
```
--> Weird bug for `embryo_control_E7`: `I0 error in FASTQ dile ... line: 0 failed to fill whole buffer`
----> There is 2 embryo_control_E7 files in BaseSpace; QC fail and not  QC failed; lets test both;
------> I indeed imported the QC fail one; let's use the QC not fail: *Success*


--> Bug for `cellranger_count_humangastruloid_UNTREATED24hr`: `Sequence and quality length mismatch: file: "/scr1/users/roulet/Akizu_Lab/002_scRNAseq/003__YAP1/input/24hgastruloidhumanUN/24hgastruloidhumanUN_S1_L001_R2_001.fastq.gz", line: 563052004`
----> There is 2 24hgastruloidhumanUN files in BaseSpace; QC fail and not  QC failed; lets test both;
------> I indeed imported the QC fail one; let's use the QC not fail: *Success*




## RNA contamination and doublet detection
- doublet detection using [scrublet](https://github.com/swolock/scrublet) **on the filtered matrix**
- ambient RNA correction using `soupX` in R before generating the Seurat object

```bash
srun --mem=500g --pty bash -l
conda deactivate # base environment needed
python3 scrublet.py [input_path] [output_path]
# Run doublet detection/scrublet sample per sample
python3 scripts/scrublet_doublets.py E7mousecontrolQCNOTfail/outs/filtered_feature_bc_matrix output/doublets/embryo_E7_control.tsv
python3 scripts/scrublet_doublets.py E7mousecYAPKO/outs/filtered_feature_bc_matrix output/doublets/embryo_E7_cYAPKO.tsv

python3 scripts/scrublet_doublets.py 24hgastruloidhumanUNQCNOTfail/outs/filtered_feature_bc_matrix output/doublets/humangastruloid_UNTREATED24hr.tsv
python3 scripts/scrublet_doublets.py 24hgastruloidhumanDASA/outs/filtered_feature_bc_matrix output/doublets/humangastruloid_DASATINIB24hr.tsv

```
Doublet detection score:
- humangastruloid_UNTREATED24hr: 42% (previous time: 0% doublet)
- humangastruloid_DASATINIB24hr: 0.2% (previous time: 34.2% doublet)
- embryo_E7_control: 7.2% (previous time: 0.1% doublet)
- embryo_E7_cYAPKO: 6.2% (previous time: 3.6%)
--> Successfully assigned doublet





# embryo E7 (second sample) analysis in Seurat


Then `conda activate scRNAseq` and go into R and filtered out **RNA contamination and start with SEURAT**.

Tuto seurat [here](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html)

**SoupX RNA cleaning**

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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

# soupX decontamination
## Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('E7mousecontrolQCNOTfail/outs') 
sc = load10X('E7mousecYAPKO/outs') 

## Assess % of conta
pdf("output/soupX/autoEstCont_E7mousecYAPKO.pdf", width=10, height=10)
pdf("output/soupX/autoEstCont_E7mousecontrol.pdf", width=10, height=10)
sc = autoEstCont(sc)
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/out_E7mousecontrol.RData")
save(out, file = "output/soupX/out_E7mousecYAPKO.RData")

```

```bash
conda activate scRNAseqV2
```

Let's 1st try to **perform clustering that look like the E775 time-point** (use same marker); display some genes and generate a **shiny app for Conchi** --> Shared to Conchi and discuss; then see if we put both exp together



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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

## Load the matrix and Create SEURAT object
load("output/soupX/out_E7mousecontrol.RData")
srat_WT_E7 <- CreateSeuratObject(counts = out, project = "WT_E7") # 32,285 features across 1,829 samples

load("output/soupX/out_E7mousecYAPKO.RData")
srat_cYAPKO_E7 <- CreateSeuratObject(counts = out, project = "cYAPKO_E7") # 32,285 features across 1,897 samples


# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
srat_WT_E7[["percent.mt"]] <- PercentageFeatureSet(srat_WT_E7, pattern = "^mt-")
srat_WT_E7[["percent.rb"]] <- PercentageFeatureSet(srat_WT_E7, pattern = "^Rp[sl]")

srat_cYAPKO_E7[["percent.mt"]] <- PercentageFeatureSet(srat_cYAPKO_E7, pattern = "^mt-")
srat_cYAPKO_E7[["percent.rb"]] <- PercentageFeatureSet(srat_cYAPKO_E7, pattern = "^Rp[sl]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/embryo_E7_control.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_WT_E7 <- AddMetaData(srat_WT_E7,doublets)
srat_WT_E7$Doublet_score <- as.numeric(srat_WT_E7$Doublet_score) # make score as numeric
head(srat_WT_E7[[]])

doublets <- read.table("output/doublets/embryo_E7_cYAPKO.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_cYAPKO_E7 <- AddMetaData(srat_cYAPKO_E7,doublets)
srat_cYAPKO_E7$Doublet_score <- as.numeric(srat_cYAPKO_E7$Doublet_score) # make score as numeric
head(srat_cYAPKO_E7[[]])



## Plot
pdf("output/seurat/VlnPlot_QC_embryo_control_E7.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QC_embryo_cYAPKO_E7.pdf", width=10, height=6)
VlnPlot(srat_cYAPKO_E7, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
pdf("output/seurat/FeatureScatter_QC_rb_mt_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_rb_mt_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "percent.rb", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_mt_doublet_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_mt_doublet_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "percent.mt", feature2 = "Doublet_score")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "nFeature_RNA", feature2 = "Doublet_score")
dev.off()




## After seeing the plot; add QC information in our seurat object
### V1 not super stringeant
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 500 & srat_WT_E7@meta.data$QC == 'Pass','Low_nFeature',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 500 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$percent.mt > 10 & srat_WT_E7@meta.data$QC == 'Pass','High_MT',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 500 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
table(srat_WT_E7[['QC']])
## 
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 500 & srat_cYAPKO_E7@meta.data$QC == 'Pass','Low_nFeature',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 500 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$percent.mt > 10 & srat_cYAPKO_E7@meta.data$QC == 'Pass','High_MT',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 500 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
table(srat_cYAPKO_E7[['QC']])


### V2 more stringeant
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC == 'Pass','Low_nFeature',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$percent.mt > 10 & srat_WT_E7@meta.data$QC == 'Pass','High_MT',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
table(srat_WT_E7[['QC']])
## 
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC == 'Pass','Low_nFeature',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$percent.mt > 10 & srat_cYAPKO_E7@meta.data$QC == 'Pass','High_MT',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
table(srat_cYAPKO_E7[['QC']])






# Quality check after QC filtering
pdf("output/seurat/VlnPlot_QCPass_embryo_control_E7.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QCPass_embryo_cYAPKO_E7.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QCPassV2_embryo_control_E7.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QCPassV2_embryo_cYAPKO_E7.pdf", width=10, height=6)

VlnPlot(subset(srat_cYAPKO_E7, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()







## subset my seurat object to only analyze the cells that pass the QC
srat_WT_E7 <- subset(srat_WT_E7, subset = QC == 'Pass')
srat_cYAPKO_E7 <- subset(srat_cYAPKO_E7, subset = QC == 'Pass')
srat_WT_E7$condition <- "WT_E7"
srat_cYAPKO_E7$condition <- "cYAPKO_E7"



mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
srat_WT_E7 <- NormalizeData(srat_WT_E7, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_WT_E7)
srat_WT_E7 <- ScaleData(srat_WT_E7, features = all.genes) # zero-centres and scales it

srat_cYAPKO_E7 <- NormalizeData(srat_cYAPKO_E7, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_cYAPKO_E7)
srat_cYAPKO_E7 <- ScaleData(srat_cYAPKO_E7, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
srat_WT_E7 <- CellCycleScoring(srat_WT_E7, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_WT_E7[[]]$Phase)
srat_cYAPKO_E7 <- CellCycleScoring(srat_cYAPKO_E7, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_cYAPKO_E7[[]]$Phase)

set.seed(42)

# elbow
## srat_WT_E7_elbow = SCTransform(srat_WT_E7, method = "glmGamPoi", ncells = 1624, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% RunPCA(npcs = 50, verbose = FALSE)


pdf("output/seurat/Elbow_srat_WT_E7QCV2.pdf", width=10, height=10)
ElbowPlot(srat_WT_E7_elbow) # 10 or 15 or 19
dev.off()


# clustering
## Optimal parameter 50 dim
srat_WT_E7 <- SCTransform(srat_WT_E7, method = "glmGamPoi", ncells = 788, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 50, verbose = FALSE)
srat_cYAPKO_E7 <- SCTransform(srat_cYAPKO_E7, method = "glmGamPoi", ncells = 862, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 50, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(srat_WT_E7 = srat_WT_E7, srat_cYAPKO_E7 = srat_cYAPKO_E7)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
embryoE7.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
embryoE7.combined.sct <- IntegrateData(anchorset = embryoE7.anchors, normalization.method = "SCT")
set.seed(42)
DefaultAssay(embryoE7.combined.sct) <- "integrated"
embryoE7.combined.sct <- RunPCA(embryoE7.combined.sct, verbose = FALSE, npcs = 50)
embryoE7.combined.sct <- RunUMAP(embryoE7.combined.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
embryoE7.combined.sct <- FindNeighbors(embryoE7.combined.sct, reduction = "pca", k.param = 5, dims = 1:50)
embryoE7.combined.sct <- FindClusters(embryoE7.combined.sct, resolution = 0.2, verbose = FALSE, algorithm = 4)
embryoE7.combined.sct$condition <- factor(embryoE7.combined.sct$condition, levels = c("WT_E7", "cYAPKO_E7")) # Reorder untreated 1st




## Optimal parameter 19 dim  _ V1
srat_WT_E7 <- SCTransform(srat_WT_E7, method = "glmGamPoi", ncells = 788, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO_E7 <- SCTransform(srat_cYAPKO_E7, method = "glmGamPoi", ncells = 862, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(srat_WT_E7 = srat_WT_E7, srat_cYAPKO_E7 = srat_cYAPKO_E7)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
embryoE7.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
embryoE7.combined.sct <- IntegrateData(anchorset = embryoE7.anchors, normalization.method = "SCT")
set.seed(42)
DefaultAssay(embryoE7.combined.sct) <- "integrated"
embryoE7.combined.sct <- RunPCA(embryoE7.combined.sct, verbose = FALSE, npcs = 19)
embryoE7.combined.sct <- RunUMAP(embryoE7.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryoE7.combined.sct <- FindNeighbors(embryoE7.combined.sct, reduction = "pca", k.param = 4, dims = 1:19)
embryoE7.combined.sct <- FindClusters(embryoE7.combined.sct, resolution = 0.07, verbose = FALSE, algorithm = 4)
embryoE7.combined.sct$condition <- factor(embryoE7.combined.sct$condition, levels = c("WT_E7", "cYAPKO_E7")) # Reorder untreated 1st


## Optimal parameter 19 dim  _ V2
srat_WT_E7 <- SCTransform(srat_WT_E7, method = "glmGamPoi", ncells = 788, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO_E7 <- SCTransform(srat_cYAPKO_E7, method = "glmGamPoi", ncells = 862, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(srat_WT_E7 = srat_WT_E7, srat_cYAPKO_E7 = srat_cYAPKO_E7)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
embryoE7.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
embryoE7.combined.sct <- IntegrateData(anchorset = embryoE7.anchors, normalization.method = "SCT")
set.seed(42)
DefaultAssay(embryoE7.combined.sct) <- "integrated"
embryoE7.combined.sct <- RunPCA(embryoE7.combined.sct, verbose = FALSE, npcs = 19)
embryoE7.combined.sct <- RunUMAP(embryoE7.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryoE7.combined.sct <- FindNeighbors(embryoE7.combined.sct, reduction = "pca", k.param = 20, dims = 1:19)
embryoE7.combined.sct <- FindClusters(embryoE7.combined.sct, resolution = 0.72, verbose = FALSE, algorithm = 4)
embryoE7.combined.sct$condition <- factor(embryoE7.combined.sct$condition, levels = c("WT_E7", "cYAPKO_E7")) # Reorder untreated 1st




## CLUSTERING testing
srat_WT_E7 <- SCTransform(srat_WT_E7, method = "glmGamPoi", ncells = 788, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO_E7 <- SCTransform(srat_cYAPKO_E7, method = "glmGamPoi", ncells = 862, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(srat_WT_E7 = srat_WT_E7, srat_cYAPKO_E7 = srat_cYAPKO_E7)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
embryoE7.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
embryoE7.combined.sct <- IntegrateData(anchorset = embryoE7.anchors, normalization.method = "SCT")
set.seed(42)
DefaultAssay(embryoE7.combined.sct) <- "integrated"
embryoE7.combined.sct <- RunPCA(embryoE7.combined.sct, verbose = FALSE, npcs = 19)
embryoE7.combined.sct <- RunUMAP(embryoE7.combined.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryoE7.combined.sct <- FindNeighbors(embryoE7.combined.sct, reduction = "pca", k.param = 20, dims = 1:19)
embryoE7.combined.sct <- FindClusters(embryoE7.combined.sct, resolution = 0.72, verbose = FALSE, algorithm = 4)
embryoE7.combined.sct$condition <- factor(embryoE7.combined.sct$condition, levels = c("WT_E7", "cYAPKO_E7")) # Reorder untreated 1st



####


pdf("output/seurat/UMAP_control_cYAPKO_E7_50dimOpt.pdf", width=10, height=4)
pdf("output/seurat/UMAP_control_cYAPKO_E7_19dimOpt.pdf", width=10, height=4)

pdf("output/seurat/UMAP_control_cYAPKO_E7_test.pdf", width=10, height=4)
pdf("output/seurat/UMAP_control_cYAPKO_E7_19dimOpt_V2.pdf", width=10, height=4)

DimPlot(embryoE7.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()




### Updated marker gene list from Alex + WINNER pro-winner
Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # 1 \  1 $ 1 + 1 %  1
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # 3 \ 4 $ 3  + 4 %  2
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # 2 \ 3 $ 8  + 3 %  7
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # 4 \ 2 $ 2 + 2 %  3
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # 5 \ 5 $ 5 +  5 %  4
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # 2 \ 3 $ 4  +  3 %  8
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # 7 \ 7 $ 6  +    %  5
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # 8 \ 8 $ 9  +    %  9
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # 11  \  $ 11 +   %  11
Caudal_Neurectoderm = c("Olig3","Hes3","Lrrc7","Cntfr") # 6  \   $ 5  +  %  4
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # 5  \   $ 7 + %  6
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # 9  \  $  10 + %  12
Intermediate_Mesoderm = c("Bik","Arg1","Meox1","Ism1") # 4 \  2  $  2 + %  mix of 3 5 6 (let's not name it as more general than 3 5 6)
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # 12  \  $  12 + %  13
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # 4 and 10 \  2  $  2 + %  10 (and 3)
Notocord = c("Shh","Noto","Vtn","Fgg") #  13 \   $  13 ( wierd separated) + %  14
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # 14  \   $  14 + %  15
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3") # between 12 and 13    + %  18


#### Conchi new list 20231214
Epiblast = c("Pou5f1", "Sox2", "Cdh1")
Primitive_Streak = c("Mixl1", "Eomes", "T")
Blood_Progenitor = c("Gata1", "Tal1", "Kdr", "Runx1", "Ets2", "Etv2")
Nascent_Mesoderm = c("Mesp1")
Endoderm =c("Sox17", "Foxa2")  # These cells will also be positive for Primitive streak markers, such as Eomes or T
ExE_Endoderm = c("Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6") # This cells will be negative or have low levels of Primitive streak markers
ExE_Ectoderm =c("Elf5", "Rhox5")


# test_newconchi list 20231214
## 50 dim
DefaultAssay(embryoE7.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_NODAL_50dim.pdf", width=10, height=12)
FeaturePlot(embryoE7.combined.sct, features = c("Nodal", "Axin2", "Sp5"), split.by = "condition", max.cutoff = 5, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Epiblast_List20231214_50dim.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Epiblast, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Primitive_Streak_List20231214_50dim.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Primitive_Streak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Blood_Progenitor_List20231214_50dim.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Blood_Progenitor, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Nascent_Mesoderm_List20231214_50dim.pdf", width=5, height=5)
FeaturePlot(embryoE7.combined.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Endoderm_List20231214_50dim.pdf", width=10, height=5)
FeaturePlot(embryoE7.combined.sct, features = Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_ExE_Endoderm_List20231214_50dim.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = ExE_Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_ExE_Ectoderm_List20231214_50dim.pdf", width=10, height=5)
FeaturePlot(embryoE7.combined.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()



## 19 dim
DefaultAssay(embryoE7.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_NODAL_19dim.pdf", width=10, height=12)
FeaturePlot(embryoE7.combined.sct, features = c("Nodal", "Axin2", "Sp5"), split.by = "condition", max.cutoff = 5, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Epiblast_List20231214_19dim.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Epiblast, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Primitive_Streak_List20231214_19dim.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Primitive_Streak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Blood_Progenitor_List20231214_19dim.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Blood_Progenitor, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Nascent_Mesoderm_List20231214_19dim.pdf", width=5, height=5)
FeaturePlot(embryoE7.combined.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Endoderm_List20231214_19dim.pdf", width=10, height=5)
FeaturePlot(embryoE7.combined.sct, features = Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_ExE_Endoderm_List20231214_19dim.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = ExE_Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_ExE_Ectoderm_List20231214_19dim.pdf", width=10, height=5)
FeaturePlot(embryoE7.combined.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()





# test
DefaultAssay(embryoE7.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA
# DefaultAssay(embryoE7.combined.sct) <- "RNA" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_NODAL_test.pdf", width=10, height=15)
FeaturePlot(embryoE7.combined.sct, features = c("Nodal", "Axin2", "Sp5"), split.by = "condition", max.cutoff = 5, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Notocord_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Notocord, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Mesenchyme_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Mesenchyme, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_ExE_Ectoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Nascent_Mesoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Blood_Progenitor_1_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Blood_Progenitor_1, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Epiblast_PrimStreak_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Epiblast_PrimStreak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Primordial_Germ_Cells_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Primordial_Germ_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Somitic_Mesoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Somitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Caudal_Mesoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Caudal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Haematodenothelial_progenitors_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Haematodenothelial_progenitors, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Paraxial_Mesoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Paraxial_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Caudal_Neurectoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Caudal_Neurectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Pharyngeal_Mesoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Pharyngeal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Blood_Progenitor_2_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Blood_Progenitor_2, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Intermediate_Mesoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Intermediate_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Surface_Ectoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Surface_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Mixed_Mesoderm_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Mixed_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Gut_test.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Gut, max.cutoff = 3, cols = c("grey", "red"))
dev.off()




# display features

pdf("output/seurat/UMAP_control_cYAPKO_E7_Phase.pdf", width=10, height=6)
DimPlot(embryoE7.combined.sct, reduction = "umap", split.by = "condition", label=TRUE, group.by = "Phase")
dev.off()
pdf("output/seurat/UMAP_control_cYAPKO_E7_nCount_RNA.pdf", width=10, height=6)
FeaturePlot(embryoE7.combined.sct, "nCount_RNA")
dev.off()
pdf("output/seurat/UMAP_control_cYAPKO_E7_percent.mt.pdf", width=10, height=6)
FeaturePlot(embryoE7.combined.sct, "percent.mt")
dev.off()
pdf("output/seurat/UMAP_control_cYAPKO_E7_nFeature_RNA.pdf", width=10, height=6)
FeaturePlot(embryoE7.combined.sct, "nFeature_RNA")
dev.off()

pdf("output/seurat/UMAP_control_cYAPKO_E7_Phase_QCV2.pdf", width=10, height=6)
DimPlot(embryoE7.combined.sct, reduction = "umap", split.by = "condition", label=TRUE, group.by = "Phase")
dev.off()
pdf("output/seurat/UMAP_control_cYAPKO_E7_nCount_RNA_QCV2.pdf", width=10, height=6)
FeaturePlot(embryoE7.combined.sct, "nCount_RNA")
dev.off()
pdf("output/seurat/UMAP_control_cYAPKO_E7_percent.mt_QCV2.pdf", width=10, height=6)
FeaturePlot(embryoE7.combined.sct, "percent.mt")
dev.off()
pdf("output/seurat/UMAP_control_cYAPKO_E7_nFeature_RNA_QCV2.pdf", width=10, height=6)
FeaturePlot(embryoE7.combined.sct, "nFeature_RNA")
dev.off()



DefaultAssay(embryoE7.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA
# QCV1
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Notocord.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Notocord, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Mesenchyme.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Mesenchyme, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_ExE_Ectoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Nascent_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Blood_Progenitor_1.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Blood_Progenitor_1, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Epiblast_PrimStreak.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Epiblast_PrimStreak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Primordial_Germ_Cells.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Primordial_Germ_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Somitic_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Somitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Caudal_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Caudal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Haematodenothelial_progenitors.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Haematodenothelial_progenitors, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Paraxial_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Paraxial_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Caudal_Neurectoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Caudal_Neurectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Pharyngeal_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Pharyngeal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Blood_Progenitor_2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Blood_Progenitor_2, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Intermediate_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Intermediate_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Surface_Ectoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Surface_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Mixed_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Mixed_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Gut.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Gut, max.cutoff = 3, cols = c("grey", "red"))
dev.off()


# QCV2
DefaultAssay(embryoE7.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Notocord_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Notocord, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Mesenchyme_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Mesenchyme, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_ExE_Ectoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Nascent_Mesoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Blood_Progenitor_1_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Blood_Progenitor_1, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Epiblast_PrimStreak_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Epiblast_PrimStreak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Primordial_Germ_Cells_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Primordial_Germ_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Somitic_Mesoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Somitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Caudal_Mesoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Caudal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Haematodenothelial_progenitors_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Haematodenothelial_progenitors, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Paraxial_Mesoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Paraxial_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Caudal_Neurectoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Caudal_Neurectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Pharyngeal_Mesoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Pharyngeal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Blood_Progenitor_2_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Blood_Progenitor_2, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Intermediate_Mesoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Intermediate_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Surface_Ectoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Surface_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Mixed_Mesoderm_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Mixed_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Gut_QCV2.pdf", width=10, height=10)
FeaturePlot(embryoE7.combined.sct, features = Gut, max.cutoff = 3, cols = c("grey", "red"))
dev.off()




### conclusion QCV1 and previous parameter clustering:
- Noto: bottom
- mesenchyme: top a tiny bit
- ExE_Ectoderm: good, cluster on the right
- Nascent_Mesdoerm: good, top
- Blood_Progenitor_1: good, top
- Epiblast_Primstreak: good, most of the cells middle
- Germ_cells: good, form small cluster  middle

- Gut: bad, lowly express/absent
- Haematodenothelial_progenitors: bad, lowly express
- Mesoderm: bad, lowly and dispersed
- Caudal_Mesoderm: bad, dispersed
- Somitic_Mesdoerm: bad, lowly express
- Surface_Ectoderm: bad, lowly express
- Paraxial_Mesdoerm: bad, dispersed
- Blood_Progenitor_2: bad, lowly express
- Caudal_Neuroectoderm: bad, lowly express
- Pharyngeal_Mesoderm: bad, everywhere
- Intermediate_Mesoderm: bad, lowly express


--> Any cluster in the non express regions? 


### conclusion V2 QCV2 and previous parameter clustering:
- Noto: bottom and up
- mesenchyme: bottom
- ExE_Ectoderm: good, left side
- Nascent_Mesdoerm: good, left side
- Blood_Progenitor_1: good, top
- Epiblast_Primstreak: good, most of the cells middle
- Germ_cells: good, form small cluster  middle
- Gut: good, right
- Caudal_Mesoderm: ok, up and top
- Surface_Ectoderm: good, middle
- Paraxial_Mesdoerm: ok, middle right
- Blood_Progenitor_2: good bottom
- Caudal_Neuroectoderm: good, top
- Pharyngeal_Mesoderm: good, bottom left
- Intermediate_Mesoderm: ok, bottom
- Mixed_Mesoderm: good, middle left

- Haematodenothelial_progenitors: bad, lowly express
- Mesoderm: bad, lowly and dispersed
- Somitic_Mesdoerm: bad, lowly dispersed





# Rename cluster
## 50 dim


new.cluster.ids <- c(
  "Epiblast",
  "Exe_Ectoderm",
  "Primitive_Streak-Nascent_Mesoderm_1",
  "Exe_Endoderm_1",
  "Primitive_Streak-Nascent_Mesoderm_2",
  "Endoderm",
  "Exe_Endoderm_2",
  "Blood_Progenitor"
)

names(new.cluster.ids) <- levels(embryoE7.combined.sct)
embryoE7.combined.sct <- RenameIdents(embryoE7.combined.sct, new.cluster.ids)

embryoE7.combined.sct$cluster.annot <- Idents(embryoE7.combined.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_control_cYAPKO_E7_label_50dim.pdf", width=12, height=6)
DimPlot(embryoE7.combined.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_control_cYAPKO_E7_label2_50dim.pdf", width=10, height=5)
DimPlot(embryoE7.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

#overlapping condition
pdf("output/seurat/UMAP_control_cYAPKO_E7_label_overlap_50dim.pdf", width=6, height=5)
DimPlot(embryoE7.combined.sct, reduction = "umap", group.by = "condition", pt.size = 0.5, cols = c("blue","red"))
dev.off()


# All in dotplot
DefaultAssay(embryoE7.combined.sct) <- "SCT"

  "Epiblast",
  "Exe_Ectoderm",
  "Primitive_Streak-Nascent_Mesoderm_1",
  "Exe_Endoderm_1",
  "Primitive_Streak-Nascent_Mesoderm_2",
  "Endoderm",
  "Exe_Endoderm_2",
  "Blood_Progenitor"


#### Conchi new list 20231214

Epiblast = c("Pou5f1", "Sox2", "Cdh1") # 
Exe_Ectoderm = c("Elf5", "Rhox5") # 
Primitive_Streak-Nascent_Mesoderm_1 = c("Mixl1", "Eomes", "T", "Mesp1") # 
Exe_Endoderm_1 = c("Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6") # 
Primitive_Streak-Nascent_Mesoderm_2 = c("Mixl1", "Eomes", "T", "Mesp1") # 
Endoderm = c("Sox17", "Foxa2") # 
Exe_Endoderm_2 = c("Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6") # 
Blood_Progenitor = c("Gata1", "Tal1", "Kdr", "Runx1", "Ets2", "Etv2") # 



all_markers <- c(
  "Pou5f1", "Sox2", "Cdh1",
  "Elf5", "Rhox5",
  "Mixl1", "Eomes", "T", "Mesp1",
  "Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6",
  "Sox17", "Foxa2",
  "Gata1", "Tal1", "Kdr", "Runx1", "Ets2", "Etv2"
)


levels(embryoE7.combined.sct) <- c(
  "Blood_Progenitor",
  "Endoderm",
  "Exe_Endoderm_2",
  "Exe_Endoderm_1",
  "Primitive_Streak-Nascent_Mesoderm_2",
  "Primitive_Streak-Nascent_Mesoderm_1",
  "Exe_Ectoderm",
  "Epiblast"
)


pdf("output/seurat/DotPlot_SCT_control_cYAPKO_E7_ident_50dim.pdf", width=12, height=3)
DotPlot(embryoE7.combined.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()



# SAve rds 
saveRDS(embryoE7.combined.sct, file = "output/seurat/embryoE7.combined.sct_50dim.rds")







## 19 dim 


new.cluster.ids <- c(
  "Epiblast",
  "Exe_Endoderm",
  "Primitive_Streak-Nascent_Mesoderm",
  "Exe_Ectoderm",
  "Endoderm",
  "Blood_Progenitor"
)


new.cluster.ids <- c(
  "Epiblast",
  "Exe_Ectoderm",
  "Primitive_Streak",
  "Exe_Endoderm_1",
  "Nascent_Mesoderm",
  "Cardiac_Mesoderm",
  "Endoderm",
  "Exe_Endoderm_2",
  "Blood_Progenitor"
)


names(new.cluster.ids) <- levels(embryoE7.combined.sct)
embryoE7.combined.sct <- RenameIdents(embryoE7.combined.sct, new.cluster.ids)

embryoE7.combined.sct$cluster.annot <- Idents(embryoE7.combined.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_control_cYAPKO_E7_label_19dim.pdf", width=12, height=5)
pdf("output/seurat/UMAP_control_cYAPKO_E7_label_19dim_V2.pdf", width=12, height=5)

DimPlot(embryoE7.combined.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_control_cYAPKO_E7_label2_19dim.pdf", width=10, height=5)
pdf("output/seurat/UMAP_control_cYAPKO_E7_label2_19dim_V2.pdf", width=10, height=5)
DimPlot(embryoE7.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.7, label.size = 4)
dev.off()



#overlapping condition
pdf("output/seurat/UMAP_control_cYAPKO_E7_label_overlap_19dim.pdf", width=6, height=4)

pdf("output/seurat/UMAP_control_cYAPKO_E7_label_overlap_19dim_V2.pdf", width=6, height=4)
DimPlot(embryoE7.combined.sct, reduction = "umap", group.by = "condition", pt.size = 0.5, cols = c("blue","red"))
dev.off()


# All in dotplot
DefaultAssay(embryoE7.combined.sct) <- "SCT"

#### Conchi new list 20231214
###### V1
Epiblast = c("Pou5f1", "Sox2", "Cdh1") # 
Exe_Endoderm = c("Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6") # 
Primitive_Streak-Nascent_Mesoderm = c("Mixl1", "Eomes", "T", "Mesp1") # 
Exe_Ectoderm = c("Elf5", "Rhox5") # 
Endoderm = c("Sox17", "Foxa2") # 
Blood_Progenitor = c("Gata1", "Tal1", "Kdr", "Runx1", "Ets2", "Etv2") # 

all_markers <- c(
  "Pou5f1", "Sox2", "Cdh1",
  "Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6",
  "Mixl1", "Eomes", "T", "Mesp1",
  "Elf5", "Rhox5",
  "Sox17", "Foxa2",
  "Gata1", "Tal1", "Kdr", "Runx1", "Ets2", "Etv2"
)
levels(embryoE7.combined.sct) <- c(
  "Blood_Progenitor",
  "Endoderm",
  "Exe_Ectoderm",
  "Primitive_Streak-Nascent_Mesoderm",
  "Exe_Endoderm",
  "Epiblast"
)

#### V2

Epiblast = c("Pou5f1", "Sox2", "Cdh1") # 
Primitive_Streak = c("Eomes", "Fgf8", "Nodal") # 
Nascent_Mesoderm = c("Mesp1") # 
Cardiac_Mesoderm = c("Hand1") # 
Endoderm = c("Sox17", "Foxa2") # 
Exe_Endoderm = c("Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6") # 
Exe_Ectoderm = c("Elf5", "Rhox5") # 
Blood_Progenitor = c("Gata1", "Tal1", "Kdr", "Runx1", "Ets2", "Etv2") # 


all_markers <- c(
  "Pou5f1", "Sox2", "Cdh1",
  "Eomes", "Fgf8", "Nodal",
  "Mesp1",
  "Hand1",
  "Sox17", "Foxa2",
  "Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6",
  "Elf5", "Rhox5",
  "Gata1", "Tal1", "Kdr", "Runx1", "Ets2", "Etv2"
)


levels(embryoE7.combined.sct) <- c(
  "Blood_Progenitor",
  "Exe_Ectoderm",
  "Exe_Endoderm_2",
  "Exe_Endoderm_1",
  "Endoderm",
  "Cardiac_Mesoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast"
)




pdf("output/seurat/DotPlot_SCT_control_cYAPKO_E7_ident_19dim.pdf", width=12, height=3)
pdf("output/seurat/DotPlot_SCT_control_cYAPKO_E7_ident_19dim_V2.pdf", width=12, height=3)

DotPlot(embryoE7.combined.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()




## previous E775 time point cell type marker genes expression

Blood_Progenitor = c("Sox18","Esam","Rhoj","Flt4","Gypa","Gata1","Cited4","Gfi1b") # 11  \  $ 11 +   %  11
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # 3 \ 4 $ 3  + 4 %  2
Epiblast_PrimStreak_Nascent_Mesoderm = c("Hesx1","Otx2","Pax2","Otx2os1","Col9a1","Mesp1","Osr1") # 1 \  1 $ 1 + 1 %  1


all_markers <- c(
  "Sox18","Esam","Rhoj","Flt4","Gypa","Gata1","Cited4","Gfi1b",
  "Tex19.1","Elf5","Wnt7b","Dnmt3l",
  "Hesx1","Otx2","Pax2","Otx2os1","Col9a1","Mesp1","Osr1"

)


levels(embryoE7.combined.sct) <- c(
  "Blood_Progenitor",
  "Endoderm",
  "Exe_Ectoderm",
  "Primitive_Streak-Nascent_Mesoderm",
  "Exe_Endoderm",
  "Epiblast"
)


pdf("output/seurat/DotPlot_SCT_control_cYAPKO_E7_ident_19dim_E775markers.pdf", width=12, height=3)
DotPlot(embryoE7.combined.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()





# SAve rds 
## saveRDS(embryoE7.combined.sct, file = "output/seurat/embryoE7.combined.sct_19dim_V2.rds")
embryoE7.combined.sct = readRDS(file = "output/seurat/embryoE7.combined.sct_19dim_V2.rds")





## Downsampling with bootstrap to compare the nb of cell per cell types

library("tidyverse")

### Identify the unique clusters
unique_clusters <- unique(Idents(embryoE7.combined.sct))

### Create empty matrices to store cell counts
control_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
cYAPKO_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(control_clusters_counts) <- unique_clusters
colnames(cYAPKO_clusters_counts) <- unique_clusters

### Loop through 100 iterations
embryoE7.combined.sct_WT <- which(embryoE7.combined.sct$orig.ident == 'WT_E7') # 788 cells
embryoE7.combined.sct_cYAPKO <- which(embryoE7.combined.sct$orig.ident == 'cYAPKO_E7') # 862 cells

for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  embryoE7.combined.sct_cYAPKO_downsample <- sample(embryoE7.combined.sct_cYAPKO, 788)
  embryoE7.combined.sct_integrated_downsample <- embryoE7.combined.sct[,c(embryoE7.combined.sct_WT, embryoE7.combined.sct_cYAPKO_downsample)]

  # Count nb of cells in each cluster
  control_clusters <- table(Idents(embryoE7.combined.sct_integrated_downsample)[embryoE7.combined.sct_integrated_downsample$condition == "WT_E7"])
  cYAPKO_clusters <- table(Idents(embryoE7.combined.sct_integrated_downsample)[embryoE7.combined.sct_integrated_downsample$condition == "cYAPKO_E7"])

  # Align the counts with the unique clusters
  control_clusters_counts[i, names(control_clusters)] <- as.numeric(control_clusters)
  cYAPKO_clusters_counts[i, names(cYAPKO_clusters)] <- as.numeric(cYAPKO_clusters)
}



### Calculate mean and standard error
mean_control_clusters <- colMeans(control_clusters_counts)
mean_cYAPKO_clusters <- colMeans(cYAPKO_clusters_counts)
std_error_cYAPKO_clusters <- apply(cYAPKO_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("WT_E7", "cYAPKO_E7")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- control_clusters_counts[j,i]
    contingency_table[1,2] <- cYAPKO_clusters_counts[j,i]
    contingency_table[2,1] <- sum(control_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(cYAPKO_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_control_clusters),
  untreated = mean_control_clusters,
  dasatinib = mean_cYAPKO_clusters,
  std_error_cYAPKO = std_error_cYAPKO_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_cYAPKO, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "WT_E7", "cYAPKO_E7"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.05, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("WT_E7", "cYAPKO_E7")) # Reorder untreated 1st
plot_data$cluster <- factor(plot_data$cluster, levels = c(  "Blood_Progenitor",
  "Exe_Ectoderm",
  "Exe_Endoderm_2",
  "Exe_Endoderm_1",
  "Endoderm",
  "Cardiac_Mesoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast")) # Reorder untreated 1st


# Plotting using ggplot2
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling10_clean_embryo_E7_19dim.pdf", width=5, height=4)
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling10_clean_embryo_E7_19dim_V2.pdf", width=5, height=4)

ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = filter(plot_data, condition == "cYAPKO_E7"),
    aes(label = significance, y = value + std_error_cYAPKO_clusters),
    vjust = -0.8,
    position = position_dodge(0.9), size = 5
  ) +
  scale_fill_manual(values = c("WT_E7" = "blue", "cYAPKO_E7" = "red")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()



# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

DefaultAssay(embryoE7.combined.sct) <- "RNA"

embryoE7.combined.sct <- NormalizeData(embryoE7.combined.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(embryoE7.combined.sct)
embryoE7.combined.sct <- ScaleData(embryoE7.combined.sct, features = all.genes) # zero-centres and scales it


## what genes change in different conditions for cells of the same type

embryoE7.combined.sct$celltype.stim <- paste(embryoE7.combined.sct$cluster.annot, embryoE7.combined.sct$condition,
    sep = "-")
Idents(embryoE7.combined.sct) <- "celltype.stim"

# use RNA corrected count for DEGs
embryoE7.combined.sct <- PrepSCTFindMarkers(embryoE7.combined.sct)



## Automation::
cell_types <- c("Blood_Progenitor",
  "Exe_Ectoderm",
  "Exe_Endoderm_2",
  "Exe_Endoderm_1",
  "Endoderm",
  "Cardiac_Mesoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast")

for (cell_type in cell_types) {
  response_name <- paste(cell_type, "cYAPKO_E7.response", sep = ".")
  ident_1 <- paste(cell_type, "-cYAPKO_E7", sep = "")
  ident_2 <- paste(cell_type, "-WT_E7", sep = "")

  response <- FindMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
  
  print(head(response, n = 15))
  
  file_name <- paste("output/seurat/", cell_type, "-cYAPKO_E7_response_19dim_V2.txt", sep = "")
  write.table(response, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE)
}



### Find all markers 
all_markers <- FindAllMarkers(embryoE7.combined.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_WT_cYAPKO_all_markers_E7_19di_V2_.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# Display the top 10 CONSERVED marker genes of each cluster
Idents(embryoE7.combined.sct) <- "cluster.annot"

## DEGs cluster versus all other
Blood_Progenitor.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = "Blood_Progenitor", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor")
Endoderm.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = "Endoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Endoderm")
Exe_Ectoderm.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = "Exe_Ectoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Exe_Ectoderm")
Exe_Endoderm_1.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = "Exe_Endoderm_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Exe_Endoderm_1")
Exe_Endoderm_2.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = "Exe_Endoderm_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Exe_Endoderm_2")
Epiblast.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = "Epiblast", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Epiblast")
Nascent_Mesoderm.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = 'Nascent_Mesoderm', grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Nascent_Mesoderm")
Primitive_Streak.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = 'Primitive_Streak', grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Primitive_Streak")
Cardiac_Mesoderm.conserved <- FindConservedMarkers(embryoE7.combined.sct, assay = "RNA", ident.1 = 'Cardiac_Mesoderm', grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Cardiac_Mesoderm")


## Combine all conserved markers into one data frame
all_conserved <- bind_rows(Blood_Progenitor.conserved, Endoderm.conserved, Exe_Ectoderm.conserved, Exe_Endoderm_1.conserved, Exe_Endoderm_2.conserved, Epiblast.conserved , Primitive_Streak.conserved, Nascent_Mesoderm.conserved, Cardiac_Mesoderm.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers_embryo_E7_19dim_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Blood_Progenitor",
  "Exe_Ectoderm",
  "Exe_Endoderm_2",
  "Exe_Endoderm_1",
  "Endoderm",
  "Cardiac_Mesoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Blood_Progenitor",
  "Exe_Ectoderm",
  "Exe_Endoderm_2",
  "Exe_Endoderm_1",
  "Endoderm",
  "Cardiac_Mesoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast")))

## Find the top 3 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Blood_Progenitor",
  "Endoderm",
  "Exe_Ectoderm",
  "Primitive_Streak-Nascent_Mesoderm",
  "Exe_Endoderm",
  "Epiblast"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 3) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Blood_Progenitor",
  "Endoderm",
  "Exe_Ectoderm",
  "Primitive_Streak-Nascent_Mesoderm",
  "Exe_Endoderm",
  "Epiblast")))



## Write the top 10 conserved markers for each cluster to a file
write.table(top10_conserved, file = "output/seurat/srat_top10_conserved_markers_embryo_E7_19dim_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Visualize the top 10/3 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(embryoE7.combined.sct) <- c("Blood_Progenitor",
  "Exe_Ectoderm",
  "Exe_Endoderm_2",
  "Exe_Endoderm_1",
  "Endoderm",
  "Cardiac_Mesoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast")

pdf("output/seurat/DotPlot_SCT_top10_conserved_embryo_E7_19dim_V2.pdf", width=21, height=4)
DotPlot(embryoE7.combined.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# save
## saveRDS(embryoE7.combined.sct, file = "output/seurat/embryoE7.combined.sct_19dim_V2.rds")
embryoE7.combined.sct <- readRDS(file = "output/seurat/embryoE7.combined.sct_19dim.rds")



## Re-calculate DEGs keeping ALL genes
embryoE7.combined.sct$celltype.stim <- paste(embryoE7.combined.sct$cluster.annot, embryoE7.combined.sct$condition,
    sep = "-")
Idents(embryoE7.combined.sct) <- "celltype.stim"

Blood_Progenitor <- FindMarkers(embryoE7.combined.sct, ident.1 = "Blood_Progenitor-cYAPKO_E7", ident.2 = "Blood_Progenitor-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Endoderm <- FindMarkers(embryoE7.combined.sct, ident.1 = "Endoderm-cYAPKO_E7", ident.2 = "Endoderm-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")     
Exe_Ectoderm <- FindMarkers(embryoE7.combined.sct, ident.1 = "Exe_Ectoderm-cYAPKO_E7", ident.2 = "Exe_Ectoderm-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")     
Primitive_Streak <- FindMarkers(embryoE7.combined.sct, ident.1 = "Primitive_Streak-cYAPKO_E7", ident.2 = "Primitive_Streak-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Nascent_Mesoderm <- FindMarkers(embryoE7.combined.sct, ident.1 = "Nascent_Mesoderm-cYAPKO_E7", ident.2 = "Nascent_Mesoderm-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
Cardiac_Mesoderm <- FindMarkers(embryoE7.combined.sct, ident.1 = "Cardiac_Mesoderm-cYAPKO_E7", ident.2 = "Cardiac_Mesoderm-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")   
Exe_Endoderm_1 <- FindMarkers(embryoE7.combined.sct, ident.1 = "Exe_Endoderm_1-cYAPKO_E7", ident.2 = "Exe_Endoderm_1-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")   
Exe_Endoderm_2 <- FindMarkers(embryoE7.combined.sct, ident.1 = "Exe_Endoderm_2-cYAPKO_E7", ident.2 = "Exe_Endoderm_2-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Epiblast <- FindMarkers(embryoE7.combined.sct, ident.1 = "Epiblast-cYAPKO_E7", ident.2 = "Epiblast-WT_E7",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")   



### save output
#### write.table(Cardiac_Mesoderm, file = "output/seurat/Cardiac_Mesoderm-cYAPKO_response_E7_19dim_allGenes_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)



# Check some genes
DefaultAssay(embryoE7.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

## Conchi geneList 20231214
pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_NODAL_19dim.pdf", width=10, height=12)
FeaturePlot(embryoE7.combined.sct, features = c("Yap1","Nodal", "Axin2", "Sp5"), split.by = "condition", max.cutoff = 5, cols = c("grey", "red"))
dev.off()


## YAP1-NODAL

YAP1_NODAL <- c("Shh", "Dmrt1", "Cited2", 'Dact2', "Tgif2", "Acvr1b", 'Smad3', 'Dand5', 'Smad2','Cer1', 'Nodal', 'Foxh1', 'Dact1', 'Tdgf1', 'Cripto3', 'Acvr1c', 'Cfc1', 'Cfc1b') 

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_19dim_YAP1_NODAL.pdf", width=10, height=70)
FeaturePlot(embryoE7.combined.sct, features = YAP1_NODAL, max.cutoff = 4, cols = c("grey", "red"), split.by = "condition")
dev.off()

## Anxa gene family


pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_ANXA_19dim.pdf", width=10, height=25)
FeaturePlot(embryoE7.combined.sct, features = c("Anxa1","Amotl2","Ccn1","Ccn2","Tagln","Anxa3","Chchd2"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()



## PRC2-related genes

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_PRC2_19dim.pdf", width=10, height=22)
FeaturePlot(embryoE7.combined.sct, features = c("Ezh2", "Kdm6b", "Gata6", "Cdx2","T"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


### Check YAP1 hippo bulk-regulated genes

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_YAP1_hippo_19dim.pdf", width=7, height=30)

FeaturePlot(embryoE7.combined.sct, features = c("Tead1", "Tead4", "Ccn2", "Ccn1","Areg","Myc","Gli2","Vim","Axl","Birc5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()



pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_YAP1_19dim.pdf", width=10, height=5)

FeaturePlot(embryoE7.combined.sct, features = c("Yap1"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()


## Wnt related

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_Wnt_19dim.pdf", width=10, height=7)

FeaturePlot(embryoE7.combined.sct, features = c("Axin2", "Sp5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()



## Ridge plot

embryoE7.combined.sct_Epiblast <- subset(embryoE7.combined.sct, idents = "Epiblast")

pdf("output/seurat/RidgePlot_SCT_control_cYAPKO_E7_NODAL_Epiblast.pdf", width=10, height=5)
RidgePlot(embryoE7.combined.sct_Epiblast, features = c('Nodal'), cols = c("blue","red"), group.by = 'condition', ncol =1)
dev.off()

pdf("output/seurat/RidgePlot_SCT_control_cYAPKO_E7_NODAL_Epiblast_V2.pdf", width=5, height=2)
RidgePlot(embryoE7.combined.sct_Epiblast, features = c('Nodal'), cols = c("blue","red"), group.by = 'condition', ncol =1)
dev.off()


```


## Pursue E7 only analysis, meeting Liz 20240419 - paper E7 embryo

Load data from where we stop, DEG has been done; RNA is scale and norm already!
- GSEA
- SCPA
- heatmap


```bash
srun --mem=500g --pty bash -l
conda activate scRNAseqV2
```

```R
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
library("gprofiler2") 

embryoE7.combined.sct <- readRDS(file = "output/seurat/embryoE7.combined.sct_19dim_V2.rds")

# GSEA plot (code2)
library("fgsea")
## Re-calculate DEGs keeping ALL genes -- ALREADY DONE UP
DefaultAssay(embryoE7.combined.sct) <- "RNA"



#### import all clsuter DEGs output :
cluster_types <- c("Blood_Progenitor", "Endoderm", "Exe_Ectoderm", "Primitive_Streak", "Nascent_Mesoderm", "Cardiac_Mesoderm", "Exe_Endoderm_1", "Exe_Endoderm_2", "Epiblast")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-cYAPKO_response_E7_19dim_allGenes_V2.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}


## load list of genes to test
Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

fgsea_sets <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt")$Genes,
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt")$Genes,
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt")$Genes,
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")$Genes
)





## Rank genes based on FC
genes <- Epiblast %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(ES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -NES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()




## plot GSEA
pdf("output/Pathway/GSEA_E7_ManualGeneLists_Vitamin_A-Epiblast.pdf", width=5, height=3)

plotEnrichment(fgsea_sets[["Vitamin_A"]],
               ranks) + labs(title="Vitamin_A-Epiblast") +
               theme_bw()
dev.off()



# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c("Blood_Progenitor", "Endoderm", "Exe_Ectoderm", "Primitive_Streak", "Nascent_Mesoderm", "Cardiac_Mesoderm", "Exe_Endoderm_1", "Exe_Endoderm_2", "Epiblast")

## Initialize an empty list to store the results for each cluster type
all_results <- list()
## Loop over each cluster type
for (cluster in cluster_types) {
  
  # Extract genes for the current cluster
  genes <- get(cluster) %>% 
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC)
  
  ranks <- deframe(genes)
  
  # Run GSEA for the current cluster
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(ES))
  
  # Extract summary table and add cluster column
  fgseaResTidy_summary = fgseaResTidy %>% 
    dplyr::select(pathway, pval, padj, ES, size, NES) %>%
    mutate(cluster = cluster) %>%
    arrange(padj) %>% 
    head()
  
  # Store results in the list
  all_results[[cluster]] <- fgseaResTidy_summary
}
## Combine results from all cluster types into one table
final_results <- bind_rows(all_results, .id = "cluster")


write.table(final_results, file = c("output/Pathway/gsea_output_E7_ManualGeneLists.txt"), sep = "\t", quote = FALSE, row.names = FALSE)



# Heatmap all GSEA


pdf("output/Pathway/heatmap_E7_gsea_padj.pdf", width=5, height=5)
ggplot(final_results, aes(x=cluster, y=pathway, fill=ES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Enrichment\nScore") +
  geom_text(aes(label=sprintf("%.2f", ES)), 
            color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()



pdf("output/Pathway/heatmap_E7_gsea_pval_greyTile.pdf", width=5, height=5)
ggplot(final_results, aes(x=cluster, y=pathway)) + 
  geom_tile(aes(fill = ifelse(pval <= 0.05, NES, NA)), color = "black") +  # Conditional fill based on significance
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Normalized Enrichment\nScore", na.value="grey") +
  # geom_text(aes(label=sprintf("%.2f", ES)), 
  #           color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change between pvalue, qvalue,p.adjust
  #           size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1

dev.off()



# SCPA (code3)


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

DefaultAssay(embryoE7.combined.sct) <- "RNA" # Recommended 

# Import manually curated gene list as a list
#### Need to folllow this format
pathways <- msigdbr("Mus musculus", "C2") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise
pathways$REACTOME_SIGNALING_BY_RETINOIC_ACID$Genes
### Convert table to list of genes
Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

pathways <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt"),
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt"),
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt"),
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")
)

pathways$Hippo$Genes


# Compare pathway activity between condition within  
cell_types <- unique(embryoE7.combined.sct$cluster.annot)
embryoE7.combined.sct_split <- SplitObject(embryoE7.combined.sct, split.by = "condition")

# SCPA comparison 
## keep all columns
scpa_out_all <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryoE7.combined.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryoE7.combined.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_all[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}
## keep only PAthway and qval column (Best for heatmap qval representation)
scpa_out <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryoE7.combined.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryoE7.combined.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}


saveRDS(scpa_out_all, file = "output/Pathway/scpa_out_eachCellTypes_ManualGeneLists_E7.rds")

console_output <- capture.output(print(scpa_out_all))
writeLines(console_output, "output/Pathway/scpa_out_console_ManualGeneLists_E7.txt")

## Save output as R object
# saveRDS(scpa_out, file = "output/Pathway/scpa_out_eachCellTypes_E7.rds")
# scpa_out_all <- readRDS("output/Pathway/scpa_out_all_eachCellTypes_E7.rds")


# plot output

# Code to save output for each cell type comparison
clusters = c(
"Blood_Progenitor", "Endoderm", "Exe_Ectoderm", "Primitive_Streak", "Nascent_Mesoderm", "Cardiac_Mesoderm", "Exe_Endoderm_1", "Exe_Endoderm_2", "Epiblast"
) # 

### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(embryoE7.combined.sct,
                       meta1 = "condition", value_meta1 = "WT_E7",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  cYAPKO <- seurat_extract(embryoE7.combined.sct,
                           meta1 = "condition", value_meta1 = "cYAPKO_E7",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_E7_ManualGeneLists_", cluster, ".txt")
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}



# BALOW CLEAN CODE; dotplot>enrichplot>heatmap>violinplot
# DOT PLOT pepresetnation
## load all the comparison for each cell type (FC qval information)
clusters = c(
"Blood_Progenitor", "Endoderm", "Exe_Ectoderm", "Primitive_Streak", "Nascent_Mesoderm", "Cardiac_Mesoderm", "Exe_Endoderm_1", "Exe_Endoderm_2", "Epiblast"
) # 

## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_E7_ManualGeneLists_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster))
# --> import the gene pathways (used table from Conchi `Pathwyas of interest*.xlsx`)
pathways = c("Hippo", "Nodal_TGF", "WNT", "Vitamin_A")
all_data_pathways = as_tibble(all_data) %>% 
  filter(Pathway %in% pathways)

# tidy the df
all_data_pathways_tidy <- all_data_pathways %>%
  filter(qval >= 1.4) 



# REFINE COLOR
custom_color <- function(fc_value){
  ifelse(fc_value >= 5, "strong_blue", 
         ifelse(fc_value > 2, "light_blue",
                ifelse(fc_value <= -5, "strong_red", 
                       ifelse(fc_value < -2, "light_red", "grey"))))
}

# Add a column for this custom color
all_data_pathways_tidy <- all_data_pathways_tidy %>%
  mutate(custom_col = sapply(FC, custom_color))

pdf("output/Pathway/dotplot_E7_ManualGeneLists_FCtresh.pdf", width=10, height=4)
ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, color = custom_col), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = c("strong_blue" = "dodgerblue4",
                                "light_blue" = "lightblue2",
                                "grey" = "grey",
                                "light_red" = "indianred1",
                                "strong_red" = "red3")) +
  theme_bw() +
  labs(size = "q-value", color = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



pdf("output/Pathway/dotplot_E7_ManualGeneLists_grey.pdf", width=10, height=4)
ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, fill = "black"), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  labs(size = "q-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# heatmap XXX to do if needed XXX (code1)

XXXX Create new lsit of genes XXX

Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

pathways <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt"),
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt"),
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt"),
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")
)

pathways$Hippo$Genes



genes_of_interest <- pathways$Hippo$Genes
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(embryo.combined.sct@assays$RNA@data)]
### extract WT and cYAPKO gene expression values from RNA assay
WT_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "WT" & embryo.combined.sct$cluster.annot == "Unknown_1"]]
cYAPKO_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO" & embryo.combined.sct$cluster.annot == "Unknown_1"]]
### mean expression values for each gene
WT_mean <- rowMeans(WT_expression)
cYAPKO_mean <- rowMeans(cYAPKO_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  WT = WT_mean,
  cYAPKO = cYAPKO_mean
) %>% 
pivot_longer(cols = c(WT, cYAPKO), names_to = "Condition", values_to = "Expression")
## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(WT_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)


# if few genes:
pdf("output/Pathway/heatmap_embryo_RNAexpression_E7_ManualGeneLists_Nodal_TGF_Blood_Progenitor_1.pdf", width=6, height=2)


library("viridis")
ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()

# if many genesL
pdf("output/Pathway/heatmap_embryo_RNAexpression_E7_ManualGeneLists_WNT_Blood_Progenitor_1.pdf", width=3, height=2)

ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") 
dev.off()



# SCPA and GSEA 
gsea_output = final_results %>% 
  rename("pathway" = "Pathway") %>%
  dplyr::select(Pathway, pval, padj, ES, cluster)
scpa_output = all_data_pathways_tidy  %>%
  dplyr::select(Pathway, qval, cluster)

scpa_gsea = gsea_output %>%
  left_join(scpa_output)

## Filter data for GSEA padj < 0.05
filtered_data <- scpa_gsea
## Assign color based on your criteria
filtered_data$is_significant <- filtered_data$padj < 0.05

pdf("output/Pathway/dotplot_scpa_gsea_E7.pdf", width=10, height=4)

ggplot(filtered_data, aes(x=cluster, y=Pathway, size=qval)) +
  geom_point(aes(color=ifelse(is_significant, ES, NA), 
                 shape=is_significant), alpha=0.7) +
  scale_color_gradient2(low="blue", high="red", midpoint=0, na.value="black", name = "GSEA Enrichment score") +
  scale_size_continuous(name="SCPA qval") +  # Title for size scale
  scale_shape_manual(values=c(`TRUE`=16, `FALSE`=16)) + # Using circle shape for both
  theme_light() +
  labs(title="SCPA and GSEA", x="Cluster", y="Pathway") +
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# show some UMAP between conditions
DefaultAssay(embryoE7.combined.sct) <- "SCT"
## Primitive Streak (code4)

pdf("output/seurat/FeaturePlot_SCT_control_cYAPKO_E7_T_19dim_V2.pdf", width=10, height=5)

FeaturePlot(embryoE7.combined.sct, features = c("T"), split.by = "condition", max.cutoff = 5, cols = c("grey", "red"))
dev.off()



# Test different Pathway collections and generate enrichment plot for each cell types (C8, C5, H)
DefaultAssay(embryoE7.combined.sct) <- "SCT" # SCT is OK https://github.com/jackbibby1/SCPA/issues/67

## import Pathways
pathways <- msigdbr("Mus musculus", "C5") %>%  # C2
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
"Blood_Progenitor", "Endoderm", "Exe_Ectoderm", "Primitive_Streak", "Nascent_Mesoderm", "Cardiac_Mesoderm", "Exe_Endoderm_1", "Exe_Endoderm_2", "Epiblast"
)
clusters = c(
"Epiblast"
)

### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(embryoE7.combined.sct,
                       meta1 = "condition", value_meta1 = "WT_E7",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  cYAPKO <- seurat_extract(embryoE7.combined.sct,
                           meta1 = "condition", value_meta1 = "cYAPKO_E7",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(cYAPKO, WT),   # list(WT, cYAPKO) = wrong direction for fc
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCT_E7_19dim_V2_SCPA_C5_", cluster, ".txt")       # CHANGE HERE PATHWAYS !!!!!!
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}
### --> WATCH OUT HERE list(WT, cYAPKO) lead to FC in wrong direction!!! for C2 file are shit (used instead the one that end with -cYPAKO_WT), for C5 that is good



### Plot to count how many pos/neg FC
SCT_SCPA_C2 <- list()
for (cluster in clusters) {
  # Create the filename based on the cluster name
  input_filename <- paste0("output/Pathway/SCT_E7_19dim_V2_SCPA_C2_", cluster, ".txt")
  
  # Import the data
  current_data <- read.delim(input_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Store the data in the list with the cluster name as the list name
  SCT_SCPA_C2[[cluster]] <- current_data
}



### Tidy the data; convert list into tibble

tidy_SCT_C2 <- map2_df(SCT_SCPA_C2, names(SCT_SCPA_C2), ~ {
    data_frame <- .x
    data_frame$cluster <- .y
    return(data_frame)
  }) %>%
  add_column(db = "C2") %>% 
  as_tibble() %>%
  mutate(FC_direction = case_when(
    FC < -0 ~ "positive",
    FC > 0 ~ "negative",
    TRUE ~ "neutral"  # for values between -1 and 1
  ))

# Plot
pdf("output/Pathway/Count_Pathway_FC_direction_SCT_E7_19dim_V2_SCPA_C2.pdf", width=6, height=3)
tidy_SCT_C2 %>%
  filter(qval > 1.4) %>%
  group_by(cluster, db, FC_direction) %>%
  summarise(n=n()) %>%
  ggplot(aes(x = FC_direction, y = n, fill = cluster)) +  # Moved fill inside aes
  geom_bar(stat = "identity", position = "dodge") +  # Add dodging
  facet_wrap(~db, scale = "free") +  # Facet by db
  labs(title = "Distribution of FC_direction per cluster",
       y = "Count", x = "FC Direction") +
  theme_bw()
dev.off()






# inverse pop1 and pop2 to have FC in good direction:
"Blood_Progenitor", "Endoderm", "Exe_Ectoderm", "Primitive_Streak", "Nascent_Mesoderm", "Cardiac_Mesoderm", "Exe_Endoderm_1", "Exe_Endoderm_2", "Epiblast"
#### Extract data for WT and cYAPKO based on current value
WT <- seurat_extract(embryoE7.combined.sct,
                      meta1 = "condition", value_meta1 = "WT_E7",
                      meta2 = "cluster.annot", value_meta2 = "Primitive_Streak")

cYAPKO <- seurat_extract(embryoE7.combined.sct,
                          meta1 = "condition", value_meta1 = "cYAPKO_E7",
                          meta2 = "cluster.annot", value_meta2 = "Primitive_Streak")

##### Compare pathways
cYAPKO_WT <- compare_pathways(samples = list(cYAPKO, WT),    # initially WT, cYPAKO
                              pathways = pathways,
                              parallel = TRUE, cores = 8)

write.table(cYAPKO_WT, file = "output/Pathway/SCPA_Primitive_Streak_SCT_E7_19dim_V2_SCPA_C2-cYAPKO_WT", sep = "\t", quote = FALSE, row.names = FALSE)



cYAPKO_WT_filter <- cYAPKO_WT %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

pdf("output/Pathway/plot_embryo_msigdb_Primitive_Streak_SCT_E7_19dim_V2_SCPA_C2-cYAPKO_WT.pdf", width=5, height=5)
ggplot(cYAPKO_WT_filter, aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = cYAPKO_WT_filter$color, stroke = 0.3) +
  xlim(-10, 10) +
  ylim(0, 5) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
dev.off()




# Rank and color FC with treshold 5
### import all C2 comparisons

SCPA_C2 <- list()
for (cluster in clusters) {
  # Create the filename based on the cluster name
  input_filename <- paste0("output/Pathway/SCPA_", cluster, "_SCT_E7_19dim_V2_SCPA_C2-cYAPKO_WT")
  
  # Import the data
  current_data <- read.delim(input_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Store the data in the list with the cluster name as the list name
  SCPA_C2[[cluster]] <- current_data
}

### convert list into tibble
tidy_SCT_C2 <- map2_df(SCPA_C2, names(SCPA_C2), ~ {
    data_frame <- .x
    data_frame$cluster <- .y
    return(data_frame)
  }) %>%
  as_tibble() %>%
  mutate(FC_direction = case_when(
    FC < -2.5 ~ "positive",
    FC > 2.5 ~ "negative",
    TRUE ~ "neutral"  # for values between -1 and 1
  ))

tidy_SCT_C2_Epiblast = 
  tidy_SCT_C2 %>%
  filter(cluster == "Epiblast", qval > 1.4) %>%
  arrange(qval) %>%
  mutate(rank = (row_number() - 1) / (n() - 1) * 100)

### plot rank


#### color scale
pdf("output/Pathway/plot_rank-Epiblast_E7_19dim_V2_SCPA_C2.pdf", width=5, height=5)
tidy_SCT_C2_Epiblast %>%
  ggplot(., aes(x = qval, y = rank)) +
  geom_jitter(aes(color = -FC, alpha = FC_direction), size=0.1, width = 0.2, height = 2) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, guide = "colourbar") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = 1, "neutral" = 0.5)) +
  labs(x = "Qval", y = "Pathway rank", size = "q-value", color = "Fold Change") +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept =95)
dev.off()







```

--> SCPA FC negative = less express in KO!; use file that finish with `-cYAPKO_WT`



# human gastruloid 24hr (second sample) analysis in Seurat


Then `conda activate scRNAseq` and go into R and filtered out **RNA contamination and start with SEURAT**.

Tuto seurat [here](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html)

**SoupX RNA cleaning**

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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

# soupX decontamination
## Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('24hgastruloidhumanUNQCNOTfail/outs') 
sc = load10X('24hgastruloidhumanDASA/outs') 

## Assess % of conta
pdf("output/soupX/autoEstCont_24hgastruloidhumanUN.pdf", width=10, height=10)
pdf("output/soupX/autoEstCont_24hgastruloidhumanDASA.pdf", width=10, height=10)
sc = autoEstCont(sc)
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/out_24hgastruloidhumanUN.RData")
save(out, file = "output/soupX/out_24hgastruloidhumanDASA.RData")

```

```bash
conda activate scRNAseqV2
```

Let's 1st try to **perform clustering that look like the E775 time-point** (use same marker); display some genes and generate a **shiny app for Conchi** --> Shared to Conchi and discuss; then see if we put both exp together



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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

## Load the matrix and Create SEURAT object
load("output/soupX/out_24hgastruloidhumanUN.RData")
srat_UNTREATED24hr <- CreateSeuratObject(counts = out, project = "UNTREATED24hr") # 32,285 features across 1,829 samples

load("output/soupX/out_24hgastruloidhumanDASA.RData")
srat_DASATINIB24hr <- CreateSeuratObject(counts = out, project = "DASATINIB24hr") # 32,285 features across 1,897 samples


# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
srat_UNTREATED24hr[["percent.mt"]] <- PercentageFeatureSet(srat_UNTREATED24hr, pattern = "^MT-")
srat_UNTREATED24hr[["percent.rb"]] <- PercentageFeatureSet(srat_UNTREATED24hr, pattern = "^RP[SL]")

srat_DASATINIB24hr[["percent.mt"]] <- PercentageFeatureSet(srat_DASATINIB24hr, pattern = "^MT-")
srat_DASATINIB24hr[["percent.rb"]] <- PercentageFeatureSet(srat_DASATINIB24hr, pattern = "^RP[SL]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/humangastruloid_UNTREATED24hr.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_UNTREATED24hr <- AddMetaData(srat_UNTREATED24hr,doublets)
srat_UNTREATED24hr$Doublet_score <- as.numeric(srat_UNTREATED24hr$Doublet_score) # make score as numeric
head(srat_UNTREATED24hr[[]])

doublets <- read.table("output/doublets/humangastruloid_DASATINIB24hr.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_DASATINIB24hr <- AddMetaData(srat_DASATINIB24hr,doublets)
srat_DASATINIB24hr$Doublet_score <- as.numeric(srat_DASATINIB24hr$Doublet_score) # make score as numeric
head(srat_DASATINIB24hr[[]])



## Plot
pdf("output/seurat/VlnPlot_QC_humangastruloid_UNTREATED24hr.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QC_humangastruloid_DASATINIB24hr.pdf", width=10, height=6)
VlnPlot(srat_DASATINIB24hr, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()



## After seeing the plot; add QC information in our seurat object
### V1 stringent
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 2500 & srat_UNTREATED24hr@meta.data$QC == 'Pass','Low_nFeature',srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 2500 & srat_UNTREATED24hr@meta.data$QC != 'Pass' & srat_UNTREATED24hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_UNTREATED24hr@meta.data$QC,sep = ','),srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$percent.mt > 15 & srat_UNTREATED24hr@meta.data$QC == 'Pass','High_MT',srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 2500 & srat_UNTREATED24hr@meta.data$QC != 'Pass' & srat_UNTREATED24hr@meta.data$QC != 'High_MT',paste('High_MT',srat_UNTREATED24hr@meta.data$QC,sep = ','),srat_UNTREATED24hr@meta.data$QC)
table(srat_UNTREATED24hr[['QC']])
## 
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 2500 & srat_DASATINIB24hr@meta.data$QC == 'Pass','Low_nFeature',srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 2500 & srat_DASATINIB24hr@meta.data$QC != 'Pass' & srat_DASATINIB24hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_DASATINIB24hr@meta.data$QC,sep = ','),srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$percent.mt > 15 & srat_DASATINIB24hr@meta.data$QC == 'Pass','High_MT',srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 2500 & srat_DASATINIB24hr@meta.data$QC != 'Pass' & srat_DASATINIB24hr@meta.data$QC != 'High_MT',paste('High_MT',srat_DASATINIB24hr@meta.data$QC,sep = ','),srat_DASATINIB24hr@meta.data$QC)
table(srat_DASATINIB24hr[['QC']])



### V2 more stringent
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 2750 & srat_UNTREATED24hr@meta.data$QC == 'Pass','Low_nFeature',srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 2750 & srat_UNTREATED24hr@meta.data$QC != 'Pass' & srat_UNTREATED24hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_UNTREATED24hr@meta.data$QC,sep = ','),srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$percent.mt > 15 & srat_UNTREATED24hr@meta.data$QC == 'Pass','High_MT',srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 2750 & srat_UNTREATED24hr@meta.data$QC != 'Pass' & srat_UNTREATED24hr@meta.data$QC != 'High_MT',paste('High_MT',srat_UNTREATED24hr@meta.data$QC,sep = ','),srat_UNTREATED24hr@meta.data$QC)
table(srat_UNTREATED24hr[['QC']])
## 
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 2750 & srat_DASATINIB24hr@meta.data$QC == 'Pass','Low_nFeature',srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 2750 & srat_DASATINIB24hr@meta.data$QC != 'Pass' & srat_DASATINIB24hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_DASATINIB24hr@meta.data$QC,sep = ','),srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$percent.mt > 15 & srat_DASATINIB24hr@meta.data$QC == 'Pass','High_MT',srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 2750 & srat_DASATINIB24hr@meta.data$QC != 'Pass' & srat_DASATINIB24hr@meta.data$QC != 'High_MT',paste('High_MT',srat_DASATINIB24hr@meta.data$QC,sep = ','),srat_DASATINIB24hr@meta.data$QC)
table(srat_DASATINIB24hr[['QC']])



### V2 more stringent
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 3000 & srat_UNTREATED24hr@meta.data$QC == 'Pass','Low_nFeature',srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 3000 & srat_UNTREATED24hr@meta.data$QC != 'Pass' & srat_UNTREATED24hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_UNTREATED24hr@meta.data$QC,sep = ','),srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$percent.mt > 15 & srat_UNTREATED24hr@meta.data$QC == 'Pass','High_MT',srat_UNTREATED24hr@meta.data$QC)
srat_UNTREATED24hr[['QC']] <- ifelse(srat_UNTREATED24hr@meta.data$nFeature_RNA < 3000 & srat_UNTREATED24hr@meta.data$QC != 'Pass' & srat_UNTREATED24hr@meta.data$QC != 'High_MT',paste('High_MT',srat_UNTREATED24hr@meta.data$QC,sep = ','),srat_UNTREATED24hr@meta.data$QC)
table(srat_UNTREATED24hr[['QC']])
## 
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 3000 & srat_DASATINIB24hr@meta.data$QC == 'Pass','Low_nFeature',srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 3000 & srat_DASATINIB24hr@meta.data$QC != 'Pass' & srat_DASATINIB24hr@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_DASATINIB24hr@meta.data$QC,sep = ','),srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$percent.mt > 15 & srat_DASATINIB24hr@meta.data$QC == 'Pass','High_MT',srat_DASATINIB24hr@meta.data$QC)
srat_DASATINIB24hr[['QC']] <- ifelse(srat_DASATINIB24hr@meta.data$nFeature_RNA < 3000 & srat_DASATINIB24hr@meta.data$QC != 'Pass' & srat_DASATINIB24hr@meta.data$QC != 'High_MT',paste('High_MT',srat_DASATINIB24hr@meta.data$QC,sep = ','),srat_DASATINIB24hr@meta.data$QC)
table(srat_DASATINIB24hr[['QC']])





# Quality check after QC filtering
pdf("output/seurat/VlnPlot_QCPassV1_humangastruloid_UNTREATED24hr.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QCPassV1_humangastruloid_DASATINIB24hr.pdf", width=10, height=6)


VlnPlot(subset(srat_DASATINIB24hr, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()







## subset my seurat object to only analyze the cells that pass the QC
srat_UNTREATED24hr <- subset(srat_UNTREATED24hr, subset = QC == 'Pass')
srat_DASATINIB24hr <- subset(srat_DASATINIB24hr, subset = QC == 'Pass')
srat_UNTREATED24hr$condition <- "UNTREATED24hr"
srat_DASATINIB24hr$condition <- "DASATINIB24hr"


s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
srat_UNTREATED24hr <- NormalizeData(srat_UNTREATED24hr, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_UNTREATED24hr)
srat_UNTREATED24hr <- ScaleData(srat_UNTREATED24hr, features = all.genes) # zero-centres and scales it

srat_DASATINIB24hr <- NormalizeData(srat_DASATINIB24hr, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_DASATINIB24hr)
srat_DASATINIB24hr <- ScaleData(srat_DASATINIB24hr, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
srat_UNTREATED24hr <- CellCycleScoring(srat_UNTREATED24hr, s.features = s.genes, g2m.features = g2m.genes)
table(srat_UNTREATED24hr[[]]$Phase)
srat_DASATINIB24hr <- CellCycleScoring(srat_DASATINIB24hr, s.features = s.genes, g2m.features = g2m.genes)
table(srat_DASATINIB24hr[[]]$Phase)

set.seed(42)

# elbow
## srat_UNTREATED24hr_elbow = SCTransform(srat_UNTREATED24hr, method = "glmGamPoi", ncells = 7420, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% RunPCA(npcs = 50, verbose = FALSE)


pdf("output/seurat/Elbow_srat_UNTREATED24hr.pdf", width=10, height=10)
ElbowPlot(srat_UNTREATED24hr_elbow) # 10 or 15 or 19
dev.off()


# clustering
set.seed(42)

### Previous optimal parameters

# Run SCTransform
## Version OK with 2000 treshold RNA
srat_UNTREATED24hr <- SCTransform(srat_UNTREATED24hr, method = "glmGamPoi", ncells = 7331, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 25, verbose = FALSE)
srat_DASATINIB24hr <- SCTransform(srat_DASATINIB24hr, method = "glmGamPoi", ncells = 6613, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 25, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(srat_UNTREATED24hr = srat_UNTREATED24hr, srat_DASATINIB24hr = srat_DASATINIB24hr)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
humangastruloid24hr.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
humangastruloid24hr.combined.sct <- IntegrateData(anchorset = humangastruloid24hr.anchors, normalization.method = "SCT")
# Perform integrated analysis (check active assay is 'integrated')

DefaultAssay(humangastruloid24hr.combined.sct) <- "integrated"
humangastruloid24hr.combined.sct <- RunPCA(humangastruloid24hr.combined.sct, verbose = FALSE, npcs = 25)
humangastruloid24hr.combined.sct <- RunUMAP(humangastruloid24hr.combined.sct, reduction = "pca", dims = 1:25, verbose = FALSE)
humangastruloid24hr.combined.sct <- FindNeighbors(humangastruloid24hr.combined.sct, reduction = "pca", k.param = 5, dims = 1:25)
humangastruloid24hr.combined.sct <- FindClusters(humangastruloid24hr.combined.sct, resolution = 0.1, verbose = FALSE, algorithm = 4)
humangastruloid24hr.combined.sct$condition <- factor(humangastruloid24hr.combined.sct$condition, levels = c("UNTREATED24hr", "DASATINIB24hr")) # Reorder untreated 1st



pdf("output/seurat/UMAP_UNTREATED24hr_DASATINIB24hr_test.pdf", width=10, height=6)

pdf("output/seurat/UMAP_UNTREATED24hr_DASATINIB24hr_25dim.pdf", width=10, height=6)

DimPlot(humangastruloid24hr.combined.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()




#### Conchi new list 20231214
Epiblast = c("POU5F1", "NANOG", "SOX2", "CDH1", "NODAL")
Primitive_Streak = c("TBXT", "EOMES", "MIXL1")
Nascent_Mesoderm = c("MESP1")
Endoderm =c("SOX17", "FOXA2")  # These cells will also be positive for Primitive streak markers, such as Eomes or T


# test_newconchi list 20231214
## 25 dim
DefaultAssay(humangastruloid24hr.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_NODAL_25dim.pdf", width=10, height=5)
FeaturePlot(humangastruloid24hr.combined.sct, features = c("NODAL"), split.by = "condition", max.cutoff = 5, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_Epiblast_List20231214_25dim.pdf", width=10, height=10)
FeaturePlot(humangastruloid24hr.combined.sct, features = Epiblast, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_Primitive_Streak_List20231214_25dim.pdf", width=5, height=5)
FeaturePlot(humangastruloid24hr.combined.sct, features = Primitive_Streak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_Nascent_Mesoderm_List20231214_25dim.pdf", width=5, height=5)
FeaturePlot(humangastruloid24hr.combined.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_Endoderm_List20231214_25dim.pdf", width=10, height=5)
FeaturePlot(humangastruloid24hr.combined.sct, features = Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()




# display features

pdf("output/seurat/UMAP_humangastruloid24hr_Phase.pdf", width=10, height=6)
DimPlot(humangastruloid24hr.combined.sct, reduction = "umap", split.by = "condition", label=TRUE, group.by = "Phase")
dev.off()
pdf("output/seurat/UMAP_humangastruloid24hr_nCount_RNA.pdf", width=10, height=6)
FeaturePlot(humangastruloid24hr.combined.sct, "nCount_RNA")
dev.off()
pdf("output/seurat/UMAP_humangastruloid24hr_percent.mt.pdf", width=10, height=6)
FeaturePlot(humangastruloid24hr.combined.sct, "percent.mt")
dev.off()
pdf("output/seurat/UMAP_humangastruloid24hr_nFeature_RNA.pdf", width=10, height=6)
FeaturePlot(humangastruloid24hr.combined.sct, "nFeature_RNA")
dev.off()






# Rename cluster
## 25 dim


new.cluster.ids <- c(
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast",
  "Endoderm",
  "Unknown"
)

names(new.cluster.ids) <- levels(humangastruloid24hr.combined.sct)
humangastruloid24hr.combined.sct <- RenameIdents(humangastruloid24hr.combined.sct, new.cluster.ids)

humangastruloid24hr.combined.sct$cluster.annot <- Idents(humangastruloid24hr.combined.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_humangastruloid24hr_label_25dim.pdf", width=12, height=6)
DimPlot(humangastruloid24hr.combined.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 5)
dev.off()

pdf("output/seurat/UMAP_humangastruloid24hr_label2_25dim.pdf", width=10, height=10)
DimPlot(humangastruloid24hr.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 5)
dev.off()



# SAve 
## saveRDS(humangastruloid24hr.combined.sct, file = "output/seurat/humangastruloid24hr.combined.sct_25dim.rds")
humangastruloid24hr.combined.sct = readRDS(file = "output/seurat/humangastruloid24hr.combined.sct_25dim.rds")





#overlapping condition
pdf("output/seurat/UMAP_humangastruloid24hr_label_overlap_25dim.pdf", width=6, height=5)
DimPlot(humangastruloid24hr.combined.sct, reduction = "umap", group.by = "condition", pt.size = 0.5, cols = c("blue","red"))
dev.off()



# All in dotplot
DefaultAssay(humangastruloid24hr.combined.sct) <- "SCT"

  "Epiblast",
  "Exe_Ectoderm",
  "Primitive_Streak-Nascent_Mesoderm_1",
  "Exe_Endoderm_1",
  "Primitive_Streak-Nascent_Mesoderm_2",
  "Endoderm",
  "Exe_Endoderm_2",
  "Blood_Progenitor"


#### Conchi new list 20231214

Epiblast = c("POU5F1", "NANOG", "SOX2", "CDH1", "NODAL")
Primitive_Streak = c("TBXT", "EOMES", "MIXL1")
Nascent_Mesoderm = c("MESP1")
Endoderm =c("SOX17", "FOXA2")  # These cells will also be positive for Primitive streak markers, such as Eomes or T



all_markers <- c(
"POU5F1", "NANOG", "SOX2", "CDH1", "NODAL",
"TBXT", "EOMES", "MIXL1",
"MESP1",
"SOX17", "FOXA2"
)


levels(humangastruloid24hr.combined.sct) <- c(
  "Unknown",
  "Endoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast"
)


pdf("output/seurat/DotPlot_humangastruloid24hr_ident_25dim.pdf", width=8, height=3)
DotPlot(humangastruloid24hr.combined.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()




## Downsampling with bootstrap to compare the nb of cell per cell types

library("tidyverse")

### Identify the unique clusters
unique_clusters <- unique(Idents(humangastruloid24hr.combined.sct))

### Create empty matrices to store cell counts
UNTREATED24hr_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
DASATINIB24hr_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(UNTREATED24hr_clusters_counts) <- unique_clusters
colnames(DASATINIB24hr_clusters_counts) <- unique_clusters

# count nb of cells in each cluster
untreated_clusters <- table(Idents(humangastruloid24hr.combined.sct)[humangastruloid24hr.combined.sct$condition == "UNTREATED24hr"])
print(untreated_clusters)
dasatinib_clusters <- table(Idents(humangastruloid24hr.combined.sct)[humangastruloid24hr.combined.sct$condition == "DASATINIB24hr"])
print(dasatinib_clusters)



### Loop through 100 iterations
humangastruloid24hr.combined.sct_UNTREATED24hr <- which(humangastruloid24hr.combined.sct$orig.ident == 'UNTREATED24hr') # 7331 cells
humangastruloid24hr.combined.sct_DASATINIB24hr <- which(humangastruloid24hr.combined.sct$orig.ident == 'DASATINIB24hr') # 6613 cells

for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  humangastruloid24hr.combined.sct_UNTREATED24hr_downsample <- sample(humangastruloid24hr.combined.sct_UNTREATED24hr, 6613)
  humangastruloid24hr.combined.sct_integrated_downsample <- humangastruloid24hr.combined.sct[,c(humangastruloid24hr.combined.sct_UNTREATED24hr_downsample, humangastruloid24hr.combined.sct_DASATINIB24hr)]

  # Count nb of cells in each cluster
  UNTREATED24hr_clusters <- table(Idents(humangastruloid24hr.combined.sct_integrated_downsample)[humangastruloid24hr.combined.sct_integrated_downsample$condition == "UNTREATED24hr"])
  DASATINIB24hr_clusters <- table(Idents(humangastruloid24hr.combined.sct_integrated_downsample)[humangastruloid24hr.combined.sct_integrated_downsample$condition == "DASATINIB24hr"])

  # Align the counts with the unique clusters
  UNTREATED24hr_clusters_counts[i, names(UNTREATED24hr_clusters)] <- as.numeric(UNTREATED24hr_clusters)
  DASATINIB24hr_clusters_counts[i, names(DASATINIB24hr_clusters)] <- as.numeric(DASATINIB24hr_clusters)
}



### Calculate mean and standard error
mean_UNTREATED24hr_clusters <- colMeans(UNTREATED24hr_clusters_counts)
mean_DASATINIB24hr_clusters <- colMeans(DASATINIB24hr_clusters_counts)
std_error_UNTREATED24hr_clusters <- apply(UNTREATED24hr_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("UNTREATED24hr", "DASATINIB24hr")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- UNTREATED24hr_clusters_counts[j,i]
    contingency_table[1,2] <- DASATINIB24hr_clusters_counts[j,i]
    contingency_table[2,1] <- sum(UNTREATED24hr_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(DASATINIB24hr_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_UNTREATED24hr_clusters),
  untreated = mean_UNTREATED24hr_clusters,
  dasatinib = mean_DASATINIB24hr_clusters,
  std_error_UNTREATED24hr = std_error_UNTREATED24hr_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_UNTREATED24hr, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "UNTREATED24hr", "DASATINIB24hr"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.05, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("UNTREATED24hr", "DASATINIB24hr")) # Reorder untreated 1st
plot_data$cluster <- factor(plot_data$cluster, levels = c(  "Unknown",
  "Endoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast")) # Reorder untreated 1st


# Plotting using ggplot2
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling100_clean_humangastruloid24hr_25dim.pdf", width=5, height=4)
ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = filter(plot_data, condition == "DASATINIB24hr"),
    aes(label = significance, y = value + std_error_UNTREATED24hr_clusters),
    vjust = -0.8,
    position = position_dodge(0.9), size = 5
  ) +
  scale_fill_manual(values = c("UNTREATED24hr" = "blue", "DASATINIB24hr" = "red")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,4000)
dev.off()




# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

DefaultAssay(humangastruloid24hr.combined.sct) <- "RNA"

humangastruloid24hr.combined.sct <- NormalizeData(humangastruloid24hr.combined.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(humangastruloid24hr.combined.sct)
humangastruloid24hr.combined.sct <- ScaleData(humangastruloid24hr.combined.sct, features = all.genes) # zero-centres and scales it


## what genes change in different conditions for cells of the same type

humangastruloid24hr.combined.sct$celltype.stim <- paste(humangastruloid24hr.combined.sct$cluster.annot, humangastruloid24hr.combined.sct$condition,
    sep = "-")
Idents(humangastruloid24hr.combined.sct) <- "celltype.stim"

# use RNA corrected count for DEGs
humangastruloid24hr.combined.sct <- PrepSCTFindMarkers(humangastruloid24hr.combined.sct)


## Automation::
cell_types <- c("Unknown",
  "Endoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast")

for (cell_type in cell_types) {
  response_name <- paste(cell_type, "DASATINIB24hr.response", sep = ".")
  ident_1 <- paste(cell_type, "-DASATINIB24hr", sep = "")
  ident_2 <- paste(cell_type, "-UNTREATED24hr", sep = "")

  response <- FindMarkers(humangastruloid24hr.combined.sct, assay = "RNA", ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
  
  print(head(response, n = 15))
  
  file_name <- paste("output/seurat/", cell_type, "-DASATINIB24hr_25dim.txt", sep = "")
  write.table(response, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE)
}



### Find all markers 
all_markers <- FindAllMarkers(humangastruloid24hr.combined.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_humangastruloid24hr_all_markers_25dim.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# Display the top 10 CONSERVED marker genes of each cluster
Idents(humangastruloid24hr.combined.sct) <- "cluster.annot"

## DEGs cluster versus all other
Unknown.conserved <- FindConservedMarkers(humangastruloid24hr.combined.sct, assay = "RNA", ident.1 = "Unknown", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknown")
Endoderm.conserved <- FindConservedMarkers(humangastruloid24hr.combined.sct, assay = "RNA", ident.1 = "Endoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Endoderm")
Nascent_Mesoderm.conserved <- FindConservedMarkers(humangastruloid24hr.combined.sct, assay = "RNA", ident.1 = "Nascent_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Nascent_Mesoderm")
Primitive_Streak.conserved <- FindConservedMarkers(humangastruloid24hr.combined.sct, assay = "RNA", ident.1 = "Primitive_Streak", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Primitive_Streak")
Epiblast.conserved <- FindConservedMarkers(humangastruloid24hr.combined.sct, assay = "RNA", ident.1 = "Epiblast", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Epiblast")


## Combine all conserved markers into one data frame
all_conserved <- bind_rows(Unknown.conserved, Endoderm.conserved, Nascent_Mesoderm.conserved, Primitive_Streak.conserved, Epiblast.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers_humangastruloid24hr_25dim.txt", sep = "\t", quote = FALSE, row.names = TRUE)

## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("Unknown",
  "Endoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("Unknown",
  "Endoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast")))



## Write the top 10 conserved markers for each cluster to a file
write.table(top10_conserved, file = "output/seurat/srat_top10_conserved_markers_humangastruloid24hr_25dim.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Visualize the top 10/3 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(humangastruloid24hr.combined.sct) <- c("Unknown",
  "Endoderm",
  "Nascent_Mesoderm",
  "Primitive_Streak",
  "Epiblast")

pdf("output/seurat/DotPlot_SCT_top10_conserved_humangastruloid24hr_25dim.pdf", width=18, height=3)
DotPlot(humangastruloid24hr.combined.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# save
## saveRDS(humangastruloid24hr.combined.sct, file = "output/seurat/humangastruloid24hr.combined.sct_25dim.rds")
humangastruloid24hr.combined.sct <- readRDS(file = "output/seurat/humangastruloid24hr.combined.sct_25dim.rds")



## Re-calculate DEGs keeping ALL genes
humangastruloid24hr.combined.sct$celltype.stim <- paste(humangastruloid24hr.combined.sct$cluster.annot, humangastruloid24hr.combined.sct$condition,
    sep = "-")
Idents(humangastruloid24hr.combined.sct) <- "celltype.stim"

Unknown <- FindMarkers(humangastruloid24hr.combined.sct, ident.1 = "Unknown-DASATINIB24hr", ident.2 = "Unknown-UNTREATED24hr",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Endoderm <- FindMarkers(humangastruloid24hr.combined.sct, ident.1 = "Endoderm-DASATINIB24hr", ident.2 = "Endoderm-UNTREATED24hr",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Nascent_Mesoderm <- FindMarkers(humangastruloid24hr.combined.sct, ident.1 = "Nascent_Mesoderm-DASATINIB24hr", ident.2 = "Nascent_Mesoderm-UNTREATED24hr",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Primitive_Streak <- FindMarkers(humangastruloid24hr.combined.sct, ident.1 = "Primitive_Streak-DASATINIB24hr", ident.2 = "Primitive_Streak-UNTREATED24hr",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Epiblast <- FindMarkers(humangastruloid24hr.combined.sct, ident.1 = "Epiblast-DASATINIB24hr", ident.2 = "Epiblast-UNTREATED24hr",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
    



### save output
#### write.table(Epiblast, file = "output/seurat/Epiblast-humangastruloid24hr_25dim_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)





# Check some genes
DefaultAssay(humangastruloid24hr.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

## Conchi geneList 20231214
pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_NODAL_25dim.pdf", width=8, height=12)
FeaturePlot(humangastruloid24hr.combined.sct, features = c("YAP1","NODAL", "AXIN2", "SP5"), split.by = "condition", max.cutoff = 5, cols = c("grey", "red"))
dev.off()


## YAP1-NODAL

YAP1_NODAL <- c("Shh", "Dmrt1", "Cited2", 'Dact2', "Tgif2", "Acvr1b", 'Smad3', 'Dand5', 'Smad2','Cer1', 'Nodal', 'Foxh1', 'Dact1', 'Tdgf1', 'Cripto3', 'Acvr1c', 'Cfc1', 'Cfc1b') 

YAP1_NODAL <- c("SHH", "DMRT1", "CITED2", 'DACT2', "TGIF2", "ACVR1B", 'SMAD3', 'DAND5', 'SMAD2','CER1', 'NODAL', 'FOXH1', 'DACT1', 'TDGF1', 'TDGF1P3', 'ACVR1C', 'CFC1', 'CFC1B') 

pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_25dim_YAP1_NODAL.pdf", width=10, height=70)
FeaturePlot(humangastruloid24hr.combined.sct, features = YAP1_NODAL, max.cutoff = 4, cols = c("grey", "red"), split.by = "condition")
dev.off()

## Anxa gene family
"Anxa1","Amotl2","Ccn1","Ccn2","Tagln","Anxa3","Chchd2"

pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_ANXA_25dim.pdf", width=10, height=25)
FeaturePlot(humangastruloid24hr.combined.sct, features = c("ANXA1","AMOTL2","CCN1","CCN2","TAGLN","ANXA3","CHCHD2"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()



## PRC2-related genes

pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_PRC2_25dim.pdf", width=10, height=22)
FeaturePlot(humangastruloid24hr.combined.sct, features = c("EZH2", "KDM6B", "GATA6", "CDX2"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


### Check YAP1 hippo bulk-regulated genes
"Tead1", "Tead4", "Ccn2", "Ccn1","Areg","Myc","Gli2","Vim","Axl","Birc5"

pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_YAP1_hippo_25dim.pdf", width=7, height=30)
FeaturePlot(humangastruloid24hr.combined.sct, features = c("TEAD1", "TEAD4", "CCN2", "CCN1","AREG","MYC","GLI2","VIM","AXL","BIRC5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()



pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_YAP1_25dim.pdf", width=10, height=5)

FeaturePlot(humangastruloid24hr.combined.sct, features = c("YAP1"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()


## Wnt related

pdf("output/seurat/FeaturePlot_SCT_humangastruloid24hr_Wnt_25dim.pdf", width=8, height=7)

FeaturePlot(humangastruloid24hr.combined.sct, features = c("AXIN2", "SP5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()


```



# Embryo E7 - E7.75 data integration


Let's integrate E7 E7.75 dataset into a single dataframe.

Can use (seurat v5 data integration method)[https://satijalab.org/seurat/articles/seurat5_integration]; they integrate together data from different sequencing platforms; same cell type. 


Let's try using our previous method, and simply add in the list of srat; all samples:
- clean data E7 and E7.75 (same parameters as before)
- integrate all samples
- separate plot per WT and KO and color cells for E7 and E7.75

*NOTE: these methods are mostly to integrate data from the same condition! So same time-point, as it tries to "force" cells that are probably the same cell type to cluster closer together : "These methods aim to identify shared cell states that are present across different datasets, even if they were collected from different individuals, experimental conditions, technologies, or even species" from (here)[https://satijalab.org/seurat/archive/v3.0/integration.html]. But let's try! As using this we identify cluster only present in KO sometime!*

```bash
conda activate scRNAseqV3
```

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
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

# E7 data

## Load the matrix and Create SEURAT object
load("output/soupX/out_E7mousecontrol.RData")
srat_WT_E7 <- CreateSeuratObject(counts = out, project = "WT_E7") # 32,285 features across 1,829 samples

load("output/soupX/out_E7mousecYAPKO.RData")
srat_cYAPKO_E7 <- CreateSeuratObject(counts = out, project = "cYAPKO_E7") # 32,285 features across 1,897 samples

# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
srat_WT_E7[["percent.mt"]] <- PercentageFeatureSet(srat_WT_E7, pattern = "^mt-")
srat_WT_E7[["percent.rb"]] <- PercentageFeatureSet(srat_WT_E7, pattern = "^Rp[sl]")

srat_cYAPKO_E7[["percent.mt"]] <- PercentageFeatureSet(srat_cYAPKO_E7, pattern = "^mt-")
srat_cYAPKO_E7[["percent.rb"]] <- PercentageFeatureSet(srat_cYAPKO_E7, pattern = "^Rp[sl]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/embryo_E7_control.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_WT_E7 <- AddMetaData(srat_WT_E7,doublets)
srat_WT_E7$Doublet_score <- as.numeric(srat_WT_E7$Doublet_score) # make score as numeric
head(srat_WT_E7[[]])

doublets <- read.table("output/doublets/embryo_E7_cYAPKO.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_cYAPKO_E7 <- AddMetaData(srat_cYAPKO_E7,doublets)
srat_cYAPKO_E7$Doublet_score <- as.numeric(srat_cYAPKO_E7$Doublet_score) # make score as numeric
head(srat_cYAPKO_E7[[]])


### V2 more stringeant
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC == 'Pass','Low_nFeature',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$percent.mt > 10 & srat_WT_E7@meta.data$QC == 'Pass','High_MT',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
table(srat_WT_E7[['QC']])
## 
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC == 'Pass','Low_nFeature',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$percent.mt > 10 & srat_cYAPKO_E7@meta.data$QC == 'Pass','High_MT',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
table(srat_cYAPKO_E7[['QC']])



## subset my seurat object to only analyze the cells that pass the QC
srat_WT_E7 <- subset(srat_WT_E7, subset = QC == 'Pass')
srat_cYAPKO_E7 <- subset(srat_cYAPKO_E7, subset = QC == 'Pass')
srat_WT_E7$condition <- "WT"
srat_cYAPKO_E7$condition <- "cYAPKO"



mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
srat_WT_E7 <- NormalizeData(srat_WT_E7, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_WT_E7)
srat_WT_E7 <- ScaleData(srat_WT_E7, features = all.genes) # zero-centres and scales it

srat_cYAPKO_E7 <- NormalizeData(srat_cYAPKO_E7, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_cYAPKO_E7)
srat_cYAPKO_E7 <- ScaleData(srat_cYAPKO_E7, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
srat_WT_E7 <- CellCycleScoring(srat_WT_E7, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_WT_E7[[]]$Phase)
srat_cYAPKO_E7 <- CellCycleScoring(srat_cYAPKO_E7, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_cYAPKO_E7[[]]$Phase)

set.seed(42)


# E7.75 data


## Load the matrix and Create SEURAT object
load("output/soupX/out_embryo_control_e775_mice.RData")
srat_WT_E775 <- CreateSeuratObject(counts = out, project = "WT") # 32,285 features across 5,568 samples

load("output/soupX/out_embryo_cYAPKO_e775_mice.RData")
srat_cYAPKO_E775 <- CreateSeuratObject(counts = out, project = "cYAPKO") # 32,285 features across 4,504 samples


# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
srat_WT_E775[["percent.mt"]] <- PercentageFeatureSet(srat_WT_E775, pattern = "^mt-")
srat_WT_E775[["percent.rb"]] <- PercentageFeatureSet(srat_WT_E775, pattern = "^Rp[sl]")

srat_cYAPKO_E775[["percent.mt"]] <- PercentageFeatureSet(srat_cYAPKO_E775, pattern = "^mt-")
srat_cYAPKO_E775[["percent.rb"]] <- PercentageFeatureSet(srat_cYAPKO_E775, pattern = "^Rp[sl]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/embryo_control.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_WT_E775 <- AddMetaData(srat_WT_E775,doublets)
srat_WT_E775$Doublet_score <- as.numeric(srat_WT_E775$Doublet_score) # make score as numeric
head(srat_WT_E775[[]])

doublets <- read.table("output/doublets/embryo_cYAPKO.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_cYAPKO_E775 <- AddMetaData(srat_cYAPKO_E775,doublets)
srat_cYAPKO_E775$Doublet_score <- as.numeric(srat_cYAPKO_E775$Doublet_score) # make score as numeric
head(srat_cYAPKO_E775[[]])


## After seeing the plot; add QC information in our seurat object
## V1 QC; optimal with vst V1; dim 19 k param 70 es 0.9 --> THE WINNER pro winner
srat_WT_E775[['QC']] <- ifelse(srat_WT_E775@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_WT_E775[['QC']] <- ifelse(srat_WT_E775@meta.data$nFeature_RNA < 2000 & srat_WT_E775@meta.data$QC == 'Pass','Low_nFeature',srat_WT_E775@meta.data$QC)
srat_WT_E775[['QC']] <- ifelse(srat_WT_E775@meta.data$nFeature_RNA < 2000 & srat_WT_E775@meta.data$QC != 'Pass' & srat_WT_E775@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_WT_E775@meta.data$QC,sep = ','),srat_WT_E775@meta.data$QC)
srat_WT_E775[['QC']] <- ifelse(srat_WT_E775@meta.data$percent.mt > 25 & srat_WT_E775@meta.data$QC == 'Pass','High_MT',srat_WT_E775@meta.data$QC)
srat_WT_E775[['QC']] <- ifelse(srat_WT_E775@meta.data$nFeature_RNA < 2000 & srat_WT_E775@meta.data$QC != 'Pass' & srat_WT_E775@meta.data$QC != 'High_MT',paste('High_MT',srat_WT_E775@meta.data$QC,sep = ','),srat_WT_E775@meta.data$QC)
table(srat_WT_E775[['QC']])
## 
srat_cYAPKO_E775[['QC']] <- ifelse(srat_cYAPKO_E775@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_cYAPKO_E775[['QC']] <- ifelse(srat_cYAPKO_E775@meta.data$nFeature_RNA < 2000 & srat_cYAPKO_E775@meta.data$QC == 'Pass','Low_nFeature',srat_cYAPKO_E775@meta.data$QC)
srat_cYAPKO_E775[['QC']] <- ifelse(srat_cYAPKO_E775@meta.data$nFeature_RNA < 2000 & srat_cYAPKO_E775@meta.data$QC != 'Pass' & srat_cYAPKO_E775@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_cYAPKO_E775@meta.data$QC,sep = ','),srat_cYAPKO_E775@meta.data$QC)
srat_cYAPKO_E775[['QC']] <- ifelse(srat_cYAPKO_E775@meta.data$percent.mt > 25 & srat_cYAPKO_E775@meta.data$QC == 'Pass','High_MT',srat_cYAPKO_E775@meta.data$QC)
srat_cYAPKO_E775[['QC']] <- ifelse(srat_cYAPKO_E775@meta.data$nFeature_RNA < 2000 & srat_cYAPKO_E775@meta.data$QC != 'Pass' & srat_cYAPKO_E775@meta.data$QC != 'High_MT',paste('High_MT',srat_cYAPKO_E775@meta.data$QC,sep = ','),srat_cYAPKO_E775@meta.data$QC)
table(srat_cYAPKO_E775[['QC']])




## subset my seurat object to only analyze the cells that pass the QC
srat_WT_E775 <- subset(srat_WT_E775, subset = QC == 'Pass')
srat_cYAPKO_E775 <- subset(srat_cYAPKO_E775, subset = QC == 'Pass')
srat_WT_E775$condition <- "WT"
srat_cYAPKO_E775$condition <- "cYAPKO"



mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
srat_WT_E775 <- NormalizeData(srat_WT_E775, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_WT_E775)
srat_WT_E775 <- ScaleData(srat_WT_E775, features = all.genes) # zero-centres and scales it

srat_cYAPKO_E775 <- NormalizeData(srat_cYAPKO_E775, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_cYAPKO_E775)
srat_cYAPKO_E775 <- ScaleData(srat_cYAPKO_E775, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
srat_WT_E775 <- CellCycleScoring(srat_WT_E775, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_WT_E775[[]]$Phase)
srat_cYAPKO_E775 <- CellCycleScoring(srat_cYAPKO_E775, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_cYAPKO_E775[[]]$Phase)

set.seed(42)




# Both E7 and E7.75 used 19 dims
srat_WT_E7 <- SCTransform(srat_WT_E7, method = "glmGamPoi", ncells = 788, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO_E7 <- SCTransform(srat_cYAPKO_E7, method = "glmGamPoi", ncells = 862, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)
srat_WT_E775 <- SCTransform(srat_WT_E775, method = "glmGamPoi", ncells = 4812, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 19, verbose = FALSE)
srat_cYAPKO_E775 <- SCTransform(srat_cYAPKO_E775, method = "glmGamPoi", ncells = 3621, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 19, verbose = FALSE)

### Add time as factor:
srat_WT_E775$time <- "E775"
srat_cYAPKO_E775$time <- "E775"
srat_WT_E7$time <- "E7"
srat_cYAPKO_E7$time <- "E7"

# Data integration (check active assay is 'SCT')
srat.list <- list(srat_WT_E7 = srat_WT_E7, srat_WT_E775 = srat_WT_E775, srat_cYAPKO_E7 = srat_cYAPKO_E7, srat_cYAPKO_E775 = srat_cYAPKO_E775)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

embryo.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
embryo.combined.E7E775.sct <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT")

set.seed(42)

## cluster all samples
DefaultAssay(embryo.combined.E7E775.sct) <- "integrated"
embryo.combined.E7E775.sct <- RunPCA(embryo.combined.E7E775.sct, verbose = FALSE, npcs = 19)
embryo.combined.E7E775.sct <- RunUMAP(embryo.combined.E7E775.sct, reduction = "pca", dims = 1:19, verbose = FALSE)
embryo.combined.E7E775.sct <- FindNeighbors(embryo.combined.E7E775.sct, reduction = "pca", k.param = 25, dims = 1:19)
embryo.combined.E7E775.sct <- FindClusters(embryo.combined.E7E775.sct, resolution = 1.5, verbose = FALSE, algorithm = 4)
embryo.combined.E7E775.sct$condition <- factor(embryo.combined.E7E775.sct$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st
pdf("output/seurat/UMAP_embryo_E7E775.pdf", width=10, height=6)
DimPlot(embryo.combined.E7E775.sct, reduction = "umap", split.by = "condition", label=TRUE)
dev.off()



pdf("output/seurat/UMAP_embryo_E7E775_time.pdf", width=6, height=6)
DimPlot(embryo.combined.E7E775.sct, reduction = "umap", group.by = "time", label=TRUE, pt.size = 0.5, shuffle = TRUE, label.size = 6, cols = c('E7' = '#78C850', 'E775' = '#FDB813') )
dev.off()

pdf("output/seurat/UMAP_embryo_E7E775_timeSplit.pdf", width=15, height=6)
DimPlot(embryo.combined.E7E775.sct, reduction = "umap", split.by = "time", label=TRUE, pt.size = 0.5, shuffle = TRUE, label.size = 2.5 )
dev.off()


# cell type assignment

### E775 markers (Updated marker gene list from Alex + WINNER pro-winner)
Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # 
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l") # 
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # 
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # 
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # 
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # 
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # 
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # 
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # 
Caudal_Neurectoderm = c("Olig3","Hes3","Lrrc7","Cntfr") # 
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # 
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # 
Intermediate_Mesoderm = c("Bik","Arg1","Meox1","Ism1") # 
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # 
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # 
Notocord = c("Shh","Noto","Vtn","Fgg") #  
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # 
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3") # 


#### E7 markers (Conchi new list 20231214)
Epiblast = c("Pou5f1", "Sox2", "Cdh1")
Primitive_Streak = c("Mixl1", "Eomes", "T")
Blood_Progenitor = c("Gata1", "Tal1", "Kdr", "Runx1", "Ets2", "Etv2")
Nascent_Mesoderm = c("Mesp1")
Endoderm =c("Sox17", "Foxa2")  # These cells will also be positive for Primitive streak markers, such as Eomes or T
ExE_Endoderm = c("Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6") # This cells will be negative or have low levels of Primitive streak markers
ExE_Ectoderm =c("Elf5", "Rhox5")

#### combine markers:
# Epiblast_PrimStreak = c("Hesx1","Otx2","Pax2","Otx2os1") # 
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l", "Rhox5")
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # 
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # 
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # 
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # 
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # 
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # 
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # 
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # 
Blood_Progenitor_E7 = c("Tal1", "Kdr", "Runx1", "Ets2", "Etv2") # removed Gata1
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # 
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # 
Notocord = c("Shh","Noto","Vtn","Fgg") #  
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # 
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3") # 
Epiblast = c("Pou5f1", "Sox2", "Cdh1")
Primitive_Streak = c("Mixl1", "Eomes", "T")
Endoderm =c("Sox17", "Foxa2")  # These cells will also be positive for Primitive streak markers, such as Eomes or T
ExE_Endoderm = c("Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6") # This cells will be negative or have low levels of Primitive streak markers
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # 

#######################
# SAVE: saveRDS(embryo.combined.E7E775.sct, file = "output/seurat/embryo.combined.E7E775.sct.rds")
# LOAD: embryo.combined.E7E775.sct <- readRDS(file = "output/seurat/embryo.combined.E7E775.sct.rds")
#######################

DefaultAssay(embryo.combined.E7E775.sct) <- "SCT" # 



# expression of marker genes
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-ExE_Ectoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = ExE_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Nascent_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Nascent_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Somitic_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Somitic_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Caudal_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Caudal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Haematodenothelial_progenitors.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Haematodenothelial_progenitors, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Paraxial_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Paraxial_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Mesenchyme.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Mesenchyme, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Blood_Progenitor_1.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Blood_Progenitor_1, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Blood_Progenitor_2.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Blood_Progenitor_2, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Blood_Progenitor_E7.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Blood_Progenitor_E7, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Surface_Ectoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Surface_Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Mixed_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Mixed_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Notocord.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Notocord, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Gut.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Gut, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Primordial_Germ_Cells.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Primordial_Germ_Cells, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Epiblast.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Epiblast, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Primitive_Streak.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Primitive_Streak, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Endoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-ExE_Endoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = ExE_Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Pharyngeal_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = Pharyngeal_Mesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775-Primitive_Streak_Cardiac_Mesoderm.pdf", width=10, height=10)
FeaturePlot(embryo.combined.E7E775.sct, features = c("Eomes", "Fgf8", "Nodal","Hand1"), max.cutoff = 3, cols = c("grey", "red"))
dev.off()



new.cluster.ids <- c(
  "Epiblast_1", 
  "Caudal_Mesoderm",
  "ExE_Ectoderm_1",
  "Somitic_Mesoderm",
  "Nascent_Mesoderm",
  "Haematodenothelial_progenitors",
  "Mesenchyme",
  "Unknown_1",
  "Mixed_Mesoderm",
  "Epiblast_2",
  "Blood_Progenitor_1",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "ExE_Endoderm",
  "Paraxial_Mesoderm",
  "ExE_Ectoderm_2",
  "Unknown_2",
  "Endoderm",
  "Pharyngeal_Mesoderm",
  "Gut",
  "Unknown_3",
  "Notocord",
  "Primordial_Germ_Cells",
  "Blood_Progenitor_3"
)


names(new.cluster.ids) <- levels(embryo.combined.E7E775.sct)
embryo.combined.E7E775.sct <- RenameIdents(embryo.combined.E7E775.sct, new.cluster.ids)

embryo.combined.E7E775.sct$cluster.annot <- Idents(embryo.combined.E7E775.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_embryo_E7E775_label.pdf", width=12, height=6)
DimPlot(embryo.combined.E7E775.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_embryo_E7E775_label_noSplit.pdf", width=10, height=6)
DimPlot(embryo.combined.E7E775.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()


#overlapping condition
pdf("output/seurat/UMAP_embryo_E7E775_label_overlap.pdf", width=6, height=5)
DimPlot(embryo.combined.E7E775.sct, reduction = "umap", group.by = "condition", pt.size = 0.000001, cols = c("blue","red"))
dev.off()



# All in dotplot
DefaultAssay(embryo.combined.E7E775.sct) <- "SCT"

Epiblast = c("Pou5f1", "Sox2", "Cdh1")
Caudal_Mesoderm = c("Wnt3a","Nkx1-2","Apln","Rxrg") # 
ExE_Ectoderm = c("Tex19.1","Elf5","Wnt7b","Dnmt3l", "Rhox5")
Somitic_Mesoderm = c("Meox1","Aldh1a2","Synm","Pcdh8") # 
Nascent_Mesoderm = c("Col9a1","Mesp1","Osr1") # 
Haematodenothelial_progenitors = c("Hoxa9","Hoxa10","Tbx4","Pitx1") # 
Mesenchyme = c("Tdo2","Adamts2","Colec11","Snai2") # 
Mixed_Mesoderm = c("Tbx6","Dll1","Hes7","Fgf17") # 
Blood_Progenitor_1 = c("Sox18","Esam","Rhoj","Flt4") # 
Surface_Ectoderm = c("Tfap2a","Npnt","Wnt4","Id4") # 
Blood_Progenitor_2 = c("Gypa","Gata1","Cited4","Gfi1b") # 
ExE_Endoderm = c("Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6") # This cells will be negative or have low levels of Primitive streak markers
Paraxial_Mesoderm = c("Ebf2","Ptgfr","Ebf3","Col23a1") # 
Endoderm =c("Sox17", "Foxa2")  # These cells will also be positive for Primitive streak markers, such as Eomes or T
Pharyngeal_Mesoderm = c("Nkx2-5","Tbx5","Mef2c","Myocd") # 
Gut = c("Pga5","Hs3cst1","Wfdc1","Islr2") # 
Notocord = c("Shh","Noto","Vtn","Fgg") #  
Primordial_Germ_Cells = c("Dppa3","Sprr2a3","Irf1","Ifitm3") # 
Primitive_Streak = c("Mixl1", "Eomes", "T")
Blood_Progenitor_E7 = c("Tal1", "Kdr", "Runx1", "Ets2", "Etv2") # removed Gata1


all_markers <- c(
  "Pou5f1", "Sox2", "Cdh1",
  "Wnt3a","Nkx1-2","Apln","Rxrg",
  "Tex19.1","Elf5","Wnt7b","Dnmt3l", "Rhox5",
  "Meox1","Aldh1a2","Synm","Pcdh8",
  "Col9a1","Mesp1","Osr1",
  "Hoxa9","Hoxa10","Tbx4","Pitx1",
  "Tdo2","Adamts2","Colec11","Snai2",
  "Tbx6","Dll1","Hes7","Fgf17",
  "Sox18","Esam","Rhoj","Flt4",
  "Tfap2a","Npnt","Wnt4","Id4",
  "Gypa","Gata1","Cited4","Gfi1b",
  "Ttr", "Dab2", "Sox7", "Fabp1", "Fabp2", "Gata4", "Cldn6",
  "Ebf2","Ptgfr","Ebf3","Col23a1",
  "Sox17", "Foxa2",
  "Nkx2-5","Tbx5","Mef2c","Myocd",
  "Pga5","Hs3cst1","Wfdc1","Islr2",
  "Shh","Noto","Vtn","Fgg",
  "Dppa3","Sprr2a3","Irf1","Ifitm3",
  "Mixl1", "Eomes", "T",
  "Tal1", "Kdr", "Runx1", "Ets2", "Etv2"
)


levels(embryo.combined.E7E775.sct) <- c(
  "Epiblast_1", 
  "Caudal_Mesoderm",
  "ExE_Ectoderm_1",
  "Somitic_Mesoderm",
  "Nascent_Mesoderm",
  "Haematodenothelial_progenitors",
  "Mesenchyme",
  "Unknown_1",
  "Mixed_Mesoderm",
  "Epiblast_2",
  "Blood_Progenitor_1",
  "Surface_Ectoderm",
  "Blood_Progenitor_2",
  "ExE_Endoderm",
  "Paraxial_Mesoderm",
  "ExE_Ectoderm_2",
  "Unknown_2",
  "Endoderm",
  "Pharyngeal_Mesoderm",
  "Gut",
  "Unknown_3",
  "Notocord",
  "Primordial_Germ_Cells",
  "Blood_Progenitor_3"
)
pdf("output/seurat/DotPlot_SCT_embryo_E7E775_ident.pdf", width=16.5, height=4.5)
DotPlot(embryo.combined.E7E775.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()



## Downsampling with bootstrap to compare the nb of cell per cell types

library("tidyverse")

### Identify the unique clusters
unique_clusters <- unique(Idents(embryo.combined.E7E775.sct))

### Create empty matrices to store cell counts
control_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
cYAPKO_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(control_clusters_counts) <- unique_clusters
colnames(cYAPKO_clusters_counts) <- unique_clusters

### Loop through 100 iterations
embryo.combined.E7E775.sct_WT <- which(embryo.combined.E7E775.sct$orig.ident == 'WT')
embryo.combined.E7E775.sct_cYAPKO <- which(embryo.combined.E7E775.sct$orig.ident == 'cYAPKO')

for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  embryo.combined.E7E775.sct_WT_downsample <- sample(embryo.combined.E7E775.sct_WT, 4483)
  embryo.combined.E7E775.sct_integrated_downsample <- embryo.combined.E7E775.sct[,c(embryo.combined.E7E775.sct_cYAPKO, embryo.combined.E7E775.sct_WT_downsample)]

  # Count nb of cells in each cluster
  control_clusters <- table(Idents(embryo.combined.E7E775.sct_integrated_downsample)[embryo.combined.E7E775.sct_integrated_downsample$condition == "WT"])
  cYAPKO_clusters <- table(Idents(embryo.combined.E7E775.sct_integrated_downsample)[embryo.combined.E7E775.sct_integrated_downsample$condition == "cYAPKO"])

  # Align the counts with the unique clusters
  control_clusters_counts[i, names(control_clusters)] <- as.numeric(control_clusters)
  cYAPKO_clusters_counts[i, names(cYAPKO_clusters)] <- as.numeric(cYAPKO_clusters)
}



### Calculate mean and standard error
mean_control_clusters <- colMeans(control_clusters_counts)
mean_cYAPKO_clusters <- colMeans(cYAPKO_clusters_counts)
std_error_WT_clusters <- apply(control_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("WT", "cYAPKO")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- control_clusters_counts[j,i]
    contingency_table[1,2] <- cYAPKO_clusters_counts[j,i]
    contingency_table[2,1] <- sum(control_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(cYAPKO_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_control_clusters),
  untreated = mean_control_clusters,
  dasatinib = mean_cYAPKO_clusters,
  std_error_WT = std_error_WT_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_WT, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "WT", "cYAPKO"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.05, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("WT", "cYAPKO")) # Reorder untreated 1st
plot_data$cluster <- factor(plot_data$cluster, levels = c(
    "Epiblast_1",
    "Epiblast_2",
    "Unknown_3" ,
    "Nascent_Mesoderm",
    "Somitic_Mesoderm",
    "Unknown_2" ,
    "Caudal_Mesoderm",
    "Pharyngeal_Mesoderm",
    "Unknown_1" ,
    "Paraxial_Mesoderm",
    "Mixed_Mesoderm",
    "Notocord",
    "Haematodenothelial_progenitors",
    "Blood_Progenitor_1",
    "Blood_Progenitor_2",
    "Blood_Progenitor_3",
    "Endoderm",
    "ExE_Endoderm",
    "Mesenchyme",
    "Surface_Ectoderm",
    "ExE_Ectoderm_1",
    "ExE_Ectoderm_2",
    "Primordial_Germ_Cells",
    "Gut"
)) # Reorder untreated 1st


# Plotting using ggplot2
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling100_embryo_E7E775.pdf", width=9, height=4)
ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = filter(plot_data, condition == "cYAPKO"),
    aes(label = significance, y = value + std_error_WT_clusters),
    vjust = -0.8,
    position = position_dodge(0.9), size = 5
  ) +
  scale_fill_manual(values = c("WT" = "blue", "cYAPKO" = "red")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

DefaultAssay(embryo.combined.E7E775.sct) <- "RNA"

embryo.combined.E7E775.sct <- NormalizeData(embryo.combined.E7E775.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(embryo.combined.E7E775.sct)
embryo.combined.E7E775.sct <- ScaleData(embryo.combined.E7E775.sct, features = all.genes) # zero-centres and scales it


## what genes change in different conditions for cells of the same type

embryo.combined.E7E775.sct$celltype.stim <- paste(embryo.combined.E7E775.sct$cluster.annot, embryo.combined.E7E775.sct$condition,
    sep = "-")
Idents(embryo.combined.E7E775.sct) <- "celltype.stim"

# use RNA corrected count for DEGs
embryo.combined.E7E775.sct <- PrepSCTFindMarkers(embryo.combined.E7E775.sct)


## Automation::
cell_types <- c(
  "Epiblast_1",
    "Epiblast_2",
    "Unknown_3" ,
    "Nascent_Mesoderm",
    "Somitic_Mesoderm",
    "Unknown_2" ,
    "Caudal_Mesoderm",
    "Pharyngeal_Mesoderm",
    "Unknown_1" ,
    "Paraxial_Mesoderm",
    "Mixed_Mesoderm",
    "Notocord",
    "Haematodenothelial_progenitors",
    "Blood_Progenitor_1",
    "Blood_Progenitor_2",
   # "Blood_Progenitor_3", Removed as not in WT
    "Endoderm",
    "ExE_Endoderm",
    "Mesenchyme",
    "Surface_Ectoderm",
    "ExE_Ectoderm_1",
    "ExE_Ectoderm_2",
    "Primordial_Germ_Cells",
    "Gut")

for (cell_type in cell_types) {
  response_name <- paste(cell_type, "cYAPKO.response", sep = ".")
  ident_1 <- paste(cell_type, "-cYAPKO", sep = "")
  ident_2 <- paste(cell_type, "-WT", sep = "")

  response <- FindMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = ident_1, ident.2 = ident_2, verbose = FALSE)
  
  print(head(response, n = 15))
  
  file_name <- paste("output/seurat/", cell_type, "-cYAPKO_response_E7E775.txt", sep = "")
  write.table(response, file = file_name, sep = "\t", quote = FALSE, row.names = TRUE)
}



### Find all markers 
all_markers <- FindAllMarkers(embryo.combined.E7E775.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_embryo_E7E775_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# Display the top 10 CONSERVED marker genes of each cluster
Idents(embryo.combined.E7E775.sct) <- "cluster.annot"

## DEGs cluster versus all other
Primordial_Germ_Cells.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Primordial_Germ_Cells", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Primordial_Germ_Cells")
Epiblast_1.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Epiblast_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Epiblast_1")
Epiblast_2.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Epiblast_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Epiblast_2")
Unknown_3.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Unknown_3", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknown_3")
Nascent_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Nascent_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Nascent_Mesoderm")
Somitic_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Somitic_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Somitic_Mesoderm")

Unknown_2.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Unknown_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknown_2")
Caudal_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Caudal_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Caudal_Mesoderm")
Pharyngeal_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Pharyngeal_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Pharyngeal_Mesoderm")
Unknown_1.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Unknown_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Unknown_1")
Paraxial_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Paraxial_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Paraxial_Mesoderm")
Mixed_Mesoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Mixed_Mesoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mixed_Mesoderm")
Notocord.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Notocord", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Notocord")
Haematodenothelial_progenitors.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Haematodenothelial_progenitors", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Haematodenothelial_progenitors")

Blood_Progenitor_1.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Blood_Progenitor_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_1")
Blood_Progenitor_2.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Blood_Progenitor_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_2")
Blood_Progenitor_3.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Blood_Progenitor_3", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Blood_Progenitor_3")
Endoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Endoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Endoderm")
ExE_Endoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "ExE_Endoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "ExE_Endoderm")
Mesenchyme.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Mesenchyme", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Mesenchyme")
Surface_Ectoderm.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Surface_Ectoderm", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Surface_Ectoderm")
ExE_Ectoderm_1.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "ExE_Ectoderm_1", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "ExE_Ectoderm_1")
ExE_Ectoderm_2.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "ExE_Ectoderm_2", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "ExE_Ectoderm_2")
Primordial_Germ_Cells.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Primordial_Germ_Cells", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Primordial_Germ_Cells")
Gut.conserved <- FindConservedMarkers(embryo.combined.E7E775.sct, assay = "RNA", ident.1 = "Gut", grouping.var = "condition", verbose = TRUE) %>% mutate(cluster = "Gut")




#######################
# SAVE: saveRDS(embryo.combined.E7E775.sct, file = "output/seurat/embryo.combined.E7E775.sct.rds")
# LOAD: embryo.combined.E7E775.sct <- readRDS(file = "output/seurat/embryo.combined.E7E775.sct.rds")
#######################





## Combine all conserved markers into one data frame
all_conserved <- bind_rows(Primordial_Germ_Cells.conserved,Epiblast_1.conserved,Epiblast_2.conserved,Unknown_3.conserved,Nascent_Mesoderm.conserved,Somitic_Mesoderm.conserved,Unknown_2.conserved,Caudal_Mesoderm.conserved,Pharyngeal_Mesoderm.conserved,Unknown_1.conserved,Paraxial_Mesoderm.conserved,Mixed_Mesoderm.conserved,Notocord.conserved,Haematodenothelial_progenitors.conserved,Blood_Progenitor_1.conserved,Blood_Progenitor_2.conserved,Blood_Progenitor_3.conserved,Endoderm.conserved,ExE_Endoderm.conserved,Mesenchyme.conserved,Surface_Ectoderm.conserved,ExE_Ectoderm_1.conserved,ExE_Ectoderm_2.conserved,Primordial_Germ_Cells.conserved,Gut.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers_embryo_E7E775.txt", sep = "\t", quote = FALSE, row.names = TRUE)




## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c(
    "Epiblast_1",
    "Epiblast_2",
    "Unknown_3" ,
    "Nascent_Mesoderm",
    "Somitic_Mesoderm",
    "Unknown_2" ,
    "Caudal_Mesoderm",
    "Pharyngeal_Mesoderm",
    "Unknown_1" ,
    "Paraxial_Mesoderm",
    "Mixed_Mesoderm",
    "Notocord",
    "Haematodenothelial_progenitors",
    "Blood_Progenitor_1",
    "Blood_Progenitor_2",
    "Blood_Progenitor_3",
    "Endoderm",
    "ExE_Endoderm",
    "Mesenchyme",
    "Surface_Ectoderm",
    "ExE_Ectoderm_1",
    "ExE_Ectoderm_2",
    "Primordial_Germ_Cells",
    "Gut"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c(
      "Epiblast_1",
    "Epiblast_2",
    "Unknown_3" ,
    "Nascent_Mesoderm",
    "Somitic_Mesoderm",
    "Unknown_2" ,
    "Caudal_Mesoderm",
    "Pharyngeal_Mesoderm",
    "Unknown_1" ,
    "Paraxial_Mesoderm",
    "Mixed_Mesoderm",
    "Notocord",
    "Haematodenothelial_progenitors",
    "Blood_Progenitor_1",
    "Blood_Progenitor_2",
    "Blood_Progenitor_3",
    "Endoderm",
    "ExE_Endoderm",
    "Mesenchyme",
    "Surface_Ectoderm",
    "ExE_Ectoderm_1",
    "ExE_Ectoderm_2",
    "Primordial_Germ_Cells",
    "Gut")))


## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c(
    "Epiblast_1",
    "Epiblast_2",
    "Unknown_3" ,
    "Nascent_Mesoderm",
    "Somitic_Mesoderm",
    "Unknown_2" ,
    "Caudal_Mesoderm",
    "Pharyngeal_Mesoderm",
    "Unknown_1" ,
    "Paraxial_Mesoderm",
    "Mixed_Mesoderm",
    "Notocord",
    "Haematodenothelial_progenitors",
    "Blood_Progenitor_1",
    "Blood_Progenitor_2",
    "Blood_Progenitor_3",
    "Endoderm",
    "ExE_Endoderm",
    "Mesenchyme",
    "Surface_Ectoderm",
    "ExE_Ectoderm_1",
    "ExE_Ectoderm_2",
    "Primordial_Germ_Cells",
    "Gut"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 10) %>% 
  ungroup() %>% 
  arrange(match(cluster, c(


## Find the top 10 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c(
    "Epiblast_1",
    "Epiblast_2",
    "Unknown_3" ,
    "Nascent_Mesoderm",
    "Somitic_Mesoderm",
    "Unknown_2" ,
    "Caudal_Mesoderm",
    "Pharyngeal_Mesoderm",
    "Unknown_1" ,
    "Paraxial_Mesoderm",
    "Mixed_Mesoderm",
    "Notocord",
    "Haematodenothelial_progenitors",
    "Blood_Progenitor_1",
    "Blood_Progenitor_2",
    "Blood_Progenitor_3",
    "Endoderm",
    "ExE_Endoderm",
    "Mesenchyme",
    "Surface_Ectoderm",
    "ExE_Ectoderm_1",
    "ExE_Ectoderm_2",
    "Primordial_Germ_Cells",
    "Gut"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 3) %>% 
  ungroup() %>% 
  arrange(match(cluster, c(
    "Epiblast_1",
    "Epiblast_2",
    "Unknown_3" ,
    "Nascent_Mesoderm",
    "Somitic_Mesoderm",
    "Unknown_2" ,
    "Caudal_Mesoderm",
    "Pharyngeal_Mesoderm",
    "Unknown_1" ,
    "Paraxial_Mesoderm",
    "Mixed_Mesoderm",
    "Notocord",
    "Haematodenothelial_progenitors",
    "Blood_Progenitor_1",
    "Blood_Progenitor_2",
    "Blood_Progenitor_3",
    "Endoderm",
    "ExE_Endoderm",
    "Mesenchyme",
    "Surface_Ectoderm",
    "ExE_Ectoderm_1",
    "ExE_Ectoderm_2",
    "Primordial_Germ_Cells",
    "Gut")))



## Write the top 10 conserved markers for each cluster to a file
write.table(top10_conserved, file = "output/seurat/srat_top10_conserved_markers_embryo_E7E775.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Visualize the top 10/3 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(embryo.combined.E7E775.sct) <- c(
    "Epiblast_1",
    "Epiblast_2",
    "Unknown_3" ,
    "Nascent_Mesoderm",
    "Somitic_Mesoderm",
    "Unknown_2" ,
    "Caudal_Mesoderm",
    "Pharyngeal_Mesoderm",
    "Unknown_1" ,
    "Paraxial_Mesoderm",
    "Mixed_Mesoderm",
    "Notocord",
    "Haematodenothelial_progenitors",
    "Blood_Progenitor_1",
    "Blood_Progenitor_2",
    "Blood_Progenitor_3",
    "Endoderm",
    "ExE_Endoderm",
    "Mesenchyme",
    "Surface_Ectoderm",
    "ExE_Ectoderm_1",
    "ExE_Ectoderm_2",
    "Primordial_Germ_Cells",
    "Gut")

pdf("output/seurat/DotPlot_SCT_top3_conserved_embryo_E7E775.pdf", width=20, height=5)
DotPlot(embryo.combined.E7E775.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


################## 
## SAVE saveRDS(embryo.combined.E7E775.sct, file = "output/seurat/embryo.combined.E7E775.sct.rds")
## LOAD embryo.combined.E7E775.sct <- readRDS(file = "output/seurat/embryo.combined.E7E775.sct.rds")
##################


# Check some genes
DefaultAssay(embryo.combined.E7E775.sct) <- "SCT" # For vizualization either use SCT or norm RNA

### Check YAP1 hippo bulk-regulated genes:
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775_YAP1_hippo.pdf", width=7, height=30)
# CTGF, CRY61, AREG, MYC, GLI2, VIM, AXL and BIRC5

FeaturePlot(embryo.combined.E7E775.sct, features = c("Tead1", "Tead4", "Ccn2", "Ccn1","Areg","Myc","Gli2","Vim","Axl","Birc5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
dev.off()



### Check YAP1 NODAL bulk-regulated genes:

YAP1_NODAL <- c("Shh", "Dmrt1", "Cited2", 'Dact2', "Tgif2", "Acvr1b", 'Smad3', 'Dand5', 'Smad2','Cer1', 'Nodal', 'Foxh1', 'Dact1', 'Tdgf1', 'Cripto3', 'Acvr1c', 'Cfc1', 'Cfc1b') 

pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775_YAP1_NODAL.pdf", width=10, height=80)
FeaturePlot(embryo.combined.E7E775.sct, features = YAP1_NODAL, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

### Check ANXA NODAL bulk-regulated genes:


pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775_ANXA.pdf", width=10, height=35)
FeaturePlot(embryo.combined.E7E775.sct, features = c("Anxa1","Amotl2","Ccn1","Ccn2","Tagln","Anxa3","Chchd2"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

## cell cycle


pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775_cellCycle.pdf", width=10, height=5)
FeaturePlot(embryo.combined.E7E775.sct,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
dev.off()  


## PRC2-related genes

pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775_PRC2.pdf", width=10, height=22)
FeaturePlot(embryo.combined.E7E775.sct, features = c("Ezh2", "Kdm6b", "Gata6", "Cdx2","T"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

## Hand1
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775_HAND1.pdf", width=10, height=5)
FeaturePlot(embryo.combined.E7E775.sct, features = c("Hand1"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()

## Yap1
pdf("output/seurat/FeaturePlot_SCT_embryo_E7E775_Yap1.pdf", width=10, height=5)
FeaturePlot(embryo.combined.E7E775.sct, features = c("Yap1"), max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()





# GSEA plot
library("fgsea")
## Re-calculate DEGs keeping ALL genes
DefaultAssay(embryo.combined.E7E775.sct) <- "RNA" # For DEGs



embryo.combined.E7E775.sct$celltype.stim <- paste(embryo.combined.E7E775.sct$cluster.annot, embryo.combined.E7E775.sct$condition,
    sep = "-")
Idents(embryo.combined.E7E775.sct) <- "celltype.stim"

Epiblast_1 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Epiblast_1-cYAPKO", ident.2 = "Epiblast_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  

Epiblast_2 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Epiblast_2-cYAPKO", ident.2 = "Epiblast_2-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  


Paraxial_Mesoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Paraxial_Mesoderm-cYAPKO", ident.2 = "Paraxial_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Nascent_Mesoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Nascent_Mesoderm-cYAPKO", ident.2 = "Nascent_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Pharyngeal_Mesoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Pharyngeal_Mesoderm-cYAPKO", ident.2 = "Pharyngeal_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
Mesenchyme <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Mesenchyme-cYAPKO", ident.2 = "Mesenchyme-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    

Mixed_Mesoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Mixed_Mesoderm-cYAPKO", ident.2 = "Mixed_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
Blood_Progenitor_1 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Blood_Progenitor_1-cYAPKO", ident.2 = "Blood_Progenitor_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  

ExE_Ectoderm_1 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "ExE_Ectoderm_1-cYAPKO", ident.2 = "ExE_Ectoderm_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Somitic_Mesoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Somitic_Mesoderm-cYAPKO", ident.2 = "Somitic_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Caudal_Mesoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Caudal_Mesoderm-cYAPKO", ident.2 = "Caudal_Mesoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Haematodenothelial_progenitors <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Haematodenothelial_progenitors-cYAPKO", ident.2 = "Haematodenothelial_progenitors-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Blood_Progenitor_2 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Blood_Progenitor_2-cYAPKO", ident.2 = "Blood_Progenitor_2-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Blood_Progenitor_3 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Blood_Progenitor_3-cYAPKO", ident.2 = "Blood_Progenitor_3-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Surface_Ectoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Surface_Ectoderm-cYAPKO", ident.2 = "Surface_Ectoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Gut <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Gut-cYAPKO", ident.2 = "Gut-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Unknown_1 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Unknown_1-cYAPKO", ident.2 = "Unknown_1-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Unknown_2 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Unknown_2-cYAPKO", ident.2 = "Unknown_2-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Unknown_3 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Unknown_3-cYAPKO", ident.2 = "Unknown_3-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Notocord <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Notocord-cYAPKO", ident.2 = "Notocord-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
Primordial_Germ_Cells <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Primordial_Germ_Cells-cYAPKO", ident.2 = "Primordial_Germ_Cells-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  


Endoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "Endoderm-cYAPKO", ident.2 = "Endoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  

ExE_Ectoderm_2 <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "ExE_Ectoderm_2-cYAPKO", ident.2 = "ExE_Ectoderm_2-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  

ExE_Endoderm <- FindMarkers(embryo.combined.E7E775.sct, ident.1 = "ExE_Endoderm-cYAPKO", ident.2 = "ExE_Endoderm-WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  



### save output
####
write.table(Epiblast_1, file = "output/seurat/Epiblast_1-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Epiblast_2, file = "output/seurat/Epiblast_2-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Unknown_3, file = "output/seurat/Unknown_3-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Nascent_Mesoderm, file = "output/seurat/Nascent_Mesoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Somitic_Mesoderm, file = "output/seurat/Somitic_Mesoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Unknown_2, file = "output/seurat/Unknown_2-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Caudal_Mesoderm, file = "output/seurat/Caudal_Mesoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Pharyngeal_Mesoderm, file = "output/seurat/Pharyngeal_Mesoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Unknown_1, file = "output/seurat/Unknown_1-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Paraxial_Mesoderm, file = "output/seurat/Paraxial_Mesoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Mixed_Mesoderm, file = "output/seurat/Mixed_Mesoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Notocord, file = "output/seurat/Notocord-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Haematodenothelial_progenitors, file = "output/seurat/Haematodenothelial_progenitors-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Blood_Progenitor_1, file = "output/seurat/Blood_Progenitor_1-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Blood_Progenitor_2, file = "output/seurat/Blood_Progenitor_2-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Endoderm, file = "output/seurat/Endoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(ExE_Endoderm, file = "output/seurat/ExE_Endoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Mesenchyme, file = "output/seurat/Mesenchyme-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Surface_Ectoderm, file = "output/seurat/Surface_Ectoderm-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(ExE_Ectoderm_1, file = "output/seurat/ExE_Ectoderm_1-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(ExE_Ectoderm_2, file = "output/seurat/ExE_Ectoderm_2-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Primordial_Germ_Cells, file = "output/seurat/Primordial_Germ_Cells-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Gut, file = "output/seurat/Gut-cYAPKO_response_E7E775_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)

##--> No blood_prog_3 because no cells in WT

#### import all clsuter DEGs output :
cluster_types <- c("Epiblast_1", "Epiblast_2", "Paraxial_Mesoderm", "Nascent_Mesoderm", "Pharyngeal_Mesoderm", "Mesenchyme", "Mixed_Mesoderm", "Blood_Progenitor_1", "ExE_Ectoderm_1", "Somitic_Mesoderm", "Caudal_Mesoderm", "Haematodenothelial_progenitors", "Blood_Progenitor_2", "Surface_Ectoderm", "Gut", "Unknown_1", "Unknown_2", "Unknown_3", "Notocord", "Primordial_Germ_Cells", "Endoderm", "ExE_Ectoderm_2", "ExE_Endoderm")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-cYAPKO_response_E7E775_allGenes.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}


## load list of genes to test
Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

fgsea_sets <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt")$Genes,
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt")$Genes,
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt")$Genes,
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")$Genes
)





## Rank genes based on FC
genes <- Haematodenothelial_progenitors %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(ES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -NES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()




## plot GSEA
pdf("output/Pathway/GSEA_E7E775_ManualGeneLists_Vitamin_A-Haematodenothelial_progenitors.pdf", width=5, height=3)

plotEnrichment(fgsea_sets[["Vitamin_A"]],
               ranks) + labs(title="Vitamin_A-Haematodenothelial_progenitors") +
               theme_bw()
dev.off()



# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c("Epiblast_1", "Epiblast_2", "Paraxial_Mesoderm", "Nascent_Mesoderm", "Pharyngeal_Mesoderm", "Mesenchyme", "Mixed_Mesoderm", "Blood_Progenitor_1", "ExE_Ectoderm_1", "Somitic_Mesoderm", "Caudal_Mesoderm", "Haematodenothelial_progenitors", "Blood_Progenitor_2", "Surface_Ectoderm", "Gut", "Unknown_1", "Unknown_2", "Unknown_3", "Notocord", "Primordial_Germ_Cells", "Endoderm", "ExE_Ectoderm_2", "ExE_Endoderm")

## Initialize an empty list to store the results for each cluster type
all_results <- list()
## Loop over each cluster type
for (cluster in cluster_types) {
  
  # Extract genes for the current cluster
  genes <- get(cluster) %>% 
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC)
  
  ranks <- deframe(genes)
  
  # Run GSEA for the current cluster
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(ES))
  
  # Extract summary table and add cluster column
  fgseaResTidy_summary = fgseaResTidy %>% 
    dplyr::select(pathway, pval, padj, ES, size, NES) %>%
    mutate(cluster = cluster) %>%
    arrange(padj) %>% 
    head()
  
  # Store results in the list
  all_results[[cluster]] <- fgseaResTidy_summary
}
## Combine results from all cluster types into one table
final_results <- bind_rows(all_results, .id = "cluster")


write.table(final_results, file = c("output/Pathway/gsea_output_E7E775_ManualGeneLists.txt"), sep = "\t", quote = FALSE, row.names = FALSE)



# Heatmap all GSEA


pdf("output/Pathway/heatmap_E7E775_gsea_padj.pdf", width=5, height=5)
ggplot(final_results, aes(x=cluster, y=pathway, fill=ES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Enrichment\nScore") +
  geom_text(aes(label=sprintf("%.2f", ES)), 
            color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()



pdf("output/Pathway/heatmap_E7E775_gsea_pval_greyTile.pdf", width=5, height=5)
ggplot(final_results, aes(x=cluster, y=pathway)) + 
  geom_tile(aes(fill = ifelse(pval <= 0.05, NES, NA)), color = "black") +  # Conditional fill based on significance
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Normalized Enrichment\nScore", na.value="grey") +
  # geom_text(aes(label=sprintf("%.2f", ES)), 
  #           color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change between pvalue, qvalue,p.adjust
  #           size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1

dev.off()



# SCPA


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

DefaultAssay(embryo.combined.E7E775.sct) <- "RNA" # Recommended 



# Import manually curated gene list as a list
#### Need to folllow this format
pathways <- msigdbr("Mus musculus", "C2") %>%
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise
pathways$REACTOME_SIGNALING_BY_RETINOIC_ACID$Genes
### Convert table to list of genes
Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

pathways <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt"),
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt"),
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt"),
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")
)

pathways$Hippo$Genes


# Compare pathway activity between condition within  
cell_types <- unique(embryo.combined.E7E775.sct$cluster.annot)
embryo.combined.E7E775.sct_split <- SplitObject(embryo.combined.E7E775.sct, split.by = "condition")

# SCPA comparison
## keep all columns
scpa_out_all <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryo.combined.E7E775.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryo.combined.E7E775.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out_all[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}
## keep only PAthway and qval column (Best for heatmap qval representation)
scpa_out <- list()
for (i in cell_types) {
  
  WT <- seurat_extract(embryo.combined.E7E775.sct_split$WT, 
                            meta1 = "cluster.annot", value_meta1 = i)
  
  cYAPKO <- seurat_extract(embryo.combined.E7E775.sct_split$cYAPKO, 
                          meta1 = "cluster.annot", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out[[i]] <- compare_pathways(list(WT, cYAPKO), pathways, parallel = TRUE, cores = 8) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments
}


saveRDS(scpa_out_all, file = "output/Pathway/scpa_out_eachCellTypes_ManualGeneLists_E7E775.rds")

console_output <- capture.output(print(scpa_out_all))
writeLines(console_output, "output/Pathway/scpa_out_console_ManualGeneLists_E7E775.txt")

## Save output as R object
# saveRDS(scpa_out, file = "output/Pathway/scpa_out_eachCellTypes_E7E775.rds")
# scpa_out_all <- readRDS("output/Pathway/scpa_out_all_eachCellTypes_E7E775.rds")


# plot output

# Code to save output for each cell type comparison
clusters = c(
"Epiblast_1", "Epiblast_2", "Paraxial_Mesoderm", "Nascent_Mesoderm", "Pharyngeal_Mesoderm", "Mesenchyme", "Mixed_Mesoderm", "Blood_Progenitor_1", "ExE_Ectoderm_1", "Somitic_Mesoderm", "Caudal_Mesoderm", "Haematodenothelial_progenitors", "Blood_Progenitor_2", "Surface_Ectoderm", "Gut", "Unknown_1", "Unknown_2", "Unknown_3", "Notocord", "Primordial_Germ_Cells", "Endoderm", "ExE_Ectoderm_2", "ExE_Endoderm"
) # Blood_Prog_3 removed

### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(embryo.combined.E7E775.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  cYAPKO <- seurat_extract(embryo.combined.E7E775.sct,
                           meta1 = "condition", value_meta1 = "cYAPKO",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, cYAPKO),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_E7E775_ManualGeneLists_", cluster, ".txt")
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}



# BALOW CLEAN CODE; dotplot>enrichplot>heatmap>violinplot
# DOT PLOT pepresetnation
## load all the comparison for each cell type (FC qval information)
clusters = c(
"Epiblast_1", "Epiblast_2", "Paraxial_Mesoderm", "Nascent_Mesoderm", "Pharyngeal_Mesoderm", "Mesenchyme", "Mixed_Mesoderm", "Blood_Progenitor_1", "ExE_Ectoderm_1", "Somitic_Mesoderm", "Caudal_Mesoderm", "Haematodenothelial_progenitors", "Blood_Progenitor_2", "Surface_Ectoderm", "Gut", "Unknown_1", "Unknown_2", "Unknown_3", "Notocord", "Primordial_Germ_Cells", "Endoderm", "ExE_Ectoderm_2", "ExE_Endoderm"
) # Blood_Prog_3 removed

## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_E7E775_ManualGeneLists_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster))
# --> import the gene pathways (used table from Conchi `Pathwyas of interest*.xlsx`)
pathways = c("Hippo", "Nodal_TGF", "WNT", "Vitamin_A")
all_data_pathways = as_tibble(all_data) %>% 
  filter(Pathway %in% pathways)

# tidy the df
all_data_pathways_tidy <- all_data_pathways %>%
  filter(qval >= 1.4) 



# REFINE COLOR
custom_color <- function(fc_value){
  ifelse(fc_value >= 5, "strong_blue", 
         ifelse(fc_value > 2, "light_blue",
                ifelse(fc_value <= -5, "strong_red", 
                       ifelse(fc_value < -2, "light_red", "grey"))))
}

# Add a column for this custom color
all_data_pathways_tidy <- all_data_pathways_tidy %>%
  mutate(custom_col = sapply(FC, custom_color))

pdf("output/Pathway/dotplot_E7E775_ManualGeneLists_FCtresh.pdf", width=10, height=4)
ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, color = custom_col), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  scale_color_manual(values = c("strong_blue" = "dodgerblue4",
                                "light_blue" = "lightblue2",
                                "grey" = "grey",
                                "light_red" = "indianred1",
                                "strong_red" = "red3")) +
  theme_bw() +
  labs(size = "q-value", color = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



pdf("output/Pathway/dotplot_E7E775_ManualGeneLists_grey.pdf", width=10, height=4)
ggplot(all_data_pathways_tidy, aes(x = cluster, y = Pathway)) + 
  geom_point(aes(size = qval, fill = "black"), pch=16, alpha=0.7) +   
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  labs(size = "q-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# heatmap XXX to do if needed XXX

Hippo = read_table(file = c("output/Pathway/Manual_geneList_hippo.txt"))
Nodal_TGF = read_table(file = c("output/Pathway/Manual_geneList_Nodal_TGF.txt"))
WNT = read_table(file = c("output/Pathway/Manual_geneList_WNT.txt"))
Vitamin_A = read_table(file = c("output/Pathway/Manual_geneList_Vitamin_A.txt"))

pathways <- list(
  Hippo = read_table(file = "output/Pathway/Manual_geneList_hippo.txt"),
  Nodal_TGF = read_table(file = "output/Pathway/Manual_geneList_Nodal_TGF.txt"),
  WNT = read_table(file = "output/Pathway/Manual_geneList_WNT.txt"),
  Vitamin_A = read_table(file = "output/Pathway/Manual_geneList_Vitamin_A.txt")
)

pathways$Hippo$Genes



genes_of_interest <- pathways$Hippo$Genes
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(embryo.combined.sct@assays$RNA@data)]
### extract WT and cYAPKO gene expression values from RNA assay
WT_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "WT" & embryo.combined.sct$cluster.annot == "Unknown_1"]]
cYAPKO_expression <- embryo.combined.sct@assays$RNA@data[genes_of_interest, colnames(embryo.combined.sct)[embryo.combined.sct$condition == "cYAPKO" & embryo.combined.sct$cluster.annot == "Unknown_1"]]
### mean expression values for each gene
WT_mean <- rowMeans(WT_expression)
cYAPKO_mean <- rowMeans(cYAPKO_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  WT = WT_mean,
  cYAPKO = cYAPKO_mean
) %>% 
pivot_longer(cols = c(WT, cYAPKO), names_to = "Condition", values_to = "Expression")
## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(WT_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)


# if few genes:
pdf("output/Pathway/heatmap_embryo_RNAexpression_E7E775_ManualGeneLists_Nodal_TGF_Blood_Progenitor_1.pdf", width=6, height=2)


library("viridis")
ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()

# if many genesL
pdf("output/Pathway/heatmap_embryo_RNAexpression_E7E775_ManualGeneLists_WNT_Blood_Progenitor_1.pdf", width=3, height=2)

ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") 
dev.off()



# SCPA and GSEA 
gsea_output = final_results %>% 
  rename("pathway" = "Pathway") %>%
  dplyr::select(Pathway, pval, padj, ES, cluster)
scpa_output = all_data_pathways_tidy  %>%
  dplyr::select(Pathway, qval, cluster)

scpa_gsea = gsea_output %>%
  left_join(scpa_output)

## Filter data for GSEA padj < 0.05
filtered_data <- scpa_gsea
## Assign color based on your criteria
filtered_data$is_significant <- filtered_data$padj < 0.05

pdf("output/Pathway/dotplot_scpa_gsea_E7E775.pdf", width=10, height=4)

ggplot(filtered_data, aes(x=cluster, y=Pathway, size=qval)) +
  geom_point(aes(color=ifelse(is_significant, ES, NA), 
                 shape=is_significant), alpha=0.7) +
  scale_color_gradient2(low="blue", high="red", midpoint=0, na.value="black", name = "GSEA Enrichment score") +
  scale_size_continuous(name="SCPA qval") +  # Title for size scale
  scale_shape_manual(values=c(`TRUE`=16, `FALSE`=16)) + # Using circle shape for both
  theme_light() +
  labs(title="SCPA and GSEA", x="Cluster", y="Pathway") +
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


```

--> Pseudotime done in a separate section 


















# Shiny app second time points and data integration

Created with [ShinyCell](https://github.com/SGDDNB/ShinyCell); and follow [shinyapps](https://www.shinyapps.io/) to put it online 

```bash
conda activate scRNAseqV2
# conda install -c anaconda hdf5
```

```R
# installation
## devtools::install_github("SGDDNB/ShinyCell")
## install.packages('rsconnect')


# Packages
library("Seurat")
library("ShinyCell")
library("rsconnect")

# Data import EMBRYO E7 50 dim
embryoE7.combined.sct <- readRDS(file = "output/seurat/embryoE7.combined.sct_50dim.rds")
DefaultAssay(embryoE7.combined.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(embryoE7.combined.sct)

makeShinyApp(embryoE7.combined.sct, scConf, gene.mapping = TRUE,
             shiny.title = "Embryo_E7_50dim_V1",
             shiny.dir = "shinyApp_embryo_E7_50dim_V1/") 

rsconnect::deployApp('shinyApp_embryo_E7_50dim_V1')



# Data import EMBRYO E7 19 dim
embryoE7.combined.sct <- readRDS(file = "output/seurat/embryoE7.combined.sct_19dim.rds")
DefaultAssay(embryoE7.combined.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(embryoE7.combined.sct)

makeShinyApp(embryoE7.combined.sct, scConf, gene.mapping = TRUE,
             shiny.title = "Embryo_E7_19dim_V1",
             shiny.dir = "shinyApp_embryo_E7_19dim_V1/") 

rsconnect::deployApp('shinyApp_embryo_E7_19dim_V1')



# Data import EMBRYO E7 19 dim
embryoE7.combined.sct <- readRDS(file = "output/seurat/embryoE7.combined.sct_19dim_V2.rds")
DefaultAssay(embryoE7.combined.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(embryoE7.combined.sct)

makeShinyApp(embryoE7.combined.sct, scConf, gene.mapping = TRUE,
             shiny.title = "Embryo_E7_19dim_V2",
             shiny.dir = "shinyApp_embryo_E7_19dim_V2/") 

rsconnect::deployApp('shinyApp_embryo_E7_19dim_V2')




# Data import humangastruloid24hr
humangastruloid24hr.combined.sct <- readRDS(file = "output/seurat/humangastruloid24hr.combined.sct_19dim.rds")
DefaultAssay(humangastruloid24hr.combined.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(humangastruloid24hr.combined.sct)

makeShinyApp(humangastruloid24hr.combined.sct, scConf, gene.mapping = TRUE,
             shiny.title = "humangastruloid24hr_25dim_V1",
             shiny.dir = "shinyApp_humangastruloid24hr_25dim_V1/") 

rsconnect::deployApp('shinyApp_humangastruloid24hr_25dim_V1')




# Data import EMBRYO E7 E775 19dim data integration
embryo.combined.E7E775.sct <- readRDS(file = "output/seurat/embryo.combined.E7E775.sct.rds")
DefaultAssay(embryo.combined.E7E775.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(embryo.combined.E7E775.sct)

makeShinyApp(embryo.combined.E7E775.sct, scConf, gene.mapping = TRUE,
             shiny.title = "Embryo_E7E775_19dim_V1",
             shiny.dir = "shinyApp_Embryo_E7E775_19dim_V1/") 

rsconnect::deployApp('shinyApp_Embryo_E7E775_19dim_V1')



```

--> Shiny apps:
- https://roulethomas.shinyapps.io/shinyApp_embryo_E7_50dim_V1/ 
- https://roulethomas.shinyapps.io/shinyApp_embryo_E7_19dim_V1/ 
- https://roulethomas.shinyapps.io/shinyapp_embryo_e7e775_19dim_v1/ 



--> On 20231218 I shared to Conchi the first embryo E7 analysis with DEGs, UMAP marker genes, ppt (`YAP1_scRNAseq_embryo_E7_TR_V1.pptx`)
--> Share to Conchi the Conserved Marker list (`srat_all_conserved_markers_embryo_filter_E7.xlsx`). To avoid confusion, I did some filtering: For each cell type; I only keep log2FC positive (= correspond to gene more highly express in this cell types) and I told her to filter per pvalue which is the max_pvalue. Like this, she will only see the highly express genes in each cluster
----> Conchi provided an updated list (check Dropbox)














