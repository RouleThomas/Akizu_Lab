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

## human gastruloid_analysis V2_integration only

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
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_DASATINIB72hr_cellCycle_V2.pdf", width=10, height=10)
VlnPlot(humangastruloid.combined.sct,features = c("S.Score","G2M.Score"), split.by = "condition") & 
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

### Loop through 100 iterations
for (i in 1:10) {
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


### V2

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
pdf("output/seurat/VlnPlot_SCT_UNTREATED72hr_DASATINIB72hr_cellCycle_V2.pdf", width=10, height=10)
VlnPlot(humangastruloid.combined.sct,features = c("S.Score","G2M.Score"), split.by = "condition") & 
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

### Loop through 100 iterations
for (i in 1:10) {
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


### V2

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

YAP1_NODAL <- c("SHH", "DMRT1", "CITED2", 'DACT2', "TGIF2", "ACVR1B", 'SMAD3', 'DAND5', 'SMAD2','CER1', 'NODAL', 'FOXH1', 'DACT1', 'TDGF1', 'TDGF1P3', 'ACVR1C', 'CFC1', 'CFC1B') 

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_NODAL_V2.pdf", width=10, height=80)
FeaturePlot(humangastruloid.combined.sct, features = YAP1_NODAL, max.cutoff = 3, cols = c("grey", "red"), split.by = "condition")
dev.off()


### Check YAP1 hippo bulk-regulated genes:
DefaultAssay(humangastruloid.combined.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_UNTREATED72hr_DASATINIB72hr_YAP1_hippo_V2.pdf", width=7, height=25)

FeaturePlot(humangastruloid.combined.sct, features = c("TEAD1", "TEAD4", "CCN2", "CCN1","AREG","MYC","GLI2","VIM","AXL","BIRC5"), split.by = "condition", max.cutoff = 3, cols = c("grey", "red"))
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

```