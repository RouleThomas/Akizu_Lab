# Project

For Conchi' project, 3D gastruloid paper (24, and 72hrs); response to reviewer, task from email 20250626 ("we are working on the Stem Cell Reports paper-> back to gastruloids"):

- Compare our untreated 72hr gastruloid with human gastruloid from this [paper](https://www.science.org/doi/10.1126/sciadv.ado1350#sec-1)
    - Data:  ArrayExpress under accession code: E-MTAB-12045



Method used from the paper: "Quality control analysis was performed with Seurat (version 4.0.5) (56) on the raw count matrices, keeping only cells with the number of expressed genes greater or equal to 1000 and the fraction of mitochondrial genes below or equal to 0.20. This resulted in a total of 11,983 cells on D0, 8590 on D2, 7751 on D3, 18,494 on D4, and 23,593 on D8. An additional filtering step was performed with SoupX to take potential environmental RNA contamination into account. After normalization with function NormalizeData in Seurat (56) and correcting for batch effects with the FindIntegrationAnchors function (reduction = “cca”), cluster analysis was performed with the FindClusters function (resolution = 0.1) on the top 20 PCA components computed from the top 2000 genes selected with SelectIntegrationFeatures. The optimal resolution value was chosen according to the output of the clustree function (version 0.4.4) (57). To improve sensitivity to small clusters, cluster analysis at each time point was refined using the concept of entropy of mixing. Briefly, we selected genes frequently detected in small cell neighborhoods by ranking them according to their entropy of mixing computed in local regions defined from the kNN graph. Then, sub-clustering was performed on clusters enriched with low-entropy genes, verified by Fisher’s test. This analysis resulted in one additional cluster on D2 and two additional clusters on D8."

# Data download

From [here](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-12045?query=E-MTAB-12045), I click on ENA ID and arrived [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB55498?show=reads)

--> Let's only download the day 3 data. Samples are annotated 1-12; the **D3 is sample 1 and 7**; corresponding to the two batches. Sample name = Library name. I selected all *FASTQ submitted files* for sample1 and 7; then `Generate Download script`

Refer to [this](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-12045/sdrf) to know sample name. I copy and past to `sample_name.txt`



```bash
# Download script
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096279/SIGAA8_S12_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096315/SIGAD8_S23_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096316/SIGAD8_S23_L002_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096315/SIGAD8_S23_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096316/SIGAD8_S23_L002_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096318/SIGAD8_S24_L002_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096317/SIGAD8_S24_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096276/SIGAA8_S10_L002_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096277/SIGAA8_S11_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096318/SIGAD8_S24_L002_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096312/SIGAD8_S21_L002_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096276/SIGAA8_S10_L002_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096311/SIGAD8_S21_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096312/SIGAD8_S21_L002_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096311/SIGAD8_S21_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096317/SIGAD8_S24_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096282/SIGAA8_S9_L002_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096282/SIGAA8_S9_L002_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096280/SIGAA8_S12_L002_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096275/SIGAA8_S10_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096275/SIGAA8_S10_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096313/SIGAD8_S22_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096278/SIGAA8_S11_L002_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096277/SIGAA8_S11_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096280/SIGAA8_S12_L002_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096278/SIGAA8_S11_L002_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096314/SIGAD8_S22_L002_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096313/SIGAD8_S22_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096281/SIGAA8_S9_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096314/SIGAD8_S22_L002_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096279/SIGAA8_S12_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10096281/SIGAA8_S9_L001_R1_001.fastq.gz
```




# Counting with cellranger count

```bash 
conda activate scRNAseq
which cellranger

# Run counting sample per sample for humangastruloid
sbatch scripts/cellranger_count_D3_hG_Sample1.sh # 46094525 ok
sbatch scripts/cellranger_count_D3_hG_Sample7.sh # 46094636 ok

```




# RNA contamination and doublet detection
- doublet detection using [scrublet](https://github.com/swolock/scrublet) **on the filtered matrix**
- ambient RNA correction using `soupX` in R before generating the Seurat object

```bash
srun --mem=500g --pty bash -l
conda deactivate # base environment needed
python3 scrublet.py [input_path] [output_path]
# Run doublet detection/scrublet sample per sample
python3 scripts/scrublet_doublets.py Sample1/outs/filtered_feature_bc_matrix output/doublets/Sample1.tsv
python3 scripts/scrublet_doublets.py Sample7/outs/filtered_feature_bc_matrix output/doublets/Sample7.tsv

```
Doublet detection score:
- Sample1: 10% expected 13% estimated
- Sample7: 






# Sample1 and Sample7 integration and clustering


Then `conda activate scRNAseq` and go into R and filtered out **RNA contamination and start with SEURAT**.

Tuto seurat [here](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html)

## SoupX RNA cleaning

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
sc = load10X('Sample1/outs') 
sc = load10X('Sample7/outs') 

## Assess % of conta
pdf("output/soupX/autoEstCont_Sample1.pdf", width=10, height=10)
pdf("output/soupX/autoEstCont_Sample7.pdf", width=10, height=10)
sc = autoEstCont(sc)
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/out_Sample1.RData")
save(out, file = "output/soupX/out_Sample7.RData")

```


## Clustering, integration and projection - version1 raw


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
load("output/soupX/out_Sample1.RData")
Sample1 <- CreateSeuratObject(counts = out, project = "Sample1") # 36,601 features across 5,303 samples

load("output/soupX/out_Sample7.RData")
Sample7 <- CreateSeuratObject(counts = out, project = "Sample7") # 36,601 features across 5,235 samples






# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
Sample1[["percent.mt"]] <- PercentageFeatureSet(Sample1, pattern = "^MT-")
Sample1[["percent.rb"]] <- PercentageFeatureSet(Sample1, pattern = "^RP[SL]")

Sample7[["percent.mt"]] <- PercentageFeatureSet(Sample7, pattern = "^MT-")
Sample7[["percent.rb"]] <- PercentageFeatureSet(Sample7, pattern = "^RP[SL]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/Sample1.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
Sample1 <- AddMetaData(Sample1,doublets)
Sample1$Doublet_score <- as.numeric(Sample1$Doublet_score) # make score as numeric
head(Sample1[[]])

doublets <- read.table("output/doublets/Sample7.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
Sample7 <- AddMetaData(Sample7,doublets)
Sample7$Doublet_score <- as.numeric(Sample7$Doublet_score) # make score as numeric
head(Sample7[[]])



## Plot
pdf("output/seurat/VlnPlot_QC_Sample1.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QC_Sample7.pdf", width=10, height=6)

VlnPlot(Sample7, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()



## After seeing the plot; add QC information in our seurat object


### Same filtering as in the paper; addind remnoval of doublet
Sample1[['QC']] <- ifelse(Sample1@meta.data$Is_doublet == 'True','Doublet','Pass')
Sample1[['QC']] <- ifelse(Sample1@meta.data$nFeature_RNA < 1000 & Sample1@meta.data$QC == 'Pass','Low_nFeature',Sample1@meta.data$QC)
Sample1[['QC']] <- ifelse(Sample1@meta.data$nFeature_RNA < 1000 & Sample1@meta.data$QC != 'Pass' & Sample1@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',Sample1@meta.data$QC,sep = ','),Sample1@meta.data$QC)
Sample1[['QC']] <- ifelse(Sample1@meta.data$percent.mt > 20 & Sample1@meta.data$QC == 'Pass','High_MT',Sample1@meta.data$QC)
Sample1[['QC']] <- ifelse(Sample1@meta.data$nFeature_RNA < 1000 & Sample1@meta.data$QC != 'Pass' & Sample1@meta.data$QC != 'High_MT',paste('High_MT',Sample1@meta.data$QC,sep = ','),Sample1@meta.data$QC)
table(Sample1[['QC']]) # 4,131
## 
Sample7[['QC']] <- ifelse(Sample7@meta.data$Is_doublet == 'True','Doublet','Pass')
Sample7[['QC']] <- ifelse(Sample7@meta.data$nFeature_RNA < 1000 & Sample7@meta.data$QC == 'Pass','Low_nFeature',Sample7@meta.data$QC)
Sample7[['QC']] <- ifelse(Sample7@meta.data$nFeature_RNA < 1000 & Sample7@meta.data$QC != 'Pass' & Sample7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',Sample7@meta.data$QC,sep = ','),Sample7@meta.data$QC)
Sample7[['QC']] <- ifelse(Sample7@meta.data$percent.mt > 20 & Sample7@meta.data$QC == 'Pass','High_MT',Sample7@meta.data$QC)
Sample7[['QC']] <- ifelse(Sample7@meta.data$nFeature_RNA < 1000 & Sample7@meta.data$QC != 'Pass' & Sample7@meta.data$QC != 'High_MT',paste('High_MT',Sample7@meta.data$QC,sep = ','),Sample7@meta.data$QC)
table(Sample7[['QC']]) # 3,545
#--> 7676 total cells for me, versus 7751 on D3 (looks good)



## subset my seurat object to only analyze the cells that pass the QC
Sample1 <- subset(Sample1, subset = QC == 'Pass')
Sample7 <- subset(Sample7, subset = QC == 'Pass')
Sample1$replicate <- "Rep1"
Sample7$replicate <- "Rep2"


s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
Sample1 <- NormalizeData(Sample1, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(Sample1)
Sample1 <- ScaleData(Sample1, features = all.genes) # zero-centres and scales it

Sample7 <- NormalizeData(Sample7, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(Sample7)
Sample7 <- ScaleData(Sample7, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
Sample1 <- CellCycleScoring(Sample1, s.features = s.genes, g2m.features = g2m.genes)
table(Sample1[[]]$Phase)
Sample7 <- CellCycleScoring(Sample7, s.features = s.genes, g2m.features = g2m.genes)
table(Sample7[[]]$Phase)

set.seed(42)




# clustering
set.seed(42)

### Previous optimal parameters

# Run SCTransform
## Version OK with 2000 treshold RNA
Sample1 <- SCTransform(Sample1, method = "glmGamPoi", ncells = 4131, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 2000) %>% 
    RunPCA(npcs = 20, verbose = FALSE)
Sample7 <- SCTransform(Sample7, method = "glmGamPoi", ncells = 3545, vars.to.regress = c("nCount_RNA","percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 2000) %>%
    RunPCA(npcs = 20, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(Sample1 = Sample1, Sample7 = Sample7)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 2000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
D3_hG.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
D3_hG <- IntegrateData(anchorset = D3_hG.anchors, normalization.method = "SCT")
# Perform integrated analysis (check active assay is 'integrated')

DefaultAssay(D3_hG) <- "integrated"
D3_hG <- RunPCA(D3_hG, verbose = FALSE, npcs = 20)
D3_hG <- RunUMAP(D3_hG, reduction = "pca", dims = 1:20, verbose = FALSE)
D3_hG <- FindNeighbors(D3_hG, reduction = "pca", k.param = 5, dims = 1:20)
D3_hG <- FindClusters(D3_hG, resolution = 0.1, verbose = FALSE, algorithm = 4)


pdf("output/seurat/UMAP_D3_hG-dim20kparam5res01-splitReplicate.pdf", width=10, height=6)
DimPlot(D3_hG, reduction = "umap", split.by = "replicate", label=TRUE)
dev.off()


pdf("output/seurat/UMAP_D3_hG-dim20kparam5res01.pdf", width=6, height=6)
DimPlot(D3_hG, reduction = "umap", label=TRUE)
dev.off()



#### List from the paper
Ectoderm= c("ZIC1", "EN2")
PS1= c("SOX17", "FOXC1", "TBX6")
Mesendoderm= c("TBXT", "MIXL1", "EOMES")
PS2= c("TBXT", "MIXL1", "WNT3", "EVX1", "NKX1.2")
EarlyMesoderm= c("TBXT", "MESP1", "MESP2", "MSGN1", "FGF8")
IntermediateMesoderm= c("PDGFRA", "OSRl", "PMP22", "TMEM88", "KDR")
ParaxialMesoderm= c("TBX6", "MESP2", "HES7", "MSGN1", "NKX1.2")
SomiticMesoderm= c("RIPPLY1", "MEOXl", "MEOX2", "UNCX", "FGFl8", "HOXC8")
CardiacMesoderm= c("TNNI1", "TNNT2", "HAND1", "BMPS", "GATA4", "GATA6")
NeuralCrest= c("SOX2", "SOX10", "TFAP2A", "TFAP2B", "TFAP2C", "CRABP1")
NeuralTube= c("SOX2", "SOX2", "SOX3", "PAX6", "PAX7", "NEUROG2", "HES5")
Neurons= c("ELAVL3", "TUBB3", "MAP2", "STMN2", "ASCL1", "POU3F1", "POU4F1")
Endoderm= c("SOX17", "FOXA2", "MNX1", "EPCAM", "APOA1")
NMP= c("SOX2", "NKX1.2", "TBXT", "HOXB8", "HOXC9", "WNT3A", "WNT8A")



DefaultAssay(D3_hG) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-Ectoderm.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = Ectoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-PS1.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = PS1, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-Mesendoderm.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = Mesendoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-PS2.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = PS2, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-EarlyMesoderm.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = EarlyMesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-IntermediateMesoderm.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = IntermediateMesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-ParaxialMesoderm.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = ParaxialMesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-SomiticMesoderm.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = SomiticMesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-CardiacMesoderm.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = CardiacMesoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-NeuralCrest.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = NeuralCrest, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-NeuralTube.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = NeuralTube, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-Neurons.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = Neurons, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-Endoderm.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = Endoderm, max.cutoff = 3, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_D3_hG-dim20-NMP.pdf", width=10, height=10)
FeaturePlot(D3_hG, features = NMP, max.cutoff = 3, cols = c("grey", "red"))
dev.off()





# Rename cluster
Cluster1= PS (TBXT, MIXL1, WNT3, EVX1)
Cluster2= ParaxialMesoderm (TBX6, MESP2, MSGN1)
Cluster3= IntermediateMesoderm (PMP22, PDGFRA, TMEM88, KDR)
Cluster4= CardiacMesoderm (TNNI1, TNNT2, HAND1, GATA4, GATA6)
Cluster5= Ectoderm (ZIC1, EN2)
Cluster6= Endoderm (FOXA2, MNX1, EPCAM)


new.cluster.ids <- c(
  "PS",
  "ParaxialMesoderm",
  "IntermediateMesoderm",
  "CardiacMesoderm",
  "Ectoderm",
  "Endoderm"
)

names(new.cluster.ids) <- levels(D3_hG)
D3_hG <- RenameIdents(D3_hG, new.cluster.ids)

D3_hG$cluster.annot <- Idents(D3_hG) # create a new slot in my seurat object


pdf("output/seurat/UMAP_UMAP_D3_hG-dim20kparam5res01_label.pdf", width=6, height=5)
DimPlot(D3_hG, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 5)
dev.off()



# SAve 
## saveRDS(D3_hG, file = "output/seurat/D3_hG-dim20kparam5res01_label.rds")
D3_hG = readRDS(file = "output/seurat/D3_hG-dim20kparam5res01_label.rds")





# All in dotplot
DefaultAssay(D3_hG) <- "SCT"

#### Gene marker

Cluster1= PS (TBXT, MIXL1, WNT3, EVX1)
Cluster2= ParaxialMesoderm (TBX6, MESP2, MSGN1)
Cluster3= IntermediateMesoderm (PMP22, PDGFRA, TMEM88, KDR)
Cluster4= CardiacMesoderm (TNNI1, TNNT2, HAND1, GATA4, GATA6)
Cluster5= Ectoderm (ZIC1, EN2)
Cluster6= Endoderm (FOXA2, MNX1, EPCAM)



all_markers <- c(
  "TBXT", "MIXL1", "WNT3", "EVX1",
  "TBX6", "MESP2", "MSGN1",
  "PMP22", "PDGFRA", "TMEM88", "KDR",
  "TNNI1", "TNNT2", "HAND1", "GATA4", "GATA6",
  "ZIC1", "EN2",
  "FOXA2", "MNX1", "EPCAM"
)


levels(D3_hG) <- c(
  "PS",
  "ParaxialMesoderm",
  "IntermediateMesoderm",
  "CardiacMesoderm",
  "Ectoderm",
  "Endoderm"
)


pdf("output/seurat/DotPlot_D3_hG-dim20kparam5res01.pdf", width=9, height=2.5)
DotPlot(D3_hG, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()


######################################################################
##### scRNAseq projection ###################################
######################################################################


# Import our data

humangastruloid_dim25kparam15res02 <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid.combined.sct_V2-dim25kparam15res02.rds")
humangastruloid_dim25kparam50res07 <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid.combined.sct_V2-dim25kparam50res07.rds")




# Do scRNAseq projection (reference humanGasutrla and query humangastruloid_dim25kparam15res02) - V1 
# Set assay
DefaultAssay(D3_hG) <- "RNA"
DefaultAssay(humangastruloid_dim25kparam15res02) <- "RNA"


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(D3_hG),
  rownames(humangastruloid_dim25kparam15res02)
)

anchors <- FindTransferAnchors(
  reference = D3_hG,
  query = humangastruloid_dim25kparam15res02,
  normalization.method = "LogNormalize",
  reduction = "pcaproject",
  dims = 1:20,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)



# Transfer labels (e.g. cluster_id from reference)
predictions <- TransferData(
  anchorset = anchors,
  refdata = D3_hG$cluster.annot,  # or sub_cluster if you prefer
  dims = 1:20,
  k.weight = 100
)

# Add predicted labels to query object
humangastruloid_dim25kparam15res02 <- AddMetaData(humangastruloid_dim25kparam15res02, metadata = predictions)

pdf("output/seurat/UMAP_D3_hG-Projected_humangastruloid_dim25kparam15res02.pdf", width=10, height=6)
DimPlot(humangastruloid_dim25kparam15res02, group.by = "predicted.id", label = TRUE)
dev.off()
#--> This is the opposite as done in the paper: we annotate our cells using the other cell type annotation


# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
D3_hG <- RunUMAP(
  D3_hG,
  dims = 1:20,
  reduction = "pca",
  return.model = TRUE  
)



anchors <- FindTransferAnchors(
  reference = D3_hG,
  query = humangastruloid_dim25kparam15res02,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:20,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = D3_hG$cluster.annot,
  dims = 1:20,
  k.weight = 100
)

humangastruloid_dim25kparam15res02 <- AddMetaData(humangastruloid_dim25kparam15res02, metadata = predictions)


# Step 3: Project query onto reference UMAP
humangastruloid_dim25kparam15res02 <- MapQuery(
  anchorset = anchors,
  reference = D3_hG,
  query = humangastruloid_dim25kparam15res02,
  refdata = list(cluster_id = "cluster.annot"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
# Panel 1: D3_hG
p1 <- DimPlot(D3_hG, reduction = "umap", group.by = "cluster.annot", label = TRUE, pt.size = 1) +
  ggtitle("D3_hG")

# Panel 2: Projected gastruloid
p2 <- DimPlot(humangastruloid_dim25kparam15res02, reduction = "ref.umap", group.by = "cluster.annot", label = TRUE, pt.size = 1) +
  ggtitle("Human gastruloid 72hr")

# Panel 3: True overlay

# Plot reference in gray
D3_hG$dummy_group <- "Reference"
p_ref <- DimPlot(D3_hG, reduction = "umap", group.by = "dummy_group", cols = "lightgray", pt.size = 1) +
  NoLegend()

# Plot query separately to extract point positions + colors
p_query <- DimPlot(humangastruloid_dim25kparam15res02, reduction = "ref.umap", group.by = "cluster.annot", pt.size = 1) +
  NoAxes() + NoLegend()

# Extract layers and combine manually
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(data = query_layer, aes(x = x, y = y), color = query_layer$colour, size = 1) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humangastruloid_dim25kparam15res02$prediction.score.max,
  group = "72hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humangastruloid_dim25kparam15res02$prediction.score.max,
  predicted_id = humangastruloid_dim25kparam15res02$predicted.id
)
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


```


--> Poor overlap with paraxial mesoderm with this dataset

## Projection - version2 clean


Let's load the .rds integrated object from version1 but re-do the projection





```bash
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
library("ggpubr")

set.seed(42)


# import D3_hG from version1 raw

D3_hG = readRDS(file = "output/seurat/D3_hG-dim20kparam5res01_label.rds")

DefaultAssay(D3_hG) <- "RNA"


# Import our data

humangastruloid_dim25kparam15res02 <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid.combined.sct_V2-dim25kparam15res02.rds")
humangastruloid_dim25kparam50res07 <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid.combined.sct_V2-dim25kparam50res07.rds")


DefaultAssay(humangastruloid_dim25kparam15res02) <- "RNA"


# Do scRNAseq projection (reference D3_hG and query humangastruloid_dim25kparam15res02) - V1 

# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
humangastruloid_dim25kparam15res02 <- RunUMAP(
  humangastruloid_dim25kparam15res02,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


# Find anchors using PCA projection
shared_features <- intersect(
  rownames(D3_hG),
  rownames(humangastruloid_dim25kparam15res02)
)



anchors <- FindTransferAnchors(
  reference = humangastruloid_dim25kparam15res02 ,
  query = D3_hG,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = humangastruloid_dim25kparam15res02$cluster.annot,
  dims = 1:25,
  k.weight = 100
)

D3_hG <- AddMetaData(D3_hG, metadata = predictions)


# Step 3: Project query onto reference UMAP
D3_hG <- MapQuery(
  anchorset = anchors,
  reference = humangastruloid_dim25kparam15res02,
  query = D3_hG,
  refdata = list(cluster_id = "cluster.annot"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
# Step 1: Create a color palette for all cluster levels
# Get all unique cluster names from both reference and query
all_clusters <- union(
  unique(humangastruloid_dim25kparam15res02$cluster.annot),
  unique(D3_hG$predicted.id)
)

# Assign colors (adjust palette as needed or use scales::hue_pal())
cluster_colors <- setNames(
  scales::hue_pal()(length(all_clusters)),
  sort(all_clusters)
)

# Panel 1: Human gastruloid
p1 <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "umap",
  group.by = "cluster.annot",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 2: Projected gastrula
p2 <- DimPlot(
  D3_hG,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("D3_hG")

# Panel 3: Overlay
# Reference in gray
humangastruloid_dim25kparam15res02$dummy_group <- "Reference"
p_ref <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Query plotted separately with correct colors
p_query <- DimPlot(
  D3_hG,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Extract and overlay
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,
    size = 1,
    show.legend = FALSE
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay-order1-version2.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = D3_hG$prediction.score.max,
  group = "D3_hG"
)

# Plot
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay-order1-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = D3_hG$prediction.score.max,
  predicted_id = D3_hG$predicted.id
)
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay-order1-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()








##############################################################
# The other way ##############################################
##############################################################
DefaultAssay(humangastruloid_dim25kparam15res02) <- "RNA"


# Step 1: Find anchors (you may have already done this)
D3_hG <- RunUMAP(
  D3_hG,
  dims = 1:20,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = D3_hG,
  query = humangastruloid_dim25kparam15res02,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:20,
  features = shared_features,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = D3_hG$cluster.annot,
  dims = 1:20,
  k.weight = 100
)

humangastruloid_dim25kparam15res02 <- AddMetaData(humangastruloid_dim25kparam15res02, metadata = predictions)


# Step 3: Project query onto reference UMAP
humangastruloid_dim25kparam15res02 <- MapQuery(
  anchorset = anchors,
  reference = D3_hG,
  query = humangastruloid_dim25kparam15res02,
  refdata = list(cluster_id = "cluster.annot"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
all_clusters <- union(
  unique(D3_hG$cluster.annot),
  unique(humangastruloid_dim25kparam15res02$predicted.id)
)

# Step 2: Assign consistent colors
cluster_colors <- setNames(scales::hue_pal()(length(all_clusters)), sort(all_clusters))
# Panel 1: Human gastrula
p1 <- DimPlot(
  D3_hG,
  reduction = "umap",
  group.by = "cluster.annot",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("D3_hG")

# Panel 2: Projected gastruloid
p2 <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "ref.umap",
  group.by = "predicted.id",
  label = TRUE,
  pt.size = 1,
  cols = cluster_colors
) + ggtitle("Human gastruloid 72hr")

# Panel 3: Overlay

# Plot reference (D3_hG) in gray
D3_hG$dummy_group <- "Reference"
p_ref <- DimPlot(
  D3_hG,
  reduction = "umap",
  group.by = "dummy_group",
  cols = "lightgray",
  pt.size = 1
) + NoLegend()

# Plot projected query separately with colors
p_query <- DimPlot(
  humangastruloid_dim25kparam15res02,
  reduction = "ref.umap",
  group.by = "predicted.id",
  pt.size = 1,
  cols = cluster_colors
) + NoAxes() + NoLegend()

# Overlay: Extract and draw query points on top of reference
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(
    data = query_layer,
    aes(x = x, y = y),
    color = query_layer$colour,  # already hex
    size = 1
  ) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4: Export
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay-version2.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()





# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humangastruloid_dim25kparam15res02$prediction.score.max,
  group = "72hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humangastruloid_dim25kparam15res02$prediction.score.max,
  predicted_id = humangastruloid_dim25kparam15res02$predicted.id
)
pdf("output/seurat/UMAP_D3_hG-reference_query_overlay-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


```





















XXX BELOW NOT MOD




# Analyzis with the processed data


Let's try to directly used the processed data from the authors


```bash
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
library("ggpubr")

set.seed(42)


# import humanGastrula .rds files

raw_reads <- readRDS("input/raw_reads.rds")
umap_info <- readRDS("input/umap.rds")



# Set cell names as rownames
rownames(raw_reads) <- umap_info$cell_name

# Convert to seurat
humanGastrula <- CreateSeuratObject(counts = t(raw_reads))  # transpose since genes are in columns

# Build UMAP matrix
umap_mat <- as.matrix(umap_info[, c("X", "X1")])
rownames(umap_mat) <- umap_info$cell_name
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")

# Attach UMAP to Seurat
humanGastrula[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAP_", assay = "RNA")

# Make sure rownames are cell names
rownames(umap_info) <- umap_info$cell_name

# Add desired metadata columns to the Seurat object
humanGastrula$cluster_id <- umap_info$cluster_id
humanGastrula$sub_cluster <- umap_info$sub_cluster


head(humanGastrula@meta.data)


# Re-perform clustering

DefaultAssay(humanGastrula) <- "RNA"

## NORMALIZE AND SCALE DATA 
humanGastrula <- NormalizeData(humanGastrula, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(humanGastrula)
humanGastrula <- ScaleData(humanGastrula, features = all.genes) # zero-centres and scales it


## Find variable features
humanGastrula = FindVariableFeatures(humanGastrula, nfeatures = 3000)

humanGastrula <- RunPCA(humanGastrula, verbose = FALSE, npcs = 25)
humanGastrula <- RunUMAP(humanGastrula, reduction = "pca", dims = 1:25, verbose = FALSE)
humanGastrula <- FindNeighbors(humanGastrula, reduction = "pca", k.param = 15, dims = 1:25)
humanGastrula <- FindClusters(humanGastrula, resolution = 0.3, verbose = FALSE, algorithm = 4)



# Plot to confirm it work

pdf("output/seurat/UMAP_humanGastrula.pdf", width=10, height=6)
DimPlot(humanGastrula, reduction = "umap", label=TRUE, group.by = "cluster_id")
dev.off()






# Import our data

humangastruloid_dim25kparam15res02 <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid.combined.sct_V2-dim25kparam15res02.rds")
humangastruloid_dim25kparam50res07 <- readRDS(file = "../003__YAP1/output/seurat/humangastruloid.combined.sct_V2-dim25kparam50res07.rds")




# Do scRNAseq projection (reference humanGasutrla and query humangastruloid_dim25kparam15res02) - V1 
# Set assay
DefaultAssay(humanGastrula) <- "RNA"
DefaultAssay(humangastruloid_dim25kparam15res02) <- "RNA"


# Find anchors using PCA projection
anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = humangastruloid_dim25kparam15res02,
  normalization.method = "LogNormalize", # use "SCT" if both objects are SCTransformed; use "LogNormalize" otherwise
  reduction = "pcaproject",
  dims = 1:25,
  k.anchor = 100,
  k.filter = 500
)

# Transfer labels (e.g. cluster_id from reference)
predictions <- TransferData(
  anchorset = anchors,
  refdata = humanGastrula$cluster_id,  # or sub_cluster if you prefer
  dims = 1:25,
  k.weight = 100
)

# Add predicted labels to query object
humangastruloid_dim25kparam15res02 <- AddMetaData(humangastruloid_dim25kparam15res02, metadata = predictions)

pdf("output/seurat/UMAP_humanGastrula-Projected_humangastruloid_dim25kparam15res02.pdf", width=10, height=6)
DimPlot(humangastruloid_dim25kparam15res02, group.by = "predicted.id", label = TRUE)
dev.off()

#--> V1 is the opposite as done in the paper: we annotate our cells using the other cell type annotation




# V2 - like in the paper

# Step 1: Find anchors (you may have already done this)
humanGastrula <- RunUMAP(
  humanGastrula,
  dims = 1:25,
  reduction = "pca",
  return.model = TRUE  
)


anchors <- FindTransferAnchors(
  reference = humanGastrula,
  query = humangastruloid_dim25kparam15res02,
  normalization.method = "LogNormalize",  # or "LogNormalize"
  reference.reduction = "pca",  
  reduction = "pcaproject",
  dims = 1:25,
  k.anchor = 100,
  k.filter = 500
)

# Step 2: Transfer labels from reference to query
predictions <- TransferData(
  anchorset = anchors,
  refdata = humanGastrula$cluster_id,
  dims = 1:25,
  k.weight = 100
)

humangastruloid_dim25kparam15res02 <- AddMetaData(humangastruloid_dim25kparam15res02, metadata = predictions)


# Step 3: Project query onto reference UMAP
humangastruloid_dim25kparam15res02 <- MapQuery(
  anchorset = anchors,
  reference = humanGastrula,
  query = humangastruloid_dim25kparam15res02,
  refdata = list(cluster_id = "cluster_id"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# Step 4: Generate UMAP plots
# Panel 1: Human gastrula
p1 <- DimPlot(humanGastrula, reduction = "umap", group.by = "cluster_id", label = TRUE, pt.size = 1) +
  ggtitle("Human gastrula")

# Panel 2: Projected gastruloid
p2 <- DimPlot(humangastruloid_dim25kparam15res02, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, pt.size = 1) +
  ggtitle("Human gastruloid 72hr")

# Panel 3: True overlay

# Plot reference in gray
humanGastrula$dummy_group <- "Reference"
p_ref <- DimPlot(humanGastrula, reduction = "umap", group.by = "dummy_group", cols = "lightgray", pt.size = 1) +
  NoLegend()

# Plot query separately to extract point positions + colors
p_query <- DimPlot(humangastruloid_dim25kparam15res02, reduction = "ref.umap", group.by = "predicted.id", pt.size = 1) +
  NoAxes() + NoLegend()

# Extract layers and combine manually
g_ref <- p_ref[[1]]
query_layer <- ggplot_build(p_query[[1]])$data[[1]]

g_overlay <- g_ref +
  geom_point(data = query_layer, aes(x = x, y = y), color = query_layer$colour, size = 1) +
  ggtitle("Overlay") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5: Export
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay.pdf", width = 30, height = 7)
(p1 | p2 | g_overlay)
dev.off()



# Prediciton score
# Extract scores and metadata
# Create a dummy timepoint to allow boxplot
df <- data.frame(
  prediction_score = humangastruloid_dim25kparam15res02$prediction.score.max,
  group = "72hr Gastruloid"
)

# Plot
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score.pdf", width = 5, height = 7)
ggplot(df, aes(x = group, y = prediction_score)) +
  geom_boxplot(fill = "white", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


df2 <- data.frame(
  prediction_score = humangastruloid_dim25kparam15res02$prediction.score.max,
  predicted_id = humangastruloid_dim25kparam15res02$predicted.id
)
pdf("output/seurat/UMAP_humanGastrula-reference_query_overlay-prediction_score_predicted_id.pdf", width = 5, height = 4)
ggplot(df2, aes(x = predicted_id, y = prediction_score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  coord_flip() +
  labs(title = "Prediction score", x = NULL, y = NULL)
dev.off()


```

--> Pretty well recapitulated!!! I generated a plot that say the prediction score for each cluster, it say how well each of the cluster is recapitulated in our data (>0.8 can be considered very high).












# Data analyzis

XXX MAYBE NOT NEEDED AS I TRY THE PROCESSED DATA


Follow the same method as we used in Conchi paper, BUT use their QC filetering parameters; we need clustering as close to what they publish.








