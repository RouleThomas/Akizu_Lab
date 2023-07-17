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





























## Generate sc feature counts for a single library
[Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) and [cell ranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count)


generate single cell feature counts for a single library:

```bash
# fasq have to be .fastq.gz; so lets zip them
gzip input/SRR8734990_* # it took >2hrs

# rename so that they respect naming convention
mv input/SRR8734990_1.fastq.gz input/SRR8734990_S1_L001_R1_001.fastq.gz
mv input/SRR8734990_2.fastq.gz input/SRR8734990_S1_L001_R2_001.fastq.gz
mv input/SRR8734990_3.fastq.gz input/SRR8734990_S1_L001_I1_001.fastq.gz

# ex of command
cellranger count --id=count \
                   --transcriptome=meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/ \
                   --sample=SRR8734990

# run into sbatch
sbatch scripts/cellranger_count.sh # 890789 # ok


```
- *NOTE: it is important to rename our files to respect cell ranger nomenclature*

--> run succesfully!







## RNA preprocessing (doublet and soupX)

Preliminary steps, starting with the **aggregated** `filtered_feature_bc_matrix` from Cellranger 10x:
- doublet detection using `scrublet`
- ambient RNA correction using `soupX`



First step is to detect doublet using [scrublet](https://github.com/swolock/scrublet). Run this in `base` conda env

```bash
srun --mem=500g --pty bash -l

python3 scrublet.py [input_path] [output_path]

python3 scripts/scrublet_doublets.py count/outs/count/filtered_feature_bc_matrix output/doublets/scrublet_50dOrga.tsv
```


--> Successfully assigned doublet
**NOTE: do not name any file `scrublet.smthg` or it will fail!**

Then `conda activate scRNAseq` and go into R and filtered out **RNA contamination and start with SEURAT**:

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

# Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('count/outs/count') 

# Assess % of conta
pdf("output/soupX/autoEstCont_50dOrga.pdf", width=10, height=10)
sc = autoEstCont(sc)
dev.off()

# Generate the corrected matrix
out = adjustCounts(sc)

# Save the matrix
save(out, file = "output/soupX/out_50dOrga.RData")

# SEURAT
srat <- CreateSeuratObject(counts = out, project = "orga50d") # 

# QUALITY CONTROL add mitochondrial and Ribosomal conta
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

# ADD doublet information (scrublet)
doublets <- read.table("output/doublets/scrublet_50dOrga.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat <- AddMetaData(srat,doublets)
srat$Doublet_score <- as.numeric(srat$Doublet_score) # make score as numeric
head(srat[[]])

# Quality check
pdf("output/seurat/VlnPlot_QC_50dOrga.pdf", width=10, height=6)
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_50dOrga.pdf", width=5, height=5)
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_50dOrga.pdf", width=5, height=5)
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAcount_rb_50dOrga.pdf", width=5, height=5)
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")
dev.off()
pdf("output/seurat/FeatureScatter_QC_rb_mt_50dOrga.pdf", width=5, height=5)
FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_50dOrga.pdf", width=5, height=5)
FeatureScatter(srat, feature1 = "nFeature_RNA", feature2 = "Doublet_score")
dev.off()

# After seeing the plot; add QC information in our seurat object
srat[['QC']] <- ifelse(srat@meta.data$Is_doublet == 'True','Doublet','Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC == 'Pass','Low_nFeature',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 20 & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT',paste('High_MT',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
table(srat[['QC']])

# Quality check after QC filtering
pdf("output/seurat/VlnPlot_QCPass_50dOrga.pdf", width=10, height=6)
VlnPlot(subset(srat, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

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