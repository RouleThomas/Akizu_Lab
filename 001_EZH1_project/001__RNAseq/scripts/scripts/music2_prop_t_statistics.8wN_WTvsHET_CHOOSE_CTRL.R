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
library("gprofiler2") 
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")
library("MuSiC2")
library("MuSiC")
library("biomaRt")


print("packages loaded")
# CHECK CONTROL FILE


## import seurat object
CHOOSE_CTRL_annot_srt <- readRDS(file = "output/MuSiC2/CHOOSE_CTRL_annot_srt.rds")
DefaultAssay(CHOOSE_CTRL_annot_srt) <- "RNA" # 



## isolate the cells of interest; control one
cells_to_keep <- WhichCells(CHOOSE_CTRL_annot_srt, expression = celltype_cl_coarse2 != "NA")

CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2 <- subset(CHOOSE_CTRL_annot_srt, cells = cells_to_keep)




## Transform scRNAseq data in ExpressionSet Class
CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.sceset = ExpressionSet(assayData = as.matrix(GetAssayData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2)), phenoData =  new("AnnotatedDataFrame",CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2@meta.data))


## Transform scRNAseq data in SingleCellExperiment (for music2 t statistics)

CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment <- as.SingleCellExperiment(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2, assay = "RNA")









# Import raw RNAseq read counts
# code for ensembl to genesymbol conv
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
genes <- X8wN_WT_R1$Geneid
#### Get the mapping from Ensembl ID to gene symbol
genes_mapped  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = genes,
                      mart = ensembl)
####

### import featureCounts output
#### 8wN WT R1
X8wN_WT_R1 <- read.delim("output/featurecounts_hg38/8wN_WT_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R1$Geneid <- gsub("\\..*", "", X8wN_WT_R1$Geneid)
X8wN_WT_R1_geneSymbol <- merge(X8wN_WT_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam')
## Calculate median as gene dupplicated with the conversino!
X8wN_WT_R1_count_summary <- X8wN_WT_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R1_count_matrix <- as.matrix(X8wN_WT_R1_count_summary$median_count)
rownames(X8wN_WT_R1_count_matrix) <- X8wN_WT_R1_count_summary$external_gene_name
colnames(X8wN_WT_R1_count_matrix) <- "8wN_WT_R1"



#### 8wN WT R2
X8wN_WT_R2 <- read.delim("output/featurecounts_hg38/8wN_WT_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R2$Geneid <- gsub("\\..*", "", X8wN_WT_R2$Geneid)
X8wN_WT_R2_geneSymbol <- merge(X8wN_WT_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R2_count_summary <- X8wN_WT_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R2_count_matrix <- as.matrix(X8wN_WT_R2_count_summary$median_count)
rownames(X8wN_WT_R2_count_matrix) <- X8wN_WT_R2_count_summary$external_gene_name
colnames(X8wN_WT_R2_count_matrix) <- "8wN_WT_R2"


#### 8wN WT R3
X8wN_WT_R3 <- read.delim("output/featurecounts_hg38/8wN_WT_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R3$Geneid <- gsub("\\..*", "", X8wN_WT_R3$Geneid)
X8wN_WT_R3_geneSymbol <- merge(X8wN_WT_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R3_count_summary <- X8wN_WT_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R3_count_matrix <- as.matrix(X8wN_WT_R3_count_summary$median_count)
rownames(X8wN_WT_R3_count_matrix) <- X8wN_WT_R3_count_summary$external_gene_name
colnames(X8wN_WT_R3_count_matrix) <- "8wN_WT_R3"


#### 8wN WT R4
X8wN_WT_R4 <- read.delim("output/featurecounts_hg38/8wN_WT_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_WT_R4$Geneid <- gsub("\\..*", "", X8wN_WT_R4$Geneid)
X8wN_WT_R4_geneSymbol <- merge(X8wN_WT_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_WT_R4_count_summary <- X8wN_WT_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_WT_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_WT_R4_count_matrix <- as.matrix(X8wN_WT_R4_count_summary$median_count)
rownames(X8wN_WT_R4_count_matrix) <- X8wN_WT_R4_count_summary$external_gene_name
colnames(X8wN_WT_R4_count_matrix) <- "8wN_WT_R4"


#### 8wN HET R1
X8wN_HET_R1 <- read.delim("output/featurecounts_hg38/8wN_HET_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R1$Geneid <- gsub("\\..*", "", X8wN_HET_R1$Geneid)
X8wN_HET_R1_geneSymbol <- merge(X8wN_HET_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R1_count_summary <- X8wN_HET_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R1_count_matrix <- as.matrix(X8wN_HET_R1_count_summary$median_count)
rownames(X8wN_HET_R1_count_matrix) <- X8wN_HET_R1_count_summary$external_gene_name
colnames(X8wN_HET_R1_count_matrix) <- "8wN_HET_R1"



#### 8wN HET R2
X8wN_HET_R2 <- read.delim("output/featurecounts_hg38/8wN_HET_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R2$Geneid <- gsub("\\..*", "", X8wN_HET_R2$Geneid)
X8wN_HET_R2_geneSymbol <- merge(X8wN_HET_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R2_count_summary <- X8wN_HET_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R2_count_matrix <- as.matrix(X8wN_HET_R2_count_summary$median_count)
rownames(X8wN_HET_R2_count_matrix) <- X8wN_HET_R2_count_summary$external_gene_name
colnames(X8wN_HET_R2_count_matrix) <- "8wN_HET_R2"


#### 8wN HET R3
X8wN_HET_R3 <- read.delim("output/featurecounts_hg38/8wN_HET_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R3$Geneid <- gsub("\\..*", "", X8wN_HET_R3$Geneid)
X8wN_HET_R3_geneSymbol <- merge(X8wN_HET_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R3_count_summary <- X8wN_HET_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R3_count_matrix <- as.matrix(X8wN_HET_R3_count_summary$median_count)
rownames(X8wN_HET_R3_count_matrix) <- X8wN_HET_R3_count_summary$external_gene_name
colnames(X8wN_HET_R3_count_matrix) <- "8wN_HET_R3"


#### 8wN HET R4
X8wN_HET_R4 <- read.delim("output/featurecounts_hg38/8wN_HET_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_HET_R4$Geneid <- gsub("\\..*", "", X8wN_HET_R4$Geneid)
X8wN_HET_R4_geneSymbol <- merge(X8wN_HET_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_HET_R4_count_summary <- X8wN_HET_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_HET_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_HET_R4_count_matrix <- as.matrix(X8wN_HET_R4_count_summary$median_count)
rownames(X8wN_HET_R4_count_matrix) <- X8wN_HET_R4_count_summary$external_gene_name
colnames(X8wN_HET_R4_count_matrix) <- "8wN_HET_R4"


#### 8wN KO R1
X8wN_KO_R1 <- read.delim("output/featurecounts_hg38/8wN_KO_R1.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R1$Geneid <- gsub("\\..*", "", X8wN_KO_R1$Geneid)
X8wN_KO_R1_geneSymbol <- merge(X8wN_KO_R1, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R1_count_summary <- X8wN_KO_R1_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R1_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R1_count_matrix <- as.matrix(X8wN_KO_R1_count_summary$median_count)
rownames(X8wN_KO_R1_count_matrix) <- X8wN_KO_R1_count_summary$external_gene_name
colnames(X8wN_KO_R1_count_matrix) <- "8wN_KO_R1"



#### 8wN KO R2
X8wN_KO_R2 <- read.delim("output/featurecounts_hg38/8wN_KO_R2.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R2$Geneid <- gsub("\\..*", "", X8wN_KO_R2$Geneid)
X8wN_KO_R2_geneSymbol <- merge(X8wN_KO_R2, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R2_count_summary <- X8wN_KO_R2_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R2_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R2_count_matrix <- as.matrix(X8wN_KO_R2_count_summary$median_count)
rownames(X8wN_KO_R2_count_summary) <- X8wN_KO_R2_count_summary$external_gene_name
colnames(X8wN_KO_R2_count_summary) <- "8wN_KO_R2"



#### 8wN KO R3
X8wN_KO_R3 <- read.delim("output/featurecounts_hg38/8wN_KO_R3.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R3$Geneid <- gsub("\\..*", "", X8wN_KO_R3$Geneid)
X8wN_KO_R3_geneSymbol <- merge(X8wN_KO_R3, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R3_count_summary <- X8wN_KO_R3_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R3_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R3_count_matrix <- as.matrix(X8wN_KO_R3_count_summary$median_count)
rownames(X8wN_KO_R3_count_summary) <- X8wN_KO_R3_count_summary$external_gene_name
colnames(X8wN_KO_R3_count_summary) <- "8wN_KO_R3"



#### 8wN KO R4
X8wN_KO_R4 <- read.delim("output/featurecounts_hg38/8wN_KO_R4.txt", header=TRUE, stringsAsFactors=FALSE, skip = 1) 
### convert enesembl ID to gene Symbol
X8wN_KO_R4$Geneid <- gsub("\\..*", "", X8wN_KO_R4$Geneid)
X8wN_KO_R4_geneSymbol <- merge(X8wN_KO_R4, genes_mapped, by.x = 'Geneid', by.y = 'ensembl_gene_id', all.x = TRUE) %>% 
  drop_na() %>%
  dplyr::select(external_gene_name, 'output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam')
### Calculate median as gene dupplicated with the conversino!
X8wN_KO_R4_count_summary <- X8wN_KO_R4_geneSymbol %>%
  group_by(external_gene_name) %>%
  summarize(median_count = median(`output.STAR_hg38.8wN_KO_R4_Aligned.sortedByCoord.out.bam`)) %>%
  ungroup() 
X8wN_KO_R4_count_matrix <- as.matrix(X8wN_KO_R4_count_summary$median_count)
rownames(X8wN_KO_R4_count_summary) <- X8wN_KO_R4_count_summary$external_gene_name
colnames(X8wN_KO_R4_count_summary) <- "8wN_KO_R4"



print("bulk RNAseq data imported")


## WT versus HET comparison
all_counts <- cbind(X8wN_WT_R1_count_matrix, X8wN_WT_R2_count_matrix, X8wN_WT_R3_count_matrix, X8wN_WT_R4_count_matrix, X8wN_HET_R1_count_matrix, X8wN_HET_R2_count_matrix, X8wN_HET_R3_count_matrix, X8wN_HET_R4_count_matrix)
### Create the phenotype data frame
pheno_data <- data.frame(sampleID = colnames(all_counts),
                         group = c("WT", "WT","WT", "WT", "HET", "HET", "HET", "HET"),
                         row.names = colnames(all_counts))
### Convert pheno_data to AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data = pheno_data)
# Make sure that your combined matrix is indeed a matrix
all_counts_matrix <- as.matrix(all_counts)
### Now create the ExpressionSet object
WT_HET_8wN.bulkeset <- ExpressionSet(assayData = all_counts_matrix,
                      phenoData = phenoData)
table(WT_HET_8wN.bulkeset$group)

bulk.control.mtx = exprs(WT_HET_8wN.bulkeset)[, WT_HET_8wN.bulkeset$group == 'WT']
bulk.case.mtx = exprs(WT_HET_8wN.bulkeset)[, WT_HET_8wN.bulkeset$group == 'HET']

# TO AVOID SUBSCRIPTY OUT OF BOUND ERROR:
### subset bulk features that are only present in scRNAseq assay
###### Assuming 'featureNames' can be used to retrieve the features from both ExpressionSets
common_features <- intersect(featureNames(WT_HET_8wN.bulkeset), featureNames(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.sceset))
# Subset the WT_HET_8wN.bulkeset to only include common features
WT_HET_8wN.bulkeset_subset <- WT_HET_8wN.bulkeset[common_features, ]


bulk.control.mtx = exprs(WT_HET_8wN.bulkeset_subset)[, WT_HET_8wN.bulkeset_subset$group == 'WT']
bulk.case.mtx = exprs(WT_HET_8wN.bulkeset_subset)[, WT_HET_8wN.bulkeset_subset$group == 'HET']

# music2 deconvolution music2_prop_t_statistics
set.seed(42)
CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment$sampleID <- rownames(colData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment))

print("running est")


est = music2_prop_t_statistics(bulk.control.mtx = bulk.control.mtx, bulk.case.mtx = bulk.case.mtx, sc.sce = CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment, clusters = 'celltype_cl_coarse2', samples = 'sampleID', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "IP", "INP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "oRG", "RG","vRG")
, n_resample=20, sample_prop=0.5,cutoff_c=0.05,cutoff_r=0.01)


print("saving est")


# save image just in case
save(est, file = "output/MuSiC2/est_8wN_WTvsHET_CTRL.RData")
