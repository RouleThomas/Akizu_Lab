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


print("importing music2 custom function")
my_music2_prop_toast_V2 <- function (bulk.control.mtx, bulk.case.mtx, sc.sce, clusters, 
          samples, select.ct, expr_low = 20, prop_r = 0.1, eps_c = 0.05, 
          eps_r = 0.01, cutoff_c = 10^(-3), cutoff_r = 10^(-3), cap = 0.3, 
          maxiter = 200, markers = NULL, ct.cov = FALSE, cell_size = NULL,
          centered = FALSE, normalize = FALSE) {
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if (length(gene.bulk) < 0.1 * min(nrow(bulk.control.mtx), 
                                    nrow(bulk.case.mtx))) {
    stop("Not enough genes for bulk data! Please check gene annotations.")
  }
  bulk.mtx = cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, 
  ])
  Pheno = data.frame(condition = factor(c(rep("control", ncol(bulk.control.mtx)), 
                                          rep("case", ncol(bulk.case.mtx))), levels = c("control", 
                                                                                        "case")))
  rownames(Pheno) = colnames(bulk.mtx)
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if (length(gene_all) < 0.2 * min(length(gene.bulk), nrow(sc.sce))) {
    stop("Not enough genes between bulk and single-cell data! Please check gene annotations.")
  }
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr >= expr_low])
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)]
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)]
  prop_control = music_prop(bulk.mtx = bulk.control, sc.sce = sc.sce, 
                            clusters = clusters, samples = samples, select.ct = select.ct, 
                            markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                            iter.max = 1000, nu = 1e-04, eps = 0.01, centered = centered, 
                            normalize = normalize, verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.mtx = bulk.case, sc.sce = sc.sce, 
                             clusters = clusters, samples = samples, select.ct = select.ct, 
                             markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                             iter.max = 1000, nu = 1e-04, eps = 0.01, centered = centered, 
                             normalize = normalize, verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control, prop_CASE)
  iter = 1
  ncell = length(select.ct)
  id_conv = NULL
  while (iter <= maxiter) {
    Y_raw = log1p(bulk.mtx)
    design = Pheno
    Prop <- prop_all[rownames(Pheno), ]
    Design_out <- makeDesign(design, Prop)
    fitted_model <- fitModel(Design_out, Y_raw)
    res_table <- csTest(fitted_model, coef = "condition", verbose = F)
    mex = apply(prop_all, 2, mean)
    lr = NULL
    for (celltype in select.ct) {
      m = mex[celltype]
      DE = res_table[[celltype]]
      pval = DE$fdr
      names(pval) = rownames(DE)
      pval = pval[names(pval) %in% exp_genel]
      if (m >= prop_r) {
        lr = c(lr, names(pval[pval <= cutoff_c & pval <= 
                                quantile(pval, prob = cap)]))
      }
      else {
        lr = c(lr, names(pval[pval <= cutoff_r & pval <= 
                                quantile(pval, prob = cap)]))
      }
    }
    lr = unique(lr)
    l = setdiff(gene_all, lr)
    sc.iter.sce = sc.sce[l, ]
    if (length(id_conv) > 0) {
      case_sample = bulk.case[, !colnames(bulk.case) %in% 
                                id_conv]
    }
    else {
      case_sample = bulk.case
    }
    prop_case = music_prop(bulk.mtx = case_sample, sc.sce = sc.iter.sce, 
                           clusters = clusters, samples = samples, select.ct = select.ct, 
                           verbose = F)$Est.prop.weighted
    prop_CASE = rbind(prop_case, prop_case_fix)
    if (length(id_conv) == 1) {
      rownames(prop_CASE) = c(rownames(prop_case), id_conv)
    }
    prop_all = rbind(prop_control, prop_CASE)
    prop_case = prop_case[rownames(prop_case_ini), ]

    replace_problematic_values <- function(x, replace_with = 0) {
      ifelse(is.na(x) | is.nan(x) | is.infinite(x), replace_with, x)
    }

    # Ensure 'prop_case' and 'prop_case_ini' do not contain problematic values
    prop_case <- replace_problematic_values(prop_case)
    prop_case_ini <- replace_problematic_values(prop_case_ini)

    # Use 'safe_division' function to avoid division by zero
    safe_division <- function(x, y) {
      ifelse(y == 0, 0, x / y)
    }

    # Calculate absolute differences and replace problematic values with zero
    pc <- abs(prop_case - prop_case_ini)
    pc <- replace_problematic_values(pc)

    # Initialize convergence matrix with ones
    conv <- matrix(1, nrow = nrow(pc), ncol = ncol(pc))

    # Update the convergence matrix based on conditions
    conv[prop_case_ini <= prop_r] <- ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
    pc[prop_case_ini > prop_r] <- safe_division(pc[prop_case_ini > prop_r], prop_case_ini[prop_case_ini > prop_r])
    conv[prop_case_ini > prop_r] <- ifelse(pc[prop_case_ini > prop_r] < eps_c, 0, 1)


                                               
    convf = apply(conv, 1, function(x) {
      all(x == 0)
    })
    all_converge = FALSE
    id_conv = c(id_conv, names(convf[convf == TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% 
                                id_conv, ]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv, 
    ]
    if (is.vector(prop_case_ini)) {
      all_converge = TRUE
      break
    }
    else if (nrow(prop_case_ini) == 0) {
      all_converge = TRUE
      break
    }
    iter = iter + 1
  }
  if (all_converge) {
    return(list(Est.prop = prop_all, convergence = TRUE, 
                n.iter = iter, DE.genes = lr))
  }
  else {
    return(list(Est.prop = prop_all, convergence = FALSE, 
                id.not.converge = rownames(prop_case_ini)))
  }
}



print("running est")


set.seed(42)
CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment$sampleID <- rownames(colData(CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment))

est = my_music2_prop_toast_V2(bulk.control.mtx = bulk.control.mtx, bulk.case.mtx = bulk.case.mtx, sc.sce = CHOOSE_CTRL_annot_srt.subset_celltype_cl_coarse2.SingleCellExperiment, clusters = 'celltype_cl_coarse2', samples = 'sampleID', select.ct = c("Astrocytes", "ccRG", "ccvRG", "CGE_IN", "CGE_LGE_IN", "IP", "INP", "L23", "L4", "L56", "L6_CThPN", "LGE_IN", "oRG", "RG","vRG"))


print("saving est")


# save image just in case
save(est, file = "output/MuSiC2/est_8wN_WTvsHET_CTRL_my_music2_prop_toast_V2.RData")
