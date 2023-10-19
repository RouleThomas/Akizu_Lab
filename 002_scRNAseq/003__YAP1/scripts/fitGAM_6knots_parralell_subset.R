#! /usr/bin/env Rscript



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
library("mgcv")

set.seed(42)
# Load
load("output/condiments/condiments_embryo_V2clust.RData")


set.seed(42)
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 8

### Define the desired lineages to keep
desired_lineages <- paste0("Lineage", 1:5)
### Extract pseudotimes and cellweights for the desired lineages from the embryo object
pseudotimes <- slingPseudotime(SlingshotDataSet(embryo), na = FALSE) [, desired_lineages]
cellweights <- slingCurveWeights(embryo)[, desired_lineages]
### Retain only those cells with non-zero weights
sub_weights <- cellweights[rowSums(cellweights != 0) > 0,]
sub_pseudotimes <- pseudotimes[rownames(pseudotimes) %in% rownames(sub_weights),]
### Subset the embryo SingleCellExperiment object to retain only relevant cells
embryo_sub <- embryo[,colnames(embryo) %in% rownames(sub_weights)]
### Extract the counts
counts <- assays(embryo_sub)$counts %>% as.matrix()

# Set up the conditions and GAM control
my_conditions <- as.factor(embryo_sub$condition)
control <- gam.control()
control$maxit <- 1000


fitgam <-  fitGAM(counts = counts,
                  pseudotime = sub_pseudotimes,
                  cellWeights = sub_weights,
                  control = control, 
                  conditions = my_conditions,
                  nknots = 6,
                  parallel=TRUE,
                  BPPARAM = BPPARAM)




condRes <- conditionTest(fitgam, l2fc = log2(2), lineages = TRUE)

condRes$padj <- p.adjust(condRes$pvalue_lineage1, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)
conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

write.table(conditionGenes, file = "output/condiments/conditionGenes_nknots6.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# save image just in case
save.image(file="output/condiments/condiments_embryo_V2clust_fitGAM.RData")

