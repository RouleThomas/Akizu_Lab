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


set.seed(42)
# Load
load("output/condiments/condiments_embryo_V1.RData")

# run DEGs
embryo <- fitGAM(counts = embryo, conditions = factor(embryo$condition), nknots = 6) # change nknots here
condRes <- conditionTest(embryo, l2fc = log2(2), lineages = TRUE)

condRes$padj <- p.adjust(condRes$pvalue_lineage1, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)
conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

write.table(conditionGenes, file = "output/condiments/conditionGenes_nknots6.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# save image just in case
save.image(file="output/condiments/condiments_embryo_V2.RData")

