#! /usr/bin/env Rscript


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



load("output/condiments/condiments_multiome_WT_Bap1KO_QCV2vC1.RData")

set.seed(42)


#### DEGs trajectory per trajectory
counts <- multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident) # identify conditions

#### Extract the pseudotimes and cell weights for the 2 lineage
pseudotimes <- slingPseudotime(RNA_WT_Bap1KO, na = FALSE) [,2]  # HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cellweights <- slingCurveWeights(RNA_WT_Bap1KO) [,2] # HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj2 <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, 
     nknots = 6,
     sce = TRUE
   )

saveRDS(traj2, file = "output/condiments/traj2_RNA_common.rds")


