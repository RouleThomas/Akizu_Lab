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



load("output/condiments/condiments_embryo_V2clust.RData")

set.seed(42)


#### DEGs trajectory per trajectory
counts <- embryo.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(embryo.combined.sct$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the 2nd lineage
pseudotimes <- slingPseudotime(embryo, na = FALSE) [,2]
cellweights <- slingCurveWeights(embryo) [,2]
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

saveRDS(traj2, file = "output/condiments/traj2.rds")


