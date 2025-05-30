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




load("output/condiments/condiments-Part_MLI2_subset_START27_END19_points100extendpc1stretch1-version5dim50kparam30res25.RData")




#### DEGs trajectory per trajectory
counts <- WT_Kcnc1_CB_integrateMerge.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(WT_Kcnc1_CB_integrateMerge.sct$condition) # identify conditions
#### Extract the pseudotimes and cell weights for the FIRST lineage
pseudotimes <- slingPseudotime(Part_MLI2_subset, na = FALSE) [,1]
cellweights <- slingCurveWeights(Part_MLI2_subset) [,1]
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

saveRDS(traj1, file = "output/condiments/traj1_Part_MLI2_subset-version5dim50kparam30res25.rds")


