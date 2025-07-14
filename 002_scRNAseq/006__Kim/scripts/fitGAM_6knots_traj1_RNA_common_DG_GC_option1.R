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



load("output/condiments/condiments-Part_DG_GC_option1_subset_STARTNSCquiesc_ENDDGGC_points100extendpc1stretch1distMetmnn-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.RData")
set.seed(42)


#### DEGs trajectory per trajectory
counts <- multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the FIRST lineage
pseudotimes <- slingPseudotime(Part_DG_GC_subset, na = FALSE) [,1]
cellweights <- slingCurveWeights(Part_DG_GC_subset) [,1]
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

saveRDS(traj1, file = "output/condiments/traj1_Part_DG_GC_option1_subset-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.rds")



