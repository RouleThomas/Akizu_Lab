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


load("output/condiments/condiments_RNA_WT_Bap1KO_Part_PyNs_RSC_UL.Rdata")


set.seed(42)


#### DEGs trajectory per trajectory
counts <- RNA_WT_Bap1KO.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(RNA_WT_Bap1KO.sct$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(Part_PyNs_RSC_UL, na = FALSE) [,3]
cellweights <- slingCurveWeights(Part_PyNs_RSC_UL) [,3]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj_Part_PyNs_RSC_UL <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, 
     nknots = 6,
     sce = TRUE
   )

saveRDS(traj_Part_PyNs_RSC_UL, file = "output/condiments/traj_Part_PyNs_RSC_UL.rds")


