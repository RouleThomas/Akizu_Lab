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



load("output/condiments/condiments_humangastruloidDASATINIB72hr.RData")

set.seed(42)


#### DEGs trajectory per trajectory
counts <- humangastruloid_DASATINIB72hr.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(humangastruloid_DASATINIB72hr.combined.sct$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(humangastruloid_DASATINIB72hr, na = FALSE) [,2]
cellweights <- slingCurveWeights(humangastruloid_DASATINIB72hr) [,2]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj2_humangastruloidDASATINIB72hr <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = NULL, 
     nknots = 6,
     sce = TRUE
   )

saveRDS(traj2_humangastruloidDASATINIB72hr, file = "output/condiments/traj2_humangastruloidDASATINIB72hr.rds")


