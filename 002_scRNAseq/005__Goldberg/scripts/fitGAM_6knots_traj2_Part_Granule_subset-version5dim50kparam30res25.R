


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



load("output/condiments/condiments-Part_Granule_subset_START26_END11_points100extendpc1stretch1-version5dim50kparam30res25.RData")
set.seed(42)



#### DEGs trajectory per trajectory
counts <- WT_Kcnc1_CB_integrateMerge.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(WT_Kcnc1_CB_integrateMerge.sct$condition) # identify conditions
#### Extract the pseudotimes and cell weights for the second lineage
pseudotimes <- slingPseudotime(Part_Granule_subset, na = FALSE) [,2]
cellweights <- slingCurveWeights(Part_Granule_subset) [,2]
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

saveRDS(traj2, file = "output/condiments/traj2_Part_Granule_subset-version5dim50kparam30res25.rds")


