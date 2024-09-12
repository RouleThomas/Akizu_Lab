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
library("cowplot")
library("scales")
library("pheatmap")
library("readr")



################### Time Course effect COMMON CONDITIONS ######################################################

## TRAJECTORY3 ##################
set.seed(42)
traj5_noCondition_humangastruloid2472hrs <- readRDS("output/condiments/traj5_noCondition_humangastruloid2472hrs.rds")



## Identify Activation point = peak (maximum expression) of each gene along the pseudotime trajectory
### Identify peak of expression (max expr) of these Time-course DEG
traj5_noCondition_humangastruloid2472hrs
#### Extract pseudotime values
pseudotime <- colData(traj5_noCondition_humangastruloid2472hrs)$crv$pseudotime
#### Extract the expression matrix
expr_matrix <- assays(traj5_noCondition_humangastruloid2472hrs)$counts
#### Ensure the pseudotime values are named with the same cell names as the expression matrix columns
names(pseudotime) <- colnames(expr_matrix)
#### Function to find the peak pseudotime for each gene (raw and smoothed)
find_max_pseudotime <- function(gene_expr, pseudotime) {
  # Raw peak pseudotime
  raw_peak_pseudotime <- pseudotime[which.max(gene_expr)]
  # Smooth gene expression using loess
  smooth_model <- loess(gene_expr ~ pseudotime)
  smooth_expr <- predict(smooth_model)
  # Smooth peak pseudotime
  smooth_peak_pseudotime <- pseudotime[which.max(smooth_expr)]
  return(list(raw_peak_pseudotime = raw_peak_pseudotime, 
              smooth_peak_pseudotime = smooth_peak_pseudotime))
}
#### Apply the function to all genes
peak_values <- apply(expr_matrix, 1, function(x) find_max_pseudotime(as.numeric(x), pseudotime))
#### Convert the results to a data frame
peak_df <- data.frame(
  gene = rownames(expr_matrix),
  raw_peak_pseudotime = sapply(peak_values, `[[`, "raw_peak_pseudotime"),
  smooth_peak_pseudotime = sapply(peak_values, `[[`, "smooth_peak_pseudotime")
) %>% as_tibble()

write.table(peak_df, file = c("output/condiments/traj5_noCondition_humangastruloid2472hrs_ActivationPoint.txt"),sep="\t", quote=FALSE, row.names=FALSE)
