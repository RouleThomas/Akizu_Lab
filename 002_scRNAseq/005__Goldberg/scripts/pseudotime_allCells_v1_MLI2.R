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
set.seed(42)



# Data import - all samples and genotype CB
WT_Kcnc1_CB_integrateMerge.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_integrateMerge-dim40kparam15res03-labelv1.rds")


DefaultAssay(WT_Kcnc1_CB_integrateMerge.sct) <- "RNA" # According to condiments workflow


# convert to SingleCellExperiment
WT_Kcnc1_CB <- as.SingleCellExperiment(WT_Kcnc1_CB_integrateMerge.sct, assay = "RNA")


# tidy
df <- bind_cols(
  as.data.frame(reducedDims(WT_Kcnc1_CB)$UMAP),
  as.data.frame(colData(WT_Kcnc1_CB)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$condition,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

########################################
## PLOT keeping all cells ##########
########################################
# Testing
WT_Kcnc1_CB <- slingshot(WT_Kcnc1_CB, reducedDim = 'UMAP',
                 clusterLabels = colData(WT_Kcnc1_CB)$cluster.annot,
                 start.clus = c("MLI2_3"), end.clus = c("MLI2_2") ,approx_points = 100, extend = 'n', stretch = 1)

#test reduceDim PCA or subset endoderm
topologyTest(SlingshotDataSet(WT_Kcnc1_CB), WT_Kcnc1_CB$condition) #  
sdss <- slingshot_conditions(SlingshotDataSet(WT_Kcnc1_CB), WT_Kcnc1_CB$condition)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "condition")
pdf("output/condiments/UMAP_trajectory_WT_Kcnc1_CB-dim40kparam15res03-START-MLI2_3-END-MLI2_2-points100extendnstretch1.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  theme_classic()
dev.off()


# Save RData
save.image(file = "output/condiments/pseudotime_allCells_v1_MLI2.RData")



