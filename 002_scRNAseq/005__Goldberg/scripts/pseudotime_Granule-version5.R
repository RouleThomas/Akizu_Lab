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
#WT_Kcnc1_CB_integrateMerge.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_integrateMerge-dim40kparam15res03-labelv1.rds")
#WT_Kcnc1_CB_integrateMerge.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_integrateMerge-version3QCversion4dim50kparam30res05-V1_numeric.rds")
#WT_Kcnc1_CB_integrateMerge.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_integrateMerge-version3QCversion4dim50kparam30res12-V1_numeric.rds")

#WT_Kcnc1_CB_integrateMerge.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_integrateMerge-version3QCversion4dim50kparam30res25-V1_numeric.rds")

WT_Kcnc1_CB_integrateMerge.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_integrateMerge-version5dim50kparam30res25-V1_numeric.rds")


DefaultAssay(WT_Kcnc1_CB_integrateMerge.sct) <- "RNA" # According to condiments workflow

# convert to SingleCellExperiment
WT_Kcnc1_CB <- as.SingleCellExperiment(WT_Kcnc1_CB_integrateMerge.sct, assay = "RNA")



## SEPARATE CELLS for each trajectory ############################

########################################################
############################ Granule ############################
########################################################

# First filter based on cell type
Part_Granule <- WT_Kcnc1_CB[, WT_Kcnc1_CB$seurat_clusters %in% c("26","21","16","10","5","12","14","6","1","36","13","8","15","11","3","34")]
table(Part_Granule$seurat_clusters) # to double check


# tidy
df <- bind_cols(
  as.data.frame(reducedDims(Part_Granule)$UMAP),
  as.data.frame(colData(Part_Granule)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
pdf("output/condiments/UMAP_WT_Kcnc1_CB-version5dim50kparam30res25-Part_Granule.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = seurat_clusters)) +
  geom_point(size = .7) +
  labs(col = "Genotype") +
  theme_classic()
dev.off()



## Second filter based on UMAP coordinate
umap_coords <- reducedDims(Part_Granule)$UMAP

# Filter conditions based on your description:
# Keep cells with UMAP_1 > -3 and UMAP_2 < 2.5

selected_cells <- umap_coords[,1] > -8 & umap_coords[,2] > 1 &  umap_coords[,1] < 5.25

# Subset your SCE object
Part_Granule_subset <- Part_Granule[, selected_cells]

# Check resulting subset
dim(Part_Granule_subset)

df <- bind_cols(
  as.data.frame(reducedDims(Part_Granule_subset)$UMAP),
  as.data.frame(colData(Part_Granule_subset)[, -3])
  ) %>%
  sample_frac(1)




## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$condition,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores






## PLOT with Separate trajectories
### Testing Area ############


##########################################


Part_Granule_subset <- slingshot(Part_Granule_subset, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_Granule_subset)$seurat_clusters,
                 start.clus = "26", end.clus = c("11") ,approx_points = 100, extend = 'pc1', stretch = 1)




#test reduceDim PCA or subset endoderm
topologyTest(SlingshotDataSet(Part_Granule_subset), Part_Granule_subset$condition) #  


sdss <- slingshot_conditions(SlingshotDataSet(Part_Granule_subset), Part_Granule_subset$condition)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "condition")



#  

pdf("output/condiments/UMAP_trajectory_separated_WT_Kcnc1_CB-version5dim50kparam30res25-Part_Granule_subset-START26_END11_points100extendpc1stretch1.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  theme_classic()
dev.off()




save.image(file="output/condiments/condiments-Part_Granule_subset_START26_END11_points100extendpc1stretch1-version5dim50kparam30res25.RData")





