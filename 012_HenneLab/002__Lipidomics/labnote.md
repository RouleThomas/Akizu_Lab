# Project

Samples and conditions:
- genotypes (3): WT, KD SNX13, KD SNX14
- condition (3): DMSO, LLOME, recovery
- fraction (2): Cell, lysosome

**Components class** include:
- TAG - Triacylglycerol (or Triglyceride) 
- CE - Cholesteryl ester 
- DAG - Diacylglycerol 
- SM - Sphingomyelin 
- LPC - Lysophosphatidylcholine 
- PC - Phosphatidylcholine 
- LPE - Lysophosphatidylethanolamine 
- PG - Phosphatidylglycerol 
- PI - Phosphatidylinositol 
- PS - Phosphatidylserine 
- PA - Phosphatidic acid 
- LPG - Lysophosphatidylglycerol 
- LPI - Lysophosphatidylinositol 
- LPS - Lysophosphatidylserine 
- BMP - Bis(monoacylglycero)phosphate 
- PE - Phosphatidylethanolamine lipids 



**Output values** available:
- Raw abundance of each component
    - Summary values of degree of saturation for each class of components


--> Run independently limma for each fraction; and run independent comparison WT vs siSNX13 and WT vs siSNX14


--> Let's manually tidy the Xcell file as `For Thomas_Lipidomics_tidy.xlsx`: 
- rename all columns name by removing any special characters and space + rename all; columns with replicate name
- Added column Class (ie. TAG, CE, ...), Size (chain length), Saturation (Unsaturated or Saturated)


Let's try to follow the LOGLIMMA method used in Proteomics; we will need to **normalize value per total lipid content for each class**.





# Testing - v1 - LOGLIMMA



```bash
conda activate deseq2
```



```R

# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("limma")
library("readxl")




# Import file
lip <- read_excel(
  path  = "input/For_Thomas_Lipidomics_tidy.xlsx",
  sheet = "tidy clean"
)

# Tidy files
sample_cols <- names(lip)[str_detect(names(lip), "^(Lysosomes|Cells)-")]

lip_tidy <- lip %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance",
    values_transform = list(abundance = as.numeric)
  ) %>%
  extract(
    col = sample,
    into = c("fraction", "condition", "genotype", "replicate"),
    regex = "^(Lysosomes|Cells)-([A-Za-z0-9]+)_([^&]+)&(\\d+)$",
    remove = FALSE
  ) 



########################################
## Normalize to sum per Class ##########
########################################


lip_tidy_norm <- lip_tidy %>%
  group_by(
    sample,
    fraction,
    condition,
    genotype,
    replicate,
    Class
  ) %>%
  mutate(
    Class_total     = sum(abundance, na.rm = TRUE),
    abundance_Class = ifelse(Class_total > 0,
                             abundance / Class_total,
                             NA_real_)
  ) %>%
  ungroup()






########################################
## Scram vs siSNX13 - Cells ####################
########################################

df <- lip_tidy_norm %>%
  filter(fraction == "Cells",
         genotype %in% c("Scram", "siSNX13")) %>%
  dplyr::select(Class, Size, Saturation, Component, sample, condition, genotype, replicate, abundance_Class)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Component, sample, abundance_Class) %>%
  pivot_wider(names_from = sample, values_from = abundance_Class) %>%
  column_to_rownames("Component") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)


## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))

stopifnot(all(coldata$sample == colnames(log2_expr)))

rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = c("DMSO","LLOME","RECOVERY"))
coldata$genotype  <- factor(coldata$genotype,  levels = c("Scram","siSNX13"))
coldata$replicate <- factor(coldata$replicate)



## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr, design_full)
fit <- eBayes(fit)

## 7) "LRT-like" global test: are interaction terms jointly != 0?
## This is the analogue of testing full vs reduced (~ genotype + condition)
int_cols <- grep("^genotypesiSNX13:condition", colnames(design_full))

interaction_global <- topTableF(
  fit,
  int_cols,
  number = Inf
) 
head(interaction_global)
#--> Signficant Component = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Component") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX13.conditionLLOME),
      abs(genotypesiSNX13.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 0)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Component

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 257




# Build ordered sample table
sample_order <- coldata %>%
  mutate(
    condition = factor(condition, levels = c("DMSO","LLOME","RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram","siSNX13"))
  ) %>%
  arrange(genotype, condition, replicate)

# Reorder matrix
log2_expr_sig_ord <- log2_expr_sig[, sample_order$sample]




expr_scaled <- t(scale(t(log2_expr_sig_ord)))

row_dist   <- dist(expr_scaled, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")




k <- 9
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Component = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Component)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Component")


n_rep <- sample_order %>%
  count(genotype, condition)


block_sizes <- sample_order %>%
  count(genotype, condition, .drop = FALSE) %>%   # add fraction here too if needed
  pull(n)

gaps_col <- cumsum(block_sizes)[-length(block_sizes)] 

row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-SNX13_Cells-9Cluster.pdf", width = 4, height = 4)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = gaps_col,
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()
#--> Clean, no outliers



out_tbl <- cluster_tbl %>%
  left_join(df %>% dplyr::select(Component, Class, Size, Saturation) %>% unique) %>%
  arrange(cluster, Component)

write_tsv(out_tbl, "output/GENELIST-LOGLIMMA-SNX13_Cells-9Cluster.tsv")



pca <- prcomp(t(log2_expr_sig), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-SNX13_Cells.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...




# Show top 10 Component per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Component, cluster), by = "Component") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX13")),
    log2_abund = log2(abundance_Class + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Component) %>%
  summarise(mean_abund = mean(abundance_Class, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Component")) %>%
  group_by(cluster, Component, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(" (", Component, ")")
  )
## Plot: one PDF per cluster (top 10 panels)
for (cl in sort(unique(plot_df$cluster))) {
  p <- plot_df %>%
    filter(cluster == cl) %>%
    ggplot(aes(x = condition, y = mean_log2, color = genotype, group = genotype)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.6) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Scram" = "black", "siSNX13" = "red")) +
    labs(
      title = paste0("Top 10 lipids — Cluster ", cl, " (Scram vs siSNX13)"),
      x = NULL,
      y = "log2(Normalized abundance + 1)",
      color = "genotype"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX13_Cells-9Cluster.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}






# Global lipid characteristics
comp_annot_all <- df %>%
  distinct(Component, Class, Size, Saturation)

prop_class_all <- comp_annot_all %>%
  count(Class, name = "n") %>%
  mutate(prop = n / sum(n))
prop_size_all <- comp_annot_all %>%
  count(Size, name = "n") %>%
  mutate(prop = n / sum(n))
prop_sat_all <- comp_annot_all %>%
  count(Saturation, name = "n") %>%
  mutate(prop = n / sum(n))

pdf("output/global_composition-Cells.pdf", width = 3, height = 4)
ggplot(prop_class_all, aes(x = "All components", y = prop, fill = Class)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global lipid class composition") +
  theme_classic()
ggplot(prop_size_all, aes(x = "All components", y = prop, fill = Size)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global size composition") +
  theme_classic()
ggplot(prop_sat_all, aes(x = "All components", y = prop, fill = Saturation)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global saturation composition") +
  theme_classic()
dev.off()


# Investigate lipid characteristics of each clusters
comp_annot_clust <- df %>%
  inner_join(cluster_tbl %>% select(Component, cluster), by = "Component") %>%
  mutate(cluster = factor(cluster, levels = sort(unique(cluster))))


prop_class <- comp_annot_clust %>%
  count(cluster, Class, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_size <- comp_annot_clust %>%
  count(cluster, Size, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_sat <- comp_annot_clust %>%
  count(cluster, Saturation, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


pdf("output/cluster_composition-LOGLIMMA-SNX13_Cells-9Cluster.pdf", width = 6, height = 4)
ggplot(prop_class, aes(x = cluster, y = prop, fill = Class)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Lipid Class composition per cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
ggplot(prop_size, aes(x = cluster, y = prop, fill = Size)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Size composition per cluster") +
  theme_classic()
ggplot(prop_sat, aes(x = cluster, y = prop, fill = Saturation)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Saturation composition per cluster") +
  theme_classic()
dev.off()






comp_annot_clust_all <- comp_annot_clust %>%
  select(Component, Class, Size, Saturation, cluster) %>%
  bind_rows(comp_annot_all)

cluster_levels <- c("Global", sort(unique(cluster_tbl$cluster)))

comp_annot_clust_all <- comp_annot_clust_all %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  mutate(
    cluster = ifelse(is.na(cluster), "Global", as.character(cluster)),
    cluster = factor(cluster, levels = c("Global", sort(unique(cluster[cluster != "Global"]))))
  )


prop_class <- comp_annot_clust_all %>%
  count(cluster, Class, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_size <- comp_annot_clust_all %>%
  count(cluster, Size, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_sat <- comp_annot_clust_all %>%
  count(cluster, Saturation, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

pdf("output/cluster_composition_andGlobal-LOGLIMMA-SNX13_Cells-9Cluster.pdf", width = 6, height = 4)
ggplot(prop_class, aes(x = cluster, y = prop, fill = Class)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Lipid Class composition (clusters vs global)") +
  theme_classic()
ggplot(prop_size, aes(x = cluster, y = prop, fill = Size)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Size composition (clusters vs global)") +
  theme_classic()
ggplot(prop_sat, aes(x = cluster, y = prop, fill = Saturation)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Saturation composition (clusters vs global)") +
  theme_classic()
dev.off()






########################################
## Scram vs siSNX13 - Lysosomes ####################
########################################

df <- lip_tidy_norm %>%
  filter(fraction == "Lysosomes",
         genotype %in% c("Scram", "siSNX13")) %>%
  dplyr::select(Class, Size, Saturation, Component, sample, condition, genotype, replicate, abundance_Class)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Component, sample, abundance_Class) %>%
  pivot_wider(names_from = sample, values_from = abundance_Class) %>%
  column_to_rownames("Component") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)


## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))

stopifnot(all(coldata$sample == colnames(log2_expr)))

rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = c("DMSO","LLOME","RECOVERY"))
coldata$genotype  <- factor(coldata$genotype,  levels = c("Scram","siSNX13"))
coldata$replicate <- factor(coldata$replicate)



## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr, design_full)
fit <- eBayes(fit)

## 7) "LRT-like" global test: are interaction terms jointly != 0?
## This is the analogue of testing full vs reduced (~ genotype + condition)
int_cols <- grep("^genotypesiSNX13:condition", colnames(design_full))

interaction_global <- topTableF(
  fit,
  int_cols,
  number = Inf
) 
head(interaction_global)
#--> Signficant Component = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Component") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX13.conditionLLOME),
      abs(genotypesiSNX13.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 0)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Component

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 257




# Build ordered sample table
sample_order <- coldata %>%
  mutate(
    condition = factor(condition, levels = c("DMSO","LLOME","RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram","siSNX13"))
  ) %>%
  arrange(genotype, condition, replicate)

# Reorder matrix
log2_expr_sig_ord <- log2_expr_sig[, sample_order$sample]




expr_scaled <- t(scale(t(log2_expr_sig_ord)))

row_dist   <- dist(expr_scaled, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")




k <- 8
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Component = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Component)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Component")


n_rep <- sample_order %>%
  count(genotype, condition)


block_sizes <- sample_order %>%
  count(genotype, condition, .drop = FALSE) %>%   # add fraction here too if needed
  pull(n)

gaps_col <- cumsum(block_sizes)[-length(block_sizes)] 

row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-SNX13_Lysosomes-8Cluster.pdf", width = 4, height = 4)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = gaps_col,
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()
#--> Clean, no outliers



out_tbl <- cluster_tbl %>%
  left_join(df %>% dplyr::select(Component, Class, Size, Saturation) %>% unique) %>%
  arrange(cluster, Component)

write_tsv(out_tbl, "output/GENELIST-LOGLIMMA-SNX13_Lysosomes-8Cluster.tsv")





pca <- prcomp(t(log2_expr_sig), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-SNX13_Lysosomes.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...



# Show top 10 Component per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Component, cluster), by = "Component") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX13")),
    log2_abund = log2(abundance_Class + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Component) %>%
  summarise(mean_abund = mean(abundance_Class, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Component")) %>%
  group_by(cluster, Component, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(" (", Component, ")")
  )
## Plot: one PDF per cluster (top 10 panels)
for (cl in sort(unique(plot_df$cluster))) {
  p <- plot_df %>%
    filter(cluster == cl) %>%
    ggplot(aes(x = condition, y = mean_log2, color = genotype, group = genotype)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.6) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Scram" = "black", "siSNX13" = "red")) +
    labs(
      title = paste0("Top 10 lipids — Cluster ", cl, " (Scram vs siSNX13)"),
      x = NULL,
      y = "log2(Normalized abundance + 1)",
      color = "genotype"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX13_Lysosomes-8Cluster.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}




# Global lipid characteristics
comp_annot_all <- df %>%
  distinct(Component, Class, Size, Saturation)

prop_class_all <- comp_annot_all %>%
  count(Class, name = "n") %>%
  mutate(prop = n / sum(n))
prop_size_all <- comp_annot_all %>%
  count(Size, name = "n") %>%
  mutate(prop = n / sum(n))
prop_sat_all <- comp_annot_all %>%
  count(Saturation, name = "n") %>%
  mutate(prop = n / sum(n))

pdf("output/global_composition-Lysosomes.pdf", width = 3, height = 4)
ggplot(prop_class_all, aes(x = "All components", y = prop, fill = Class)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global lipid class composition") +
  theme_classic()
ggplot(prop_size_all, aes(x = "All components", y = prop, fill = Size)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global size composition") +
  theme_classic()
ggplot(prop_sat_all, aes(x = "All components", y = prop, fill = Saturation)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global saturation composition") +
  theme_classic()
dev.off()


# Investigate lipid characteristics of each clusters
comp_annot_clust <- df %>%
  inner_join(cluster_tbl %>% select(Component, cluster), by = "Component") %>%
  mutate(cluster = factor(cluster, levels = sort(unique(cluster))))


prop_class <- comp_annot_clust %>%
  count(cluster, Class, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_size <- comp_annot_clust %>%
  count(cluster, Size, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_sat <- comp_annot_clust %>%
  count(cluster, Saturation, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


pdf("output/cluster_composition-LOGLIMMA-SNX13_Lysosomes-8Cluster.pdf", width = 6, height = 4)
ggplot(prop_class, aes(x = cluster, y = prop, fill = Class)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Lipid Class composition per cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
ggplot(prop_size, aes(x = cluster, y = prop, fill = Size)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Size composition per cluster") +
  theme_classic()
ggplot(prop_sat, aes(x = cluster, y = prop, fill = Saturation)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Saturation composition per cluster") +
  theme_classic()
dev.off()




comp_annot_clust_all <- comp_annot_clust %>%
  select(Component, Class, Size, Saturation, cluster) %>%
  bind_rows(comp_annot_all)

cluster_levels <- c("Global", sort(unique(cluster_tbl$cluster)))

comp_annot_clust_all <- comp_annot_clust_all %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  mutate(
    cluster = ifelse(is.na(cluster), "Global", as.character(cluster)),
    cluster = factor(cluster, levels = c("Global", sort(unique(cluster[cluster != "Global"]))))
  )


prop_class <- comp_annot_clust_all %>%
  count(cluster, Class, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_size <- comp_annot_clust_all %>%
  count(cluster, Size, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_sat <- comp_annot_clust_all %>%
  count(cluster, Saturation, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

pdf("output/cluster_composition_andGlobal-LOGLIMMA-SNX13_Lysosomes-8Cluster.pdf", width = 6, height = 4)
ggplot(prop_class, aes(x = cluster, y = prop, fill = Class)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Lipid Class composition (clusters vs global)") +
  theme_classic()
ggplot(prop_size, aes(x = cluster, y = prop, fill = Size)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Size composition (clusters vs global)") +
  theme_classic()
ggplot(prop_sat, aes(x = cluster, y = prop, fill = Saturation)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Saturation composition (clusters vs global)") +
  theme_classic()
dev.off()



########################################
## Scram vs siSNX14 - Cells ####################
########################################

df <- lip_tidy_norm %>%
  filter(fraction == "Cells",
         genotype %in% c("Scram", "siSNX14")) %>%
  dplyr::select(Class, Size, Saturation, Component, sample, condition, genotype, replicate, abundance_Class)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Component, sample, abundance_Class) %>%
  pivot_wider(names_from = sample, values_from = abundance_Class) %>%
  column_to_rownames("Component") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)


## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))

stopifnot(all(coldata$sample == colnames(log2_expr)))

rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = c("DMSO","LLOME","RECOVERY"))
coldata$genotype  <- factor(coldata$genotype,  levels = c("Scram","siSNX14"))
coldata$replicate <- factor(coldata$replicate)



## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr, design_full)
fit <- eBayes(fit)

## 7) "LRT-like" global test: are interaction terms jointly != 0?
## This is the analogue of testing full vs reduced (~ genotype + condition)
int_cols <- grep("^genotypesiSNX14:condition", colnames(design_full))

interaction_global <- topTableF(
  fit,
  int_cols,
  number = Inf
) 
head(interaction_global)
#--> Signficant Component = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Component") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX14.conditionLLOME),
      abs(genotypesiSNX14.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 0)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Component

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 258




# Build ordered sample table
sample_order <- coldata %>%
  mutate(
    condition = factor(condition, levels = c("DMSO","LLOME","RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram","siSNX14"))
  ) %>%
  arrange(genotype, condition, replicate)

# Reorder matrix
log2_expr_sig_ord <- log2_expr_sig[, sample_order$sample]




expr_scaled <- t(scale(t(log2_expr_sig_ord)))

row_dist   <- dist(expr_scaled, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")




k <- 8
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Component = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Component)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Component")


n_rep <- sample_order %>%
  count(genotype, condition)


block_sizes <- sample_order %>%
  count(genotype, condition, .drop = FALSE) %>%   # add fraction here too if needed
  pull(n)

gaps_col <- cumsum(block_sizes)[-length(block_sizes)] 

row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-SNX14_Cells-8Cluster.pdf", width = 4, height = 4)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = gaps_col,
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()
#--> Clean, no outliers

out_tbl <- cluster_tbl %>%
  left_join(df %>% dplyr::select(Component, Class, Size, Saturation) %>% unique) %>%
  arrange(cluster, Component)

write_tsv(out_tbl, "output/GENELIST-LOGLIMMA-SNX14_Cells-8Cluster.tsv")


pca <- prcomp(t(log2_expr_sig), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-SNX14_Cells.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...



# Show top 10 Component per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Component, cluster), by = "Component") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX14")),
    log2_abund = log2(abundance_Class + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Component) %>%
  summarise(mean_abund = mean(abundance_Class, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Component")) %>%
  group_by(cluster, Component, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(" (", Component, ")")
  )
## Plot: one PDF per cluster (top 10 panels)
for (cl in sort(unique(plot_df$cluster))) {
  p <- plot_df %>%
    filter(cluster == cl) %>%
    ggplot(aes(x = condition, y = mean_log2, color = genotype, group = genotype)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.6) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Scram" = "black", "siSNX14" = "red")) +
    labs(
      title = paste0("Top 10 lipids — Cluster ", cl, " (Scram vs siSNX14)"),
      x = NULL,
      y = "log2(Normalized abundance + 1)",
      color = "genotype"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX14_Cells-8Cluster.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}






# Global lipid characteristics
comp_annot_all <- df %>%
  distinct(Component, Class, Size, Saturation)

prop_class_all <- comp_annot_all %>%
  count(Class, name = "n") %>%
  mutate(prop = n / sum(n))
prop_size_all <- comp_annot_all %>%
  count(Size, name = "n") %>%
  mutate(prop = n / sum(n))
prop_sat_all <- comp_annot_all %>%
  count(Saturation, name = "n") %>%
  mutate(prop = n / sum(n))

pdf("output/global_composition-Cells.pdf", width = 3, height = 4)
ggplot(prop_class_all, aes(x = "All components", y = prop, fill = Class)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global lipid class composition") +
  theme_classic()
ggplot(prop_size_all, aes(x = "All components", y = prop, fill = Size)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global size composition") +
  theme_classic()
ggplot(prop_sat_all, aes(x = "All components", y = prop, fill = Saturation)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global saturation composition") +
  theme_classic()
dev.off()


# Investigate lipid characteristics of each clusters
comp_annot_clust <- df %>%
  inner_join(cluster_tbl %>% select(Component, cluster), by = "Component") %>%
  mutate(cluster = factor(cluster, levels = sort(unique(cluster))))


prop_class <- comp_annot_clust %>%
  count(cluster, Class, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_size <- comp_annot_clust %>%
  count(cluster, Size, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_sat <- comp_annot_clust %>%
  count(cluster, Saturation, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


pdf("output/cluster_composition-LOGLIMMA-SNX14_Cells-8Cluster.pdf", width = 6, height = 4)
ggplot(prop_class, aes(x = cluster, y = prop, fill = Class)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Lipid Class composition per cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
ggplot(prop_size, aes(x = cluster, y = prop, fill = Size)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Size composition per cluster") +
  theme_classic()
ggplot(prop_sat, aes(x = cluster, y = prop, fill = Saturation)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Saturation composition per cluster") +
  theme_classic()
dev.off()





comp_annot_clust_all <- comp_annot_clust %>%
  select(Component, Class, Size, Saturation, cluster) %>%
  bind_rows(comp_annot_all)

cluster_levels <- c("Global", sort(unique(cluster_tbl$cluster)))

comp_annot_clust_all <- comp_annot_clust_all %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  mutate(
    cluster = ifelse(is.na(cluster), "Global", as.character(cluster)),
    cluster = factor(cluster, levels = c("Global", sort(unique(cluster[cluster != "Global"]))))
  )


prop_class <- comp_annot_clust_all %>%
  count(cluster, Class, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_size <- comp_annot_clust_all %>%
  count(cluster, Size, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_sat <- comp_annot_clust_all %>%
  count(cluster, Saturation, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

pdf("output/cluster_composition_andGlobal-LOGLIMMA-SNX14_Cells-8Cluster.pdf", width = 6, height = 4)
ggplot(prop_class, aes(x = cluster, y = prop, fill = Class)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Lipid Class composition (clusters vs global)") +
  theme_classic()
ggplot(prop_size, aes(x = cluster, y = prop, fill = Size)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Size composition (clusters vs global)") +
  theme_classic()
ggplot(prop_sat, aes(x = cluster, y = prop, fill = Saturation)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Saturation composition (clusters vs global)") +
  theme_classic()
dev.off()





########################################
## Scram vs siSNX14 - Lysosomes ####################
########################################

df <- lip_tidy_norm %>%
  filter(fraction == "Lysosomes",
         genotype %in% c("Scram", "siSNX14")) %>%
  dplyr::select(Class, Size, Saturation, Component, sample, condition, genotype, replicate, abundance_Class)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Component, sample, abundance_Class) %>%
  pivot_wider(names_from = sample, values_from = abundance_Class) %>%
  column_to_rownames("Component") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)


## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))

stopifnot(all(coldata$sample == colnames(log2_expr)))

rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = c("DMSO","LLOME","RECOVERY"))
coldata$genotype  <- factor(coldata$genotype,  levels = c("Scram","siSNX14"))
coldata$replicate <- factor(coldata$replicate)



## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr, design_full)
fit <- eBayes(fit)

## 7) "LRT-like" global test: are interaction terms jointly != 0?
## This is the analogue of testing full vs reduced (~ genotype + condition)
int_cols <- grep("^genotypesiSNX14:condition", colnames(design_full))

interaction_global <- topTableF(
  fit,
  int_cols,
  number = Inf
) 
head(interaction_global)
#--> Signficant Component = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Component") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX14.conditionLLOME),
      abs(genotypesiSNX14.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 0)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Component

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 261




# Build ordered sample table
sample_order <- coldata %>%
  mutate(
    condition = factor(condition, levels = c("DMSO","LLOME","RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram","siSNX14"))
  ) %>%
  arrange(genotype, condition, replicate)

# Reorder matrix
log2_expr_sig_ord <- log2_expr_sig[, sample_order$sample]




expr_scaled <- t(scale(t(log2_expr_sig_ord)))

row_dist   <- dist(expr_scaled, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")




k <- 8
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Component = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Component)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Component")


n_rep <- sample_order %>%
  count(genotype, condition)


block_sizes <- sample_order %>%
  count(genotype, condition, .drop = FALSE) %>%   # add fraction here too if needed
  pull(n)

gaps_col <- cumsum(block_sizes)[-length(block_sizes)] 

row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-SNX14_Lysosomes-8Cluster.pdf", width = 4, height = 4)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = gaps_col,
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()
#--> Clean, no outliers



out_tbl <- cluster_tbl %>%
  left_join(df %>% dplyr::select(Component, Class, Size, Saturation) %>% unique) %>%
  arrange(cluster, Component)

write_tsv(out_tbl, "output/GENELIST-LOGLIMMA-SNX14_Lysosomes-8Cluster.tsv")



pca <- prcomp(t(log2_expr_sig), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-SNX14_Lysosomes.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...



# Show top 10 Component per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Component, cluster), by = "Component") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX14")),
    log2_abund = log2(abundance_Class + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Component) %>%
  summarise(mean_abund = mean(abundance_Class, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Component")) %>%
  group_by(cluster, Component, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(" (", Component, ")")
  )
## Plot: one PDF per cluster (top 10 panels)
for (cl in sort(unique(plot_df$cluster))) {
  p <- plot_df %>%
    filter(cluster == cl) %>%
    ggplot(aes(x = condition, y = mean_log2, color = genotype, group = genotype)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.6) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Scram" = "black", "siSNX14" = "red")) +
    labs(
      title = paste0("Top 10 lipids — Cluster ", cl, " (Scram vs siSNX14)"),
      x = NULL,
      y = "log2(Normalized abundance + 1)",
      color = "genotype"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX14_Lysosomes-8Cluster.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}







# Global lipid characteristics
comp_annot_all <- df %>%
  distinct(Component, Class, Size, Saturation)

prop_class_all <- comp_annot_all %>%
  count(Class, name = "n") %>%
  mutate(prop = n / sum(n))
prop_size_all <- comp_annot_all %>%
  count(Size, name = "n") %>%
  mutate(prop = n / sum(n))
prop_sat_all <- comp_annot_all %>%
  count(Saturation, name = "n") %>%
  mutate(prop = n / sum(n))

pdf("output/global_composition-Lysosomes.pdf", width = 3, height = 4)
ggplot(prop_class_all, aes(x = "All components", y = prop, fill = Class)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global lipid class composition") +
  theme_classic()
ggplot(prop_size_all, aes(x = "All components", y = prop, fill = Size)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global size composition") +
  theme_classic()
ggplot(prop_sat_all, aes(x = "All components", y = prop, fill = Saturation)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of components", title = "Global saturation composition") +
  theme_classic()
dev.off()




# Investigate lipid characteristics of each clusters
comp_annot_clust <- df %>%
  inner_join(cluster_tbl %>% select(Component, cluster), by = "Component") %>%
  mutate(cluster = factor(cluster, levels = sort(unique(cluster))))


prop_class <- comp_annot_clust %>%
  count(cluster, Class, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_size <- comp_annot_clust %>%
  count(cluster, Size, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_sat <- comp_annot_clust %>%
  count(cluster, Saturation, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


pdf("output/cluster_composition-LOGLIMMA-SNX14_Lysosomes-8Cluster.pdf", width = 6, height = 4)
ggplot(prop_class, aes(x = cluster, y = prop, fill = Class)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Lipid Class composition per cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
ggplot(prop_size, aes(x = cluster, y = prop, fill = Size)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Size composition per cluster") +
  theme_classic()
ggplot(prop_sat, aes(x = cluster, y = prop, fill = Saturation)) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components", title = "Saturation composition per cluster") +
  theme_classic()
dev.off()






comp_annot_clust_all <- comp_annot_clust %>%
  select(Component, Class, Size, Saturation, cluster) %>%
  bind_rows(comp_annot_all)

cluster_levels <- c("Global", sort(unique(cluster_tbl$cluster)))

comp_annot_clust_all <- comp_annot_clust_all %>%
  mutate(cluster = factor(cluster, levels = cluster_levels)) %>%
  mutate(
    cluster = ifelse(is.na(cluster), "Global", as.character(cluster)),
    cluster = factor(cluster, levels = c("Global", sort(unique(cluster[cluster != "Global"]))))
  )


prop_class <- comp_annot_clust_all %>%
  count(cluster, Class, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_size <- comp_annot_clust_all %>%
  count(cluster, Size, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_sat <- comp_annot_clust_all %>%
  count(cluster, Saturation, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

pdf("output/cluster_composition_andGlobal-LOGLIMMA-SNX14_Lysosomes-8Cluster.pdf", width = 6, height = 4)
ggplot(prop_class, aes(x = cluster, y = prop, fill = Class)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Lipid Class composition (clusters vs global)") +
  theme_classic()
ggplot(prop_size, aes(x = cluster, y = prop, fill = Size)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Size composition (clusters vs global)") +
  theme_classic()
ggplot(prop_sat, aes(x = cluster, y = prop, fill = Saturation)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Proportion of components",
       title = "Saturation composition (clusters vs global)") +
  theme_classic()
dev.off()





```


--> This version works but it is confusing to have both Lysosomes and Cell fraction separated.. Let's try integrate both values into a single one next.









# Testing - v2 - ratio lysosome / cell fractions

Let's try to integrate Lysosomes and Cells fraction by dividing lysosome/cell numbers for each sample. As positive control, we should see the most enriched in lysosome being 40:5 (ie. PA(40:5) or DAG(40:5))



```bash
conda activate deseq2
```



```R

# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("limma")
library("readxl")




# Import file
lip <- read_excel(
  path  = "input/For_Thomas_Lipidomics_tidy.xlsx",
  sheet = "tidy clean"
)

# Tidy files
sample_cols <- names(lip)[str_detect(names(lip), "^(Lysosomes|Cells)-")]

lip_tidy <- lip %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance",
    values_transform = list(abundance = as.numeric)
  ) %>%
  extract(
    col = sample,
    into = c("fraction", "condition", "genotype", "replicate"),
    regex = "^(Lysosomes|Cells)-([A-Za-z0-9]+)_([^&]+)&(\\d+)$",
    remove = FALSE
  ) 


############################################################
## 1) Normalize within Class (your current approach)
############################################################

lip_tidy_norm <- lip_tidy %>%
  group_by(sample, fraction, condition, genotype, replicate, Class) %>%
  mutate(
    Class_total     = sum(abundance, na.rm = TRUE),
    abundance_Class = ifelse(Class_total > 0, abundance / Class_total, NA_real_)
  ) %>%
  ungroup()


############################################################
## 2) Compute replicate-level log2 ratio (Lys/Cells), then mean ± SEM
############################################################

lip_ratio_rep <- lip_tidy_norm %>%
  # safety: collapse any duplicates within replicate/fraction
  group_by(Component, Class, Size, Saturation, condition, genotype, replicate, fraction) %>%
  summarise(abundance_Class = mean(abundance_Class, na.rm = TRUE), .groups = "drop") %>%
  filter(fraction %in% c("Lysosomes", "Cells")) %>%
  pivot_wider(names_from = fraction, values_from = abundance_Class) %>%
  mutate(
    ratio_Lys_over_Cells = log2((Lysosomes + 1) / (Cells + 1))
  )

lip_tidy_norm_ratio <- lip_ratio_rep %>%
  group_by(Component, Class, Size, Saturation, condition, genotype) %>%
  summarise(
    mean_ratio = mean(ratio_Lys_over_Cells, na.rm = TRUE),
    sd_ratio   = sd(ratio_Lys_over_Cells, na.rm = TRUE),
    n_rep      = sum(!is.na(ratio_Lys_over_Cells)),
    sem_ratio  = sd_ratio / sqrt(n_rep),
    ymin = mean_ratio - sem_ratio,
    ymax = mean_ratio + sem_ratio,
    .groups = "drop"
  )

# Quick view
lip_tidy_norm_ratio %>%
  select(Component, condition, genotype, mean_ratio, sem_ratio, n_rep) %>%
  arrange(desc(mean_ratio))


#####################################################################
## QC plot: Scram DMSO — TOP 20 and BOTTOM 20 with error bars
#####################################################################

ratio_ranked <- lip_tidy_norm_ratio %>%
  filter(genotype == "Scram", condition == "DMSO") %>%
  filter(is.finite(mean_ratio), is.finite(ymin), is.finite(ymax)) %>%
  arrange(desc(mean_ratio))


# Top 10 per lipid Class
top10_per_class <- ratio_ranked %>%
  group_by(Class) %>%
  slice_max(order_by = mean_ratio, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  # Create a unique per-class label so ordering works independently in each facet
  mutate(Component_class = paste0(Component, "___", Class))

# Build factor levels that are ordered high -> low within each Class
level_tbl <- top10_per_class %>%
  arrange(Class, desc(mean_ratio)) %>%
  distinct(Class, Component_class)

top10_per_class <- top10_per_class %>%
  mutate(
    Component_class = factor(Component_class, levels = level_tbl$Component_class)
  )

p_top10_class <- ggplot(top10_per_class, aes(x = Component_class, y = mean_ratio)) +
  geom_col() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.25, linewidth = 0.6) +
  coord_flip() +
  facet_wrap(~ Class, scales = "free_y") +
  scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +  # show only Component name
  labs(
    title = "Top 10 Components per Lipid Class\nLysosomes / Cells ratio (Scram, DMSO)",
    x = NULL,
    y = "log2((Lysosomes + 1) / (Cells + 1)) ± SEM"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

pdf("output/ratio_Lysosomes_over_Cells_top10_perClass_withSEM.pdf",
    width = 10, height = 6)
print(p_top10_class)
dev.off()






```


Lysosome and whole-cell fractions prepared and measured in separate experiments → raw intensities are not directly comparable:
- Ratio lysosome vs cell cannot be computed from raw values (we need a spike-in / common quantitative anchor).
- Available normalization is within lipid class, which informs redistribution within a class, not absolute lysosomal enrichment: Ratios based on class-normalized values reflect within-class composition changes only and must be interpreted as such.

--> Easier to focus on Lysosome fraction only; identify interesting lipid changes and check whehter they also occur in Cell fraction; if not, these are likely speficic of lysosome: even easier: forget the cell fraction!


