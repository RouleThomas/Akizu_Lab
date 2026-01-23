# Project

Samples and conditions:
- genotypes (3): WT, KD SNX13, KD SNX14
- condition (3): DMSO, LLOME, recovery

Lysosome fraction only.

Output values available:
- Raw abundance from MS
- Relative abundance (ie. each protein expressed as a fraction of the total proteome) 
- Relative abundance - empty vector (ie. empty vector value subtracted from relative abundance)


--> Let's start with Relative abundance.


# Testing - v1

File is `For Thomas_Proteomics.xlsx`; let's:
- open in R
- work with `Normalized to Sum` columns

Samples names are:
- DMSO Scram1
- DMSO siSNX13-1
- DMSO siSNX14-1
- DMSO empty pulldown1
- LLOMe scram1
- LLOMe siSNX13-1
- LLOMe siSNX14-1
- LLOMe EMPTY1
- 4hr Recovery Scram1
- Recovery siSNX13-1
- Recovery siSNX14-1
- Recovery EMPTY1

--> Let's manually tidy the Xcell file as `For Thomas_Proteomics_tidy.xlsx`: I simply rename all columns name by removing any special characters and space + rename all; columns with replicate name.




```bash
conda activate deseq2
```


```R
# Load packages
library("tidyverse")
library("readxl")
library("pheatmap")

# Import file
prot <- read_excel(
  path  = "input/For_Thomas_Proteomics_tidy.xlsx",
  sheet = "Proteins"
)

# Tidy files
sample_cols <- names(prot)[str_detect(names(prot), "^(Raw|Norm)-")]

prot_tidy <- prot %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance",
    values_transform = list(abundance = as.numeric)
  ) %>%
  extract(
    col = sample,
    into = c("type", "condition", "genotype", "replicate"),
    regex = "^(Raw|Norm)-([A-Za-z0-9]+)_([^&]+)&(\\d+)$",
    remove = FALSE
  )

## QC
prot_tidy %>% count(type, condition, genotype, replicate) %>% arrange(type, condition, genotype, replicate)


##############################
# HIGH Protein only ##########
##############################

# --- 1) filter: High only + normalized only ---
df <- prot_tidy %>%
  filter(`Protein_FDR` == "High",
         type == "Norm")

# --- 2) average replicates per genotype x condition (one value per cell in heatmap) ---
df_sum <- df %>%
  group_by(`Accession`, `GeneSymbol`, genotype, condition) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop")

# --- 3) enforce desired column order (GENOTYPE blocks, then time course) ---
genotype_levels  <- c("Scram", "siSNX13", "siSNX14", "empty")
condition_levels <- c("DMSO", "LLOME", "RECOVERY")

desired_col_order <- as.vector(
  unlist(lapply(genotype_levels, function(g)
    paste(g, condition_levels, sep = "_")
  ))
)

df_sum <- df_sum %>%
  mutate(
    genotype  = factor(genotype, levels = genotype_levels),
    condition = factor(condition, levels = condition_levels),
    col_id    = paste(as.character(genotype), as.character(condition), sep = "_")
  )

# --- 4) build matrix: rows=protein, cols=genotype_condition ---
mat_df <- df_sum %>%
  select(Accession, GeneSymbol, col_id, abundance) %>%
  pivot_wider(names_from = col_id, values_from = abundance)

# rownames: GeneSymbol__Accession (unique + readable)
mat <- mat_df %>%
  mutate(
    row_name = ifelse(
      is.na(GeneSymbol) | GeneSymbol == "",
      Accession,
      paste0(GeneSymbol, "__", Accession)
    )
  ) %>%
  select(-Accession, -GeneSymbol) %>%
  as.data.frame()

rownames(mat) <- mat$row_name
mat$row_name <- NULL

# keep only desired columns that exist, in the EXACT order you want
col_order_present <- intersect(desired_col_order, colnames(mat))
mat <- mat[, col_order_present, drop = FALSE]

# --- 5) log-transform (helps proteomics scale; NA stays NA) ---
mat_log <- log2(mat + 1)

# --- 6) order rows by Scram_DMSO high -> low ---
if ("Scram_DMSO" %in% colnames(mat_log)) {
  mat_log <- mat_log[order(mat_log[, "Scram_DMSO"], decreasing = TRUE, na.last = TRUE), ]
}

# --- 7) plot heatmap: NO clustering ---
pdf("output/pheatmap-Norm_High-allNoCluster.pdf", width = 4, height = 4)
pheatmap(
  mat_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  fontsize_col = 10
)
dev.off()






##################################################
# HIGH Protein AND detected (no NaN) only #########
##################################################

# --- 1) filter: High only + normalized only ---
df <- prot_tidy %>%
  filter(`Protein_FDR` == "High",
         type == "Norm")

# --- 2) average replicates per genotype x condition AND REMOVE any proteins with NaN values
df_sum <- df %>%
  group_by(`Accession`, `GeneSymbol`, genotype, condition) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Accession) %>%
  filter(!any(is.na(abundance))) %>%
  ungroup()

# --- 3) enforce desired column order (GENOTYPE blocks, then time course) ---
genotype_levels  <- c("Scram", "siSNX13", "siSNX14", "empty")
condition_levels <- c("DMSO", "LLOME", "RECOVERY")

desired_col_order <- as.vector(
  unlist(lapply(genotype_levels, function(g)
    paste(g, condition_levels, sep = "_")
  ))
)

df_sum <- df_sum %>%
  mutate(
    genotype  = factor(genotype, levels = genotype_levels),
    condition = factor(condition, levels = condition_levels),
    col_id    = paste(as.character(genotype), as.character(condition), sep = "_")
  )

# --- 4) build matrix: rows=protein, cols=genotype_condition ---
mat_df <- df_sum %>%
  select(Accession, GeneSymbol, col_id, abundance) %>%
  pivot_wider(names_from = col_id, values_from = abundance)

# rownames: GeneSymbol__Accession (unique + readable)
mat <- mat_df %>%
  mutate(
    row_name = ifelse(
      is.na(GeneSymbol) | GeneSymbol == "",
      Accession,
      paste0(GeneSymbol, "__", Accession)
    )
  ) %>%
  select(-Accession, -GeneSymbol) %>%
  as.data.frame()

rownames(mat) <- mat$row_name
mat$row_name <- NULL

# keep only desired columns that exist, in the EXACT order you want
col_order_present <- intersect(desired_col_order, colnames(mat))
mat <- mat[, col_order_present, drop = FALSE]

# --- 5) log-transform (helps proteomics scale; NA stays NA) ---
mat_log <- log2(mat + 1)

# --- 6) order rows by Scram_DMSO high -> low ---
if ("Scram_DMSO" %in% colnames(mat_log)) {
  mat_log <- mat_log[order(mat_log[, "Scram_DMSO"], decreasing = TRUE, na.last = TRUE), ]
}

# --- 7) plot heatmap: NO clustering ---
pdf("output/pheatmap-Norm_High_noNA-allNoCluster.pdf", width = 4, height = 4)
pheatmap(
  mat_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  fontsize_col = 10
)
dev.off()


# Clustering relative to scram


# helper to subtract Scram per condition
delta_mat <- cbind(
  SNX13_DMSO     = mat_log[,"siSNX13_DMSO"]     - mat_log[,"Scram_DMSO"],
  SNX13_LLOME    = mat_log[,"siSNX13_LLOME"]    - mat_log[,"Scram_LLOME"],
  SNX13_RECOVERY = mat_log[,"siSNX13_RECOVERY"] - mat_log[,"Scram_RECOVERY"],

  SNX14_DMSO     = mat_log[,"siSNX14_DMSO"]     - mat_log[,"Scram_DMSO"],
  SNX14_LLOME    = mat_log[,"siSNX14_LLOME"]    - mat_log[,"Scram_LLOME"],
  SNX14_RECOVERY = mat_log[,"siSNX14_RECOVERY"] - mat_log[,"Scram_RECOVERY"]
)

delta_scaled <- t(scale(t(delta_mat)))


pdf("output/pheatmap-Norm_High_noNA-DeltaScram-Cluster.pdf", width = 4, height = 4)
pheatmap(
  delta_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Genotype-specific deviation from Scram"
)
dev.off()

#--> Not sure looks great hard to interpret






##################################################
# HIGH Protein AND detected (no NaN) only ########
# Scram vs SNX13 #################################
##################################################

# --- 1) filter: High only + normalized only ---
df <- prot_tidy %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX13"))

# --- 2) average replicates per genotype x condition AND REMOVE any proteins with NaN values
df_sum <- df %>%
  group_by(`Accession`, `GeneSymbol`, genotype, condition) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Accession) %>%
  filter(!any(is.na(abundance))) %>%
  ungroup()

# --- 3) enforce desired column order (GENOTYPE blocks, then time course) ---
genotype_levels  <- c("Scram", "siSNX13", "siSNX14", "empty")
condition_levels <- c("DMSO", "LLOME", "RECOVERY")

desired_col_order <- as.vector(
  unlist(lapply(genotype_levels, function(g)
    paste(g, condition_levels, sep = "_")
  ))
)

df_sum <- df_sum %>%
  mutate(
    genotype  = factor(genotype, levels = genotype_levels),
    condition = factor(condition, levels = condition_levels),
    col_id    = paste(as.character(genotype), as.character(condition), sep = "_")
  )

# --- 4) build matrix: rows=protein, cols=genotype_condition ---
mat_df <- df_sum %>%
  select(Accession, GeneSymbol, col_id, abundance) %>%
  pivot_wider(names_from = col_id, values_from = abundance)

# rownames: GeneSymbol__Accession (unique + readable)
mat <- mat_df %>%
  mutate(
    row_name = ifelse(
      is.na(GeneSymbol) | GeneSymbol == "",
      Accession,
      paste0(GeneSymbol, "__", Accession)
    )
  ) %>%
  select(-Accession, -GeneSymbol) %>%
  as.data.frame()

rownames(mat) <- mat$row_name
mat$row_name <- NULL

# keep only desired columns that exist, in the EXACT order you want
col_order_present <- intersect(desired_col_order, colnames(mat))
mat <- mat[, col_order_present, drop = FALSE]

# --- 5) log-transform (helps proteomics scale; NA stays NA) ---
mat_log <- log10(mat + 1)

# --- 6) order rows by Scram_DMSO high -> low ---
if ("Scram_DMSO" %in% colnames(mat_log)) {
  mat_log <- mat_log[order(mat_log[, "Scram_DMSO"], decreasing = TRUE, na.last = TRUE), ]
}

# --- 7) plot heatmap: NO clustering ---
pdf("output/pheatmap-Norm_High_noNA-scram_SNX13_NoCluster.pdf", width = 4, height = 4)
pheatmap(
  mat_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  fontsize_col = 10,
  color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(101),
)
dev.off()

## Identify changes siSNX13 vs Scram
delta_mat <- cbind(
  DMSO     = mat_log[, "siSNX13_DMSO"]     - mat_log[, "Scram_DMSO"],
  LLOME    = mat_log[, "siSNX13_LLOME"]    - mat_log[, "Scram_LLOME"],
  RECOVERY = mat_log[, "siSNX13_RECOVERY"] - mat_log[, "Scram_RECOVERY"]
)
rownames(delta_mat) <- rownames(mat_log)
delta_scaled <- t(scale(t(delta_mat)))
hc <- hclust(dist(delta_scaled), method = "ward.D2")


## clustering
clusters <- cutree(hc, k = 4)  
mat_log_ord <- mat_log[hc$order, ]

cluster_df <- data.frame(
  Cluster = factor(clusters[hc$order])
)
rownames(cluster_df) <- rownames(mat_log_ord)

pdf("output/pheatmap-Norm_High_noNA-scram_SNX13_4Cluster.pdf", width = 3, height = 4)
pheatmap(
  mat_log_ord,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = cluster_df,
  show_rownames = FALSE,
  fontsize_col = 10,
  color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(101),
)
dev.off()


## Prioritize strong changes ##########
pattern_strength <- rowSums(abs(delta_mat))
### Combine clustering + strength into one table
cluster_df <- data.frame(
  protein = rownames(delta_mat),
  cluster = clusters,
  pattern_strength = pattern_strength
)
### Rank proteins within each cluster by strength
cluster_df_ranked <- cluster_df %>%
  group_by(cluster) %>%
  arrange(desc(pattern_strength)) %>%
  mutate(rank_in_cluster = row_number()) %>%
  ungroup()

ord <- order(cluster_df$cluster, -cluster_df$pattern_strength)
mat_log_ord2 <- mat_log[ord, ]

cluster_annot <- data.frame(
  Cluster = factor(cluster_df$cluster[ord])
)
rownames(cluster_annot) <- rownames(mat_log_ord2)


pdf("output/pheatmap-Norm_High_noNA-scram_SNX13_4Cluster-strength.pdf", width = 3, height = 4)
pheatmap(
  mat_log_ord2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = cluster_annot,
  show_rownames = FALSE,
  fontsize_col = 10,
  color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(101)
)
dev.off()






##################################################
# HIGH Protein AND detected (no NaN) only ########
# Scram vs SNX13  ################################
# Filter high changes - top 20% ############################
##################################################


pattern_strength <- rowSums(abs(delta_mat))

## Top 20% strongest changes
cutoff <- quantile(pattern_strength, 0.80)
keep <- pattern_strength >= cutoff


## Filter
delta_mat_filt    <- delta_mat[keep, ]
mat_log_filt      <- mat_log[keep, ]
pattern_strength_filt <- pattern_strength[keep]


delta_scaled_filt <- t(scale(t(delta_mat_filt)))
hc_filt <- hclust(dist(delta_scaled_filt), method = "ward.D2")
clusters_filt <- cutree(hc_filt, k = 10)


ord <- order(clusters_filt, -pattern_strength_filt)
mat_log_ord <- mat_log_filt[ord, ]


cluster_annot <- data.frame(
  Cluster = factor(clusters_filt[ord])
)
rownames(cluster_annot) <- rownames(mat_log_ord)



pdf("output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-Top20percent.pdf", width = 3, height = 4)
pheatmap(
  t(scale(t(mat_log_ord))),  # row-scaled DISPLAY
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = cluster_annot,
  show_rownames = FALSE
)
dev.off()


## Save gene list

cluster_tidy <- data.frame(
  row_name = rownames(mat_log_ord),
  cluster  = clusters_filt[ord],
  stringsAsFactors = FALSE
) %>%
  separate(
    row_name,
    into = c("GeneSymbol", "Accession"),
    sep = "__",
    remove = TRUE,
    fill = "right"
  ) %>%
  mutate(
    GeneSymbol = ifelse(is.na(GeneSymbol) | GeneSymbol == "",
                        Accession,
                        GeneSymbol)
  ) 



write.table(
  cluster_tidy,
  file = "output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-Top20percent.csv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)




## Plot top 10 proteins per cluster (to double check their profile in agreement)

# 1) Build a tidy table of cluster membership in your current ordering
#    (rownames(mat_log_ord) are "GeneSymbol__Accession")
# -----------------------------
cluster_map <- data.frame(
  row_name = rownames(mat_log_ord),
  cluster  = clusters_filt[ord],
  stringsAsFactors = FALSE
) %>%
  separate(row_name, into = c("GeneSymbol", "Accession"), sep = "__", fill = "right", remove = TRUE) %>%
  mutate(
    GeneSymbol = ifelse(is.na(GeneSymbol) | GeneSymbol == "", Accession, GeneSymbol)
  )

# -----------------------------
# 2) Pick top 10 proteins per cluster (based on your current order already)
# -----------------------------
top10_per_cluster <- cluster_map %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()

hits <- top10_per_cluster$Accession

# -----------------------------
# 3) Pull replicate-level data (NOT df_sum) and compute mean ± SEM
#    IMPORTANT: prot_tidy must have replicate column
# -----------------------------
df_rep <- prot_tidy %>%
  filter(
    Protein_FDR == "High",
    type == "Norm",
    genotype %in% c("Scram", "siSNX13"),
    Accession %in% hits
  ) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX13"))
  ) %>%
  left_join(top10_per_cluster %>% select(Accession, cluster),
            by = c("Accession"))

# mean + SEM across replicates
df_plot <- df_rep %>%
  group_by(cluster, Accession, GeneSymbol, genotype, condition) %>%
  summarise(
    mean_ab = mean(abundance, na.rm = TRUE),
    sd_ab   = sd(abundance, na.rm = TRUE),
    n       = sum(!is.na(abundance)),
    sem_ab  = sd_ab / sqrt(n),
    .groups = "drop"
  )

# -----------------------------
# 4) Plot: one PDF per cluster, facets = 10 proteins
# -----------------------------
dir.create("output/cluster_lineplots", showWarnings = FALSE, recursive = TRUE)

for (k in sort(unique(df_plot$cluster))) {

  df_k <- df_plot %>% filter(cluster == k) %>%
    mutate(panel = paste0(GeneSymbol, " (", Accession, ")"))

  p <- ggplot(df_k, aes(x = condition, y = mean_ab, group = genotype, color = genotype)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_ab - sem_ab, ymax = mean_ab + sem_ab),
                  width = 0.15, linewidth = 0.5) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    labs(
      title = paste0("Top 10 proteins — Cluster ", k, " (Scram vs siSNX13)"),
      x = NULL,
      y = "Normalized abundance (mean ± SEM)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 8)
    )

  ggsave(
    filename = sprintf("output/cluster_lineplots/Cluster_%02d_top10_timecourse_Scram_vs_SNX13-10Cluster-Top20percent.pdf", k),
    plot = p,
    width = 7,
    height = 6
  )
}






##################################################
# HIGH Protein AND detected (no NaN) only ########
# Scram vs SNX13  ################################
# Filter high changes - FC 0.58 ############################
##################################################


pattern_strength <- rowSums(abs(delta_mat))

## at least one timepoint with |log2FC| ≥ 0.58

keep <- apply(abs(delta_mat), 1, max) >= 0.58

## Filter
delta_mat_filt    <- delta_mat[keep, ]
mat_log_filt      <- mat_log[keep, ]
pattern_strength_filt <- pattern_strength[keep]


delta_scaled_filt <- t(scale(t(delta_mat_filt)))
hc_filt <- hclust(dist(delta_scaled_filt), method = "ward.D2")
clusters_filt <- cutree(hc_filt, k = 15)


ord <- order(clusters_filt, -pattern_strength_filt)
mat_log_ord <- mat_log_filt[ord, ]


cluster_annot <- data.frame(
  Cluster = factor(clusters_filt[ord])
)
rownames(cluster_annot) <- rownames(mat_log_ord)



pdf("output/pheatmap-Norm_High_noNA-scram_SNX13_15Cluster-FC058.pdf", width = 3, height = 4)
pheatmap(
  t(scale(t(mat_log_ord))),  # row-scaled DISPLAY
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = cluster_annot,
  show_rownames = FALSE
)
dev.off()

## Save gene list

cluster_tidy <- data.frame(
  row_name = rownames(mat_log_ord),
  cluster  = clusters_filt[ord],
  stringsAsFactors = FALSE
) %>%
  separate(
    row_name,
    into = c("GeneSymbol", "Accession"),
    sep = "__",
    remove = TRUE,
    fill = "right"
  ) %>%
  mutate(
    GeneSymbol = ifelse(is.na(GeneSymbol) | GeneSymbol == "",
                        Accession,
                        GeneSymbol)
  ) 

write.table(
  cluster_tidy,
  file = "output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-FC058.csv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## Plot top 10 proteins per cluster (to double check their profile in agreement)

# 1) Build a tidy table of cluster membership in your current ordering
#    (rownames(mat_log_ord) are "GeneSymbol__Accession")
# -----------------------------
cluster_map <- data.frame(
  row_name = rownames(mat_log_ord),
  cluster  = clusters_filt[ord],
  stringsAsFactors = FALSE
) %>%
  separate(row_name, into = c("GeneSymbol", "Accession"), sep = "__", fill = "right", remove = TRUE) %>%
  mutate(
    GeneSymbol = ifelse(is.na(GeneSymbol) | GeneSymbol == "", Accession, GeneSymbol)
  )

# -----------------------------
# 2) Pick top 10 proteins per cluster (based on your current order already)
# -----------------------------
top10_per_cluster <- cluster_map %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()

hits <- top10_per_cluster$Accession

# -----------------------------
# 3) Pull replicate-level data (NOT df_sum) and compute mean ± SEM
#    IMPORTANT: prot_tidy must have replicate column
# -----------------------------
df_rep <- prot_tidy %>%
  filter(
    Protein_FDR == "High",
    type == "Norm",
    genotype %in% c("Scram", "siSNX13"),
    Accession %in% hits
  ) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX13"))
  ) %>%
  left_join(top10_per_cluster %>% select(Accession, cluster),
            by = c("Accession"))

# mean + SEM across replicates
df_plot <- df_rep %>%
  group_by(cluster, Accession, GeneSymbol, genotype, condition) %>%
  summarise(
    mean_ab = mean(abundance, na.rm = TRUE),
    sd_ab   = sd(abundance, na.rm = TRUE),
    n       = sum(!is.na(abundance)),
    sem_ab  = sd_ab / sqrt(n),
    .groups = "drop"
  )

# -----------------------------
# 4) Plot: one PDF per cluster, facets = 10 proteins
# -----------------------------
dir.create("output/cluster_lineplots", showWarnings = FALSE, recursive = TRUE)

for (k in sort(unique(df_plot$cluster))) {

  df_k <- df_plot %>% filter(cluster == k) %>%
    mutate(panel = paste0(GeneSymbol, " (", Accession, ")"))

  p <- ggplot(df_k, aes(x = condition, y = mean_ab, group = genotype, color = genotype)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_ab - sem_ab, ymax = mean_ab + sem_ab),
                  width = 0.15, linewidth = 0.5) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    labs(
      title = paste0("Top 10 proteins — Cluster ", k, " (Scram vs siSNX13)"),
      x = NULL,
      y = "Normalized abundance (mean ± SEM)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 8)
    )

  ggsave(
    filename = sprintf("output/cluster_lineplots/Cluster_%02d_top10_timecourse_Scram_vs_SNX13-10Cluster-FC058.pdf", k),
    plot = p,
    width = 7,
    height = 6
  )
}


## Test abundance values directly ##
mat_disp <- mat_log_ord
# cap by global quantiles (makes heatmap readable without changing ordering)
lo <- quantile(mat_disp, 0.02, na.rm = TRUE)
hi <- quantile(mat_disp, 0.98, na.rm = TRUE)
mat_disp <- pmin(pmax(mat_disp, lo), hi)
pdf("output/pheatmap-Norm_High_noNA-scram_SNX13_6Cluster-FC058_Data.pdf", width = 3.5, height = 4.5)
pheatmap(
  mat_disp,                 # <-- RAW log2 values
  cluster_rows = FALSE,      # keep our custom ordering
  cluster_cols = FALSE,
  annotation_row = cluster_annot,
  show_rownames = FALSE,
  fontsize_col = 9
)
dev.off()
#--> NOT great...










##################################################
# HIGH Protein AND detected (no NaN) only ########
# Scram vs SNX13  ################################
# Filter high changes - FC 1 ############################
##################################################


pattern_strength <- rowSums(abs(delta_mat))

## at least one timepoint with |log2FC| ≥ 1

keep <- apply(abs(delta_mat), 1, max) >= 1

## Filter
delta_mat_filt    <- delta_mat[keep, ]
mat_log_filt      <- mat_log[keep, ]
pattern_strength_filt <- pattern_strength[keep]


delta_scaled_filt <- t(scale(t(delta_mat_filt)))
hc_filt <- hclust(dist(delta_scaled_filt), method = "ward.D2")
clusters_filt <- cutree(hc_filt, k = 10)


ord <- order(clusters_filt, -pattern_strength_filt)
mat_log_ord <- mat_log_filt[ord, ]


cluster_annot <- data.frame(
  Cluster = factor(clusters_filt[ord])
)
rownames(cluster_annot) <- rownames(mat_log_ord)



pdf("output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-FC1.pdf", width = 3, height = 4)
pheatmap(
  t(scale(t(mat_log_ord))),  # row-scaled DISPLAY
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = cluster_annot,
  show_rownames = FALSE
)
dev.off()


## Save gene list

cluster_tidy <- data.frame(
  row_name = rownames(mat_log_ord),
  cluster  = clusters_filt[ord],
  stringsAsFactors = FALSE
) %>%
  separate(
    row_name,
    into = c("GeneSymbol", "Accession"),
    sep = "__",
    remove = TRUE,
    fill = "right"
  ) %>%
  mutate(
    GeneSymbol = ifelse(is.na(GeneSymbol) | GeneSymbol == "",
                        Accession,
                        GeneSymbol)
  ) 


write.table(
  cluster_tidy,
  file = "output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-FC1.csv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#--> NOT enough proteins with FC 1


```




Best so far:
- Only keep `Accession == High` and **remove any Accession if values is NaN** in one condition
- Work with Scram vs SNX13; or Scram vs SNX14 and filter to keep delta FC > 0.58



Pipeline:
- Used normalized-to-sum proteomics data
- Proteins with missing values in any condition were removed.
- Only proteins annotated with High protein FDR confidence were retained.
- Abundance values were log-transformed (log2(x + 1)).
- Proteins showing a genotype-dependent change were selected by requiring an ABS(L2FC) ≥ 0.58 in at least one condition (siSNX13 relative to Scram).
- Selected proteins were clustered based on the temporal pattern of their genotype-specific log2 fold-changes.
- For visualization, protein abundance values were row-scaled (Z-scored) to emphasize relative expression patterns across conditions.


--> Generate plot of unique genes to check if that respect the heamtap: yes, kind of.. Maybe increase clustering





# Testing - v2 - LOGLIMMA

Let's try to use bulk RNAseq based method for this.. From [this](https://www.reddit.com/r/bioinformatics/comments/12umj5l/mass_spectrometry_raw_count_data_and_analysis_a/) forum, it seems people used limma DE for protein...


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
prot <- read_excel(
  path  = "input/For_Thomas_Proteomics_tidy.xlsx",
  sheet = "Proteins"
)

# Tidy files
sample_cols <- names(prot)[str_detect(names(prot), "^(Raw|Norm)-")]

prot_tidy <- prot %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance",
    values_transform = list(abundance = as.numeric)
  ) %>%
  extract(
    col = sample,
    into = c("type", "condition", "genotype", "replicate"),
    regex = "^(Raw|Norm)-([A-Za-z0-9]+)_([^&]+)&(\\d+)$",
    remove = FALSE
  ) 


########################################
## All samples PCA ####################
########################################


df <- prot_tidy %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX13", "siSNX14")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)




## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)
## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))
## 5) Strict filter: remove any protein with ≥1 NA
keep <- complete.cases(log2_expr)
log2_expr_f <- log2_expr[keep, , drop = FALSE]

pca <- prcomp(t(log2_expr_f), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-allSamples.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#



########################################
## Mutants DMSO PCA ####################
########################################


df <- prot_tidy %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("siSNX13", "siSNX14"),
         condition %in% c("DMSO")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)




## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)
## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))
## 5) Strict filter: remove any protein with ≥1 NA
keep <- complete.cases(log2_expr)
log2_expr_f <- log2_expr[keep, , drop = FALSE]

pca <- prcomp(t(log2_expr_f), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-MutatntsDMSO.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...









########################################
## Scram vs siSNX13 ####################
########################################

df <- prot_tidy %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX13")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
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

## 5) Strict filter: remove any protein with ≥1 NA
keep <- complete.cases(log2_expr)
log2_expr_f <- log2_expr[keep, , drop = FALSE]

## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr_f, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX13.conditionLLOME),
      abs(genotypesiSNX13.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 1)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Accession

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 1055




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




k <- 6
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")

pdf("output/pheatmap-LOGLIMMA-SNX13-6Cluster.pdf", width = 4, height = 4)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = cumsum(table(sample_order$genotype, sample_order$condition)[1,]),
         main = "Interaction proteins – replicate-aware clustering")
dev.off()
#--> Some replicate seems outliers...




pca <- prcomp(t(log2_expr_sig), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-SNX13.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...






############################################################
## Scram vs siSNX13 - without outlier (siSNX13 DMSO Rep1) ####################
############################################################

df <- prot_tidy %>%
  filter(!(condition == "DMSO" & genotype == "siSNX13" & replicate == "1")) %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX13")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
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

## 5) Strict filter: remove any protein with ≥1 NA
keep <- complete.cases(log2_expr)
log2_expr_f <- log2_expr[keep, , drop = FALSE]

## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr_f, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX13.conditionLLOME),
      abs(genotypesiSNX13.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 1)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Accession

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 1055




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




k <- 6
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")

pdf("output/pheatmap-LOGLIMMA-SNX13filter-6Cluster.pdf", width = 4, height = 4)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = cumsum(table(sample_order$genotype, sample_order$condition)[1,]),
         main = "Interaction proteins – replicate-aware clustering")
dev.off()
#--> Some replicate seems outliers...




pca <- prcomp(t(log2_expr_sig), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-SNX13filter-6Cluster.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...








########################################
## Scram vs siSNX14 ####################
########################################

df <- prot_tidy %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX14")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
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

## 5) Strict filter: remove any protein with ≥1 NA
keep <- complete.cases(log2_expr)
log2_expr_f <- log2_expr[keep, , drop = FALSE]

## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr_f, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX14.conditionLLOME),
      abs(genotypesiSNX14.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 1)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Accession

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   



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




k <- 6
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")

pdf("output/pheatmap-LOGLIMMA-SNX14-6Cluster.pdf", width = 4, height = 4)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = cumsum(table(sample_order$genotype, sample_order$condition)[1,]),
         main = "Interaction proteins – replicate-aware clustering")
dev.off()
#--> Some replicate seems outliers...

pca <- prcomp(t(log2_expr_sig), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-SNX14-6Cluster.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...





############################################################
## Scram vs siSNX14 - without outlier (siSNX14 DMSO Rep3) ####################
############################################################

df <- prot_tidy %>%
  filter(!(condition == "DMSO" & genotype == "siSNX14" & replicate == "3")) %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX14")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
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

## 5) Strict filter: remove any protein with ≥1 NA
keep <- complete.cases(log2_expr)
log2_expr_f <- log2_expr[keep, , drop = FALSE]

## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr_f, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX14.conditionLLOME),
      abs(genotypesiSNX14.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 1)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Accession

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 1055




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




k <- 6
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")

pdf("output/pheatmap-LOGLIMMA-SNX14filter-6Cluster.pdf", width = 4, height = 4)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = cumsum(table(sample_order$genotype, sample_order$condition)[1,]),
         main = "Interaction proteins – replicate-aware clustering")
dev.off()
#--> Some replicate seems outliers...




pca <- prcomp(t(log2_expr_sig), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  coldata
)
pdf("output/pca-LOGLIMMA-SNX14filter-6Cluster.pdf", width = 4, height = 4)
ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = replicate), vjust = -1) +
  theme_bw()
dev.off()
#--> PCA is not too bad...




```


--> I tested using Raw (non sum normalized) and it was bad: So **sum normalized should be used.** However looking at signal overall I noticed **two outliers samples**: 
- siSNX13 DMSO Rep1
- siSNX14 DMSO Rep3

I think better to remove them for the analysis...





# LOGLIMMA

Let's follow `# Testing - v2 - LOGLIMMA` with the log limma method:
- Log transform
- Construct linear model with interaction (~genotype + condition)
- Filter signficant changes (FC > 1)
- Scale signal per row; per accession
- Unbiased clustering 
- Heatmap vizualization

--> Let's remove the outliers samples


In this version any proteins with 1 NA value in one sample was removed! 


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
prot <- read_excel(
  path  = "input/For_Thomas_Proteomics_tidy.xlsx",
  sheet = "Proteins"
)

# Tidy abundance output file
sample_cols <- names(prot)[str_detect(names(prot), "^(Raw|Norm)-")]

prot_tidy <- prot %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance",
    values_transform = list(abundance = as.numeric)
  ) %>%
  extract(
    col = sample,
    into = c("type", "condition", "genotype", "replicate"),
    regex = "^(Raw|Norm)-([A-Za-z0-9]+)_([^&]+)&(\\d+)$",
    remove = FALSE
  ) 


############################################################
## Scram vs siSNX13 - without outlier (siSNX13 DMSO Rep1) ####################
############################################################

df <- prot_tidy %>%
  filter(!(condition == "DMSO" & genotype == "siSNX13" & replicate == "1")) %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX13")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
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

## 5) Strict filter: remove any protein with ≥1 NA
keep <- complete.cases(log2_expr)
log2_expr_f <- log2_expr[keep, , drop = FALSE]


####################################################
## Export the proteins I removed: ####
removed_accessions <- rownames(log2_expr)[!keep]
# Keep ALL columns/rows from df for those proteins
df_removed <- df %>%
  filter(Accession %in% removed_accessions) %>%
  arrange(Accession, genotype, condition, replicate, sample)
write.table(
  df_removed,
  file = "output/removed_proteins_rows_with_NA-SNX13filter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
######################################################


## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr_f, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX13.conditionLLOME),
      abs(genotypesiSNX13.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 1)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Accession

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 1055




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




k <- 12
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")


n_rep <- sample_order %>%
  count(genotype, condition)
block_sizes <- n_rep$n


row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-SNX13filter-12Cluster.pdf", width = 4, height = 5)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = gaps_col,
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()

out_tbl <- cluster_tbl %>%
  select(Accession, GeneSymbol, cluster) %>%
  arrange(cluster, Accession)

write_tsv(out_tbl, "output/GENELIST-LOGLIMMA-SNX13filter-12Cluster.tsv")



# Show top 10 genes per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Accession, GeneSymbol, cluster), by = "Accession") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX13")),
    log2_abund = log2(abundance + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Accession, GeneSymbol.x) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Accession", "GeneSymbol.x")) %>%
  group_by(cluster, Accession, GeneSymbol.x, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(GeneSymbol.x, " (", Accession, ")")
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
      title = paste0("Top 10 proteins — Cluster ", cl, " (Scram vs siSNX13)"),
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
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX13filter-12Cluster.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}







############################################################
## Scram vs siSNX14 - without outlier (siSNX14 DMSO Rep3) ####################
############################################################

df <- prot_tidy %>%
  filter(!(condition == "DMSO" & genotype == "siSNX14" & replicate == "3")) %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX14")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
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

## 5) Strict filter: remove any protein with ≥1 NA
keep <- complete.cases(log2_expr)
log2_expr_f <- log2_expr[keep, , drop = FALSE]




####################################################
## Export the proteins I removed: ####
removed_accessions <- rownames(log2_expr)[!keep]
# Keep ALL columns/rows from df for those proteins
df_removed <- df %>%
  filter(Accession %in% removed_accessions) %>%
  arrange(Accession, genotype, condition, replicate, sample)
write.table(
  df_removed,
  file = "output/removed_proteins_rows_with_NA-SNX14filter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
######################################################




## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr_f, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  



interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX14.conditionLLOME),
      abs(genotypesiSNX14.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 1)   # <-- threshold here


sig_acc <- interaction_tbl_filt$Accession

log2_expr_sig <- log2_expr[rownames(log2_expr) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_sig)   # should be 1055




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




k <- 10
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")


n_rep <- sample_order %>%
  count(genotype, condition)
block_sizes <- n_rep$n


row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-SNX14filter-10Cluster.pdf", width = 4, height = 5)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = gaps_col,
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()

out_tbl <- cluster_tbl %>%
  select(Accession, GeneSymbol, cluster) %>%
  arrange(cluster, Accession)

write_tsv(out_tbl, "output/GENELIST-LOGLIMMA-SNX14filter-10Cluster.tsv")




# Show top 10 genes per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Accession, GeneSymbol, cluster), by = "Accession") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX14")),
    log2_abund = log2(abundance + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Accession, GeneSymbol.x) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Accession", "GeneSymbol.x")) %>%
  group_by(cluster, Accession, GeneSymbol.x, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(GeneSymbol.x, " (", Accession, ")")
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
      title = paste0("Top 10 proteins — Cluster ", cl, " (Scram vs siSNX14)"),
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
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX14filter-10Cluster.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}

```


--> Seems to work great, individual gene plot in agreement with the heamtap
  --> Although filtering NA value proteins might be a bit rough, see next version





  
# LOGLIMMA - NA values filtering - v1

Let's follow previous version `# LOGLIMMA`; just improve the filtering of NA values:
- Keep proteins that have at least 2 replicates in at least 1 condition and assing 0 (or an arbitrary number close to 0) to the NA values in other conditions



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
prot <- read_excel(
  path  = "input/For_Thomas_Proteomics_tidy.xlsx",
  sheet = "Proteins"
)

# Tidy abundance output file
sample_cols <- names(prot)[str_detect(names(prot), "^(Raw|Norm)-")]

prot_tidy <- prot %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance",
    values_transform = list(abundance = as.numeric)
  ) %>%
  extract(
    col = sample,
    into = c("type", "condition", "genotype", "replicate"),
    regex = "^(Raw|Norm)-([A-Za-z0-9]+)_([^&]+)&(\\d+)$",
    remove = FALSE
  ) 


############################################################
## Scram vs siSNX13 - without outlier (siSNX13 DMSO Rep1) - ####################
############################################################

df <- prot_tidy %>%
  filter(!(condition == "DMSO" & genotype == "siSNX13" & replicate == "1")) %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX13")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)


## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))

stopifnot(all(coldata$sample == colnames(log2_expr)))

coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = c("DMSO","LLOME","RECOVERY"))
coldata$genotype  <- factor(coldata$genotype,  levels = c("Scram","siSNX13"))
coldata$replicate <- factor(coldata$replicate)

## log2_expr is a matrix with rownames = Accession, colnames = sample
## coldata has rownames = sample, and columns condition, genotype, replicate

group <- interaction(coldata$condition, coldata$genotype, drop = TRUE)

# For each protein, count non-NA per group
ok_by_group <- sapply(levels(group), function(g) {
  cols <- rownames(coldata)[group == g]
  rowSums(!is.na(log2_expr[, cols, drop = FALSE])) >= 2
})

keep <- apply(ok_by_group, 1, any)
log2_expr_keep <- log2_expr[keep, , drop = FALSE]


log2_expr_imp <- log2_expr_keep

for (j in seq_len(ncol(log2_expr_imp))) {
  x <- log2_expr_imp[, j]
  min_x <- min(x, na.rm = TRUE)
  # put missing slightly below the observed minimum
  x[is.na(x)] <- min_x - 1
  log2_expr_imp[, j] <- x
}



# Table that say which Accession got NA values and whtehr it has been change with another value
imp_info <- tibble(
  Accession = rownames(log2_expr_keep),
  n_imputed = rowSums(is.na(log2_expr_keep)),
  any_imputed = n_imputed > 0
)









####################################################
## Export the proteins I removed (and why)
####################################################
# 1) Removed accessions
removed_accessions <- rownames(log2_expr)[!keep]
# 2) Per-protein non-NA counts per (condition x genotype) group
group <- interaction(coldata$condition, coldata$genotype, drop = TRUE)
non_na_counts <- sapply(levels(group), function(g) {
  cols <- rownames(coldata)[group == g]
  rowSums(!is.na(log2_expr[, cols, drop = FALSE]))
})
non_na_counts <- as.data.frame(non_na_counts)
non_na_counts$Accession <- rownames(log2_expr)
# max replicates observed in ANY group (your keep criterion was >=2)
non_na_counts$max_non_na_any_group <- apply(non_na_counts[levels(group)], 1, max)
# 3) Build a small “removed summary” table
removed_summary <- non_na_counts %>%
  dplyr::filter(Accession %in% removed_accessions) %>%
  dplyr::arrange(desc(max_non_na_any_group), Accession)
write.table(
  removed_summary,
  file = "output/removed_proteins_summary-SNX13filterNAvaluesFilter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# 4) Export the full long df rows for removed proteins (all samples/conditions)
df_removed <- df %>%
  dplyr::filter(Accession %in% removed_accessions) %>%
  dplyr::arrange(Accession, genotype, condition, replicate, sample)
write.table(
  df_removed,
  file = "output/removed_proteins_longrows-SNX13filterNAvaluesFilter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
####################################################





## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr_imp, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  

############################
## max_abs_interaction >= 2  ####

interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX13.conditionLLOME),
      abs(genotypesiSNX13.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 2)   # <-- threshold here adj.P.Val < 0.05, max_abs_interaction >= 1


sig_acc <- interaction_tbl_filt$Accession

log2_expr_imp_sig <- log2_expr_imp[rownames(log2_expr_imp) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_imp_sig)   # should be 3429/1815/1023




# Build ordered sample table
sample_order <- coldata %>%
  mutate(
    condition = factor(condition, levels = c("DMSO","LLOME","RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram","siSNX13"))
  ) %>%
  arrange(genotype, condition, replicate)

# Reorder matrix
log2_expr_imp_sig_ord <- log2_expr_imp_sig[, sample_order$sample]




expr_scaled <- t(scale(t(log2_expr_imp_sig_ord)))

row_dist   <- dist(expr_scaled, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")




k <- 7
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")


n_rep <- sample_order %>%
  count(genotype, condition)
block_sizes <- n_rep$n


row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-SNX13filterNAvaluesFilter-7Cluster-max_abs_interaction2.pdf", width = 4, height = 5)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = head(cumsum(sample_order %>% count(genotype, condition) %>% arrange(genotype, condition) %>% pull(n)), -1),
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()


cluster_tbl <- cluster_tbl %>%
  left_join(imp_info, by = "Accession") %>%
  mutate(
    n_imputed = replace_na(n_imputed, 0L),
    any_imputed = replace_na(any_imputed, FALSE)
  )
  
write_tsv(cluster_tbl, "output/GENELIST-LOGLIMMA-SNX13filterNAvaluesFilter-7Cluster-max_abs_interaction2.tsv")



# Show top 10 genes per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Accession, GeneSymbol, cluster), by = "Accession") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX13")),
    log2_abund = log2(abundance + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Accession, GeneSymbol.x) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Accession", "GeneSymbol.x")) %>%
  group_by(cluster, Accession, GeneSymbol.x, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(GeneSymbol.x, " (", Accession, ")")
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
      title = paste0("Top 10 proteins — Cluster ", cl, " (Scram vs siSNX13)"),
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
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX13filterNAvaluesFilter-7Cluster-max_abs_interaction2.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}
































############################################################
## Scram vs siSNX14 - without outlier (siSNX14 DMSO Rep3) ####################
############################################################

df <- prot_tidy %>%
  filter(!(condition == "DMSO" & genotype == "siSNX14" & replicate == "3")) %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX14")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)


## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))

stopifnot(all(coldata$sample == colnames(log2_expr)))

coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = c("DMSO","LLOME","RECOVERY"))
coldata$genotype  <- factor(coldata$genotype,  levels = c("Scram","siSNX14"))
coldata$replicate <- factor(coldata$replicate)

## log2_expr is a matrix with rownames = Accession, colnames = sample
## coldata has rownames = sample, and columns condition, genotype, replicate

group <- interaction(coldata$condition, coldata$genotype, drop = TRUE)

# For each protein, count non-NA per group
ok_by_group <- sapply(levels(group), function(g) {
  cols <- rownames(coldata)[group == g]
  rowSums(!is.na(log2_expr[, cols, drop = FALSE])) >= 2
})

keep <- apply(ok_by_group, 1, any)
log2_expr_keep <- log2_expr[keep, , drop = FALSE]


log2_expr_imp <- log2_expr_keep

for (j in seq_len(ncol(log2_expr_imp))) {
  x <- log2_expr_imp[, j]
  min_x <- min(x, na.rm = TRUE)
  # put missing slightly below the observed minimum
  x[is.na(x)] <- min_x - 1
  log2_expr_imp[, j] <- x
}



# Table that say which Accession got NA values and whtehr it has been change with another value
imp_info <- tibble(
  Accession = rownames(log2_expr_keep),
  n_imputed = rowSums(is.na(log2_expr_keep)),
  any_imputed = n_imputed > 0
)









####################################################
## Export the proteins I removed (and why)
####################################################
# 1) Removed accessions
removed_accessions <- rownames(log2_expr)[!keep]
# 2) Per-protein non-NA counts per (condition x genotype) group
group <- interaction(coldata$condition, coldata$genotype, drop = TRUE)
non_na_counts <- sapply(levels(group), function(g) {
  cols <- rownames(coldata)[group == g]
  rowSums(!is.na(log2_expr[, cols, drop = FALSE]))
})
non_na_counts <- as.data.frame(non_na_counts)
non_na_counts$Accession <- rownames(log2_expr)
# max replicates observed in ANY group (your keep criterion was >=2)
non_na_counts$max_non_na_any_group <- apply(non_na_counts[levels(group)], 1, max)
# 3) Build a small “removed summary” table
removed_summary <- non_na_counts %>%
  dplyr::filter(Accession %in% removed_accessions) %>%
  dplyr::arrange(desc(max_non_na_any_group), Accession)
write.table(
  removed_summary,
  file = "output/removed_proteins_summary-SNX14filterNAvaluesFilter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# 4) Export the full long df rows for removed proteins (all samples/conditions)
df_removed <- df %>%
  dplyr::filter(Accession %in% removed_accessions) %>%
  dplyr::arrange(Accession, genotype, condition, replicate, sample)
write.table(
  df_removed,
  file = "output/removed_proteins_longrows-SNX14filterNAvaluesFilter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
####################################################





## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata)

fit <- lmFit(log2_expr_imp, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  

############################
## max_abs_interaction >= 2  ####

interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX14.conditionLLOME),
      abs(genotypesiSNX14.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 2)   # <-- threshold here adj.P.Val < 0.05, max_abs_interaction >= 1


sig_acc <- interaction_tbl_filt$Accession

log2_expr_imp_sig <- log2_expr_imp[rownames(log2_expr_imp) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_imp_sig)   # should be 3205/1959/1129




# Build ordered sample table
sample_order <- coldata %>%
  mutate(
    condition = factor(condition, levels = c("DMSO","LLOME","RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram","siSNX13"))
  ) %>%
  arrange(genotype, condition, replicate)

# Reorder matrix
log2_expr_imp_sig_ord <- log2_expr_imp_sig[, sample_order$sample]




expr_scaled <- t(scale(t(log2_expr_imp_sig_ord)))

row_dist   <- dist(expr_scaled, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")




k <- 8
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")


n_rep <- sample_order %>%
  count(genotype, condition)
block_sizes <- n_rep$n


row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-SNX14filterNAvaluesFilter-8Cluster-max_abs_interaction2.pdf", width = 4, height = 5)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = head(cumsum(sample_order %>% count(genotype, condition) %>% arrange(genotype, condition) %>% pull(n)), -1),
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()


cluster_tbl <- cluster_tbl %>%
  left_join(imp_info, by = "Accession") %>%
  mutate(
    n_imputed = replace_na(n_imputed, 0L),
    any_imputed = replace_na(any_imputed, FALSE)
  )
  
write_tsv(cluster_tbl, "output/GENELIST-LOGLIMMA-SNX14filterNAvaluesFilter-8Cluster-max_abs_interaction2.tsv")



# Show top 10 genes per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Accession, GeneSymbol, cluster), by = "Accession") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX14")),
    log2_abund = log2(abundance + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Accession, GeneSymbol.x) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Accession", "GeneSymbol.x")) %>%
  group_by(cluster, Accession, GeneSymbol.x, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(GeneSymbol.x, " (", Accession, ")")
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
      title = paste0("Top 10 proteins — Cluster ", cl, " (Scram vs siSNX14)"),
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
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX14filterNAvaluesFilter-8Cluster-max_abs_interaction2.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}





```


--> GOOD! But lets perform the filtering/updating of NA values directly with all genotypes and condition, to avoid having a different set of proteins compared between each of my different comparions (Scramble vs siSNX13, Scramble vs siSNX14)


















































# LOGLIMMA - NA values filtering - v2

- remove outlier samples
- Filter the NA values and generate a new table with normalize abundance
- Filter sample of interest for comparison
  - Scramble vs siSNX13
  - Scramble vs siSNX14
  - DMSO: Scramble vs siSNX13 vs siSNX14




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
prot <- read_excel(
  path  = "input/For_Thomas_Proteomics_tidy.xlsx",
  sheet = "Proteins"
)

# Tidy abundance output file
sample_cols <- names(prot)[str_detect(names(prot), "^(Raw|Norm)-")]

prot_tidy <- prot %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance",
    values_transform = list(abundance = as.numeric)
  ) %>%
  extract(
    col = sample,
    into = c("type", "condition", "genotype", "replicate"),
    regex = "^(Raw|Norm)-([A-Za-z0-9]+)_([^&]+)&(\\d+)$",
    remove = FALSE
  ) 


####################################################################################################
## Filter outlier samples (siSNX13 DMSO Rep1, siSNX14 DMSO Rep3)  #######
####################################################################################################

df <- prot_tidy %>%
  filter(!(condition == "DMSO" & genotype == "siSNX13" & replicate == "1")) %>%
  filter(!(condition == "DMSO" & genotype == "siSNX14" & replicate == "3")) %>%
  filter(`Protein_FDR` == "High",
         type == "Norm",
         genotype %in% c("Scram", "siSNX13", "siSNX14")) %>%
  dplyr::select(Accession, GeneSymbol, sample, condition, genotype, replicate, abundance)


## 2) Expression matrix: proteins (rows) x samples (cols)
expr_mat <- df %>%
  select(Accession, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("Accession") %>%
  as.matrix()


log2_expr <- log2(expr_mat + 1)


## 4) Sample metadata (like your colData)
coldata <- df %>%
  distinct(sample, condition, genotype, replicate) %>%
  arrange(match(sample, colnames(log2_expr)))

stopifnot(all(coldata$sample == colnames(log2_expr)))

coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = c("DMSO","LLOME","RECOVERY"))
coldata$genotype  <- factor(coldata$genotype,  levels = c("Scram","siSNX13","siSNX14"))
coldata$replicate <- factor(coldata$replicate)

## log2_expr is a matrix with rownames = Accession, colnames = sample
## coldata has rownames = sample, and columns condition, genotype, replicate


####################################################################################################
## Filter and/or transform NA values  #######
####################################################################################################


group <- interaction(coldata$condition, coldata$genotype, drop = TRUE)

# For each protein, count non-NA per group
ok_by_group <- sapply(levels(group), function(g) {
  cols <- rownames(coldata)[group == g]
  rowSums(!is.na(log2_expr[, cols, drop = FALSE])) >= 2
})

keep <- apply(ok_by_group, 1, any)
log2_expr_keep <- log2_expr[keep, , drop = FALSE]


log2_expr_imp <- log2_expr_keep

for (j in seq_len(ncol(log2_expr_imp))) {
  x <- log2_expr_imp[, j]
  min_x <- min(x, na.rm = TRUE)
  # put missing slightly below the observed minimum
  x[is.na(x)] <- min_x - 1
  log2_expr_imp[, j] <- x
}


# Table that say which Accession got NA values and whtehr it has been change with another value
imp_info <- tibble(
  Accession = rownames(log2_expr_keep),
  n_imputed = rowSums(is.na(log2_expr_keep)),
  any_imputed = n_imputed > 0
)




####################################################
## Export the proteins I removed (and why)
####################################################
# 1) Removed accessions
removed_accessions <- rownames(log2_expr)[!keep]
# 2) Per-protein non-NA counts per (condition x genotype) group
group <- interaction(coldata$condition, coldata$genotype, drop = TRUE)
non_na_counts <- sapply(levels(group), function(g) {
  cols <- rownames(coldata)[group == g]
  rowSums(!is.na(log2_expr[, cols, drop = FALSE]))
})
non_na_counts <- as.data.frame(non_na_counts)
non_na_counts$Accession <- rownames(log2_expr)
# max replicates observed in ANY group (your keep criterion was >=2)
non_na_counts$max_non_na_any_group <- apply(non_na_counts[levels(group)], 1, max)
# 3) Build a small “removed summary” table
removed_summary <- non_na_counts %>%
  dplyr::filter(Accession %in% removed_accessions) %>%
  dplyr::arrange(desc(max_non_na_any_group), Accession)
write.table(
  removed_summary,
  file = "output/removed_proteins_summary-AllNAvaluesFilter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# 4) Export the full long df rows for removed proteins (all samples/conditions)
df_removed <- df %>%
  dplyr::filter(Accession %in% removed_accessions) %>%
  dplyr::arrange(Accession, genotype, condition, replicate, sample)
write.table(
  df_removed,
  file = "output/removed_proteins_longrows-AllNAvaluesFilter.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
####################################################






####################################################
## Scramble vs siSNX13 #############################
####################################################

coldata__Scramble_vs_siSNX13 = coldata %>%
  filter(genotype %in% c("Scram", "siSNX13")) %>%
  mutate(
    genotype  = factor(genotype, levels = c("Scram","siSNX13")),
    condition = factor(condition, levels = c("DMSO","LLOME","RECOVERY"))
  )

log2_expr_imp__Scramble_vs_siSNX13 <- log2_expr_imp[
  , rownames(coldata__Scramble_vs_siSNX13), drop = FALSE
]





## 6) Model: full model = genotype + condition + genotype:condition
design_full <- model.matrix(~ genotype * condition, data = coldata__Scramble_vs_siSNX13)

fit <- lmFit(log2_expr_imp__Scramble_vs_siSNX13, design_full)
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
#--> Signficant Accession = different between genotype during the time-course condition

interaction_tbl <- interaction_global %>%
  rownames_to_column(var = "Accession") %>%
  as_tibble() 
  

############################
## max_abs_interaction >= 2  ####

interaction_tbl_filt <- interaction_tbl %>%
  mutate(
    max_abs_interaction = pmax(
      abs(genotypesiSNX13.conditionLLOME),
      abs(genotypesiSNX13.conditionRECOVERY),
      na.rm = TRUE
    )
  ) %>%
  filter(adj.P.Val < 0.05, max_abs_interaction >= 2)   # <-- threshold here adj.P.Val < 0.05, max_abs_interaction >= 1


sig_acc <- interaction_tbl_filt$Accession

log2_expr_imp__Scramble_vs_siSNX13_sig <- log2_expr_imp__Scramble_vs_siSNX13[rownames(log2_expr_imp__Scramble_vs_siSNX13) %in% sig_acc, , drop = FALSE]
nrow(log2_expr_imp__Scramble_vs_siSNX13_sig)   # should be 1837




# Build ordered sample table
sample_order <- coldata__Scramble_vs_siSNX13 %>%
  mutate(
    condition = factor(condition, levels = c("DMSO","LLOME","RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram","siSNX13"))
  ) %>%
  arrange(genotype, condition, replicate)

# Reorder matrix
log2_expr_imp__Scramble_vs_siSNX13_sig_ord <- log2_expr_imp__Scramble_vs_siSNX13_sig[, sample_order$sample]




expr_scaled <- t(scale(t(log2_expr_imp__Scramble_vs_siSNX13_sig_ord)))

row_dist   <- dist(expr_scaled, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")




k <- 12
row_clusters <- cutree(row_hclust, k = k)


cluster_tbl <- tibble(
  Accession = rownames(expr_scaled),
  cluster   = row_clusters
)

# add gene symbols
gene_map <- df %>% distinct(Accession, GeneSymbol)
cluster_tbl <- cluster_tbl %>% left_join(gene_map, by = "Accession")


n_rep <- sample_order %>%
  count(genotype, condition)
block_sizes <- n_rep$n


row_annot <- data.frame(Cluster = factor(row_clusters))
rownames(row_annot) <- rownames(expr_scaled)

clust_cols <- setNames(RColorBrewer::brewer.pal(max(3, k), "Set3")[1:k], levels(row_annot$Cluster))

ann_colors <- list(Cluster = clust_cols)




pdf("output/pheatmap-LOGLIMMA-AllNAvaluesFilter-Scramble_vs_siSNX13-12Cluster-max_abs_interaction2.pdf", width = 4, height = 5)
pheatmap(expr_scaled,
         cluster_rows = row_hclust,
         cluster_cols = FALSE,
         cutree_rows = k,
         show_rownames = FALSE,
         gaps_col = head(cumsum(sample_order %>% count(genotype, condition) %>% arrange(genotype, condition) %>% pull(n)), -1),
         annotation_row = row_annot,
         annotation_colors = ann_colors )
dev.off()


cluster_tbl <- cluster_tbl %>%
  left_join(imp_info, by = "Accession") %>%
  mutate(
    n_imputed = replace_na(n_imputed, 0L),
    any_imputed = replace_na(any_imputed, FALSE)
  )
  
write_tsv(cluster_tbl, "output/GENELIST-LOGLIMMA-AllNAvaluesFilter-Scramble_vs_siSNX13-12Cluster-max_abs_interaction2.tsv")



XXXY HERE pick the best clustering and pursue! !!!!


# Show top 10 genes per cluster
## Make a tidy table with cluster + log2(norm+1) values (replicate-level)
prot_tidy_for_plot <- df %>%
  left_join(cluster_tbl %>% select(Accession, GeneSymbol, cluster), by = "Accession") %>%
  filter(!is.na(cluster)) %>%
  mutate(
    condition = factor(condition, levels = c("DMSO", "LLOME", "RECOVERY")),
    genotype  = factor(genotype, levels = c("Scram", "siSNX13")),
    log2_abund = log2(abundance + 1)
  )
top10_per_cluster <- prot_tidy_for_plot %>%
  group_by(cluster, Accession, GeneSymbol.x) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = mean_abund, n = 10, with_ties = FALSE) %>%
  ungroup()
## Summarise mean ± SEM for plotting (keeps bio reps)
plot_df <- prot_tidy_for_plot %>%
  semi_join(top10_per_cluster, by = c("cluster", "Accession", "GeneSymbol.x")) %>%
  group_by(cluster, Accession, GeneSymbol.x, genotype, condition) %>%
  summarise(
    mean_log2 = mean(log2_abund, na.rm = TRUE),
    sem_log2  = sd(log2_abund, na.rm = TRUE) / sqrt(sum(!is.na(log2_abund))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_log2 - sem_log2,
    ymax = mean_log2 + sem_log2,
    panel_label = paste0(GeneSymbol.x, " (", Accession, ")")
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
      title = paste0("Top 10 proteins — Cluster ", cl, " (Scram vs siSNX13)"),
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
    filename = sprintf("output/top10_cluster_%02d-LOGLIMMA-SNX13filterNAvaluesFilter-7Cluster-max_abs_interaction2.pdf", cl),
    plot = p,
    width = 5, height = 6
  )
}






















```

















# Functional analysis 


## Test


We will use clusterProfile package. Tutorial [here](https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html).


Let's test functional analyses (GO BP; KEGG) on the `pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-FC058.csv` version


**IMPORTANT NOTE: When doing GO, do NOT set a universe (background list of genes) it perform better!**


```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")

## Read GTF file
gtf_file <- "../../Master/meta/gencode.v47.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


### GeneSymbol list of signif DEG qval 0.05 FC 0.58
output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-FC058.csv




############ scram_SNX13_10Cluster-FC058 ############
scram_SNX13_10Cluster <- read_tsv("output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-FC058.csv")

ego <- enrichGO(gene = as.character(scram_SNX13_10Cluster$GeneSymbol), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)


pdf("output/GO/dotplot_BP-upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()



## Loop
df <- read_tsv("output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-FC058.csv") %>%
  dplyr::select(GeneSymbol, cluster) %>%
  unique() %>%
  mutate(cluster = as.integer(cluster)) %>%      
  filter(!is.na(cluster)) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(cluster, GeneSymbol, .keep_all = TRUE)

clusters <- sort(unique(df$cluster))
print(clusters)

out_pdf <- "output/GO/dotplot_GOBP-Norm_High_noNA-scram_SNX13_10Cluster.pdf"
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)

for (cl in clusters) {
  genes <- df %>%
    filter(cluster == cl) %>%
    pull(GeneSymbol) %>%
    unique()
  ego <- enrichGO(
    gene          = genes,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("Cluster ", cl, " (n genes = ", length(genes), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No significant GO BP terms (p<0.05)", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)
  } else {
    p <- dotplot(ego, showCategory = 10) +
      ggtitle(paste0("GO BP — Cluster ", cl, " (n genes = ", length(genes), ")")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()















entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(PSC_up$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()

pdf("output/GO/dotplot_KEGG-upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ekegg, showCategory=10)
dev.off()







```





## LOGLIMMA




Let's test functional analyses (GO BP; KEGG) on the `# LOGLIMMA` version






**IMPORTANT NOTE: When doing GO, do NOT set a universe (background list of genes) it perform better!**


```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")

## Read GTF file
gtf_file <- "../../Master/meta/gencode.v47.annotation.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


### GeneSymbol list of signif LOGLIMMA signif

output/GENELIST-LOGLIMMA-SNX13filter-12Cluster.tsv
output/GENELIST-LOGLIMMA-SNX14filter-10Cluster.tsv




############ SNX13filter-12Cluster ############
SNX13filter_12Cluster <- read_tsv("output/GENELIST-LOGLIMMA-SNX13filter-12Cluster.tsv")


## GO BP ################
ego <- enrichGO(gene = as.character(scram_SNX13_10Cluster$GeneSymbol), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)


pdf("output/GO/dotplot_BP-upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()


## Loop
df <- read_tsv("output/GENELIST-LOGLIMMA-SNX13filter-12Cluster.tsv") %>%
  dplyr::select(GeneSymbol, cluster) %>%
  unique() %>%
  mutate(cluster = as.integer(cluster)) %>%      
  filter(!is.na(cluster)) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(cluster, GeneSymbol, .keep_all = TRUE)

clusters <- sort(unique(df$cluster))
print(clusters)

out_pdf <- "output/GO/dotplot_GOBP-LOGLIMMA-SNX13filter-12Cluster.pdf"
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)

for (cl in clusters) {
  genes <- df %>%
    filter(cluster == cl) %>%
    pull(GeneSymbol) %>%
    unique()
  ego <- enrichGO(
    gene          = genes,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("Cluster ", cl, " (n genes = ", length(genes), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No significant GO BP terms (p<0.05)", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)
  } else {
    p <- dotplot(ego, showCategory = 10) +
      ggtitle(paste0("GO BP — Cluster ", cl, " (n genes = ", length(genes), ")")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()




## KEGG ################



entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(PSC_up$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()


# loop
df <- read_tsv("output/GENELIST-LOGLIMMA-SNX13filter-12Cluster.tsv") %>%
  dplyr::select(GeneSymbol, cluster) %>%
  unique() %>%
  mutate(cluster = as.integer(cluster)) %>%
  filter(!is.na(cluster)) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(cluster, GeneSymbol, .keep_all = TRUE)

clusters <- sort(unique(df$cluster))
print(clusters)

out_pdf <- "output/GO/dotplot_KEGG-LOGLIMMA-SNX13filter-12Cluster.pdf"
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)

for (cl in clusters) {

  ## symbols for this cluster
  genes_sym <- df %>%
    filter(cluster == cl) %>%
    pull(GeneSymbol) %>%
    unique()

  ## SYMBOL -> ENTREZID (KEGG uses Entrez IDs)
  entrez <- mapIds(
    x        = org.Hs.eg.db,
    keys     = genes_sym,
    keytype  = "SYMBOL",
    column   = "ENTREZID",
    multiVals = "first"
  )

  entrez_genes <- unique(na.omit(as.character(entrez)))

  ## KEGG enrichment
  ekegg <- enrichKEGG(
    gene          = entrez_genes,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  ## optional: convert Entrez IDs in result back to readable gene symbols
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  ## plot one page per cluster
  if (length(entrez_genes) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("KEGG — Cluster ", cl, " (n genes = ", length(genes_sym), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No Entrez IDs after SYMBOL->ENTREZ conversion", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)
  } else if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("KEGG — Cluster ", cl,
                     " (n genes = ", length(genes_sym),
                     "; n Entrez = ", length(entrez_genes), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No significant KEGG pathways (p<0.05)", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)
  } else {
    p <- dotplot(ekegg, showCategory = 20) +
      ggtitle(paste0("KEGG — Cluster ", cl,
                     " (n genes = ", length(genes_sym),
                     "; n Entrez = ", length(entrez_genes), ")")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()













############ SNX14filter-12Cluster ############
SNX14filter_10Cluster <- read_tsv("output/GENELIST-LOGLIMMA-SNX14filter-10Cluster.tsv")


## GO BP ################
ego <- enrichGO(gene = as.character(SNX14filter_10Cluster$GeneSymbol), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)


pdf("output/GO/dotplot_BP-upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top10.pdf", width=5, height=4)
dotplot(ego, showCategory=10)
dev.off()


## Loop
df <- read_tsv("output/GENELIST-LOGLIMMA-SNX14filter-10Cluster.tsv") %>%
  dplyr::select(GeneSymbol, cluster) %>%
  unique() %>%
  mutate(cluster = as.integer(cluster)) %>%      
  filter(!is.na(cluster)) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(cluster, GeneSymbol, .keep_all = TRUE)

clusters <- sort(unique(df$cluster))
print(clusters)

out_pdf <- "output/GO/dotplot_GOBP-LOGLIMMA-SNX14filter-10Cluster.pdf"
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)

for (cl in clusters) {
  genes <- df %>%
    filter(cluster == cl) %>%
    pull(GeneSymbol) %>%
    unique()
  ego <- enrichGO(
    gene          = genes,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("Cluster ", cl, " (n genes = ", length(genes), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No significant GO BP terms (p<0.05)", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)
  } else {
    p <- dotplot(ego, showCategory = 10) +
      ggtitle(paste0("GO BP — Cluster ", cl, " (n genes = ", length(genes), ")")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()




## KEGG ################



entrez_genes <- as.character( mapIds(org.Hs.eg.db, as.character(PSC_up$gene_name), 'ENTREZID', 'SYMBOL') )

ekegg <- enrichKEGG(gene = entrez_genes, 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)
                
pdf("output/GO/dotplot_KEGG-upregulated_q05fc058_PSC_Hypo_vs_Norm-featurecounts_multi-top20.pdf", width=7, height=7)
dotplot(ekegg, showCategory=20)
dev.off()


# loop
df <- read_tsv("output/GENELIST-LOGLIMMA-SNX14filter-10Cluster.tsv") %>%
  dplyr::select(GeneSymbol, cluster) %>%
  unique() %>%
  mutate(cluster = as.integer(cluster)) %>%
  filter(!is.na(cluster)) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(cluster, GeneSymbol, .keep_all = TRUE)

clusters <- sort(unique(df$cluster))
print(clusters)

out_pdf <- "output/GO/dotplot_KEGG-LOGLIMMA-SNX14filter-10Cluster.pdf"
pdf(out_pdf, width = 7, height = 5, onefile = TRUE)

for (cl in clusters) {

  ## symbols for this cluster
  genes_sym <- df %>%
    filter(cluster == cl) %>%
    pull(GeneSymbol) %>%
    unique()
  ## SYMBOL -> ENTREZID (KEGG uses Entrez IDs)
  entrez <- mapIds(
    x        = org.Hs.eg.db,
    keys     = genes_sym,
    keytype  = "SYMBOL",
    column   = "ENTREZID",
    multiVals = "first"
  )
  entrez_genes <- unique(na.omit(as.character(entrez)))
  ## KEGG enrichment
  ekegg <- enrichKEGG(
    gene          = entrez_genes,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  ## optional: convert Entrez IDs in result back to readable gene symbols
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  ## plot one page per cluster
  if (length(entrez_genes) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("KEGG — Cluster ", cl, " (n genes = ", length(genes_sym), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No Entrez IDs after SYMBOL->ENTREZ conversion", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)
  } else if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(paste0("KEGG — Cluster ", cl,
                     " (n genes = ", length(genes_sym),
                     "; n Entrez = ", length(entrez_genes), ")")) +
      annotate("text", x = 0, y = 0,
               label = "No significant KEGG pathways (p<0.05)", size = 5) +
      xlim(-1, 1) + ylim(-1, 1)
    print(p)
  } else {
    p <- dotplot(ekegg, showCategory = 20) +
      ggtitle(paste0("KEGG — Cluster ", cl,
                     " (n genes = ", length(genes_sym),
                     "; n Entrez = ", length(entrez_genes), ")")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}
dev.off()




```








