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


# Testing

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
clusters_filt <- cutree(hc_filt, k = 10)


ord <- order(clusters_filt, -pattern_strength_filt)
mat_log_ord <- mat_log_filt[ord, ]


cluster_annot <- data.frame(
  Cluster = factor(clusters_filt[ord])
)
rownames(cluster_annot) <- rownames(mat_log_ord)



pdf("output/pheatmap-Norm_High_noNA-scram_SNX13_10Cluster-FC058.pdf", width = 3, height = 4)
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







# Functional analysis 



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






