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



XXXY HERE!!





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


```


