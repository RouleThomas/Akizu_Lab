# Project

**H9 cell lines**
- ESC native:
    - WT: 3 Bio Rep (A1-3)
    - KO: 3 Bio Rep (A4-6)
    - KOEF1aEZH1: 3 Bio Rep (A7-9)


--> Directional mRNA library preparation (poly A enrichment), NovaSeq X Plus Series (PE150)




**Objectives:**
- Put together with CutRun 018


Novogene Input	Sample Name	
A1	R1 EZH1 WT	ESC_WT_R1
A2	R1 EZH1 KO	ESC_KO_R1
A3	R1 EZH1 KO (50 ng dox)	R1 EZH1 OEKO    ESC_OEKO_R1
A4	R2 EZH1 WT	    ESC_WT_R2
A5	R2 EZH1 KO	    ESC_KO_R2
A6	R2 EZH1 KO (50 ng dox)	R2 EZH1 OEKO    ESC_OEKO_R2
A7	R3 EZH1 WT	    ESC_WT_R3
A8	R3 EZH1 KO	    ESC_KO_R3
A9	R3 EZH1 KO (50 ng dox)	R3 EZH1 OEKO    ESC_OEKO_R3





# Pipeline
- Download data (wget)
- Rename files
- FastQC (fastqc)
- Trimming (fastp)
- Count with Kallisto (better than Salmon for switch variation analysis, and better than featureCounts; best to do transcriptome level coutning and then to gene
)

--> Detail of the overall pipeline in `Meeting_20230919_draft.xlsx` 

# Download / import data


```bash
# Following email instructions
module load lftp
lftp -c 'set sftp:auto-confirm yes;set net:max-retries 20;open sftp://X202SC25078875-Z01-F001:wfcxcay4@usftp23.novogene.com; mirror --verbose --use-pget-n=8 -c'

# Copy all .fz.gz data into input/ folder
rsync -av --include '*/' --include '*.fq.gz' --exclude '*' usftp23.novogene.com/ input/ # copy from usftp23 folder to input
find input/ -mindepth 2 -type f -exec mv -t input/ {} + # mv files from their folder to input/ folder
find input/ -type d -empty -delete # delete empty directory

```

--> All good, files created in `usftp23.novogene.com/`




# Rename file

Renamed manually as only 8 samples

--> All good 



# Fastp cleaning

```bash
sbatch scripts/fastp.sh # 50212025 ok
```

## mapping fastp trim

```bash
sbatch --dependency=afterany:50212025 scripts/STAR_mapping_fastp.sh # 50212217 ok


# Convert alignment to bigwig
conda activate deeptools

sbatch scripts/STAR_TPM_bw.sh # 50944526 xxx
```

-->  ok







# Count with Kallisto

Follow instruction [here](https://pachterlab.github.io/kallisto/download.html)

## count with Kallisto

```bash
conda activate kallisto


## run in sbatch
sbatch scripts/kallisto_count_gtf.sh # 50212654 FAIL SHOULD USE STRANDED!; 50302968 ok

# Convert pseudoalignment to bigwig
conda activate deeptools

sbatch scripts/TPM_bw.sh # 50314231 ok


# Calculate median
conda activate BedToBigwig

sbatch scripts/bigwigmerge_TPM_bw.sh # 50944579 ok
```

- *NOTE: Added `--rf-stranded --genomebam` options for strandness and pseudobam alignemt generation*


### deepTool bigwig QC


### Raw bigwig


```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_TPM.sh # 50503924 xxx



############################################
# Plot ESC ###########
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_TPM.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_WT_R1 ESC_WT_R2 ESC_WT_R3 ESC_KO_R1 ESC_KO_R2 ESC_KO_R3 ESC_OEKO_R1 ESC_OEKO_R2 ESC_OEKO_R3 \
    --colors black black black red red red blue blue blue \
    --markers 's' 'o' '>' 's' 'o' '>' 's' 'o' '>' \
    -o output/bigwig/multiBigwigSummary_TPM_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_TPM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_WT_R1 ESC_WT_R2 ESC_WT_R3 ESC_KO_R1 ESC_KO_R2 ESC_KO_R3 ESC_OEKO_R1 ESC_OEKO_R2 ESC_OEKO_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_TPM_heatmap.pdf

#################################
```

--> WT cluster well, KO and OEKO kind of overlap..







## Gene quantification with txImport and DESEQ2



Follow [this](https://nbisweden.github.io/workshop-RNAseq/2011/lab_kallisto.html#2_Quantification)

and [this](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering) for DESEQ2 part

```bash
# Extract all transcriptnames (1st) and genenames (4th) from  GTF and write to a file.  
## Keep gene version information 
awk 'BEGIN{OFS=","; print "TXNAME,GENEID"}
     $3=="transcript"{
       match($0,/transcript_id "([^"]+)"/,t);
       match($0,/gene_id "([^"]+)"/,g);
       if(t[1]!="" && g[1]!=""){
         tx=t[1]; gn=g[1];
         sub(/\.[0-9]+$/,"",tx);
         sub(/\.[0-9]+$/,"",gn);
         print tx, gn;
       }
     }' ../../Master/meta/gencode.v47.annotation.gtf \
| sort -u > ../../Master/meta/gencode.v47.annotation.tx2gene.csv

## Remove gene version information 
echo "TXNAME,GENEID" > ../../Master/meta/gencode.v47.annotation.tx2gene.csv

awk '$3=="transcript"{
       match($0,/transcript_id "([^"]+)"/,t);
       match($0,/gene_id "([^"]+)"/,g);
       if(t[1]!="" && g[1]!="") print t[1] "," g[1];
     }' OFS="," ../../Master/meta/gencode.v47.annotation.gtf \
| sort -u >> ../../Master/meta/gencode.v47.annotation.tx2gene.csv


conda activate deseq2
```

Go in R to create metadata file; and convert transcript to gene count

Metadata file format as: SampleName, SampleID, No, Model, Day, Group, Replicate


```R
# packages
library("tidyverse")
library("dplyr") # data wrangling
library("ggplot2") # plotting
library("DESeq2") # rna-seq
library("tximport") # importing kalisto transcript counts to geneLevels
library("readr") # Fast readr of files.
library("rhdf5") # read/convert kalisto output files.  

library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")

set.seed(42)

########################################
## WT KO OEKO - gene level count ############################
########################################


# Create metadata
samples <- c(
  "ESC_WT_R1","ESC_KO_R1","ESC_OEKO_R1",
  "ESC_WT_R2","ESC_KO_R2","ESC_OEKO_R2",
  "ESC_WT_R3","ESC_KO_R3","ESC_OEKO_R3"
)
mr <- data.frame(
  SampleID   = samples,
  No         = 1:length(samples),
  Model      = "ESC",
  Day        = "ESC",
  Group      = sub("ESC_","", sub("_R[0-9]","", samples)),
  Replicate  = sub(".*_R","", samples),
  row.names  = samples         # <-- set SampleName as rownames
)
mr


# List all abundance.tsv files
files <- setNames(
  file.path("output/kallisto", paste0(samples, "_quant"), "abundance.tsv"),
  samples
)
# Name the files with your sample names
names(files) <- rownames(mr)
files



# Convert transcript to gene ID
tx2gene <- read_csv("../../Master/meta/gencode.v47.annotation.tx2gene.csv")

txi <- tximport(
  files,
  type = "none",
  tx2gene = tx2gene,
  ignoreAfterBar = TRUE,
  importer = function(f) readr::read_tsv(f, show_col_types = FALSE),
  txIdCol      = "target_id",
  countsCol    = "est_counts",
  abundanceCol = "tpm",
  lengthCol    = "eff_length"
)

tpm_long <- txi$abundance %>%                 # <- matrix: genes x samples (TPM)
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(-GeneID, names_to = "sample", values_to = "TPM") %>%
  mutate(
    replicate    = str_extract(sample, "R\\d+"),
    genotype_raw = str_match(sample, "_(WT|KO|OE\\w*)_")[,2],  # WT / KO / OE*
    genotype     = case_when(
      genotype_raw %in% c("WT","KO") ~ genotype_raw,
      str_detect(genotype_raw, "^OE") ~ "OE",
      TRUE ~ genotype_raw
    )
  ) %>%
  dplyr::select(GeneID, genotype, replicate, TPM) %>%
  arrange(GeneID, genotype, replicate)

tpm_long



# Add geneSymbol
tpm_long <- tpm_long %>%
  mutate(ENSG = sub("\\..*$", "", GeneID))


# 2. Map ENSG → GeneSymbol with org.Hs.eg.db
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = tpm_long$ENSG,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# 3. Add GeneSymbol column
tpm_long$GeneSymbol <- gene_symbols[tpm_long$ENSG]

# 4. Reorder columns (put GeneSymbol after GeneID)
tpm_long <- tpm_long %>%
  relocate(GeneSymbol, .after = GeneID)


# Save gene level count #######################################
write.table(tpm_long, file = "output/deseq2/txi_Kallisto-GeneLevel-TPM.txt", sep = "\t", quote = FALSE, row.names = FALSE) # 


## plot single gene TPM  ################################################
df_gene <- tpm_long %>%
  filter(GeneSymbol == "EZH1") %>%
  mutate(log2TPM = log2(TPM + 1),
         genotype = factor(genotype, levels = c("WT","KO","OE")))  # adjust/order if needed

# choose the pairwise comparisons you want to test
comparisons <- list(c("WT","KO"), c("WT","OE"), c("KO","OE"))

pdf("output/deseq2/barplot_TPM-WTvsKOvsOEKO-EZH1.pdf", width = 3, height =4)
ggbarplot(
  df_gene,
  x = "genotype",
  y = "log2TPM",
  add = c("mean_se", "jitter"),     # bars = mean, error bars = SE, overlay points
  add.params = list(size = 1, alpha = 0.7), # point styling
  fill = "genotype",                 # color bars by group (optional)
  palette = c("WT" = "black",
            "KO" = "red",
            "OE" = "blue"),
  width = 0.7
) +
  labs(
    title = "EZH1",
    x = "Genotype",
    y = "log2(TPM + 1)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  stat_compare_means(
    comparisons = comparisons,
    method = "t.test",               # or "wilcox.test" if non-parametric
    label = "p.signif",
    hide.ns = TRUE
  ) 
dev.off()



##############################











## Plot H3K27me3 GAIN - lost genes - HEATMAP
genes_gain <- read.table(
  "../018__CutRun_DOX_V1/output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt",
  header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE,
  col.names = "GeneSymbol"
) %>% 
  dplyr::pull(GeneSymbol) %>% unique()

## 2) Average TPM per genotype (WT/KO) for those genes
avg_tpm <- tpm_long %>%
  filter(GeneSymbol %in% genes_gain, genotype %in% c("WT","KO")) %>%
  group_by(GeneSymbol, genotype) %>%
  summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop")

mat <- avg_tpm %>%
  tidyr::pivot_wider(names_from = genotype, values_from = TPM) %>%
  # drop genes with missing GeneSymbol
  filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
  # order by WT expression (high → low)
  arrange(desc(WT)) %>%
  # enforce WT first, KO second
  dplyr::select(GeneSymbol, WT, KO)

# convert to matrix
mat <- mat %>%
  tibble::column_to_rownames("GeneSymbol") %>%
  as.matrix()

# replace any NAs with 0
mat[is.na(mat)] <- 0

## Plot heatmap (TPM color scale; I like log10+1 to compress range)
pdf("output/deseq2/heatmap_TPM-H3K27me3-WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain.pdf",
    width = 2, height = max(4, nrow(mat) * 0.18))
pheatmap(
  log2(mat + 1),
  color = colorRampPalette(c("white","orange","red"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "TPM (log2+1)",
  border_color = NA
)
dev.off()


avg_tpm$genotype <-
  factor(avg_tpm$genotype,
         c("WT", "KO"))

pdf("output/deseq2/boxplot_TPM-H3K27me3-WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain.pdf", width = 4, height =4)
ggplot(avg_tpm, aes(x = genotype, y = log2(TPM+1), fill = genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +        # boxplot, light fill
  geom_jitter(aes(color = genotype), width = 0.2, size = 2, alpha = 0.5) +  # dots
  scale_fill_manual(values = c("WT" = "grey70", "KO" = "red")) +
  scale_color_manual(values = c("WT" = "grey40", "KO" = "darkred")) +
  theme_classic(base_size = 14) +
  labs(
    title = "H3K27me3-gain genes",
    x = NULL, y = "TPM (log2+1)"
  )
dev.off()






## Plot H3K27me3 LOST - lost genes - HEATMAP 
genes_lost <- read.table(
  "../018__CutRun_DOX_V1/output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt",
  header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE,
  col.names = "GeneSymbol"
) %>% 
  dplyr::pull(GeneSymbol) %>% unique()

## 2) Average TPM per genotype (WT/KO) for those genes
avg_tpm <- tpm_long %>%
  filter(GeneSymbol %in% genes_lost, genotype %in% c("WT","KO")) %>%
  group_by(GeneSymbol, genotype) %>%
  summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop")

mat <- avg_tpm %>%
  tidyr::pivot_wider(names_from = genotype, values_from = TPM) %>%
  # drop genes with missing GeneSymbol
  filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
  # order by WT expression (high → low)
  arrange(desc(WT)) %>%
  # enforce WT first, KO second
  dplyr::select(GeneSymbol, WT, KO)

# convert to matrix
mat <- mat %>%
  tibble::column_to_rownames("GeneSymbol") %>%
  as.matrix()

# replace any NAs with 0
mat[is.na(mat)] <- 0

## 4) Plot heatmap (TPM color scale; I like log10+1 to compress range)
pdf("output/deseq2/heatmap_TPM-H3K27me3-WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost.pdf",
    width = 2, height = max(4, nrow(mat) * 0.18))
pheatmap(
  log2(mat + 1),
  color = colorRampPalette(c("white","orange","red"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "TPM (log2+1)",
  border_color = NA
)
dev.off()


avg_tpm$genotype <-
  factor(avg_tpm$genotype,
         c("WT", "KO"))

pdf("output/deseq2/boxplot_TPM-H3K27me3-WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost.pdf", width = 4, height =4)
ggplot(avg_tpm, aes(x = genotype, y = log2(TPM+1), fill = genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +        # boxplot, light fill
  geom_jitter(aes(color = genotype), width = 0.2, size = 2, alpha = 0.5) +  # dots
  scale_fill_manual(values = c("WT" = "grey70", "KO" = "red")) +
  scale_color_manual(values = c("WT" = "grey40", "KO" = "darkred")) +
  theme_classic(base_size = 14) +
  labs(
    title = "H3K27me3-gain genes",
    x = NULL, y = "TPM (log2+1)"
  )
dev.off()





################################################################


########################################
## WT VS KO ############################
########################################


# Create metadata
samples <- c(
  "ESC_WT_R1","ESC_KO_R1",
  "ESC_WT_R2","ESC_KO_R2",
  "ESC_WT_R3","ESC_KO_R3"
)
mr <- data.frame(
  SampleID   = samples,
  No         = 1:length(samples),
  Model      = "ESC",
  Day        = "ESC",
  Group      = sub("ESC_","", sub("_R[0-9]","", samples)),
  Replicate  = sub(".*_R","", samples),
  row.names  = samples         # <-- set SampleName as rownames
)
mr


# List all abundance.tsv files
files <- setNames(
  file.path("output/kallisto", paste0(samples, "_quant"), "abundance.tsv"),
  samples
)
# Name the files with your sample names
names(files) <- rownames(mr)
files



# Convert transcript to gene ID
tx2gene <- read_csv("../../Master/meta/gencode.v47.annotation.tx2gene.csv")

txi <- tximport(
  files,
  type = "none",
  tx2gene = tx2gene,
  ignoreAfterBar = TRUE,
  importer = function(f) readr::read_tsv(f, show_col_types = FALSE),
  txIdCol      = "target_id",
  countsCol    = "est_counts",
  abundanceCol = "tpm",
  lengthCol    = "eff_length"
)



# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)

# Normalize IDs (remove version suffixes in both sets)
tx_genes      <- rownames(txi$counts)
tx_genes_nov  <- sub("\\..*$", "", tx_genes)                  # ENSG00000123456.3 -> ENSG00000123456
xy_genes_nov  <- sub("\\..*$", "", genes_X_Y$ensembl_gene_id) # usually no version, but safe

# Build keep mask and filter all txi matrices
keep <- !(tx_genes_nov %in% xy_genes_nov)

txi$counts    <- txi$counts[keep, , drop = FALSE]
txi$abundance <- txi$abundance[keep, , drop = FALSE]
txi$length    <- txi$length[keep, , drop = FALSE]

cat("Removed", sum(!keep), "X/Y genes; kept", sum(keep), "genes.\n")




# DGE with DESEQ2
dds <- DESeqDataSetFromTximport(txi, mr, ~Group)

# Pre-filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Assign ref
dds$Group <- relevel(dds$Group, ref = "WT")

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, name="Group_KO_vs_WT")

# shrinkage
resultsNames(dds)
res <- lfcShrink(dds, coef="Group_KO_vs_WT", type="apeglm")
res

# explore output DEG
res05 <- results(dds, alpha=0.05)
summary(res05)






# Identify DEGs and count them

## padj 0.05 FC 0.25 ##################################
res_df <- res %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.25 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.25 & res_df$padj == TRUE, na.rm = TRUE)

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.25 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.25 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.25)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.25)'

pdf("output/deseq2/plotVolcano_res_q05fc025_ESC_KO_vs_ESC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.25,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2/res_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.25 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)






## padj 0.05 log2FC 0.58 (~50% increase) ##################################
res_df <- res %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.58 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'

pdf("output/deseq2/plotVolcano_res_q05fc058_ESC_KO_vs_ESC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2/res_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.58 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.58 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc058_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc058_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)






######################################
# GAIN H3K27me3 from `001*/018*` ###################
######################################

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols



res_Gain= read.table("../018__CutRun_DOX_V1/output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE, col.names = "GeneSymbol") %>%
  as_tibble() %>%
  inner_join(as_tibble(res), by = "GeneSymbol") %>% unique() %>% dplyr::filter(GeneSymbol != "NA")


res_df <- res_Gain %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.25 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.25 & res_df$padj == TRUE, na.rm = TRUE)


# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.25 & res_Gain$padj < 5e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.25 & res_Gain$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.25)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.25)'

pdf("output/deseq2/plotVolcano_res_Gain_q05fc025_ESC_KO_vs_ESC_WT.pdf", width=4, height=5)    
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.25,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 75, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 75, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
#write.table(res, file = "output/deseq2/res_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.25 & res_Gain$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.25 & res_Gain$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_ESC_KO_vs_ESC_WT-res_Gain.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_ESC_KO_vs_ESC_WT-res_Gain.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



######################################
# LOST H3K27me3 from `001*/018*` ###################
######################################

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols



res_Lost= read.table("../018__CutRun_DOX_V1/output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE, col.names = "GeneSymbol") %>%
  as_tibble() %>%
  inner_join(as_tibble(res), by = "GeneSymbol") %>% unique() %>% dplyr::filter(GeneSymbol != "NA")


res_df <- res_Lost %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.25 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.25 & res_df$padj == TRUE, na.rm = TRUE)


# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_Lost$log2FoldChange < -0.25 & res_Lost$padj < 5e-2, 'Sky Blue',
    ifelse(res_Lost$log2FoldChange > 0.25 & res_Lost$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.25)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.25)'

pdf("output/deseq2/plotVolcano_res_Lost_q05fc025_ESC_KO_vs_ESC_WT.pdf", width=4, height=5)    
EnhancedVolcano(res_Lost,
  lab = res_Lost$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.25,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 75, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 75, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
#write.table(res, file = "output/deseq2/res_ESC_KO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange > 0.25 & res_Gain$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_Gain[!is.na(res_Gain$log2FoldChange) & !is.na(res_Gain$padj) & res_Gain$log2FoldChange < -0.25 & res_Gain$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_ESC_KO_vs_ESC_WT-res_Gain.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_ESC_KO_vs_ESC_WT-res_Gain.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


















########################################
## WT VS OEKO ############################
########################################


# Create metadata
samples <- c(
  "ESC_WT_R1","ESC_OEKO_R1",
  "ESC_WT_R2","ESC_OEKO_R2",
  "ESC_WT_R3","ESC_OEKO_R3"
)
mr <- data.frame(
  SampleID   = samples,
  No         = 1:length(samples),
  Model      = "ESC",
  Day        = "ESC",
  Group      = sub("ESC_","", sub("_R[0-9]","", samples)),
  Replicate  = sub(".*_R","", samples),
  row.names  = samples         # <-- set SampleName as rownames
)
mr


# List all abundance.tsv files
files <- setNames(
  file.path("output/kallisto", paste0(samples, "_quant"), "abundance.tsv"),
  samples
)
# Name the files with your sample names
names(files) <- rownames(mr)
files



# Convert transcript to gene ID
tx2gene <- read_csv("../../Master/meta/gencode.v47.annotation.tx2gene.csv")

txi <- tximport(
  files,
  type = "none",
  tx2gene = tx2gene,
  ignoreAfterBar = TRUE,
  importer = function(f) readr::read_tsv(f, show_col_types = FALSE),
  txIdCol      = "target_id",
  countsCol    = "est_counts",
  abundanceCol = "tpm",
  lengthCol    = "eff_length"
)



# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)

# Normalize IDs (remove version suffixes in both sets)
tx_genes      <- rownames(txi$counts)
tx_genes_nov  <- sub("\\..*$", "", tx_genes)                  # ENSG00000123456.3 -> ENSG00000123456
xy_genes_nov  <- sub("\\..*$", "", genes_X_Y$ensembl_gene_id) # usually no version, but safe

# Build keep mask and filter all txi matrices
keep <- !(tx_genes_nov %in% xy_genes_nov)

txi$counts    <- txi$counts[keep, , drop = FALSE]
txi$abundance <- txi$abundance[keep, , drop = FALSE]
txi$length    <- txi$length[keep, , drop = FALSE]

cat("Removed", sum(!keep), "X/Y genes; kept", sum(keep), "genes.\n")




# DGE with DESEQ2
dds <- DESeqDataSetFromTximport(txi, mr, ~Group)

# Pre-filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Assign ref
dds$Group <- relevel(dds$Group, ref = "WT")

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, name="Group_OEKO_vs_WT")

# shrinkage
resultsNames(dds)
res <- lfcShrink(dds, coef="Group_OEKO_vs_WT", type="apeglm")
res

# explore output DEG
res05 <- results(dds, alpha=0.05)
summary(res05)






# Identify DEGs and count them

## padj 0.05 FC 0.25 ##################################
res_df <- res %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.25 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.25 & res_df$padj == TRUE, na.rm = TRUE)

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.25 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.25 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.25)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.25)'

pdf("output/deseq2/plotVolcano_res_q05fc025_ESC_OEKO_vs_ESC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'OEKO vs WT, ESC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.25,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2/res_ESC_OEKO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.25 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_ESC_OEKO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_ESC_OEKO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)










## padj 0.05 FC 0.58 ##################################
res_df <- res %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.58 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'

pdf("output/deseq2/plotVolcano_res_q05fc058_ESC_OEKO_vs_ESC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'OEKO vs WT, ESC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 2.0,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24) ) +
  annotate("text", x = 3, y = 140, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 140, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2/res_ESC_OEKO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.58 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.58 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc058_ESC_OEKO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc058_ESC_OEKO_vs_ESC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


```








# Shiny app


generate the `tpm_all_sample.txt` file and then go into `001_EZH1*/001__RNAseq` to create shiny app V2 including these

```R
library("tidyverse")
library("biomaRt")

# Import TPM (done in Kallisto)
long_data <- read.delim("output/deseq2/txi_Kallisto-GeneLevel-TPM.txt", 
                 header = TRUE, 
                 sep = "\t", 
                 stringsAsFactors = FALSE) %>%
              as_tibble()


long_data_log2tpm = long_data %>%
  mutate(log2tpm = log2(TPM + 1)) %>%
  dplyr::rename("Genotype" = "genotype",
                "Replicate" = "replicate",
                "external_gene_name" = "GeneSymbol") %>%
  add_column(Tissue = "ESC") %>%
  dplyr::select(external_gene_name, Tissue, Genotype, Replicate, TPM, log2tpm)

## Save
write.table(long_data_log2tpm, file = c("output/deseq2/long_data_log2tpm_Akoto001019.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


```


--> Next here: `001_EZH1*/001__RNAseq` (`## Shiny app V3; including Ciceri RNAseq neuron diff dataset + Akoto 001/015 RNAseq`)




















