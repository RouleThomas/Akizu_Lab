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


## Convert alignment to bigwig
conda activate deeptools
sbatch scripts/STAR_TPM_bw.sh # 50944526 ok

## Calculate median
conda activate BedToBigwig
sbatch scripts/bigwigmerge_STAR_TPM_bw.sh # 51671403 ok




# Remove the X chromomsoem
sbatch scripts/STAR_removeXchr.sh # 51676961 ok

## Convert alignment to bigwig
sbatch --dependency=afterany:51676961 scripts/STAR_noXchr_TPM_bw.sh # 51676963 ok



```

-->  ok

--> X chr removed using `samtools view`




### deepTool bigwig QC - STAR mapping


### Raw bigwig


```bash
conda activate deeptools

###################################
# Include X chr #####################
###################################
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_STAR_TPM.sh # 51072872 ok


############################################
# Plot ESC ###########
## PCA
plotPCA -in output/bigwig_STAR/multiBigwigSummary_STAR_TPM.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_WT_R1 ESC_WT_R2 ESC_WT_R3 ESC_KO_R1 ESC_KO_R2 ESC_KO_R3 ESC_OEKO_R1 ESC_OEKO_R2 ESC_OEKO_R3 \
    --colors black black black red red red blue blue blue \
    --markers 's' 'o' '>' 's' 'o' '>' 's' 'o' '>' \
    -o output/bigwig_STAR/multiBigwigSummary_STAR_TPM_plotPCA.pdf \
    --plotWidth 7

## Heatmap
plotCorrelation \
    -in output/bigwig_STAR/multiBigwigSummary_STAR_TPM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_WT_R1 ESC_WT_R2 ESC_WT_R3 ESC_KO_R1 ESC_KO_R2 ESC_KO_R3 ESC_OEKO_R1 ESC_OEKO_R2 ESC_OEKO_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_STAR/multiBigwigSummary_STAR_TPM_heatmap.pdf

#################################




###################################
# Without X chr #####################
###################################

sbatch scripts/multiBigwigSummary_STAR_TPM_noXchr.sh # 51799981 xxx


############################################
# Plot ESC ###########
## PCA
plotPCA -in output/bigwig_STAR/multiBigwigSummary_STAR_noXchr_TPM.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_WT_R1 ESC_WT_R2 ESC_WT_R3 ESC_KO_R1 ESC_KO_R2 ESC_KO_R3 ESC_OEKO_R1 ESC_OEKO_R2 ESC_OEKO_R3 \
    --colors black black black red red red blue blue blue \
    --markers 's' 'o' '>' 's' 'o' '>' 's' 'o' '>' \
    -o output/bigwig_STAR/multiBigwigSummary_STAR_noXchr_TPM_plotPCA.pdf \
    --plotWidth 7

## Heatmap
plotCorrelation \
    -in output/bigwig_STAR/multiBigwigSummary_STAR_noXchr_TPM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_WT_R1 ESC_WT_R2 ESC_WT_R3 ESC_KO_R1 ESC_KO_R2 ESC_KO_R3 ESC_OEKO_R1 ESC_OEKO_R2 ESC_OEKO_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_STAR/multiBigwigSummary_STAR_noXchr_TPM_heatmap.pdf

#################################




```

--> **STAR plot PCA cluster the same as Kallisto one**; meaning the bad clustering is not due to the new Kallisto method used.

--> With or without X chr does not change things a lot at the RNA level; PCA looks very similar





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


### deepTool bigwig QC - Kallisto mapping


### Raw bigwig


```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_TPM.sh # 50503924 ok



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
# GAIN H3K27me3 from `001*/018*` DIFFREPS ###################
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
# LOST H3K27me3 from `001*/018*` DIFFREPS ###################
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







######################################
# GAIN H3K27me3 from `001*/018*` DESEQ2 MACS2 ###################
######################################
res <- read.table("output/deseq2/res_ESC_KO_vs_ESC_WT.txt",
                  sep = "\t", 
                  header = TRUE, 
                  row.names = 1, 
                  stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column(var = "gene") %>%
  as_tibble()



gene_Gain <- read.csv(
  "../018__CutRun_DOX_V1/output/edgeR/DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  # 1) keep significant peaks and assign effect direction
  filter(padj < 0.05) %>%
  mutate(effect_dir = case_when(
    log2FoldChange >=  0.58 ~ "pos",
    log2FoldChange <= -0.58 ~ "neg",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(effect_dir)) %>%
  # 2) determine gene-level direction
  group_by(geneSymbol) %>%
  summarise(
    n_pos = sum(effect_dir == "pos"),
    n_neg = sum(effect_dir == "neg"),
    direction = case_when(
      n_pos > 0 & n_neg > 0 ~ "Mix",
      n_pos > 0             ~ "Gain",
      n_neg > 0             ~ "Lost",
      TRUE                  ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  # 3) keep only Gain genes (exclude Mix)
  filter(direction == "Gain") %>%
  dplyr::select(geneSymbol) %>%
  unique() %>%
  dplyr::rename("GeneSymbol" = "geneSymbol")
  


res_Gain = res %>%
  inner_join(gene_Gain)



res_df <- res_Gain %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)



keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.58 & res_Gain$padj < 5e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.58 & res_Gain$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.25)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.25)'

pdf("output/deseq2/plotVolcano_res_Gain-DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3-q05fc058_ESC_KO_vs_ESC_WT.pdf", width=4, height=5)    
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
  annotate("text", x = 3, y = 15, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 15, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()









## WHY SOME GENES ARE MISSING between gain and RNAseq ########################
gene_Gain <- read.csv(
  "../018__CutRun_DOX_V1/output/edgeR/DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  # 1) keep significant peaks and assign effect direction
  filter(padj < 0.05) %>%
  mutate(effect_dir = case_when(
    log2FoldChange >=  0.58 ~ "pos",
    log2FoldChange <= -0.58 ~ "neg",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(effect_dir)) %>%
  # 2) determine gene-level direction
  group_by(geneSymbol) %>%
  summarise(
    n_pos = sum(effect_dir == "pos"),
    n_neg = sum(effect_dir == "neg"),
    direction = case_when(
      n_pos > 0 & n_neg > 0 ~ "Mix",
      n_pos > 0             ~ "Gain",
      n_neg > 0             ~ "Lost",
      TRUE                  ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  # 3) keep only Gain genes (exclude Mix)
  filter(direction == "Gain") %>%
  dplyr::select(geneSymbol) %>%
  unique() %>%
  dplyr::rename("GeneSymbol" = "geneSymbol")

## Seprate gene presnet and absent after fusion RNAseq and CutRun
res_Gain_MISSED = gene_Gain %>%
  anti_join(res, by = "GeneSymbol")

res_Gain_PRESENT = gene_Gain %>%
  inner_join(res, by = "GeneSymbol")


## Load TPM info of each genes
tpm_long <- read.csv(
  "output/deseq2/txi_Kallisto-GeneLevel-TPM.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

avg_tpm <- tpm_long %>%
  group_by(GeneSymbol, genotype) %>%
  summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop")


mat <- avg_tpm %>%
  tidyr::pivot_wider(names_from = genotype, values_from = TPM) %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
  arrange(desc(WT)) %>%
  dplyr::select(GeneSymbol, WT, KO) %>%
  tibble::column_to_rownames("GeneSymbol") %>%
  as.matrix()

## 5) Replace any NAs with 0
mat[is.na(mat)] <- 0

## 6) Plot heatmap (log2+1 TPM)
# 0) Inputs: tpm_long, res_Gain_MISSED (GeneSymbol), res_Gain_PRESENT (GeneSymbol)
genes_missing <- unique(res_Gain_MISSED$GeneSymbol)
genes_present <- unique(res_Gain_PRESENT$GeneSymbol)

# 1) Helper to build WT/KO matrix for a gene list
make_mat <- function(genes) {
  avg_tpm <- tpm_long %>%
    dplyr::filter(GeneSymbol %in% genes, genotype %in% c("WT", "KO")) %>%
    group_by(GeneSymbol, genotype) %>%
    summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop")

  mat <- avg_tpm %>%
    pivot_wider(names_from = genotype, values_from = TPM) %>%
    mutate(WT = coalesce(WT, 0), KO = coalesce(KO, 0)) %>%
    filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
    arrange(desc(WT)) %>%
    dplyr::select(GeneSymbol, WT, KO) %>%
    tibble::column_to_rownames("GeneSymbol") %>%
    as.matrix()

  mat[is.na(mat)] <- 0
  mat
}

mat_missing <- make_mat(genes_missing)
mat_present <- make_mat(genes_present)

# 2) Shared color scale across both panels
vmin  <- min(log2(c(mat_missing, mat_present) + 1), na.rm = TRUE)
vmax  <- max(log2(c(mat_missing, mat_present) + 1), na.rm = TRUE)
ncols <- 100
cols  <- colorRampPalette(c("white","orange","red"))(ncols)
brks  <- seq(vmin, vmax, length.out = ncols + 1)

# 3) Draw side-by-side into one PDF
pdf("output/deseq2/heatmap_TPM-H3K27me3-WTvsKO_missing_vs_present_Gain.pdf",
    width = 6,
    height = max(4, max(nrow(mat_missing), nrow(mat_present)) * 0.18))

p1 <- pheatmap(
  log2(mat_missing + 1),
  color = cols, breaks = brks,
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "Missing in RNA (Gain genes)",
  border_color = NA, silent = TRUE
)

p2 <- pheatmap(
  log2(mat_present + 1),
  color = cols, breaks = brks,
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "Present in RNA (Gain genes)",
  border_color = NA, silent = TRUE
)

gridExtra::grid.arrange(grobs = list(p1$gtable, p2$gtable), ncol = 2)
dev.off()




# Helper: average TPM per gene per genotype (one dot per gene)
avg_by_gene <- function(genes, status_label) {
  tpm_long %>%
    filter(GeneSymbol %in% genes, genotype %in% c("WT","KO")) %>%
    group_by(GeneSymbol, genotype) %>%
    summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(status = status_label)
}

df_missing  <- avg_by_gene(genes_missing, "Missing in RNA (Gain)")
df_present  <- avg_by_gene(genes_present, "Present in RNA (Gain)")
df_plot <- bind_rows(df_missing, df_present) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT","KO")),
    status   = factor(status, levels = c("Missing in RNA (Gain)", "Present in RNA (Gain)"))
  )

# Plot
pdf("output/deseq2/boxplot_TPM-H3K27me3-WTvsKO_missing_vs_present_Gain.pdf",
    width = 5, height = 4.5)

ggplot(df_plot, aes(x = genotype, y = log2(TPM + 1), fill = genotype, color = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.25) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
             size = 1.6, alpha = 0.75) +
  facet_wrap(~ status, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = c(WT = "black", KO = "red")) +
  scale_color_manual(values = c(WT = "black", KO = "red")) +
  labs(x = NULL, y = "TPM (log2 + 1)", title = "WT vs KO — Gain genes (Missing vs Present in RNA)") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dev.off()









# --- 2) Create one jitter value PER GENE (per facet) and reuse for both WT & KO
# Segments data (one row per gene/facet), with WT/KO y-values and shared jitter

second <- "KO"          # <-- change to "OEKO" if that's your column value
levels2 <- c("WT", second)

# -- inputs assumed: tpm_long, res_Gain_MISSED (GeneSymbol), res_Gain_PRESENT (GeneSymbol)
genes_missing <- unique(res_Gain_MISSED$GeneSymbol)
genes_present <- unique(res_Gain_PRESENT$GeneSymbol)

avg_by_gene <- function(genes, status_label) {
  tpm_long %>%
    filter(GeneSymbol %in% genes, genotype %in% levels2) %>%
    group_by(GeneSymbol, genotype) %>%
    summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(status = status_label)
}


df_missing <- avg_by_gene(genes_missing, "Missing in RNA (Gain)")
df_present <- avg_by_gene(genes_present, "Present in RNA (Gain)")

df_plot <- bind_rows(df_missing, df_present) %>%
  mutate(
    genotype = factor(genotype, levels = levels2),
    status   = factor(status, levels = c("Missing in RNA (Gain)", "Present in RNA (Gain)")),
    y        = log2(TPM + 1)   # <-- THIS is the 'y' column
  )

# one jitter per gene (per facet), reused for both genotypes
df_segs <- df_plot %>%
  dplyr::select(GeneSymbol, status, genotype, y) %>%
  pivot_wider(names_from = genotype, values_from = y) %>%
  dplyr::filter(!is.na(.data[[levels2[1]]]) & !is.na(.data[[levels2[2]]])) %>%
  mutate(
    j    = rnorm(n(), 0, 0.12),
    x_1  = 1 + j,
    x_2  = 2 + j
  )

df_pts <- df_plot %>%
  left_join(df_segs %>% dplyr::select(GeneSymbol, status, j), by = c("GeneSymbol","status")) %>%
  mutate(xj = as.numeric(genotype) + j)




df_plot <- bind_rows(df_missing, df_present) %>%
  mutate(
    genotype = factor(genotype, levels = levels2),
    status   = factor(status, levels = c("Missing in RNA (Gain)", "Present in RNA (Gain)")),
    y        = log2(TPM + 1)   # <-- THIS is the 'y' column
  )


df_segs <- df_plot %>%
  dplyr::select(GeneSymbol, status, genotype, y) %>%
  pivot_wider(names_from = genotype, values_from = y) %>%
  filter(!is.na(WT) & !is.na(KO)) %>%
  mutate(
    j    = rnorm(n(), mean = 0, sd = 0.12),   # shared jitter per gene
    x_wt = 1 + j,
    x_ko = 2 + j
  )

# Points data with same jitter
df_pts <- df_plot %>%
  left_join(df_segs %>% dplyr::select(GeneSymbol, status, j), by = c("GeneSymbol","status")) %>%
  mutate(xj = as.numeric(genotype) + j)

# --- 3) Plot: boxplots + paired points + arrowed segments
pdf("output/deseq2/boxplot_TPM-H3K27me3-WTvsKO_missing_vs_present_Gain_arrows.pdf",
    width = 5, height = 4.5)

ggplot() +
  # Boxplots (no jitter here; just the two categories)
  geom_boxplot(
    data = df_pts,
    aes(x = genotype, y = y, fill = genotype, color = genotype),
    outlier.shape = NA, width = 0.55, alpha = 0.25
  ) +
  # Arrows WT ->     aes(x = x_wt, xend = x_ko, y = WT, yend = KO),
# per gene (uses shared jitter so lines go through the points)
  geom_segment(
    data = df_segs,
    aes(x = x_wt, xend = x_ko, y = WT, yend = KO),
    inherit.aes = FALSE,
    alpha = 0.35, linewidth = 0.3,
    arrow = arrow(length = unit(2, "mm"), type = "closed")
  ) +
  # Points (one per gene per genotype)
  geom_point(
    data = df_pts,
    aes(x = xj, y = y, color = genotype),
    size = 1.6, alpha = 0.85
  ) +
  facet_wrap(~ status, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = c(WT = "black", KO = "red")) +
  scale_color_manual(values = c(WT = "black", KO = "red")) +
  labs(x = NULL, y = "TPM (log2 + 1)",
       title = "WT vs KO — Gain genes (Missing vs Present in RNA)\nBoxplots, per-gene dots, and WT→KO arrows") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dev.off()





pdf("output/deseq2/boxplot_TPM-H3K27me3-WTvsKO_missing_and_present_Gain_arrows.pdf",
    width = 2, height = 3)

ggplot() +
  # Boxplots (no jitter here; just the two categories)
  geom_boxplot(
    data = df_pts,
    aes(x = genotype, y = y, fill = genotype, color = genotype),
    outlier.shape = NA, width = 0.55, alpha = 0.25
  ) +
  # Arrows WT ->     aes(x = x_wt, xend = x_ko, y = WT, yend = KO),
# per gene (uses shared jitter so lines go through the points)
  geom_segment(
    data = df_segs,
    aes(x = x_wt, xend = x_ko, y = WT, yend = KO),
    inherit.aes = FALSE,
    alpha = 0.35, linewidth = 0.2,
    arrow = arrow(length = unit(2, "mm"), type = "closed")
  ) +
  # Points (one per gene per genotype)
  geom_point(
    data = df_pts,
    aes(x = xj, y = y, color = genotype),
    size = 0.5, alpha = 0.5
  ) +
  scale_fill_manual(values = c(WT = "black", KO = "red")) +
  scale_color_manual(values = c(WT = "black", KO = "red")) +
  labs(x = NULL, y = "TPM (log2 + 1)",
       title = "WT vs KO — Gain genes (Missing vs Present in RNA)\nBoxplots, per-gene dots, and WT→KO arrows") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dev.off()


########################








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














######################################
# GAIN H3K27me3 from `001*/018*` DESEQ2 MACS2 ###################
######################################
res <- read.table("output/deseq2/res_ESC_OEKO_vs_ESC_WT.txt",
                  sep = "\t", 
                  header = TRUE, 
                  row.names = 1, 
                  stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column(var = "gene") %>%
  as_tibble()



gene_Gain <- read.csv(
  "../018__CutRun_DOX_V1/output/edgeR/DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  # 1) keep significant peaks and assign effect direction
  filter(padj < 0.05) %>%
  mutate(effect_dir = case_when(
    log2FoldChange >=  0.58 ~ "pos",
    log2FoldChange <= -0.58 ~ "neg",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(effect_dir)) %>%
  # 2) determine gene-level direction
  group_by(geneSymbol) %>%
  summarise(
    n_pos = sum(effect_dir == "pos"),
    n_neg = sum(effect_dir == "neg"),
    direction = case_when(
      n_pos > 0 & n_neg > 0 ~ "Mix",
      n_pos > 0             ~ "Gain",
      n_neg > 0             ~ "Lost",
      TRUE                  ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  # 3) keep only Gain genes (exclude Mix)
  filter(direction == "Gain") %>%
  dplyr::select(geneSymbol) %>%
  unique() %>%
  dplyr::rename("GeneSymbol" = "geneSymbol")
  


res_Gain = res %>%
  inner_join(gene_Gain)



res_df <- res_Gain %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)



keyvals <- ifelse(
  res_Gain$log2FoldChange < -0.58 & res_Gain$padj < 5e-2, 'Sky Blue',
    ifelse(res_Gain$log2FoldChange > 0.58 & res_Gain$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.25)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.25)'

pdf("output/deseq2/plotVolcano_res_Gain-DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3-q05fc058_ESC_OEKO_vs_ESC_WT.pdf", width=4, height=5)    
EnhancedVolcano(res_Gain,
  lab = res_Gain$GeneSymbol,
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
  annotate("text", x = 3, y = 15, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 15, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()









## WHY SOME GENES ARE MISSING between gain and RNAseq ########################
gene_Gain <- read.csv(
  "../018__CutRun_DOX_V1/output/edgeR/DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  # 1) keep significant peaks and assign effect direction
  filter(padj < 0.05) %>%
  mutate(effect_dir = case_when(
    log2FoldChange >=  0.58 ~ "pos",
    log2FoldChange <= -0.58 ~ "neg",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(effect_dir)) %>%
  # 2) determine gene-level direction
  group_by(geneSymbol) %>%
  summarise(
    n_pos = sum(effect_dir == "pos"),
    n_neg = sum(effect_dir == "neg"),
    direction = case_when(
      n_pos > 0 & n_neg > 0 ~ "Mix",
      n_pos > 0             ~ "Gain",
      n_neg > 0             ~ "Lost",
      TRUE                  ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  # 3) keep only Gain genes (exclude Mix)
  filter(direction == "Gain") %>%
  dplyr::select(geneSymbol) %>%
  unique() %>%
  dplyr::rename("GeneSymbol" = "geneSymbol")

## Seprate gene presnet and absent after fusion RNAseq and CutRun
res_Gain_MISSED = gene_Gain %>%
  anti_join(res, by = "GeneSymbol")

res_Gain_PRESENT = gene_Gain %>%
  inner_join(res, by = "GeneSymbol")


## Load TPM info of each genes
tpm_long <- read.csv(
  "output/deseq2/txi_Kallisto-GeneLevel-TPM.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

avg_tpm <- tpm_long %>%
  group_by(GeneSymbol, genotype) %>%
  summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop")


mat <- avg_tpm %>%
  tidyr::pivot_wider(names_from = genotype, values_from = TPM) %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
  arrange(desc(WT)) %>%
  dplyr::select(GeneSymbol, WT, OE) %>%
  tibble::column_to_rownames("GeneSymbol") %>%
  as.matrix()

## 5) Replace any NAs with 0
mat[is.na(mat)] <- 0

## 6) Plot heatmap (log2+1 TPM)
# 0) Inputs: tpm_long, res_Gain_MISSED (GeneSymbol), res_Gain_PRESENT (GeneSymbol)
genes_missing <- unique(res_Gain_MISSED$GeneSymbol)
genes_present <- unique(res_Gain_PRESENT$GeneSymbol)

# 1) Helper to build WT/OE matrix for a gene list
make_mat <- function(genes) {
  avg_tpm <- tpm_long %>%
    dplyr::filter(GeneSymbol %in% genes, genotype %in% c("WT", "OE")) %>%
    group_by(GeneSymbol, genotype) %>%
    summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop")

  mat <- avg_tpm %>%
    pivot_wider(names_from = genotype, values_from = TPM) %>%
    mutate(WT = coalesce(WT, 0), OE = coalesce(OE, 0)) %>%
    filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
    arrange(desc(WT)) %>%
    dplyr::select(GeneSymbol, WT, OE) %>%
    tibble::column_to_rownames("GeneSymbol") %>%
    as.matrix()

  mat[is.na(mat)] <- 0
  mat
}

mat_missing <- make_mat(genes_missing)
mat_present <- make_mat(genes_present)

# 2) Shared color scale across both panels
vmin  <- min(log2(c(mat_missing, mat_present) + 1), na.rm = TRUE)
vmax  <- max(log2(c(mat_missing, mat_present) + 1), na.rm = TRUE)
ncols <- 100
cols  <- colorRampPalette(c("white","orange","red"))(ncols)
brks  <- seq(vmin, vmax, length.out = ncols + 1)

# 3) Draw side-by-side into one PDF
pdf("output/deseq2/heatmap_TPM-H3K27me3-WTvsOEKO_missing_vs_present_Gain.pdf",
    width = 6,
    height = max(4, max(nrow(mat_missing), nrow(mat_present)) * 0.18))

p1 <- pheatmap(
  log2(mat_missing + 1),
  color = cols, breaks = brks,
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "Missing in RNA (Gain genes)",
  border_color = NA, silent = TRUE
)

p2 <- pheatmap(
  log2(mat_present + 1),
  color = cols, breaks = brks,
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "Present in RNA (Gain genes)",
  border_color = NA, silent = TRUE
)

gridExtra::grid.arrange(grobs = list(p1$gtable, p2$gtable), ncol = 2)
dev.off()




# Helper: average TPM per gene per genotype (one dot per gene)
avg_by_gene <- function(genes, status_label) {
  tpm_long %>%
    filter(GeneSymbol %in% genes, genotype %in% c("WT","OE")) %>%
    group_by(GeneSymbol, genotype) %>%
    summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(status = status_label)
}

df_missing  <- avg_by_gene(genes_missing, "Missing in RNA (Gain)")
df_present  <- avg_by_gene(genes_present, "Present in RNA (Gain)")
df_plot <- bind_rows(df_missing, df_present) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT","OE")),
    status   = factor(status, levels = c("Missing in RNA (Gain)", "Present in RNA (Gain)"))
  )

# Plot
pdf("output/deseq2/boxplot_TPM-H3K27me3-WTvsOEKO_missing_vs_present_Gain.pdf",
    width = 5, height = 4.5)

ggplot(df_plot, aes(x = genotype, y = log2(TPM + 1), fill = genotype, color = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.25) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
             size = 1.6, alpha = 0.75) +
  facet_wrap(~ status, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = c(WT = "black", OE = "blue")) +
  scale_color_manual(values = c(WT = "black", OE = "blue")) +
  labs(x = NULL, y = "TPM (log2 + 1)", title = "WT vs OEKO — Gain genes (Missing vs Present in RNA)") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dev.off()










# --- 2) Create one jitter value PER GENE (per facet) and reuse for both WT & KO
# Segments data (one row per gene/facet), with WT/OE y-values and shared jitter

second <- "OE"          # <-- change to "OEKO" if that's your column value
levels2 <- c("WT", second)

# -- inputs assumed: tpm_long, res_Gain_MISSED (GeneSymbol), res_Gain_PRESENT (GeneSymbol)
genes_missing <- unique(res_Gain_MISSED$GeneSymbol)
genes_present <- unique(res_Gain_PRESENT$GeneSymbol)

avg_by_gene <- function(genes, status_label) {
  tpm_long %>%
    filter(GeneSymbol %in% genes, genotype %in% levels2) %>%
    group_by(GeneSymbol, genotype) %>%
    summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(status = status_label)
}


df_missing <- avg_by_gene(genes_missing, "Missing in RNA (Gain)")
df_present <- avg_by_gene(genes_present, "Present in RNA (Gain)")

df_plot <- bind_rows(df_missing, df_present) %>%
  mutate(
    genotype = factor(genotype, levels = levels2),
    status   = factor(status, levels = c("Missing in RNA (Gain)", "Present in RNA (Gain)")),
    y        = log2(TPM + 1)   # <-- THIS is the 'y' column
  )

# one jitter per gene (per facet), reused for both genotypes
df_segs <- df_plot %>%
  dplyr::select(GeneSymbol, status, genotype, y) %>%
  pivot_wider(names_from = genotype, values_from = y) %>%
  dplyr::filter(!is.na(.data[[levels2[1]]]) & !is.na(.data[[levels2[2]]])) %>%
  mutate(
    j    = rnorm(n(), 0, 0.12),
    x_1  = 1 + j,
    x_2  = 2 + j
  )

df_pts <- df_plot %>%
  left_join(df_segs %>% dplyr::select(GeneSymbol, status, j), by = c("GeneSymbol","status")) %>%
  mutate(xj = as.numeric(genotype) + j)




df_plot <- bind_rows(df_missing, df_present) %>%
  mutate(
    genotype = factor(genotype, levels = levels2),
    status   = factor(status, levels = c("Missing in RNA (Gain)", "Present in RNA (Gain)")),
    y        = log2(TPM + 1)   # <-- THIS is the 'y' column
  )


df_segs <- df_plot %>%
  dplyr::select(GeneSymbol, status, genotype, y) %>%
  pivot_wider(names_from = genotype, values_from = y) %>%
  filter(!is.na(WT) & !is.na(OE)) %>%
  mutate(
    j    = rnorm(n(), mean = 0, sd = 0.12),   # shared jitter per gene
    x_wt = 1 + j,
    x_oe = 2 + j
  )

# Points data with same jitter
df_pts <- df_plot %>%
  left_join(df_segs %>% dplyr::select(GeneSymbol, status, j), by = c("GeneSymbol","status")) %>%
  mutate(xj = as.numeric(genotype) + j)

# --- 3) Plot: boxplots + paired points + arrowed segments
pdf("output/deseq2/boxplot_TPM-H3K27me3-WTvsOEKO_missing_vs_present_Gain_arrows.pdf",
    width = 5, height = 4.5)

ggplot() +
  # Boxplots (no jitter here; just the two categories)
  geom_boxplot(
    data = df_pts,
    aes(x = genotype, y = y, fill = genotype, color = genotype),
    outlier.shape = NA, width = 0.55, alpha = 0.25
  ) +
  # Arrows WT -> OE per gene (uses shared jitter so lines go through the points)
  geom_segment(
    data = df_segs,
    aes(x = x_wt, xend = x_oe, y = WT, yend = OE),
    inherit.aes = FALSE,
    alpha = 0.35, linewidth = 0.3,
    arrow = arrow(length = unit(2, "mm"), type = "closed")
  ) +
  # Points (one per gene per genotype)
  geom_point(
    data = df_pts,
    aes(x = xj, y = y, color = genotype),
    size = 1.6, alpha = 0.85
  ) +
  facet_wrap(~ status, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = c(WT = "black", OE = "blue")) +
  scale_color_manual(values = c(WT = "black", OE = "blue")) +
  labs(x = NULL, y = "TPM (log2 + 1)",
       title = "WT vs OE — Gain genes (Missing vs Present in RNA)\nBoxplots, per-gene dots, and WT→KO arrows") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dev.off()



pdf("output/deseq2/boxplot_TPM-H3K27me3-OEWTvsKO_missing_and_present_Gain_arrows.pdf",
    width = 2, height = 3)

ggplot() +
  # Boxplots (no jitter here; just the two categories)
  geom_boxplot(
    data = df_pts,
    aes(x = genotype, y = y, fill = genotype, color = genotype),
    outlier.shape = NA, width = 0.55, alpha = 0.25
  ) +
  # Arrows WT ->     aes(x = x_wt, xend = x_oe, y = WT, yend = KO),
# per gene (uses shared jitter so lines go through the points)
  geom_segment(
    data = df_segs,
    aes(x = x_wt, xend = x_oe, y = WT, yend = OE),
    inherit.aes = FALSE,
    alpha = 0.35, linewidth = 0.2,
    arrow = arrow(length = unit(2, "mm"), type = "closed")
  ) +
  # Points (one per gene per genotype)
  geom_point(
    data = df_pts,
    aes(x = xj, y = y, color = genotype),
    size = 0.5, alpha = 0.5
  ) +
  scale_fill_manual(values = c(WT = "black", OE = "blue")) +
  scale_color_manual(values = c(WT = "black", OE = "blue")) +
  labs(x = NULL, y = "TPM (log2 + 1)",
       title = "WT vs OE — Gain genes (Missing vs Present in RNA)\nBoxplots, per-gene dots, and WT→OE arrows") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dev.off()




pdf("output/deseq2/boxplot_TPM-H3K27me3-OEWTvsKO_missing_and_present_Gain_arrows.pdf",
    width = 2, height = 3)

ggplot() +
  # Boxplots (no jitter here; just the two categories)
  geom_boxplot(
    data = df_pts,
    aes(x = genotype, y = y, fill = genotype, color = genotype),
    outlier.shape = NA, width = 0.55, alpha = 0.25
  ) +
  # Arrows WT ->     aes(x = x_wt, xend = x_oe, y = WT, yend = KO),
# per gene (uses shared jitter so lines go through the points)
  geom_segment(
    data = df_segs,
    aes(x = x_wt, xend = x_oe, y = WT, yend = OE),
    inherit.aes = FALSE,
    alpha = 0.35, linewidth = 0.2,
    arrow = arrow(length = unit(2, "mm"), type = "closed")
  ) +
  # Points (one per gene per genotype)
  geom_point(
    data = df_pts,
    aes(x = xj, y = y, color = genotype),
    size = 0.5, alpha = 0.5
  ) +
  scale_fill_manual(values = c(WT = "black", OE = "blue")) +
  scale_color_manual(values = c(WT = "black", OE = "blue")) +
  labs(x = NULL, y = "TPM (log2 + 1)",
       title = "WT vs OE — Gain genes (Missing vs Present in RNA)\nBoxplots, per-gene dots, and WT→OE arrows") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

dev.off()

########################



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







# Differential alternative mRNA splicing


## Run IsoformSwitchAnalyzeR prerequisets Part1 - webserver CPC2, PFAM, IUPRED2A, SIGNALP

After running `isoformSwitchAnalysisPart1()` from  `## IsoformSwitchAnalyzeR usage`


**Run these in Webserver** :
- coding potential with [CPC2](https://cpc2.gao-lab.org/), I put `output/IsoformSwitchAnalyzeR/isoformSwitchAnalyzeR_isoform_nt.fasta`;
  --> `IsoformSwitchAnalyzeR_kallisto/result_cpc2.txt`
- protein domain with [PFAM](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan), manually on the cluster, see below:
  --> OK
- Prediction of Intrinsically Unstructured Proteins with [IUPred2](https://iupred2a.elte.hu/); V3 do not support multi FASTA file with Default: `IUPred2 long disorder (default)`;
  --> `IsoformSwitchAnalyzeR_kallisto/isoformSwitchAnalyzeR_isoform_AA_complete.result`
- Prediction of signal peptide with [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) with option Eukarya/short output/fast; and save the *Prediction summary*
  - --> `IsoformSwitchAnalyzeR_kallisto/prediction_results_SignalIP6.txt` --> LOOK GOOD


**PFAM reformatting**:

```bash
# Load packages and modules
conda activate deseq2
module load HMMER


# Go to the software
cd ../../Master/software
cd pfam_scan

# Run it
./pfam_scan.py ../../../001_EZH1_Project/019__RNAseq_ESC_V1/output/IsoformSwitchAnalyzeR_kallisto/isoformSwitchAnalyzeR_isoform_AA_complete.fasta ../pfamdb/ -out ../../../001_EZH1_Project/019__RNAseq_ESC_V1/output/pfam/pfam_results_kallisto.txt


# Re-format the output
cd ../../../001_EZH1_Project/019__RNAseq_ESC_V1
python3 scripts/reformat_pfam_kallisto.py
```



## Run IsoformSwitchAnalyzeR Part1 and part2




```bash
conda activate IsoformSwitchAnalyzeRv5
```

```R
# packages
library("IsoformSwitchAnalyzeR")
library("rtracklayer")


set.seed(42)


# Kallisto ####################################################

# Importing the Data
salmonQuant <- importIsoformExpression(
    parentDir = "output/kallisto/")

# metadata file
myDesign = data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('.*_(WT|KO|OEKO)_.*', '\\1', colnames(salmonQuant$abundance)[-1])
)


# 
aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "../../Master/meta/gencode.v47.chr_patch_hapl_scaff.annotation.gtf", # gencode.v47.annotation.gtf gencode.v47.chr_patch_hapl_scaff.annotation.gtf
    isoformNtFasta       = "../../Master/meta/salmon/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    fixStringTieAnnotationProblem = TRUE,
    showProgress = FALSE
)
summary(aSwitchList)


SwitchList <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist   = aSwitchList,
    pathToOutput = 'output/IsoformSwitchAnalyzeR_kallisto',
    outputSequences      = TRUE, # change to TRUE whan analyzing your own data 
    prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)


#--> Run in WebServer the CPC2, PFAM, IUPRED2A, SIGNALP

analysSwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList, 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "output/IsoformSwitchAnalyzeR_kallisto/result_cpc2.txt",
  pathToPFAMresultFile      = "output/pfam/pfam_results_kallisto_reformat.txt",
  pathToIUPred2AresultFile  = "output/IsoformSwitchAnalyzeR_kallisto/isoformSwitchAnalyzeR_isoform_AA_complete.result",
  pathToSignalPresultFile   = "output/IsoformSwitchAnalyzeR_kallisto/prediction_results_SignalIP6.txt",
  outputPlots               = TRUE
)

## SAVE IMAGE R SESSION
# save.image("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto.RData")
# load("output/IsoformSwitchAnalyzeR_kallisto/IsoformSwitchAnalyzeR_kallisto.RData")
##


## Generate plot for a gene - All panels
LSAMP
### WT vs KO
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/switchPlot-WTKO-ECEL1.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(analysSwitchList, gene= "ECEL1", condition1= "WT", condition2= "KO", reverseMinus = FALSE)
dev.off()
### WT vs OEKO
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/switchPlot-WTOEKO-HTRA1.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(analysSwitchList, gene= "HTRA1", condition1= "WT", condition2= "OEKO", reverseMinus = FALSE)
dev.off()


## Only gene and isoform expression
### WT vs KO
pdf("output/IsoformSwitchAnalyzeR_kallisto/switchPlotGeneExp-WTKO-ECEL1.pdf", width=3, height=5)
switchPlotGeneExp(
  switchAnalyzeRlist = analysSwitchList,
  gene = "ECEL1",
  condition1 = "WT",
  condition2 = "KO"
)  
dev.off()
pdf("output/IsoformSwitchAnalyzeR_kallisto/switchPlotIsoExp-WTKO-ECEL1.pdf", width=5, height=5)
switchPlotIsoExp(
  switchAnalyzeRlist = analysSwitchList,
  gene = "ECEL1",
  condition1 = "WT",
  condition2 = "KO"
) +
  scale_fill_manual(values = c("WT" = "black", "KO" = "red")) +
  scale_color_manual(values = c("WT" = "black", "KO" = "red"))
dev.off()


### WT vs OEKO

pdf("output/IsoformSwitchAnalyzeR_kallisto/switchPlotGeneExp-WTOEKO-HTRA1.pdf", width=3, height=5)
switchPlotGeneExp(
  switchAnalyzeRlist = analysSwitchList,
  gene = "HTRA1",
  condition1 = "WT",
  condition2 = "OEKO"
) 
dev.off()
pdf("output/IsoformSwitchAnalyzeR_kallisto/switchPlotIsoExp-WTOEKO-HTRA1.pdf", width=5, height=5)
switchPlotIsoExp(
  switchAnalyzeRlist = analysSwitchList,
  gene = "HTRA1",
  condition1 = "WT",
  condition2 = "OEKO"
) +
  scale_fill_manual(values = c("WT" = "black", "OEKO" = "blue")) +
  scale_color_manual(values = c("WT" = "black", "OEKO" = "blue"))
dev.off()






# Genome-wide Summaries
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/extractSwitchOverlap.pdf', onefile = FALSE, height=6, width = 9)
extractSwitchOverlap(
    analysSwitchList,
    filterForConsequences=TRUE,
    plotIsoforms = FALSE
)
dev.off()


pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/extractConsequenceSummary.pdf', onefile = FALSE, height=6, width = 9)
extractConsequenceSummary(
    analysSwitchList,
    consequencesToAnalyze='all',
    plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
dev.off()
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/extractConsequenceSummary_Genes.pdf', onefile = FALSE, height=6, width = 9)
extractConsequenceSummary(
    analysSwitchList,
    consequencesToAnalyze='all',
    plotGenes = TRUE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
dev.off()


# Consequence Enrichment Analysis
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/extractConsequenceEnrichment.pdf', onefile = FALSE, height=6, width = 12)
extractConsequenceEnrichment(
    analysSwitchList,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()



# Splicing Enrichment Analysis

pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/extractSplicingEnrichment.pdf', onefile = FALSE, height=6, width = 12)
extractSplicingEnrichment(
    analysSwitchList,
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()

#Overview Plots
## Volcano like plot

pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/Overview_Plots.pdf', onefile = FALSE, height=3, width = 6)
ggplot(data=analysSwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_1) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off()

pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/Overview_Plots3.pdf', onefile = FALSE, height=6, width = 6)
ggplot(data=analysSwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off()

## count nb of isoforms:
significant_isoforms <- analysSwitchList$isoformFeatures %>%
  filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05)
###  Count for WT vs KO
WT_vs_KO <- significant_isoforms %>%
  filter(condition_1 == "KO" & condition_2 == "WT") %>%
  nrow()
### Count for WT vs OEKO
WT_vs_OEKO <- significant_isoforms %>%
  filter(condition_1 == "OEKO" & condition_2 == "WT") %>%
  nrow()
### Print the counts
cat("Number of significant isoform switches (WT vs KO):", WT_vs_KO, "\n")
cat("Number of significant isoform switches (WT vs OEKO):", WT_vs_OEKO, "\n")

## Import GTF to have gene name chromosome 
gtf <- import("../../Master/meta/gencode.v47.chr_patch_hapl_scaff.annotation.gtf")
## Convert to data frame
gtf_df <- as.data.frame(gtf)
## Keep only gene entries
gene_df <- gtf_df %>%
  filter(type == "gene") %>%
  dplyr::select(seqnames, gene_name) %>%
  distinct() %>%
  rename(chromosome = seqnames) %>%
  as_tibble()

# Save output list of genes
## Signif only without X chr
significant_isoforms__KO_geneSymbol =  significant_isoforms %>%
  filter(condition_1 == "KO" & condition_2 == "WT") %>%
  dplyr::select(gene_name) %>% unique() %>% left_join(gene_df) %>% as_tibble() %>% dplyr::filter(chromosome != "chrX") %>%   dplyr::select(gene_name)
write.table(
  significant_isoforms__KO_geneSymbol,
  file = "output/IsoformSwitchAnalyzeR_kallisto/significant_isoforms_dIF01qval05__KO_geneSymbol.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
## Signif only without X chr AND with consequence
significant_isoforms__KO_geneSymbol =  significant_isoforms %>%
  filter(condition_1 == "KO" & condition_2 == "WT" & switchConsequencesGene == TRUE) %>%
  dplyr::select(gene_name) %>% unique() %>% left_join(gene_df) %>% as_tibble() %>% dplyr::filter(chromosome != "chrX") %>%   dplyr::select(gene_name)
write.table(
  significant_isoforms__KO_geneSymbol,
  file = "output/IsoformSwitchAnalyzeR_kallisto/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE__KO_geneSymbol.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


## Signif only without X chr
significant_isoforms__OEKO_geneSymbol =  significant_isoforms %>%
  filter(condition_1 == "OEKO" & condition_2 == "WT") %>%
  dplyr::select(gene_name) %>% unique() %>% left_join(gene_df) %>% as_tibble() %>% dplyr::filter(chromosome != "chrX")  %>%   dplyr::select(gene_name)
write.table(
  significant_isoforms__OEKO_geneSymbol,
  file = "output/IsoformSwitchAnalyzeR_kallisto/significant_isoforms_dIF01qval05__OEKO_geneSymbol.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
## Signif only without X chr AND with consequence
significant_isoforms__OEKOKO_geneSymbol =  significant_isoforms %>%
  filter(condition_1 == "OEKO" & condition_2 == "WT" & switchConsequencesGene == TRUE) %>%
  dplyr::select(gene_name) %>% unique() %>% left_join(gene_df) %>% as_tibble() %>% dplyr::filter(chromosome != "chrX") %>%   dplyr::select(gene_name)
write.table(
  significant_isoforms__OEKOKO_geneSymbol,
  file = "output/IsoformSwitchAnalyzeR_kallisto/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE__OEKO_geneSymbol.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)




### Switch vs Gene changes:
pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/Overview_Plots2.pdf', onefile = FALSE, height=3, width = 6)
ggplot(data=analysSwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) + 
    facet_wrap(~ condition_1) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw()
dev.off()






```


--> There should be ways to investigate the list of isoforms more carefully but I m just gonna do a GO..
  --> **Gene list of interest**: Gene with isforom switch and a functional consequence (loss of coding potential... etc..) `output/IsoformSwitchAnalyzeR_kallisto/significant_isoforms_dIF01qval05switchConsequencesGeneTRUE__[OEKO or KO]_geneSymbol.txt`
    --> Check in `001*/018*` whether these genes are bound with EZH1 (Venn overlap using consens peak does ot show high overlap)


--> **IsoformSwitchAnalyzeR evaluates changes within each isoform, so the X chromosome doesn’t bias the overall analysis**. I can simply exclude isoform switches detected on chrX.



More info [here](https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html)



# Investigate EZH1 KO deletion



Let's confirm that EZH1 KO deletion 8bp in exon7 has occured.




## Using BCF


### bcftools installation

Follow guidelines from [this](https://samtools.github.io/bcftools/howtos/install.html)

```bash

cd ../../Master/software
mkdir bcftools

cd bcftools

# download bcf tools
git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/samtools/bcftools.git
cd bcftools

make


# Add shortcut to launch bcftools

nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software/bcftools/bcftools
# Restart terminal


module load SAMtools
```

--> TO USE bcftools, need first to run  `module load SAMtools` asnd then type `bcftools`


### bcftools usage



```bash
module load SAMtools


# Quick investagion
# (A) Do your BAMs actually have "chr17"?
samtools idxstats output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam | head

# (B) What’s the coverage in/near that exon?
samtools depth -r chr17:42720200-42720520 output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam \
  | awk '{sum+=$3; n++} END{print "avg_depth:", (n?sum/n:0), "positions:", n}'

#--> YEs, there is enough coverage


# (C) Do any reads show a deletion (look for D in CIGAR) ?
samtools view output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam chr17:42720200-42720520 \
  | awk '$6 ~ /[0-9]+D/ {print $3,$4,$6}' | head
  samtools view output/STAR/fastp/ESC_WT_R2_Aligned.sortedByCoord.out.bam chr17:42720200-42720520 \
  | awk '$6 ~ /[0-9]+D/ {print $3,$4,$6}' | head
  samtools view output/STAR/fastp/ESC_WT_R3_Aligned.sortedByCoord.out.bam chr17:42720200-42720520 \
  | awk '$6 ~ /[0-9]+D/ {print $3,$4,$6}' | head

samtools view output/STAR/fastp/ESC_KO_R1_Aligned.sortedByCoord.out.bam chr17:42720200-42720520 \
  | awk '$6 ~ /[0-9]+D/ {print $3,$4,$6}' | head
samtools view output/STAR/fastp/ESC_KO_R2_Aligned.sortedByCoord.out.bam chr17:42720200-42720520 \
  | awk '$6 ~ /[0-9]+D/ {print $3,$4,$6}' | head
samtools view output/STAR/fastp/ESC_KO_R3_Aligned.sortedByCoord.out.bam chr17:42720200-42720520 \
  | awk '$6 ~ /[0-9]+D/ {print $3,$4,$6}' | head

#--> YES, we see 16D in the CIGAR, meaning there is 16bp deletion detected in KO samples!!!





# Quantification 

####################################################
## 001/019 samples #################################
####################################################


REGION="chr17:42720200-42720520"   # your exon7 window

printf "sample\ttotal_reads\tdel16_reads\tfraction\n" > output/variant_calling/exon7_16D_cigar_counts.tsv
for bam in output/STAR/fastp/*_Aligned.sortedByCoord.out.bam; do
  sample=$(basename "$bam" _Aligned.sortedByCoord.out.bam)

  # total reads overlapping region (primary, properly mapped-ish; tweak filters if needed)
  total=$(samtools view -F 0x904 -q 20 "$bam" "$REGION" | wc -l)

  # reads whose CIGAR contains a 16-nt deletion anywhere within the region
  del16=$(samtools view -F 0x904 -q 20 "$bam" "$REGION" \
          | awk '$6 ~ /16D/ {c++} END{print c+0}')

  frac=0; if [ "$total" -gt 0 ]; then frac=$(awk -v a="$del16" -v b="$total" 'BEGIN{printf "%.4f", a/b}'); fi
  printf "%s\t%s\t%s\t%s\n" "$sample" "$total" "$del16" "$frac" >> output/variant_calling/exon7_16D_cigar_counts.tsv
done
column -t output/variant_calling/exon7_16D_cigar_counts.tsv

#--> This show total read in the window and read with the deletion; but it include reads that does not overlap the deletion



####################################################
## 001/001 samples #################################
####################################################


REGION="chr17:42720200-42720520"   # your exon7 window

printf "sample\ttotal_reads\tdel8_reads\tfraction\n" > output/variant_calling/exon7_8D_cigar_counts-001001.tsv
for bam in ../001__RNAseq/output/STAR_hg38/ESC*_Aligned.sortedByCoord.out.bam; do
  sample=$(basename "$bam" _Aligned.sortedByCoord.out.bam)

  # total reads overlapping region (primary, properly mapped-ish; tweak filters if needed)
  total=$(samtools view -F 0x904 -q 20 "$bam" "$REGION" | wc -l)

  # reads whose CIGAR contains a 16-nt deletion anywhere within the region
  del16=$(samtools view -F 0x904 -q 20 "$bam" "$REGION" \
          | awk '$6 ~ /8D/ {c++} END{print c+0}')

  frac=0; if [ "$total" -gt 0 ]; then frac=$(awk -v a="$del16" -v b="$total" 'BEGIN{printf "%.4f", a/b}'); fi
  printf "%s\t%s\t%s\t%s\n" "$sample" "$total" "$del16" "$frac" >> output/variant_calling/exon7_8D_cigar_counts-001001.tsv
done
column -t output/variant_calling/exon7_8D_cigar_counts-001001.tsv


#--> This show total read in the window and read with the deletion; but it include reads that does not overlap the deletion





###################################################################################################################
# Only include reads that overlap with the deletion ###############################################################
###################################################################################################################

####################################################
## 001/019 samples #################################
####################################################

# --- set your event ---
CHR=chr17
DEL_START=42720375     # exact start of the deletion (leftmost deleted base, 1-based)
LEN=16                 # use 8 for your 8-bp sample
PAD=5

OUT=output/variant_calling/exon7_${LEN}D_exact.tsv
printf "sample\tspan_reads\tdel${LEN}_reads\tfraction\n" > "$OUT"

for bam in output/STAR/fastp/*_Aligned.sortedByCoord.out.bam; do
  s=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
  samtools view -F 0x904 -q 20 "$bam" "${CHR}:$((DEL_START-PAD))-$((DEL_START+PAD))" \
  | awk -v S="$DEL_START" -v L="$LEN" -v sample="$s" '
    function covers_start(cigar,pos,   re,tok,len,op){
      re="[0-9]+[MIDNSHP=X]"
      while (match(cigar,re)){
        tok=substr(cigar,RSTART,RLENGTH)
        len=substr(tok,1,length(tok)-1)+0
        op=substr(tok,length(tok),1)
        # count reads spanning S by M/=/X OR by a deletion that covers S
        if ((op=="M"||op=="="||op=="X") && S>=pos && S<pos+len) return 1
        if (op=="D" && S>=pos && S<pos+len) return 1
        if (op=="M"||op=="="||op=="X"||op=="D"||op=="N") pos+=len
        cigar=substr(cigar,RSTART+RLENGTH)
      } return 0
    }
    function has_exact_del(cigar,pos,   re,tok,len,op){
      re="[0-9]+[MIDNSHP=X]"
      while (match(cigar,re)){
        tok=substr(cigar,RSTART,RLENGTH)
        len=substr(tok,1,length(tok)-1)+0
        op=substr(tok,length(tok),1)
        if (op=="D" && len==L && pos==S) return 1
        if (op=="M"||op=="="||op=="X"||op=="D"||op=="N") pos+=len
        cigar=substr(cigar,RSTART+RLENGTH)
      } return 0
    }
    {
      q=$1; pos=$4; cig=$6
      if (covers_start(cig,pos)) seen[q]=1
      if (has_exact_del(cig,pos)) del[q]=1
    }
    END{
      for (q in seen) tot++
      for (q in del)  alt++
      frac=(tot?alt/tot:0)
      printf "%s\t%d\t%d\t%.4f\n", sample, tot, alt, frac
    }' >> "$OUT"
done

column -t "$OUT"






####################################################
## 001/015 samples #################################
####################################################

# --- set your event ---
CHR=chr17
DEL_START=42720375     # exact start of the deletion (leftmost deleted base, 1-based)
LEN=16                 # use 8 for your 8-bp sample
PAD=5

OUT=output/variant_calling/exon7_${LEN}D_exact-001015.tsv
printf "sample\tspan_reads\tdel${LEN}_reads\tfraction\n" > "$OUT"

for bam in ../015__RNAseq_PSC/output/STAR/fastp/*_Aligned.sortedByCoord.out.bam; do
  s=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
  samtools view -F 0x904 -q 20 "$bam" "${CHR}:$((DEL_START-PAD))-$((DEL_START+PAD))" \
  | awk -v S="$DEL_START" -v L="$LEN" -v sample="$s" '
    function covers_start(cigar,pos,   re,tok,len,op){
      re="[0-9]+[MIDNSHP=X]"
      while (match(cigar,re)){
        tok=substr(cigar,RSTART,RLENGTH)
        len=substr(tok,1,length(tok)-1)+0
        op=substr(tok,length(tok),1)
        # count reads spanning S by M/=/X OR by a deletion that covers S
        if ((op=="M"||op=="="||op=="X") && S>=pos && S<pos+len) return 1
        if (op=="D" && S>=pos && S<pos+len) return 1
        if (op=="M"||op=="="||op=="X"||op=="D"||op=="N") pos+=len
        cigar=substr(cigar,RSTART+RLENGTH)
      } return 0
    }
    function has_exact_del(cigar,pos,   re,tok,len,op){
      re="[0-9]+[MIDNSHP=X]"
      while (match(cigar,re)){
        tok=substr(cigar,RSTART,RLENGTH)
        len=substr(tok,1,length(tok)-1)+0
        op=substr(tok,length(tok),1)
        if (op=="D" && len==L && pos==S) return 1
        if (op=="M"||op=="="||op=="X"||op=="D"||op=="N") pos+=len
        cigar=substr(cigar,RSTART+RLENGTH)
      } return 0
    }
    {
      q=$1; pos=$4; cig=$6
      if (covers_start(cig,pos)) seen[q]=1
      if (has_exact_del(cig,pos)) del[q]=1
    }
    END{
      for (q in seen) tot++
      for (q in del)  alt++
      frac=(tot?alt/tot:0)
      printf "%s\t%d\t%d\t%.4f\n", sample, tot, alt, frac
    }' >> "$OUT"
done

column -t "$OUT"





####################################################
## 001/001 samples #################################
####################################################

CHR=chr17
DEL_START=42720376   # <-- exact start of the 16-bp deletion from IGV 42720374 42720391
LEN=8               
PAD=2                # tiny fetch padding


OUT=output/variant_calling/exon7_${LEN}D_exact-001001.tsv
printf "sample\tspan_reads\tdel${LEN}_reads\tfraction\n" > "$OUT"

for bam in ../001__RNAseq/output/STAR_hg38/ESC*_Aligned.sortedByCoord.out.bam; do
  s=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
  samtools view -F 0x904 -q 20 "$bam" "${CHR}:$((DEL_START-PAD))-$((DEL_START+PAD))" \
  | awk -v S="$DEL_START" -v L="$LEN" -v sample="$s" '
    function covers_start(cigar,pos,   re,tok,len,op){
      re="[0-9]+[MIDNSHP=X]"
      while (match(cigar,re)){
        tok=substr(cigar,RSTART,RLENGTH)
        len=substr(tok,1,length(tok)-1)+0
        op=substr(tok,length(tok),1)
        # count reads spanning S by M/=/X OR by a deletion that covers S
        if ((op=="M"||op=="="||op=="X") && S>=pos && S<pos+len) return 1
        if (op=="D" && S>=pos && S<pos+len) return 1
        if (op=="M"||op=="="||op=="X"||op=="D"||op=="N") pos+=len
        cigar=substr(cigar,RSTART+RLENGTH)
      } return 0
    }
    function has_exact_del(cigar,pos,   re,tok,len,op){
      re="[0-9]+[MIDNSHP=X]"
      while (match(cigar,re)){
        tok=substr(cigar,RSTART,RLENGTH)
        len=substr(tok,1,length(tok)-1)+0
        op=substr(tok,length(tok),1)
        if (op=="D" && len==L && pos==S) return 1
        if (op=="M"||op=="="||op=="X"||op=="D"||op=="N") pos+=len
        cigar=substr(cigar,RSTART+RLENGTH)
      } return 0
    }
    {
      q=$1; pos=$4; cig=$6
      if (covers_start(cig,pos)) seen[q]=1
      if (has_exact_del(cig,pos)) del[q]=1
    }
    END{
      for (q in seen) tot++
      for (q in del)  alt++
      frac=(tot?alt/tot:0)
      printf "%s\t%d\t%d\t%.4f\n", sample, tot, alt, frac
    }' >> "$OUT"
done

column -t "$OUT"







```


--




