# Project

**H9 cell lines**
- PSC native:
    - WT: 3 Bio Rep (A1-3)
    - KO: 3 Bio Rep (A4-6)
    - KOEF1aEZH1: 3 Bio Rep (A7-9)


--> Directional mRNA library preparation (poly A enrichment),NovaSeq X Plus Series (PE150)




**Objectives:**
- Put together with CutRun







# Pipeline
- Download data (wget)
- Rename files
- FastQC (fastqc)
- Trimming (fastp)
- Histone content (R)
- Mapping (bowtie2)
- Spike-in scaling (DiffBind)
- Bigwig generation (deepTools)
- peak calling (MACS2)
- peak assignment to gene (ChIPseeker)

--> Detail of the overall pipeline in `Meeting_20230919_draft.xlsx` 

# Download / import data


```bash
# Following email instructions
wget -r -b -c --user=X202SC24100552-Z01-F001 --password=f40sgxjy ftp://usftp21.novogene.com:21/


# Copy all .fz.gz data into input_raw/ folder
rsync -av --include '*/' --include '*.fq.gz' --exclude '*' usftp21.novogene.com/01.RawData/ input_raw/ # copy from usftp21 folder to input_raw
find input_raw/ -mindepth 2 -type f -exec mv -t input_raw/ {} + # mv files from their folder to input_raw/ folder
find input_raw/ -type d -empty -delete # delete empty directory

```

--> All good, files created in `usftp21.novogene.com/`




# Rename file

Renamed manually as only 8 samples



```bash
cp input_raw_Novogene/*.gz input/
```

--> All good 



# Fastp cleaning

```bash
sbatch scripts/fastp.sh # 28206738 ok
```



## mapping fastp trim

```bash
sbatch --dependency=afterany:28206738 scripts/STAR_mapping_fastp.sh # 28206965 ok; except KO_R2 fail.

# KO_R2 in interactive
## Re-generate fastp file
x=(
    "PSC_KO_R2"
    )

for x in "${x[@]}"; do
    fastp -i input_raw/${x}_1.fq.gz -I input_raw/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done

## Re-Run STAR
module load STAR/2.7.3a-GCC-9.3.0
module load SAMtools/1.16.1-GCC-11.3.0

x=(
    "PSC_KO_R2"
)

for x in "${x[@]}"; do
	STAR --genomeDir ../../Master/meta/STAR_hg38/ \
		--runThreadN 6 \
		--readFilesCommand zcat \
		--readFilesIn output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix output/STAR/fastp/${x}_
    samtools index output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done
# --> All good
```

--> Around 97% uniq aligned reads

--> KO_R2 bug: After *started mapping*: `ReadAlignChunk_processChunks.cpp:177:processChunks EXITING because of FATAL ERROR in input reads: unknown file format: the read ID should start with @ or >` 
  --> After inspection of the line with : `zcat output/fastp/PSC_KO_R2_1.fq.gz | awk 'NR % 4 == 1 && $0 !~ /^@/' > problematic_lines.txt` it also output : `gzip: invalid compressed data` indicating the **PSC_KO_R2_1.fq.gz** is corrupted: let's re-generate fastp files in interactive,.
  --> All good, fastp bugged, re-run fastp and STAR and that worked!



# Count with featureCounts





Count on gene features with parameter
```bash
conda activate featurecounts

# slight test
## -s for stranded
featureCounts -p -C -O -M --fraction -s 2 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf \
	-o output/test/PSC_WT_R1.txt output/STAR/fastp/PSC_WT_R1_Aligned.sortedByCoord.out.bam



# all samples:
sbatch scripts/featurecounts.sh # 28754176 ok 



```
test with `PSC_WT_R1`
--> Around 902% of succesfully assigned alignments with `-p -C -O -M --fraction -s 2` parameters...
--> all good!!

--> All sampels ~90% uniq aligned reads

## Calculate TPM and RPKM


Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
conda activate deseq2
# Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_TPM.sh # 28772752 xxx
# mv all output to output/tpm or rpkm folder
mv output/featurecounts/*tpm* output/tpm/
mv output/featurecounts/*rpkm* output/rpkm/
```

--> All good. 


# Shiny app


generate the `tpm_all_sample.txt` file and then go into `001_EZH1*/001__RNAseq` to create shiny app V2 including these

```R
library("tidyverse")
library("biomaRt")

# Display some genes in TPM: 
# ---> The code below is not perfect; issue at the geneSymbol conversin; to troubleshoot later; but it work
#### Generate TPM for ALL samples
#### collect all samples ID
samples <- c(  "PSC_WT_R1", "PSC_WT_R2", "PSC_WT_R3", "PSC_KO_R1", "PSC_KO_R2", "PSC_KO_R3", "PSC_KOEF1aEZH1_R1", "PSC_KOEF1aEZH1_R2",  "PSC_KOEF1aEZH1_R3")

## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/tpm/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR.")) %>%
    rename(!!sample := starts_with("output.STAR."))
}

## Merge all dataframe into a single one
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")
write.csv(tpm_all_sample, file="output/tpm/tpm_all_sample_Akoto.txt")
### If need to import: tpm_all_sample <- read_csv("output/tpm/tpm_all_sample_Akoto.txt") %>% dplyr::select(-"...1") #To import



## add geneSymbol
tpm_all_sample$Geneid <- gsub("\\..*", "", tpm_all_sample$Geneid) # remove Ensembl gene id version

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## Convert Ensembl gene IDs to gene symbols
tpm_all_genesymbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = tpm_all_sample$Geneid,
                     mart = ensembl)

tpm_all_withGenesymbols = tpm_all_sample %>% 
  dplyr::rename("ensembl_gene_id" = "Geneid") %>%
  left_join(tpm_all_genesymbols)


## Convert data from wide to long keep only geneSymbol
long_data <- tpm_all_withGenesymbols %>% 
  dplyr::select(-ensembl_gene_id) %>%
  drop_na() %>%
  tidyr::pivot_longer(-external_gene_name, names_to = "condition", values_to = "TPM") %>% 
  tidyr::separate(condition, into = c("Tissue", "Genotype", "Replicate"), sep = "_")

long_data_log2tpm = long_data %>%
  mutate(log2tpm = log2(TPM + 1))
## Save
write.table(long_data_log2tpm, file = c("output/tpm/long_data_log2tpm_Akoto.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


```


--> Next here: `001_EZH1*/001__RNAseq` (`## Shiny app V3; including Ciceri RNAseq neuron diff dataset + Akoto 001/015 RNAseq`)


# bigwig 
## Generate Bigwig coverage files

Let's generate **TPM coverage**:

```bash
conda activate deeptools
# run time-per-time:
sbatch scripts/TPM_bw.sh # 28775122 ok
```


Let's merge the bigwig into 1 file with wiggletools (will do average of bigwig signal and not sum, many options see [github](https://github.com/Ensembl/WiggleTools)):



**Run wiggletools:**
```bash
conda activate BedToBigwig
sbatch scripts/bigwigmerge_TPM.sh # 28990349 ok
```
*NOTE: bigwig are merge into 1 bedgraph which is then converted into 1 bigwig (wiggletools cannot output bigwig directly so need to pass by bedgraph or wiggle in between)*







## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all.sh # 28990409 ok

# Plot all
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 KOEF1aEZH1_R1 KOEF1aEZH1_R2 KOEF1aEZH1_R3 \
    --colors black black black red red red blue blue blue \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 KOEF1aEZH1_R1 KOEF1aEZH1_R2 KOEF1aEZH1_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf

```




# DEGs with deseq2 (and CutRun integration)

**IMPORTANT NOTE: Here it is advisable to REMOVE all genes from chromosome X and Y BEFORE doing the DEGs analysis (X chromosome re-activation occurs in some samples, notably these with more cell passage; in our case, the HET and KO)**
--> It is good to do this on the count matrix see [here](https://support.bioconductor.org/p/119932/)
### 'one-by-one' comparison
Comparison WT vs mutant:
- PSC KO vs WT
- PSC HET vs WT



### PSC KO vs WT

```bash
conda activate deseq2
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")


# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("PSC_WT_R1", "PSC_WT_R2" ,"PSC_WT_R3" ,"PSC_KO_R1" ,"PSC_KO_R2", "PSC_KO_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  dplyr::select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")




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

pdf("output/deseq2/plotVolcano_res_q05fc025_PSC_KO_vs_PSC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC',
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
write.table(res, file = "output/deseq2/res_PSC_KO_vs_PSC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.25 & res$padj < 5e-2, ]
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.25 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_PSC_KO_vs_PSC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_PSC_KO_vs_PSC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





## padj 0.05 FC 0.55 ##################################
res_df <- res %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.5 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.5 & res_df$padj == TRUE, na.rm = TRUE)

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.5)'

pdf("output/deseq2/plotVolcano_res_q05fc05_PSC_KO_vs_PSC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
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
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc05_PSC_KO_vs_PSC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc05_PSC_KO_vs_PSC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```






### PSC KOEF1aEZH1 vs WT

```bash
conda activate deseq2
```
Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")
library("org.Hs.eg.db")
library("biomaRt")

library("RColorBrewer")
library("pheatmap")
library("AnnotationDbi")


# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("PSC_WT_R1", "PSC_WT_R2" ,"PSC_WT_R3" ,"PSC_KOEF1aEZH1_R1" ,"PSC_KOEF1aEZH1_R2", "PSC_KOEF1aEZH1_R3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL

# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Not including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  dplyr::select(-replicate) %>%
  bind_cols(data.frame(samples))
### Including replicate
coldata_raw <- data.frame(samples) %>%
  separate(samples, into = c("time", "genotype", "replicate"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix),
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_WT", type="apeglm")




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

pdf("output/deseq2/plotVolcano_res_q05fc025_PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs WT, PSC',
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
write.table(res, file = "output/deseq2/res_PSC_KOEF1aEZH1_vs_PSC_WT.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.25 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_PSC_KOEF1aEZH1_vs_PSC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_PSC_KOEF1aEZH1_vs_PSC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





## padj 0.05 FC 0.55 ##################################
res_df <- res %>% as.data.frame() %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.5 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.5 & res_df$padj == TRUE, na.rm = TRUE)

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Hs.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.5 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.5 & res$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.5)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.5)'

pdf("output/deseq2/plotVolcano_res_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, PSC',
  pCutoff = 5e-2,         #
  FCcutoff = 0.5,
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
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.5 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.5 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```


- NOTE: *gene with padj of NA* is normal; happen for genes with outlier counts, having low mean expression values, and having all 0 counts; explain [here](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA)




## Generate GTF files of DEGs genes


```bash
# Generate gtf file from gene list:
output/deseq2/upregulated_q05fc05_PSC_KO_vs_PSC_WT.txt
output/deseq2/downregulated_q05fc05_PSC_KO_vs_PSC_WT.txt

output/deseq2/upregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.txt
output/deseq2/downregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.txt

### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/deseq2/upregulated_q05fc05_PSC_KO_vs_PSC_WT.txt > output/deseq2/upregulated_q05fc05_PSC_KO_vs_PSC_WT_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/deseq2/downregulated_q05fc05_PSC_KO_vs_PSC_WT.txt > output/deseq2/downregulated_q05fc05_PSC_KO_vs_PSC_WT_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/deseq2/upregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.txt > output/deseq2/upregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/deseq2/downregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.txt > output/deseq2/downregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT_as_gtf_geneSymbol.txt

## Filter the gtf
grep -Ff output/deseq2/upregulated_q05fc05_PSC_KO_vs_PSC_WT_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_upregulated_q05fc05_PSC_KO_vs_PSC_WT.gtf
grep -Ff output/deseq2/downregulated_q05fc05_PSC_KO_vs_PSC_WT_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_downregulated_q05fc05_PSC_KO_vs_PSC_WT.gtf

grep -Ff output/deseq2/upregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_upregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.gtf
grep -Ff output/deseq2/downregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_downregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_WT.gtf

```

--> Will be used to generate deepTool plots at `001*/016*` integration to check RNA fit well with CutRun




