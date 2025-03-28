# Project

Function of AMP2 in mice.
- 3 organs: CB cerebellum, CT cortex, HP hippocampus

# Reference genome

Reference genome and annotation downloaded from [ENCODE](https://www.encodeproject.org/data-standards/reference-sequences/):
- mm10 XY reference genome (ENCODE3 used only one reference genome for analysis)
- mm10 GENCODE VM21 merged annotations gtf file 

--> Genome Files in `/Akizu_Lab/Master/meta_mice`

# FASTQC on raw

```bash
sbatch scripts/fastqc_raw_CB.sh # 4611595 ok
sbatch scripts/fastqc_raw_CT.sh # 4611594 ok
sbatch scripts/fastqc_raw_HP.sh # 4611593 ok
```

--> QC of the raw pretty good! No overrepresented sequences; would worth mapp them in addition to the fastp-trimmed one

# Quality control with FASTP (trim)

```bash
sbatch scripts/fastp_CB.sh # 4611647 ok
sbatch scripts/fastp_CT.sh # 4611648 ok
sbatch scripts/fastp_HP.sh # 4611649 ok
```

Run fastqc on fastp-trimmed files

```bash
sbatch scripts/fastqc_fastp_CB.sh # 4626730 ok
sbatch scripts/fastqc_fastp_CT.sh # 4626731 ok
sbatch scripts/fastqc_fastp_HP.sh # 4626732 ok
```

--> QC loos good

# Mapping with STAR

## Genome indexation

```bash
sbatch scripts/STAR_index_mm10.sh # 4626872 ok
```

## mapping raw without trimming

```bash
sbatch --dependency=afterany:4626872 scripts/STAR_mapping_raw_CB.sh # 4626915 ok
sbatch --dependency=afterany:4626872 scripts/STAR_mapping_raw_CT.sh # 4626916 ok
sbatch --dependency=afterany:4626872 scripts/STAR_mapping_raw_HP.sh # 4626917 ok
```

--> Let's compil the number of uniquely mapped reads for all files (add it in the Google Drive `RNAseq_infos.xlsx` file)

```bash
# Print nb of uniq map reads for raw mapping
for file in output/STAR/raw/*Log.final.out; do
    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $file | awk '{print $NF}')
    echo "$file: Number of uniquely mapped reads: $uniquely_mapped_reads"
done > output/STAR/raw/uniq_map_reads_counts.txt
```


## mapping fastp trim

```bash
sbatch --dependency=afterany:4626872 scripts/STAR_mapping_fastp_CB.sh # 4626931 ok
sbatch --dependency=afterany:4626872 scripts/STAR_mapping_fastp_CT.sh # 4626932 ok
sbatch --dependency=afterany:4626872 scripts/STAR_mapping_fastp_HP.sh # 4626933 ok
```

--> Let's compil the number of uniquely mapped reads for all files (add it in the Google Drive `RNAseq_infos.xlsx` file)

```bash
# Print nb of uniq map reads for raw mapping
for file in output/STAR/fastp/*Log.final.out; do
    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $file | awk '{print $NF}')
    echo "$file: Number of uniquely mapped reads: $uniquely_mapped_reads"
done > output/STAR/fastp/uniq_map_reads_counts.txt
```

--> Overall the fastp is better; more uniquely aligned reads


# Count with featureCounts


Count on gene features with parameter
```bash
conda activate featurecounts

# all samples:
sbatch scripts/featurecounts_CB.sh # 4633070' delete by mistake, relaunch as 463448
sbatch scripts/featurecounts_CT.sh # 4633091 ok
sbatch scripts/featurecounts_HP.sh # 4633092 ok
```

--> More than 80% of succesfully assigned alignments
----> Very few DEGs, with this version; try increase the read mapping rate


```bash
conda activate featurecounts

# all samples:
sbatch scripts/featurecounts_CB_relax1.sh # 4634454
sbatch scripts/featurecounts_CT_relax1.sh # 4634463
sbatch scripts/featurecounts_HP_relax1.sh # 4634470
```
--> More than 90% of succesfully assigned alignments
----> Does not improve the DEGs, same results

Let's try count with the raw reads:

```bash
conda activate featurecounts

# all samples:
sbatch scripts/featurecounts_CB_raw.sh # 4634015
sbatch scripts/featurecounts_CT_raw.sh # run if CB_raw increase nb of DEgs
sbatch scripts/featurecounts_HP_raw.sh # run if CB_raw increase nb of DEgs
```
--> More than 80% of succesfully assigned alignments
----> Does not improve the DEGs, same results

--> Let's use the default version (not the raw, not the relax1); it is the one I am the most confident about


# DEGs with deseq2


**IMPORTANT NOTE: Here it is advisable to REMOVE all genes from chromosome X and Y BEFORE doing the DEGs analysis (X chromosome re-activation occurs in some samples, notably these with more cell passage; in our case, the HET and KO)**


## PCA and clustering with deseq2 hg38

In R; Import all counts and combine into one matrix, deseq2 dataframe (dds)
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("CB14_Het", "CB42_Het", "CB43_Het", "CB20_KO", "CB38_KO" ,"CB41_KO", "CT14_Het", "CT42_Het", "CT43_Het", "CT20_KO" ,"CT38_KO" ,"CT41_KO" ,"HP14_Het" ,"HP42_Het", "HP43_Het" ,"HP20_KO", "HP38_KO", "HP41_KO")

samples <- c("CB14_Het", "CB42_Het", "CB43_Het", "CB20_KO", "CB38_KO" ,"CB41_KO")

samples <- c("HP14_Het" ,"HP42_Het", "HP43_Het" ,"HP20_KO", "HP38_KO", "HP41_KO")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

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
counts_all_matrix = make_matrix(select(counts_all, -Geneid), pull(counts_all, Geneid)) 

## Create colData file that describe all our samples
### Including replicate
coldata_raw <- data.frame(samples) %>%
  mutate(tissue = substr(samples, 1, 2),             # separate the 1st two character
         samples = substr(samples, 3, nchar(samples))) %>% # separate the 1st two character
  separate(samples, into = c("replicate", "genotype"), sep = "_") %>%
  bind_cols(data.frame(samples))

## transform df into matrix
coldata = make_matrix(select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet

#### Tissue and genotype
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype + tissue)
### Genotype only
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype )


# Data normalization
vsd <- vst(dds, blind=TRUE)
rld <- rlog(dds, blind=TRUE) # last 5min

# Vizualization for quality metrics
## Heatmap of the sample-to-sample distances
### vsd 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$tissue, vsd$genotype, vsd$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("output/deseq2/heatmap_cluster_vsd.pdf", width=5, height=6)
pdf("output/deseq2/heatmap_cluster_vsd_HP.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

### rld
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$tissue, rld$genotype, rld$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("output/deseq2/heatmap_cluster_rld.pdf", width=5, height=6)
pdf("output/deseq2/heatmap_cluster_rld_HP.pdf", width=5, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


## PCA
### vsd 
pdf("output/deseq2/PCA_vsd.pdf", width=10, height=10)
pdf("output/deseq2/PCA_vsd_HP.pdf", width=10, height=10)

pcaData <- plotPCA(vsd, intgroup=c("tissue", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=tissue, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()
dev.off()

### rld 
pdf("output/deseq2/PCA_rld.pdf", width=10, height=10)
pdf("output/deseq2/PCA_rld_HP.pdf", width=10, height=10)

pcaData <- plotPCA(rld, intgroup=c("tissue", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=tissue, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()
dev.off()
```



## 'one-by-one' comparison
Comparison tisse per tissue:
- Het vs KO in CB
- Het vs KO in CT
- Het vs KO in HP

*NOTE: For the relax1 featureCounts version, number are with decimal as I used the `fraction` option, so I need to round them to make it work with deseq2 at the dds step.*

###  Het vs KO in CB

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("CB14_Het", "CB42_Het", "CB43_Het", "CB20_KO", "CB38_KO" ,"CB41_KO")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}


for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_relax1/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts_raw/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/raw/")) %>%
    rename(!!sample := starts_with("output/STAR/raw/"))
}


# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
counts_all_matrix = make_matrix(dplyr::select(counts_all, -Geneid), pull(counts_all, Geneid)) 

counts_all_matrix = make_matrix(dplyr::select(counts_all_filtered, -Geneid), pull(counts_all_filtered, Geneid)) 

## Create colData file that describe all our samples
### Including replicate
coldata_raw <- data.frame(samples) %>%
  mutate(tissue = substr(samples, 1, 2),             # separate the 1st two character
         samples = substr(samples, 3, nchar(samples))) %>% # separate the 1st two character
  separate(samples, into = c("replicate", "genotype"), sep = "_") %>%
  bind_cols(data.frame(samples))


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts_all_matrix), # HERE  NEED TO ROUND!!
                              colData = coldata,
                              design= ~ genotype)

dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix, 
                              colData = coldata,
                              design= ~ genotype)


# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "Het")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_Het", type="apeglm")

## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

###### save output
res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2/filtered_CB_KO_vs_CB_Het_lfcnormal.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
######



keyvals <- ifelse(
  res$log2FoldChange < 0 & res$pvalue < 10e-6, 'Sky Blue',
    ifelse(res$log2FoldChange > 0 & res$pvalue < 10e-6, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; FC > 0)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; FC < 0)'

pdf("output/deseq2/plotVolcano_res_CB_KO_vs_CB_Het.pdf", width=7, height=6)    
pdf("output/deseq2/relax1_plotVolcano_res_CB_KO_vs_CB_Het.pdf", width=7, height=6)    
pdf("output/deseq2/raw_plotVolcano_res_CB_KO_vs_CB_Het.pdf", width=7, height=6)    

EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'Het vs KO, CB',
  pCutoff = 10e-6,
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1)  + 
  theme_bw() 
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.1 & res$pvalue < 10e-6)
downregulated_genes <- sum(res$log2FoldChange < -0.1 & res$pvalue < 10e-6)

# Save as gene list for GO analysis:
### Complete table with GeneSymbol
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.1 & res$pvalue < 10e-6, ]
#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.1 & res$pvalue < 10e-6, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_CB_KO_vs_CB_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_CB_KO_vs_CB_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```





###  Het vs KO in CT

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("CT14_Het", "CT42_Het", "CT43_Het", "CT20_KO", "CT38_KO" ,"CT41_KO")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}


# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(samples) %>%
  mutate(tissue = substr(samples, 1, 2),             # separate the 1st two character
         samples = substr(samples, 3, nchar(samples))) %>% # separate the 1st two character
  separate(samples, into = c("replicate", "genotype"), sep = "_") %>%
  bind_cols(data.frame(samples))


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix, 
                              colData = coldata,
                              design= ~ genotype)


# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "Het")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_Het", type="apeglm")



## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

###### save output
res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2/filtered_CT_KO_vs_CT_Het.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
######



keyvals <- ifelse(
  res$log2FoldChange < 0 & res$pvalue < 10e-6, 'Sky Blue',
    ifelse(res$log2FoldChange > 0 & res$pvalue < 10e-6, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; FC > 0)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; FC < 0)'

pdf("output/deseq2/plotVolcano_res_CT_KO_vs_CT_Het.pdf", width=7, height=6)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'Het vs KO, CB',
  pCutoff = 10e-6,
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1)  + 
  theme_bw() 
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.1 & res$pvalue < 10e-6)
downregulated_genes <- sum(res$log2FoldChange < -0.1 & res$pvalue < 10e-6)

# Save as gene list for GO analysis:
### Complete table with GeneSymbol
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.1 & res$pvalue < 10e-6, ]
#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.1 & res$pvalue < 10e-6, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_CT_KO_vs_CT_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_CT_KO_vs_CT_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```





###  Het vs KO in HP

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("HP14_Het", "HP42_Het", "HP43_Het", "HP20_KO", "HP38_KO" ,"HP41_KO")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}


# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(samples) %>%
  mutate(tissue = substr(samples, 1, 2),             # separate the 1st two character
         samples = substr(samples, 3, nchar(samples))) %>% # separate the 1st two character
  separate(samples, into = c("replicate", "genotype"), sep = "_") %>%
  bind_cols(data.frame(samples))


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix, 
                              colData = coldata,
                              design= ~ genotype)


# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "Het")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KO_vs_Het", type="apeglm")



## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$GeneSymbol <- gene_symbols

###### save output
res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2/filtered_HP_KO_vs_HP_Het.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
######



keyvals <- ifelse(
  res$log2FoldChange < 0 & res$pvalue < 10e-6, 'Sky Blue',
    ifelse(res$log2FoldChange > 0 & res$pvalue < 10e-6, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; FC > 0)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; FC < 0)'

pdf("output/deseq2/plotVolcano_res_HP_KO_vs_HP_Het.pdf", width=7, height=6)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'Het vs KO, CB',
  pCutoff = 10e-6,
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1)  + 
  theme_bw() 
dev.off()

upregulated_genes <- sum(res$log2FoldChange > 0.1 & res$pvalue < 10e-6)
downregulated_genes <- sum(res$log2FoldChange < -0.1 & res$pvalue < 10e-6)

# Save as gene list for GO analysis:
### Complete table with GeneSymbol
#### Filter for up-regulated genes
upregulated <- res[res$log2FoldChange > 0.1 & res$pvalue < 10e-6, ]
#### Filter for down-regulated genes
downregulated <- res[res$log2FoldChange < -0.1 & res$pvalue < 10e-6, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_HP_KO_vs_HP_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_HP_KO_vs_HP_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```



## Alternative method to increase the nb of DEGs (without lfc shrinkage)

--> Labeled as *res05*

--> Generate NaiaraPlot for Volcano plot with label genes from the **Excell file (`Genes for volcano plot labels.xlsx`)**:
`c("Mndal", "Lilrb4a", "Gfap", "Ly86", "Cd68", "Ptprc", "Trem2", "Mef2a", "Hexb", "Pros1", "Cst7", "C1qb", "Cd14", "Csf1", "C4b", "Fcer1g", "Lyz2", "Tmem119", "Selplg", "Vsir", "P2ry12")`
AND with **GO:2000427 (positive regulation of apoptotic cell clearance: )**: `c("Abca7", "C2", "C3", "Ccl2", "Cd300lf", "Trem2")`

###  Het vs KO in CB

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("CB14_Het", "CB42_Het", "CB43_Het", "CB20_KO", "CB38_KO" ,"CB41_KO")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(samples) %>%
  mutate(tissue = substr(samples, 1, 2),             # separate the 1st two character
         samples = substr(samples, 3, nchar(samples))) %>% # separate the 1st two character
  separate(samples, into = c("replicate", "genotype"), sep = "_") %>%
  bind_cols(data.frame(samples))


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix, 
                              colData = coldata,
                              design= ~ genotype)


# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "Het")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res05 <- results(dds, alpha=0.05)
summary(res05)


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res05)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res05$GeneSymbol <- gene_symbols

###### save output
res05 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2/filtered_CB_KO_vs_CB_Het_res05.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
res05 = read.table("output/deseq2/filtered_CB_KO_vs_CB_Het_res05.txt", header = TRUE, sep = "\t")
######


keyvals <- ifelse(
  res05$log2FoldChange < 0 & res05$padj < 5e-2, 'Sky Blue',
    ifelse(res05$log2FoldChange > 0 & res05$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; FC > 0)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; FC < 0)'

pdf("output/deseq2/plotVolcano_res05_CB_KO_vs_CB_Het.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Het vs KO, CB',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1)  + 
  theme_bw() 
dev.off()

## NaiaraPlot
### Select genes only if significant
significant_genes <- res05[res05$padj < 5e-2, ]
selected_significant_genes <- intersect(significant_genes$GeneSymbol, c("Mndal", "Lilrb4a", "Gfap", "Ly86", "Cd68", "Ptprc", "Trem2", "Mef2a", "Hexb", "Pros1", "Cst7", "C1qb", "Cd14", "Csf1", "C4b", "Fcer1g", "Lyz2", "Tmem119", "Selplg", "Vsir", "P2ry12", "Abca7", "C2", "C3", "Ccl2", "Cd300lf", "Trem2"))

pdf("output/deseq2/plotVolcano_res05_CB_KO_vs_CB_Het_GenesLabeled.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = selected_significant_genes,
  title = 'Het vs KO, CB',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 5,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  arrowheads = FALSE,
  max.overlaps = 100)  + 
  theme_bw() +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))
dev.off()

pdf("output/deseq2/plotVolcano_res05_CB_KO_vs_CB_Het_GenesLabeled_v2.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = selected_significant_genes,
  title = 'Het vs KO, CB',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 3.5,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  arrowheads = FALSE,
  max.overlaps = 100)  + 
  theme_bw() +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))
dev.off()



upregulated_genes <- sum(res05$log2FoldChange > 0.1 & res05$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res05$log2FoldChange < -0.1 & res05$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange > 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail

### Complete table with GeneSymbol
#### Filter for up-regulated genes
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange > 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail
upregulated <- filtered_res05[filtered_res05$log2FoldChange > 0.1 & filtered_res05$padj < 5e-2, ]
#### Filter for down-regulated genes
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange < 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail
downregulated <- filtered_res05[filtered_res05$log2FoldChange < -0.1 & filtered_res05$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```






###  Het vs KO in CT

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("CT14_Het", "CT42_Het", "CT43_Het", "CT20_KO", "CT38_KO" ,"CT41_KO")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(samples) %>%
  mutate(tissue = substr(samples, 1, 2),             # separate the 1st two character
         samples = substr(samples, 3, nchar(samples))) %>% # separate the 1st two character
  separate(samples, into = c("replicate", "genotype"), sep = "_") %>%
  bind_cols(data.frame(samples))


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix, 
                              colData = coldata,
                              design= ~ genotype)


# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "Het")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res05 <- results(dds, alpha=0.05)
summary(res05)


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res05)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res05$GeneSymbol <- gene_symbols

###### save output
res05 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2/filtered_CT_KO_vs_CT_Het_res05.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
res05 = read.table("output/deseq2/filtered_CT_KO_vs_CT_Het_res05.txt", header = TRUE, sep = "\t")
######



keyvals <- ifelse(
  res05$log2FoldChange < 0 & res05$padj < 5e-2, 'Sky Blue',
    ifelse(res05$log2FoldChange > 0 & res05$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; FC > 0)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; FC < 0)'

pdf("output/deseq2/plotVolcano_res05_CT_KO_vs_CT_Het.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Het vs KO, CT',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1)  + 
  theme_bw() 
dev.off()

## NaiaraPlot
### Select genes only if significant
significant_genes <- res05[res05$padj < 5e-2, ]
selected_significant_genes <- intersect(significant_genes$GeneSymbol, c("Mndal", "Lilrb4a", "Gfap", "Ly86", "Cd68", "Ptprc", "Trem2", "Mef2a", "Hexb", "Pros1", "Cst7", "C1qb", "Cd14", "Csf1", "C4b", "Fcer1g", "Lyz2", "Tmem119", "Selplg", "Vsir", "P2ry12", "Abca7", "C2", "C3", "Ccl2", "Cd300lf", "Trem2"))

pdf("output/deseq2/plotVolcano_res05_CT_KO_vs_CT_Het_GenesLabeled.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = selected_significant_genes,
  title = 'Het vs KO, CT',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 5,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  arrowheads = FALSE,
  max.overlaps = 100)  + 
  theme_bw() +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))
dev.off()
pdf("output/deseq2/plotVolcano_res05_CT_KO_vs_CT_Het_GenesLabeled_v2.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = selected_significant_genes,
  title = 'Het vs KO, CT',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 3.5,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  arrowheads = FALSE,
  max.overlaps = 100)  + 
  theme_bw() +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))
dev.off()

upregulated_genes <- sum(res05$log2FoldChange > 0.1 & res05$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res05$log2FoldChange < -0.1 & res05$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange > 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail

### Complete table with GeneSymbol
#### Filter for up-regulated genes
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange > 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail
upregulated <- filtered_res05[filtered_res05$log2FoldChange > 0.1 & filtered_res05$padj < 5e-2, ]
#### Filter for down-regulated genes
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange < 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail
downregulated <- filtered_res05[filtered_res05$log2FoldChange < -0.1 & filtered_res05$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```







###  Het vs KO in HP

Go in R
```R
# Load packages
library("DESeq2")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("apeglm")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("biomaRt")

# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("HP14_Het", "HP42_Het", "HP43_Het", "HP20_KO", "HP38_KO" ,"HP41_KO")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/fastp/")) %>%
    rename(!!sample := starts_with("output/STAR/fastp/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

# Remove X and Y chromosome genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
### Including replicate
coldata_raw <- data.frame(samples) %>%
  mutate(tissue = substr(samples, 1, 2),             # separate the 1st two character
         samples = substr(samples, 3, nchar(samples))) %>% # separate the 1st two character
  separate(samples, into = c("replicate", "genotype"), sep = "_") %>%
  bind_cols(data.frame(samples))


## transform df into matrix
coldata = make_matrix(dplyr::select(coldata_raw, -samples), pull(coldata_raw, samples))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix, 
                              colData = coldata,
                              design= ~ genotype)


# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "Het")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res05 <- results(dds, alpha=0.05)
summary(res05)


## Plot-volcano
### GeneSymbol ID
gene_ids <- rownames(res05)
stripped_gene_ids <- sub("\\..*", "", gene_ids)
gene_symbols <- mapIds(org.Mm.eg.db, keys = stripped_gene_ids,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res05$GeneSymbol <- gene_symbols

###### save output
res05 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  write.table(file = "output/deseq2/filtered_HP_KO_vs_HP_Het_res05.txt", 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
res05 = read.table("output/deseq2/filtered_HP_KO_vs_HP_Het_res05.txt", header = TRUE, sep = "\t")
######



keyvals <- ifelse(
  res05$log2FoldChange < 0 & res05$padj < 5e-2, 'Sky Blue',
    ifelse(res05$log2FoldChange > 0 & res05$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; FC > 0)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; FC < 0)'

pdf("output/deseq2/plotVolcano_res05_HP_KO_vs_HP_Het.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Het vs KO, CT',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 1)  + 
  theme_bw() 
dev.off()

## NaiaraPlot
### Select genes only if significant
significant_genes <- res05[res05$padj < 5e-2, ]
selected_significant_genes <- intersect(significant_genes$GeneSymbol, c("Mndal", "Lilrb4a", "Gfap", "Ly86", "Cd68", "Ptprc", "Trem2", "Mef2a", "Hexb", "Pros1", "Cst7", "C1qb", "Cd14", "Csf1", "C4b", "Fcer1g", "Lyz2", "Tmem119", "Selplg", "Vsir", "P2ry12", "Abca7", "C2", "C3", "Ccl2", "Cd300lf", "Trem2"))

pdf("output/deseq2/plotVolcano_res05_HP_KO_vs_HP_Het_GenesLabeled.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = selected_significant_genes,
  title = 'Het vs KO, HP',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 5,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  arrowheads = FALSE,
  max.overlaps = 100)  + 
  theme_bw() +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))
dev.off()
pdf("output/deseq2/plotVolcano_res05_HP_KO_vs_HP_Het_GenesLabeled_v2.pdf", width=7, height=6)  
EnhancedVolcano(res05,
  lab = res05$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = selected_significant_genes,
  title = 'Het vs KO, HP',
  pCutoff = 5e-2,         # I can use this pCutoff !!
  FCcutoff = 0,
  pointSize = 1.0,
  labSize = 3.5,
  colCustom = keyvals,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  arrowheads = FALSE,
  max.overlaps = 100)  + 
  theme_bw() +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))
dev.off()


upregulated_genes <- sum(res05$log2FoldChange > 0.1 & res05$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res05$log2FoldChange < -0.1 & res05$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange > 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail

### Complete table with GeneSymbol
#### Filter for up-regulated genes
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange > 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail
upregulated <- filtered_res05[filtered_res05$log2FoldChange > 0.1 & filtered_res05$padj < 5e-2, ]
#### Filter for down-regulated genes
filtered_res05 <- res05[!is.na(res05$log2FoldChange) & !is.na(res05$padj) &
                        res05$log2FoldChange < 0.1 & res05$padj < 5e-2, ]  # REMOVE NA VALUE in eirhter lFC of padj; otherwise will fail
downregulated <- filtered_res05[filtered_res05$log2FoldChange < -0.1 & filtered_res05$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

--> Even though that is not optimal, let's use the non shrinkage method to identify our DEGs and perform functional analysis.
----> Very low DEGs are detected when using the shrinkage method...


# Functional analysis with enrichR

**IMPOPRTANT NOTE: Run the reading and processing ONE BY ONE !!! Otherwise, lead to bug!!!!**

```R
# libr
library("tidyverse")
library("enrichR")

# Define databases for enrichment
dbs <- c("KEGG_2019_Mouse") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt
output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt

output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt
output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt

output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt
output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$KEGG_2019_Mouse
down <- edown$KEGG_2019_Mouse
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)



# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_KEGG_2019_Mouse_CB_KO_vs_CB_Het.pdf", width=8, height=2)
pdf("output/GO/enrichR_KEGG_2019_Mouse_CT_KO_vs_CT_Het.pdf", width=8, height=2)
pdf("output/GO/enrichR_KEGG_2019_Mouse_HP_KO_vs_HP_Het.pdf", width=12, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_CB_KO_vs_CB_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_CT_KO_vs_CT_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KEGG_2019_Mouse_HP_KO_vs_HP_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)








# Define databases for enrichment
dbs <- c("WikiPathways_2019_Mouse") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt
output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt

output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt
output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt

output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt
output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$WikiPathways_2019_Mouse
down <- edown$WikiPathways_2019_Mouse
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_CB_KO_vs_CB_Het.pdf", width=12, height=2)
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_CT_KO_vs_CT_Het.pdf", width=12, height=2.5)
pdf("output/GO/enrichR_WikiPathways_2019_Mouse_HP_KO_vs_HP_Het.pdf", width=12, height=4)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_CB_KO_vs_CB_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_CT_KO_vs_CT_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_WikiPathways_2019_Mouse_HP_KO_vs_HP_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)








# Define databases for enrichment
dbs <- c("KOMP2_Mouse_Phenotypes_2022") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt
output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt

output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt
output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt

output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt
output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$KOMP2_Mouse_Phenotypes_2022
down <- edown$KOMP2_Mouse_Phenotypes_2022
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_KOMP2_Mouse_Phenotypes_2022_CB_KO_vs_CB_Het.pdf", width=12, height=2)
pdf("output/GO/enrichR_KOMP2_Mouse_Phenotypes_2022_CT_KO_vs_CT_Het.pdf", width=12, height=2.5)
pdf("output/GO/enrichR_KOMP2_Mouse_Phenotypes_2022_HP_KO_vs_HP_Het.pdf", width=12, height=4)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_KOMP2_Mouse_Phenotypes_2022_CB_KO_vs_CB_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KOMP2_Mouse_Phenotypes_2022_CT_KO_vs_CT_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_KOMP2_Mouse_Phenotypes_2022_HP_KO_vs_HP_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)









# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt
output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt

output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt
output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt

output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt
output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25) ##  Adjust if you don't want the top 5
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_CB_KO_vs_CB_Het.pdf", width=15, height=8)
pdf("output/GO/enrichR_GO_Biological_Process_2023_CT_KO_vs_CT_Het.pdf", width=15, height=9)
pdf("output/GO/enrichR_GO_Biological_Process_2023_HP_KO_vs_HP_Het.pdf", width=15, height=9)

pdf("output/GO/enrichR_GO_Biological_Process_2023_CT_KO_vs_CT_Het_FigV3.pdf", width=15, height=9)
pdf("output/GO/enrichR_GO_Biological_Process_2023_CB_KO_vs_CB_Het_FigV3.pdf", width=15, height=9)
pdf("output/GO/enrichR_GO_Biological_Process_2023_HP_KO_vs_HP_Het_FigV3.pdf", width=15, height=9)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 25, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 40)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_CB_KO_vs_CB_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_CT_KO_vs_CT_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_HP_KO_vs_HP_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment -- adjusted p value re-order
dbs <- c("GO_Biological_Process_2023") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt
output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt

output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt
output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt

output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt
output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Adjusted.P.value, decreasing = FALSE), ], 5) ##  Adjust if you don't want the top 5
down <- head(down[order(down$Adjusted.P.value, decreasing = FALSE), ], 5) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)



# Plotting with enhanced aesthetics

pdf("output/GO/enrichR_GO_Biological_Process_2023_CT_KO_vs_CT_Het_FigV4.pdf", width=20, height=9)
pdf("output/GO/enrichR_GO_Biological_Process_2023_CB_KO_vs_CB_Het_FigV4.pdf", width=15, height=9)
pdf("output/GO/enrichR_GO_Biological_Process_2023_HP_KO_vs_HP_Het_FigV4.pdf", width=15, height=9)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 10, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 40)
  )
dev.off()




# Define databases for enrichment
dbs <- c("GO_Cellular_Component_2023") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt
output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt

output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt
output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt

output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt
output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Cellular_Component_2023
down <- edown$GO_Cellular_Component_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Cellular_Component_2023_CB_KO_vs_CB_Het.pdf", width=8, height=2)
pdf("output/GO/enrichR_GO_Cellular_Component_2023_CT_KO_vs_CT_Het.pdf", width=8, height=2)
pdf("output/GO/enrichR_GO_Cellular_Component_2023_HP_KO_vs_HP_Het.pdf", width=15, height=6)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Cellular_Component_2023_CB_KO_vs_CB_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Cellular_Component_2023_CT_KO_vs_CT_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Cellular_Component_2023_HP_KO_vs_HP_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("GO_Molecular_Function_2023") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt
output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt

output/deseq2/downregulated_res05_CT_KO_vs_CT_Het.txt
output/deseq2/upregulated_res05_CT_KO_vs_CT_Het.txt

output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt
output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/downregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_res05_HP_KO_vs_HP_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Molecular_Function_2023
down <- edown$GO_Molecular_Function_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Molecular_Function_2023_CB_KO_vs_CB_Het.pdf", width=15, height=2)
pdf("output/GO/enrichR_GO_Molecular_Function_2023_CT_KO_vs_CT_Het.pdf", width=12, height=2)
pdf("output/GO/enrichR_GO_Molecular_Function_2023_HP_KO_vs_HP_Het.pdf", width=12, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Molecular_Function_2023_CB_KO_vs_CB_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Molecular_Function_2023_CT_KO_vs_CT_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Molecular_Function_2023_HP_KO_vs_HP_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
```
--> HDSigDB_Mouse_2021 tested and not clear; and Mouse_Gene_Atlas tested too and give like organs..


--> KEGG; HP many genes related to disease neurodegenerescence!
--> GO BP; show ATP-related stuff for HP! GO CC mitochondria-related






# Bigwig coverage files

Let's generate **TPM coverage**:

```bash
conda activate deeptools
# run time-per-time:
sbatch scripts/TPM_bw_CB.sh # 4635227 ok
sbatch scripts/TPM_bw_CT.sh # 4635235 ok
sbatch scripts/TPM_bw_HP.sh # 4635238 ok
```

Let's merge the bigwig into 1 file with wiggletools (will do average of bigwig signal and not sum, many options see [github](https://github.com/Ensembl/WiggleTools)):

--> Files looks good, replicates looks comparable

**Run wiggletools:**
```bash
conda activate BedToBigwig
sbatch --dependency=afterany:4635227:4635235:4635238 scripts/bigwigmerge_TPM.sh # 4635565v
```
*NOTE: bigwig are merge into 1 bedgraph which is then converted into 1 bigwig (wiggletools cannot output bigwig directly so need to pass by bedgraph or wiggle in between)*



# TPM gene expression
Code to check expression level of some genes (TPM)

## Calculate TPM and RPKM

Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
conda activate deseq2
# Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_TPM.sh # 5376123
# mv all output to output/tpm or rpkm folder
mv output/featurecounts/*tpm* output/tpm/
mv output/featurecounts/*rpkm* output/rpkm_hg38/
```

All good. 

If needed to **display gene with TPM**:

```R
# library
library("biomaRt")
library("tidyverse")
# Plot with TPM instead of baseMean (Naiara plot)
## import tpm
#### Generate TPM for ALL samples
#### collect all samples ID
samples <- c("HP14_Het", "HP42_Het", "HP43_Het", "HP20_KO", "HP38_KO", "HP41_KO",
   "CT14_Het", "CT42_Het", "CT43_Het", "CT20_KO", "CT38_KO", "CT41_KO",
   "CB14_Het", "CB42_Het", "CB43_Het", "CB20_KO", "CB38_KO", "CB41_KO")

## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/tpm/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR.")) %>%
    rename(!!sample := starts_with("output.STAR."))
}

## Merge all dataframe into a single one
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")
# write.csv(tpm_all_sample, file="../001__RNAseq/output/tpm/tpm_all_sample.txt")
### If need to import: tpm_all_sample <- read_csv("../001__RNAseq/output/tpm/tpm_all_sample.txt") %>% dplyr::select(-("...1"))#To import

# plot some genes
tpm_all_sample_tidy <- tpm_all_sample %>%
  gather(key = 'variable', value = 'tpm', -Geneid) %>%
  mutate(tissue = substr(variable, 1, 2),             # separate the 1st two character
         variable = substr(variable, 3, nchar(variable))) %>% # separate the 1st two character
  separate(variable, into = c("replicate", "genotype"), sep = "_") %>%
  rename(gene = Geneid)


tpm_all_sample_tidy$gene <- gsub("\\..*", "", tpm_all_sample_tidy$gene)

## convert gene Ensembl to symbol 
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Convert Ensembl gene IDs to gene symbols
genesymbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = tpm_all_sample_tidy$gene,
                     mart = ensembl)

# Merge gene symbols to your dataframe
tpm_all_sample_tidy <- left_join(tpm_all_sample_tidy, genesymbols, 
                                 by = c("gene" = "ensembl_gene_id"))

# genes 
c("Cd68", "Tlr2", "Trem2") # 


plot_data <- tpm_all_sample_tidy %>%
  unique() %>%
  filter(external_gene_name %in% c("Cd68", "Tlr2", "Trem2")) %>%
  group_by(gene, genotype, tissue,external_gene_name) %>%
  summarise(mean_log2tpm = mean(log2(tpm + 1)),
            se_log2tpm = sd(log2(tpm + 1)) / sqrt(n())) %>%
  ungroup()


plot_data$tissue <-
  factor(plot_data$tissue,
         c("HP", "CT","CB"))

        
# Plot
pdf("output/tpm/tpm__Cd68_Tlr2_Trem2.pdf", width=8, height=4)
ggplot(plot_data, aes(x = external_gene_name, y = mean_log2tpm, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = mean_log2tpm - se_log2tpm, ymax = mean_log2tpm + se_log2tpm),
    width = 0.25,
    position = position_dodge(width = 0.9)
  ) +
  theme_bw() +
  ylab("log2(TPM + 1)") +
  xlab("Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~tissue)
dev.off()


# ADD STAT with ggpubr

library("ggpubr")
my_comparisons <- list( c("Het", "KO") )

c("Cd68", "Tlr2", "Trem2")

tpm_all_sample_tidy$tissue <-
  factor(tpm_all_sample_tidy$tissue,
         c("HP", "CT","CB"))


## Plot
pdf("output/tpm/tpm_Cd68_stat.pdf", width=5, height=4)
pdf("output/tpm/tpm_Tlr2_stat.pdf", width=5, height=4)
pdf("output/tpm/tpm_Trem2_stat.pdf", width=5, height=4)

tpm_all_sample_tidy %>%
  filter(external_gene_name %in% c("Trem2") ) %>% 
  unique() %>%
  mutate(TPM = log2(tpm + 1) ) %>%
    ggboxplot(., x = "genotype", y = "TPM",
                 fill = "genotype",
                 palette = c("grey","red3")) +
      # Add the statistical comparisons
      stat_compare_means(comparisons = my_comparisons, 
                        method = "t.test", 
                        aes(group = genotype)) +
      theme_bw() +
      facet_wrap(~tissue) +
      ylab("log2(TPM + 1)")
dev.off()


# Display expression in heatmap

### genes 

tpm_all_sample_tidy_clean = tpm_all_sample_tidy %>%
  unique()

#### save/import output 
# write.table(tpm_all_sample_tidy_clean, "output/tpm/tpm_all_sample_tidy_clean.txt", sep="\t", row.names=FALSE, quote=FALSE)
# tpm_all_sample_tidy_clean =  as.tibble(read.table("output/tpm/tpm_all_sample_tidy_clean.txt", 
                                        header = TRUE, 
                                        sep = "\t", 
                                        quote = "", 
                                        stringsAsFactors = FALSE) )



c("Impdh1", "Impdh2", "Ampd1", "Ampd2", "Ampd3", "Prps1", "Ppat", "Gart", "Pfas", "Paics", "Adsl", "Atic", "Hprt", "Ada") # 
c("Abca7", "C4a", "Ccl2", "Cd300lf", "Cfb", "Or4c58", "Trem2")
c("Mndal", "Lilrb4a","Gfap","Ly86","Cd68","Ptprc","Trem2","Mef2a","Hexb","Pros1","Cst7","C1qb","Cd14","Csf1","C4b","Fcer1g","Lyz2") # microglia genes V1
c("C3ar1","Cd68","Ctss","Fcer1g","Hexb","Lair1","Ly86","Lyz2","Mef2a","Pabpc1","Sparc","Timp2","Trem2") # microglia genes V2
c("C3ar1","Cd68","Ctss","Fcer1g","Hexb","Lair1","Ly86","Lyz2","Mef2a","Pabpc1","Timp2","Trem2") # microglia genes V3


plot_data <- tpm_all_sample_tidy_clean %>%
  unique() %>%
  filter(external_gene_name %in% c("C3ar1","Cd68","Ctss","Fcer1g","Hexb","Lair1","Ly86","Lyz2","Mef2a","Pabpc1","Timp2","Trem2")) %>%
  group_by(gene, genotype, tissue,external_gene_name) %>%
  summarise(mean_log2tpm = mean(log2(tpm + 1)),
            se_log2tpm = sd(log2(tpm + 1)) / sqrt(n())) %>%
  ungroup() %>%
  dplyr::select(external_gene_name, genotype, tissue, mean_log2tpm) %>%
  unite(sample, genotype, tissue, sep = "_") %>%
  mutate(sample = factor(sample, levels = c("Het_HP", "KO_HP", "Het_CT", "KO_CT", "Het_CB", "KO_CB")))


## Re-order based on Het_HP expr
### Calculate the mean expression for Het_HP sample
Het_HP_expression <- plot_data %>%
  filter(sample == "Het_HP") %>%
  arrange(mean_log2tpm) %>%
  pull(external_gene_name)


# Reorder the gene factor levels based on the Het_HP expression
plot_data$external_gene_name <- factor(plot_data$external_gene_name, levels = Het_HP_expression)

## heatmap

pdf("output/tpm/heatmap_geneList1.pdf", width=5, height=4)
pdf("output/tpm/heatmap_positiveRegulationOfApoptoticCellClearance.pdf", width=5, height=4)
pdf("output/tpm/heatmap_microgliaGenes_V1.pdf", width=5, height=4)
pdf("output/tpm/heatmap_microgliaGenes_V1_formatWide.pdf", width=6, height=3)
pdf("output/tpm/heatmap_microgliaGenes_V1_formatWide2.pdf", width=6, height=4)


pdf("output/tpm/heatmap_microgliaGenes_V2.pdf", width=5, height=4)
pdf("output/tpm/heatmap_microgliaGenes_V2_formatWide.pdf", width=6, height=3)
pdf("output/tpm/heatmap_microgliaGenes_V2_formatWide2.pdf", width=6, height=4)

pdf("output/tpm/heatmap_microgliaGenes_V3.pdf", width=5, height=4)
pdf("output/tpm/heatmap_microgliaGenes_V3_formatWide.pdf", width=6, height=3)
pdf("output/tpm/heatmap_microgliaGenes_V3_formatWide2.pdf", width=6, height=4)

ggplot(plot_data, aes(x = sample, y = external_gene_name, fill = mean_log2tpm)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 4) +    # mid 4 for geneList1; 2.25 for *CellClearance 4.5 for micrlogia genes; V2 5.5
  labs(x = "Sample", y = "Gene", fill = "Expression (log2 TPM)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



## meeting 20240529
### add pvalue / log2fc information (file to use are `deseq2/filtered*res05`)

HP_deseq2 <- as_tibble(read.table("output/deseq2/filtered_HP_KO_vs_HP_Het_res05.txt", header = TRUE, sep = "\t")) %>%
  rename("external_gene_name" = "GeneSymbol") %>%
  inner_join(plot_data) %>%
  dplyr::select(external_gene_name,log2FoldChange,padj) %>%
  unique() %>%
  mutate(sample = list(c("Het_HP", "KO_HP"))) %>%
  unnest(sample)
CT_deseq2 <- as_tibble(read.table("output/deseq2/filtered_CT_KO_vs_CT_Het_res05.txt", header = TRUE, sep = "\t")) %>%
  rename("external_gene_name" = "GeneSymbol") %>%
  inner_join(plot_data) %>%
  dplyr::select(external_gene_name,log2FoldChange,padj) %>%
  unique() %>%
  mutate(sample = list(c("Het_CT", "KO_CT"))) %>%
  unnest(sample)
CB_deseq2 <- as_tibble(read.table("output/deseq2/filtered_CB_KO_vs_CB_Het_res05.txt", header = TRUE, sep = "\t")) %>%
  rename("external_gene_name" = "GeneSymbol") %>%
  inner_join(plot_data) %>%
  dplyr::select(external_gene_name,log2FoldChange,padj) %>%
  unique() %>%
  mutate(sample = list(c("Het_CB", "KO_CB"))) %>%
  unnest(sample)

deseq2 = HP_deseq2 %>%
  bind_rows(CT_deseq2) %>%
  bind_rows(CB_deseq2)

### combine tpm and log2fc
## version1
plot_data_deseq2 = plot_data %>%
  left_join(deseq2) %>%
  mutate(
    log2FoldChange = replace_na(log2FoldChange, 0),
    padj = replace_na(padj, 1),
    neg_log10_padj = -log10(padj) # -log10 padj comonly used; so that high number bigger dot
  )
## version2 - small dot for non signif (Preferred)
plot_data_deseq2 <- plot_data %>%
  left_join(deseq2) %>%
  mutate(
    log2FoldChange = replace_na(log2FoldChange, 0),
    padj = replace_na(padj, 1),
    neg_log10_padj = ifelse(padj > 0.05, 0.1, -log10(padj)) # Set very small size for non-significant padj
  )


## Re-order based on sample
plot_data_deseq2$sample <-
  factor(plot_data_deseq2$sample,
         c("Het_HP", "KO_HP", "Het_CT", "KO_CT", "Het_CB", "KO_CB"))

## Re-order based on Het_HP expr
### Calculate the mean expression for Het_HP sample
Het_HP_expression <- plot_data %>%
  filter(sample == "Het_HP") %>%
  arrange(mean_log2tpm) %>%
  pull(external_gene_name)

plot_data_deseq2$external_gene_name <- factor(plot_data_deseq2$external_gene_name, levels = Het_HP_expression)




# Generate the dot plot
pdf("output/tpm/dotplot_microgliaGenes_V3.pdf", width = 5, height = 4)

ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = neg_log10_padj, color = mean_log2tpm)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 4) +
  labs(x = "Sample", y = "Gene", color = "Expression (log2 TPM)", size = "-log10(padj)") + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


# Custom function to format the size legend
custom_size_scale <- scale_size_continuous(name = "-log10(padj)", 
                                           breaks = c(0.1, 2, 3, 4), 
                                           labels = c("<1.30103", "2", "3", "4"),
                                           range = c(0.5, 4))  # Adjust the minimum size to make it visible

pdf("output/tpm/dotplot_microgliaGenes_V4.pdf", width=5, height=4)
ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = neg_log10_padj, color = mean_log2tpm)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 4) +
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", color = "Expression (log2 TPM)", size = "-log10(padj)")  +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf("output/tpm/dotplot_microgliaGenes_V5.pdf", width=5, height=4)
ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = neg_log10_padj, color = mean_log2tpm)) +
  scale_color_gradientn(colors = c("white", "#FFCCCC", "#FF6666", "#FF3333", "#990000", "#660000"),  # Enhanced red gradient ending with dark red
                        values = scales::rescale(c(0, 3, 4, 4.5, 7, 8)),   # Adjust values to ensure gradient coverage
                        name = "Expression (log2 TPM)") +
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", size = "-log10(padj)") +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  )
dev.off()

pdf("output/tpm/dotplot_microgliaGenes_V6.pdf", width=5, height=4)
ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = neg_log10_padj, color = mean_log2tpm)) +
  scale_color_gradientn(colors = c("white", "#D3D3D3", "#A9A9A9", "#696969", "#4d4d4d", "#333333"),  # White to dark grey gradient
                        values = scales::rescale(c(0, 3, 4, 4.5, 7, 8)),  # Adjust values to ensure gradient coverage
                        name = "Expression (log2 TPM)") +
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", size = "-log10(padj)") +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  )
dev.off()



library("viridis")


pdf("output/tpm/dotplot_microgliaGenes_V7.pdf", width=5, height=4)
ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = neg_log10_padj, color = mean_log2tpm)) +
  scale_color_viridis(name = "Expression (log2 TPM)", option = "inferno") +  # Using magma palette
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", size = "-log10(padj)") +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  )
dev.off()


pdf("output/tpm/dotplot_microgliaGenes_V8.pdf", width=5, height=4)
ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = neg_log10_padj, color = mean_log2tpm)) +
  scale_color_viridis(name = "Expression (log2 TPM)", option = "mako", direction = -1) +  # rocket, inferno, mako
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", size = "-log10(padj)") +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  )
dev.off()


pdf("output/tpm/dotplot_microgliaGenes_V9.pdf", width=5, height=4)
ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = neg_log10_padj, color = mean_log2tpm)) +
  scale_color_gradientn(colors = c("white", "yellow", "#FFCC00", "#FF9900", "#FF6600", "#FF3300", "#CC0000", "#990000"),  # Custom heat gradient
                        values = scales::rescale(c(0, 2, 3, 4, 5, 6, 7, 8)),  # Adjust values to ensure gradient coverage
                        name = "Expression (log2 TPM)") +
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", size = "-log10(padj)") +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  )
dev.off()



## dot size at TPM; color as pvalue, grey color for non significant


pdf("output/tpm/dotplot_microgliaGenes_V10.pdf", width=5, height=4)
ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = mean_log2tpm, color = neg_log10_padj)) +
  scale_color_gradientn(colors = c("#F0F0F0", "#F0F0F0", "#F0F0F0", "#FF3333", "#990000", "#660000"),  # Enhanced red gradient ending with dark red
                        values = scales::rescale(c(0, 1.30103, 2, 3, 4, 6, 8)),   # Adjust values to ensure gradient coverage
                        name = "-log10(padj)") +
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", size = "mean_log2tpm") +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  )
dev.off()



plot_data_deseq2 <- plot_data_deseq2 %>%
  group_by(external_gene_name, sample) %>%
  mutate(median_tpm = median(mean_log2tpm, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(significance = ifelse(padj > 0.05, "Not Significant", "Significant"),
         neg_log10_padj = -log10(padj))

# Custom function to format the size legend
custom_size_scale <- scale_size_continuous(name = "Mean log2 TPM", 
                                           range = c(0.5, 4))  # Adjust the minimum size to make it visible

pdf("output/tpm/dotplot_microgliaGenes_V11.pdf", width=5, height=4)
ggplot(plot_data_deseq2, aes(x = sample, y = external_gene_name)) +
  geom_point(aes(size = mean_log2tpm, color = ifelse(neg_log10_padj < 1.30103, "#F0F0F0", "#FF3333"))) +
  scale_color_identity(name = "-log10(padj)", guide = "legend") +
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", size = "Mean log2 TPM") +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  )
dev.off()



## log2fc and padj only
HP_deseq2 <- as_tibble(read.table("output/deseq2/filtered_HP_KO_vs_HP_Het_res05.txt", header = TRUE, sep = "\t")) %>%
  rename("external_gene_name" = "GeneSymbol") %>%
  inner_join(plot_data) %>%
  dplyr::select(external_gene_name, log2FoldChange, padj) %>%
  unique() %>%
  add_column(tissue = "HP")

CT_deseq2 <- as_tibble(read.table("output/deseq2/filtered_CT_KO_vs_CT_Het_res05.txt", header = TRUE, sep = "\t")) %>%
  rename("external_gene_name" = "GeneSymbol") %>%
  inner_join(plot_data) %>%
  dplyr::select(external_gene_name, log2FoldChange, padj) %>%
  unique() %>%
  add_column(tissue = "CT")

CB_deseq2 <- as_tibble(read.table("output/deseq2/filtered_CB_KO_vs_CB_Het_res05.txt", header = TRUE, sep = "\t")) %>%
  rename("external_gene_name" = "GeneSymbol") %>%
  inner_join(plot_data) %>%
  dplyr::select(external_gene_name, log2FoldChange, padj) %>%
  unique() %>%
  add_column(tissue = "CB")

deseq2 = HP_deseq2 %>%
  bind_rows(CT_deseq2) %>%
  bind_rows(CB_deseq2)

deseq2_clean = deseq2 %>%
  mutate(
    log2FoldChange = replace_na(log2FoldChange, 0),
    padj = replace_na(padj, 1),
    neg_log10_padj = ifelse(padj > 0.05, 0.1, -log10(padj))  # Set very small size for non-significant padj
  )

deseq2_clean$tissue <- factor(deseq2_clean$tissue, levels = c("HP", "CT", "CB"))

# Custom function to format the size legend
custom_size_scale <- scale_size_continuous(name = "-log10(padj)", 
                                           range = c(0.1, 4),  # Adjust the minimum size to make it very small for non-significant values
                                           breaks = c(0.1, 1, 2, 3, 4),
                                           labels = c("<1.30103", "1", "2", "3", "4"))

# Plot
pdf("output/tpm/dotplot_microgliaGenes_V12.pdf", width = 3, height = 3)
ggplot(deseq2_clean, aes(x = tissue, y = external_gene_name)) +
  geom_point(aes(size = neg_log10_padj, color = log2FoldChange)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "log2FoldChange") +
  custom_size_scale + # Apply the custom size scale
  labs(x = "Sample", y = "Gene", size = "-log10(padj)") +
  theme_minimal() + # Change to a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),  # Remove legend background
    legend.key = element_blank()  # Remove legend key background
  ) +
  guides(size = guide_legend(override.aes = list(size = c(0.1, 1, 2, 3, 4))))
dev.off()


## boxplot
### all points!

plot_data <- tpm_all_sample_tidy_clean %>%
  unique() %>%
  filter(external_gene_name %in% c("Mndal", "Lilrb4a","Gfap","Ly86","Cd68","Ptprc","Trem2","Mef2a","Hexb","Pros1","Cst7","C1qb","Cd14","Csf1","C4b","Fcer1g","Lyz2")) %>%
  group_by(gene, genotype, tissue,external_gene_name) %>%
  summarise(tpm = log2(tpm + 1)  ) %>%
  ungroup() %>%
  dplyr::select(external_gene_name, genotype, tissue, tpm) %>%
  unite(sample, genotype, tissue, sep = "_") %>%
  mutate(sample = factor(sample, levels = c("Het_HP", "KO_HP", "Het_CT", "KO_CT", "Het_CB", "KO_CB")))


pdf("output/tpm/boxplot_microgliaGenes_V1_allPoints.pdf", width=5, height=4)
ggplot(plot_data, aes(x = sample, y = tpm)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, color = "black", alpha = 0.5) +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



### mean (BETTER)
plot_data <- tpm_all_sample_tidy_clean %>%
  unique() %>%
  filter(external_gene_name %in% c("Mndal", "Lilrb4a","Gfap","Ly86","Cd68","Ptprc","Trem2","Mef2a","Hexb","Pros1","Cst7","C1qb","Cd14","Csf1","C4b","Fcer1g","Lyz2")) %>%
  group_by(gene, genotype, tissue,external_gene_name) %>%
  summarise(mean_log2tpm = mean(log2(tpm + 1) ) ) %>%
  ungroup() %>%
  dplyr::select(external_gene_name, genotype, tissue, mean_log2tpm) %>%
  unite(sample, genotype, tissue, sep = "_") %>%
  mutate(sample = factor(sample, levels = c("Het_HP", "KO_HP", "Het_CT", "KO_CT", "Het_CB", "KO_CB")))


pdf("output/tpm/boxplot_microgliaGenes_V1.pdf", width=5, height=4)
ggplot(plot_data, aes(x = sample, y = mean_log2tpm)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

### median (LESS GOOD)
plot_data <- tpm_all_sample_tidy_clean %>%
  unique() %>%
  filter(external_gene_name %in% c("Mndal", "Lilrb4a","Gfap","Ly86","Cd68","Ptprc","Trem2","Mef2a","Hexb","Pros1","Cst7","C1qb","Cd14","Csf1","C4b","Fcer1g","Lyz2")) %>%
  group_by(gene, genotype, tissue,external_gene_name) %>%
  summarise(median_log2tpm = median(log2(tpm + 1) ) ) %>%
  ungroup() %>%
  dplyr::select(external_gene_name, genotype, tissue, median_log2tpm) %>%
  unite(sample, genotype, tissue, sep = "_") %>%
  mutate(sample = factor(sample, levels = c("Het_HP", "KO_HP", "Het_CT", "KO_CT", "Het_CB", "KO_CB")))


pdf("output/tpm/boxplot_microgliaGenes_V1_median.pdf", width=6, height=4)
ggplot(plot_data, aes(x = sample, y = median_log2tpm)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_bw()
dev.off()




## boxplot with ggpubr and statistics
## AllPoints (BETTER)
plot_data <- tpm_all_sample_tidy_clean %>%
  unique() %>%
  filter(external_gene_name %in% c("C3ar1","Cd68","Ctss","Fcer1g","Hexb","Lair1","Ly86","Lyz2","Mef2a","Pabpc1","Timp2","Trem2")) %>%
  group_by(gene, genotype, tissue,external_gene_name) %>%
  summarise(tpm = log2(tpm + 1)  ) %>%
  ungroup() %>%
  dplyr::select(external_gene_name, genotype, tissue, tpm) %>%
  unite(sample, genotype, tissue, sep = "_") %>%
  mutate(sample = factor(sample, levels = c("Het_HP", "KO_HP", "Het_CT", "KO_CT", "Het_CB", "KO_CB")))

pdf("output/tpm/boxplot_microgliaGenes_V3_allPoints_stat.pdf", width=6, height=4)
ggboxplot(plot_data, x = "sample", y = "tpm",
  add.params = list(size = 1, alpha = 0.5),
      fill = "sample", palette = c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"),  add = "jitter") + theme_classic() +
  stat_compare_means(comparisons = list( c("Het_HP", "KO_HP"), c("Het_CT", "KO_CT"), c("Het_CB", "KO_CB") ))  # Add pairwise comparisons p-value
dev.off()

## Mean
plot_data <- tpm_all_sample_tidy_clean %>%
  unique() %>%
  filter(external_gene_name %in% c("C3ar1","Cd68","Ctss","Fcer1g","Hexb","Lair1","Ly86","Lyz2","Mef2a","Pabpc1","Timp2","Trem2")) %>%
  group_by(gene, genotype, tissue,external_gene_name) %>%
  summarise(mean_log2tpm = mean(log2(tpm + 1) ) ) %>%
  ungroup() %>%
  dplyr::select(external_gene_name, genotype, tissue, mean_log2tpm) %>%
  unite(sample, genotype, tissue, sep = "_") %>%
  mutate(sample = factor(sample, levels = c("Het_HP", "KO_HP", "Het_CT", "KO_CT", "Het_CB", "KO_CB")))


pdf("output/tpm/boxplot_microgliaGenes_V3_mean_stat.pdf", width=6, height=4)
ggboxplot(plot_data, x = "sample", y = "mean_log2tpm",
  add.params = list(size = 1, alpha = 0.5),
      fill = "sample", palette = c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"),  add = "jitter") + theme_classic() +
  stat_compare_means(comparisons = list( c("Het_HP", "KO_HP"), c("Het_CT", "KO_CT"), c("Het_CB", "KO_CB") ))  # Add pairwise comparisons p-value
dev.off()


# adjust pvalue

p.adjust(c(0.037, 0.12, 0.076), method = "bonferroni")

p.adjust(c(0.037, 0.12, 0.076), method = "BH")


```








# GSEA


Let's do GSEA analysis with the gene set lists from Naiara; include:
- B_cells
- DAM_microglia
- Granulocytes_1
- Granulocytes_2
- Immature_B_cells
- Microglia_1
- Microglia_2
- Microglia_3
- Monocyte
- Perivascular_MF
- T_NK_cells



```R
# Packages
library("tidyverse")
library("clusterProfiler")
library("msigdbr") # BiocManager::install("msigdbr")
library("org.Mm.eg.db")
library("enrichplot") # for gseaplot2()
library("pheatmap")

# 


# import DEGs
HP = read.table("output/deseq2/filtered_HP_KO_vs_HP_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
HP_geneSymbol = HP %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene
CB = read.table("output/deseq2/filtered_CB_KO_vs_CB_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
CB_geneSymbol = CB %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene
CT = read.table("output/deseq2/filtered_CT_KO_vs_CT_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
CT_geneSymbol = CT %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene
# import gene signature marker lists
## example for 1
B_cells = read.table("output/gsea/B_cells.txt", header = FALSE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::rename("gene" = "V1") %>%
  add_column(cellName = "B_cells")
## as a function
read_cell_file <- function(cell_name) {
  filepath <- paste0("output/gsea/", cell_name, ".txt")
  df <- read.table(filepath, header = FALSE, sep = "\t") %>%
    as_tibble() %>%
    dplyr::rename("gene" = "V1") %>%
    add_column(cellName = cell_name)
  
  return(df)
} 
cell_names <- c("B_cells", "DAM_microglia", "Granulocytes_1", "Granulocytes_2", "Immature_B_cells", "Microglia_1", "Microglia_2", "Microglia_3", "Monocyte", "Perivascular_MF", "T_NK_cells")
all_data <- lapply(cell_names, read_cell_file) %>%
  bind_rows()


# Order our DEG
## Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- CT_geneSymbol$log2FoldChange  ### CHAGNE HERE DATA!!!!!!!
names(lfc_vector) <- CT_geneSymbol$GeneSymbol ### CHAGNE HERE DATA!!!!!!!
## We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
### Set the seed so our results are reproducible:
set.seed(42)


# run GSEA
## without pvalue
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = all_data %>% dplyr::select(cellName,gene), # Need to be in that order...
)



gsea_result_df <- data.frame(gsea_results@result)
# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_CT_complete.tsv"
  )
)

# plots
c("B_cells", "DAM_microglia", "Granulocytes_1", "Granulocytes_2", "Immature_B_cells", "Microglia_1", "Microglia_2", "Microglia_3", "Monocyte", "Perivascular_MF", "T_NK_cells")


pdf("output/gsea/CT_T_NK_cells.pdf", width=10, height=8)

enrichplot::gseaplot(
  gsea_results,
  geneSetID = "T_NK_cells",
  title = "T_NK_cells",
  color.line = "#0d76ff"
)
dev.off()



### Set up the heatmap
c("B_cells", "Immature_B_cells","T_NK_cells") # Immune Cells Lymphoid 
c("Microglia_1", "Microglia_2", "Microglia_3", "DAM_microglia") # Microglia
c("Granulocytes_1", "Granulocytes_2" , "Monocyte", "Perivascular_MF") # Immune Cells Myeloid 

gsea_result_df_CT <- gsea_result_df
gsea_result_df_CB <- gsea_result_df
gsea_result_df_HP <- gsea_result_df


desired_ids <- c(
"B_cells", "Immature_B_cells","T_NK_cells","Microglia_1", "Microglia_2", "Microglia_3", "DAM_microglia","Granulocytes_1", "Granulocytes_2" , "Monocyte", "Perivascular_MF"
)

gsea_result_df_CT_filt = as_tibble(gsea_result_df_CT) %>%
  dplyr::select(ID, enrichmentScore, qvalue) %>% # change btween pvalue, qvalue,p.adjust
  add_column(tissue = "CT")
gsea_result_df_CB_filt = as_tibble(gsea_result_df_CB) %>%
  dplyr::select(ID, enrichmentScore, qvalue) %>%
  add_column(tissue = "CB")
gsea_result_df_HP_filt = as_tibble(gsea_result_df_HP) %>%
  dplyr::select(ID, enrichmentScore, qvalue) %>%
  add_column(tissue = "HP")

gsea_result_df_tidy = gsea_result_df_CT_filt %>%
  bind_rows(gsea_result_df_CB_filt) %>%
  bind_rows(gsea_result_df_HP_filt)
# Filter the data for desired IDs

filtered_data <- gsea_result_df_tidy %>%
  filter(ID %in% desired_ids)

filtered_data$ID = factor(filtered_data$ID, c("B_cells", "Immature_B_cells","T_NK_cells","Microglia_1", "Microglia_2", "Microglia_3", "DAM_microglia","Granulocytes_1", "Granulocytes_2" , "Monocyte", "Perivascular_MF"
))

pdf("output/gsea/heatmap_gsea_qvalue.pdf", width=3, height=3)
ggplot(filtered_data, aes(x=ID, y=tissue, fill=enrichmentScore)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Enrichment\nScore") +
  geom_text(aes(label=sprintf("%.2f", enrichmentScore)), 
            color = ifelse(filtered_data$qvalue <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()


```

Let's GSEA with updated Naiara gene list (20231120) and test:
- log2fc ranking
- log2fc and pvalue ranking (combined_score = log2FoldChange * -log10(pvalue) )




```R
# Packages
library("tidyverse")
library("clusterProfiler")
library("msigdbr") # BiocManager::install("msigdbr")
library("org.Mm.eg.db")
library("enrichplot") # for gseaplot2()
library("pheatmap")

# import DEGs
HP = read.table("output/deseq2/filtered_HP_KO_vs_HP_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
HP_geneSymbol = HP %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene
CB = read.table("output/deseq2/filtered_CB_KO_vs_CB_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
CB_geneSymbol = CB %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene
CT = read.table("output/deseq2/filtered_CT_KO_vs_CT_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
CT_geneSymbol = CT %>%
  filter(!is.na(GeneSymbol)) # filter to keep only the geneSymbol gene


# import gene signature lists
## top 100
read_cell_file <- function(cell_name) {
  filepath <- paste0("output/gsea/", cell_name, "_v2.txt")
  df <- read.table(filepath, header = FALSE, sep = "\t") %>%
    as_tibble() %>%
    dplyr::rename("gene" = "V1") %>%
    add_column(cellName = cell_name)
  
  return(df)
} 
cell_names <- c("DAM_microglia", "homeostatic_microglia")
all_data <- lapply(cell_names, read_cell_file) %>%
  bind_rows()

## top 50
read_cell_file <- function(cell_name) {
  filepath <- paste0("output/gsea/", cell_name, "_top50_v2.txt")
  df <- read.table(filepath, header = FALSE, sep = "\t") %>%
    as_tibble() %>%
    dplyr::rename("gene" = "V1") %>%
    add_column(cellName = cell_name)
  
  return(df)
} 
cell_names <- c("DAM_microglia", "homeostatic_microglia")
all_data <- lapply(cell_names, read_cell_file) %>%
  bind_rows()


# METHOD 1 log2fc ordering only

## Order our DEG
### Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- CT_geneSymbol$log2FoldChange  ### CHAGNE HERE DATA!!!!!!!
names(lfc_vector) <- CT_geneSymbol$GeneSymbol ### CHAGNE HERE DATA!!!!!!!
### We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
### Set the seed so our results are reproducible:
set.seed(42)


# run GSEA
## without pvalue
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = all_data %>% dplyr::select(cellName,gene), # Need to be in that order...
)



gsea_result_df <- data.frame(gsea_results@result)
# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_CT_top50_v2_complete.tsv"
  )
)

# plots
c("DAM_microglia", "homeostatic_microglia")


pdf("output/gsea/CT_DAM_microglia_top50_v2.pdf", width=10, height=8)

enrichplot::gseaplot(
  gsea_results,
  geneSetID = "DAM_microglia",
  title = "DAM_microglia",
  color.line = "#0d76ff"
)
dev.off()



### Set up the heatmap
c("DAM_microglia", "homeostatic_microglia")

gsea_result_df_CT <- gsea_result_df
gsea_result_df_CB <- gsea_result_df
gsea_result_df_HP <- gsea_result_df


desired_ids <- c("DAM_microglia", "homeostatic_microglia")

gsea_result_df_CT_filt = as_tibble(gsea_result_df_CT) %>%
  dplyr::select(ID, NES, qvalue) %>% # change btween pvalue, qvalue,p.adjust
  add_column(tissue = "CT")
gsea_result_df_CB_filt = as_tibble(gsea_result_df_CB) %>%
  dplyr::select(ID, NES, qvalue) %>%
  add_column(tissue = "CB")
gsea_result_df_HP_filt = as_tibble(gsea_result_df_HP) %>%
  dplyr::select(ID, NES, qvalue) %>%
  add_column(tissue = "HP")

gsea_result_df_tidy = gsea_result_df_CT_filt %>%
  bind_rows(gsea_result_df_CB_filt) %>%
  bind_rows(gsea_result_df_HP_filt)
# Filter the data for desired IDs

filtered_data <- gsea_result_df_tidy %>%
  filter(ID %in% desired_ids)

filtered_data$ID = factor(filtered_data$ID, c("homeostatic_microglia", "DAM_microglia"))

pdf("output/gsea/heatmap_gsea_qvalue_v2.pdf", width=3.5, height=3.5)
pdf("output/gsea/heatmap_gsea_qvalue_top50_v2.pdf", width=3.5, height=3.5)

ggplot(filtered_data, aes(x=ID, y=tissue, fill=NES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="NES") +
  geom_text(aes(label=sprintf("%.2f", NES)), 
            color = ifelse(filtered_data$qvalue <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=4) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()





# METHOD2 order log2fc and pvalue

HP = read.table("output/deseq2/filtered_HP_KO_vs_HP_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
HP_geneSymbol = HP %>%
  filter(!is.na(GeneSymbol)) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue))
CB = read.table("output/deseq2/filtered_CB_KO_vs_CB_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
CB_geneSymbol = CB %>%
  filter(!is.na(GeneSymbol)) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue))
CT = read.table("output/deseq2/filtered_CT_KO_vs_CT_Het_res05.txt", header = TRUE, sep = "\t") %>%
  as_tibble() 
CT_geneSymbol = CT %>%
  filter(!is.na(GeneSymbol)) %>%
  mutate(combined_score = log2FoldChange * -log10(pvalue))





## Order our DEG
### Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- HP_geneSymbol$combined_score  ### CHAGNE HERE DATA!!!!!!!
names(lfc_vector) <- HP_geneSymbol$GeneSymbol ### CHAGNE HERE DATA!!!!!!!
### We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)
### Set the seed so our results are reproducible:
set.seed(42)


# run GSEA
## without pvalue
gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 1,
  maxGSSize = 5000,
  pvalueCutoff = 1,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = all_data %>% dplyr::select(cellName,gene), # Need to be in that order...
)



gsea_result_df <- data.frame(gsea_results@result)
# Save output
readr::write_tsv(
  gsea_result_df,
  file.path("output/gsea/gsea_results_HP_top50_v2_complete_combinedScore.tsv"
  )
)

# plots
c("DAM_microglia", "homeostatic_microglia")


pdf("output/gsea/HP_homeostatic_microglia_top50_v2_combinedScore.pdf", width=10, height=8)

enrichplot::gseaplot(
  gsea_results,
  geneSetID = "homeostatic_microglia",
  title = "homeostatic_microglia",
  color.line = "#0d76ff"
)
dev.off()



### Set up the heatmap
c("DAM_microglia", "homeostatic_microglia")

gsea_result_df_CT <- gsea_result_df
gsea_result_df_CB <- gsea_result_df
gsea_result_df_HP <- gsea_result_df


desired_ids <- c("DAM_microglia", "homeostatic_microglia")

gsea_result_df_CT_filt = as_tibble(gsea_result_df_CT) %>%
  dplyr::select(ID, NES, p.adjust) %>% # change btween pvalue, qvalue,p.adjust
  add_column(tissue = "CT")
gsea_result_df_CB_filt = as_tibble(gsea_result_df_CB) %>%
  dplyr::select(ID, NES, p.adjust) %>%
  add_column(tissue = "CB")
gsea_result_df_HP_filt = as_tibble(gsea_result_df_HP) %>%
  dplyr::select(ID, NES, p.adjust) %>%
  add_column(tissue = "HP")

gsea_result_df_tidy = gsea_result_df_CT_filt %>%
  bind_rows(gsea_result_df_CB_filt) %>%
  bind_rows(gsea_result_df_HP_filt)
# Filter the data for desired IDs

filtered_data <- gsea_result_df_tidy %>%
  filter(ID %in% desired_ids)

filtered_data$ID = factor(filtered_data$ID, c("homeostatic_microglia", "DAM_microglia"))

pdf("output/gsea/heatmap_gsea_padjust_v2_combinedScore.pdf", width=3.5, height=3.5)
pdf("output/gsea/heatmap_gsea_padjust_top50_v2_combinedScore.pdf", width=3.5, height=3.5)

ggplot(filtered_data, aes(x=ID, y=tissue, fill=NES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="NES") +
  geom_text(aes(label=sprintf("%.2f", NES)), 
            color = ifelse(filtered_data$p.adjust <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=4) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()



```










# Brain deconvolution

Let's try the [BrainDeconvShinny](https://voineagulab.shinyapps.io/BrainDeconvShiny/)

```R
#### collect all samples ID
samples <- c("HP14_Het", "HP42_Het", "HP43_Het", "HP20_KO", "HP38_KO", "HP41_KO",
   "CT14_Het", "CT42_Het", "CT43_Het", "CT20_KO", "CT38_KO", "CT41_KO",
   "CB14_Het", "CB42_Het", "CB43_Het", "CB20_KO", "CB38_KO", "CB41_KO")



## Make a loop for importing all rpkm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/rpkm/", sample, "_rpkm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR.")) %>%
    rename(!!sample := starts_with("output.STAR."))
}

## Merge all dataframe into a single one
rpkm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")

rpkm_all_sample_tidy <- rpkm_all_sample %>%
  gather(key = 'variable', value = 'rpkm', -Geneid) %>%
  mutate(tissue = substr(variable, 1, 2),             # separate the 1st two character
         variable = substr(variable, 3, nchar(variable))) %>% # separate the 1st two character
  separate(variable, into = c("replicate", "genotype"), sep = "_") %>%
  rename(gene = Geneid)

rpkm_all_sample_tidy$gene <- gsub("\\..*", "", rpkm_all_sample_tidy$gene) # remove Ensembl gene id version



## Save sample per sample
rpkm_all_sample_tidy_HP_Het_BrainDeconvShiny = rpkm_all_sample_tidy %>%
  filter(tissue == "HP", genotype == "KO") %>%
  dplyr::select(gene,replicate,rpkm) %>%
  group_by(gene, replicate) %>%
  summarise(rpkm = mean(rpkm, na.rm = TRUE)) %>%  # some gene id are dupplicated as we remove version id; need do mean...
  spread(key = replicate, value = rpkm)

write.table(rpkm_all_sample_tidy_HP_Het_BrainDeconvShiny, file = "output/rpkm/rpkm_all_sample_tidy_HP_KO_BrainDeconvShiny.txt", sep = "\t", quote = FALSE, row.names = FALSE)

```



# Venn diagram detected gene in each tissue

Let's identify specifically/uniquely expressed genes in each tissue; in Het 1st:
- Import tpm table
- add a detected column (> 0.5 or 1 = yes)
- Export the detected gene list for each tissue
- Venn diagram [online](https://bioinformatics.psb.ugent.be/webtools/Venn/)





```R
#### collect all samples ID
samples <- c("HP14_Het", "HP42_Het", "HP43_Het", "HP20_KO", "HP38_KO", "HP41_KO",
   "CT14_Het", "CT42_Het", "CT43_Het", "CT20_KO", "CT38_KO", "CT41_KO",
   "CB14_Het", "CB42_Het", "CB43_Het", "CB20_KO", "CB38_KO", "CB41_KO")



## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/tpm/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR.")) %>%
    rename(!!sample := starts_with("output.STAR."))
}

## Merge all dataframe into a single one
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")

tpm_all_sample_tidy <- tpm_all_sample %>%
  gather(key = 'variable', value = 'tpm', -Geneid) %>%
  mutate(tissue = substr(variable, 1, 2),             # separate the 1st two character
         variable = substr(variable, 3, nchar(variable))) %>% # separate the 1st two character
  separate(variable, into = c("replicate", "genotype"), sep = "_") %>%
  rename(gene = Geneid)

tpm_all_sample_tidy$gene <- gsub("\\..*", "", tpm_all_sample_tidy$gene) # remove Ensembl gene id version

## Calculate median for each sample
tpm_all_sample_tidy_median = tpm_all_sample_tidy %>%
  group_by(gene, genotype, tissue) %>%
  summarise(median = median(tpm))  

tpm_all_sample_tidy_median_detected = tpm_all_sample_tidy_median %>%
  mutate(detected = ifelse(median > 0.5, "yes", "no"))


## Save detected gene list

write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "HP", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected05_HP_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "CB", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected05_CB_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "CT", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected05_CT_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "HP", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected1_HP_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "CB", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected1_CB_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "CT", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected1_CT_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "HP", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected5_HP_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "CB", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected5_CB_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "Het", tissue == "CT", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected5_CT_Het.txt", sep = "\t", quote = FALSE, row.names = FALSE)


write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "HP", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected05_HP_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "CB", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected05_CB_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "CT", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected05_CT_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "HP", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected1_HP_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "CB", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected1_CB_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "CT", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected1_CT_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "HP", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected5_HP_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "CB", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected5_CB_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tpm_all_sample_tidy_median_detected %>% filter(genotype == "KO", tissue == "CT", detected == "yes") %>% ungroup() %>% dplyr::select(gene) %>% unique(), file = "output/tpm/detected5_CT_KO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

# Functional analysis dotplots

I did venn diagram to identify the HP specific genes and HP-gentoype-specific one. Let's use the 0.5 tpm treshold detection:
- HP specific genes in Het; 934, `output/GO/HP_HET_934.txt`
- HP specific genes in KO; 993,  `output/GO/HP_KO_993.txt`
- Het-specific HP-specific genes; 808, `output/GO/HP_HETspe_808`
- KO-specofoc HP-specific genes; 867, `output/GO/HP_KOspe_867`



```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Mm.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")


## Read GTF file
gtf_file <- "../../Master/meta/ENCFF159KBI.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()


# Import genes_cluster list 
HP_HET_934 = read_table(file = "output/GO/HP_HET_934.txt", col_names = FALSE)
HP_KO_993 = read_table(file = "output/GO/HP_KO_993.txt", col_names = FALSE)
HP_HETspe_808 = read_table(file = "output/GO/HP_HETspe_808.txt", col_names = FALSE)
HP_KOspe_867 = read_table(file = "output/GO/HP_KOspe_867.txt", col_names = FALSE)



## Run GO enrichment analysis 

ego <- enrichGO(gene = as.character(HP_KOspe_867$X1), 
                keyType = "ENSEMBL",     # Use ENSEMBL if want to use ENSG000XX format
                OrgDb = org.Mm.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
###

## Save GO analyses
GO_summary <- data.frame(ego)

write.csv(GO_summary, "output/GO/BP_HP_KOspe_867.csv")
write.csv(GO_summary, "output/GO/MF_HP_KOspe_867.csv")
write.csv(GO_summary, "output/GO/CC_HP_KOspe_867.csv")

## Vizualization

pdf("output/GO/dotplot_MF_HP_KOspe_867.pdf", width=6, height=3)
dotplot(ego, showCategory=20)
dev.off()
pdf("output/GO/emapplot_MF_HP_KOspe_867.pdf", width=12, height=14)
emapplot(pairwise_termsim(ego), showCategory = 20)
dev.off()
```





# Venn diagram with DEGs

## With all DEGs together (up and down)

- Collect all DEGs tissue per tissue (Het vs KO): *on xls*
- Venn diagram tissue per tissue: *on Venn webtool*
- GO on tissue specific DEGs genes: *Files transfered to `output/deseq2/Venn_DEG_*spe_*.txt`*


--> Below code modified to show only 1 set of genes (no up and down, only 1 set)


```R
# packages

library("tidyverse")
library("enrichR")


# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 

### GeneSymbol list of DEGs per tissue
output/deseq2/Venn_DEG_HPspe_426.txt
output/deseq2/Venn_DEG_CTspe_118.txt
output/deseq2/Venn_DEG_CBspe_268.txt

# IF starting with geneSymbol

## Read and preprocess data for DEGs genes
gene_names_up <- read.csv("output/deseq2/Venn_DEG_CBspe_268.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
up$type <- "up"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 50) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)

# Combine the two dataframes
gos <- up
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
new_order <- up_pathways
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
## gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_HPspe_426.pdf", width=20, height=8)
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_CTspe_118.pdf", width=10, height=2)
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_CBspe_268.pdf", width=20, height=6)



ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 12, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Green")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_HPspe_426.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_CTspe_118.txt", sep="\t", row.names=FALSE, quote=FALSE)






# Define databases for enrichment
dbs <- c("MSigDB_Hallmark_2020") # 

### GeneSymbol list of DEGs per tissue
output/deseq2/Venn_DEG_HPspe_426.txt
output/deseq2/Venn_DEG_CTspe_118.txt
output/deseq2/Venn_DEG_CBspe_268.txt

# IF starting with geneSymbol

## Read and preprocess data for DEGs genes
gene_names_up <- read.csv("output/deseq2/Venn_DEG_CBspe_268.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$MSigDB_Hallmark_2020
up$type <- "up"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 50) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)

# Combine the two dataframes
gos <- up
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
new_order <- up_pathways
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
## gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_CBspe_268.pdf", width=20, height=6)


pdf("output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_HPspe_426.pdf", width=15, height=8)
pdf("output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_CTspe_118.pdf", width=10, height=3.5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 12, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Green")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_HPspe_426.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_CTspe_118.txt", sep="\t", row.names=FALSE, quote=FALSE)





```



## With up and down DEGs separated

- Collect **up** DEGs tissue per tissue (Het vs KO)
- Venn diagram tissue per tissue
- GO on tissue specific DEGs genes
--> Repeat with **down** DEGs


```R
# packages
library("tidyverse")
library("enrichR")



# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/Venn_DEG_Up_HPspe_216.txt
output/deseq2/Venn_DEG_Down_HPspe_211.txt

output/deseq2/Venn_DEG_Up_CTspe_44.txt
output/deseq2/Venn_DEG_Down_CTspe_74.txt

output/deseq2/Venn_DEG_Up_CBspe_157.txt
output/deseq2/Venn_DEG_Down_CBspe_113.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/Venn_DEG_Down_CBspe_113.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/Venn_DEG_Up_CBspe_157.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25) ##  Adjust if you don't want the top 5
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
## gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_Down_Up_HPspe.pdf", width=20, height=8)
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_Down_Up_CBspe.pdf", width=20, height=8)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 12, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_Down_Up_HPspe.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_Down_Up_CBspe.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("MSigDB_Hallmark_2020") # 


### GeneSymbol list of signif up/down genes in each genotypes
output/deseq2/Venn_DEG_Up_HPspe_216.txt
output/deseq2/Venn_DEG_Down_HPspe_211.txt

output/deseq2/Venn_DEG_Up_CTspe_44.txt
output/deseq2/Venn_DEG_Down_CTspe_74.txt

output/deseq2/Venn_DEG_Up_CBspe_157.txt
output/deseq2/Venn_DEG_Down_CBspe_113.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/deseq2/Venn_DEG_Down_HPspe_211.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/Venn_DEG_Up_HPspe_216.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$MSigDB_Hallmark_2020
down <- edown$MSigDB_Hallmark_2020
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25) ##  Adjust if you don't want the top 5
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
## gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics

pdf("output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_Down_Up_CBspe.pdf", width=15, height=2)
pdf("output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_Down_Up_CTspe.pdf", width=15, height=4)
pdf("output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_Down_Up_HPspe.pdf", width=15, height=6)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 12, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG pathways") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_Venn_DEG_Down_Up_HPspe.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_Down_Up_CBspe.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_MSigDB_Hallmark_2020_Venn_DEG_Down_Up_CTspe.txt", sep="\t", row.names=FALSE, quote=FALSE)




```




# Upload files to GEO

Go [here](https://www.ncbi.nlm.nih.gov/geo/info/seq.html); and follow instructions in `Transfer Files`. Connect to my personal space (`uploads/thomasroule@orcid_A787EGG4`) and transfer files.

- Create a clean `GEO` folder with all `*fq.gz` and `*bigwig` (re-name file so that they have same prefix; only extension differ)
- Fill in the `seq_template.xlsx` (`Metada` and `MD5` sheet notably)
- submit files

```bash
# do file integrity check with md5
md5sum * | awk '{print $2 "\t" $1}' > md5sums.txt




module load lftp

# connect to ftp
lftp -u geoftp,inAlwokhodAbnib5 ftp-private.ncbi.nlm.nih.gov # geoftp = username; inAlwokhodAbnib5 = pwd
cd uploads/thomasroule@orcid_A787EGG4

mirror -R ../001__RNAseq/geo_sub_RNAseq_AMPD2/
```

--> **FOLDER `geo_sub_RNAseq_AMPD2` DELETED TO FREE SPACE.** 
