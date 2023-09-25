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
gene_names_down <- read.csv("output/deseq2/downregulated_res05_CB_KO_vs_CB_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/deseq2/upregulated_res05_CB_KO_vs_CB_Het.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
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
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)


# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_CB_KO_vs_CB_Het.pdf", width=15, height=8)
pdf("output/GO/enrichR_GO_Biological_Process_2023_CT_KO_vs_CT_Het.pdf", width=16, height=6)
pdf("output/GO/enrichR_GO_Biological_Process_2023_HP_KO_vs_HP_Het.pdf", width=15, height=6)

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
    axis.text.x = element_text(size = 20)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_CB_KO_vs_CB_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_CT_KO_vs_CT_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_HP_KO_vs_HP_Het.txt", sep="\t", row.names=FALSE, quote=FALSE)







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
write.csv(tpm_all_sample, file="../001__RNAseq/output/tpm/tpm_all_sample.txt")
### If need to import: tpm_all_sample <- read_csv("../001__RNAseq/output/tpm/tpm_all_sample.txt") #To import

XXXX

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
XXXX


library("ggpubr")
my_comparisons <- list( c("Het", "KO") )

# Plot
pdf("output/tpm/tpm__Cd68_Tlr2_Trem2_stat.pdf", width=8, height=4)
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
  facet_wrap(~tissue) +
  # Add the statistical comparisons
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label.y = 5,       # Adjust this to position the p-value labels correctly
                     aes(group = genotype))
dev.off()









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

















