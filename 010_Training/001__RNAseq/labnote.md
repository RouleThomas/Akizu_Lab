# Overview

RNAseq data analysis - two conditions (ie. WT vs mutants; untreated vs treated...):
- Data cleaning (FASTP, FASTQC)
- Mapping (STAR)
- Counting (featureCounts)
- DGE (DESeq2)



Example will be on three bio rep of 2 conditions, **WT and KO in human**;
    --> sample name: `WT_Rep1, WT_Rep2, ..., KO_Rep1, ....`

*Paired end Directional RNAseq* --> If data is something else (undirectional, single end; please read docs to adapt each code)


# Data cleaning


Install [fastp](https://github.com/OpenGene/fastp).

```bash
# Download in Master/Software/
wget https://opengene.org/fastp/fastp # download fastp
chmod a+x ./fastp # install fastp
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software



# Restart terminal to apply changes!!


# Test installation:
fastp --help
```

Run fastp:

```bash
# create fastp output firectory
mkdir output/fastp

# Create script
nano scripts/fastp.sh

# Run cleaning
sbatch scripts/fastp.sh #   

# OPTIONAL - Run FASTQC on raw and clean fastq
## create output folders
mkdir output/fastqc/raw
mkdir output/fastqc/fastp

sbatch scripts/fastqc_raw.sh #   
sbatch scripts/fastqc_fastp.sh #   
```

--> Script will output clean reads to `output/fastp`


# Data mapping

## Prerequisite - installation

Download human genome:
- [ENCODE](https://www.encodeproject.org/data-standards/reference-sequences/) --> Follow guideline to download

--> Transfer genome FASTA file to `/master/meta` : For this example:
- genome used is: `/master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta`
- annotation used is: `/master/meta/gencode.v47.annotation.gtf`

```bash
# Create output folder
mkdir /Master/meta/STAR_hg38
# Create scripts
nano scripts/STAR_index_hg38.sh
# Run script
## Index human genome for use with STAR
sbatch scripts/STAR_index_hg38.sh #  
```


## Run mapping

```bash
# Create output folder
mkdir output/STAR
# Create scripts
nano scripts/STAR_mapping.sh
# Run script
sbatch scripts/STAR_mapping.sh
```

--> Script will output mapped reads to `output/STAR`


# Count with featureCounts

Let's count the number of reads within gene, using our `/master/meta/gencode.v47.annotation.gtf` annotation

## Prerequisite - installation

Create a **featurecounts; conda environment**

```bash
conda create -c bioconda -n featurecounts subread
conda activate featurecounts
```

## Run featureCounts

```bash
# activate conda environment
conda activate featurecounts

# Create output folder
mkdir output/featurecounts

# run counting
sbatch scripts/featurecounts.sh # 
```

--> Script will output mapped reads to `output/featurecounts`


# Generate Bigwig coverage files

## Prerequisite - installation

Create a **deeptools; conda environment**
```bash
conda create -n deeptools -c bioconda deeptools # 
```

## Run bamCoverage

Let's generate **bigwig file for vizualization in IGV**:

```bash
# Activate conda environment
conda activate deeptools

# Create output folder
mkdir output/bigwig

# run 
sbatch scripts/TPM_bw.sh # 
```

--> Script will output bigwig coverage files to `output/bigwig`


To generate **median tracks**; first create a **BedToBigwig conda environment**
```bash
conda create -n BedToBigwig
conda install -c bioconda bedtools
conda install -c bioconda ucsc-bedgraphtobigwig
conda install -n ucsc openssl=1.0 
```

Then generate bigwig files:
```bash
# Activate conda environment
conda activate BedToBigwig

# run 
sbatch scripts/bigwigmerge_STAR_TPM_bw.sh # 
```



# Calculate TPM and RPKM

Use custom R script `RPKM_TPM_featurecounts.R` to generate TPM and RPKM counts from featureCounts output:
```bash
# Create output folder
mkdir output/tpm
mkdir output/rpkm

# Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_to_TPMRPKM.sh # 

# mv all output TPM and RPKM quantification to output/tpm or rpkm folder
mv output/featurecounts_hg38/*tpm* output/tpm/
mv output/featurecounts_hg38/*rpkm* output/rpkm/
```
--> Script will output TPM and RPKM counts to `output/tpm`



# DEGs with DESeq2 

Let's do a WT vs KO comparison and identify DEGs using DESEQ2 in R.

--> Prior starting, I would recommend following this **tutorial for learning [R language](https://r4ds.hadley.nz/intro.html)**


## Prerequisite - installation

Create a **deseq2; conda environment**


```bash
conda create -n deseq2 -c conda-forge r-base=4.2.2

# activate conda env
conda activate deseq2

# install extra packages within conda env
conda install -c conda-forge r-ragg 
conda install -c conda-forge r-xml

# Load modules to avoid tidyverse installation zlib issue
#module load zlib

# Load ressources/"computation power"
srun --mem=50g --pty bash -l
```

Open R and install deseq2 (**To open R press `R` and to leave it press `CTRL+D`; NEVER save your workspace**) **--> Install each package one by one! Check error if any; PACKAGE INSTALLATION NEEDS TO BE DONE ONLY ONCE (then, just load your packages)!**:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")


install.packages("tidyverse")
install.packages("pheatmap")
BiocManager::install("apeglm")
install.packages("factoextra")
BiocManager::install("rtracklayer")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("org.Hs.eg.db")
```

- *NOTE: After installing a package in R; **R ask whether you want to update old packages; with choice all/some/none --> Always select none(ie. type n)**  to avoid any issue of package compatibilities*



## If installation failed

### OPTION1 - Install DESEQ2 using conda

```bash
# Activate conda env
conda activate deseq2

# Install DESEQ2 using conda
conda install bioconda::bioconductor-deseq2
```

Open R and try loading DESEQ2

```R
library("DESeq2")


```




### OPTION2 - install from scratch DESEQ2 with new R version
If you encounter error installing DESEQ2 try installing the more recent R version from [DESEQ2 from Bioconductor](https://bioconductor.org/packages/release/bioc/html/DESeq2.html):
```bash
# Create conda env
conda create -n deseq2_v1 -c conda-forge r-base=4.5 # HERE type the last version of R as indicated in Bioconductor

# Activate environment
conda activate deseq2_v1
```

Then go to R and try installing DESeq2:
```R
# Install DESEQ2
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```


### OPTION3 - Use CHOP pre-make environment

Simply load module with R and DESeq2 installed:

```bash
module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1
```

```R
library("DESeq2")
library("tidyverse")
```


--> Good luck :) !


### OPTION4 - Work on Respublica with Rstudio

```bash
conda create -n deseq2_v7 -c conda-forge -c bioconda r-base=4.4 r-ragg r-xml bioconductor-deseq2
conda activate deseq2_v7
conda install -c conda-forge r-tidyverse
conda install -c conda-forge -c bioconda bioconductor-rtracklayer
```


### OPTION5 - Work on Respublica with Rstudio

R4.4.0 work with DESEQ2 and tidyverse







## DESeq2 WT vs KO


The below script will be run in R, in interactive mode; *recommended to add ressources using for example `srun --mem=100g --pty bash -l`*. 

The script will:
- Collect gene symbol name from annotation GTF file
- Use featureCounts counts output and perform DGEs with DESEQ2
- Optional step included to filter out X and Y chromosomes
- Generate volcano plot of DEGs
- output list of DEGs


```bash
# Activate conda environment
conda activate deseq2

# Create output folder
mkdir output/deseq2
```

Before working in R, load computational ressources:
```bash
# Load ressources/"computation power"
srun --mem=50g --pty bash -l

# make sure you are located in your working directory
cd 001*/001*
```


Then open R by typing `R`:

```R
# Load packages (each time you open a new session)
library("rtracklayer")
library("DESeq2")
library("tidyverse")
library("EnhancedVolcano")
library("apeglm")

# --> Make sure all these packages are installed; install as needed using instal.packages() or bioconductor(): you only need to install them once.


# Make sure you are in your working directory
getwd()

# import annotation GTF file to recover gene name, as featurecounts show gene ID ENSG....
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf") # change name accordingly
## Extract geneId and geneSymbol
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
## Rename columns
colnames(gene_table) <- c("geneId", "geneSymbol")

########################################################################
### If import() function is not working; use this workaround: #########
gtf <- read_tsv(
  "../../Master/meta/gencode.v47.annotation.gtf",
  col_names = c(
    "seqname", "source", "feature",
    "start", "end", "score",
    "strand", "frame", "attribute"
  ),
  comment = "#",
  show_col_types = FALSE
)
gene_table <- gtf %>%
  filter(feature == "gene") %>%   # keep only gene entries
  transmute(
    geneId = str_extract(attribute, 'gene_id "[^"]+"') %>%
      str_replace_all('gene_id "|"', ""),
    geneSymbol = str_extract(attribute, 'gene_name "[^"]+"') %>%
      str_replace_all('gene_name "|"', "")
  ) %>%
  distinct() %>%
  drop_na() %>%
  as_tibble()
########################################################################


# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("WT_Rep1", "WT_Rep2", "WT_Rep3", "KO_Rep1", "KO_Rep2", "KO_Rep3")

## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
}

# Merge all dataframe into a single one
counts_all <- reduce(sample_data, full_join, by = "Geneid")

########## OPTIONAL - Remove X and Y chromosome genes #########
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_X_Y <- getBM(attributes = c("ensembl_gene_id"),
                   filters = "chromosome_name",
                   values = c("X", "Y"),
                   mart = ensembl)
counts_all$stripped_geneid <- sub("\\..*", "", counts_all$Geneid)
counts_all_filtered <- counts_all %>%
  filter(!stripped_geneid %in% genes_X_Y$ensembl_gene_id)
counts_all_filtered$stripped_geneid <- NULL
#################################################################


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

# Create colData meta file that describe all our samples
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
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable to avoid extreme FC values
resultsNames(dds) # Here print the output you obtain in the below "coef= ..."
res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm") # install apeglm if needed with `BiocManager::install("apeglm")`

# Add geneSymbol
res_tibble = as_tibble(rownames_to_column(as.data.frame(res), var = "geneId")) %>%
  left_join(gene_table)

# Save the output of the DGEs anlaysis (include gene names, log2FC, pvalue...)
write_tsv(
  res_tibble,
  "output/deseq2/res_with_geneSymbol.txt" # update folder accordingly
)



# Identify DEGs and count them

## padj 0.05 log2FC 0.58 ##################################
res_df <- res_tibble %>% dplyr::select("baseMean", "log2FoldChange", "padj") %>% mutate(padj = ifelse(padj <= 0.05, TRUE, FALSE))
n_upregulated <- sum(res_df$log2FoldChange > 0.58 & res_df$padj == TRUE, na.rm = TRUE)
n_downregulated <- sum(res_df$log2FoldChange < -0.58 & res_df$padj == TRUE, na.rm = TRUE)



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, 'Sky Blue',
    ifelse(res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, 'Orange',
      'grey'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'


pdf("output/deseq2/plotVolcano_res_q05fc058_ESC_KO_vs_ESC_WT_STAR.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res_tibble$geneSymbol,
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
### Complete table with geneSymbol
write.table(res_tibble, file = "output/deseq2/res_ESC_KO_vs_ESC_WT-STAR.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/deseq2/upregulated_q05fc058_ESC_KO_vs_ESC_WT-STAR.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/deseq2/downregulated_q05fc058_ESC_KO_vs_ESC_WT-STAR.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


```

--> Script will output plots and gene list to `output/deseq2`

Once you finish working, you can exit the job and leave ressource by typing `exit`





