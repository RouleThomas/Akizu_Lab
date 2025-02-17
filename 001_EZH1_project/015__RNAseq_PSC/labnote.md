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

# Count with Salmon

## install Salmon

Let's follow [this](https://combine-lab.github.io/salmon/getting_started/)



```bash
conda create -n salmon salmon
# --> bug salmon: error while loading shared libraries: libboost_thread.so.1.60.0: cannot open shared object file: No such file or directory 

conda create -n salmon1 -c conda-forge -c bioconda salmon=1.3.0

```

- *NOTE: using conda, it bug 1st, issue discussed [here](https://github.com/COMBINE-lab/salmon/issues/565); installed v1.3.0 make it work!*

--> all good


## count with salmon



Counting is perform on the transcriptome, 1st downlaod and import the transcriptome to `../../Master/meta/salmon` . Download frmop [ensembl](https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/)


```bash
conda activate salmon1

# build transcriptome index
salmon index -t ../../Master/meta/salmon/Homo_sapiens.GRCh38.cdna.all.fa.gz -i ../../Master/meta/salmon/Homo_sapiens

# perform counting
## light testing
salmon quant -i ../../Master/meta/salmon/Homo_sapiens -l A \
         -1 output/fastp/PSC_WT_R1_1.fq.gz \
         -2 output/fastp/PSC_WT_R1_2.fq.gz \
         -p 8 --validateMappings -o salmon/PSC_WT_R1_quant

## run in sbatch
sbatch scripts/salmon_count_1.sh # 32405964 ok
sbatch scripts/salmon_count_2.sh # 32406155 ok
sbatch scripts/salmon_count_3.sh # 32406159 ok

## Run with gtf; recommended when using isoformSwitchAnalyzeR - Canonical chr
sbatch scripts/salmon_count_gtf.sh # 32547743 ok

## Run with gtf; recommended when using isoformSwitchAnalyzeR - All chr
sbatch scripts/salmon_count_gtf_allchr.sh # 32551178 xxx
```
- *NOTE: index adapted for 75bp or longer with default `-k 31`, good for us as 150bp*

--> ~>90% mapping; for gtf version too


# Count with Kallisto

Follow instruction [here](https://pachterlab.github.io/kallisto/download.html)

## install Kallisto


```bash
conda create -n kallisto -c bioconda -c conda-forge kallisto
```

## count with Kallisto

```bash
conda activate kallisto

# build transcriptome index
kallisto index -i ../../Master/meta/kallisto/transcripts.idx ../../Master/meta/salmon/Homo_sapiens.GRCh38.cdna.all.fa.gz # seems it cannot use a gtf

# perform counting
## light testing
kallisto quant -i ../../Master/meta/kallisto/transcripts.idx -o output/kallisto -b 100 output/fastp/PSC_WT_R1_1.fq.gz output/fastp/PSC_WT_R1_2.fq.gz -g ../../Master/meta/gencode.v47.annotation.gtf

## run in sbatch
sbatch scripts/kallisto_count_gtf_1.sh # 32559439 ok
sbatch scripts/kallisto_count_gtf_2.sh # 32559511 ok
sbatch scripts/kallisto_count_gtf_3.sh # 32559512 ok

```



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

- Sample from this *015/Jasmine* experiment (`multiBigwigSummary_all`)
- Sample between this *015/Jasmine* and *001/Carolina* experiment (multiBigwigSummary_001015)


```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all.sh # 28990409 ok
sbatch scripts/multiBigwigSummary_001015.sh # 35427778 ok
sbatch scripts/multiBigwigSummary_001015_WTKO.sh # 35429121 ok


# Plot all
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 KOEF1aEZH1_R1 KOEF1aEZH1_R2 KOEF1aEZH1_R3 \
    --colors black black black red red red blue blue blue \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf
plotPCA -in output/bigwig/multiBigwigSummary_001015.npz \
    --transpose \
    --ntop 0 \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 KOEF1aEZH1_R1 KOEF1aEZH1_R2 KOEF1aEZH1_R3 WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 \
    --colors black black black red red red blue blue blue grey grey grey maroon maroon maroon \
    -o output/bigwig/multiBigwigSummary_001015_plotPCA.pdf
plotPCA -in output/bigwig/multiBigwigSummary_001015_WTKO.npz \
    --transpose \
    --ntop 0 \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 \
    --colors black black black red red red grey grey grey maroon maroon maroon \
    -o output/bigwig/multiBigwigSummary_001015_WTKO_plotPCA.pdf


## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 KOEF1aEZH1_R1 KOEF1aEZH1_R2 KOEF1aEZH1_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_001015.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 KOEF1aEZH1_R1 KOEF1aEZH1_R2 KOEF1aEZH1_R3 WT_R1 WT_R2 WT_R3 KO_R1 KO_R2 KO_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_001015_heatmap.pdf


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

  
pdf("output/deseq2/plotVolcano_res_q05fc025_PSC_KO_vs_PSC_WT_EZH1.pdf", width=7, height=8)  
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  selectLab = "EZH1",
  drawConnectors = TRUE,
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


pdf("output/deseq2/plotVolcano_res_q05fc025_PSC_KOEF1aEZH1_vs_PSC_WT_EZH1.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  selectLab = "EZH1",
  drawConnectors = TRUE,
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




### PSC KOEF1aEZH1 vs KO

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
samples <- c("PSC_KOEF1aEZH1_R1" ,"PSC_KOEF1aEZH1_R2", "PSC_KOEF1aEZH1_R3", "PSC_KO_R1", "PSC_KO_R2" ,"PSC_KO_R3")

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
dds$genotype <- relevel(dds$genotype, ref = "KO")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_KOEF1aEZH1_vs_KO", type="apeglm")




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

pdf("output/deseq2/plotVolcano_res_q05fc025_PSC_KOEF1aEZH1_vs_PSC_KO.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs KO, PSC',
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
  annotate("text", x = 3, y = 40, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 40, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2/res_PSC_KOEF1aEZH1_vs_PSC_KO.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.25 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_PSC_KOEF1aEZH1_vs_PSC_KO.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_PSC_KOEF1aEZH1_vs_PSC_KO.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





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

pdf("output/deseq2/plotVolcano_res_q05fc05_PSC_KOEF1aEZH1_vs_PSC_KO.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs KO, PSC',
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
  annotate("text", x = 3, y = 40, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 40, 
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
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_KO.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc05_PSC_KOEF1aEZH1_vs_PSC_KO.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```






### PSC WT(015) vs WT(001)

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



# SAMPLE 015 #########################
# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("PSC_WT_R1", "PSC_WT_R2" ,"PSC_WT_R3")
## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
}
# Merge all dataframe into a single one
counts_all_015 <- reduce(sample_data, full_join, by = "Geneid")
# SAMPLE 001 #########################
samples <- c("ESC_WT_R1", "ESC_WT_R2" ,"ESC_WT_R3")
## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}
# Merge all dataframe into a single one
counts_all_001 <- reduce(sample_data, full_join, by = "Geneid")
# MERGE 001 015
counts_all = counts_all_015 %>%
  left_join(counts_all_001)
# RENAME WT
colnames(counts_all) <- colnames(counts_all) %>%
  gsub("^PSC_WT", "PSC_WTS", .) %>%  # Replace WT with WTS for PSC_
  gsub("^ESC_WT", "ESC_H9", .)      # Replace WT with H9 for ESC_
samples <- c("PSC_WTS_R1", "PSC_WTS_R2" ,"PSC_WTS_R3", "ESC_H9_R1", "ESC_H9_R2" ,"ESC_H9_R3")


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
dds$genotype <- relevel(dds$genotype, ref = "H9")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below
res <- lfcShrink(dds, coef="genotype_WTS_vs_H9", type="apeglm")


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

pdf("output/deseq2/plotVolcano_res_q05fc025_PSC_WTS_vs_PSC_H9.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'WTS vs H9, PSC',
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
  annotate("text", x = 3, y = 300, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 300, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2/res_PSC_WTS_vs_PSC_H9.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.25 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_PSC_WTS_vs_PSC_H9.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_PSC_WTS_vs_PSC_H9.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





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

pdf("output/deseq2/plotVolcano_res_q05fc05_PSC_WTS_vs_PSC_H9.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'WTS vs H9, PSC',
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
  annotate("text", x = 3, y = 280, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 280, 
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
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc05_PSC_WTS_vs_PSC_H9.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc05_PSC_WTS_vs_PSC_H9.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```





### PSC KO(015) vs WT(001)

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


# SAMPLE 001 #########################
samples <- c("ESC_WT_R1", "ESC_WT_R2" ,"ESC_WT_R3")
## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}
# Merge all dataframe into a single one
counts_all_001 <- reduce(sample_data, full_join, by = "Geneid")

# SAMPLE 015 #########################
# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("PSC_KO_R1", "PSC_KO_R2" ,"PSC_KO_R3")
## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
}
# Merge all dataframe into a single one
counts_all_015 <- reduce(sample_data, full_join, by = "Geneid")

# MERGE 001 015
counts_all = counts_all_001  %>%
  left_join(counts_all_015)

samples <- c("ESC_WT_R1", "ESC_WT_R2" ,"ESC_WT_R3", "PSC_KO_R1", "PSC_KO_R2" ,"PSC_KO_R3")


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

pdf("output/deseq2/plotVolcano_res_q05fc025_PSC_KO015_vs_PSC_WT001.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs H9, PSC',
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
  annotate("text", x = 3, y = 300, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 300, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2/res_PSC_KO015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.25 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_PSC_KO015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_PSC_KO015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





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

pdf("output/deseq2/plotVolcano_res_q05fc05_PSC_KO015_vs_PSC_WT001.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs H9, PSC',
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
  annotate("text", x = 3, y = 300, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 300, 
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
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc05_PSC_KO015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc05_PSC_KO015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```




### PSC KOEF1aEZH1(015) vs WT(001)

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


# SAMPLE 001 #########################
samples <- c("ESC_WT_R1", "ESC_WT_R2" ,"ESC_WT_R3")
## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("../001__RNAseq/output/featurecounts_hg38/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR_hg38/")) %>%
    rename(!!sample := starts_with("output/STAR_hg38/"))
}
# Merge all dataframe into a single one
counts_all_001 <- reduce(sample_data, full_join, by = "Geneid")

# SAMPLE 015 #########################
# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("PSC_KOEF1aEZH1_R1", "PSC_KOEF1aEZH1_R2" ,"PSC_KOEF1aEZH1_R3")
## Make a loop for importing all featurecounts data and keep only ID and count column
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/featurecounts/", sample, ".txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
    dplyr::select(Geneid, starts_with("output/STAR/")) %>%
    rename(!!sample := starts_with("output/STAR/"))
}
# Merge all dataframe into a single one
counts_all_015 <- reduce(sample_data, full_join, by = "Geneid")

# MERGE 001 015
counts_all = counts_all_001  %>%
  left_join(counts_all_015)

samples <- c("ESC_WT_R1", "ESC_WT_R2" ,"ESC_WT_R3", "PSC_KOEF1aEZH1_R1", "PSC_KOEF1aEZH1_R2" ,"PSC_KOEF1aEZH1_R3")


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

pdf("output/deseq2/plotVolcano_res_q05fc025_PSC_KOEF1aEZH1015_vs_PSC_WT001.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs H9, PSC',
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
  annotate("text", x = 3, y = 300, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 300, 
           label = paste(n_downregulated), hjust = 0, size = 6, color = "darkred")
dev.off()



# Save as gene list for GO analysis:
### Complete table with GeneSymbol
write.table(res, file = "output/deseq2/res_PSC_KOEF1aEZH1015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, row.names = TRUE) # that is without X and Y chr genes
### GO EntrezID Up and Down
#### Filter for up-regulated genes
upregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 0.25 & res$padj < 5e-2, ]

#### Filter for down-regulated genes
downregulated <- res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -0.25 & res$padj < 5e-2, ]
#### Save
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc025_PSC_KOEF1aEZH1015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc025_PSC_KOEF1aEZH1015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)





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

pdf("output/deseq2/plotVolcano_res_q05fc05_PSC_KOEF1aEZH1015_vs_PSC_WT001.pdf", width=7, height=8)    
EnhancedVolcano(res,
  lab = res$GeneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KOEF1aEZH1 vs H9, PSC',
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
  annotate("text", x = 3, y = 300, 
           label = paste(n_upregulated), hjust = 1, size = 6, color = "darkred") +
  annotate("text", x = -3, y = 300, 
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
write.table(upregulated$GeneSymbol, file = "output/deseq2/upregulated_q05fc05_PSC_KOEF1aEZH1015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$GeneSymbol, file = "output/deseq2/downregulated_q05fc05_PSC_KOEF1aEZH1015_vs_PSC_WT001.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```










--> WTS (`015/Jasmine`) vs H9 (`001/Carolina`) are VERY different, tons of DEGs!!!


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



# Investigate DEGs 001015 (Jasmine) and 001001 (Carolina)

Related to `015__RNAseq_PSC_v2_meeting20250122.pptx`; check the genomic distribution of the genes similarly deregulated in KO and KOEF1aEZH1. If they are all positioned within the same genomic region that could suggest CNV occuring during mutant generations



```bash
conda activate deseq2
```

```R
# Packages
library("ggplot2")
library("dplyr")
library("readr")
library("tidyverse")

# Import the gene list file
gene_list <- read.table("output/deseq2/vennDiagram_overlapUpKOUpKOEF1aEZH1_256genes.txt", header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()
gene_list <- read.table("output/deseq2/vennDiagram_overlapDownKODownKOEF1aEZH1_143genes.txt", header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()
# Import the BED file
gene_bed <- read.table("../../Master/meta/ENCFF159KBI_geneSymbol.bed", header = FALSE, stringsAsFactors = FALSE, sep = "\t") %>% as_tibble()
colnames(gene_bed) <- c("chromosome", "start", "end", "gene", "strand")

# Merge the files based on gene symbol
merged_data <- gene_list %>%
  inner_join(gene_bed, by = "gene")

# Prepare data for plotting
# Order chromosomes correctly (e.g., chr1, chr2, ..., chrX)
merged_data$chromosome <- factor(merged_data$chromosome, levels = unique(gene_bed$chromosome))

# Generate the plot
pdf("output/deseq2/genomicLocation-vennDiagram_overlapDownKODownKOEF1aEZH1_143genes.pdf", width=7, height=6)    
ggplot(merged_data, aes(x = chromosome, y = start)) +
  geom_point() +
  labs(
    title = "Distribution of Genes Along the Genome",
    x = "Chromosome",
    y = "Genomic Position (Coordinate)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()

```


--> DEGs do not localized at same place within the genome; likely no CNV happeniong in KO and KOEF1aEZH1




# Differential alternative mRNA splicing

Let's follow method employed in [Jhanji et al 2024](https://www.biorxiv.org/content/10.1101/2024.12.02.625500v1.full.pdf); briefly:
- Isoform switching between WT and KO with *DEXSeq* method implemented in the *IsoformSwitchAnalyzeR* package (identifies bins with differential isoform switching between conditions, quantified as the difference in isoform fraction (dIF), calculated as IF_mutant - IF_control)
- Significant isoform switching was defined by a dIF > 0.05 and a false discovery rate (FDR) < 0.05


## IsoformSwitchAnalyzeR installation in R

```bash
conda activate deseq2V3
```


```R
BiocManager::install("IsoformSwitchAnalyzeR")
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR", build_vignettes = TRUE)
#--> Fail with pfamAnalyzerR after


```
--> Was good but then lead to a bug. Try reinstall dplyr or the dev version but fail; lets install R 4.3.0 and install the dev version from this

Let's try to copy scRNAseq environment who is R4.3.0



```bash
conda create --name IsoformSwitchAnalyzeR --clone scRNAseq
conda activate IsoformSwitchAnalyzeR
```



```R
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install(version='devel')
BiocManager::install("IsoformSwitchAnalyzeR")
```
--> Fail R 4.3.0 required


```bash
conda create -n IsoformSwitchAnalyzeRv1 -c conda-forge -c bioconda r-base=4.3.0
conda activate IsoformSwitchAnalyzeRv1
```


```R

if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
if (!requireNamespace("pfamAnalyzeR", quietly = TRUE)){
    devtools::install_github("kvittingseerup/pfamAnalyzeR")
}
#--> textshaping configuration failed
if (!requireNamespace("devtools", quietly = TRUE)){
install.packages("devtools")
}
devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR", build_vignettes = TRUE)
#--> pkgdown fail
```
--> Need R 4.5 devel..

Try R 4.4.2 (last version)


```bash
conda create -n IsoformSwitchAnalyzeRv2 -c conda-forge -c bioconda r-base=4.4.2
conda activate IsoformSwitchAnalyzeRv2
```

```R
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
if (!requireNamespace("pfamAnalyzeR", quietly = TRUE)){
    devtools::install_github("kvittingseerup/pfamAnalyzeR")
}
#--> Error dependied miniUI, pkgdown, roxygen2, rversions, urlchecker
install.packages("miniUI")
#-->fail
```

Try install developer version of R, follow [this](https://cran.r-project.org/doc/manuals/r-devel/R-admin.html)


```bash
conda create -n IsoformSwitchAnalyzeRv3 
conda activate IsoformSwitchAnalyzeRv3
```

Download [here](https://cran.r-project.org/src/base-prerelease/)  	`R-devel.tar.gz` and trasnfer to `/Master/software`

```bash
cd ../../Master/software
#  Extract the Source File
tar -xzvf R-devel.tar.gz
cd R-devel

./configure 
#--> Error with X11
./configure LIBnn=lib
#--> Add no-X option
./configure --prefix=/scr1/users/roulet/Akizu_Lab/Master/software/R-devel-install --with-x=no --enable-R-shlib LIBnn=lib
#-->error libcurl
conda install -c conda-forge libcurl
./configure --prefix=/scr1/users/roulet/Akizu_Lab/Master/software/R-devel-install --with-x=no --enable-R-shlib LIBnn=lib
#--> look good
make
# --> work!

```

Now switch to R, for that type: `/scr1/users/roulet/Akizu_Lab/Master/software/R-devel/bin/R` to run *R developer* version `R4.5.0`

```R
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
if (!requireNamespace("pfamAnalyzeR", quietly = TRUE)){
    devtools::install_github("kvittingseerup/pfamAnalyzeR")
}
#--> bug with pkgdown; lets try instlling 
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install(version='devel')
BiocManager::install("IsoformSwitchAnalyzeR")
#--> IsoformSwitchAnalyzeR nota avail for Bioconductor 3.21
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR", build_vignettes = TRUE)
#--> error textshaping
# Back to conda and install manually `conda install -c conda-forge harfbuzz fribidi cairo pango libpng libtiff freetype`
install.packages("textshaping", dependencies = TRUE)
install.packages("ragg", dependencies = TRUE)
```
--> Still fail...



Let's try to copy `ChIPseqSpikeInFree conda env` which have devtools and R/3.6.1 installed...

```bash
conda create --name IsoformSwitchAnalyzeRv4 --clone ChIPseqSpikeInFree
#--> Bug, conda env corrupted! 
conda activate IsoformSwitchAnalyzeRv4
```

Let's try to copy  `Signac_Pando` which have devtools and R/4.3.3 installed... 

```bash
conda create --name IsoformSwitchAnalyzeRv5 --clone Signac_Pando

conda activate IsoformSwitchAnalyzeRv5
```

```R
devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR", build_vignettes = TRUE)
#--> fail

BiocManager::install("IsoformSwitchAnalyzeR")
```
--> Work!!


If fail, try find conda env or docker to use IsoformSwitchAnalyzeR...






Let's try install through [conda](https://anaconda.org/bioconda/bioconductor-isoformswitchanalyzer)



```bash
conda create -n IsoformSwitchAnalyzeRv6 -c bioconda -c conda-forge bioconductor-isoformswitchanalyzer
#--> Fail, try another one
conda create -n IsoformSwitchAnalyzeRv6 -c bioconda -c conda-forge bioconda/label/cf201901::bioconductor-isoformswitchanalyzer
#-->

```









## IsoformSwitchAnalyzeR usage

I follow tutorial from [here](https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html).

We need to re-do the quantification, let's use Salmon or Kallisto. Also it seems doing transcript level quantification and afterwards summarizing to gene level is much better! so let's do this!





```bash
conda activate IsoformSwitchAnalyzeRv5
```

```R
# packages
library("IsoformSwitchAnalyzeR")


# Salmon ####################################################


# Importing the Data
salmonQuant <- importIsoformExpression(
    parentDir = "output/salmon_gtf_allchr/")

# metadata file
myDesign = data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('.*_(WT|KO|KOEF1aEZH1)_.*', '\\1', colnames(salmonQuant$abundance)[-1])
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
    pathToOutput = 'output/IsoformSwitchAnalyzeR',
    outputSequences      = TRUE, # change to TRUE whan analyzing your own data 
    prepareForWebServers = TRUE  # change to TRUE if you will use webservers for external sequence analysis
)

#--> Run in WebServer the CPC2, PFAM, IUPRED2A, SIGNALP

analysSwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList, 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "output/IsoformSwitchAnalyzeR/result_cpc2.txt",
  pathToPFAMresultFile      = "output/pfam/pfam_results_reformat.txt",
  pathToIUPred2AresultFile  = "output/IsoformSwitchAnalyzeR/isoformSwitchAnalyzeR_isoform_AA_complete.result",
  pathToSignalPresultFile   = "output/IsoformSwitchAnalyzeR/prediction_results_SignalIP6.txt",
  outputPlots               = TRUE
)

## SAVE IMAGE R SESSION
# save.image("IsoformSwitchAnalyzeR_v1.RData")
# load("IsoformSwitchAnalyzeR_v1.RData")
##


## Generate plot for a gene

pdf(file = '<outoutDirAndFileName>.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(switchAnalyzeRlist, gene='<gene_name>')
dev.off()




# Genome-wide Summaries
pdf(file = 'output/IsoformSwitchAnalyzeR/extractSwitchOverlap.pdf', onefile = FALSE, height=6, width = 9)
extractSwitchOverlap(
    analysSwitchList,
    filterForConsequences=TRUE,
    plotIsoforms = FALSE
)
dev.off()


pdf(file = 'output/IsoformSwitchAnalyzeR/extractConsequenceSummary.pdf', onefile = FALSE, height=6, width = 9)
extractConsequenceSummary(
    analysSwitchList,
    consequencesToAnalyze='all',
    plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)
dev.off()

# Consequence Enrichment Analysis
pdf(file = 'output/IsoformSwitchAnalyzeR/extractConsequenceEnrichment.pdf', onefile = FALSE, height=6, width = 12)
extractConsequenceEnrichment(
    analysSwitchList,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()



# Splicing Enrichment Analysis

pdf(file = 'output/IsoformSwitchAnalyzeR/extractSplicingEnrichment.pdf', onefile = FALSE, height=6, width = 12)
extractSplicingEnrichment(
    analysSwitchList,
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off()

#Overview Plots
## Volcano like plot

pdf(file = 'output/IsoformSwitchAnalyzeR/Overview_Plots.pdf', onefile = FALSE, height=3, width = 6)
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

pdf(file = 'output/IsoformSwitchAnalyzeR/Overview_Plots3.pdf', onefile = FALSE, height=6, width = 6)
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
### Count for WT vs KOEF1aEZH1
WT_vs_KOEF1aEZH1 <- significant_isoforms %>%
  filter(condition_1 == "KOEF1aEZH1" & condition_2 == "WT") %>%
  nrow()
### Print the counts
cat("Number of significant isoform switches (WT vs KO):", WT_vs_KO, "\n")
cat("Number of significant isoform switches (WT vs KOEF1aEZH1):", WT_vs_KOEF1aEZH1, "\n")


### Switch vs Gene changes:
pdf(file = 'output/IsoformSwitchAnalyzeR/Overview_Plots2.pdf', onefile = FALSE, height=3, width = 6)
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









# Kallisto ####################################################

# Importing the Data
salmonQuant <- importIsoformExpression(
    parentDir = "output/kallisto/")

# metadata file
myDesign = data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('.*_(WT|KO|KOEF1aEZH1)_.*', '\\1', colnames(salmonQuant$abundance)[-1])
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
# save.image("IsoformSwitchAnalyzeR_v2_kallisto.RData")
# load("IsoformSwitchAnalyzeR_v2_kallisto.RData")
##


## Generate plot for a gene

pdf(file = 'output/IsoformSwitchAnalyzeR_kallisto/switchPlot_WTKO_CRBN.pdf', onefile = FALSE, height=6, width = 9)
switchPlot(analysSwitchList, gene= "CRBN", condition1= "WT", condition2= "KOEF1aEZH1")
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
### Count for WT vs KOEF1aEZH1
WT_vs_KOEF1aEZH1 <- significant_isoforms %>%
  filter(condition_1 == "KOEF1aEZH1" & condition_2 == "WT") %>%
  nrow()
### Print the counts
cat("Number of significant isoform switches (WT vs KO):", WT_vs_KO, "\n")
cat("Number of significant isoform switches (WT vs KOEF1aEZH1):", WT_vs_KOEF1aEZH1, "\n")


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



- NOTE: Issue with importRdata(); error in dplyr::full_join(), seems need v 4.3.0 of R.... Issue discussed [here](https://github.com/kvittingseerup/IsoformSwitchAnalyzeR/issues/189). or install a more recent version of `IsoformSwitchAnalyzeR` discussed [here](https://github.com/kvittingseerup/IsoformSwitchAnalyzeR/issues/168); 
  - also discuss to try re-installing dplyr, so I did `install.packages("dplyr")` and `library(dplyr)` and it fail
  - version of my IsoformSwitchAnalyzeR is 1.20.0


- NOTE: For the gtf, they recommend using ENCODE (not Ensembl) and with haplotypes (all non canonical chr); I downloaded from [here](https://www.gencodegenes.org/human/); and transfer to Master/meta/ ; file to use= `gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz`


- NOTE: Got a warning indicated5719 isoforms only found in the annotation and that could cause discrepancies; tested Salmon quantification using gtf (canonical chr, and all chr); and same Warning message. Let's ignore as they also say it could be trigger because I used Salmon
  --> **Let's use Salmon quantification and IsoformSwitchAnalyzeR with the gtf with all chr (canonical and non canonical)**
    --> **Kallisto do not show error message, better be sage than sorry so let's use kallisto!**



--> FASTA Files at `output/IsoformSwitchAnalyzeR`, can be used *to perform the external analysis tools (Pfam (for prediction of protein domains), SignalP (for prediction of signal peptides), CPAT or CPC2 (for prediction of coding potential, since they perform similar analysis so it is only necessary to run one of them)), Prediction of sub-cellular location (DeepLoc2) and isoform toplogy (DeepTmHmm), IUPred2A or NetSurfP-2 (for prediction of Intrinsically Disordered Regions (IDR) )*


Let's run these in Webserver :
- coding potential with [CPC2](https://cpc2.gao-lab.org/), I put `output/IsoformSwitchAnalyzeR/isoformSwitchAnalyzeR_isoform_nt.fasta`;
  - salmon: `IsoformSwitchAnalyzeR/result_cpc2.txt` --> LOOK GOOD
  - kallisto: `IsoformSwitchAnalyzeR_kallisto/result_cpc2.txt`
- protein domain with [PFAM](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan), run in code below at `### PFAM`;
  - salmon: `pfam_results.txt` --> LOOK GOOD; but need reformat (see below custom python script)
  - kallisto: `pfam_results_kallisto.txt` --> LOOK GOOD; but need reformat (see below custom python script)
- Prediction of Intrinsically Unstructured Proteins with [IUPred2](https://iupred2a.elte.hu/); V3 do not support multi FASTA file;
  - salmon: `IsoformSwitchAnalyzeR/isoformSwitchAnalyzeR_isoform_AA_complete.result` --> LOOK GOOD
  - kallisto: `IsoformSwitchAnalyzeR_kallisto/isoformSwitchAnalyzeR_isoform_AA_complete.result`
- Prediction of signal peptide with [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) with option Eukarya/short output/fast; and save the *Prediction summary*
  - salmon: `IsoformSwitchAnalyzeR/prediction_results_SignalIP6.txt` --> LOOK GOOD
  - kallisto: `IsoformSwitchAnalyzeR_kallisto/prediction_results_SignalIP6.txt` --> LOOK GOOD



### PFAM

Follow instructions [here](https://github.com/aziele/pfam_scan).
--> In the end I ran `pfam_scan.py` and reformat the output with a custom python script so that it match the output of the `pfam_scan.pl` command (which I couldn't use due to perl dependecies).

```bash
conda activate deseq2
module load HMMER


# download and prepare database
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

mkdir pfamdb
gunzip -c Pfam-A.hmm.dat.gz > pfamdb/Pfam-A.hmm.dat
gunzip -c Pfam-A.hmm.gz > pfamdb/Pfam-A.hmm
rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz

hmmpress pfamdb/Pfam-A.hmm

# install PFAM
cd ../../Master/software
git clone https://github.com/aziele/pfam_scan
cd pfam_scan
./pfam_scan.py --help

# run PFAM
cd pfam_scan
## salmon
./pfam_scan.py ../../../001_EZH1_Project/015__RNAseq_PSC/output/IsoformSwitchAnalyzeR/isoformSwitchAnalyzeR_isoform_AA_complete.fasta ../pfamdb/ -out ../../../001_EZH1_Project/015__RNAseq_PSC/output/pfam/pfam_results.txt
## kallisto
./pfam_scan.py ../../../001_EZH1_Project/015__RNAseq_PSC/output/IsoformSwitchAnalyzeR_kallisto/isoformSwitchAnalyzeR_isoform_AA_complete.fasta ../pfamdb/ -out ../../../001_EZH1_Project/015__RNAseq_PSC/output/pfam/pfam_results_kallisto.txt
```


--> The output is weird, not as in the tutorial; lets try using the `pfam_scan.pl` [script](https://github.com/SMRUCC/GCModeller/blob/master/src/interops/scripts/PfamScan/PfamScan/pfam_scan.pl). Please see after, need reformating, a **custom python script has been created.**


```bash
cd ../../Master/software
nano pfam_scan.pl

# make script exectubale
chmod +x pfam_scan.pl
# run it
./pfam_scan.pl -fasta ../../001_EZH1_Project/015__RNAseq_PSC/output/IsoformSwitchAnalyzeR/isoformSwitchAnalyzeR_isoform_AA_complete.fasta -dir pfamdb/ -outfile ../../001_EZH1_Project/015__RNAseq_PSC/output/pfam/pfam_results_script.txt

```

--> Issue running this, so instead let's use a **custom python script to reformat my `./pfam_scan.py`**

```bash
# salmon
nano scripts/reformat_pfam.py
python3 scripts/reformat_pfam.py
# kallisto
nano scripts/reformat_pfam_kallisto.py
python3 scripts/reformat_pfam_kallisto.py
```
--> Works!

