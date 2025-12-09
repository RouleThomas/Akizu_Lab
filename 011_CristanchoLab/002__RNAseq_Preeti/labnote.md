# Project

Two cell lines (ReNcell and hiPSCs) in two conditions (24h Normoxia vs 24h Hypoxia) in 4 Bio reps.


--> XXX Directional mRNA library preparation (poly A enrichment), NovaSeq X Plus Series (PE150)


Key 	Sample name 	Belong experiment
A1	ReN_Nor_24h	Exp 3 (Pooled)
A2	ReN_Nor_24h	Exp4 (Pooled)
A3	ReN_Nor_24h	Exp5 (Pooled)
A4	ReN_Nor_24h	Exp6 (Pooled)
A5	ReN_Hyp_24h	Exp 3 (Pooled)
A6	ReN_Hyp_24h	Exp4 (Pooled)
A7	ReN_Hyp_24h	Exp5 (Pooled)
A8	ReN_Hyp_24h	Exp6 (Pooled)
B1	hiPSCs_Nor_24h 	Exp 3 (Pooled)
B2	hiPSCs_Nor_24h 	Exp4 (Pooled)
B3	hiPSCs_Nor_24h 	Exp5 (Pooled)
B4	hiPSCs_Nor_24h 	Exp6 (Pooled)
B5	hiPSCs_Hyp_24h 	Exp 3 (Pooled)
B6	hiPSCs_Hyp_24h 	Exp4 (Pooled)
B7	hiPSCs_Hyp_24h 	Exp5 (Pooled)
B8	hiPSCs_Hyp_24h 	Exp6 (Pooled)





# Data access

Access Cristancho lab folder from CHOP computer: `/mnt/isilon/cristancho_data/Preeti/nsc_rnaseq/downloaded_data/01.RawData` 
Sample name: `/mnt/isilon/cristancho_data/Preeti/nsc_rnaseq/nsc_rna_seq_sample details.xlsx` 






# Pipeline
- Copy data to my `/scr1` working env
- Rename files
- FastQC (fastqc)
- Trimming (fastp)
- Count with featureCounts
- DEGs with DESEQ2





# Cp / import data


```bash
# Data copied from Preeti folder to input_raw/

## cp all data to input folder
cp input_raw/*/*.fq.gz input/

#--> Manually rename all files 
```

--> All good, files succesfully renamed according to `/mnt/isilon/cristancho_data/Preeti/nsc_rnaseq/nsc_rna_seq_sample details.xlsx` with following nomenclature: `[PSC or ReN]_[Norm or Hypo]_[Rep1 - 4]_1.fq`


A1	ReN_Norm_Rep1   ReN_Nor_24h 
A2	ReN_Norm_Rep2   ReN_Nor_24h
A3	ReN_Norm_Rep3   ReN_Nor_24h
A4	ReN_Norm_Rep4   ReN_Nor_24h
A5	ReN_Hypo_Rep1   ReN_Hyp_24h
A6	ReN_Hypo_Rep2   ReN_Hyp_24h
A7	ReN_Hypo_Rep3   ReN_Hyp_24h
A8	ReN_Hypo_Rep4   ReN_Hyp_24h
B1	PSC_Norm_Rep1   hiPSCs_Nor_24h 
B2	PSC_Norm_Rep2   hiPSCs_Nor_24h 
B3	PSC_Norm_Rep3   hiPSCs_Nor_24h 
B4	PSC_Norm_Rep4   hiPSCs_Nor_24h 
B5	PSC_Hypo_Rep1   hiPSCs_Hyp_24h 
B6	PSC_Hypo_Rep2   hiPSCs_Hyp_24h 
B7	PSC_Hypo_Rep3   hiPSCs_Hyp_24h 
B8	PSC_Hypo_Rep4   hiPSCs_Hyp_24h 




# Fastp cleaning

```bash
sbatch scripts/fastp.sh # 60456895 ok
```

# Fastqc 

Let's run fastqc to check why so many reads assigned to no features in featurecounts...


```bash
# Raw fastq
sbatch scripts/fastqc_raw.sh # 60548567 xxx


# fastp-clean fastq
sbatch scripts/fastqc_fastp.sh # 60548570 xxx
```





# STAR mapping fastp trim

```bash
sbatch --dependency=afterany:60456895 scripts/STAR_mapping_fastp.sh # 60457138 ok


## Convert alignment to bigwig
conda activate deeptools
sbatch --dependency=afterany:60457138 scripts/STAR_TPM_bw.sh # 60457411 ok

## Calculate median
conda activate BedToBigwig
sbatch --dependency=afterany:60457411 scripts/bigwigmerge_STAR_TPM_bw.sh # 60458041 ok

```

--> All good





# deepTool bigwig QC - STAR mapping


## TPM bigwig

XXXY HERE DO BIGWIG PCA sbathc job ran already !!!!!!!!!


```bash
conda activate deeptools

###################################
# Include X chr - All samples #####################
###################################
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_STAR_TPM.sh # 60553000 xxx

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
# Include X chr - PSC samples #####################
###################################
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_STAR_TPM-PSC.sh # 60553025 xxx

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
# Include X chr - ReN samples #####################
###################################
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_STAR_TPM-ReN.sh # 60553031 xxx

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



```

--> XXX








# Count with featureCounts


Data look stranded: To confirm with Preeti.

AMPD stringenT: featureCounts -p -C -O \
AMPD relax: featureCounts -p -C -O -M --fraction \


```bash
conda activate featurecounts

# slight test
## -s 2 for stranded
featureCounts -p -C -O -M --fraction -s 2 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 49% assigned

## -s 1 for stranded
featureCounts -p -C -O -M --fraction -s 1 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 16% assigned

## unstranded, not counting multimapped reads
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 50% assigned! mostly not assigned due to multimapping; seems to be stranded in the end

## unstranded, counting multimapped reads
featureCounts -p -C -O -M --fraction \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 60%! Still not great...; the rest fall into unassigned no features

## unstranded, counting multimapped reads, testing previous gene annotations
featureCounts -p -C -O -M --fraction \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 60%! issue does not come from using the new gene annotation; the rest fall into unassigned no features

## count on gene (REMOVE -O), not counting multimapped reads 
featureCounts -p -C -s 2 -t gene -g gene_id \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 48%!

## count on gene (REMOVE -O), counting multimapped reads
featureCounts -p -C -M --fraction -s 2 -t gene -g gene_id \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/ReN_Norm_Rep1.txt output/STAR/fastp/ReN_Norm_Rep1_Aligned.sortedByCoord.out.bam
#--> 66%!




# all samples:
sbatch scripts/featurecounts_unstranded.sh # 60549123 ok
#--> 50-65% uniquely aligned reads - NO, data is stranded
sbatch scripts/featurecounts_multi_unstranded.sh # 60549134 ok
#--> 60-75% uniquely aligned reads - NO, data is stranded
sbatch scripts/featurecounts.sh # 60547513 ok
#--> 50-60% uniquely aligned reads - NO, many multimapped reads!

## The two options below are good:
sbatch scripts/featurecounts_multi.sh # 60547516 ok
#--> 55-65% uniquely aligned reads - YES: output/featurecounts_multi
sbatch scripts/featurecounts_multi_gene.sh # 61302605 xxx
#--> xxx-xxx% uniquely aligned reads - YES: output/featurecounts_multi_gene



```

test with `ReN_Norm_Rep1`
- ~50%, 16% alignment with stranded paramters `-s 2` and `-s 1`, respectively + looking at IGV bam, **files are stranded**
- Many mapping to *unassigned features*, and *multimapped reads*.. Let's generate two versions:
  - one counting multimapped reads (`*_multi`), and another one removing multimapped reads, improve a bit but not so much
- I noticed checking bigwig on IGV that many reads fall within introns, should explain the high *unassigned features*
  - For this let's do a version where I **count on gene, not exon**; **when counting on gene level do NOT use` -O`**; indeed this count reads in gene that overlap twice, but it will lead to artificial correlation and issue with DESEQ2 when counting on gene body...


--> Seems **count on exon, stranded, and count multimapped reads is the best approach**: `output/featurecounts_multi` (count on gene as also been generated at `output/featurecounts_multi_gene`)






## Calculate TPM and RPKM


Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
conda activate deseq2
# output/featurecounts_multi
## Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_TPM.sh # 61302936 ok
## mv all output to output/tpm or rpkm folder
mv output/featurecounts_multi/*tpm* output/tpm_featurecounts_multi/
mv output/featurecounts_multi/*rpkm* output/rpkm_featurecounts_multi/


# output/featurecounts_multi_gene
## Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_multi_gene_TPM.sh # 61303337
## mv all output to output/tpm or rpkm folder
mv output/featurecounts_multi_gene/*tpm* output/tpm_featurecounts_multi_gene/
mv output/featurecounts_multi_gene/*rpkm* output/rpkm_featurecounts_multi_gene/
```

--> All good. 





XXXY HERE !!!


# DEGs with deseq2 (featurecounts)

**IMPORTANT NOTE: Here it is advisable to REMOVE all genes from chromosome X and Y BEFORE doing the DEGs analysis (X chromosome re-activation occurs in some samples, notably these with more cell passage; in our case, the HET and KO)**
--> It is good to do this on the count matrix see [here](https://support.bioconductor.org/p/119932/)
### 'one-by-one' comparison
Comparison WT vs mutant:
- ESc KO vs WT
- ESC OEKO vs WT



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
library("rtracklayer")


# import GTF for gene name
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
## Extract geneId and geneSymbol
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
## Rename columns
colnames(gene_table) <- c("geneId", "geneSymbol")





# import featurecounts output and keep only gene ID and counts
## collect all samples ID
samples <- c("ESC_WT_R1", "ESC_WT_R2" ,"ESC_WT_R3" ,"ESC_KO_R1" ,"ESC_KO_R2", "ESC_KO_R3")

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


# Add geneSymbol
res_tibble = as_tibble(rownames_to_column(as.data.frame(res), var = "geneId")) %>%
  left_join(gene_table)


# Identify DEGs and count them

## padj 0.05 FC 0.58 ##################################
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




