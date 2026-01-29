--> This folder is the **clean + reproducible** version of: `001_CristanchoLab/002__RNAseq_Preeti` (original exploration)

# Data Overview

- **Data type**: Bulk RNA-seq
- **Library**: Directional mRNA (poly(A) enrichment)
- **Platform**: Illumina NovaSeq X Plus (paired-end, 150 bp)
- **Factors**:
  - **Cell type**: ReN, PSC
  - **Condition**: Normoxia (24h), Hypoxia (24h)
  - **Replicates**: 4 biological replicates per group (total n = 16)



# Data processing

## Sample Renaming

| Key | Original sample name | New sample name |
| --: | -------------------- | --------------- |
|  A1 | ReN_Nor_24h_Exp3     | ReN_Norm_Rep1   |
|  A2 | ReN_Nor_24h_Exp4     | ReN_Norm_Rep2   |
|  A3 | ReN_Nor_24h_Exp5     | ReN_Norm_Rep3   |
|  A4 | ReN_Nor_24h_Exp6     | ReN_Norm_Rep4   |
|  A5 | ReN_Hyp_24h_Exp3     | ReN_Hypo_Rep1   |
|  A6 | ReN_Hyp_24h_Exp4     | ReN_Hypo_Rep2   |
|  A7 | ReN_Hyp_24h_Exp5     | ReN_Hypo_Rep3   |
|  A8 | ReN_Hyp_24h_Exp6     | ReN_Hypo_Rep4   |
|  B1 | hiPSCs_Nor_24h_Exp3  | PSC_Norm_Rep1   |
|  B2 | hiPSCs_Nor_24h_Exp4  | PSC_Norm_Rep2   |
|  B3 | hiPSCs_Nor_24h_Exp5  | PSC_Norm_Rep3   |
|  B4 | hiPSCs_Nor_24h_Exp6  | PSC_Norm_Rep4   |
|  B5 | hiPSCs_Hyp_24h_Exp3  | PSC_Hypo_Rep1   |
|  B6 | hiPSCs_Hyp_24h_Exp4  | PSC_Hypo_Rep2   |
|  B7 | hiPSCs_Hyp_24h_Exp5  | PSC_Hypo_Rep3   |
|  B8 | hiPSCs_Hyp_24h_Exp6  | PSC_Hypo_Rep4   |



## Read pre-processing & Quality filtering

Reads were trimmed using **fastp (VERSION)** with default adapter detection and quality filtering parameters.

The following command was run for each (paired-end) sample:
```bash
fastp \
  -i <sample>_1.fq.gz \
  -I <sample>_2.fq.gz \
  -o <sample>_1.fq.gz \
  -O <sample>_2.fq.gz \
  -j <sample>.json \
  -h <sample>.html
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/fastp.sh`


## Read mapping

Reads were aligned to the human reference genome (**hg38**, *GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta*) using **STAR** (v2.7.3a).  
Resulting BAM files were sorted by coordinate and indexed using **SAMtools** (v1.16.1).

The following command was run for each (paired-end) sample:
```bash
STAR \
  --genomeDir <genome_index> \
  --readFilesIn <sample>_1.fq.gz <sample>_2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate

samtools index <sample>_Aligned.sortedByCoord.out.bam
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/STAR_mapping_fastp.sh`


## Gene-level quantification

Gene-level read counts were generated using **featureCounts** (v2.0.1) with the **GENCODE v47** gene annotation.

The following command was run for each sample:
```bash
featureCounts -p -C -O -M --fraction -s 2
	-a <annotation.gtf> \
	-o <sample>.txt <sample>_Aligned.sortedByCoord.out.bam
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/featurecounts_multi.sh`



## Isoform-level quantification


Isoform-level read counts were generated using **kallisto** (v0.44.0) with the **GENCODE v47** gene annotation.


The following command was run for each sample:
```bash
kallisto quant \
  -i <transcripts.idx> \
  -o <sample>_quant \
  -b 100 \
  -g <annotation.gtf> \
  --rf-stranded \
  --genomebam \
  --chromosomes GRCh38_chrom_sizes.tab \
  <sample>_1.fq.gz <sample>_2.fq.gz
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/kallisto_count_gtf.sh`




## Transcript abundance estimation (TPM/RPKM)

Gene-level TPM and RPKM values were computed from the featureCounts output using a custom R script `RPKM_TPM_featurecounts.R`. For each input featureCounts table, the script converts raw counts to TPM/RPKM using the feature length column (Length) and outputs two tab-separated files: `*_tpm.txt` and `*_rpkm.txt`.

The following command was run for each sample:
```bash
Rscript scripts/RPKM_TPM_featurecounts.R <featureCounts_output>.txt <output_prefix>
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/scripts/featurecounts_TPM.sh`



## Coverage track & QC

Genome-wide coverage tracks (.bigWig) were generated from coordinate-sorted BAM files using **bamCoverage** (v3.5.1). Coverage was normalized using BPM (bins per million mapped reads) with single–base resolution.

The following command was run for each sample:
```bash
bamCoverage \
  --bam <sample>_Aligned.sortedByCoord.out.bam \
  --outFileName <sample>.bw \
  --outFileFormat bigwig \
  --normalizeUsing BPM \
  --binSize 1
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/bigwigmerge_STAR_TPM_bw.sh `

For each biological condition, replicate bigWig files were aggregated by computing the median coverage signal across using **wiggletools** (v1.0.0). Median coverage tracks were first generated in bedGraph format, subsequently sorted using **bedtools** (v2.30.0), and finally converted back to bigWig format using **bedGraphToBigWig** (v4).

```bash
# Compute median coverage across replicates
wiggletools write_bg <condition>_median.bedGraph \
  median <rep1>.bw <rep2>.bw <rep3>.bw <rep4>.bw
# Sort bedGraph
bedtools sort -i <condition>_median.bedGraph > <condition>_median.sorted.bedGraph
# Convert bedGraph to bigWig
bedGraphToBigWig \
  <condition>_median.sorted.bedGraph \
  GRCh38_chrom_sizes.tab \
  <condition>_median.bw
```

--> Full script: `011_CristanchoLab/002_RNAseq_Preeti/scripts/bigwigmerge_STAR_TPM_bw.sh`






# Data analysis

## Gene expression (TPM)

Let's display transcript abundance in TPM for some genes of interest in **R** (v4.2.2):


```R
# Load packages
library("tidyverse")
library("rtracklayer")
library("ggpubr")

set.seed(42) # set seed for reproducibility

# Import gene annotation file to convert gene ID to gene symbol
gtf <- import("../../Master/meta/gencode.v47.annotation.gtf")
gene_table <- mcols(gtf) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct() %>%
  as_tibble()
colnames(gene_table) <- c("Geneid", "geneSymbol")


# Import TPM counts for each sample and merge them in a unique table
## Collect samples IDs
samples <- c("PSC_Norm_Rep1", "PSC_Norm_Rep2" ,"PSC_Norm_Rep3" ,"PSC_Norm_Rep4" ,"PSC_Hypo_Rep1", "PSC_Hypo_Rep2", "PSC_Hypo_Rep3", "PSC_Hypo_Rep4", "ReN_Norm_Rep1", "ReN_Norm_Rep2" ,"ReN_Norm_Rep3" ,"ReN_Norm_Rep4" ,"ReN_Hypo_Rep1", "ReN_Hypo_Rep2", "ReN_Hypo_Rep3", "ReN_Hypo_Rep4")

## Loop over samples and import TPM
sample_data <- list()
for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/tpm_featurecounts_multi/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    dplyr::select(Geneid, starts_with("output.STAR.")) %>%
    dplyr::rename(!!sample := starts_with("output.STAR."))
}

## Merge all samples into a single table and add gene symbol
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")  %>%
  left_join(gene_table)
###### --- Export --- ######
write.table(tpm_all_sample,  file = "output/tpm_featurecounts_multi/tpm_all_sample.tsv",  sep = "\t",  row.names = FALSE,  quote = FALSE)
###### -------------- ######


# Re-shape TPM counts table
## Shape to long format and extract sample metada from the column names
tpm_all_sample_tidy <- tpm_all_sample %>%
  pivot_longer(
    cols = -c(Geneid, geneSymbol),
    names_to = "variable",
    values_to = "tpm"
  ) %>%
  separate(
    variable,
    into = c("tissue", "condition", "replicate"),
    sep = "_"
  ) %>%
  dplyr::rename(
    gene = Geneid
  )

## Plot example: EPHB2 (log2(TPM+1)), Norm vs Hypo, for each tissue
my_comparisons <- list( c("Norm", "Hypo") )
pdf("output/tpm_featurecounts_multi/tpm-EPHB2.pdf", width=4, height=3)
tpm_all_sample_tidy %>%
  filter(geneSymbol %in% c("EPHB2") ) %>% 
  unique() %>%
  mutate(TPM = log2(tpm + 1) ) %>%
    ggboxplot(., x = "condition", y = "TPM",
                 fill = "condition",
                 palette = c("blue","red")) +
      # Add the statistical comparisons
      stat_compare_means(comparisons = my_comparisons, 
                        method = "t.test", 
                        aes(group = condition)) +
      theme_bw() +
      facet_wrap(~tissue) +
      ylab("log2(TPM + 1)") +
      ggtitle("EPHB2")
dev.off()


## Replicate-level TPM expression of hypoxia marker genes
pdf("output/tpm_featurecounts_multi/tpm-hypoxia_genes-dots_by_rep.pdf", width=5, height=15)
tpm_all_sample_tidy %>%
  filter(geneSymbol %in% c( "VEGFA",  "ADM",  "EGLN3",  "BNIP3",  "BNIP3L",  "CA9",  "CA12",  "LDHA",  "SLC2A1",  "ENO1",  "PFKP",  "HK2",  "ALDOA",  "PDK1")) %>%
  distinct(geneSymbol, tissue, condition, replicate, tpm, .keep_all = TRUE) %>%
  mutate(
    condition = factor(condition, levels = c("Norm", "Hypo")),
    replicate = factor(replicate),
    TPM = log2(tpm + 1)
  ) %>%
  ggplot(aes(x = condition, y = TPM, color = replicate)) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    size = 2.2,
    alpha = 0.9
  ) +
  theme_bw() +
  facet_grid(geneSymbol ~ tissue, scales = "free_y") +
  ylab("log2(TPM + 1)") +
  xlab("") +
  guides(color = guide_legend(title = "Bio Rep"))

dev.off()
```

**OUTPUT FILES**:
- Gene-level TPM expression matrix for all samples: `output/tpm_featurecounts_multi/tpm_all_sample.tsv`
- Gene-specific TPM boxplots: `output/tpm_featurecounts_multi/tpm-[GENE OF INTEREST]`
- Replicate-level TPM expression of hypoxia marker genes: `output/tpm_featurecounts_multi/tpm-hypoxia_genes-dots_by_rep.pdf`


**CONCLUSION**:

--> Due to ReN cells showing far less DEGs than PSC; we checked whether some replicates were more or less sensitive to the Hypoxia condition by checking hypoxia marker genes and label replicates: No replicate effect detected; all Hypoxia marker genes respond similarly between replicates.

Because ReN cells exhibited fewer differentially expressed genes than hiPSCs, replicate-level TPM expression of established hypoxia marker genes was examined using dot plots to assess potential variability in hypoxia sensitivity. No replicate-specific effects were observed: hypoxia markers showed consistent TPM expression patterns across replicates, indicating that the reduced number of DEGs in ReN cells is unlikely to be driven by replicate-level variability.



## Differential Gene Expression 



## Alternative splicing 



## Functional enrichment 





