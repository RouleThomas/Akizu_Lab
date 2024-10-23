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
sbatch scripts/fastp.sh # 28206738 xxx
```



## mapping fastp trim

```bash
sbatch --dependency=afterany:28206738 scripts/STAR_mapping_fastp.sh # 28206965 xxx
```

--> XXX Around 80-90% uniq aligned reads



# Count with featureCounts


XXXY WATCHE OUT DATA IS STRADNED!!!





Count on gene features with parameter
```bash
conda activate featurecounts

# all samples:
sbatch --dependency=afterany:13367330 scripts/featurecounts.sh # 13368068 bad wrong gtf
sbatch scripts/featurecounts.sh # 13504690 ok

sbatch scripts/featurecounts_stranded.sh # 13522388 ok

featureCounts -p -C -O -M --fraction \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf \
	-o output/test/50dN_WT_R1.txt output/STAR/fastp/50dN_WT_R1_Aligned.sortedByCoord.out.bam


```
test with `50dN_WT_R1`
--> Around 50% of succesfully assigned alignments with `-p -C -O` parameters...
--> data is stranded, so `-s 1` = 2% and `-s 2` = 58.1% ; same as the unstranded
--> counting multimapped reads with `-M --fraction` = 73%
----> Let's use the default method; as in `001_EZH1*/001__RNAseq`
------> So overall 50-75% alignment; that is OK


## Calculate TPM and RPKM


Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
conda activate deseq2
# Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_TPM.sh # 13526347 ok
# mv all output to output/tpm or rpkm folder
mv output/featurecounts/*tpm* output/tpm/
mv output/featurecounts/*rpkm* output/rpkm/
```

All good. 


# Shiny app


generate the `tpm_all_sample.txt` file and then go into `001_EZH1*/001__RNAseq` to create shiny app V2 including these

```R
library("tidyverse")


# Display some genes in TPM: 
# ---> The code below is not perfect; issue at the geneSymbol conversin; to troubleshoot later; but it work
#### Generate TPM for ALL samples
#### collect all samples ID
samples <- c( "100dN_WT_R1", "75dN_WT_R1", "50dN_WT_R1", "25dN_WT_R1", "NPC_WT_R1", "ESC_WT_R1","100dN_WT_R2", "75dN_WT_R2",  "50dN_WT_R2",    "25dN_WT_R2","NPC_WT_R2", "ESC_WT_R2", "100dN_WT_R3",  "75dN_WT_R3", "50dN_WT_R3",  "25dN_WT_R3",  "NPC_WT_R3",  "ESC_WT_R3")

## Make a loop for importing all tpm data and keep only ID and count column
sample_data <- list()

for (sample in samples) {
  sample_data[[sample]] <- read_delim(paste0("output/tpm/", sample, "_tpm.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    select(Geneid, starts_with("output.STAR.")) %>%
    rename(!!sample := starts_with("output.STAR."))
}

## Merge all dataframe into a single one
tpm_all_sample <- purrr::reduce(sample_data, full_join, by = "Geneid")
write.csv(tpm_all_sample, file="output/tpm/tpm_all_sample_Ciceri.txt")
### If need to import: tpm_all_sample <- read_csv("output/tpm/tpm_all_sample_Ciceri.txt") %>% dplyr::select(-"...1") #To import



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
write.table(long_data_log2tpm, file = c("output/tpm/long_data_log2tpm_Ciceri.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


```


--> Next here: `001_EZH1*/001__RNAseq` (`## Shiny app V2; including Ciceri RNAseq neuron diff dataset`)



# Generate Bigwig coverage files

Let's generate **TPM coverage**:

```bash
conda activate deeptools
# run time-per-time:
sbatch scripts/TPM_bw.sh # 15126648 ok
```


Let's merge the bigwig into 1 file with wiggletools (will do average of bigwig signal and not sum, many options see [github](https://github.com/Ensembl/WiggleTools)):


**Run wiggletools:**
```bash
conda activate BedToBigwig
sbatch --dependency=afterany:15126648 scripts/bigwigmerge_TPM.sh # 15127640 ok
```
*NOTE: bigwig are merge into 1 bedgraph which is then converted into 1 bigwig (wiggletools cannot output bigwig directly so need to pass by bedgraph or wiggle in between)*





