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



