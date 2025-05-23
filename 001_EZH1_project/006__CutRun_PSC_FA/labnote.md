# Project

1 condition, 3 genotypes, plenty of AB
- PSC (ESC); WT, KO and KOEF1aEZH1 (KO EZH1 with ef1a strong promoter EZH1-HA); AB: EZH1cs, EZH2, H3K27me3, SUZ12 and IGG


**Objectives:**
- Check EZH1 is working without overexpressing
- Compare level of H3K27me3 between WT and KO


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

Go [there](http://data-deliver.novogene.com/login/X202SC23092052-Z01-F001) and enter credetnial: (check email Novogen)

I created a `nano url.txt` with all link and used `wget -i url.txt` to download them all (1 link per raw); then `mv input_raw_Novogene/*fq.gz input` .

# Rename file

I created a tab separated file with current / new file names (keeping the .fq.gz sufix) and then:

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_map.txt
```

--> All good

# Fastp cleaning

```bash
sbatch scripts/fastp_1.sh # 6677035 ok
sbatch scripts/fastp_2.sh # 6677038 ok
sbatch scripts/fastp_3.sh # 6677039 ok
```


# FastQC

**Raw:**
```bash
sbatch scripts/fastqc_1.sh # 6677168 ok
sbatch scripts/fastqc_2.sh # 6677169 ok
sbatch scripts/fastqc_3.sh # 6677171 ok
```

--> all good

**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:6677035 scripts/fastqc_fastp_1.sh # 6677362 ok
sbatch --dependency=afterany:6677038 scripts/fastqc_fastp_2.sh # 6677368 ok
sbatch --dependency=afterany:6677039 scripts/fastqc_fastp_3.sh # 6677369 ok
```

--> all good


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch scripts/bowtie2_1.sh # 6678343 ok
sbatch scripts/bowtie2_2.sh # 6678369 ok
sbatch scripts/bowtie2_3.sh # 6678383 ok
```

--> Looks good


Mapping on E coli --> TO DO LATER! 

```bash
conda activate bowtie2

sbatch scripts/bowtie2_MG1655_1.sh # 15829778 ok
sbatch scripts/bowtie2_MG1655_2.sh # 15829779 ok
sbatch scripts/bowtie2_MG1655_3.sh # 15829786 ok

```

--> between 0.5 - 2% uniquely aligned reads (not a lot..; previously `005__CutRun` 10% (in `003__CutRun` was less than 1%) )





## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-6678343.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_6678343.txt

for file in slurm-6678369.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_6678369.txt

for file in slurm-6678383.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_6678383.txt
```

Add these values to `/home/roulet/001_EZH1_project/006__CutRun_PSC_FA/samples_006.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >60% input reads as been uniquely mapped to the genome (75% non uniq)



# Calculate histone content
**XX Not done; not necessary as only H3K27me3... XX**


# Samtools and read filtering

--> See `METHOD GOOD TO FOLLOW` in `003__CutRun` labnote

## Marking dupplicates
```bash
conda activate bowtie2

sbatch scripts/samtools_1.sh # 6712195 ok
sbatch scripts/samtools_2.sh # 6712214 ok
sbatch scripts/samtools_3.sh # 6712231 ok
```

## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.

```bash
conda activate bowtie2

sbatch --dependency=afterany:6712195 scripts/samtools_unique_1.sh # 6712296 ok
sbatch --dependency=afterany:6712214 scripts/samtools_unique_2.sh # 6712305 ok
sbatch --dependency=afterany:6712231 scripts/samtools_unique_3.sh # 6712321 ok
```

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch --dependency=afterany:15829778 scripts/samtools_MG1655_unique_1.sh # 15829889 ok
sbatch --dependency=afterany:15829779 scripts/samtools_MG1655_unique_2.sh # 15829890 ok
sbatch --dependency=afterany:15829786 scripts/samtools_MG1655_unique_3.sh # 15829891 ok

# count the nb of unique reads
sbatch scripts/samtools_MG1655_unique_count.sh # 29772149 xxx

```

--> More information on this step in the `005__CutRun` labnote


While working on `016__integration_PSC_V1`, noticed that bam files are corrupted; transformed into `*.pgbam`; let's re-run samtools_unique for all files:


```bash
conda activate bowtie2

sbatch scripts/samtools_unique_correct.sh # 29565424 ok
```

--> all good.



# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:6712296 scripts/bamtobigwig_unique_1.sh # 6712485 ok
sbatch --dependency=afterany:6712305 scripts/bamtobigwig_unique_2.sh # 6712486 ok
sbatch --dependency=afterany:6712321 scripts/bamtobigwig_unique_3.sh # 6712487 ok

sbatch scripts/bamtobigwig_unique_1_missing.sh # 6731623 ok

# raw non unique bigwig
sbatch scripts/bamtobigwig_1.sh # 13165557 ok
sbatch scripts/bamtobigwig_2.sh # 13165558 ok
sbatch scripts/bamtobigwig_3.sh # 13165559 ok
```



- KOEF1aEZH1
*Pass*: SUZ12, EZH1cs, EZH2, H3K27me3
*Failed*: HA
- WT
*Pass*: SUZ12, EZH2, H3K27me3
*Failed*: HA, EZH1cs (better than earlier, more peaks, but still poorly)
- KO
*Pass*: H3K27me3
*Failed*: HA, SUZ12, EZH1cs, EZH2





## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch --dependency=afterany:6731623 scripts/multiBigwigSummary_all.sh # 6731623


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KOEF1aEZH1_SUZ12 PSC_KOEF1aEZH1_EZH2 PSC_KOEF1aEZH1_HA PSC_KOEF1aEZH1_EZH1cs PSC_KOEF1aEZH1_H3K27me3 PSC_KOEF1aEZH1_IGG PSC_WT_SUZ12 PSC_WT_EZH2 PSC_WT_HA PSC_WT_EZH1cs PSC_WT_H3K27me3 PSC_WT_IGG PSC_KO_SUZ12 PSC_KO_EZH2 PSC_KO_HA PSC_KO_EZH1cs PSC_KO_H3K27me3 PSC_KO_IGG \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf


## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KOEF1aEZH1_SUZ12 PSC_KOEF1aEZH1_EZH2 PSC_KOEF1aEZH1_HA PSC_KOEF1aEZH1_EZH1cs PSC_KOEF1aEZH1_H3K27me3 PSC_KOEF1aEZH1_IGG PSC_WT_SUZ12 PSC_WT_EZH2 PSC_WT_HA PSC_WT_EZH1cs PSC_WT_H3K27me3 PSC_WT_IGG PSC_KO_SUZ12 PSC_KO_EZH2 PSC_KO_HA PSC_KO_EZH1cs PSC_KO_H3K27me3 PSC_KO_IGG \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf


```






# MACS2 peak calling on bam unique

--> IGG samples used as control

--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_raw_1.sh # 6749578 ok
sbatch scripts/macs2_raw_2.sh # 6749582 ok
sbatch scripts/macs2_raw_3.sh # 6749583 ok
```



Then keep only the significant peaks (re-run the script to test different qvalue cutoff) and remove peaks overlapping with blacklist regions. MACS2 column9 output is -log10(qvalue) format so if we want 0.05; 
- q0.05: `q value = -log10(0.05) = 1.30103`
- q0.01 = 2
- q0.005 = 2.30103
- q0.001 = 3
- q0.0001 = 4
- q0.00001 = 5

```bash
conda activate bowtie2 # for bedtools
sbatch scripts/macs2_raw_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive

# quick command to print median size of peak within a bed
awk '{print $3-$2}' your_bed_file.bed | sort -n | awk 'BEGIN {c=0; sum=0;} {a[c++]=$1; sum+=$1;} END {if (c%2) print a[int(c/2)]; else print (a[c/2-1]+a[c/2])/2;}'
```

**Optimal qvalue** according to IGV:
- PSC_KOEF1aEZH1_SUZ12: 1.30103 (2.3 more true peak)
- PSC_KOEF1aEZH1_EZH2: 1.30103
- PSC_KOEF1aEZH1_EZH1cs: 1.30103
- PSC_KOEF1aEZH1_H3K27me3: 3
- PSC_WT_SUZ12: 1.30103 (2.3 more true peak)
- PSC_WT_EZH2: 1.30103
- PSC_WT_EZH1cs: FAIL
- PSC_WT_H3K27me3: 1.30103 (many true peak surprinsgly!)
- PSC_KO_SUZ12: FAIL
- PSC_KO_EZH2: FAIL
- PSC_KO_EZH1cs: FAIL
- PSC_KO_H3K27me3: 1.30103 (many true peak surprinsgly!)



## Assign peak to genes for NPC and PSC - raw macs2


```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")


# Import macs2 peaks
## KOEF
SUZ12 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_SUZ12_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
EZH2 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_EZH2_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
EZH1 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_EZH1cs_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3 = as_tibble(read.table('output/macs2/broad_blacklist_qval3/PSC_KOEF1aEZH1_H3K27me3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 


## Tidy peaks #-->> Re-Run from here with different qvalue!!
SUZ12_gr = makeGRangesFromDataFrame(SUZ12,keep.extra.columns=TRUE)
EZH2_gr = makeGRangesFromDataFrame(EZH2,keep.extra.columns=TRUE)
EZH1_gr = makeGRangesFromDataFrame(EZH1,keep.extra.columns=TRUE)
H3K27me3_gr = makeGRangesFromDataFrame(H3K27me3,keep.extra.columns=TRUE)
gr_list <- list(SUZ12=SUZ12_gr, EZH2=EZH2_gr, EZH1=EZH1_gr,  H3K27me3=H3K27me3_gr)
## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
SUZ12_annot <- as.data.frame(peakAnnoList[["SUZ12"]]@anno)
EZH2_annot <- as.data.frame(peakAnnoList[["EZH2"]]@anno)
EZH1_annot <- as.data.frame(peakAnnoList[["EZH1"]]@anno)
H3K27me3_annot <- as.data.frame(peakAnnoList[["H3K27me3"]]@anno)
## Convert entrez gene IDs to gene symbols
SUZ12_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = SUZ12_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
SUZ12_annot$gene <- mapIds(org.Hs.eg.db, keys = SUZ12_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
EZH2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH2_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
EZH1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH1_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(SUZ12_annot, file="output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_SUZ12_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(EZH2_annot, file="output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_EZH2_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(EZH1_annot, file="output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_annot, file="output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_H3K27me3_qval3.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
SUZ12_annot_promoterAnd5 = tibble(SUZ12_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
EZH2_annot_promoterAnd5 = tibble(EZH2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
EZH1_annot_promoterAnd5 = tibble(EZH1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_annot_promoterAnd5 = tibble(H3K27me3_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))    

### Save output gene lists
SUZ12_annot_promoterAnd5_geneSymbol = SUZ12_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
EZH2_annot_promoterAnd5_geneSymbol = EZH2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
EZH1_annot_promoterAnd5_geneSymbol = EZH1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_annot_promoterAnd5_geneSymbol = H3K27me3_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()   

write.table(SUZ12_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_SUZ12_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(EZH2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(EZH1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_H3K27me3_qval3_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    




## WT  

SUZ12 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_WT_SUZ12_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
EZH2 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_WT_EZH2_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_WT_H3K27me3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

SUZ12 = as_tibble(read.table('output/macs2/broad_blacklist_qval2.30103/PSC_WT_SUZ12_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
EZH2 = as_tibble(read.table('output/macs2/broad_blacklist_qval2.30103/PSC_WT_EZH2_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3 = as_tibble(read.table('output/macs2/broad_blacklist_qval2.30103/PSC_WT_H3K27me3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

## Tidy peaks #-->> Re-Run from here with different qvalue!!
SUZ12_gr = makeGRangesFromDataFrame(SUZ12,keep.extra.columns=TRUE)
EZH2_gr = makeGRangesFromDataFrame(EZH2,keep.extra.columns=TRUE)
H3K27me3_gr = makeGRangesFromDataFrame(H3K27me3,keep.extra.columns=TRUE)
gr_list <- list(SUZ12=SUZ12_gr, EZH2=EZH2_gr, H3K27me3=H3K27me3_gr)
## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
SUZ12_annot <- as.data.frame(peakAnnoList[["SUZ12"]]@anno)
EZH2_annot <- as.data.frame(peakAnnoList[["EZH2"]]@anno)
H3K27me3_annot <- as.data.frame(peakAnnoList[["H3K27me3"]]@anno)
## Convert entrez gene IDs to gene symbols
SUZ12_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = SUZ12_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
SUZ12_annot$gene <- mapIds(org.Hs.eg.db, keys = SUZ12_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
EZH2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH2_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(SUZ12_annot, file="output/ChIPseeker/annotation_macs2_PSC_WT_SUZ12_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(EZH2_annot, file="output/ChIPseeker/annotation_macs2_PSC_WT_EZH2_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_annot, file="output/ChIPseeker/annotation_macs2_PSC_WT_H3K27me3_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
SUZ12_annot_promoterAnd5 = tibble(SUZ12_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
EZH2_annot_promoterAnd5 = tibble(EZH2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_annot_promoterAnd5 = tibble(H3K27me3_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))    

### Save output gene lists
SUZ12_annot_promoterAnd5_geneSymbol = SUZ12_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
EZH2_annot_promoterAnd5_geneSymbol = EZH2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_annot_promoterAnd5_geneSymbol = H3K27me3_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()   

write.table(SUZ12_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_PSC_WT_SUZ12_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(EZH2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_PSC_WT_EZH2_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_PSC_WT_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    



## KO



H3K27me3 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_KO_H3K27me3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3 = as_tibble(read.table('output/macs2/broad_blacklist_qval2.30103/PSC_KO_H3K27me3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

## Tidy peaks #-->> Re-Run from here with different qvalue!!

H3K27me3_gr = makeGRangesFromDataFrame(H3K27me3,keep.extra.columns=TRUE)
gr_list <- list(H3K27me3=H3K27me3_gr)
## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame

H3K27me3_annot <- as.data.frame(peakAnnoList[["H3K27me3"]]@anno)
## Convert entrez gene IDs to gene symbols

H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table

write.table(H3K27me3_annot, file="output/ChIPseeker/annotation_macs2_PSC_KO_H3K27me3_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!

H3K27me3_annot_promoterAnd5 = tibble(H3K27me3_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))    

### Save output gene lists

H3K27me3_annot_promoterAnd5_geneSymbol = H3K27me3_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()   


write.table(H3K27me3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_PSC_KO_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    

   
```


# Deeptools plot in KOEF1aEZH1

Check co-localization of EZH1cs, EZH2, SUZ12, H3K27me3 in KOEF1aEZH1 using the unique bigiwg raw


```bash
# all genes
sbatch scripts/matrix_TSS_10kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_raw_allGenes.sh # 16290946 ok
sbatch scripts/matrix_TSS_5kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_raw_allGenes.sh # 16290958 ok
sbatch scripts/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_raw_allGenes.sh # 16290920 ok

# macs2 EZH1cs peaks
sbatch scripts/matrix_TSS_10kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_raw_EZH1csPeaks_macs2q1.30103.sh # 16291108 ok
sbatch scripts/matrix_TSS_5kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_raw_EZH1cspeaks_macs2q1.30103.sh # 16291163 ok
sbatch scripts/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_raw_EZH1cspeaks_macs2q1.30103.sh # 16291207 ok

# macs2 EZH2 peaks
sbatch scripts/matrix_TSS_10kb_EZH2_EZH1cs_SUZ12_H3K27me3_IGG_raw_EZH2Peaks_macs2q1.30103.sh # 16346993 ok
sbatch scripts/matrix_TSS_5kb_EZH2_EZH1cs_SUZ12_H3K27me3_IGG_raw_EZH2Peaks_macs2q1.30103.sh # 16347015 ok
sbatch scripts/matrix_TSS_2kb_EZH2_EZH1cs_SUZ12_H3K27me3_IGG_raw_EZH2Peaks_macs2q1.30103.sh # 16347040 ok

# macs2 SUZ12 peaks
sbatch scripts/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103.sh # 17892267 ok


# Comparison KOEF1aEZH1 with WT (from 005__CutRun)
## all genes
sbatch scripts/matrix_TSS_10kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_allGenes.sh # 16347385 ok
sbatch scripts/matrix_TSS_5kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_allGenes.sh # 16347404 ok
sbatch scripts/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_allGenes.sh # 16347417 ok


## EZH1 peaks
sbatch scripts/matrix_TSS_10kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.sh # 16347554 ok
sbatch scripts/matrix_TSS_5kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.sh # 16347591 fail; 16380981 ok
sbatch scripts/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.sh # 16347612 fail; 16381048 ok


# Comparison KOEF1aEZH1 with WT (from 006__CutRun this CutRun) _ THOR (TMM; no spike in)
sbatch scripts/matrix_TSS_5kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.sh # 18177870 ok
sbatch scripts/matrix_TSS_5kb_EZH1cs_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.sh # 19041320 ok
sbatch scripts/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.sh # 19041369 ok
sbatch scripts/matrix_TSS_5kb_SUZ12_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.sh # 19041458 ok
sbatch scripts/matrix_TSS_5kb_H3K27me3_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.sh # 19041477 ok


# Comparison KOEF1aEZH1 with WT (from 006__CutRun this CutRun) _ THOR (spike in TMM for EZH2, SUZ12, H3K27me3; for EZH1; TMM solely)
sbatch scripts/matrix_TSS_5kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT006_THORandTHORTMM_EZH1csPeaks_macs2q1.30103.sh # 18255674 ok


```

*NOTE: EZH1cs and EZH2 optimal MACS2 qval 1.30103 with 1,271 peaks*

--> Clear co-localization of EZH1, EZH2, SUZ12, and H3K27me3 ; especially when looking at the EZH1cs or EZH2 peaks

--> KOEF1aEZH1 and WT look the same for 



# Deeptools plot for comparison EZH1 EZH2 co-localization

Let's check whether EZH1 and EZH2 co-localize; identify EZH1 and EZH2 solely bound regions and check H3K27me3 levels.

- Venn diagram online of the genes assigned with EZH1 and EZH2 peaks (file named `Venn_overlap_EZH1EZH2__*.txt`)
- Exctract and import to cluster these list of genes
- Generate new gtf (EZH1; EZH1-EZH2; EZH2)
- Check EZH1, EZH2, and H3K27me3 deeptool profile in these 3 groups of genes





```bash

# Generate gtf from gene Symbol list
perl -p -i -e 's/\r$//' output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH1only.txt  # THIS TO CONVERT windowns to UNIX; as .txt from windows...
perl -p -i -e 's/\r$//' output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH2only.txt  # THIS TO CONVERT windowns to UNIX; as .txt from windows...
perl -p -i -e 's/\r$//' output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH1andEZH2.txt  # THIS TO CONVERT windowns to UNIX; as .txt from windows...

### Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH1only.txt > output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH1only_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH2only.txt > output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH2only_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH1andEZH2.txt > output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH1andEZH2_as_gtf_geneSymbol.txt


### Filter the gtf
grep -Ff output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH1only_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI-EZH1EZH2__EZH1only.gtf
grep -Ff output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH2only_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI-EZH1EZH2__EZH2only.gtf
grep -Ff output/ChIPseeker/Venn_overlap_EZH1EZH2__EZH1andEZH2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI-EZH1EZH2__EZH1andEZH2.gtf

```


Now deepTools

```bash
conda activate deeptools

sbatch scripts/matrix_TSS_10kb_bigwig_unique-EZH1EZH2_PSC.sh # 7075199 ok


```

--> EZH1 and EZH2 seems to always co-localize.







# Deeptools plot for comparison EZH1 CutRun in native (`005`) and FA (`006`) conditions in PSC KOEF1aEZH1

## Gene comparison
- Collect all genes bound with EZH1cs (=860 genes); according to macs2 qval 1.3 (`output/ChIPseeker/annot_macs2_KOEF1aEZH1_EZH1cs_qval1.30103_promoterAnd5_geneSymbol.txt`) in 005 and 006 = `annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006.txt`
- Generate gtf of these genes
- deeptool plots with bigwig EZH1cs `005` and `006`



```bash

# Generate gtf from gene Symbol list
perl -p -i -e 's/\r$//' output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006.txt  # THIS TO CONVERT windowns to UNIX; as .txt from windows...


### Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006.txt > output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006_as_gtf_geneSymbol.txt

### Filter the gtf
grep -Ff output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI-annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006.gtf


```


Now deepTools with *bigiwg_unique* files

```bash
conda activate deeptools

# all 860 genes
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006.sh # 11447101 ok


```

--> **FA for EZH1cs CutRun in PSC works better than native**; the signal to noise ratio is higher


## Peak comparison

- Collect native peaks and check EZH1 profile for native and FA
- Same starting with FA peaks


```bash
# peak in 005 native
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103-CutRun005.sh # 11455500 ok

# peak in 006 FA
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103-CutRun006.sh # 11455846 ok

```



# Deeptools plot for comparison EZH1 CutRun in WT vs KOEF1EZH1


Is the signal in WT, true or not? --> Collect true EZH1cs signal in PSC_KOEF1EZH1 in FA condition (better than native) and check signal in WT
- Overlap= signal is true; play with IGG to clean it
- No overlap= what we see is mostly ‘noise’; no need computational trick


Now deepTools with *bigiwg_unique* files

```bash
conda activate deeptools

# all EZH1cs peaks in KOEF1aEZH1
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103-CutRun005and006_KOEF1aEZH1andWT.sh # 11457719 ok

# all 860 genes bound with EZH1cs
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006_KOEF1aEZH1andWT.sh # 11457352 ok



```


--> There is overlap so the very low EZH1cs signal in WT is true signal as it overlap with KOEF1EZH1 signal.
----> There might be computational tricks to clean out low and noisy EZH1cs signal in the WT...?





# Deeptools plot for comparison EZH1 CutRun in WT vs KOEF1EZH1


Is the signal in WT, true or not? --> Collect true EZH1cs signal in PSC_KOEF1EZH1 in FA condition (better than native) and check signal in WT
- Overlap= signal is true; play with IGG to clean it
- No overlap= what we see is mostly ‘noise’; no need computational trick


Now deepTools with *bigiwg_unique* files

```bash
conda activate deeptools

# all EZH1cs peaks in KOEF1aEZH1
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103-CutRun005and006_KOEF1aEZH1andWT.sh # 11457719 ok

# all 860 genes bound with EZH1cs
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1_EZH1_qval1.30103_promoterAnd5_geneSymbol-CutRun005and006_KOEF1aEZH1andWT.sh # 11457352 ok



```


--> There is overlap so the very low EZH1cs signal in WT is true signal as it overlap with KOEF1EZH1 signal.
----> There might be computational tricks to clean out low and noisy EZH1cs signal in the WT...?





# Deeptools plot for comparison EZH2 CutRun in WT vs KOEF1EZH1

## Gene comparison

To check EZH2 AB specificity? Not clear to me:

*“High overlap could suggest that EZH2 recognize EZH1”*… Why?? Because many EZH1 in KOEF1EZH1 so EZH2 could recognize them? High overlap could also suggest EZH2 binding is not affected by EZH1 overexpression, right?

- Collect all 1,437 genes that are bound with EZH2 in either WT and KOEF1aEZH1 and check EZH2 profile (file = `output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1andWT_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt`)
- generate gtf of these genes


```bash

# Generate gtf from gene Symbol list
perl -p -i -e 's/\r$//' output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1andWT_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt  # THIS TO CONVERT windowns to UNIX; as .txt from windows...


### Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1andWT_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1andWT_EZH2_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt

### Filter the gtf
grep -Ff output/ChIPseeker/annotation_macs2_PSC_KOEF1aEZH1andWT_EZH2_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI-annotation_macs2_PSC_KOEF1aEZH1andWT_EZH2_qval1.30103_promoterAnd5_geneSymbol.gtf


```




Now deepTools with *bigiwg_unique* files

```bash
conda activate deeptools

# all 1,437 genes bound with EZH2
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1andWT_EZH2_qval1.30103_promoterAnd5_geneSymbol.sh # 11460412 ok

```


--> very high overlap!



## Peak comparison

- Collect WT peaks and check EZH2 profile for WT and KOEF1EZH1
- Same starting with KOEF1EZH1 peaks



```bash
# peak in WT
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_WT_EZH2_qval1.30103-WTandKOEF1aEZH1.sh # 11460783 ok

# peak in KOEF1aEZH1
sbatch scripts/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_KOEF1aEZH1_EZH2_qval1.30103-WTandKOEF1aEZH1.sh # 11460843 ok
```

--> very high overlap!





# Spike in factor

Let's do the analysis for H3K27me3 only; **compare WT vs KOEF1aEZH1 and WT vs KO**

Test 2 spikein normalization method (histone and Ecoli)


## Calculate histone content

--> This histone content will be used to generate a scaling factor which will be used to histone-scaled our library size. The calcul/method to follow is from `003__CutRun/output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt`

**Pipeline:**
- Count the histone barcode on the clean reads
- Calculate SF (group by sample (replicate) and AB and calculate the total nb of reads. Then proportion of reads = nb read in sample / total reads. SF = min(proportion) / sample proportion)

--> Paste metric values in `samples_006.xlsx` in Google Drive

## Count the histone barcode on the clean reads



```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp.sh # 15836536 ok

```


--> It output the nb of reads found for each histone; then simply copy paste to the excell file `output/spikein/SpikeIn_QC_fastp_006.xlsx` in GoogleDrive

- `PSC_WT_H3K27me3`: enriched in H3K27me3
- `PSC_KO_H3K27me3`: enriched in H3K27me3
- `PSC_KOEF1aEZH1_H3K27me3`: enriched in H3K27me3



## histone spike in factor


```R
# package
library("tidyverse")
library("readxl")
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_006.xlsx") 

## H3K27me3 
spikein_H3K27me3 = spikein %>%
    filter(Target == "H3K27me3",
    sample_ID %in% c("PSC_KO_H3K27me3", "PSC_KOEF1aEZH1_H3K27me3","PSC_WT_H3K27me3")) %>%
    group_by(sample_ID, AB) %>%
    summarise(aligned=sum(counts))
# Total reads per IP
spikein_H3K27me3_total = spikein_H3K27me3 %>%
    ungroup() %>%
    group_by(AB) %>%
    mutate(total = sum(aligned)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_H3K27me3_read_prop = spikein_H3K27me3 %>%
    left_join(spikein_H3K27me3_total) %>%
    mutate(read_prop = aligned / total)
spikein_H3K27me3_read_prop_min = spikein_H3K27me3_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_H3K27me3_scaling_factor = spikein_H3K27me3_read_prop %>%
    left_join(spikein_H3K27me3_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_H3K27me3_scaling_factor, file="output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt", sep="\t", quote=FALSE, row.names=FALSE)





```

--> All good, histone SF closely similar to the MG1655 ones




### Quality control plot

Then look at the xlsx file from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) to generate quality control plot. Use R cluster for vizualization (file is `spikein_QC.xlsx` in Google Drive), file in `output/spikein`.
```R
# package
library("tidyverse")
library("readxl")
# import df adn tidy to remove AB used in sample_ID
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_006.xlsx") %>%
  separate(sample_ID, into = c("type", "condition", "tag"), sep = "_") %>%
  mutate(sample_ID = paste(type, condition, sep = "_")) %>%
  select(-type, -condition, -tag, -tissue)


# NPC
# data processing
spikein_sum_Barcode_read <- spikein %>%
	select(-Barcode) %>%
  	group_by(sample_ID, Target, AB) %>% 
	summarise(sum_read=sum(counts)) %>%
	unique() 

spikein_sum_Total_read <- spikein %>%
	select(-Barcode, -Target) %>%
  	group_by(sample_ID, AB) %>% 
	summarise(total_read=sum(counts)) %>%
	unique() 	

spikein_all <- spikein_sum_Barcode_read %>%
	left_join(spikein_sum_Total_read) %>%
	mutate(target_norm = (sum_read/total_read)*100)

## Histone scaling for H3K27me3
spikein_all_scale = spikein_all %>%
  group_by(sample_ID) %>%
  # Find the target_norm value when Target is H3K27me3 and AB is H3K27me3
  mutate(scaling_factor = ifelse(Target == "H3K27me3" & AB == "H3K27me3", target_norm, NA)) %>%
  # Fill the scaling_factor column with the appropriate value within each group
  fill(scaling_factor, .direction = "downup") %>%
  # Scale the target_norm values
  mutate(scaled_target_norm = target_norm / scaling_factor * 100) %>%
  # Remove the scaling_factor column
  select(-scaling_factor) %>%
  # Ungroup the data
  ungroup()
# Plot
pdf("output/spikein/QC_histone_spike_in_H3K27me3.pdf", width = 10, height = 4)
spikein_all_scale %>%
    filter(
           AB %in% c("H3K27me3", "IGG")) %>%
        ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
        geom_col(position = "dodge") +
        facet_wrap(~sample_ID, nrow=1) +
        geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()




```


--> All good H3K27me3 enriched 



# Ecoli scaling factor (copy from `008__CutRun`)
## Mapping E coli

**Do it for H3K27me3 for WT, KOEF1aEZH1, KO**

- Map our reads to the E. coli genome using same parameters as for human.
- Count the number of aligned reads to the spike-in control sequences for each sample `samtools view -S -F 4 -c sample.sam > sample_spikein_count.txt`
- Do the math for scaling factor, same method as when using histone spike-in

```bash
# count nb of reads aligned to genome

samtools view -S -F 4 -c output/spikein/PSC_KOEF1aEZH1_H3K27me3_MG1655.sam > output/spikein/PSC_KOEF1aEZH1_H3K27me3-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KOEF1aEZH1_HA_MG1655.sam > output/spikein/PSC_KOEF1aEZH1_HA-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KOEF1aEZH1_EZH1cs_MG1655.sam > output/spikein/PSC_KOEF1aEZH1_EZH1cs-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KOEF1aEZH1_EZH2_MG1655.sam > output/spikein/PSC_KOEF1aEZH1_EZH2-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KOEF1aEZH1_SUZ12_MG1655.sam > output/spikein/PSC_KOEF1aEZH1_SUZ12-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KOEF1aEZH1_IGG_MG1655.sam > output/spikein/PSC_KOEF1aEZH1_IGG-spikein_count.txt

samtools view -S -F 4 -c output/spikein/PSC_KO_H3K27me3_MG1655.sam > output/spikein/PSC_KO_H3K27me3-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KO_HA_MG1655.sam > output/spikein/PSC_KO_HA-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KO_EZH1cs_MG1655.sam > output/spikein/PSC_KO_EZH1cs-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KO_EZH2_MG1655.sam > output/spikein/PSC_KO_EZH2-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KO_SUZ12_MG1655.sam > output/spikein/PSC_KO_SUZ12-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_KO_IGG_MG1655.sam > output/spikein/PSC_KO_IGG-spikein_count.txt

samtools view -S -F 4 -c output/spikein/PSC_WT_H3K27me3_MG1655.sam > output/spikein/PSC_WT_H3K27me3-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_WT_HA_MG1655.sam > output/spikein/PSC_WT_HA-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_WT_EZH1cs_MG1655.sam > output/spikein/PSC_WT_EZH1cs-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_WT_EZH2_MG1655.sam > output/spikein/PSC_WT_EZH2-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_WT_SUZ12_MG1655.sam > output/spikein/PSC_WT_SUZ12-spikein_count.txt
samtools view -S -F 4 -c output/spikein/PSC_WT_IGG_MG1655.sam > output/spikein/PSC_WT_IGG-spikein_count.txt
```


--> There is some uniq mapped reads, around xxx% (in `003__CutRun` was less than 1%)

Now calculate SF in R, as for histone SF:


```R
# package
library("tidyverse")
library("readxl")
library("ggpubr")

# SF H3K27me3
spikein <- read_excel("output/spikein/SpikeIn_MG1655_006.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "H3K27me3")
# Total reads per IP
spikein_H3K27me3_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_H3K27me3_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_H3K27me3_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)

# SF EZH1cs
spikein <- read_excel("output/spikein/SpikeIn_MG1655_006.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "EZH1cs")
# Total reads per IP
spikein_EZH1cs_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_EZH1cs_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_EZH1cs_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)



# SF EZH2
spikein <- read_excel("output/spikein/SpikeIn_MG1655_006.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "EZH2")
# Total reads per IP
spikein_EZH2_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_EZH2_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_EZH2_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)


# SF HA
spikein <- read_excel("output/spikein/SpikeIn_MG1655_006.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "HA")
# Total reads per IP
spikein_HA_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_HA_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_HA_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)

# SF SUZ12
spikein <- read_excel("output/spikein/SpikeIn_MG1655_006.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "SUZ12")
# Total reads per IP
spikein_SUZ12_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_SUZ12_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_SUZ12_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)


```

-->  histone vs MG1655 SF; same direction
----> GOOD!!




# Spike in scaling
## With MG1655 spike in for PTM CutRun


--> Let's use MG1655 as default method fopr spike in normalization! Seems more accurate and can be used for all AB!

**Using our scaling factor, let's estimate the 'new' library size** and provide it to `dba.normalize(library = c(1000, 12000))` = Like that our library size will be change taking into account our scaling factor! **Then we can normalize with library-size, RLE or TMM**... (issue discussed [here](https://support.bioconductor.org/p/9147040/)) 


### Adjust library size with MG1655 scaling factor and apply normalization
Total number of reads is our library size (used samtools flagstat to double check) :

`samtools flagstat output/bowtie2/*.dupmark.sorted.bam` used to obtain library size (first value=library size)
--> Values save in GoogleDrive `006__*/samples_006.xlsx`. Histone-norm-library-size = library-size * SF. Using the non-reciprocal scaling factor, we increase the library-size; the more histone enriched, the more library size is increased, thus the more signal will decrease.

Now let's use these new histone-scaled library size and normalize with library-size,TMM or RLE. Let's use the **unique bam files** together with the **unique bam MACS2 raw files (xlsx, not the bed with pre-filtered qvalue)**

***Key points:***
- **Let's do 1 DiffBind per AB (H3K27me3, H3K4me3,...) and tissue (PSC, NPC); otherwise the TMM normalization may take all, unrelated, samples into account!** --> Files are `meta_sample_macs2raw_unique*.txt`
- **For the non-histone CutRun, I will use the library size non histone scaled in DiffBind to collect TMM normalized SF**; I tested with and without specifying library size; and it does not change a lot the SF... Let's better use the one RiP method w/o providing the library size! Should provide BETTER correction


**--> MG1655 SF used:**

```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 

# ONE PER ONE
## PSC_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_H3K27me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_H3K27me3.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(37338241,7794352,38565592), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC_H3K27me3.txt")



## PSC_EZH1
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_EZH1cs.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## --> Fail cannot as no peak


## PSC_EZH1 WT vs KOEF1aEZH1
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_EZH1cs_WTvsKOEF1aEZH1.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## --> Fail cannot as no peak



## PSC_EZH2 WT vs KOEF1aEZH1
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_EZH2_WTvsKOEF1aEZH1.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_EZH2_WTvsKOEF1aEZH1.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_EZH2_WTvsKOEF1aEZH1.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_EZH2_WTvsKOEF1aEZH1.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_EZH2_WTvsKOEF1aEZH1.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(22168042,34263744), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC_EZH2_WTvsKOEF1aEZH1.txt")




## PSC_SUZ12 WT vs KOEF1aEZH1
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_SUZ12_WTvsKOEF1aEZH1.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_SUZ12_WTvsKOEF1aEZH1.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_SUZ12_WTvsKOEF1aEZH1.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_SUZ12_WTvsKOEF1aEZH1.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_SUZ12_WTvsKOEF1aEZH1.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(56199091,108075686), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC_SUZ12_WTvsKOEF1aEZH1.txt")

```

--> WT and KOEFaEZH1 cluster more together, KO is apart

--> EZH1 DiffBind TMM norm did not work as no macs2 peak are available...



# THOR

Let's use THOR, notably to have IGG scaled bigwig...!

Comparison to do; PSC WT vs KO and WT vs KOEF1aEZH1:
- H3K27me3



--> SF to use in THOR are the **reciprocal of MG1655_DiffBind_TMM**
--> Configs file created manually as `output/THOR/PSC_WTvsKO_H3K27me3.config`

--> Lets also try to use the DiffBind spike in BAM method (similarly use the reciprocal from diffBind)

--> Lets also try to use **regular TMM normalization method from THOR** (no Ecoli spike in scaling factor)


*THOR is very buggy to make it work I need to temporaly change where to look for libraries lol.. So cannot use nano anymore for example...*

*Follow these parameters: `WTvsHET_unique_Keepdup` (perform best in previous CutRun)*

```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge

# E coli DiffBind TMM scaling factor
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3.sh # 15945828 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3.sh # 15945937 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2.sh # 18253957 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12.sh # 18253958 ok

# Default THOR TMM normalization (no E coli spike in norm)
sbatch scripts/THOR_PSC_WTvsKO_H3K27me3_TMM.sh # 18083826 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_TMM.sh # 18084031 ok
sbatch scripts/THOR_PSC_WTvsKO_EZH1cs_TMM.sh # 18084234 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH1cs_TMM.sh # 18084235 ok
sbatch scripts/THOR_PSC_WTvsKO_EZH2_TMM.sh # 18084400; 18177095 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_EZH2_TMM.sh # 18084676 ok
sbatch scripts/THOR_PSC_WTvsKO_SUZ12_TMM.sh # 18084680 ok
sbatch scripts/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_TMM.sh # 18084683 ok
```

--> the KO samples for EZH2, EZH1?, SUZ12 are too shitty so the WT vs KO comparison CANNOT be used.

--> WT vs KOEF1aEZH1 are good to use





## Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks


```R

# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# H3K27me3 WTvsKO
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3/PSCWTvsKOH3K27me3-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  mutate(FC = (count_KO) / (count_WT))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_PSC_WTvsKO_H3K27me3/log2FC_qval30.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KO_qval30") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 100) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval100.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 30) %>%
  group_by(X6) %>%
  summarise(n = n())



# H3K27me3 WTvsKOEF1aEZH1
diffpeaks <- read_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/PSCWTvsKOEF1aEZH1H3K27me3-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1) / (count_WT))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/log2FC_qval30.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WT vs KOEF1aEZH1_qval30") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3/THOR_qval50.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())

```

- *NOTE: FC positive = less in KO; negative = more in KO*

**Optimal qvalue:**
--> *H3K27me3*; qval 30 (qval50 ok also,more stringeant)





## deeptools plot THOR bigwig
- all genes with peak in WT and or KO (macs2 peaks; prefered used qval2.3 as in `CutRun__009`)
- all diff. peaks WT vs KO (THOR; preferred used qval30 as in `CutRun__009`) 

Let's see WT vs KO in PSC (in neurons (`003__CutRun` and `007__CutRun`) and NPC (`006__CutRun` and `008__CutRun`) we have more H3K27me3 in KO)

And check WT vs KOEF1aEZH1; would expect increased H3K27me3 too?


```bash

# allGenes
sbatch scripts/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_allGenes.sh # 15955393 ok
sbatch scripts/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_allGenes.sh # 15955612 ok

# diff peaks
sbatch scripts/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q30_peak.sh # 15959030 ok
sbatch scripts/matrix_TSS_10kb_WTvsKO_H3K27me3_THOR_q50_peak.sh # 15959141 ok

sbatch scripts/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q30_peak.sh # 15959410 ok
sbatch scripts/matrix_TSS_10kb_WTvsKOEF1aEZH1_H3K27me3_THOR_q50_peak.sh # 15959511 ok

```

Generate gtf file from gene list; start with gene with peak in promoter (qval macs2 2.3):

```bash
# isolate all the genes bound with H3K27me3 in WT and or KO
cat output/ChIPseeker/annotation_macs2_PSC_WT_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt output/ChIPseeker/annotation_macs2_PSC_KO_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_macs2_PSC_WTKO_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt
cat output/ChIPseeker/annotation_macs2_PSC_WT_H3K27me3_qval1.30103_promoterAnd5_geneSymbol.txt output/ChIPseeker/annotation_macs2_PSC_KO_H3K27me3_qval1.30103_promoterAnd5_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_macs2_PSC_WTKO_H3K27me3_qval1.30103_promoterAnd5_geneSymbol.txt
### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
## Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_macs2_PSC_WTKO_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_WTKO_H3K27me3_qval2.30103_Promoter_5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_macs2_PSC_WTKO_H3K27me3_qval1.30103_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_WTKO_H3K27me3_qval1.30103_Promoter_5_as_gtf_geneSymbol.txt
## Filter the gtf
grep -Ff output/ChIPseeker/annotation_WTKO_H3K27me3_qval2.30103_Promoter_5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_macs2_H3K27me3_WTKO_qval2.30103_Promoter_5.gtf
grep -Ff output/ChIPseeker/annotation_WTKO_H3K27me3_qval1.30103_Promoter_5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_macs2_H3K27me3_WTKO_qval1.30103_Promoter_5.gtf


# deeptool plots
sbatch scripts/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q1.30103.sh # 15963978 ok
sbatch scripts/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103.sh # 15964044 ok



```

- *NOTE: qval1.3 was defined as optimal for macs2*
- *NOTE: qval30 (or 05) was defined as optimal for THOR*






## Assign peak to genes for NPC and PSC - THOR peaks


Isolate positive and negative THOR peaks
```bash
# positive negative peaks
## qval 30
awk -F'\t' '$14 > 1' output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval30.bed > output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval30_positive.bed
awk -F'\t' '$14 < 1' output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval30.bed > output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval30_negative.bed

## qval 50
awk -F'\t' '$14 > 1' output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval50.bed > output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval50_positive.bed
awk -F'\t' '$14 < 1' output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval50.bed > output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval50_negative.bed
```

```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")


# Import macs2 peaks
## WTvsKO _ q30
KO_gain = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval30_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
KO_lost = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval30_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)       
## WTvsKO _ q50
KO_gain = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval50_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
KO_lost = as_tibble(read.table('output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval50_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)       


## Tidy peaks #-->> Re-Run from here with different qvalue!!
KO_gain_gr = makeGRangesFromDataFrame(KO_gain,keep.extra.columns=TRUE)
KO_lost_gr = makeGRangesFromDataFrame(KO_lost,keep.extra.columns=TRUE)
gr_list <- list(KO_gain=KO_gain_gr, KO_lost=KO_lost_gr)
## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
KO_gain_annot <- as.data.frame(peakAnnoList[["KO_gain"]]@anno)
KO_lost_annot <- as.data.frame(peakAnnoList[["KO_lost"]]@anno)

## Convert entrez gene IDs to gene symbols
KO_gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KO_gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KO_gain_annot$gene <- mapIds(org.Hs.eg.db, keys = KO_gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KO_lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KO_lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KO_lost_annot$gene <- mapIds(org.Hs.eg.db, keys = KO_lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(KO_gain_annot, file="output/ChIPseeker/annotation_THOR_KO_gain_annot_qval50.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(KO_lost_annot, file="output/ChIPseeker/annotation_THOR_KO_lost_annot_qval50.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
KO_gain_annot_promoterAnd5 = tibble(KO_gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
KO_lost_annot_promoterAnd5 = tibble(KO_lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
KO_gain_annot_promoterAnd5_geneSymbol = KO_gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
KO_lost_annot_promoterAnd5_geneSymbol = KO_lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(KO_gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_THOR_KO_gain_qval50_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(KO_lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_THOR_KO_lost_qval50_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


   
```




# THOR with ChIPseeker

Generate xls file with gene name and THOR peak metrics (Naiara Slack task 20240314):
- import THOR bed output (correct qvalue) `output/THOR/THOR*/THOR_qval*.bed` (q50)
- import ChIPseeker output `output/ChIPseeker/annotation_THOR*.txt`
- combine THOR and ChIPseeker with Peak name
- output all metrics in `THOR` folder! Then filter them in xls to keep only the relevant ones


```R
library("tidyverse")

# import files
THOR = read_tsv("output/THOR/THOR_PSC_WTvsKO_H3K27me3/THOR_qval50.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character())) %>%
       rename("name" = "X4", "qval" = "X13", "FC" = "X14", "WT_count" = "X11", "KO_count" = "X12") %>%
       dplyr::select(name, WT_count, KO_count, FC, qval)

chipseeker_gain = read_tsv("output/ChIPseeker/annotation_THOR_KO_gain_annot_qval50.txt",
                      col_names = TRUE, trim_ws = TRUE) %>%
       dplyr::select(seqnames, start, end, width, name, annotation, geneChr, geneStart, geneEnd, geneLength, geneStrand, geneId, transcriptId, distanceToTSS, geneSymbol, gene) %>%
       add_column(peak = "gain")
chipseeker_lost = read_tsv("output/ChIPseeker/annotation_THOR_KO_lost_annot_qval50.txt",
                      col_names = TRUE, trim_ws = TRUE) %>%
       dplyr::select(seqnames, start, end, width, name, annotation, geneChr, geneStart, geneEnd, geneLength, geneStrand, geneId, transcriptId, distanceToTSS, geneSymbol, gene) %>%
       add_column(peak = "lost")
chipseeker = chipseeker_gain %>%
    bind_rows(chipseeker_lost)
# combine and filter

THOR_chipseeker = THOR %>% left_join(chipseeker)


write.table(THOR_chipseeker, file="output/THOR/THOR_PSC_WTvsKO_H3K27me3/THORq50_chipseeker-PSC_WTvsKO.txt", sep="\t", quote=FALSE, row.names=FALSE)



```





# Quantify signal around TSS 

--> This is to be used for Conchi project gastrulation paper `008*/001*` at `## Test Activation point with H3K27me3 (from 001*/009*) - xxx`


- Generate a bed file around TSS of each gene (250bp, 500bp, 1kb up/downstream); file already generated at `001_*/003__*/meta/ENCFF159KBI_gene_1kbTSS_sorted.bed` 
- Quantify EZH2 read density around the TSS with [bin.bw](https://rdrr.io/github/jmonlong/PopSV/man/bin.bw.html) 
- Represent data in R `ggplot`

--> I need to generate a file with geneSymbol names and signal quantification. If several signal value per gene, take the highest one; because there could be transcript not express; so better to take the higher value; which should show the regulated one? **Let's generate median value of all transcript per gene; and maximum value of all transcript per gene.**


```bash
conda activate binBw_v2
```

```R
library("PopSV")
library("tidyverse")
library("Rsamtools")
library("ggpubr")
library("data.table")
library("biomaRt")

# H3K27me3 #####################################################


## WT H3K27me3_in Peaks 1kb
bwFile <- "output/THOR/THOR_PSC_WTvsKO_H3K27me3/PSCWTvsKOH3K27me3-s1-rep0.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_1kbTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_H3K27me3_1kb")
#### Collect output and gene information
WT_H3K27me3 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_H3K27me3_1kb.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- unique(WT_H3K27me3$gene)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)

WT_H3K27me3_geneSymbol = WT_H3K27me3 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()

write.table(WT_H3K27me3_geneSymbol, file = "output/binBw/WT_H3K27me3_1kbTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)




## WT H3K27me3_in Peaks 500bp
bwFile <- "output/THOR/THOR_PSC_WTvsKO_H3K27me3/PSCWTvsKOH3K27me3-s1-rep0.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_500bpTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_H3K27me3_500bp")
#### Collect output and gene information
WT_H3K27me3 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_H3K27me3_500bp.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- unique(WT_H3K27me3$gene)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)

WT_H3K27me3_geneSymbol = WT_H3K27me3 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()

write.table(WT_H3K27me3_geneSymbol, file = "output/binBw/WT_H3K27me3_500bpTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)




## WT H3K27me3_in Peaks 250bp
bwFile <- "output/THOR/THOR_PSC_WTvsKO_H3K27me3/PSCWTvsKOH3K27me3-s1-rep0.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_250bpTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_H3K27me3_250bp")
#### Collect output and gene information
WT_H3K27me3 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_H3K27me3_250bp.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- unique(WT_H3K27me3$gene)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)

WT_H3K27me3_geneSymbol = WT_H3K27me3 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()

write.table(WT_H3K27me3_geneSymbol, file = "output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)



```

--> Works well. File generated `output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt` containing max and median EZH2 signal for all genes.




# Correlation signal H3K27me3 signal TSS with RNA

Let's verify that our quantification of H3K27me3 is in agreement with gene expression.

```bash
conda activate deseq2
```



```R
# packages
library("tidyverse")
library("ggpubr")
library("pheatmap")
library("scales")
library("RColorBrewer")
library("corrplot")


# import files
## H3K27me3
WT_H3K27me3_1kbTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_1kbTSS_geneSymbol.txt")
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_500bpTSS_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt")
## gene with H3K27me3 peaks
gene_H3K27me3 = read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/ChIPseeker/annotation_macs2_PSC_WT_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt", col_name = FALSE) %>%
    dplyr::rename("geneSymbol"="X1")


## RNA 
tpm = read_tsv("../../001_EZH1_Project/001__RNAseq/output/tpm_hg38/tpm_all_sample_geneSymbol.txt") %>%
    dplyr::select("external_gene_name", "ESC_WT_R1","ESC_WT_R2","ESC_WT_R3")  %>%
    rowwise() %>%
    mutate(median_tpm = median(c_across(starts_with("ESC_WT_")))) %>%
    ungroup() %>%
    select(external_gene_name, median_tpm) %>%
    dplyr::rename("geneSymbol"="external_gene_name")


WT_H3K27me3_1kbTSS_geneSymbol_tpm = WT_H3K27me3_1kbTSS_geneSymbol %>%
    left_join(tpm) %>%
    drop_na() %>%
    unique()

# correlation scatterplot
pdf("output/binBw/corr_H3K27me3_RNA_WT.pdf", width=5, height=4)
WT_H3K27me3_1kbTSS_geneSymbol_tpm %>%
    mutate(median_tpm = log2(median_tpm+1),
           bc_max = log2(bc_max+1)) %>%
    filter(  ) %>%  # bc_max>0 , median_tpm <5000
ggplot(., aes(x = median_tpm, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 9000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


# heatmap RNA 1st column
## All express genes
heatmap_data <- WT_H3K27me3_1kbTSS_geneSymbol_tpm %>%
    mutate(median_tpm = log2(median_tpm+1),
           bc_max = log2(bc_max)) %>%
    filter( median_tpm > 0) %>%
    select(median_tpm, bc_max) %>%
    arrange(desc(median_tpm))
heatmap_data_scaled <- heatmap_data %>%
    mutate(
        median_tpm = rescale(median_tpm, to = c(0, 1)),
        bc_max = rescale(bc_max, to = c(0, 1))
    )
## Convert to matrix for pheatmap
pdf("output/binBw/heatmap_RNA_H3K27me3_WT.pdf", width=1, height=4)
pheatmap(
  as.matrix(heatmap_data_scaled),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = c(seq(0, 0.45, length.out = 40), seq(0.46, 0.54, length.out = 20), seq(0.55, 1, length.out = 40))
)
dev.off()

## All genes with H3K27me3
heatmap_data <- gene_H3K27me3 %>%
    left_join(WT_H3K27me3_1kbTSS_geneSymbol_tpm) %>%
    drop_na() %>%
    mutate(median_tpm = log2(median_tpm+1),
           bc_max = log2(bc_max)) %>%
    filter( median_tpm > 0) %>%
    select(median_tpm, bc_max) %>%
    arrange(desc(median_tpm))
heatmap_data_scaled <- heatmap_data %>%
    mutate(
        median_tpm = rescale(median_tpm, to = c(0, 1)),
        bc_max = rescale(bc_max, to = c(0, 1))
    )
## Convert to matrix for pheatmap
pdf("output/binBw/heatmap_RNA_H3K27me3_WT_gene_H3K27me3.pdf", width=1, height=4)
pheatmap(
  as.matrix(heatmap_data_scaled),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = c(seq(0, 0.65, length.out = 40), seq(0.66, 0.74, length.out = 20), seq(0.75, 1, length.out = 40))
)
dev.off()




# heatmap H3K27me3 1st column

heatmap_data <- WT_H3K27me3_1kbTSS_geneSymbol_tpm %>%
    mutate(median_tpm = log2(median_tpm+1),
           bc_max = log2(bc_max+1)) %>%
    filter( median_tpm > 0) %>%
    select(bc_max, median_tpm ) %>%
    arrange(desc(bc_max))

heatmap_data_scaled <- heatmap_data %>%
    mutate(
        median_tpm = rescale(median_tpm, to = c(0, 1)),
        bc_max = rescale(bc_max, to = c(0, 1))
    )

## Convert to matrix for pheatmap

pdf("output/binBw/heatmap_H3K27me3_RNA_WT.pdf", width=1, height=4)
pheatmap(
  as.matrix(heatmap_data_scaled),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = c(seq(0, 0.45, length.out = 40), seq(0.46, 0.54, length.out = 20), seq(0.55, 1, length.out = 40))
)
dev.off()



## Scatter plot with regression line
pdf("output/binBw/scatter_H3K27me3_RNA_WT.pdf", width = 5, height = 4)
WT_H3K27me3_1kbTSS_geneSymbol_tpm %>%
    filter( median_tpm > 0) %>%
ggplot(., aes(x = median_tpm, y = bc_max)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue") +
  scale_x_log10() +  # Use log scale for x-axis if needed
  scale_y_log10() +  # Use log scale for y-axis if needed
  theme_bw() +
  labs(x = "log2(tpm+1)", y = "log2(bc_max)")
dev.off()


## Correlation matrix

# Compute correlation matrix
corrMatrix_data = WT_H3K27me3_1kbTSS_geneSymbol_tpm %>%
    mutate(median_tpm = log2(median_tpm+1),
           bc_max = log2(bc_max+1)) %>%
    filter( median_tpm > 0)


cor_matrix <- cor(corrMatrix_data %>% select(median_tpm, bc_max))

# Correlation matrix heatmap
pdf("output/binBw/correlation_matrix_H3K27me3_RNA_WT.pdf", width = 4, height = 4)
corrplot(cor_matrix, method = "color", col = colorRampPalette(brewer.pal(11, "RdBu"))(200), 
 mar = c(0, 0, 1, 0), addCoef.col = "black")
dev.off()


```


--> Slight anticorrelation, seems it is working, but not extremely striking...
    --> Better when keeping only gene express (tpm >0)
    --> Selecting only genes with peak in WT is



