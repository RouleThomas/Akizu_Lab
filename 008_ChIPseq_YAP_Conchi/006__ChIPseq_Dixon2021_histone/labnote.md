# Project and goals 

Re-analysis of ChIPseq dataset from (Dixon et al)[10.1126/science.abd0875].

- Focus on histone marks. Detail in `gastrulation paper/Figure 4/*pptx`

# Download data


- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject [PRJNA715662](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169207)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh #  xxx
```

--> WT genotype differ:
- H3K27ac on HUES8 WT hESCs
- H3K4me1 on HUES8 WT hESCs
- H3K27me3 on H1 WT
- H3K36me3 on H1 WT

XXXY















--> Not clear what are the inputs with QSER1 files, so let's use the 2 rep of H9 inputs for all samples


## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename_008002.txt`, `rename_008002_QSER1.txt` and `rename_008002_EZH2.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_008002.txt

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_008002_QSER1.txt

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_008002_EZH2.txt
```

--> All good 




# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 18390922 ok
sbatch scripts/fastp_QSER1.sh # 18544351 ok
sbatch scripts/fastp_EZH2.sh # 18660122 xxx
```



# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:18390922 scripts/bowtie2_raw.sh # 18391156 ok
sbatch --dependency=afterany:18544351 scripts/bowtie2_QSER1.sh # 18544412 ok
sbatch --dependency=afterany:18660122 scripts/bowtie2_EZH2.sh # 18660130 xxx

```


--> Looks good; overall ~70% uniquely aligned reads
----> Seems less uniquel mapped reads than us but they sequence FAR more depth (~20m reads vs 5 for us)


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-18391156.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_18391156.txt
for file in slurm-18544412.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_18544412.txt

for file in slurm-18660130.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_18660130.txt
```

Add these values to `/home/roulet/008_ChIPseq_YAP_Conchi/002__ChIPseq_Dixon2021/samples_002.xlsx`\
Then in R; see `/home/roulet/008_ChIPseq_YAP_Conchi/008_ChIPseq_YAP_Conchi.R`.

--> Overall >70% input reads as been uniquely mapped to the genome (90% non uniq)



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.


```bash
conda activate bowtie2

sbatch --dependency=afterany:18391156 scripts/samtools_unique_raw.sh # 18391579 ok
sbatch --dependency=afterany:18544412 scripts/samtools_unique_QSER1.sh # 18544524 ok
sbatch --dependency=afterany:18660130 scripts/samtools_unique_EZH2.sh # 18660179 xxx

```


# Generate bigwig coverage files
## Raw and DiffBindTMM bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:18391579 scripts/bamtobigwig_unique_raw.sh # 18392081 ok
sbatch --dependency=afterany:18544524 scripts/bamtobigwig_unique_QSER1.sh # 18544590 ok
sbatch --dependency=afterany:18660179 scripts/bamtobigwig_unique_EZH2.sh # 18660185 ok

# Bigwig with DiffBind TMM scaling factor
sbatch scripts/bamtobigwig_unique_hESC_QSER1_H3K27me3.sh # 18760307 ok

```

--> **DiffBind_TMM SF works great** at putting replicate together!!

PASS: H3K27me3, H3K4me3, EZH2, QSER1FLAG; DNMT3A and B a bit messy


Generate median tracks:
```bash
conda activate BedToBigwig
# raw unique bigwig
sbatch scripts/bigwigmerge_unique_raw.sh # 18543995 fail
sbatch --dependency=afterany:18544590 scripts/bigwigmerge_unique_QSER1.sh # 18544696 fail

# rerun unique and QSER1 together
sbatch scripts/bigwigmerge_unique_raw.sh # 18717741 ok
```

*NOTE: merging raw bigiwg probably not super smart; not seq depth normalized, so to take with caution! But to show presence of signal in our case, that's ok!*


## Pearson correlation heatmap on bigwig signals


```bash
conda activate deeptools

### Generate compile bigwig (.npz) files _ all raw bigwig
sbatch scripts/multiBigwigSummary_all.sh # 18661462 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_DNMT3A_R1 hESC_WT_DNMT3A_R2 hESC_WT_DNMT3B_R1 hESC_WT_DNMT3B_R2 hESC_WT_H3K27me3_R1 hESC_WT_H3K27me3_R2 hESC_WT_H3K4me3_R1 hESC_WT_H3K4me3_R2 hESC_WT_EZH2_R1 hESC_WT_input_R1 hESC_WT_input_R2 \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf


## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_DNMT3A_R1 hESC_WT_DNMT3A_R2 hESC_WT_DNMT3B_R1 hESC_WT_DNMT3B_R2 hESC_WT_H3K27me3_R1 hESC_WT_H3K27me3_R2 hESC_WT_H3K4me3_R1 hESC_WT_H3K4me3_R2 hESC_WT_EZH2_R1 hESC_WT_input_R1 hESC_WT_input_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf
```

--> Replicates/samples are very clean; cluster well together as expected



# MACS2 peak calling on bam unique



--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
sbatch scripts/macs2_broad_all.sh # 18601898 ok
sbatch scripts/macs2_broad_EZH2.sh # 18661592 ok

# pool
sbatch scripts/macs2_broad_all_pool.sh # 18601982 fail; rerun 18603726 ok 
```

- *NOTE: I forget added `*pool` suffix at the name file at job-18601982, and forget to run in broad; rerun.*



```bash
conda activate bowtie2 # for bedtools
# sbatch scripts/macs2_raw_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive
sbatch scripts/macs2_raw_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive


# quick command to print median size of peak within a bed
awk '{print $3-$2}' your_bed_file.bed | sort -n | awk 'BEGIN {c=0; sum=0;} {a[c++]=$1; sum+=$1;} END {if (c%2) print a[int(c/2)]; else print (a[c/2-1]+a[c/2])/2;}'
```

Then keep only the significant peaks (re-run the script to test different qvalue cutoff) and remove peaks overlapping with blacklist regions. MACS2 column9 output is -log10(qvalue) format so if we want 0.05; 
- q0.05: `q value = -log10(0.05) = 1.30103`
- q0.01 = 2
- q0.005 = 2.30103
- q0.001 = 3
- q0.0001 = 4
- q0.00001 = 5


**Optimal qvalue** according to IGV:
- hESC_WT_DNMT3A: 3 noisy (4 even better)
- hESC_WT_DNMT3B: 3 noisy (4 even better)
- hESC_WT_H3K27me3: 2.30103 
- hESC_WT_H3K4me3:  2.30103 
- hESC_WT_QSER1FLAG: 3
- hESC_WT_EZH2: 1.30103



# ChIPseeker peak gene assignment

## From optimal qval bed files peaks
Let's assign **peak to genes from MACS2 peak**:

**Optimal qvalue** according to IGV:
- hESC_WT_DNMT3A: 3 noisy (4 even better)
- hESC_WT_DNMT3B: 3 noisy (4 even better)
- hESC_WT_H3K27me3: 2.30103 
- hESC_WT_H3K4me3:  2.30103 
- hESC_WT_QSER1FLAG: 3
- hESC_WT_EZH2: 1.30103



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
## optimal qval
hESC_WT_DNMT3A = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/hESC_WT_DNMT3A_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_DNMT3B = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/hESC_WT_DNMT3B_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_H3K27me3 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/hESC_WT_H3K27me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_H3K4me3 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/hESC_WT_H3K4me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1FLAG = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_EZH2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_EZH2_R1_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 




# Tidy peaks 
## H3K27me3
hESC_WT_DNMT3A_gr = makeGRangesFromDataFrame(hESC_WT_DNMT3A,keep.extra.columns=TRUE)
hESC_WT_DNMT3B_gr = makeGRangesFromDataFrame(hESC_WT_DNMT3B,keep.extra.columns=TRUE)
hESC_WT_H3K27me3_gr = makeGRangesFromDataFrame(hESC_WT_H3K27me3,keep.extra.columns=TRUE)
hESC_WT_H3K4me3_gr = makeGRangesFromDataFrame(hESC_WT_H3K4me3,keep.extra.columns=TRUE)
hESC_WT_QSER1FLAG_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG,keep.extra.columns=TRUE)
hESC_WT_EZH2_gr = makeGRangesFromDataFrame(hESC_WT_EZH2,keep.extra.columns=TRUE)

gr_list <- list(hESC_WT_DNMT3A=hESC_WT_DNMT3A_gr, hESC_WT_DNMT3B=hESC_WT_DNMT3B_gr,  hESC_WT_H3K27me3=hESC_WT_H3K27me3_gr, hESC_WT_H3K4me3=hESC_WT_H3K4me3_gr,hESC_WT_QSER1FLAG=hESC_WT_QSER1FLAG_gr, hESC_WT_EZH2=hESC_WT_EZH2_gr)


# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_hESC_WT_qvalOPtimal.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_hESC_WT_qvalOPtimal.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
hESC_WT_DNMT3A_annot <- as.data.frame(peakAnnoList[["hESC_WT_DNMT3A"]]@anno)
hESC_WT_DNMT3B_annot <- as.data.frame(peakAnnoList[["hESC_WT_DNMT3B"]]@anno)
hESC_WT_H3K27me3_annot <- as.data.frame(peakAnnoList[["hESC_WT_H3K27me3"]]@anno)
hESC_WT_H3K4me3_annot <- as.data.frame(peakAnnoList[["hESC_WT_H3K4me3"]]@anno)
hESC_WT_QSER1FLAG_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG"]]@anno)
hESC_WT_EZH2_annot <- as.data.frame(peakAnnoList[["hESC_WT_EZH2"]]@anno)

## Convert entrez gene IDs to gene symbols
hESC_WT_DNMT3A_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_DNMT3A_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_DNMT3A_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_DNMT3A_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_DNMT3B_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_DNMT3B_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_DNMT3B_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_DNMT3B_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_H3K4me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_H3K4me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_H3K4me3_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_H3K4me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_EZH2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_EZH2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(hESC_WT_DNMT3A_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_DNMT3A_pool_qval3.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_DNMT3B_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_DNMT3B_pool_qval3.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_H3K27me3_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_H3K27me3_pool_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_H3K4me3_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_H3K4me3_pool_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1FLAG_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval3.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_EZH2_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE




## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_DNMT3A_annot_promoterAnd5 = tibble(hESC_WT_DNMT3A_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_DNMT3B_annot_promoterAnd5 = tibble(hESC_WT_DNMT3B_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_H3K27me3_annot_promoterAnd5 = tibble(hESC_WT_H3K27me3_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_H3K4me3_annot_promoterAnd5 = tibble(hESC_WT_H3K4me3_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1FLAG_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_EZH2_annot_promoterAnd5 = tibble(hESC_WT_EZH2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
hESC_WT_DNMT3A_annot_promoterAnd5_geneSymbol = hESC_WT_DNMT3A_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_DNMT3B_annot_promoterAnd5_geneSymbol = hESC_WT_DNMT3B_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_H3K27me3_annot_promoterAnd5_geneSymbol = hESC_WT_H3K27me3_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_H3K4me3_annot_promoterAnd5_geneSymbol = hESC_WT_H3K4me3_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_annot_promoterAnd5_geneSymbol = hESC_WT_EZH2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(hESC_WT_DNMT3A_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_DNMT3A_pool_qval3_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_DNMT3B_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_DNMT3B_pool_qval3_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_H3K27me3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_H3K27me3_pool_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_H3K4me3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_H3K4me3_pool_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1FLAG_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval3_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_EZH2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)



## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_DNMT3A_annot_noIntergenic = tibble(hESC_WT_DNMT3A_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_DNMT3B_annot_noIntergenic = tibble(hESC_WT_DNMT3B_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_H3K27me3_annot_noIntergenic = tibble(hESC_WT_H3K27me3_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_H3K4me3_annot_noIntergenic = tibble(hESC_WT_H3K4me3_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1FLAG_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_EZH2_annot_noIntergenic = tibble(hESC_WT_EZH2_annot) %>%
    filter(annotation != c("Distal Intergenic"))


### Save output gene lists
hESC_WT_DNMT3A_annot_noIntergenic_geneSymbol = hESC_WT_DNMT3A_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_DNMT3B_annot_noIntergenic_geneSymbol = hESC_WT_DNMT3B_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_H3K27me3_annot_noIntergenic_geneSymbol = hESC_WT_H3K27me3_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_H3K4me3_annot_noIntergenic_geneSymbol = hESC_WT_H3K4me3_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_annot_noIntergenic_geneSymbol = hESC_WT_EZH2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(hESC_WT_DNMT3A_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_DNMT3A_pool_qval3_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_DNMT3B_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_DNMT3B_pool_qval3_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_H3K27me3_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_H3K27me3_pool_qval2.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_H3K4me3_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_H3K4me3_pool_qval2.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1FLAG_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval3_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_EZH2_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)





## QSER1FLAG qvalue testing
hESC_WT_QSER1FLAG_qval4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval4/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1FLAG_qval5 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval5/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1FLAG_qval6 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval6/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1FLAG_qval7 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval7/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1FLAG_qval10 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval10/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1FLAG_qval12 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval12/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1FLAG_qval15 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval15/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1FLAG_qval20 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval20/hESC_WT_QSER1FLAG_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 


# Tidy peaks 
hESC_WT_QSER1FLAG_qval4_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG_qval4,keep.extra.columns=TRUE)
hESC_WT_QSER1FLAG_qval5_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG_qval5,keep.extra.columns=TRUE)
hESC_WT_QSER1FLAG_qval6_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG_qval6,keep.extra.columns=TRUE)
hESC_WT_QSER1FLAG_qval7_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG_qval7,keep.extra.columns=TRUE)
hESC_WT_QSER1FLAG_qval10_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG_qval10,keep.extra.columns=TRUE)
hESC_WT_QSER1FLAG_qval12_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG_qval12,keep.extra.columns=TRUE)
hESC_WT_QSER1FLAG_qval15_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG_qval15,keep.extra.columns=TRUE)
hESC_WT_QSER1FLAG_qval20_gr = makeGRangesFromDataFrame(hESC_WT_QSER1FLAG_qval20,keep.extra.columns=TRUE)

gr_list <- list(hESC_WT_QSER1FLAG_qval4=hESC_WT_QSER1FLAG_qval4_gr, hESC_WT_QSER1FLAG_qval5=hESC_WT_QSER1FLAG_qval5_gr,  hESC_WT_QSER1FLAG_qval6=hESC_WT_QSER1FLAG_qval6_gr, hESC_WT_QSER1FLAG_qval7=hESC_WT_QSER1FLAG_qval7_gr,hESC_WT_QSER1FLAG_qval10=hESC_WT_QSER1FLAG_qval10_gr,hESC_WT_QSER1FLAG_qval12=hESC_WT_QSER1FLAG_qval12_gr,hESC_WT_QSER1FLAG_qval15=hESC_WT_QSER1FLAG_qval15_gr,hESC_WT_QSER1FLAG_qval20=hESC_WT_QSER1FLAG_qval20_gr)


# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
#pdf("output/ChIPseeker/plotAnnoBar_hESC_WT_qvalOPtimal.pdf", width = 8, height = 3)
pdf("output/ChIPseeker/plotAnnoBar_hESC_WT_QSER1FLAG_qvalTest.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
#pdf("output/ChIPseeker/plotDistToTSS_hESC_WT_qvalOPtimal.pdf", width = 8, height = 3)
pdf("output/ChIPseeker/plotDistToTSS_hESC_WT_QSER1FLAG_qvalTest.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
hESC_WT_QSER1FLAG_qval4_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG_qval4"]]@anno)
hESC_WT_QSER1FLAG_qval5_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG_qval5"]]@anno)
hESC_WT_QSER1FLAG_qval6_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG_qval6"]]@anno)
hESC_WT_QSER1FLAG_qval7_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG_qval7"]]@anno)
hESC_WT_QSER1FLAG_qval10_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG_qval10"]]@anno)
hESC_WT_QSER1FLAG_qval12_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG_qval12"]]@anno)
hESC_WT_QSER1FLAG_qval15_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG_qval15"]]@anno)
hESC_WT_QSER1FLAG_qval20_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1FLAG_qval20"]]@anno)


## Convert entrez gene IDs to gene symbols
hESC_WT_QSER1FLAG_qval4_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval4_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval4_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval4_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval5_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval5_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval5_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval5_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval6_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval6_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval6_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval6_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval7_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval7_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval7_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval7_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval10_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval10_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval10_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval10_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval12_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval12_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval12_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval12_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval15_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval15_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval15_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval15_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval20_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval20_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1FLAG_qval20_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1FLAG_qval20_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(hESC_WT_QSER1FLAG_qval4_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval4.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1FLAG_qval5_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval5.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1FLAG_qval6_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval6.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1FLAG_qval7_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval7.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1FLAG_qval10_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval10.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1FLAG_qval12_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval12.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1FLAG_qval15_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval15.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1FLAG_qval20_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval20.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_QSER1FLAG_qval4_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_qval4_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1FLAG_qval5_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_qval5_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1FLAG_qval6_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_qval6_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1FLAG_qval7_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_qval7_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1FLAG_qval10_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_qval10_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1FLAG_qval12_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_qval12_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1FLAG_qval15_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_qval15_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1FLAG_qval20_annot_promoterAnd5 = tibble(hESC_WT_QSER1FLAG_qval20_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))



### Save output gene lists
hESC_WT_QSER1FLAG_qval4_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_qval4_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval5_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_qval5_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval6_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_qval6_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval7_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_qval7_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval10_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_qval10_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval12_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_qval12_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval15_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_qval15_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval20_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1FLAG_qval20_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(hESC_WT_QSER1FLAG_qval4_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval4_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1FLAG_qval5_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval5_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)           
write.table(hESC_WT_QSER1FLAG_qval6_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval6_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)      
write.table(hESC_WT_QSER1FLAG_qval7_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval7_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   
write.table(hESC_WT_QSER1FLAG_qval10_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval10_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   
write.table(hESC_WT_QSER1FLAG_qval12_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval12_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   
write.table(hESC_WT_QSER1FLAG_qval15_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval15_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   
write.table(hESC_WT_QSER1FLAG_qval20_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval20_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   


## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_QSER1FLAG_qval4_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_qval4_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1FLAG_qval5_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_qval5_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1FLAG_qval6_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_qval6_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1FLAG_qval7_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_qval7_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1FLAG_qval10_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_qval10_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1FLAG_qval12_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_qval12_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1FLAG_qval15_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_qval15_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1FLAG_qval20_annot_noIntergenic = tibble(hESC_WT_QSER1FLAG_qval20_annot) %>%
    filter(annotation != c("Distal Intergenic"))



### Save output gene lists
hESC_WT_QSER1FLAG_qval4_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_qval4_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval5_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_qval5_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval6_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_qval6_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval7_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_qval7_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval10_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_qval10_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval12_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_qval12_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval15_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_qval15_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1FLAG_qval20_annot_noIntergenic_geneSymbol = hESC_WT_QSER1FLAG_qval20_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(hESC_WT_QSER1FLAG_qval4_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval4_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1FLAG_qval5_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval5_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)           
write.table(hESC_WT_QSER1FLAG_qval6_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval6_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)      
write.table(hESC_WT_QSER1FLAG_qval7_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval7_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   
write.table(hESC_WT_QSER1FLAG_qval10_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval10_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   
write.table(hESC_WT_QSER1FLAG_qval12_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval12_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   
write.table(hESC_WT_QSER1FLAG_qval15_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval15_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   
write.table(hESC_WT_QSER1FLAG_qval20_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval20_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)   





## EZH2 qvalue testing
hESC_WT_EZH2_qval2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/hESC_WT_EZH2_R1_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_EZH2_qval3 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/hESC_WT_EZH2_R1_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_EZH2_qval4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval4/hESC_WT_EZH2_R1_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
    
# Tidy peaks 
hESC_WT_EZH2_qval2_gr = makeGRangesFromDataFrame(hESC_WT_EZH2_qval2,keep.extra.columns=TRUE)
hESC_WT_EZH2_qval3_gr = makeGRangesFromDataFrame(hESC_WT_EZH2_qval3,keep.extra.columns=TRUE)
hESC_WT_EZH2_qval4_gr = makeGRangesFromDataFrame(hESC_WT_EZH2_qval4,keep.extra.columns=TRUE)

gr_list <- list(hESC_WT_EZH2_qval2=hESC_WT_EZH2_qval2_gr, hESC_WT_EZH2_qval3=hESC_WT_EZH2_qval3_gr,  hESC_WT_EZH2_qval4=hESC_WT_EZH2_qval4_gr)


# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
#pdf("output/ChIPseeker/plotAnnoBar_hESC_WT_qvalOPtimal.pdf", width = 8, height = 3)
pdf("output/ChIPseeker/plotAnnoBar_hESC_WT_EZH2_qvalTest.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
#pdf("output/ChIPseeker/plotDistToTSS_hESC_WT_qvalOPtimal.pdf", width = 8, height = 3)
pdf("output/ChIPseeker/plotDistToTSS_hESC_WT_EZH2_qvalTest.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
hESC_WT_EZH2_qval2_annot <- as.data.frame(peakAnnoList[["hESC_WT_EZH2_qval2"]]@anno)
hESC_WT_EZH2_qval3_annot <- as.data.frame(peakAnnoList[["hESC_WT_EZH2_qval3"]]@anno)
hESC_WT_EZH2_qval4_annot <- as.data.frame(peakAnnoList[["hESC_WT_EZH2_qval4"]]@anno)


## Convert entrez gene IDs to gene symbols
hESC_WT_EZH2_qval2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_qval2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_EZH2_qval2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_qval2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_EZH2_qval3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_qval3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_EZH2_qval3_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_qval3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_EZH2_qval4_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_qval4_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_EZH2_qval4_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_qval4_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(hESC_WT_EZH2_qval2_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_EZH2_qval3_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval3.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_EZH2_qval4_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval4.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_EZH2_qval2_annot_promoterAnd5 = tibble(hESC_WT_EZH2_qval2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_EZH2_qval3_annot_promoterAnd5 = tibble(hESC_WT_EZH2_qval3_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_EZH2_qval4_annot_promoterAnd5 = tibble(hESC_WT_EZH2_qval4_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
    



### Save output gene lists
hESC_WT_EZH2_qval2_annot_promoterAnd5_geneSymbol = hESC_WT_EZH2_qval2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_qval3_annot_promoterAnd5_geneSymbol = hESC_WT_EZH2_qval3_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_qval4_annot_promoterAnd5_geneSymbol = hESC_WT_EZH2_qval4_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
    

write.table(hESC_WT_EZH2_qval2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_EZH2_qval3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval3_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)           
write.table(hESC_WT_EZH2_qval4_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval4_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)  
            

## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_EZH2_qval2_annot_noIntergenic = tibble(hESC_WT_EZH2_qval2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_EZH2_qval3_annot_noIntergenic = tibble(hESC_WT_EZH2_qval3_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_EZH2_qval4_annot_noIntergenic = tibble(hESC_WT_EZH2_qval4_annot) %>%
    filter(annotation != c("Distal Intergenic"))
    


### Save output gene lists
hESC_WT_EZH2_qval2_annot_noIntergenic_geneSymbol = hESC_WT_EZH2_qval2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_qval3_annot_noIntergenic_geneSymbol = hESC_WT_EZH2_qval3_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_qval4_annot_noIntergenic_geneSymbol = hESC_WT_EZH2_qval4_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
    
write.table(hESC_WT_EZH2_qval2_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval2.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_EZH2_qval3_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval3_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)           
write.table(hESC_WT_EZH2_qval4_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval4_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE) 
            


```


--> So many genes for QSER1FLAG! ~20k, and there is ~30k total genes in human... Tested more stringeant qval.
----> qval 15-20 show good peak-gene feature profile; similar to what we have in `008001` with endogeneous QSER1; also ~10k genes, as we have in `008001`

--> Overlap pretty good with our `008001`, notably using QSERFLAG qval10-20!

--> Make deepTool plot of genes QSER1 bound in Conchi `008001` vs Dixon `008002`: Done in `008001` labnote






# DiffBind TMM normalization (for sample without spikein)


Generate DiffBind TMM scaling factor and apply them to the bam. That will generate clean seq + complexity depth normalized bigwigs. --> **Exact same process as E coli spike in norm, just no library size correction!**

--> Apply Reciprocal SF to bamtobgiwg!! (see `001003` for detail at `### Histone-spike-in and TMM-norm scaled-bigwig OR DiffBind-ScaleFactor`)


- hESC QSER1FLAG WT
- hESC H3K27me3 WT


```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 


## hESC QSER1FLAG for WT 
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1FLAG.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1FLAG.RData")
load("output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1FLAG.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_QSER1FLAG.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_QSER1FLAG.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_TMM_SF = dba.normalize(sample_count_blackgreylist_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_TMM_unique_SF_hESC_QSER1FLAG.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_QSER1FLAG_blackgreylist_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_QSER1FLAG_blackgreylist_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()




## hESC H3K27me3 for WT 
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_hESC_H3K27me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/meta_sample_macs2raw_unique_hESC_H3K27me3.RData")
load("output/DiffBind/meta_sample_macs2raw_unique_hESC_H3K27me3.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_TMM_SF = dba.normalize(sample_count_blackgreylist_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_TMM_unique_SF_hESC_H3K27me3.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_H3K27me3_blackgreylist_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_H3K27me3_blackgreylist_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()


```

Can we **use different AB in the same DiffBind correction??** Will the SF provided be the same as AB per AB??
--> NOT the same number; need to be run AB per AB!!


















