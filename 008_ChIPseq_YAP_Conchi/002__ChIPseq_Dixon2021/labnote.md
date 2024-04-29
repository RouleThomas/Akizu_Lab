# Project and goals 

Re-analysis of CutRun dataset from (Dixon et al)[10.1126/science.abd0875].

- Focus on H3K4me3, H3K27me3 and DNMT3 (check if co-localization with EZH2 and H3K27me3 from `008*/001*`)
--> 2 rep WT H1 hESC


# Download data


- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject  PRJNA631028 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 18384845 xxx

```



## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename_008002.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_008002.txt
```

--> All good 




# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 18390922 xxx
```



# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:18390922 scripts/bowtie2_raw.sh # 18391156 xxx
```

XXXXXXXXXXXXXXXXXXXXXXX BELOw XXXXXXXXXXXXXXXXXXX

--> Looks good; overall ~30-80% uniquely aligned reads
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
```

Add these values to `/home/roulet/008_ChIPseq_YAP_Conchi/002__ChIPseq_Dixon2021/samples_002.xlsx`\
Then in R; see `/home/roulet/008_ChIPseq_YAP_Conchi/008_ChIPseq_YAP_Conchi.R`.

--> Overall >60% input reads as been uniquely mapped to the genome (90% non uniq)



XXXXXXXXXXXXXXXXXXXXXXX Up XXXXXXXXXXXXXXXXXXX

## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.


```bash
conda activate bowtie2

sbatch --dependency=afterany:18391156 scripts/samtools_unique_raw.sh # 18391579 xxx
```


# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:18391579 scripts/bamtobigwig_unique_raw.sh # 18392081 xxx

```

XXXXXXXXXXXXXXXXXXXXXXX BELOw  ALL XXXXXXXXXXXXXXXXXXX




- NPC
PASS: H3K4m3 (rep very diff.), H3K9me3 (a bit noisy), H3K27me3
FAIL: H3K27ac (very low signal and noisy, seems R2 work better)
- 53dN
PASS: H3K4me3, H3K9me3 (a bit noisy), H3K27ac, H3K27me3
FAIL: *H3K9me3* could be there





--> The failed one, are also failed in the bigwig Ciceri files...

--> Compare with H3K27ac from (ENCODE)[https://www.encodeproject.org/experiments/ENCSR799SRL/]


```bash
# gunzip the file

# Convert wig to bigwig
srun --mem=500g --pty bash -l

## install wigtobigwig
conda activate BedToBigwig
conda install bioconda::ucsc-wigtobigwig # fail
conda install bioconda/label/cf201901::ucsc-wigtobigwig  # fail
#### --> fail create a new conda env
conda create -n wigtobigwig -c bioconda ucsc-wigtobigwig
conda activate wigtobigwig

## convert wig to bigwig
wigToBigWig output/bigwig_hg19/GSM767343_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.SK504.wig ../../Master/meta/hg19.chrom.sizes output/bigwig_hg19/GSM767343_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.SK504.bw
wigToBigWig output/bigwig_hg19/GSM818031_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.AK220.wig ../../Master/meta/hg19.chrom.sizes output/bigwig_hg19/GSM818031_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.AK220.bw
wigToBigWig output/bigwig_hg19/GSM896162_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.AK319.wig ../../Master/meta/hg19.chrom.sizes output/bigwig_hg19/GSM896162_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.AK319.bw


```
NOTE: hg19 chrom size copy from [ucsc](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/)





Generate median tracks:
```bash
conda activate BedToBigwig
# raw unique bigwig
sbatch scripts/bigwigmerge_unique_NPC.sh # 17509895 ok
sbatch scripts/bigwigmerge_unique_53dN.sh # 17509979 ok
```

*NOTE: merging raw bigiwg probably not super smart; not seq depth normalized, so to take with caution! But to show presence of signal in our case, that's ok!*


## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools

# Generate compile bigwig (.npz) files _ hg38 Akizu analysis
sbatch scripts/multiBigwigSummary_NPC.sh # 15280856 ok
sbatch scripts/multiBigwigSummary_53dN.sh # 15280864 ok

# Generate compile bigwig (.npz) files _ hg19 Ciceri analysis
sbatch scripts/multiBigwigSummary_NPC_Ciceri.sh # 15280901 ok
sbatch scripts/multiBigwigSummary_53dN_Ciceri.sh # 15281092 ok


# Akizu with Ciceri _ NPC
sbatch scripts/multiBigwigSummary_NPC_CutRun001008_Ciceri.sh # 15290030 ok



```

**Good to use**:
- *NPC*: H3K27me3, H3K4me3, H3K9me3, IGG
- *53dN*: H3K27me3, H3K9me3

**Bad to use**:
- *NPC*: H3K27ac: no signal! Like IGG
- *53dN*: AB mix between IGG, H3K4me3, H3K27ac

--> Same observation between Akizu and Ciceri analysis(hg38 hg19)

--> The *good to use* ones nicely correlate with our data in NPC WT `CutRun__005008` (`CutRun__009`)



# MACS2 peak calling on bam unique



--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad_53dN.sh # 15401244 ok
sbatch scripts/macs2_broad_NPC.sh # 15401308 ok

sbatch scripts/macs2_broad_53dN_noIGG.sh # 15401429 ok

# pool
sbatch scripts/macs2_broad_53dN_pool.sh # 17503351 ok
sbatch scripts/macs2_broad_NPC_pool.sh # 17503352 ok


```

--> H3K27ac in NPC show 3 peak in R1! And a ~7k peaks in R2. R2 is better, but still very ugly and noisy!

--> 53dN IGG vs not using IGG: almost the same, so let's better use IGG



```bash
conda activate bowtie2 # for bedtools
# sbatch scripts/macs2_raw_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive
sbatch scripts/macs2_raw_peak_signif_pool.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive


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
- 50dN_WT_H3K4me3: 2.30103 (4 more true peaks)
- 50dN_WT_H3K27me3: 2.30103 (4 more true peaks); data is of poor quality...
- *other samples not checked*




# ChIPseeker peak gene assignment

## From optimal qval bed files peaks
Let's assign **peak to genes from MACS2 peak**:

**Optimal qvalue** according to IGV:
- 50dN_WT_H3K4me3: 2.30103 (4 more true peaks)
- 50dN_WT_H3K27me3: 2.30103 (4 more true peaks); data is of poor quality...
- *other samples not checked*



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
## H3K27me3 53dN
H3K27me3_53dN_pool_qval2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/53dN_WT_H3K27me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K4me3_53dN_pool_qval2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/53dN_WT_H3K4me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K27me3_53dN_pool_qval4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval4/53dN_WT_H3K27me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
H3K4me3_53dN_pool_qval4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval4/53dN_WT_H3K4me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

# Tidy peaks 
## H3K27me3
H3K27me3_53dN_pool_qval2_gr = makeGRangesFromDataFrame(H3K27me3_53dN_pool_qval2,keep.extra.columns=TRUE)
H3K4me3_53dN_pool_qval2_gr = makeGRangesFromDataFrame(H3K4me3_53dN_pool_qval2,keep.extra.columns=TRUE)
H3K27me3_53dN_pool_qval4_gr = makeGRangesFromDataFrame(H3K27me3_53dN_pool_qval4,keep.extra.columns=TRUE)
H3K4me3_53dN_pool_qval4_gr = makeGRangesFromDataFrame(H3K4me3_53dN_pool_qval4,keep.extra.columns=TRUE)
gr_list <- list(H3K27me3_53dN_pool_qval2=H3K27me3_53dN_pool_qval2_gr, H3K4me3_53dN_pool_qval2=H3K4me3_53dN_pool_qval2_gr,  H3K27me3_53dN_pool_qval4=H3K27me3_53dN_pool_qval4_gr, H3K4me3_53dN_pool_qval4=H3K4me3_53dN_pool_qval4_gr)


# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_53dN_bivalent_pool.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_53dN_bivalent_pool.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
H3K27me3_53dN_pool_qval2_annot <- as.data.frame(peakAnnoList[["H3K27me3_53dN_pool_qval2"]]@anno)
H3K4me3_53dN_pool_qval2_annot <- as.data.frame(peakAnnoList[["H3K4me3_53dN_pool_qval2"]]@anno)
H3K27me3_53dN_pool_qval4_annot <- as.data.frame(peakAnnoList[["H3K27me3_53dN_pool_qval4"]]@anno)
H3K4me3_53dN_pool_qval4_annot <- as.data.frame(peakAnnoList[["H3K4me3_53dN_pool_qval4"]]@anno)


## Convert entrez gene IDs to gene symbols
H3K27me3_53dN_pool_qval2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_53dN_pool_qval2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_53dN_pool_qval2_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_53dN_pool_qval2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_53dN_pool_qval2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_53dN_pool_qval2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_53dN_pool_qval2_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_53dN_pool_qval2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_53dN_pool_qval4_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_53dN_pool_qval4_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_53dN_pool_qval4_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_53dN_pool_qval4_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_53dN_pool_qval4_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_53dN_pool_qval4_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_53dN_pool_qval4_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_53dN_pool_qval4_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



## Save output table
write.table(H3K27me3_53dN_pool_qval2_annot, file="output/ChIPseeker/annotation_macs2_H3K27me3_53dN_pool_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_53dN_pool_qval2_annot, file="output/ChIPseeker/annotation_macs2_H3K4me3_53dN_pool_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_53dN_pool_qval4_annot, file="output/ChIPseeker/annotation_macs2_H3K27me3_53dN_pool_qval4.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_53dN_pool_qval4_annot, file="output/ChIPseeker/annotation_macs2_H3K4me3_53dN_pool_qval4.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K27me3_53dN_pool_qval2_annot_promoterAnd5 = tibble(H3K27me3_53dN_pool_qval2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_53dN_pool_qval2_annot_promoterAnd5 = tibble(H3K4me3_53dN_pool_qval2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_53dN_pool_qval4_annot_promoterAnd5 = tibble(H3K27me3_53dN_pool_qval4_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_53dN_pool_qval4_annot_promoterAnd5 = tibble(H3K4me3_53dN_pool_qval4_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
H3K27me3_53dN_pool_qval2_annot_promoterAnd5_geneSymbol = H3K27me3_53dN_pool_qval2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_53dN_pool_qval2_annot_promoterAnd5_geneSymbol = H3K4me3_53dN_pool_qval2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_53dN_pool_qval4_annot_promoterAnd5_geneSymbol = H3K27me3_53dN_pool_qval4_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_53dN_pool_qval4_annot_promoterAnd5_geneSymbol = H3K4me3_53dN_pool_qval4_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(H3K27me3_53dN_pool_qval2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K27me3_53dN_pool_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_53dN_pool_qval2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K4me3_53dN_pool_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_53dN_pool_qval4_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K27me3_53dN_pool_qval4_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_53dN_pool_qval4_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_H3K4me3_53dN_pool_qval4_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```

