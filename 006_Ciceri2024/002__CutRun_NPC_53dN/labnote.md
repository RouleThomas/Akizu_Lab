# Project and goals 

Re-analysis of CutRun dataset from (Ciceri et al)[https://www.nature.com/articles/s41586-023-06984-8].

- Identify histone -mod bound genes in WT and compare with our data


# Download data


- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject PRJNA803355 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 14679607 ok

```



## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename_002.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_002.txt
```

--> All good 





# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 14681036 ok
```



# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:14681036 scripts/bowtie2_NPC.sh # 14681230 ok
sbatch --dependency=afterany:14681036 scripts/bowtie2_53dN.sh # 14681246 ok
```

--> Looks good; overall ~30-80% uniquely aligned reads
----> Seems less uniquel mapped reads than us but they sequence FAR more depth (~20m reads vs 5 for us)


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-14681230.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_14681230.txt

for file in slurm-14681246.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_14681246.txt


```

Add these values to `/home/roulet/006_Ciceri2024/002__CutRun_NPC_53dN/samples_002.xlsx`\
Then in R; see `/home/roulet/006_Ciceri2024/006_Ciceri2024.R`.

--> Overall >60% input reads as been uniquely mapped to the genome (90% non uniq)



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.


```bash
conda activate bowtie2

sbatch --dependency=afterany:14681230 scripts/samtools_unique_NPC.sh # 14681286 NODE failure
sbatch --dependency=afterany:14681246 scripts/samtools_unique_53dN.sh # 14681396 NODE failure


sbatch scripts/samtools_unique_NPC_1.sh # 15101660 fail; 15116487 ok
sbatch scripts/samtools_unique_NPC_2.sh # 15275393 ok
sbatch scripts/samtools_unique_53dN_1.sh # 15101878 ok

```


# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:14681286 scripts/bamtobigwig_unique_NPC.sh # 14681462 node failure
sbatch --dependency=afterany:14681396 scripts/bamtobigwig_unique_53dN.sh # 14681480 node failure; 

sbatch --dependency=afterany:15101878 scripts/bamtobigwig_unique_53dN_1.sh # 15102289 ok
sbatch scripts/bamtobigwig_unique_NPC_1.sh # 15102387 fail; 15116725 fail; 15137560 ok
sbatch --dependency=afterany:15275393 scripts/bamtobigwig_unique_NPC_2.sh # 15275543 ok



```


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



# deepTool plots


On all genes


```bash
conda activate deeptools

# All genes all histone marks
## NPC
sbatch scripts/matrix_TSS_10kb_NPC_raw_allGenes.sh # 15402417 ok
sbatch scripts/matrix_TSS_5kb_NPC_H3K27ac_raw_allGenes.sh # 15402511 ok
sbatch scripts/matrix_TSS_10kb_NPC_H3K27me3_H3K4me3_raw_allGenes.sh # 15402602 ok


## 53dN
sbatch scripts/matrix_TSS_10kb_53dN_raw_allGenes.sh # 15402437 ok
sbatch scripts/matrix_TSS_5kb_53dN_H3K27ac_raw_allGenes.sh # 15402531 ok
sbatch scripts/matrix_TSS_10kb_53dN_H3K27me3_H3K4me3_raw_allGenes.sh # 15402648 ok


# Akizu and Ciceri H3K27me3, H3K4me3, IGG
sbatch scripts/matrix_TSS_10kb_H3K27me3_H3K4me3_CutRun001008_Ciceri_raw_allGenes.sh # 15432941 ok


```


--> signal is very poor for H3K27ac in NPC as compared to 53dN, notably for R1, almost like IGG...








# Pilot grant 20240204

Patient1 has EZH1 GOF mutation + KMT2A (H3K4 methyltransferase) KO
Hypothesis: Gene that gain H3K27me3 in patient1 are because loss of H3K4m3 so more room for H3K27me3 expansion

Isolate bivalent genes in WT: 
- Macs2 peak calling in the 2 pool replicate (qval 2.3 and 4) 
- ChIPseeker gene peak assignment
- Filter –in promoter/TSS peaks `annotation_macs2_H3K*me3_WT_pool_qval2.30103_promoterAnd5_geneSymbol`
- Venn diagram of peak enriched genes


**Metrics qval 2.3**:
- WT_H3K27me3:
    - 10,259 peaks
    - 4,055 genes
- WT_H3K4me3
    - 18,711 peaks
    - 12,726 genes
--> 2,560 bivalent genes

**Metrics qval 4**:
- WT_H3K27me3:
    - 958 peaks
    - 747 genes
- WT_H3K4me3
    - 7,248 peaks
    - 6,313 genes
--> 379 bivalent genes


## Generate deeptool plots for the H3K4me3, H3K27me3 and bivalent genes


```bash
# Generate gtf file from gene list:

### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval4.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval4_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval4.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval4_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval4.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval4_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval2.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval2_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval2.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval2_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval2.txt > output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval2_as_gtf_geneSymbol.txt




## Filter the gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval4_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval4.gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval4_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval4.gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval4_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval4.gtf

grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval2.gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval2.gtf
grep -Ff output/ChIPseeker/Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval2.gtf

# deeptool plots
## qval 2.3
sbatch scripts/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval2.sh # 17512685 ok

## qval 4
sbatch scripts/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4.sh # 17512702 ok
```


--> The Venn overlap gene filtering worked great; H3K4me3/H3K27me3/bivalent genes are clearly identified at q2.3. But some genes show no signal, even when decreasing the zMax scale; let's try more stringeant qvalues
---->  With more stringeat qvalues less good

--> Optimal qvalue is 2.3 




## enrichR functional analysis

On the **bivalent** - H3K4me3 and H3K27me3 (53dN) at qval 2.3 with the **genes that gain H3K27me3 in GOF/HET** (`003__CutRun`)




--> Below code modified to show only 1 set of genes (no up and down, only 1 set)


```R
# packages
library("tidyverse")
library("enrichR")


# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 

### GeneSymbol list of DEGs per tissue
output/ChIPseeker/Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt



# IF starting with geneSymbol

## Read and preprocess data for DEGs genes
gene_names_up <- read.csv("output/ChIPseeker/Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
up$type <- "up"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 50) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)

# Combine the two dataframes
gos <- up
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
new_order <- up_pathways
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
## gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.pdf", width=20, height=13)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 12, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Purple")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("KEGG_2021_Human") # 

### GeneSymbol list of DEGs per tissue
output/ChIPseeker/Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt


# IF starting with geneSymbol

## Read and preprocess data for DEGs genes
gene_names_up <- read.csv("output/ChIPseeker/Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)




# Extracting KEGG data and assigning types
up <- eup$KEGG_2021_Human
up$type <- "up"

# Get top enriched terms and sort by Combined.Score 
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 50) ##  Adjust if you don't want the top 5



# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)

# Combine the two dataframes
gos <- up
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
new_order <- up_pathways
gos$Term <- factor(gos$Term, levels = new_order)


# extract the top 5 rows (p adj ordered)
## gos <- head(gos, n = 5)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_KEGG_2021_Human_Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.pdf", width=20, height=12)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.8) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 12, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="Sky Blue", "up"="Purple")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 30)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_Venn_overlap_WTbivalent_003GainHETTHORq15__bivalentOnly.txt", sep="\t", row.names=FALSE, quote=FALSE)

```















