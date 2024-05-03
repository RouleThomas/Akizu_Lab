# Overview

Collaboration with Estaras lab for ChIPseq analysis:
- NR2F2, YAP, EZH2, QSER1 in CPCs and hESCs (hESC-derived cardiac progenitors)


ChIPs and input sequenced separately:
- Chip in Basespace
- inputs in dropbox. [Inputs](https://www.dropbox.com/t/5CJhDwoT2KBPREhG) for CPCs (2 inputs; sample5 = CPC untreated (minus RA) and sample6 = RA treated (plus RA) --> NR2F2 and YAP) and hESCs. and [input](https://www.dropbox.com/t/qUo4hA4bsNk6bK8J) for hESC (untreated, under selfrenewal conditions)

Analysis:
- **CPC**: untreated vs treated with YAP1, TEAD4, NR2F2; with 1 input untreated and 1 input treated; (3 bio rep for all) --> *NOTE: 1 bio rep for YAP and TEAD4 are PE... All other files are SE, so I will treat them as everyon is PE using Read1*
- **hESC** WT vs YAPKO with EZH2, QSER1, DVL2; with 1 input for WT only; (2 bio rep for EZH2 and QSER1. 1 Bio rep for DVL2.)



- *NOTE: data downloaded in local from basespace/dropbox and then imported into the cluster*
- *NOTE: data is single end and paired end; so all treated as PE*


Objective:
- provide bigwig; is there peak


# Data renaming

Let's see all files we have and rename; and confirm with Conchi it is all good (file = `rename_008001.xlsx`)

--> The files from Basespace are easy to rename. The one from the rawData.zip need to be concatenated 1st...

## Basespace files
I created a tab separated file with current (`sample_name.txt`) / new file names (keeping the .fq.gz sufix) and then:

**make sure to convert the `rename_map.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input_raw

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_map.txt
```

--> All good 


## rawData files

--> Combine files from multiple lanes.
----> I add to manually rename CE10 as name was: `20979572019.gz`... WTF!! --> `CE9_S10_L001_R1_001.fastq.gz`

Concatenate fastq discuss [here](https://www.biostars.org/p/317385/): `cat string_L001_sampleID_R1.fastq.gz string_L002_sampleID_R1.fastq.gz  > string_sampleID_R1.fastq.gz`.

--> Several of our samples dispatched in 4 lanes (`L001` to `L004`; **the concatenated are names `L14`; all unique are `L4`**)
----> input sample: `input_raw_Novogene/` output to `input/`

```bash
sbatch scripts/concatenate_CE56.sh # 17683162 ok
sbatch scripts/concatenate_CE78.sh # 17683171 ok
sbatch scripts/concatenate_CE910.sh # 17683182 ok

```

--> Only the Read1 file has been rename; like they were SE!



**make sure to convert the `rename_map2.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input_raw

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_map2.txt
```

--> All good 



# Fastp cleaning

```bash
sbatch scripts/fastp_CPC.sh # 17691945 ok
sbatch scripts/fastp_hESC.sh # 17691956 ok
```


# FastQC


**raw**
```bash
sbatch scripts/fastqc_CPC_raw.sh # 17691968 ok
sbatch scripts/fastqc_hESC_raw.sh # 17691973 ok
```


**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:17691945 scripts/fastqc_CPC_fastp.sh # 17691992 ok
sbatch --dependency=afterany:17691956 scripts/fastqc_hESC_fastp.sh # 17691996 ok
```

--> all good  


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)
--> *NOTE: I removed `--no-mixed --dovetail` and `-U` (instead of `-r`) for the fastq path as PE options* 


```bash
conda activate bowtie2

sbatch --dependency=afterany:17691945 scripts/bowtie2_CPC.sh # 17692156 xxx
sbatch scripts/bowtie2_hESC.sh # 17692160 ok

```



--> Looks good; around 70% uniq aligned reads (95% total) for hESC XXX and CPC XXX




## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-17692160.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_17692160.txt

XXXX TO RUN WHEN 17692156 FINISH XXX
for file in slurm-17692156.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_17692156.txt

```

Add these values to `/home/roulet/008_ChIPseq_YAP_Conchi/001__ChIPseq_V1/samples_008001.xlsx`\
Then in R; see `/home/roulet/008_ChIPseq_YAP_Conchi/ChIPseq_YAP.R`.

--> Overall > XXX % input reads as been uniquely mapped to the genome (XXX % non uniq)





## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.

```bash
conda activate bowtie2

sbatch scripts/samtools_unique_hESC.sh # 17725506 ok
sbatch --dependency=afterany:17692156 scripts/samtools_unique_CPC.sh # 17725583 xxx

```

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch scripts/samtools_MG1655_unique_1.sh # 9162457
sbatch scripts/samtools_MG1655_unique_2.sh # 9162461
sbatch scripts/samtools_MG1655_unique_3.sh # 9162467
```

--> More information on this step in the `005__CutRun` labnote

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


# Generate bigwig coverage files
## Raw bigwig

Let's use *ChIPQC* (as [greenscreen](https://github.com/sklasfeld/GreenscreenProject/blob/main/TUTORIAL.pdf) paper), to estimate fragment size and provide it to each respective sample. --> **NEEDED as we are in SE**


Let's install it in a new conda env `deseq2V3` in R `BiocManager::install("ChIPQC")`
- nano `scripts/ChIPQC.R` from [greenscreen github](https://github.com/sklasfeld/GreenscreenProject/blob/main/scripts/ChIPQC.R)
- create csv table `meta/sampleSheet.csv` describing each sample; from [greenscreen github](https://github.com/sklasfeld/GreenscreenProject/blob/main/meta/noMaskReads_Inputs_sampleSheet.csv)

Run in R; followed this [workshop](https://nbisweden.github.io/workshop-archive/workshop-ChIP-seq/2018-11-07/labs/lab-chipqc.html)


```bash
conda activate deseq2V3
```

```R
library("DiffBind")
library("ChIPQC")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")


#	reading in the sample information (metadata)
samples = read.csv("meta/sampleSheet_hESC.csv", sep="\t")
samples = read.csv("meta/sampleSheet_CPC.csv", sep="\t")

#	inspecting the metadata
samples

#	creating an object containing data
res=dba(sampleSheet=samples, config=data.frame(RunParallel=FALSE))

# inspecting the object
res

#	performing quality control
resqc = ChIPQC(res,annotation="hg38", config=data.frame(RunParallel=TRUE))

#	creating the quality control report in html format
ChIPQCreport(resqc)

```

--> A `ChIPQCreport` folder is created in current wd; I moved it to `ouput`


Then generate bigwig with the corresponding fragment size for each sample:



Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads. **HERE single end so need to estimate fragment size!! Extend to `FragL-ReadL`**


```bash
conda activate deeptools

# bigwig with extendReads 50 (Default ~250bp fragment length  150bp reads +100 bp)
sbatch --dependency=afterany:17725583 scripts/bamtobigwig_unique_extendReads100_CPC.sh # 17762131 ok
sbatch scripts/bamtobigwig_unique_extendReads100_hESC.sh # 17762114 ok

# bigwig with extendReads from CHIPQC
sbatch scripts/bamtobigwig_unique_extendReads_hESC.sh # 17770516 ok
sbatch scripts/bamtobigwig_unique_extendReads_CPC.sh # 17775530 ok


# bigwig with extendReads from CHIPQC and RPGC normalized (seq depth comparison)
sbatch scripts/bamtobigwig_unique_extendReads_RPGC_hESC.sh # 17786752 ok
sbatch scripts/bamtobigwig_unique_extendReads_RPGC_CPC.sh # 17786754 ok

```


--> All work in hESC

--> CPC; YAP1 and TEAD4 only R3 work (thats ok it is their control sample! And they already have YAP1 sequenced), NR2F2 R1 and R2 work (not R3), 

--> *RPGC* bigwig make replicate very heterogeneous... = *BAD*; Use instead the **extendReads** = **GOOD**


## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_extendReads_hESC.sh # 17863824 ok
sbatch scripts/multiBigwigSummary_CPC.sh # 17863829 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_extendReads_hESC.npz \
    --transpose \
    --ntop 0 \
    --labels hESC_WT_DVL2_R1 hESC_YAPKO_DVL2_R1 hESC_WT_EZH2_R1 hESC_YAPKO_EZH2_R1 hESC_WT_EZH2_R2 hESC_YAPKO_EZH2_R2 hESC_WT_QSER1_R1 hESC_YAPKO_QSER1_R1 hESC_WT_QSER1_R2 hESC_YAPKO_QSER1_R2 hESC_WT_input_R1 \
    -o output/bigwig/multiBigwigSummary_extendReads_hESC_plotPCA.pdf
plotPCA -in output/bigwig/multiBigwigSummary_extendReads_CPC.npz \
    --transpose \
    --ntop 0 \
    --labels CPC_untreated_YAP1_R1 CPC_RA_YAP1_R1 CPC_untreated_YAP1_R2 CPC_RA_YAP1_R2 CPC_untreated_TEAD4_R1 CPC_RA_TEAD4_R1 CPC_untreated_TEAD4_R2 CPC_RA_TEAD4_R2 CPC_untreated_NR2F2_R1 CPC_RA_NR2F2_R1 CPC_untreated_NR2F2_R2 CPC_RA_NR2F2_R2 CPC_RA_NR2F2_R3 CPC_untreated_NR2F2_R3 CPC_untreated_input_R3 CPC_RA_input_R3 CPC_untreated_YAP1_R3 CPC_RA_YAP1_R3 CPC_untreated_TEAD4_R3 CPC_RA_TEAD4_R3 \
    -o output/bigwig/multiBigwigSummary_extendReads_CPC_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_extendReads_hESC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels hESC_WT_DVL2_R1 hESC_YAPKO_DVL2_R1 hESC_WT_EZH2_R1 hESC_YAPKO_EZH2_R1 hESC_WT_EZH2_R2 hESC_YAPKO_EZH2_R2 hESC_WT_QSER1_R1 hESC_YAPKO_QSER1_R1 hESC_WT_QSER1_R2 hESC_YAPKO_QSER1_R2 hESC_WT_input_R1 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_extendReads_hESC_heatmap.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_extendReads_CPC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels CPC_untreated_YAP1_R1 CPC_RA_YAP1_R1 CPC_untreated_YAP1_R2 CPC_RA_YAP1_R2 CPC_untreated_TEAD4_R1 CPC_RA_TEAD4_R1 CPC_untreated_TEAD4_R2 CPC_RA_TEAD4_R2 CPC_untreated_NR2F2_R1 CPC_RA_NR2F2_R1 CPC_untreated_NR2F2_R2 CPC_RA_NR2F2_R2 CPC_RA_NR2F2_R3 CPC_untreated_NR2F2_R3 CPC_untreated_input_R3 CPC_RA_input_R3 CPC_untreated_YAP1_R3 CPC_RA_YAP1_R3 CPC_untreated_TEAD4_R3 CPC_RA_TEAD4_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_extendReads_CPC_heatmap.pdf

```

--> two big groups: H3K27me3 IP versus the other
----> Seems only H3K27me3 IP has worked here




# MACS2 peak calling on bam unique

--> input samples used as control

--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad` and `narrow`**

- *NOTE: as SE; I removed `-f BAMPE` option*

```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad_hESC.sh # 18364489 ok
sbatch scripts/macs2_broad_CPC.sh # 18365107 ok

sbatch scripts/macs2_narrow_hESC.sh # 18365200 ok
sbatch scripts/macs2_narrow_CPC.sh # 18365397 ok

# pool
sbatch scripts/macs2_broad_hESC_pool.sh # 18482872 ok
sbatch scripts/macs2_broad_CPC_pool.sh # 18482945 ok

sbatch scripts/macs2_narrow_hESC_pool.sh # 18482993 ok
sbatch scripts/macs2_narrow_CPC_pool.sh # 18483025 ok



```


--> For pool, only done on samples that worked: 
- For the YAP1, Rep3 work (but not Rep1 and Rep2)
- For TEAD4, Rep3 work (but not Rep1 and Rep2)
- For NR2F2, Rep1 and Rep2 work (but not Rep3)
- For EZH2, Rep1 and Rep2 work
- For QSER1, Rep1 and Rep2 work
- For DVL2, Rep1 work (only 1 Rep for this sample)





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
- CPC_YAP1, untreated and RA (Rep3 only): qval 1.30103 #9526  and #19659
- CPC_TEAD4, untreated and RA (Rep3 only): qval 1.30103 (more false positive maybe 2.3 if too many peaks) #12478  and #24804
- CPC_NR2F2, untreated and RA, (use pool; Rep1 and Rep2): qval 2.30103 #33501  and #48276
- hESC_EZH2, WT and KO (use pool; Rep1 and Rep2): qval 1.30103 #4288  and #4599
- hESC_QSER1 WT and KO (use pool; Rep1 and Rep2): qval 1.30103 #14170  and #14978
- hESC_DVL2 WT and KO (Rep1 only): qval 1.30103 (YAPKO work badly!) #232 and #2


--> broad identified more peaks! So let's use broad (broad simply combine [several narrow peaks into 1](https://www.biostars.org/p/245407/))

# ChIPseeker - macs2


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
## CPC qval optimal
CPC_untreated_NR2F2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/CPC_untreated_NR2F2_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
CPC_RA_NR2F2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/CPC_RA_NR2F2_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

CPC_untreated_TEAD4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/CPC_untreated_TEAD4_R3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
CPC_RA_TEAD4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/CPC_RA_TEAD4_R3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)    
    
CPC_untreated_YAP1 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/CPC_untreated_YAP1_R3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
CPC_RA_YAP1 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/CPC_RA_YAP1_R3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)    

## CPC qval3
CPC_untreated_NR2F2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/CPC_untreated_NR2F2_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
CPC_RA_NR2F2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/CPC_RA_NR2F2_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

CPC_untreated_TEAD4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/CPC_untreated_TEAD4_R3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
CPC_RA_TEAD4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/CPC_RA_TEAD4_R3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)    
    
CPC_untreated_YAP1 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/CPC_untreated_YAP1_R3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
CPC_RA_YAP1 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval3/CPC_RA_YAP1_R3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)    

## Tidy peaks #-->> Re-Run from here with different qvalue!!
CPC_untreated_NR2F2_gr = makeGRangesFromDataFrame(CPC_untreated_NR2F2,keep.extra.columns=TRUE)
CPC_RA_NR2F2_gr = makeGRangesFromDataFrame(CPC_RA_NR2F2,keep.extra.columns=TRUE)
CPC_untreated_TEAD4_gr = makeGRangesFromDataFrame(CPC_untreated_TEAD4,keep.extra.columns=TRUE)
CPC_RA_TEAD4_gr = makeGRangesFromDataFrame(CPC_RA_TEAD4,keep.extra.columns=TRUE)
CPC_untreated_YAP1_gr = makeGRangesFromDataFrame(CPC_untreated_YAP1,keep.extra.columns=TRUE)
CPC_RA_YAP1_gr = makeGRangesFromDataFrame(CPC_RA_YAP1,keep.extra.columns=TRUE)

gr_list <- list(CPC_untreated_NR2F2=CPC_untreated_NR2F2_gr, CPC_RA_NR2F2=CPC_RA_NR2F2_gr, CPC_untreated_TEAD4=CPC_untreated_TEAD4_gr,  CPC_RA_TEAD4=CPC_RA_TEAD4_gr, CPC_untreated_YAP1= CPC_untreated_YAP1_gr,CPC_RA_YAP1= CPC_RA_YAP1_gr)



## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here

                       
### Barplot
pdf("output/ChIPseeker/annotation_barplot_CPC.pdf", width=14, height=5)
pdf("output/ChIPseeker/annotation_barplot_CPC_qval3.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()



## Get annotation data frame
CPC_untreated_NR2F2_annot <- as.data.frame(peakAnnoList[["CPC_untreated_NR2F2"]]@anno)
CPC_RA_NR2F2_annot <- as.data.frame(peakAnnoList[["CPC_RA_NR2F2"]]@anno)
CPC_untreated_TEAD4_annot <- as.data.frame(peakAnnoList[["CPC_untreated_TEAD4"]]@anno)
CPC_RA_TEAD4_annot <- as.data.frame(peakAnnoList[["CPC_RA_TEAD4"]]@anno)
CPC_untreated_YAP1_annot <- as.data.frame(peakAnnoList[["CPC_untreated_YAP1"]]@anno)
CPC_RA_YAP1_annot <- as.data.frame(peakAnnoList[["CPC_RA_YAP1"]]@anno)

## Convert entrez gene IDs to gene symbols
CPC_untreated_NR2F2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = CPC_untreated_NR2F2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
CPC_untreated_NR2F2_annot$gene <- mapIds(org.Hs.eg.db, keys = CPC_untreated_NR2F2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
CPC_RA_NR2F2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = CPC_RA_NR2F2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
CPC_RA_NR2F2_annot$gene <- mapIds(org.Hs.eg.db, keys = CPC_RA_NR2F2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
CPC_untreated_TEAD4_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = CPC_untreated_TEAD4_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
CPC_untreated_TEAD4_annot$gene <- mapIds(org.Hs.eg.db, keys = CPC_untreated_TEAD4_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
CPC_RA_TEAD4_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = CPC_RA_TEAD4_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
CPC_RA_TEAD4_annot$gene <- mapIds(org.Hs.eg.db, keys = CPC_RA_TEAD4_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
CPC_untreated_YAP1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = CPC_untreated_YAP1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
CPC_untreated_YAP1_annot$gene <- mapIds(org.Hs.eg.db, keys = CPC_untreated_YAP1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
CPC_RA_YAP1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = CPC_RA_YAP1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
CPC_RA_YAP1_annot$gene <- mapIds(org.Hs.eg.db, keys = CPC_RA_YAP1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(CPC_untreated_NR2F2_annot, file="output/ChIPseeker/annotation_macs2_CPC_untreated_NR2F2_annot_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(CPC_RA_NR2F2_annot, file="output/ChIPseeker/annotation_macs2_CPC_RA_NR2F2_annot_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(CPC_untreated_TEAD4_annot, file="output/ChIPseeker/annotation_macs2_CPC_untreated_TEAD4_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(CPC_RA_TEAD4_annot, file="output/ChIPseeker/annotation_macs2_CPC_RA_TEAD4_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(CPC_untreated_YAP1_annot, file="output/ChIPseeker/annotation_macs2_CPC_untreated_YAP1_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(CPC_RA_YAP1_annot, file="output/ChIPseeker/annotation_macs2_CPC_RA_YAP1_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
CPC_untreated_NR2F2_annot_promoterAnd5 = tibble(CPC_untreated_NR2F2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
CPC_RA_NR2F2_annot_promoterAnd5 = tibble(CPC_RA_NR2F2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
CPC_untreated_TEAD4_annot_promoterAnd5 = tibble(CPC_untreated_TEAD4_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
CPC_RA_TEAD4_annot_promoterAnd5 = tibble(CPC_RA_TEAD4_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))    
CPC_untreated_YAP1_annot_promoterAnd5 = tibble(CPC_untreated_YAP1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
CPC_RA_YAP1_annot_promoterAnd5 = tibble(CPC_RA_YAP1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) 

### Save output gene lists
CPC_untreated_NR2F2_annot_promoterAnd5_geneSymbol = CPC_untreated_NR2F2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
CPC_RA_NR2F2_annot_promoterAnd5_geneSymbol = CPC_RA_NR2F2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
CPC_untreated_TEAD4_annot_promoterAnd5_geneSymbol = CPC_untreated_TEAD4_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
CPC_RA_TEAD4_annot_promoterAnd5_geneSymbol = CPC_RA_TEAD4_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()   
CPC_untreated_YAP1_annot_promoterAnd5_geneSymbol = CPC_untreated_YAP1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
CPC_RA_YAP1_annot_promoterAnd5_geneSymbol = CPC_RA_YAP1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique() 


write.table(CPC_untreated_NR2F2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_untreated_NR2F2_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(CPC_RA_NR2F2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_RA_NR2F2_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(CPC_untreated_TEAD4_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_untreated_TEAD4_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(CPC_RA_TEAD4_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_RA_TEAD4_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    

write.table(CPC_untreated_YAP1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_untreated_YAP1_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(CPC_RA_YAP1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_RA_YAP1_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    





## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
CPC_untreated_NR2F2_annot_noIntergenic = tibble(CPC_untreated_NR2F2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
CPC_RA_NR2F2_annot_noIntergenic = tibble(CPC_RA_NR2F2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
CPC_untreated_TEAD4_annot_noIntergenic = tibble(CPC_untreated_TEAD4_annot) %>%
    filter(annotation != c("Distal Intergenic"))
CPC_RA_TEAD4_annot_noIntergenic = tibble(CPC_RA_TEAD4_annot) %>%
    filter(annotation != c("Distal Intergenic"))
CPC_untreated_YAP1_annot_noIntergenic = tibble(CPC_untreated_YAP1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
CPC_RA_YAP1_annot_noIntergenic = tibble(CPC_RA_YAP1_annot) %>%
    filter(annotation != c("Distal Intergenic"))

### Save output gene lists
CPC_untreated_NR2F2_annot_noIntergenic_geneSymbol = CPC_untreated_NR2F2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
CPC_RA_NR2F2_annot_noIntergenic_geneSymbol = CPC_RA_NR2F2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
CPC_untreated_TEAD4_annot_noIntergenic_geneSymbol = CPC_untreated_TEAD4_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
CPC_RA_TEAD4_annot_noIntergenic_geneSymbol = CPC_RA_TEAD4_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()   
CPC_untreated_YAP1_annot_noIntergenic_geneSymbol = CPC_untreated_YAP1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
CPC_RA_YAP1_annot_noIntergenic_geneSymbol = CPC_RA_YAP1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique() 


write.table(CPC_untreated_NR2F2_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_untreated_NR2F2_qval2.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(CPC_RA_NR2F2_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_RA_NR2F2_qval2.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(CPC_untreated_TEAD4_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_untreated_TEAD4_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(CPC_RA_TEAD4_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_RA_TEAD4_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    

write.table(CPC_untreated_YAP1_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_untreated_YAP1_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(CPC_RA_YAP1_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_CPC_RA_YAP1_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    








## hESC
hESC_WT_EZH2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_EZH2_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_EZH2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/hESC_YAPKO_EZH2_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

hESC_WT_QSER1 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_QSER1_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_QSER1 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/hESC_YAPKO_QSER1_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

hESC_WT_DVL2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_DVL2_R1_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_DVL2 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/hESC_YAPKO_DVL2_R1_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 



## Tidy peaks #-->> Re-Run from here with different qvalue!!
hESC_WT_EZH2_gr = makeGRangesFromDataFrame(hESC_WT_EZH2,keep.extra.columns=TRUE)
hESC_YAPKO_EZH2_gr = makeGRangesFromDataFrame(hESC_YAPKO_EZH2,keep.extra.columns=TRUE)
hESC_WT_QSER1_gr = makeGRangesFromDataFrame(hESC_WT_QSER1,keep.extra.columns=TRUE)
hESC_YAPKO_QSER1_gr = makeGRangesFromDataFrame(hESC_YAPKO_QSER1,keep.extra.columns=TRUE)
hESC_WT_DVL2_gr = makeGRangesFromDataFrame(hESC_WT_DVL2,keep.extra.columns=TRUE)
hESC_YAPKO_DVL2_gr = makeGRangesFromDataFrame(hESC_YAPKO_DVL2,keep.extra.columns=TRUE)

gr_list <- list(hESC_WT_EZH2=hESC_WT_EZH2_gr, hESC_YAPKO_EZH2=hESC_YAPKO_EZH2_gr, hESC_WT_QSER1=hESC_WT_QSER1_gr,  hESC_YAPKO_QSER1=hESC_YAPKO_QSER1_gr, hESC_WT_DVL2= hESC_WT_DVL2_gr) # hESC_YAPKO_DVL2_gr removed as only 2 peak



## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here

                       
### Barplot
pdf("output/ChIPseeker/annotation_barplot_hESC.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()



## Get annotation data frame
hESC_WT_EZH2_annot <- as.data.frame(peakAnnoList[["hESC_WT_EZH2"]]@anno)
hESC_YAPKO_EZH2_annot <- as.data.frame(peakAnnoList[["hESC_YAPKO_EZH2"]]@anno)
hESC_WT_QSER1_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1"]]@anno)
hESC_YAPKO_QSER1_annot <- as.data.frame(peakAnnoList[["hESC_YAPKO_QSER1"]]@anno)
hESC_WT_DVL2_annot <- as.data.frame(peakAnnoList[["hESC_WT_DVL2"]]@anno)

## Convert entrez gene IDs to gene symbols
hESC_WT_EZH2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_EZH2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_YAPKO_EZH2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_EZH2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_YAPKO_EZH2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_EZH2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_YAPKO_QSER1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_QSER1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_YAPKO_QSER1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_QSER1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_DVL2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_DVL2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_DVL2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_DVL2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(hESC_WT_EZH2_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_YAPKO_EZH2_annot, file="output/ChIPseeker/annotation_macs2_hESC_YAPKO_EZH2_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_QSER1_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_YAPKO_QSER1_annot, file="output/ChIPseeker/annotation_macs2_hESC_YAPKO_QSER1_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_DVL2_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_DVL2_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_EZH2_annot_promoterAnd5 = tibble(hESC_WT_EZH2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_YAPKO_EZH2_annot_promoterAnd5 = tibble(hESC_YAPKO_EZH2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1_annot_promoterAnd5 = tibble(hESC_WT_QSER1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_YAPKO_QSER1_annot_promoterAnd5 = tibble(hESC_YAPKO_QSER1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))    
hESC_WT_DVL2_annot_promoterAnd5 = tibble(hESC_WT_DVL2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
hESC_WT_EZH2_annot_promoterAnd5_geneSymbol = hESC_WT_EZH2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_annot_promoterAnd5_geneSymbol = hESC_YAPKO_EZH2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_QSER1_annot_promoterAnd5_geneSymbol = hESC_YAPKO_QSER1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()   
hESC_WT_DVL2_annot_promoterAnd5_geneSymbol = hESC_WT_DVL2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()



write.table(hESC_WT_EZH2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_YAPKO_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_QSER1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_YAPKO_QSER1_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    
write.table(hESC_WT_DVL2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_DVL2_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
            



## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_EZH2_annot_noIntergenic = tibble(hESC_WT_EZH2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_YAPKO_EZH2_annot_noIntergenic = tibble(hESC_YAPKO_EZH2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1_annot_noIntergenic = tibble(hESC_WT_QSER1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_YAPKO_QSER1_annot_noIntergenic = tibble(hESC_YAPKO_QSER1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_DVL2_annot_noIntergenic = tibble(hESC_WT_DVL2_annot) %>%
    filter(annotation != c("Distal Intergenic"))


### Save output gene lists
hESC_WT_EZH2_annot_noIntergenic_geneSymbol = hESC_WT_EZH2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_annot_noIntergenic_geneSymbol = hESC_YAPKO_EZH2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1_annot_noIntergenic_geneSymbol = hESC_WT_QSER1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_QSER1_annot_noIntergenic_geneSymbol = hESC_YAPKO_QSER1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()   
hESC_WT_DVL2_annot_noIntergenic_geneSymbol = hESC_WT_DVL2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()



write.table(hESC_WT_EZH2_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_YAPKO_EZH2_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_QSER1_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_YAPKO_QSER1_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    
write.table(hESC_WT_DVL2_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_DVL2_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
            



```

--> CPC is a bit weird; as untreated vs RA as different binding profile; not clear whether it is technical issue or not. But higher qval (1.3 to 3) give same results

--> hESC is all good.







# THOR diff peaks

Let's use THOR, notably to have IGG scaled bigwig...!

Comparison to do:
- CPC = untreated vs RA; YAP1/TEAD4 (use R3), NR2F2 (use R1 and R2)
- hESC = WT vs YAPKO; DVL2 (only 1 bio rep R1), EZH2, QSER1 (use R1 and R2)


--> Configs file created manually as:
- `output/THOR/CPC_YAP1_untreatedvsRA.config`, `output/THOR/CPC_TEAD4_untreatedvsRA.config`, `output/THOR/CPC_NR2F2_untreatedvsRA.config` 
- `output/THOR/hESC_DVL2_WTvsYAPKO.config`, `output/THOR/hESC_EZH2_WTvsYAPKO.config`, `output/THOR/hESC_QSER1_WTvsYAPKO.config` 


*NOTE: Default TMM normalization applied, as no Spike in*

## Run THOR

*THOR is very buggy to make it work I need to temporaly change where to look for libraries lol.. So cannot use nano anymore for example...*

*Follow these parameters: `WTvsHET_unique_Keepdup` (perform best in previous CutRun)*

```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge

# AB per AB (Default TMM norm)
sbatch scripts/THOR_CPC_YAP1_untreatedvsRA.sh # 17867547 ok
sbatch scripts/THOR_CPC_TEAD4_untreatedvsRA.sh # 17867550 ok
sbatch scripts/THOR_CPC_NR2F2_untreatedvsRA.sh # 17867553 ok

sbatch scripts/THOR_hESC_DVL2_WTvsYAPKO.sh # 17867580 fail; 17869560 ok
sbatch scripts/THOR_hESC_EZH2_WTvsYAPKO.sh # 17867582 ok
sbatch scripts/THOR_hESC_QSER1_WTvsYAPKO.sh # 17867584 ok
```
--> NODAL signif lost EZH2 in KO !!! YEAH!



Generate median tracks:
```bash
conda activate BedToBigwig

sbatch scripts/bigwigmerge_THOR_CPC.sh # 17869306 ok
sbatch --dependency=afterany:17869560 scripts/bigwigmerge_THOR_hESC.sh # 17869572 xxx

```


--> THOR bigwig tracks looks good! NODAL loose EZH2 in YAPKO



## Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks


```R

# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# WTvsKO
diffpeaks <- read_tsv("output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_50dN_H3K27me3_WTvsKO/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("50dN_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_50dN_H3K27me3_WTvsKO/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO_qval25") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 100) %>%
  write_tsv("output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval100.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())




```

- *NOTE: FC positive = less in KO; negative = more in KO*

**Optimal qvalue:**
--> *H3K27me3*; qval 10 looks great!

--> In agreement with `003__CutRun`; in KO overall same number of gain and lost regions; and in KOEF1aEZH1 much more gain of H3K27me3 (act like the HET)

Isolate positive and negative THOR peaks to display deepTool plots

```bash
# positive negative peaks
## qval 10
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10_negative.bed


```




# deepTool plots


Venn diagram of peak-genes TED4, QSER1 and EZH2 has been generated.

- Check whether QSER1 flank EZH2 (confirm [Dixon2021 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8185639/))
- Check TED4 vs QSER1 (and EZH2)

--> use THOR bigwig for QSER1 and EZH2 and use raw unique bigwig for TEAD4 (from `008003`)


XXXX BELOW !!! 


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



