# Overview

Collaboration with Estaras lab for ChIPseq analysis:
- YAP1 and TEAD4 in pluripotency

Downlaod from this [link](https://www.dropbox.com/t/yxGparZoHZ9ynV5k) --> seems no input?? Use hESC input from `008001`




Analysis:
- compare with `008001` data



- *NOTE: data downloaded in local from basespace/dropbox and then imported into the cluster*
- *NOTE: data is single end*


Objective:
- integrate with hESC 008001, check overlap EZH2, QSER1 and TEAD4 (eg. the non-overlapping profile), QSER1 prevent histone methylation (work together with TET1)


# Data renaming

Rename manually as only 2 file
TEAD4_hESC_S33_L003_R1_001.fastq.gz = hESC_WT_TEAD4_R1.fq.gz
YAP_hESC_S30_L003_R1_001.fastq.gz = hESC_WT_YAP1_R1.fq.gz

--> All good 

--> input samples copy from `008001`



# Fastp cleaning

```bash
sbatch scripts/fastp_raw.sh # 18442302 ok
```



# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)
--> *NOTE: I removed `--no-mixed --dovetail` and `-U` (instead of `-r`) for the fastq path as PE options* 


```bash
conda activate bowtie2

sbatch --dependency=afterany:18442302 scripts/bowtie2_raw.sh # 18442369 ok
```



--> Looks good; around 75% uniq aligned reads (95% total) for TEAD4 and YAP1




## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-18442369.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_18442369.txt

```

Add these values to `/home/roulet/008_ChIPseq_YAP_Conchi/003__ChIPseq_pluripotency/samples_008003.xlsx`\
Then in R; see `/home/roulet/008_ChIPseq_YAP_Conchi/ChIPseq_YAP.R`.

--> Overall > XXX % input reads as been uniquely mapped to the genome (XXX % non uniq)




## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.

```bash
conda activate bowtie2

sbatch --dependency=afterany:18442369 scripts/samtools_unique_raw.sh # 18442478 ok

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

--> A `ChIPQCreport` folder is created in current wd; I moved it to `output`


Then generate bigwig with the corresponding fragment size for each sample:



Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads. **HERE single end so need to estimate fragment size!! Extend to `FragL-ReadL`**


```bash
conda activate deeptools

# bigwig with extendReads from CHIPQC
sbatch scripts/bamtobigwig_unique_extendReads_raw.sh # 18481263 ok
```





# MACS2 peak calling on bam unique

--> input samples used as control

--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad` and `narrow`**

- *NOTE: as SE; I removed `-f BAMPE` option*

```bash
conda activate macs2

sbatch scripts/macs2_broad_raw.sh # 18502057 ok

sbatch scripts/macs2_narrow_raw.sh # 18502059 ok
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
- hESC_YAP1, WT (Rep1 only): qval 1.30103 #1807
- hESC_TEAD4, WT (Rep1 only): qval 1.30103 #3998



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
## hESC qval optimal
hESC_WT_YAP1 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/hESC_WT_YAP1_R1_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_TEAD4 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval2.30103/hESC_WT_TEAD4_R1_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 



## Tidy peaks #-->> Re-Run from here with different qvalue!!
hESC_WT_YAP1_gr = makeGRangesFromDataFrame(hESC_WT_YAP1,keep.extra.columns=TRUE)
hESC_WT_TEAD4_gr = makeGRangesFromDataFrame(hESC_WT_TEAD4,keep.extra.columns=TRUE)

gr_list <- list(hESC_WT_YAP1=hESC_WT_YAP1_gr, hESC_WT_TEAD4=hESC_WT_TEAD4_gr)



## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here

                       
### Barplot
pdf("output/ChIPseeker/annotation_barplot_hESC.pdf", width=14, height=2)
plotAnnoBar(peakAnnoList)
dev.off()



## Get annotation data frame
hESC_WT_YAP1_annot <- as.data.frame(peakAnnoList[["hESC_WT_YAP1"]]@anno)
hESC_WT_TEAD4_annot <- as.data.frame(peakAnnoList[["hESC_WT_TEAD4"]]@anno)

## Convert entrez gene IDs to gene symbols
hESC_WT_YAP1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_YAP1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_YAP1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_YAP1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_TEAD4_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_TEAD4_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_TEAD4_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_TEAD4_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(hESC_WT_YAP1_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(hESC_WT_TEAD4_annot, file="output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_YAP1_annot_promoterAnd5 = tibble(hESC_WT_YAP1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_TEAD4_annot_promoterAnd5 = tibble(hESC_WT_TEAD4_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
    
### Save output gene lists
hESC_WT_YAP1_annot_promoterAnd5_geneSymbol = hESC_WT_YAP1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_TEAD4_annot_promoterAnd5_geneSymbol = hESC_WT_TEAD4_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
    

write.table(hESC_WT_YAP1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_TEAD4_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
            





## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_YAP1_annot_noIntergenic = tibble(hESC_WT_YAP1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_TEAD4_annot_noIntergenic = tibble(hESC_WT_TEAD4_annot) %>%
    filter(annotation != c("Distal Intergenic"))
    

### Save output gene lists
hESC_WT_YAP1_annot_noIntergenic_geneSymbol = hESC_WT_YAP1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_TEAD4_annot_noIntergenic_geneSymbol = hESC_WT_TEAD4_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(hESC_WT_YAP1_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_TEAD4_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_qval1.30103_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
            



## Keep signal everywhere !!!!!!!!!!!!!!!!!!!
hESC_WT_YAP1_annot_all = tibble(hESC_WT_YAP1_annot) 
hESC_WT_TEAD4_annot_all = tibble(hESC_WT_TEAD4_annot)
    
### Save output gene lists
hESC_WT_YAP1_annot_all_geneSymbol = hESC_WT_YAP1_annot_all %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_TEAD4_annot_all_geneSymbol = hESC_WT_TEAD4_annot_all %>%
    dplyr::select(geneSymbol) %>%
    unique()
    

write.table(hESC_WT_YAP1_annot_all_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_qval1.30103_all_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_TEAD4_annot_all_geneSymbol, file = "output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_qval1.30103_all_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
            









```

--> gene feature barplot show YAP1 and TEAD4 mostly in gene body

--> Export gene list and perform Venn diagram

--> Found that very few EZH2 diff bound genes are bound with YAP1. Maybe because YAP1 in intergenic/gene body region. So instead; assing all peak to the nearest TSS; whatever their distance and see if improve colocalization.


