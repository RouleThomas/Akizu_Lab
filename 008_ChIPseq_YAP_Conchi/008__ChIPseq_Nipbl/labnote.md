# Project and goal

Part of gastrulation/QSER1 paper with Conchi Estaras. Notably emboR revision.

--> See `002*/003*/gastrulation paper or QSER1 paper/Revision1 emboR/For Thomas For Revisions/Reviewers comments EMBOR Thomas highlight.docx` for detail: We need to check the overlap between YAP1, QSER1, YAP1:QSER1 co-bound peaks with cohesin protein.




Pipeline (Follow `008*/003*` for all; and `008*/001*` for homer peak calling):
- Download fastq from GEO (past paper from Conchi) with [ChIPseq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1579363) of *Nipbl (component of cohesin loading complex)*
- Clean/map
- Identify peaks using HOMER
- Check for overlap with `008*/001*` and `008*/003*` ChIPs


# Data download from NCBI

Seems that there is only one bio rep. Lets download IP and input:
- [IP Nipbl](https://www.ncbi.nlm.nih.gov/sra?term=SRX833420), 1 Run= SRR1745511
- [input](https://www.ncbi.nlm.nih.gov/sra?term=SRX1036445), 2 Runs= SRR2037028, SRR2037029; I selected both for FASTQ download at [this](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2037028&display=download) page; the SRR2037028 should contain both SRR2037028 and SRR2037029--> XXX Yes that's correct.


--> ChIPs are SE. Data stored in `008*/008*/input/` and transfered to HPC cluster.


# Rename files

`SRR1745511.fastq.gz` into `Nipbl.fq.gz`
`SRR2037028_SRR2037029.fastq.gz` into `input.fq.gz`



# Fastp cleaning

```bash
sbatch scripts/fastp.sh # 40428506 ok
```



# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)
--> *NOTE: I removed `--no-mixed --dovetail` and `-U` (instead of `-r`) for the fastq path as PE options* 


```bash
conda activate bowtie2

sbatch --dependency=afterany:40428506 scripts/bowtie2.sh # 40428600 ok
```


 
--> Looks ok; around 70% uniq aligned reads (97% total) 




## Quality control metrics

XXX NOT RUN BELOW XXX

Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-40428600.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_40428600.txt

```

Add these values to `/home/roulet/008_ChIPseq_YAP_Conchi/008__ChIPseq_Nipbl/samples_008008.xlsx`\
Then in R; see `/home/roulet/008_ChIPseq_YAP_Conchi/ChIPseq_YAP.R`.

--> Overall > XXX % input reads as been uniquely mapped to the genome (XXX % non uniq)

XXX NOT RUN UP XXX



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.

```bash
conda activate bowtie2

sbatch --dependency=afterany:40428600 scripts/samtools_unique_raw.sh # 40428881 ok
```



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
    - input FragL=179, ReadL=51; FragL-ReadL= 128
    - Nipbl FragL=189, ReadL=51; FragL-ReadL= 138


Then generate bigwig with the corresponding fragment size for each sample:



Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads. **HERE single end so need to estimate fragment size!! Extend to `FragL-ReadL`**


```bash
conda activate deeptools

# bigwig with extendReads from CHIPQC
sbatch scripts/bamtobigwig_unique_extendReads_raw.sh # 40532753 xxx
```






## Run HOMER peak caller

From method; use HOMER findPeaks command with input control; Option `-style factor` (for TF; otherwise `-style histone` was used)

- Convert .bam to tagDirectory for [homer](http://homer.ucsd.edu/homer/ngs/peaks.html). Method said only uniquely aligned reads where used, so let's use our `*unique.dupmark.sorted.bam` files
- Call peaks (Method differ if [simplicate](http://homer.ucsd.edu/homer/ngs/peaks.html) or [replicate](http://homer.ucsd.edu/homer/ngs/peaksReplicates.html))





```bash
conda activate homer_deseq2_V1
module load SAMtools/1.16.1-GCC-11.3.0

# Create tagDirectory to be used by homer (FAST, so run in interactive)
makeTagDirectory output/homer/input output/bowtie2/input.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/Nipbl output/bowtie2/Nipbl.unique.dupmark.sorted.bam 
#--> By default homer run in singleend


# PEAK CALLING
## Call peaks with one bio rep (simplicate)
findPeaks output/homer/Nipbl -style factor -o auto -i output/homer/input

## Convert .txt to .bed
### simplicate
pos2bed.pl output/homer/Nipbl/peaks.txt > output/homer/Nipbl/peaks.bed
```



--> All good, peak and bigwig looking good



XXXY HERE double check peak and bigiwg and then do peak annoatiton


# ChIPseeker - homer



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

# SIMPLICATE SAMPLE #####################

# Import homer peaks
# Convert .txt to .bed

hESC_WT_EZH2_R1 = as_tibble(read.table("output/homer/hESC_WT_EZH2_R1/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_EZH2_R2 = as_tibble(read.table("output/homer/hESC_WT_EZH2_R2/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_EZH2_R1 = as_tibble(read.table("output/homer/hESC_YAPKO_EZH2_R1/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_EZH2_R2 = as_tibble(read.table("output/homer/hESC_YAPKO_EZH2_R2/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1_R1 = as_tibble(read.table("output/homer/hESC_WT_QSER1_R1/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1_R2 = as_tibble(read.table("output/homer/hESC_WT_QSER1_R2/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_QSER1_R1 = as_tibble(read.table("output/homer/hESC_YAPKO_QSER1_R1/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_QSER1_R2 = as_tibble(read.table("output/homer/hESC_YAPKO_QSER1_R2/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
## 008*/003*
hESC_WT_TEAD4_R1 = as_tibble(read.table("../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_YAP1_R1 = as_tibble(read.table("../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 


## Tidy peaks 
hESC_WT_EZH2_R1_gr = makeGRangesFromDataFrame(hESC_WT_EZH2_R1,keep.extra.columns=TRUE)
hESC_WT_EZH2_R2_gr = makeGRangesFromDataFrame(hESC_WT_EZH2_R2,keep.extra.columns=TRUE)
hESC_YAPKO_EZH2_R1_gr = makeGRangesFromDataFrame(hESC_YAPKO_EZH2_R1,keep.extra.columns=TRUE)
hESC_YAPKO_EZH2_R2_gr = makeGRangesFromDataFrame(hESC_YAPKO_EZH2_R2,keep.extra.columns=TRUE)
hESC_WT_QSER1_R1_gr = makeGRangesFromDataFrame(hESC_WT_QSER1_R1,keep.extra.columns=TRUE)
hESC_WT_QSER1_R2_gr = makeGRangesFromDataFrame(hESC_WT_QSER1_R2,keep.extra.columns=TRUE)
hESC_YAPKO_QSER1_R1_gr = makeGRangesFromDataFrame(hESC_YAPKO_QSER1_R1,keep.extra.columns=TRUE)
hESC_YAPKO_QSER1_R2_gr = makeGRangesFromDataFrame(hESC_YAPKO_QSER1_R2,keep.extra.columns=TRUE)
hESC_WT_TEAD4_R1_gr = makeGRangesFromDataFrame(hESC_WT_TEAD4_R1,keep.extra.columns=TRUE)
hESC_WT_YAP1_R1_gr = makeGRangesFromDataFrame(hESC_WT_YAP1_R1,keep.extra.columns=TRUE)

gr_list <- list(hESC_WT_EZH2_R1=hESC_WT_EZH2_R1_gr, hESC_WT_EZH2_R2=hESC_WT_EZH2_R2_gr, hESC_YAPKO_EZH2_R1=hESC_YAPKO_EZH2_R1_gr,    hESC_YAPKO_EZH2_R2 = hESC_YAPKO_EZH2_R2_gr, hESC_WT_QSER1_R1 = hESC_WT_QSER1_R1_gr, hESC_WT_QSER1_R2 = hESC_WT_QSER1_R2_gr, hESC_YAPKO_QSER1_R1 = hESC_YAPKO_QSER1_R1_gr, hESC_YAPKO_QSER1_R2 = hESC_YAPKO_QSER1_R2_gr, hESC_WT_TEAD4_R1 = hESC_WT_TEAD4_R1_gr, hESC_WT_YAP1_R1 = hESC_WT_YAP1_R1_gr)

## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here

                       
### Barplot
pdf("output/ChIPseeker/annotation_barplot_hESC_homer.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()




## Get annotation data frame
hESC_WT_EZH2_R1_annot <- as.data.frame(peakAnnoList[["hESC_WT_EZH2_R1"]]@anno)
hESC_WT_EZH2_R2_annot <- as.data.frame(peakAnnoList[["hESC_WT_EZH2_R2"]]@anno)
hESC_YAPKO_EZH2_R1_annot <- as.data.frame(peakAnnoList[["hESC_YAPKO_EZH2_R1"]]@anno)
hESC_YAPKO_EZH2_R2_annot <- as.data.frame(peakAnnoList[["hESC_YAPKO_EZH2_R2"]]@anno)
hESC_WT_QSER1_R1_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1_R1"]]@anno)
hESC_WT_QSER1_R2_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1_R2"]]@anno)
hESC_YAPKO_QSER1_R1_annot <- as.data.frame(peakAnnoList[["hESC_YAPKO_QSER1_R1"]]@anno)
hESC_YAPKO_QSER1_R2_annot <- as.data.frame(peakAnnoList[["hESC_YAPKO_QSER1_R2"]]@anno)
hESC_WT_TEAD4_R1_annot <- as.data.frame(peakAnnoList[["hESC_WT_TEAD4_R1"]]@anno)
hESC_WT_YAP1_R1_annot <- as.data.frame(peakAnnoList[["hESC_WT_YAP1_R1"]]@anno)


## Convert entrez gene IDs to gene symbols
hESC_WT_EZH2_R1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_R1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_EZH2_R1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_R1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_EZH2_R2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_R2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_EZH2_R2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_R2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_YAPKO_EZH2_R1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_EZH2_R1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_YAPKO_EZH2_R1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_EZH2_R1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_YAPKO_EZH2_R2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_EZH2_R2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_YAPKO_EZH2_R2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_EZH2_R2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1_R1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1_R1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1_R1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1_R1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1_R2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1_R2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1_R2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1_R2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_YAPKO_QSER1_R1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_QSER1_R1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_YAPKO_QSER1_R1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_QSER1_R1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_YAPKO_QSER1_R2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_QSER1_R2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_YAPKO_QSER1_R2_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_QSER1_R2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_TEAD4_R1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_TEAD4_R1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_TEAD4_R1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_TEAD4_R1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_YAP1_R1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_YAP1_R1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_YAP1_R1_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_YAP1_R1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



## Save output table
write.table(hESC_WT_EZH2_R1_annot, file="output/ChIPseeker/annotation_homer_hESC_WT_EZH2_R1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_WT_EZH2_R2_annot, file="output/ChIPseeker/annotation_homer_hESC_WT_EZH2_R2_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_YAPKO_EZH2_R1_annot, file="output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_R1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_YAPKO_EZH2_R2_annot, file="output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_R2_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_WT_QSER1_R1_annot, file="output/ChIPseeker/annotation_homer_hESC_WT_QSER1_R1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_WT_QSER1_R2_annot, file="output/ChIPseeker/annotation_homer_hESC_WT_QSER1_R2_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_YAPKO_QSER1_R1_annot, file="output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_R1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_YAPKO_QSER1_R2_annot, file="output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_R2_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_WT_TEAD4_R1_annot, file="output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_WT_YAP1_R1_annot, file="output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot.txt", sep="\t", quote=F, row.names=F) 




## Keep all ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!

hESC_WT_EZH2_R1_annot_geneSymbol = hESC_WT_EZH2_R1_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_R2_annot_geneSymbol = hESC_WT_EZH2_R2_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_R1_annot_geneSymbol = hESC_YAPKO_EZH2_R1_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_R2_annot_geneSymbol = hESC_YAPKO_EZH2_R2_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1_R1_annot_geneSymbol = hESC_WT_QSER1_R1_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1_R2_annot_geneSymbol = hESC_WT_QSER1_R2_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_QSER1_R1_annot_geneSymbol = hESC_YAPKO_QSER1_R1_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_QSER1_R2_annot_geneSymbol = hESC_YAPKO_QSER1_R2_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_TEAD4_R1_annot_geneSymbol = hESC_WT_TEAD4_R1_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_YAP1_R1_annot_geneSymbol = hESC_WT_YAP1_R1_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()




write.table(hESC_WT_EZH2_R1_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_R1_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_EZH2_R2_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_R2_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_R1_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_R1_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_R2_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_R2_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1_R1_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_R1_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1_R2_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_R2_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_QSER1_R1_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_R1_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_QSER1_R2_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_R2_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_TEAD4_R1_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_YAP1_R1_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)








## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_EZH2_R1_annot_promoterAnd5 = tibble(hESC_WT_EZH2_R1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_EZH2_R2_annot_promoterAnd5 = tibble(hESC_WT_EZH2_R2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_YAPKO_EZH2_R1_annot_promoterAnd5 = tibble(hESC_YAPKO_EZH2_R1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_YAPKO_EZH2_R2_annot_promoterAnd5 = tibble(hESC_YAPKO_EZH2_R2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1_R1_annot_promoterAnd5 = tibble(hESC_WT_QSER1_R1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1_R2_annot_promoterAnd5 = tibble(hESC_WT_QSER1_R2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_YAPKO_QSER1_R1_annot_promoterAnd5 = tibble(hESC_YAPKO_QSER1_R1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_YAPKO_QSER1_R2_annot_promoterAnd5 = tibble(hESC_YAPKO_QSER1_R2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_TEAD4_R1_annot_promoterAnd5 = tibble(hESC_WT_TEAD4_R1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_YAP1_R1_annot_promoterAnd5 = tibble(hESC_WT_YAP1_R1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))




### Save output gene lists
hESC_WT_EZH2_R1_annot_promoterAnd5_geneSymbol = hESC_WT_EZH2_R1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_R2_annot_promoterAnd5_geneSymbol = hESC_WT_EZH2_R2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_R1_annot_promoterAnd5_geneSymbol = hESC_YAPKO_EZH2_R1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_R2_annot_promoterAnd5_geneSymbol = hESC_YAPKO_EZH2_R2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1_R1_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1_R1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1_R2_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1_R2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_QSER1_R1_annot_promoterAnd5_geneSymbol = hESC_YAPKO_QSER1_R1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_QSER1_R2_annot_promoterAnd5_geneSymbol = hESC_YAPKO_QSER1_R2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_TEAD4_R1_annot_promoterAnd5_geneSymbol = hESC_WT_TEAD4_R1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_YAP1_R1_annot_promoterAnd5_geneSymbol = hESC_WT_YAP1_R1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()




write.table(hESC_WT_EZH2_R1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_R1_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_EZH2_R2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_R2_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_R1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_R1_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_R2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_R2_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1_R1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_R1_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1_R2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_R2_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_QSER1_R1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_R1_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_QSER1_R2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_R2_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_TEAD4_R1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_YAP1_R1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)









## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_EZH2_R1_annot_noIntergenic = tibble(hESC_WT_EZH2_R1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_EZH2_R2_annot_noIntergenic = tibble(hESC_WT_EZH2_R2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_YAPKO_EZH2_R1_annot_noIntergenic = tibble(hESC_YAPKO_EZH2_R1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_YAPKO_EZH2_R2_annot_noIntergenic = tibble(hESC_YAPKO_EZH2_R2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1_R1_annot_noIntergenic = tibble(hESC_WT_QSER1_R1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1_R2_annot_noIntergenic = tibble(hESC_WT_QSER1_R2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_YAPKO_QSER1_R1_annot_noIntergenic = tibble(hESC_YAPKO_QSER1_R1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_YAPKO_QSER1_R2_annot_noIntergenic = tibble(hESC_YAPKO_QSER1_R2_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_TEAD4_R1_annot_noIntergenic = tibble(hESC_WT_TEAD4_R1_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_YAP1_R1_annot_noIntergenic = tibble(hESC_WT_YAP1_R1_annot) %>%
    filter(annotation != c("Distal Intergenic"))





### Save output gene lists
hESC_WT_EZH2_R1_annot_noIntergenic_geneSymbol = hESC_WT_EZH2_R1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_EZH2_R2_annot_noIntergenic_geneSymbol = hESC_WT_EZH2_R2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()    
hESC_YAPKO_EZH2_R1_annot_noIntergenic_geneSymbol = hESC_YAPKO_EZH2_R1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  
hESC_YAPKO_EZH2_R2_annot_noIntergenic_noIntergenic_geneSymbol = hESC_YAPKO_EZH2_R2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  
hESC_WT_QSER1_R1_annot_noIntergenic_noIntergenic_geneSymbol = hESC_WT_QSER1_R1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  
hESC_WT_QSER1_R2_annot_noIntergenic_noIntergenic_geneSymbol = hESC_WT_QSER1_R2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  
hESC_YAPKO_QSER1_R1_annot_noIntergenic_noIntergenic_geneSymbol = hESC_YAPKO_QSER1_R1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  
hESC_YAPKO_QSER1_R2_annot_noIntergenic_noIntergenic_geneSymbol = hESC_YAPKO_QSER1_R2_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  
hESC_WT_TEAD4_R1_annot_noIntergenic_noIntergenic_geneSymbol = hESC_WT_TEAD4_R1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  
hESC_WT_YAP1_R1_annot_noIntergenic_noIntergenic_geneSymbol = hESC_WT_YAP1_R1_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  




write.table(hESC_WT_EZH2_R1_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_R1_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_EZH2_R2_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_R2_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)        
write.table(hESC_YAPKO_EZH2_R1_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_R1_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)        
write.table(hESC_YAPKO_EZH2_R2_annot_noIntergenic_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_R2_annot_noIntergenic_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)     
write.table(hESC_WT_QSER1_R1_annot_noIntergenic_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_R1_annot_noIntergenic_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)    
write.table(hESC_WT_QSER1_R2_annot_noIntergenic_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_R2_annot_noIntergenic_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)  
write.table(hESC_YAPKO_QSER1_R1_annot_noIntergenic_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_R1_annot_noIntergenic_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)  
write.table(hESC_YAPKO_QSER1_R2_annot_noIntergenic_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_R2_annot_noIntergenic_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)  
write.table(hESC_WT_TEAD4_R1_annot_noIntergenic_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot_noIntergenic_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)  
write.table(hESC_WT_YAP1_R1_annot_noIntergenic_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot_noIntergenic_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)  







# BIO REP SAMPLE #####################

# Import homer peaks
# Convert .txt to .bed
hESC_WT_EZH2_pool = as_tibble(read.table("output/homer/hESC_WT_EZH2_outputPeaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_EZH2_pool = as_tibble(read.table("output/homer/hESC_YAPKO_EZH2_outputPeaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_QSER1_pool = as_tibble(read.table("output/homer/hESC_WT_QSER1_outputPeaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_YAPKO_QSER1_pool = as_tibble(read.table("output/homer/hESC_YAPKO_QSER1_outputPeaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 




## Tidy peaks 
hESC_WT_EZH2_pool_gr = makeGRangesFromDataFrame(hESC_WT_EZH2_pool,keep.extra.columns=TRUE)
hESC_YAPKO_EZH2_pool_gr = makeGRangesFromDataFrame(hESC_YAPKO_EZH2_pool,keep.extra.columns=TRUE)
hESC_WT_QSER1_pool_gr = makeGRangesFromDataFrame(hESC_WT_QSER1_pool,keep.extra.columns=TRUE)
hESC_YAPKO_QSER1_pool_gr = makeGRangesFromDataFrame(hESC_YAPKO_QSER1_pool,keep.extra.columns=TRUE)

gr_list <- list(hESC_WT_EZH2_pool=hESC_WT_EZH2_pool_gr, hESC_YAPKO_EZH2_pool=hESC_YAPKO_EZH2_pool_gr, hESC_WT_QSER1_pool=hESC_WT_QSER1_pool_gr,    hESC_YAPKO_QSER1_pool = hESC_YAPKO_QSER1_pool_gr)

## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
               
### Barplot
pdf("output/ChIPseeker/annotation_barplot_hESC_pool_homer.pdf", width=14, height=3)
plotAnnoBar(peakAnnoList)
dev.off()

## Get annotation data frame
hESC_WT_EZH2_pool_annot <- as.data.frame(peakAnnoList[["hESC_WT_EZH2_pool"]]@anno)
hESC_YAPKO_EZH2_pool_annot <- as.data.frame(peakAnnoList[["hESC_YAPKO_EZH2_pool"]]@anno)
hESC_WT_QSER1_pool_annot <- as.data.frame(peakAnnoList[["hESC_WT_QSER1_pool"]]@anno)
hESC_YAPKO_QSER1_pool_annot <- as.data.frame(peakAnnoList[["hESC_YAPKO_QSER1_pool"]]@anno)

## Convert entrez gene IDs to gene symbols
hESC_WT_EZH2_pool_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_pool_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_EZH2_pool_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_EZH2_pool_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_YAPKO_EZH2_pool_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_EZH2_pool_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_YAPKO_EZH2_pool_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_EZH2_pool_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_WT_QSER1_pool_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1_pool_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_WT_QSER1_pool_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_WT_QSER1_pool_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
hESC_YAPKO_QSER1_pool_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_QSER1_pool_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
hESC_YAPKO_QSER1_pool_annot$gene <- mapIds(org.Hs.eg.db, keys = hESC_YAPKO_QSER1_pool_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(hESC_WT_EZH2_pool_annot, file="output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_YAPKO_EZH2_pool_annot, file="output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_pool_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_WT_QSER1_pool_annot, file="output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(hESC_YAPKO_QSER1_pool_annot, file="output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_pool_annot.txt", sep="\t", quote=F, row.names=F) 


## Keep all ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!

hESC_WT_EZH2_pool_annot_geneSymbol = hESC_WT_EZH2_pool_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_pool_annot_geneSymbol = hESC_YAPKO_EZH2_pool_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1_pool_annot_geneSymbol = hESC_WT_QSER1_pool_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_QSER1_pool_annot_geneSymbol = hESC_YAPKO_QSER1_pool_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(hESC_WT_EZH2_pool_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_pool_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_pool_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1_pool_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_QSER1_pool_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_pool_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_EZH2_pool_annot_promoterAnd5 = tibble(hESC_WT_EZH2_pool_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_YAPKO_EZH2_pool_annot_promoterAnd5 = tibble(hESC_YAPKO_EZH2_pool_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_WT_QSER1_pool_annot_promoterAnd5 = tibble(hESC_WT_QSER1_pool_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
hESC_YAPKO_QSER1_pool_annot_promoterAnd5 = tibble(hESC_YAPKO_QSER1_pool_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
    
### Save output gene lists
hESC_WT_EZH2_pool_annot_promoterAnd5_geneSymbol = hESC_WT_EZH2_pool_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_pool_annot_promoterAnd5_geneSymbol = hESC_YAPKO_EZH2_pool_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_WT_QSER1_pool_annot_promoterAnd5_geneSymbol = hESC_WT_QSER1_pool_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_QSER1_pool_annot_promoterAnd5_geneSymbol = hESC_YAPKO_QSER1_pool_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(hESC_WT_EZH2_pool_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_pool_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_pool_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_WT_QSER1_pool_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_QSER1_pool_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_pool_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
hESC_WT_EZH2_pool_annot_noIntergenic = tibble(hESC_WT_EZH2_pool_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_YAPKO_EZH2_pool_annot_noIntergenic = tibble(hESC_YAPKO_EZH2_pool_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_WT_QSER1_pool_annot_noIntergenic = tibble(hESC_WT_QSER1_pool_annot) %>%
    filter(annotation != c("Distal Intergenic"))
hESC_YAPKO_QSER1_pool_annot_noIntergenic = tibble(hESC_YAPKO_QSER1_pool_annot) %>%
    filter(annotation != c("Distal Intergenic"))

### Save output gene lists
hESC_WT_EZH2_pool_annot_noIntergenic_geneSymbol = hESC_WT_EZH2_pool_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
hESC_YAPKO_EZH2_pool_annot_noIntergenic_geneSymbol = hESC_YAPKO_EZH2_pool_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()    
hESC_WT_QSER1_pool_annot_noIntergenic_geneSymbol = hESC_WT_QSER1_pool_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  
hESC_YAPKO_QSER1_pool_annot_noIntergenic_geneSymbol = hESC_YAPKO_QSER1_pool_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()  

write.table(hESC_WT_EZH2_pool_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(hESC_YAPKO_EZH2_pool_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_EZH2_pool_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)        
write.table(hESC_WT_QSER1_pool_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)        
write.table(hESC_YAPKO_QSER1_pool_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_hESC_YAPKO_QSER1_pool_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)     
            



# Plot annot peak location with WT samples of interest

# Convert .txt to .bed
hESC_WT_QSER1_pool = as_tibble(read.table("output/homer/hESC_WT_QSER1_outputPeaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_EZH2_pool = as_tibble(read.table("output/homer/hESC_WT_EZH2_outputPeaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)     
hESC_WT_TEAD4_R1 = as_tibble(read.table("../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
hESC_WT_YAP1_R1 = as_tibble(read.table("../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
    


## Tidy peaks 
hESC_WT_QSER1_pool_gr = makeGRangesFromDataFrame(hESC_WT_QSER1_pool,keep.extra.columns=TRUE)
hESC_WT_EZH2_pool_gr = makeGRangesFromDataFrame(hESC_WT_EZH2_pool,keep.extra.columns=TRUE)
hESC_WT_TEAD4_R1_gr = makeGRangesFromDataFrame(hESC_WT_TEAD4_R1,keep.extra.columns=TRUE)
hESC_WT_YAP1_R1_gr = makeGRangesFromDataFrame(hESC_WT_YAP1_R1,keep.extra.columns=TRUE)


gr_list <- list(hESC_WT_QSER1_pool=hESC_WT_QSER1_pool_gr, hESC_WT_EZH2_pool=hESC_WT_EZH2_pool_gr, hESC_WT_TEAD4_R1= hESC_WT_TEAD4_R1_gr, hESC_WT_YAP1_R1=hESC_WT_YAP1_R1_gr)

## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here

                       
### Barplot
pdf("output/ChIPseeker/annotation_barplot_hESC_WT_homer.pdf", width=14, height=3)
plotAnnoBar(peakAnnoList)
dev.off()



```


--> All good (noIntergenic version was used for VennDiagram in the paper)




