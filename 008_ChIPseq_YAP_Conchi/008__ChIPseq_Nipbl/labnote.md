# Project and goal

Part of gastrulation/QSER1 paper with Conchi Estaras. Notably emboR revision.

--> See `002*/003*/gastrulation paper or QSER1 paper/Revision1 emboR/For Thomas For Revisions/Reviewers comments EMBOR Thomas highlight.docx` for detail: We need to check the overlap between YAP1, QSER1, YAP1:QSER1 co-bound peaks with cohesin protein.




Pipeline (Follow `008*/003*` for all; and `008*/001*` for homer peak calling):
- Download fastq from GEO (past paper from Conchi) with [ChIPseq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1579363) of *Nipbl (component of cohesin loading complex)*
- Clean/map
- Identify peaks using HOMER
- Check for overlap with `008*/001*` and `008*/003*` ChIPs


--> Also for Grant QSER1, analysis of Ser5P-RNAPII untreated from the same GEO (see email 09/03/2025 *"The proposal will investigate the role of QSER1 in transcriptional RNAPII pausing and pause-release of developmental genes through condensate formation"*

)

# Data download from NCBI

Seems that there is only one bio rep. Lets download IP and input:
- [IP Nipbl](https://www.ncbi.nlm.nih.gov/sra?term=SRX833420), 1 Run= SRR1745511
- [input](https://www.ncbi.nlm.nih.gov/sra?term=SRX1036445), 2 Runs= SRR2037028, SRR2037029; I selected both for FASTQ download at [this](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR2037028&display=download) page; the SRR2037028 should contain both SRR2037028 and SRR2037029--> Yes that's correct.


--> ChIPs are SE. Data stored in `008*/008*/input/` and transfered to HPC cluster.


# Rename files

- `SRR1745511.fastq.gz` into `Nipbl.fq.gz`
- `SRR2037028_SRR2037029.fastq.gz` into `input.fq.gz`
- `SRR1745499.fastq.gz` into `Ser5P_RNAPII.fq.gz`

XXXY HERE prusue pol2

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
sbatch scripts/bamtobigwig_unique_extendReads_raw.sh # 40532753 ok
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

Nipbl = as_tibble(read.table("output/homer/Nipbl/peaks.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
    


## Tidy peaks 
Nipbl_gr = makeGRangesFromDataFrame(Nipbl,keep.extra.columns=TRUE)


gr_list <- list(Nipbl=Nipbl_gr)

## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here

                       
### Barplot
pdf("output/ChIPseeker/annotation_barplot_hESC_Nipbl.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()




## Get annotation data frame
Nipbl_annot <- as.data.frame(peakAnnoList[["Nipbl"]]@anno)


## Convert entrez gene IDs to gene symbols
Nipbl_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Nipbl_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Nipbl_annot$gene <- mapIds(org.Hs.eg.db, keys = Nipbl_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



## Save output table
write.table(Nipbl_annot, file="output/ChIPseeker/annotation_homer_hESC_Nipbl_annot.txt", sep="\t", quote=F, row.names=F)



## Keep all ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!

Nipbl_annot_geneSymbol = Nipbl_annot %>%
    dplyr::select(geneSymbol) %>%
    unique()
    



write.table(Nipbl_annot_geneSymbol, file = "output/ChIPseeker/annotation_homer_Nipbl_annot_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
            



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
Nipbl_annot_geneSymbol_promoterAnd5 = tibble(Nipbl_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
    

### Save output gene lists
Nipbl_annot_geneSymbol_promoterAnd5_geneSymbol = Nipbl_annot_geneSymbol_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
    


write.table(Nipbl_annot_geneSymbol_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_homer_Nipbl_annot_geneSymbol_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
            


## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
Nipbl_annot_noIntergenic = tibble(Nipbl_annot) %>%
    filter(annotation != c("Distal Intergenic"))
    


### Save output gene lists
Nipbl_annot_noIntergenic_geneSymbol = Nipbl_annot_noIntergenic %>%
    dplyr::select(geneSymbol) %>%
    unique()
    

write.table(Nipbl_annot_noIntergenic_geneSymbol, file = "output/ChIPseeker/annotation_homer_Nipbl_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
            
            


```


--> All good **(noIntergenic version was used for VennDiagram in the paper)**





# deepTools


Signal in NIPBL peaks for NIPBL, YAP1 (`008*/003*`), QSER1 (`008*/001*`)


```bash
conda activate deeptools


sbatch scripts/matrix_5kb_NIPBLYAP1QSER1_NIPBLpeaks.sh # 40544977 ok
sbatch scripts/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks.sh # 40544978 ok

sbatch scripts/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks-QSER1NIPBL.sh # interactive

```

--> All good, show postive correlation/overlapping between NIPBL, YAP1, and QSER1



# Overlap between QSER1 (008001) and NIPBL peaks - emboR?QSER1? revision tasks, email 5/20/2025





output/annotation_homer_hESC_WT_QSER1_pool_annot.bed



```bash
conda activate BedToBigwig

bedtools intersect -v -a ../007__ENCODE_hESC_histone/output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/homer/Nipbl/peaks.bed | wc -l # 9505 do NOT overlap with NIPBL
bedtools intersect -wa -a ../007__ENCODE_hESC_histone/output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/homer/Nipbl/peaks.bed | wc -l # 2957

bedtools intersect -v -a output/homer/Nipbl/peaks.bed -b ../007__ENCODE_hESC_histone/output/annotation_homer_hESC_WT_QSER1_pool_annot.bed | wc -l # 18351


# 

bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_noHeader.bed -b output/homer/Nipbl/peaks_noHeader.bed | uniq | wc -l # 2956



```


Let's extend the peak size of QSER1 and NIPBL of 200bp up and down; and 500bp up and down and check overlap again.



```bash
# Extend peak size

bedtools slop -i output/homer/Nipbl/peaks_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > output/homer/Nipbl/peaks_noHeader_extend200bp.bed
bedtools slop -i output/homer/Nipbl/peaks_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > output/homer/Nipbl/peaks_noHeader_extend500bp.bed
bedtools slop -i output/homer/Nipbl/peaks_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > output/homer/Nipbl/peaks_noHeader_extend1kp.bed
bedtools slop -i output/homer/Nipbl/peaks_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > output/homer/Nipbl/peaks_noHeader_extend2kp.bed
bedtools slop -i output/homer/Nipbl/peaks_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 5000 > output/homer/Nipbl/peaks_noHeader_extend5kp.bed


## Files extended
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend1kp.bed
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend2kp.bed
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend5kp.bed

output/homer/Nipbl/peaks_noHeader_extend200bp.bed
output/homer/Nipbl/peaks_noHeader_extend500bp.bed
output/homer/Nipbl/peaks_noHeader_extend1kp.bed
output/homer/Nipbl/peaks_noHeader_extend2kp.bed
output/homer/Nipbl/peaks_noHeader_extend5kp.bed

## overlap
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed -b output/homer/Nipbl/peaks_noHeader_extend200bp.bed | uniq | wc -l # 3733
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed -b output/homer/Nipbl/peaks_noHeader_extend500bp.bed | uniq | wc -l # 4219
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend1kp.bed -b output/homer/Nipbl/peaks_noHeader_extend1kp.bed | uniq | wc -l # 4496
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend2kp.bed -b output/homer/Nipbl/peaks_noHeader_extend2kp.bed | uniq | wc -l # 4838
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend5kp.bed -b output/homer/Nipbl/peaks_noHeader_extend5kp.bed | uniq | wc -l # 5908

```



Then overlap with only QSER1 intergenic peaks (ie. removing QSER1 promoter and TSS peaks)

XXXY



