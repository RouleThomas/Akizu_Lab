# Project and goals 

H9 (WA09) at d28
Inhibitor: DOT1Linh, EHMTinh, EZH2inh
IP: H3K27me3, H3K4me3

Re-analysis of CutRun dataset from (Ciceri et al)[https://www.nature.com/articles/s41586-023-06984-8].

EZH2 inhibitors should also target EZH1. So could be considered as an EZH1 and EZH2 inhibitors.
- Identify genes that lose H3K27me3 with inhibitor treatment = putative EZH2 target
- Check whether putative EZH2 target are the same genes that gain H3K27me3 in our EZH1 KO (if yes, would suggest our genes that gain H3K27me3 indeed gain it because of increase EZH2 activity)



# Download data




- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject PRJNA803355 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 16188011 ok

```



## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename_003.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_003.txt
```

--> All good 





# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 16197773 ok

sbatch scripts/fastp_raw_missing.sh # 16220591 ok

```
*ERROR: `EZH2inh_H3K4me3_R1` sequence and quality have different length!*
```ruby
@SRR26553062.12321957 A00333:479:HM2NVDSX3:2:2117:28230:31876/2
TCCTGCAGGGAGCTGGTGCCAGCCGACAGCCGCGCCAGGGCCGCTCCGGG
```
--> Re-launch as `scripts/fastp_raw_missing.sh`; weird the file renaiming is good...!

Here is how solve the [issue](https://github.com/OpenGene/fastp/issues/340):
- Unzipped and opened the fastq (used text editor).
- ctrl+f to sequence number giving the error - turns out was in the same line as sequence before it (i.e., they - are all supposed to start on their own line, but for some reason this one read was running up against the one directly preceding - where there should have been a space/"return" key pressed, there was not).
- Put in a space/pressed "return", putting that read on it's own line.
- Saved, re-zipped and ran.


Try re-download `EZH2inh_H3K4me3_R1`:


```bash
# download
cd tmp/
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/062/SRR26553062/SRR26553062_1.fastq.gz -o SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR265/062/SRR26553062/SRR26553062_2.fastq.gz -o SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER_2.fastq.gz

# fastp
fastp -i SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER_1.fastq.gz -I SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER_2.fastq.gz \
    -o SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER_1.fq.gz -O SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER_2.fq.gz \
    -j SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER -h SRR26553062_GSM7869111_CnR_EZH2inh_H3K4me3_rep1_Homo_sapiens_OTHER

# --> mv and rename manually file to `output/fastp`

```
--> IT WORK!! Files has been corrupted upon downloading!!!




# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch scripts/bowtie2_1.sh # 16221819 ok
sbatch scripts/bowtie2_2.sh # 16221821 ok
```

-->  XXX Looks good; overall ~30-80% uniquely aligned reads XXX
----> Seems less uniquel mapped reads than us but they sequence FAR more depth (~20m reads vs 5 for us)


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-16221819.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_16221819.txt

for file in slurm-16221821.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_16221821.txt


```

Add these values to `/home/roulet/006_Ciceri2024/003__CutRun_EpigeneticInhibitors/samples_003.xlsx`\
Then in R; see `/home/roulet/006_Ciceri2024/006_Ciceri2024.R`.

--> Overall 65-80% input reads as been uniquely mapped to the genome (90% non uniq)



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.


```bash
conda activate bowtie2

sbatch --dependency=afterany:16221819 scripts/samtools_unique_1.sh # 16221854 ok
sbatch --dependency=afterany:16221821 scripts/samtools_unique_2.sh # 16221855 ok


```


# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.



```bash
conda activate deeptools

sbatch --dependency=afterany:16221854 scripts/bamtobigwig_unique_1.sh # 16221861 ok
sbatch --dependency=afterany:16221855 scripts/bamtobigwig_unique_2.sh # 16221862 ok

```


- EZH2inh
PASS: R1, R2
FAIL: NA
- DMSO
PASS: R1, R2
FAIL: NA
- EHMTinh
PASS: XXX
FAIL: XXX
- DOT1Linh
PASS: XXX
FAIL: XXX





## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools

# Generate compile bigwig (.npz) files _ hg38 Akizu analysis
sbatch scripts/multiBigwigSummary_all.sh # 17148596 ok

sbatch scripts/multiBigwigSummary_EZH2inh_H3K27me3noAB.sh # 17148950 xxx

```

--> AB clearly clustered together

--> DMSO and EZH2inh samples cluster well as expected



# THOR without spike in, TMM default norm

No spikein in Ciceri data, so let's use THOR with default (TMM) normalization: EZH2inh vs DMSO


```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge

# Default TMM method
sbatch scripts/THOR_EZH2inh_H3K27me3.sh # 17151035 ok


```

- *NOTE: no DMSO noAB rep2, so I used R1 for both replicates*



Generate median tracks:
```bash
conda activate BedToBigwig
# Default TMM method
sbatch --dependency=afterany:17151035 scripts/bigwigmerge_THOR_EZH2inh_H3K27me3.sh # 17151091 ok

```




## Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks


```R

# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# H3K27me3 TMM default norm wthout spike in
diffpeaks <- read_tsv("output/THOR/THOR_EZH2inh_H3K27me3/EZH2inhH3K27me3-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_DMSO", "count_EZH2inh", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_DMSO, into = c("count_DMSO_1","count_DMSO_2"), sep = ":", convert = TRUE) %>%
  separate(count_EZH2inh, into = c("count_EZH2inh_1","count_EZH2inh_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_EZH2inh_1+count_EZH2inh_2) / (count_DMSO_1+count_DMSO_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_EZH2inh_H3K27me3/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("28dN_DMSO vs EZH2inh") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_EZH2inh_H3K27me3/log2FC_qval50.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("28dN_DMSO vs EZH2inh_qval50") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_EZH2inh_H3K27me3/THOR_qval50.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())




```

**Optimal qvalue:**
--> *H3K27me3*; qval 50 for TMM default no spikein






# deepTool plots


check whether genes that gain H3K27me3 in KO (in 50dN; `007__CutRun`) are EZH2-specific = the one that lose H3K27me3 with EZH2 inhibitor:
- use gain lost H3K27me3 region WT vs KO (`007__CutRun`) *THORq20*
- use gain lost EZH2inh region; *THORq50*


```bash
conda activate deeptools

# deeptool plot on PEAKS
## CutRun__007 gain/lost THOR q20 WT vs KO H3K27me3 region
sbatch scripts/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks.sh # 17164413 xxx



## gain lost EZH2 inh region
### isolate gain / lost peaks
awk -F'\t' '$16 > 1' output/THOR/THOR_EZH2inh_H3K27me3/THOR_qval50.bed > output/THOR/THOR_EZH2inh_H3K27me3/THOR_qval50_gain.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_EZH2inh_H3K27me3/THOR_qval50.bed > output/THOR/THOR_EZH2inh_H3K27me3/THOR_qval50_lost.bed

### deeptool plots
sbatch scripts/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_THOR_q50_peaks.sh # 17164635 xxx



# deeptool plot on GENES
## CutRun__007 gain/lost THOR q20 WT vs KO H3K27me3 region
XXX Need generate gtf XXX sbatch scripts/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_gene.sh #  xxx

## gain lost EZH2 inh region

XXX NEED ASSIGN EZH2 diff peak to genes
XXX Need generate gtf XXX 

```


--> xxx





# ChIPseeker peak gene assignment

## From THOR diff bound peaks
Let's assign **peak to genes from THORs peak**:

**Optimal qvalue** according to IGV:
- 28dN_H3K27me3 DMSO vs EZH2inh; qval 50


--> Assign peak to genes for 28dN:

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


# Import THOR peaks
# TMM defalt
## H3K27me3 _ q50
H3K27me3_q50_pos = as_tibble(read.table('output/THOR/THOR_EZH2inh_H3K27me3/THOR_qval50.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC >1) 
H3K27me3_q50_neg = as_tibble(read.table('output/THOR/THOR_EZH2inh_H3K27me3/THOR_qval50.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, FC=V16) %>%
    filter(FC <1) 
    
# Tidy peaks 
## H3K27me3
H3K27me3_q50_pos_gr = makeGRangesFromDataFrame(H3K27me3_q50_pos,keep.extra.columns=TRUE)
H3K27me3_q50_neg_gr = makeGRangesFromDataFrame(H3K27me3_q50_neg,keep.extra.columns=TRUE)
gr_list <- list(H3K27me3_q50_pos=H3K27me3_q50_pos_gr, H3K27me3_q50_neg=H3K27me3_q50_neg_gr)


# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_THOR_H3K27me3_EZH2inh.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_THOR_H3K27me3_EZH2inh.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame

H3K27me3_q50_pos_annot <- as.data.frame(peakAnnoList[["H3K27me3_q50_pos"]]@anno)
H3K27me3_q50_neg_annot <- as.data.frame(peakAnnoList[["H3K27me3_q50_neg"]]@anno)



## Convert entrez gene IDs to gene symbols

H3K27me3_q50_pos_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_pos_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q50_pos_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_pos_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K27me3_q50_neg_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_neg_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_q50_neg_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_q50_neg_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



## Save output table
write.table(H3K27me3_q50_pos_annot, file="output/ChIPseeker/annotation_THOR_H3K27me3_EZH2inh_q50_pos.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K27me3_q50_neg_annot, file="output/ChIPseeker/annotation_THOR_H3K27me3_EZH2inh_q50_neg.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K27me3_q50_pos_annot_promoterAnd5 = tibble(H3K27me3_q50_pos_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K27me3_q50_neg_annot_promoterAnd5 = tibble(H3K27me3_q50_neg_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))



### Save output gene lists
H3K27me3_q50_pos_annot_promoterAnd5_geneSymbol = H3K27me3_q50_pos_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_q50_neg_annot_promoterAnd5_geneSymbol = H3K27me3_q50_neg_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(H3K27me3_q50_pos_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_EZH2inh_H3K27me3_q50_pos_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_q50_neg_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THOR_EZH2inh_H3K27me3_q50_neg_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


```

--> 4,198 genes gain and 25 genes lost H3K27me3 in resp to EZH2inh (THORq50)



XXXXXXXXXXXXX CHUI AL below not mod XXXXXXXXXXXXX


# MACS2 peak calling on bam unique



--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad_53dN.sh # 15401244 ok
sbatch scripts/macs2_broad_NPC.sh # 15401308 ok

sbatch scripts/macs2_broad_53dN_noIGG.sh # 15401429 ok
```

--> H3K27ac in NPC show 3 peak in R1! And a ~7k peaks in R2. R2 is better, but still very ugly and noisy!

--> 53dN IGG vs not using IGG: almost the same, so let's better use IGG






XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX below not mod

```bash
conda activate bowtie2 # for bedtools
sbatch scripts/macs2_raw_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive
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
- 50dN_KOEF1aEZH1_H3K27me3: 1.30103 (2.3 more true peaks)
- 50dN_KO_H3K27me3: 1.30103 (2.3 more true peaks)
- 50dN_WTQ731E_H3K27me3: 1.30103 (2.3 more true peaks)

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX








