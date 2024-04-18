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
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads. **HERE single end so need to estaimate fragment size!! Extend to `FragL-ReadL`**


```bash
conda activate deeptools

# bigwig with extendReads 50 (Default ~250bp fragment length  150bp reads +100 bp)
sbatch --dependency=afterany:17725583 scripts/bamtobigwig_unique_extendReads100_CPC.sh # 17762131 xxx
sbatch scripts/bamtobigwig_unique_extendReads100_hESC.sh # 17762114 xxx

# bigwig with extendReads from CHIPQC
sbatch scripts/bamtobigwig_unique_extendReads_hESC.sh # 17770516 xxx
sbatch scripts/bamtobigwig_unique_extendReads_CPC.sh #  xxx



```

XXXXXXXXXXXX

- KOEF1aEZH1
*Pass*: H3K27me3
*Failed*: EZH1cs, EZH2, SUZ12
- KO
*Pass*: H3K27me3
*Failed*: EZH1cs, EZH2, SUZ12
- WTQ731E
*Pass*: H3K27me3
*Failed*: EZH1cs, EZH2, SUZ12
- WT (PSC)
*Pass*: NA
*Failed*: EZH1cs and H3K27me1




## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_50dN.sh # 9064423 ok 



# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels 50dN_KOEF1aEZH1_EZH1cs_R1 50dN_KOEF1aEZH1_EZH1cs_R2 50dN_KOEF1aEZH1_EZH2_R1 50dN_KOEF1aEZH1_EZH2_R2 50dN_KOEF1aEZH1_H3K27me3_R1 50dN_KOEF1aEZH1_H3K27me3_R2 50dN_KOEF1aEZH1_IGG_R1 50dN_KOEF1aEZH1_SUZ12_R1 50dN_KOEF1aEZH1_SUZ12_R2 50dN_KO_EZH1cs_R1 50dN_KO_EZH1cs_R2 50dN_KO_EZH2_R1 50dN_KO_EZH2_R2 50dN_KO_H3K27me3_R1 50dN_KO_H3K27me3_R2 50dN_KO_IGG_R1 50dN_KO_IGG_R2 50dN_KO_SUZ12_R1 50dN_KO_SUZ12_R2 50dN_WTQ731E_EZH1cs_R1 50dN_WTQ731E_EZH1cs_R2 50dN_WTQ731E_EZH2_R1 50dN_WTQ731E_EZH2_R2 50dN_WTQ731E_H3K27me3_R1 50dN_WTQ731E_H3K27me3_R2 50dN_WTQ731E_H3K27me3_R3 50dN_WTQ731E_IGG_R1 50dN_WTQ731E_IGG_R2 50dN_WTQ731E_SUZ12_R1 50dN_WTQ731E_SUZ12_R2 PSC_WT_EZH1cs_01FA PSC_WT_EZH1cs_1FA PSC_WT_H3K27me1_01FA PSC_WT_H3K27me1_1FA \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf



## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels 50dN_KOEF1aEZH1_EZH1cs_R1 50dN_KOEF1aEZH1_EZH1cs_R2 50dN_KOEF1aEZH1_EZH2_R1 50dN_KOEF1aEZH1_EZH2_R2 50dN_KOEF1aEZH1_H3K27me3_R1 50dN_KOEF1aEZH1_H3K27me3_R2 50dN_KOEF1aEZH1_IGG_R1 50dN_KOEF1aEZH1_SUZ12_R1 50dN_KOEF1aEZH1_SUZ12_R2 50dN_KO_EZH1cs_R1 50dN_KO_EZH1cs_R2 50dN_KO_EZH2_R1 50dN_KO_EZH2_R2 50dN_KO_H3K27me3_R1 50dN_KO_H3K27me3_R2 50dN_KO_IGG_R1 50dN_KO_IGG_R2 50dN_KO_SUZ12_R1 50dN_KO_SUZ12_R2 50dN_WTQ731E_EZH1cs_R1 50dN_WTQ731E_EZH1cs_R2 50dN_WTQ731E_EZH2_R1 50dN_WTQ731E_EZH2_R2 50dN_WTQ731E_H3K27me3_R1 50dN_WTQ731E_H3K27me3_R2 50dN_WTQ731E_H3K27me3_R3 50dN_WTQ731E_IGG_R1 50dN_WTQ731E_IGG_R2 50dN_WTQ731E_SUZ12_R1 50dN_WTQ731E_SUZ12_R2 PSC_WT_EZH1cs_01FA PSC_WT_EZH1cs_1FA PSC_WT_H3K27me1_01FA PSC_WT_H3K27me1_1FA \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf


```

--> two big groups: H3K27me3 IP versus the other
----> Seems only H3K27me3 IP has worked here




# MACS2 peak calling on bam unique
