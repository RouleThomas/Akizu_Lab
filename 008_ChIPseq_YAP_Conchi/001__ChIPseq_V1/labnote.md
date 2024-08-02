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
## Raw and DiffBindTMM bigwig

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


# Bigwig with DiffBind TMM scaling factor
sbatch scripts/bamtobigwig_unique_extendReads_hESC_QSER1_EZH2.sh # 18757417 ok

```


--> All work in hESC

--> CPC; YAP1 and TEAD4 only R3 work (thats ok it is their control sample! And they already have YAP1 sequenced), NR2F2 R1 and R2 work (not R3), 

--> *RPGC* bigwig make replicate very heterogeneous... = *BAD*; Use instead the **extendReads** = **GOOD**

--> For QSER1/EZH2 WT, seems **THOR and DiffBind_TMM are very similar**


## Pearson correlation heatmap on bigwig signals

- `008002` current samples
- `008002` with `008001` for QSER1 sample comparison (raw and DiffBind_TMM bigwigs)
- `008002` with `008001` for QSER1 sample comparison (raw and DiffBind_TMM bigwigs) with QSER1 and EZH2


```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_extendReads_hESC.sh # 17863824 ok
sbatch scripts/multiBigwigSummary_CPC.sh # 17863829 ok

sbatch scripts/multiBigwigSummary_008001008002_QSER1.sh # 18718911 ok
sbatch scripts/multiBigwigSummary_008001008002_QSER1_geneOnly.sh # 18719262 ok

sbatch scripts/multiBigwigSummary_008001008002_QSER1_DiffBindTMM.sh # 18868031 ok
sbatch scripts/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_bin2kb.sh # 18868238 xxx

sbatch scripts/multiBigwigSummary_008001008002_QSER1EZH2_DiffBindTMM_bin2kb.sh # 18868395 xxx


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

plotPCA -in output/bigwig/multiBigwigSummary_008001008002_QSER1.npz \
    --transpose \
    --ntop 0 \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --markers o o x x o o x \
    --colors "black" "black" "grey" "grey" "darkblue" "darkblue" "grey" \
    -o output/bigwig/multiBigwigSummary_008001008002_QSER1_plotPCA.pdf
plotPCA -in output/bigwig/multiBigwigSummary_008001008002_QSER1_geneOnly.npz \
    --transpose \
    --ntop 0 \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --markers o o x x o o x \
    --colors "black" "black" "grey" "grey" "darkblue" "darkblue" "grey" \
    -o output/bigwig/multiBigwigSummary_008001008002_QSER1_geneOnly_plotPCA.pdf
plotPCA -in output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM.npz \
    --transpose \
    --ntop 0 \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --markers o o x x o o x \
    --colors "black" "black" "grey" "grey" "darkblue" "darkblue" "grey" \
    -o output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_plotPCA.pdf
plotPCA -in output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_bin2kb.npz \
    --transpose \
    --ntop 0 \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --markers o o x x o o x \
    --colors "black" "black" "grey" "grey" "darkblue" "darkblue" "grey" \
    -o output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_bin2kb_plotPCA.pdf

plotPCA -in output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_bin2kb.npz \
    --transpose \
    --ntop 0 \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --markers o o x x o o x \
    --colors "black" "black" "grey" "grey" "darkblue" "darkblue" "grey" \
    --log2 \
    -o output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_bin2kb_plotPCAlog2.pdf

plotPCA -in output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1EZH2_DiffBindTMM_bin2kb.npz \
    --transpose \
    --ntop 0 \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_EZH2_R1 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 hESC_WT_EZH2_R1 hESC_WT_EZH2_R2 \
    --markers o o x x o o o x o o \
    --colors "black" "black" "grey" "grey" "red" "darkblue" "darkblue" "grey" "darkred" "darkred" \
    -o output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1EZH2_DiffBindTMM_bin2kb_plotPCA.pdf



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
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_008001008002_QSER1.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_008001008002_QSER1_heatmap.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_008001008002_QSER1_geneOnly.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_008001008002_QSER1_geneOnly_heatmap.pdf
plotCorrelation \
    -in output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_heatmap.pdf
plotCorrelation \
    -in output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_bin2kb.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_bin2kb_heatmap.pdf

plotCorrelation \
    -in output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1EZH2_DiffBindTMM_bin2kb.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels hESC_WT_QSER1FLAG_R1 hESC_WT_QSER1FLAG_R2 hESC_WT_inputFLAG_R1 hESC_WT_inputFLAG_R2 hESC_WT_EZH2_R1 hESC_WT_QSER1_R1 hESC_WT_QSER1_R2 hESC_WT_input_R1 hESC_WT_EZH2_R1 hESC_WT_EZH2_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1EZH2_DiffBindTMM_bin2kb_heatmap.pdf


```

--> two big groups: H3K27me3 IP versus the other
----> Seems only H3K27me3 IP has worked here

--> `*_bin2kb` include the following modifications at the `multiBigwigSummary`: `--blackListFileName /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --chromosomesToSkip chrX chrY chrM --binSize 2000`


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





# ChIPseeker - THOR peaks

- Annotate EZH2 diff bound peaks + Check whether these genes are bound with YAP1 (in `001003`)



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

## EZH2 hESC WT vs YAPKO q5
EZH2_pos = as_tibble(read.table('output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval5_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
EZH2_neg = as_tibble(read.table('output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval5_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

 ## EZH2 hESC WT vs YAPKO q4
EZH2_pos = as_tibble(read.table('output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
EZH2_neg = as_tibble(read.table('output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 

## Tidy peaks #-->> Re-Run from here with different qvalue!!
EZH2_pos_gr = makeGRangesFromDataFrame(EZH2_pos,keep.extra.columns=TRUE)
EZH2_neg_gr = makeGRangesFromDataFrame(EZH2_neg,keep.extra.columns=TRUE)


gr_list <- list(EZH2_pos=EZH2_pos_gr, EZH2_neg=EZH2_neg_gr)



## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here


### Barplot
pdf("output/ChIPseeker/annotation_barplot_THORq5_EZH2.pdf", width=14, height=5)
pdf("output/ChIPseeker/annotation_barplot_THORq4_EZH2.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()



## Get annotation data frame
EZH2_pos_annot <- as.data.frame(peakAnnoList[["EZH2_pos"]]@anno)
EZH2_neg_annot <- as.data.frame(peakAnnoList[["EZH2_neg"]]@anno)



## Convert entrez gene IDs to gene symbols
EZH2_pos_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH2_pos_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH2_pos_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH2_pos_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
EZH2_neg_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH2_neg_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH2_neg_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH2_neg_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



## Save output table
write.table(EZH2_pos_annot, file="output/ChIPseeker/annotation_THORq4_EZH2_pos_annot.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(EZH2_neg_annot, file="output/ChIPseeker/annotation_THORq4_EZH2_neg_annot.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
EZH2_pos_annot_promoterAnd5 = tibble(EZH2_pos_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
EZH2_neg_annot_promoterAnd5 = tibble(EZH2_neg_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
EZH2_pos_annot_promoterAnd5_geneSymbol = EZH2_pos_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
EZH2_neg_annot_promoterAnd5_geneSymbol = EZH2_neg_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(EZH2_pos_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(EZH2_neg_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


# Check wether EZH2 diff bound genes are YAP1 bound 001003
## import EZH2 diff bound genes
EZH2_pos_annot_promoterAnd5_geneSymbol = read.table("output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt", header = FALSE, sep = "\t") %>%
    add_column(EZH2 = "gain") %>%
    dplyr::rename("geneSymbol" ="V1")

EZH2_neg_annot_promoterAnd5_geneSymbol = read.table("output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt", header = FALSE, sep = "\t") %>%
    add_column(EZH2 = "lost") %>%
    dplyr::rename("geneSymbol" ="V1")


EZH2_posNeg_annot_promoterAnd5_geneSymbol = EZH2_pos_annot_promoterAnd5_geneSymbol %>%
    bind_rows(EZH2_neg_annot_promoterAnd5_geneSymbol) %>%
    as_tibble()

## import YAP1 bound genes in WT 008003
YAP1_annot_noIntergenic_geneSymbol = read.table("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_qval1.30103_noIntergenic_geneSymbol.txt", header = FALSE, sep = "\t") %>%
    dplyr::rename("geneSymbol" ="V1") %>%
    as_tibble()
YAP1_annot_all_geneSymbol = read.table("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_qval1.30103_all_geneSymbol.txt", header = FALSE, sep = "\t") %>%
    dplyr::rename("geneSymbol" ="V1") %>%
    as_tibble()


## Add a column YAP1, yes (YAP1 is binding) or no (YAP1 is not binding)
EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1binding = EZH2_posNeg_annot_promoterAnd5_geneSymbol %>%
  mutate(YAP1 = ifelse(geneSymbol %in% YAP1_annot_noIntergenic_geneSymbol$geneSymbol, "yes", "no"))

EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1binding_all = EZH2_posNeg_annot_promoterAnd5_geneSymbol %>%
  mutate(YAP1 = ifelse(geneSymbol %in% YAP1_annot_all_geneSymbol$geneSymbol, "yes", "no"))

write.table(EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1binding_all, file = "output/ChIPseeker/EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1binding_all.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)
```

--> *EZH2*; THORq5, 325/156 gene gain/lost

--> *EZH2*; THORq4, 406/216 gene gain/lost

--> Very few genes with differential EZH2 binding upon YAP1KO are bound with YAP1 when using YAP1 in promoter or 5.






# Functional analysis with enrichR

- **Up and down reg genes**: THOR diff bound genes 
- **Unique list of genes**: QSER1, EZH2 (YAP1, TEAD4) bound genes in WT (macs2)
- **Up and down reg genes**: THOR diff bound genes for EZH2 and bound with YAP1 in WT


- NOTE: I did NOT using webtool Venn diagram I isolated the specific genes that gain / lost H3K*me3 (eg.; like: `output/ChIPseeker/annotation_THOR_H3K*me3_q*_pos_promoterAnd5_geneSymbol_Venndiagram*.txt`)


**IMPOPRTANT NOTE: Run the reading and processing ONE BY ONE !!! Otherwise, lead to bug!!!!**

```R
# library
library("tidyverse")
library("enrichR")
library("ggrepel")

# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 

### GeneSymbol list of signif gain/lost EZH2 in WT vs KO
output/ChIPseeker/annotation_THORq5_EZH2_neg_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq5_EZH2_pos_annot_promoterAnd5_geneSymbol.txt

output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt

output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2TEAD4.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1.txt
output/ChIPseeker/none.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"
# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)
# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics 

pdf("output/GO/enrichR_GO_Biological_Process_2023_THORq5_EZH2.pdf", width=8, height=10)
pdf("output/GO/enrichR_GO_Biological_Process_2023_THORq4_EZH2.pdf", width=8, height=10)

pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_overlap_hESC_WT_QSER1EZH2.pdf", width=8, height=7)
pdf("output/GO/enrichR_GO_Biological_Process_2023_Venn_overlap_hESC_WT_QSER1EZH2YAP1.pdf", width=8, height=7)

pdf("output/GO/enrichR_GO_Biological_Process_2023_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.pdf", width=10, height=8)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="EZH2",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "GO_Biological_Process_2023") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("GO_Molecular_Function_2023") # 

### GeneSymbol list of signif gain/lost EZH2 in WT vs KO
output/ChIPseeker/annotation_THORq5_EZH2_neg_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq5_EZH2_pos_annot_promoterAnd5_geneSymbol.txt

output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt


output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2TEAD4.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1.txt
output/ChIPseeker/none.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)
# Extracting KEGG data and assigning types
up <- eup$GO_Molecular_Function_2023
down <- edown$GO_Molecular_Function_2023
up$type <- "up"
down$type <- "down"
# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)
# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)
# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string


## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics

pdf("output/GO/enrichR_GO_Molecular_Function_2023_THORq4_EZH2.pdf", width=8, height=6)
pdf("output/GO/enrichR_GO_Molecular_Function_2023_Venn_overlap_hESC_WT_QSER1EZH2YAP1.pdf", width=8, height=3)

pdf("output/GO/enrichR_GO_Molecular_Function_2023_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.pdf", width=10, height=6)

pdf("output/GO/enrichR_GO_Molecular_Function_2023_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.pdf", width=10, height=6)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="EZH2",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "GO_Molecular_Function_2023") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()

## save output
write.table(gos, "output/GO/enrichR_GO_Molecular_Function_2023_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Define databases for enrichment
dbs <- c("GO_Cellular_Component_2023") # 

### GeneSymbol list of signif gain/lost EZH2 in WT vs KO
output/ChIPseeker/annotation_THORq5_EZH2_neg_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq5_EZH2_pos_annot_promoterAnd5_geneSymbol.txt

output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt

output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2TEAD4.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Cellular_Component_2023
down <- edown$GO_Cellular_Component_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)



# Plotting with enhanced aesthetics

pdf("output/GO/enrichR_GO_Cellular_Component_2023_THORq4_EZH2.pdf", width=8, height=2)
pdf("output/GO/enrichR_GO_Cellular_Component_2023_Venn_overlap_hESC_WT_QSER1EZH2YAP1.pdf", width=8, height=3)

pdf("output/GO/enrichR_GO_Cellular_Component_2023_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.pdf", width=8, height=3)
pdf("output/GO/enrichR_GO_Cellular_Component_2023_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.pdf", width=8, height=3)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="EZH2",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "GO_Cellular_Component_2023") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Cellular_Component_2023_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment
dbs <- c("KEGG_2021_Human") # 

### GeneSymbol list of signif gain/lost EZH2 in WT vs KO
output/ChIPseeker/annotation_THORq5_EZH2_neg_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq5_EZH2_pos_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt

output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2TEAD4.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$KEGG_2021_Human
down <- edown$KEGG_2021_Human
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 20)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 20)
# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)
# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)
# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)


# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics

pdf("output/GO/enrichR_KEGG_2021_Human_THORq4_EZH2.pdf", width=8, height=3)
pdf("output/GO/enrichR_KEGG_2021_Human_Venn_overlap_hESC_WT_QSER1EZH2TEAD4.pdf", width=8, height=3)
pdf("output/GO/enrichR_KEGG_2021_Human_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.pdf", width=8, height=2)

pdf("output/GO/enrichR_KEGG_2021_Human_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.pdf", width=8, height=2)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name="EZH2",   # H3K27me3  H3K4me3
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "KEGG_2021_Human") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt", sep="\t", row.names=FALSE, quote=FALSE)


```



# Functional analysis with enrichGO (single list of genes dotplot)


We will use clusterProfile package. Tutorial [here](https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html).

Let's do a test of the pipeline with genes from cluster4 amd cluster14 from the rlog counts. Our background list will be all genes tested for differential expression.

**IMPORTANT NOTE: When doing GO, do NOT set a universe (background list of genes) it perform better!**

```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")

# GO


# Genes that gain H3K27me3 in NPC (009)
## Files
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2TEAD4.txt
output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1.txt

list = read_csv("output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1.txt", col_names = "gene_name") # CHANGE HERE!!!!!!!!
  
ego <- enrichGO(gene = as.character(list$gene_name), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "CC",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
# pdf("output/GO/dotplot_CC_Venn_overlap_hESC_WT_QSER1EZH2_top20.pdf", width=7, height=7) # CHANGE HERE!!!!!!!!
# pdf("output/GO/dotplot_CC_Venn_overlap_hESC_WT_QSER1EZH2TEAD4_top20.pdf", width=7, height=7) # CHANGE HERE!!!!!!!!
pdf("output/GO/dotplot_CC_Venn_overlap_hESC_WT_QSER1EZH2YAP1_top20.pdf", width=7, height=7) # CHANGE HERE!!!!!!!!

dotplot(ego, showCategory=20)
dev.off()

write.csv(data.frame(ego) , file = "output/GO/dotplot_CC_Venn_overlap_hESC_WT_QSER1EZH2YAP1_top20.txt", row.names=FALSE)

```








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

# THOR
sbatch scripts/bigwigmerge_THOR_CPC.sh # 17869306 ok
sbatch --dependency=afterany:17869560 scripts/bigwigmerge_THOR_hESC.sh # 17869572 ok


# raw bigwig
sbatch scripts/bigwigmerge_raw_CPC.sh # 18717584 ok
sbatch scripts/bigwigmerge_raw_hESC.sh # 18717671 xxx


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

# WTvsKO QSER1
diffpeaks <- read_tsv("output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_hESC_QSER1_WTvsYAPKO/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT vs KO") +
  theme_bw()
dev.off()

## --> Very few diff peaks!


# WTvsKO EZH2
diffpeaks <- read_tsv("output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_hESC_EZH2_WTvsYAPKO/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT vs KO") +
  theme_bw()
dev.off()


pdf("output/THOR/THOR_hESC_EZH2_WTvsYAPKO/log2FC_qval4.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 4) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT vs KO_qval5") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 4) %>%
  write_tsv("output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 5) %>%
  group_by(X6) %>%
  summarise(n = n())
```

- *NOTE: FC positive = less in KO; negative = more in KO*

**Optimal qvalue:**
--> *EZH2*; qval 5 looks great!; but qval4 better as we have NODAL :) !!





Isolate positive and negative THOR peaks to display deepTool plots

```bash
# positive negative peaks
## qval 4
awk -F'\t' '$16 > 1' output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4.bed > output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4.bed > output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_negative.bed

## qval 5
awk -F'\t' '$16 > 1' output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval5.bed > output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval5_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval5.bed > output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval5_negative.bed
```




# deepTool plots

## Non overlaping config

Venn diagram of peak-genes TED4, QSER1 and EZH2 has been generated.

- Check whether QSER1 flank EZH2 (confirm [Dixon2021 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8185639/))
- Check TED4 vs QSER1 (and EZH2)
--> Venn diagram generated online and gene list exported manually (`ChIPseeker/Venn_overlap_hESC_WT_*.txt`)


**--> use THOR bigwig for QSER1 and EZH2 and use raw unique bigwig for TEAD4 (from `008003`)**

Generate gtf file from gene list; 
- all genes
- gene with peak in promoter (qval macs2 2.3)
- QSER1 comparison Conchi vs Dixon (`008002`)
- EZH2 comparison Conchi vs Dixon (`008002`)
- QSER1, EZH2, YAP1 
- gene with differential EZH2 binding (THORq4)

```bash
### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
## Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_hESC_WT_QSER1TEAD4.txt > output/ChIPseeker/Venn_overlap_hESC_WT_QSER1TEAD4_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2TEAD4.txt > output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2TEAD4_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2.txt > output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_qval1.30103_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_qval1.30103_noIntergenic_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval10_noIntergenic_geneSymbol.txt > ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval10_noIntergenic_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval15_noIntergenic_geneSymbol.txt > ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval15_noIntergenic_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval20_noIntergenic_geneSymbol.txt > ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval20_noIntergenic_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval1.30103_promoterAnd5_geneSymbol.txt > ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval3_promoterAnd5_geneSymbol.txt > ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval3_promoterAnd5_as_gtf_geneSymbol.txt


sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1.txt > output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_as_gtf_geneSymbol.txt
cat output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt | sort | uniq | sed 's/\r$//; s/.*/gene_name "&"/' > output/ChIPseeker/annotation_THORq4_EZH2_posneg_annot_promoterAnd5_as_gtf_geneSymbol.txt # put together pos and neg



## Filter the gtf
grep -Ff output/ChIPseeker/Venn_overlap_hESC_WT_QSER1TEAD4_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_hESC_WT_QSER1TEAD4.gtf
grep -Ff output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_hESC_WT_QSER1EZH2.gtf
grep -Ff output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2TEAD4_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_hESC_WT_QSER1EZH2TEAD4.gtf

grep -Ff output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_qval1.30103_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_macs2_hESC_WT_QSER1_qval1.30103_noIntergenic.gtf
grep -Ff output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_macs2_hESC_WT_EZH2_qval1.30103_promoterAnd5.gtf

grep -Ff ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval10_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > ../002__ChIPseq_Dixon2021/meta/ENCFF159KBI_macs2_hESC_WT_QSER1FLAG_pool_qval10_noIntergenic.gtf
grep -Ff ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval15_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > ../002__ChIPseq_Dixon2021/meta/ENCFF159KBI_macs2_hESC_WT_QSER1FLAG_pool_qval15_noIntergenic.gtf
grep -Ff ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_QSER1FLAG_pool_qval20_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > ../002__ChIPseq_Dixon2021/meta/ENCFF159KBI_macs2_hESC_WT_QSER1FLAG_pool_qval20_noIntergenic.gtf

grep -Ff ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > ../002__ChIPseq_Dixon2021/meta/ENCFF159KBI_macs2_hESC_WT_EZH2_R1_qval1.30103_promoterAnd5.gtf
grep -Ff ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_R1_qval3_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > ../002__ChIPseq_Dixon2021/meta/ENCFF159KBI_macs2_hESC_WT_EZH2_R1_qval3_promoterAnd5.gtf

grep -Ff output/ChIPseeker/Venn_overlap_hESC_WT_QSER1EZH2YAP1_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gtf


grep -Ff output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_THORq4_EZH2_pos_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_THORq4_EZH2_neg_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_THORq4_EZH2_posneg_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_THORq4_EZH2_posneg_annot_promoterAnd5.gtf




# deeptool plots
## all genes
sbatch scripts/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.sh # 18879644 ok
sbatch scripts/matrix_gene_1kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.sh # 19802287 ok
sbatch scripts/matrix_gene_500bp_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.sh # 19823389 ok
sbatch scripts/matrix_gene_500bp_rawBigwig_QSER1EZH2YAP1TEAD4DVL2input_WT_allGenes.sh # 19956770 ok



## all genes EZH2, 5mC 008004
sbatch scripts/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.sh # 19836877 ok


## macs2 008001
sbatch scripts/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1TEAD4.sh # 18526196 ok
sbatch scripts/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2.sh # 18526206 ok
sbatch scripts/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2TEAD4.sh # 18526209 ok
sbatch scripts/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2only.sh # 18877131 ok
sbatch scripts/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.sh # 18879428 ok



## macs2 QSER1 008001 vs 008002 with output/bigwig raw files
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008001noIntergenic.sh # 18718029 ok
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval10.sh # 18718073 ok
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15.sh # 18718096 ok
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval20.sh # 18718131 ok


## macs2 EZH2 008001 vs 008002 with output/bigwig raw files
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008001promoterAnd5.sh # 18869271 ok
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103.sh # 18869286 ok
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval3.sh # 18869293 ok



```

--> QSER1 before TSS; EZH2 right after TSS

--> QSER1 and QSER1FLAG (`008002`) seems to well co-localize; same signal whatever list of gene used...

--> EZH1 Conchi and Dixon (`008002`) seems to well co-localize; same signal whatever list of gene used...



## EZH2 gain lost WT vs KO

```bash
# deeptool plots
## WT vs YAPKO; EZH2 signal
sbatch scripts/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegative_peak.sh # 18526815 ok
sbatch scripts/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.sh # 18542691 ok
sbatch scripts/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2gene.sh # 22311569 ok



## WT vs YAPKO; EZH2 gain lost with EZH2 and QSER1 signal in WT vs YAPKO
sbatch scripts/matrix_TSS_5kb_THOREZH2QSER1_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.sh # 19826509 ok
sbatch scripts/matrix_TSS_5kb_THOREZH2QSER1_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2gene.sh # 19943017 ok


## WT vs YAPKO; EZH2 gain lost with EZH2 and m5C signal in WT (`008004`)
sbatch scripts/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.sh # 19831082 ok
sbatch scripts/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.sh # 19835895 ok
sbatch scripts/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.sh # 19835895 ok



## Genes with EZH2 binding changes; EZH2, 5mC 008004
sbatch scripts/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene.sh # 19942671 ok
sbatch scripts/matrix_TSS_5kb_THOREZH25mCENCSR744UGB_hESCWTvsYAPKO_positiveNegativeTHORq4Separated_EZH2gene.sh # 20253104 ok

```

- *NOTE: qval4 include NODAL*

--> gain lost EZH2 regions are all bound with QSER1. No correlation with QSER1 binding changes. 

--> `_sorted` file is needed for heatmap in R `## 20240720 ppt task`

## 20240720 ppt task

Lets' generate a heatmap of the gene that gain and lost EZH2 WT vs cYAP1KO; and label some genes (gene label cannot be done with deeptools). --> `ComplexHeatmap` in R

Matrix need to be sorted; like in the plot generated by deeptools `plotHeatmap`. 

```bash
# sort the matrix
conda activate deeptools

computeMatrixOperations sort -m output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2gene.gz -R meta/ENCFF159KBI_THORq4_EZH2_pos_annot_promoterAnd5.gtf meta/ENCFF159KBI_THORq4_EZH2_neg_annot_promoterAnd5.gtf -o output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2gene_sorted.gz

conda activate deseq2
```

```R
# Load required libraries
library("ComplexHeatmap")
library("circlize")
library("data.table")
library("tidyverse")

# Read the gzipped matrix data directly
matrix_data <- fread("zcat output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2gene_sorted.gz", header = TRUE, data.table = FALSE)

# Keep the ENST gene information as row names and subset the numerical data for the heatmap
row_names <- matrix_data[, 4]
heatmap_data <- as.matrix(matrix_data[, -(1:6)])

# Replace NA values with 0 or any suitable value
heatmap_data[is.na(heatmap_data)] <- 0

# Set the row names
rownames(heatmap_data) <- row_names

# Example gene names and their positions
# Replace this with actual gene names of interest
genes_of_interest <- c("PAX7", "PAX3", "NEUROG2-AS1", "CYP26A1", "IRX3", "OTX2", "LHX8", "NEUROD1", "NKX2-4", "ZIC1", "ZIC4", "HES2", "WLS", "GDF7", "FOXA2", "TBX1", "SOX17", "FGF4", "KLF4", "WNT3", "NODAL", "CDX2", "GATA5")

# Create a logical vector indicating rows that are genes of interest
gene_labels <- ifelse(rownames(heatmap_data) %in% genes_of_interest, rownames(heatmap_data), NA)


# Draw the heatmap 
pdf("output/deeptools/pheatmap_matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2gene.pdf", width=5, height=6)
pheatmap(heatmap_data,
         cluster_rows = FALSE,  # Disable row clustering to speed up
         cluster_cols = FALSE,   # Enable column clustering
         show_rownames = TRUE,  # Show row names
         show_colnames = FALSE,  # Remove column names (X-axis labels)
         annotation_row = data.frame(Gene = gene_labels))  # Annotate specific genes
dev.off()

```
--> Work, but heatmap look very bad.. Deeptool plot heatmap are much better. Let's troubleshoot to add gene on the current heatmap

--> Let's just collect row number of each gene; or simply which one is the more gain/lost (FC value of THOR peak)


```R
# packages
library("tidyverse")
library("tidyverse")


# import
## peak to gene
output/ChIPseeker/annotation_THORq4_EZH2_pos_annot.txt
output/ChIPseeker/annotation_THORq4_EZH2_neg_annot.txt

## THOR peak FC value
output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_positive.bed
output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_negative.bed

# POSITIVE
# isolate genes / peak
genes_of_interest <- c("PAX7", "PAX3", "NEUROG2-AS1", "CYP26A1", "IRX3", "OTX2", "LHX8", "NEUROD1", "NKX2-4", "ZIC1", "ZIC4")
annotation_THORq4_EZH2_pos_annot = as_tibble( read_tsv("output/ChIPseeker/annotation_THORq4_EZH2_pos_annot.txt") ) %>%
    filter(geneSymbol %in% c("PAX7", "PAX3", "NEUROG2-AS1", "CYP26A1", "IRX3", "OTX2", "LHX8", "NEUROD1", "NKX2-4", "ZIC1", "ZIC4") ) %>%
    dplyr::select(geneSymbol, annotation, distanceToTSS, name)
## Keep the peak closest to the TSS for each gene
annotation_THORq4_EZH2_pos_annot_peakTSS = annotation_THORq4_EZH2_pos_annot %>%
  group_by(geneSymbol) %>%
  filter(abs(distanceToTSS) == min(abs(distanceToTSS))) %>%
  ungroup()
## Collect FC information of the selected peaks
THOR_qval4_positive = as_tibble(read_tsv("output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_positive.bed", col_names = FALSE)) %>%
    dplyr::rename( "name"= "X4",  "FC" = "X16") %>%
    dplyr::select(name, FC) %>%
    inner_join(annotation_THORq4_EZH2_pos_annot_peakTSS) %>%
    arrange(desc(FC))


# NEGATIVE
# isolate genes / peak
genes_of_interest <- c("HES2", "WLS", "GDF7", "FOXA2", "TBX1", "SOX17", "FGF4", "KLF4", "WNT3", "NODAL", "CDX2", "GATA5")
annotation_THORq4_EZH2_neg_annot = as_tibble( read_tsv("output/ChIPseeker/annotation_THORq4_EZH2_neg_annot.txt") ) %>%
    filter(geneSymbol %in% c("HES2", "WLS", "GDF7", "FOXA2", "TBX1", "SOX17", "FGF4", "KLF4", "WNT3", "NODAL", "CDX2", "GATA5") ) %>%
    dplyr::select(geneSymbol, annotation, distanceToTSS, name)
## Keep the peak closest to the TSS for each gene
annotation_THORq4_EZH2_neg_annot_peakTSS = annotation_THORq4_EZH2_neg_annot %>%
  group_by(geneSymbol) %>%
  filter(abs(distanceToTSS) == min(abs(distanceToTSS))) %>%
  ungroup()
## Collect FC information of the selected peaks
THOR_qval4_negative = as_tibble(read_tsv("output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_negative.bed", col_names = FALSE)) %>%
    dplyr::rename( "name"= "X4",  "FC" = "X16") %>%
    dplyr::select(name, FC) %>%
    inner_join(annotation_THORq4_EZH2_neg_annot_peakTSS) %>%
    arrange((FC))




```


# enhancer (with deeptool plot)

Let's look at ChIPseq signal in the enhancer regions from *GeneHancer data*:
- Heatmap / profile; WT ChIPseq signal (hopefully YAP1 will be bound to some enhancer)
- Collect YAP1 bound regions (macs2) directly overlapping with enhancer regions
    - Collect enhancer ID and all associated target genes
    - Check whether these genes are EZH2 diff bound


**Heatmap / profile; WT ChIPseq signal (hopefully YAP1 will be bound to some enhancer)**

```bash

sbatch scripts/matrix_TSS_5kb_rawBigwig_YAP1EZH2QSER1TEAD4DVL2_WT_GeneHancer.sh # 19945005 ok
sbatch scripts/matrix_TSS_5kb_rawBigwig_YAP1EZH2IGG_WT_GeneHancer.sh # 19944996 ok

```

--> YAP1 and other mark not specifically bound to enhancer; seems random, input signal is similar


**Collect YAP1 bound regions (macs2) directly overlapping with enhancer regions**


```bash
bedtools intersect -wa -a meta/GeneHancer_v5.20.bed -b ../003__ChIPseq_pluripotency/output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_YAP1_R1_peaks.broadPeak > meta/GeneHancer_v5_hESC_WT_YAP1_R1_peaks.bed

bedtools intersect -wa -a meta/GeneHancer_v5.20.bed -b output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_EZH2_pool_peaks.broadPeak > meta/GeneHancer_v5_hESC_WT_EZH2_pool_peaks.bed

# To check how many YAP1 (or other) peak overlap with enhancer
bedtools intersect -wa -a ../003__ChIPseq_pluripotency/output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_YAP1_R1_peaks.broadPeak -b meta/GeneHancer_v5.20.bed  > meta/hESC_WT_YAP1_R1_peaks_GeneHancer_v5.bed 

bedtools intersect -wa -a output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_EZH2_pool_peaks.broadPeak -b meta/GeneHancer_v5.20.bed > meta/hESC_WT_EZH2_pool_peaks_GeneHancer_v5.bed 

bedtools intersect -wa -a output/macs2/broad/broad_blacklist_qval1.30103/hESC_WT_QSER1_pool_peaks.broadPeak -b meta/GeneHancer_v5.20.bed > meta/hESC_WT_QSER1_pool_peaks_GeneHancer_v5.bed 

```
--> Work, `meta/GeneHancer_v5_hESC_WT_YAP1_R1_peaks.bed` contains enhancer that overlap with YAP1 peak; 1,746 unique enhancers 
----> Nearly all YAP1 peaks overlap with enhancer (1,772/1,798 YAP1 peak overlap with enhancer!)
----> EZH2, a bit less (3,918/4,288 EZH2 peak overlap with enhancer)
----> QSER1, a bit less (13,416/14,170 EZH2 peak overlap with enhancer)



Let's now modify the `meta/GeneHancer_v5.20.gff` to create a file with column enhancer_id and corresponding target gene with scores.




```bash
awk -F'\t' '{
    # Extract the attributes column
    attr = $9;
    # Initialize variables
    genehancer_id = "";
    connected_gene = "";
    score = "";
    
    # Split attributes by semicolon
    n = split(attr, attrs, ";");
    for (i = 1; i <= n; i++) {
        # Trim leading and trailing whitespace
        gsub(/^ +| +$/, "", attrs[i]);
        # Extract genehancer_id
        if (attrs[i] ~ /genehancer_id=/) {
            genehancer_id = substr(attrs[i], index(attrs[i], "=") + 1);
        }
        # Extract connected_gene
        if (attrs[i] ~ /connected_gene=/) {
            connected_gene = substr(attrs[i], index(attrs[i], "=") + 1);
        }
        # Extract score
        if (attrs[i] ~ /score=/) {
            score = substr(attrs[i], index(attrs[i], "=") + 1);
        }
    }
    # Print the extracted values
    if (genehancer_id != "" && connected_gene != "" && score != "") {
        print genehancer_id "\t" connected_gene "\t" score;
    }
}' meta/GeneHancer_v5.20.gff > meta/GeneHancer_v5.20_tidy.txt
```

Let's now go in R to check our YAP1 target enhacner gene and whether these genes experienced EZH2 binding changes.

```R
# packages
library("tidyverse")

# Files
meta/GeneHancer_v5.20_tidy.txt # enhancer gene target
meta/GeneHancer_v5_hESC_WT_YAP1_R1_peaks.bed # enhancer bound with YAP1
EZH2_pos_annot_promoterAnd5_geneSymbol # gene gain EZH2
EZH2_neg_annot_promoterAnd5_geneSymbol # gene lost EZH2

# import
## EZH2 gain lost
EZH2_pos_annot_promoterAnd5_geneSymbol = read.table("output/ChIPseeker/annotation_THORq4_EZH2_pos_annot_promoterAnd5_geneSymbol.txt", header = FALSE, sep = "\t") %>%
    add_column(EZH2 = "gain") %>%
    dplyr::rename("geneSymbol" ="V1")
EZH2_neg_annot_promoterAnd5_geneSymbol = read.table("output/ChIPseeker/annotation_THORq4_EZH2_neg_annot_promoterAnd5_geneSymbol.txt", header = FALSE, sep = "\t") %>%
    add_column(EZH2 = "lost") %>%
    dplyr::rename("geneSymbol" ="V1")

EZH2_posNeg_annot_promoterAnd5_geneSymbol = EZH2_pos_annot_promoterAnd5_geneSymbol %>%
    bind_rows(EZH2_neg_annot_promoterAnd5_geneSymbol) %>%
    as_tibble()

## enhancer
GeneHancer_v5 = as_tibble(read.table("meta/GeneHancer_v5.20_tidy.txt", header = FALSE, sep = "\t") )  %>%
    dplyr::rename("genehancer_id" ="V1", "geneSymbol" = "V2", "score" = "V3")

GeneHancer_v5_hESC_WT_YAP1_R1_peaks = as_tibble(read.table("meta/GeneHancer_v5_hESC_WT_YAP1_R1_peaks.bed", header = FALSE, sep = "\t") )  %>%
    dplyr::rename("genehancer_id" = "V4") %>%
    dplyr::select(genehancer_id) %>%
    unique()


## isolate gene target by enhacner where YAP1 is bound
GeneHancer_v5_hESC_WT_YAP1_R1_peaks_geneSymbol = GeneHancer_v5_hESC_WT_YAP1_R1_peaks %>% 
    left_join(GeneHancer_v5) %>%
    dplyr::select(geneSymbol) %>%
    unique()



## Check wether these genes are among the EZH2 diff bound genes

## Add a column YAP1, yes (YAP1 is binding) or no (YAP1 is not binding)
EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1Enhancerbinding = EZH2_posNeg_annot_promoterAnd5_geneSymbol %>%
  mutate(YAP1_enhancer = ifelse(geneSymbol %in% GeneHancer_v5_hESC_WT_YAP1_R1_peaks_geneSymbol$geneSymbol, "yes", "no"))

## also add the column for regular binding in promoter
YAP1_annot_noIntergenic_geneSymbol = read.table("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_qval1.30103_noIntergenic_geneSymbol.txt", header = FALSE, sep = "\t") %>%
    dplyr::rename("geneSymbol" ="V1") %>%
    as_tibble()

EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1binding = EZH2_posNeg_annot_promoterAnd5_geneSymbol %>%
  mutate(YAP1_promoter = ifelse(geneSymbol %in% YAP1_annot_noIntergenic_geneSymbol$geneSymbol, "yes", "no"))

EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding = EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1Enhancerbinding %>%
    left_join(EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1binding)


write.table(EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding, file = "output/ChIPseeker/EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)


# count genes with YAP1 binding in enhancer or promoter and EZH2 changes

EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding = as_tibble(read.table("output/ChIPseeker/EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding.txt", 
                        header = TRUE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE) )


EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding %>%
  filter(YAP1_promoter == "yes" | YAP1_enhancer == "yes")

EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding = 
EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding %>%
  filter(YAP1_promoter == "yes" | YAP1_enhancer == "yes",
         EZH2 == "gain") %>%
  dplyr::select(geneSymbol)

EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding = 
EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding %>%
  filter(YAP1_promoter == "yes" | YAP1_enhancer == "yes",
         EZH2 == "lost")%>%
  dplyr::select(geneSymbol)

## save output list of genes for functional analysis
write.table(EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding, file = "output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)

write.table(EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding, file = "output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)



EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding = 
EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding %>%
  filter(YAP1_promoter == "yes" | YAP1_enhancer == "yes",
         EZH2 == "lost")%>%
  dplyr::select(geneSymbol)

## save output list of genes for functional analysis

EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding = 
EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding %>%
  filter(YAP1_promoter == "yes",
         EZH2 == "gain") %>%
  dplyr::select(geneSymbol)

EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding = 
EZH2_posNeg_annot_promoterAnd5_geneSymbol_YAP1EnhancerPromoterbinding %>%
  filter(YAP1_promoter == "yes",
         EZH2 == "lost")%>%
  dplyr::select(geneSymbol)

  
write.table(EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding, file = "output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)

write.table(EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding, file = "output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE)


## command to find out which enhancer-YAP1 bound target which genes 

GeneHancer_v5_gene = GeneHancer_v5 %>%
    filter(geneSymbol == "SIX2") # !!!!!!!!!!!!!  CHANGE GENE OF INTEREST HERE  !!!!!!!!!!!!!!!!!

GeneHancer_v5_gene %>%
    inner_join(GeneHancer_v5_hESC_WT_YAP1_R1_peaks)


```

--> Not much genes with EZH2 binding changes and bound with enhancer where YAP1 is binding








# DiffBind TMM normalization (for sample without spikein)


Generate DiffBind TMM scaling factor and apply them to the bam. That will generate clean seq + complexity depth normalized bigwigs. --> **Exact same process as E coli spike in norm, just no library size correction!**

--> Apply Reciprocal SF to bamtobgiwg!! (see `001003` for detail at `### Histone-spike-in and TMM-norm scaled-bigwig OR DiffBind-ScaleFactor`)


- hESC QSER1 WT and YAPKO
- hESC EZH2 WT and YAPKO


```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 


## hESC QSER1 and EZH2 for WT and YAPKO (NOT GOOD!!)
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1EZH2.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1EZH2.RData")
load("output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1EZH2.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_QSER1EZH2.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_QSER1EZH2.pdf", width=14, height=20) 
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
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_TMM_unique_SF_hESC_QSER1EZH2.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_QSER1EZH2_blackgreylist_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_QSER1EZH2_blackgreylist_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()



## hESC QSER1 for WT and YAPKO
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1.RData")
load("output/DiffBind/meta_sample_macs2raw_unique_hESC_QSER1.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_QSER1.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_QSER1.pdf", width=14, height=20) 
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
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_TMM_unique_SF_hESC_QSER1.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_QSER1_blackgreylist_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_QSER1_blackgreylist_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()



## hESC EZH2 for WT and YAPKO
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_hESC_EZH2.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/meta_sample_macs2raw_unique_hESC_EZH2.RData")
load("output/DiffBind/meta_sample_macs2raw_unique_hESC_EZH2.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_EZH2.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_EZH2.pdf", width=14, height=20) 
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
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_TMM_unique_SF_hESC_EZH2.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_hESC_EZH2_blackgreylist_TMM.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_hESC_EZH2_blackgreylist_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

```

Can we **use different AB in the same DiffBind correction??** Will the SF provided be the same as AB per AB??
--> NOT the same number; need to be run AB per AB!!




# Quantify signal around TSS 


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

# EZH2 #####################################################

## WT EZH2_in Peaks 1kb
bwFile <- "output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_1kbTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2_1kb")
#### Collect output and gene information
WT_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_EZH2_1kb.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- unique(WT_EZH2$gene)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)

WT_EZH2_geneSymbol = WT_EZH2 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()

write.table(WT_EZH2_geneSymbol, file = "output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)



## YAPKO EZH2_in Peaks 1kb
bwFile <- "output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s2_median.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_1kbTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/YAPKO_EZH2")
#### Collect output and gene information
YAPKO_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/YAPKO_EZH2.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()
YAPKO_EZH2_geneSymbol = YAPKO_EZH2 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()
write.table(YAPKO_EZH2_geneSymbol, file = "output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)






## WT EZH2_in Peaks 500bp
bwFile <- "output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_500bpTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2_500bp")
#### Collect output and gene information
WT_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_EZH2_500bp.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()
WT_EZH2_geneSymbol = WT_EZH2 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()
write.table(WT_EZH2_geneSymbol, file = "output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)



## YAPKO EZH2_in Peaks 500bp
bwFile <- "output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s2_median.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_500bpTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/YAPKO_EZH2_500bp")
#### Collect output and gene information
YAPKO_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/YAPKO_EZH2_500bp.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()
YAPKO_EZH2_geneSymbol = YAPKO_EZH2 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()
write.table(YAPKO_EZH2_geneSymbol, file = "output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)








## WT EZH2_in Peaks 250bp
bwFile <- "output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_250bpTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_EZH2_250bp")
#### Collect output and gene information
WT_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_EZH2_250bp.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()
WT_EZH2_geneSymbol = WT_EZH2 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()
write.table(WT_EZH2_geneSymbol, file = "output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)



## YAPKO EZH2_in Peaks 250bp
bwFile <- "output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s2_median.bw"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_250bpTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/YAPKO_EZH2_250bp")
#### Collect output and gene information
YAPKO_EZH2 <- as_tibble(fread(cmd = "gunzip -c output/binBw/YAPKO_EZH2_250bp.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()
YAPKO_EZH2_geneSymbol = YAPKO_EZH2 %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc),
           bc_median = median(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_median, bc_max) %>%
    unique()
write.table(YAPKO_EZH2_geneSymbol, file = "output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)





```

--> Works well. File generated `output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt` containing max and median EZH2 signal for all genes.


# Correlation EZH2 signal with Pseudotime 002*/003*

## Test pseudotime_start_end_association() - FAIL


Related to *20240720_meeting* Let's check whether correlation between:
- EZH2 signal in TSS (this labnote), 1kb/500bp/250bp
- Pseudotime of trajectory2 (`002*/003*` at `### Time Course effect ##`)
- If correlation, compare WT vs KO

```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj2 <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_traj2_noCondition_humangastruloid72hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 

# WT #########################
WT_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt")
WT_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt")
WT_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj2_WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj2 %>%
    left_join(WT_EZH2_1kbTSS_geneSymbol)
pseudotime_traj2_WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj2 %>%
    left_join(WT_EZH2_500bpTSS_geneSymbol)
pseudotime_traj2_WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj2 %>%
    left_join(WT_EZH2_250bpTSS_geneSymbol)

# correlation scale values to -1 +1 
scale_to_range <- function(x) {
  return((x - min(x)) / (max(x) - min(x)) * 2 - 1)
}
# correlation raw
pseudotime_traj2fdr005__WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj2_WT_EZH2_1kbTSS_geneSymbol %>%
    filter(fdr<0.05, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median))  %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))
pseudotime_traj2fdr001__WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj2_WT_EZH2_1kbTSS_geneSymbol %>%
    filter(fdr<0.01, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median))  %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))
pseudotime_traj2fdr005__WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj2_WT_EZH2_500bpTSS_geneSymbol %>%
    filter(fdr<0.05, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median))  %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))
pseudotime_traj2fdr001__WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj2_WT_EZH2_500bpTSS_geneSymbol %>%
    filter(fdr<0.01, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median)) %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))

pseudotime_traj2fdr005__WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj2_WT_EZH2_250bpTSS_geneSymbol %>%
    filter(fdr<0.05, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median))  %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))
pseudotime_traj2fdr001__WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj2_WT_EZH2_250bpTSS_geneSymbol %>%
    filter(fdr<0.01, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median)) %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))

# pdf("output/binBw/corr_pseudotime_traj2fdr001__WT_EZH2_1kbTSS_geneSymbol_bc_median.pdf", width=5, height=4)
# pdf("output/binBw/corr_pseudotime_traj2fdr001__WT_EZH2_500bpTSS_geneSymbol_bc_median.pdf", width=5, height=4)
# pdf("output/binBw/corr_pseudotime_traj2fdr001__WT_EZH2_250bpTSS_geneSymbol_bc_median.pdf", width=5, height=4)

pdf("output/binBw/corr_pseudotime_traj2fdr005__WT_EZH2_1kbTSS_geneSymbol_bc_median_scaled.pdf", width=5, height=4)
ggplot(pseudotime_traj2fdr005__WT_EZH2_1kbTSS_geneSymbol, aes(x = logFClineage1_scaled, y = bc_median_scaled)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()





# YAPKO ###################


YAPKO_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt")
YAPKO_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt")
YAPKO_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj2_YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj2 %>%
    left_join(YAPKO_EZH2_1kbTSS_geneSymbol)
pseudotime_traj2_YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj2 %>%
    left_join(YAPKO_EZH2_500bpTSS_geneSymbol)
pseudotime_traj2_YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj2 %>%
    left_join(YAPKO_EZH2_250bpTSS_geneSymbol)


## correlation scale values to -1 +1 
scale_to_range <- function(x) {
  return((x - min(x)) / (max(x) - min(x)) * 2 - 1)
}
## correlation raw
pseudotime_traj2fdr005__YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj2_YAPKO_EZH2_1kbTSS_geneSymbol %>%
    filter(fdr<0.05, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median))  %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))
pseudotime_traj2fdr001__YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj2_YAPKO_EZH2_1kbTSS_geneSymbol %>%
    filter(fdr<0.01, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median))  %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))
pseudotime_traj2fdr005__YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj2_YAPKO_EZH2_500bpTSS_geneSymbol %>%
    filter(fdr<0.05, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median))  %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))
pseudotime_traj2fdr001__YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj2_YAPKO_EZH2_500bpTSS_geneSymbol %>%
    filter(fdr<0.01, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median)) %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))

pseudotime_traj2fdr005__YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj2_YAPKO_EZH2_250bpTSS_geneSymbol %>%
    filter(fdr<0.05, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median))  %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))
pseudotime_traj2fdr001__YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj2_YAPKO_EZH2_250bpTSS_geneSymbol %>%
    filter(fdr<0.01, logFClineage1>0, !is.na(logFClineage1), !is.na(bc_median)) %>%
  mutate(logFClineage1_scaled = scale_to_range(logFClineage1),
         bc_median_scaled = scale_to_range(bc_median))

# pdf("output/binBw/corr_pseudotime_traj2fdr001__YAPKO_EZH2_1kbTSS_geneSymbol_bc_median.pdf", width=5, height=4)
# pdf("output/binBw/corr_pseudotime_traj2fdr001__YAPKO_EZH2_500bpTSS_geneSymbol_bc_median.pdf", width=5, height=4)
# pdf("output/binBw/corr_pseudotime_traj2fdr001__YAPKO_EZH2_250bpTSS_geneSymbol_bc_median.pdf", width=5, height=4)

pdf("output/binBw/corr_pseudotime_traj2fdr005__YAPKO_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
ggplot(pseudotime_traj2fdr005__YAPKO_EZH2_1kbTSS_geneSymbol, aes(x = logFClineage1_scaled, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()







pdf("output/binBw/corr_pseudotime_test.pdf", width=5, height=4)
ggplot(pseudotime_traj2fdr005__WT_EZH2_1kbTSS_geneSymbol %>% filter(logFClineage1 > 0), aes(x = logFClineage1, y = bc_median)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 1, label.y = 0, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()




```


--> Scaling the values between -1 +1 does not change corr profile...

--> Using Start vs End `pseudotime_start_end_association()` is not working, poor correlation. Likely because the logFC from pseudotime_start_end_association() does not reflect activation time point... Rather difference between Start and End of the trajectory....





## Test Activation point with EZH2 - FAIL



Related to *20240720_meeting* Let's check whether correlation between:
- EZH2 signal in TSS (this labnote), 1kb/500bp/250bp
- Pseudotime Activation point of trajectory2 (`002*/003*` at `## Identify Activation point`) or `pseudotime_start_end_association()`
    - Also check only conserved cell type marker genes (top50, top100, top500, signif only..) (`002*/003*`-->`output/seurat/srat_all_conserved_markers_V2corr.txt`) *_corr is just removing of 1st column name*
    - Also check only genes with EZH2 peak in WT (`output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt`)
- If correlation, compare WT vs KO

```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj2_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_noCondition_humangastruloid72hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj2_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj2_noCondition_humangastruloid72hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj2_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_traj2_noCondition_humangastruloid72hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_StartEnd" = "fdr")  %>% 
    dplyr::select(geneSymbol, fdr_StartEnd, logFClineage1)
EZH2_peaks = read_tsv("output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt", col_names = FALSE) %>%
    dplyr::rename("geneSymbol" = "X1")
## cell type marker genes (very badly formated)
cellType_marker = read_tsv("../../002_scRNAseq/003__YAP1/output/seurat/srat_all_conserved_markers_V2corr.txt") %>%
    tidyr::separate(gene, into = c("gene", "trash"), sep = "\\.\\.\\.") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    dplyr::select(cluster, geneSymbol, max_pval)
###



pseudotime_traj2_peak_DEG_StartEnd = pseudotime_traj2_peak %>%
    left_join(pseudotime_traj2_DEG) %>%
    left_join(pseudotime_traj2_StartEnd)

# WT #########################
WT_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt")
WT_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt")
WT_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
# YAPKO #########################
YAPKO_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt")
YAPKO_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt")
YAPKO_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj2_peak_YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))

## signal EZH2 vs DEG timecourse

# pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdr0__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_DEGstartEndfdr0001__WT_EZH2_1kbTSS_geneSymbol_bc_median.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd <0.05,
           logFClineage1 > 0,
           smooth_peak_pseudotime>-1,
           bc_max>0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_logFClineage1Over5___WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol %>% 
    filter(
           logFClineage1 > 5) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrStartEnd05logFClineage1Over05___WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrDEG05logFClineage1Over1___WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol %>% 
    filter(
           logFClineage1 > 1, 
           fdr_DEG <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 500, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrDEG05logFClineage1Over1___YAPKO_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_YAPKO_EZH2_1kbTSS_geneSymbol %>% 
    filter(
           logFClineage1 > 1, 
           fdr_DEG <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 500, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()





pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_logFClineage1Over1_EZH2peaks___WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
EZH2_peaks %>% 
    left_join(pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol) %>% 
    filter(
           logFClineage1 > 1) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_EZH2_peaks__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
EZH2_peaks %>% 
    left_join(pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol) %>%
    filter(logFClineage1 > 0, fdr_StartEnd <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_test.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol %>% 
    filter( bc_max> 200 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()




## signal EZH2 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 200) %>%
  ungroup()
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_500__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_200__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_100 %>% 
    filter( ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


#### 250bp TSS
cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 200) %>%
  ungroup()
pseudotime_traj2_peak_WT_EZH2_250bpTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_EZH2_250bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_500__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_200__WT_EZH2_250bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_250bpTSS_geneSymbol_cellType_marker_100 %>% 
    filter(  ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



#### 1kb TSS
cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 5000) %>%
  ungroup()
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_500__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_5000__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol_cellType_marker_100 %>% 
    filter( logFClineage1 > 0 , fdr_DEG <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()






### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_signif05__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()




#### 1kb
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_signif05__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol_cellType_marker_signif %>% 
    filter(  ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_signif05_EZH2_peaks__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
EZH2_peaks %>% 
    left_join(pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol_cellType_marker_signif) %>% 
    filter(  ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


## Random testing
pdf("output/binBw/corr_pseudotime_test.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG <0.05,
           logFClineage1 > 0,
           smooth_peak_pseudotime>-1,
           bc_max>0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 1000) %>%
  ungroup()
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

pdf("output/binBw/corr_pseudotime_test.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol_cellType_marker_100 %>% 
    filter( bc_max > 20 ,logFClineage1 > 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


```

--> No correlation observed, despite many filtering tested (only DEG, only TC-DEG, only StartEnd positive, only H3K27me3 signal, combined...)...

--> smooth_peak and H3K27me3 signal metrics seems correct. Now need to select on which gene perform the correlation
    --> Try to use the cell type marker genes (top cell type marker genes of each of our cluster): no correlatio
    --. Try use genes with peak in EZH2: no correlation

--> A lot of test has been done: **no biologically relevant correlation** found (only <0.1...)
    --> Test same with H3K27me3 signal.. EZH2 might be transient
    
--> Only strong pos corr with `logFClineage1 > 1, fdr_DEG <0.05` 1 kbTSS; R 0.42


- **May** improve correlation:
    - Use top or significant cell type marker; use top signif works better; top 200 gave R 0.34 with 140 unique genes
    - Filter logFClineage1 > 0 , to select only genes that increase upon the pseudotime, *not always good!*
    - use bc_max and 1kb TSS (to maximise chance of collecting EZH2 peak)




## Test Activation point with H3K27me3 (from 001*/009*) - xxx

The level of H3K27me3 around TSS has been calculated in `001*/006*` at `# Quantify signal around TSS`. Only 1 Bio Rep for 1st slight test. (`output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt`)

--> If encouraging, could get another Rep from `001*/005*`

xx test select gene peak and increa dc >1


```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj2_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_noCondition_humangastruloid72hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj2_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj2_noCondition_humangastruloid72hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj2_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_traj2_noCondition_humangastruloid72hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_StartEnd" = "fdr")  %>% 
    dplyr::select(geneSymbol, fdr_StartEnd, logFClineage1)
H3K27me3_peaks = read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/ChIPseeker/annotation_macs2_PSC_WT_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt", col_names = FALSE) %>%
    dplyr::rename("geneSymbol" = "X1")

## cell type marker genes (very badly formated)
cellType_marker = read_tsv("../../002_scRNAseq/003__YAP1/output/seurat/srat_all_conserved_markers_V2corr.txt") %>%
    tidyr::separate(gene, into = c("gene", "trash"), sep = "\\.\\.\\.") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    dplyr::select(cluster, geneSymbol, max_pval)
###



pseudotime_traj2_peak_DEG_StartEnd = pseudotime_traj2_peak %>%
    left_join(pseudotime_traj2_DEG) %>%
    left_join(pseudotime_traj2_StartEnd)

# WT #########################
WT_H3K27me3_1kbTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_1kbTSS_geneSymbol.txt")
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_500bpTSS_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt")
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj2_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))


## signal H3K27me3 vs DEG timecourse

# pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdr0__WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_DEGstartEndfdr0001__WT_H3K27me3_1kbTSS_geneSymbol_bc_median.pdf", width=5, height=4)

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth__WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_logFClineageOVer1__WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol %>% 
    filter(logFClineage1 > 1 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 9000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdr_DEG05_logFClineage1pos_WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrDEG05logFClineage1Over1__WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol %>%
    filter(logFClineage1 > 1 , fdr_DEG <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 9000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrDEG05logFClineage1Over1_H3K27me3peaks__WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
H3K27me3_peaks %>% 
    left_join(pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol) %>%
    filter(logFClineage1 > 1 , fdr_DEG <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 9000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()





## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 500) %>%
  ungroup()
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_500__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_500__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_100 %>% 
    filter( ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


#### 250bp TSS
cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 200) %>%
  ungroup()
pseudotime_traj2_peak_WT_H3K27me3_250bpTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_250bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_500__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_200__WT_H3K27me3_250bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_250bpTSS_geneSymbol_cellType_marker_100 %>% 
    filter(  ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



#### 1kb TSS
cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 200) %>%
  ungroup()
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_500__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_200__WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol_cellType_marker_100 %>% 
    filter(   ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

#logFClineage1 > 0.5 , fdr_DEG <0.05




### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_signif05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()




#### 1kb
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_signif05__WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol_cellType_marker_signif %>% 
    filter(  ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_signif05_H3K27me3_peaks__WT_H3K27me3_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
H3K27me3_peaks %>% 
    left_join(pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol_cellType_marker_signif) %>% 
    filter(  ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


## Random testing
pdf("output/binBw/corr_pseudotime_test.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG <0.05,
           logFClineage1 > 0,
           smooth_peak_pseudotime>-1,
           bc_max>0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 1000) %>%
  ungroup()
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

pdf("output/binBw/corr_pseudotime_test.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_1kbTSS_geneSymbol_cellType_marker_100 %>% 
    filter( bc_max > 20 ,logFClineage1 > 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


```


--> No correlation... Except when `logFClineage1 > 1`; increase even more when filtering `bc_max > 100`



