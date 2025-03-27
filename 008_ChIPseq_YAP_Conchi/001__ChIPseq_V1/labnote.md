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





# HOMER peak caller

Let's use the same method as [Conchi previously used](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2635694) to call for peak (notably to identify the enhancer of NODAL for YAP1 ChIPseq)

## Install HOMER 4.11

```bash
conda create -n homer -c bioconda homer # command can be launch from anywhere (directory and node)
#--> WORK! But R deseq2 not present, cannot use all homer tools...

# Try install homer together with deseq2 (to use replicate peak calling)
conda create -n homer_deseq2 -c bioconda homer bioconductor-deseq2 # needed when replicate...
#--> Fail, try clone deseq2 then install homer

conda create --name homer_deseq2 --clone deseq2
conda install bioconda::homer
#--> Fail, try double anaconda

# try stingle step homer and deseq2
conda create -n homer_deseq2 -c bioconda -c conda-forge wget r-essentials bioconductor-deseq2 bioconductor-edger homer
#--> Fail

# try re install homer and then deseq2 in two step
conda create -n homer_deseq2_V1 -c bioconda homer

conda activate homer_deseq2_V1
conda install bioconda::bioconductor-deseq2
#--> Work!

```

## Run HOMER peak caller

From method; use HOMER findPeaks command with input control; Option `-style factor` (for TF; otherwise `-style histone` was used)

- Convert .bam to tagDirectory for [homer](http://homer.ucsd.edu/homer/ngs/peaks.html). Method said only uniquely aligned reads where used, so let's use our `*unique.dupmark.sorted.bam` files
- Call peaks (Method differ if [simplicate](http://homer.ucsd.edu/homer/ngs/peaks.html) or [replicate](http://homer.ucsd.edu/homer/ngs/peaksReplicates.html))





```bash
conda activate homer_deseq2_V1
module load SAMtools/1.16.1-GCC-11.3.0

# Create tagDirectory to be used by homer (FAST, so run in interactive)
makeTagDirectory output/homer/hESC_WT_EZH2_R1 output/bowtie2/hESC_WT_EZH2_R1.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/hESC_WT_EZH2_R2 output/bowtie2/hESC_WT_EZH2_R2.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/hESC_WT_input_R1 output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/hESC_WT_QSER1_R1 output/bowtie2/hESC_WT_QSER1_R1.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/hESC_WT_QSER1_R2 output/bowtie2/hESC_WT_QSER1_R2.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/hESC_YAPKO_EZH2_R1 output/bowtie2/hESC_YAPKO_EZH2_R1.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/hESC_YAPKO_EZH2_R2 output/bowtie2/hESC_YAPKO_EZH2_R2.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/hESC_YAPKO_QSER1_R1 output/bowtie2/hESC_YAPKO_QSER1_R1.unique.dupmark.sorted.bam 
makeTagDirectory output/homer/hESC_YAPKO_QSER1_R2 output/bowtie2/hESC_YAPKO_QSER1_R2.unique.dupmark.sorted.bam 

## 008*/003*
makeTagDirectory ../003__ChIPseq_pluripotency/output/homer/hESC_WT_input_R1 ../003__ChIPseq_pluripotency/output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam 
makeTagDirectory ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1 ../003__ChIPseq_pluripotency/output/bowtie2/hESC_WT_TEAD4_R1.unique.dupmark.sorted.bam 
makeTagDirectory ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1 ../003__ChIPseq_pluripotency/output/bowtie2/hESC_WT_YAP1_R1.unique.dupmark.sorted.bam 

#--> By default homer run in singleend


# PEAK CALLING
## Call peaks with multiple bio rep
getDifferentialPeaksReplicates.pl -all -style factor -t output/homer/hESC_WT_EZH2_R1 output/homer/hESC_WT_EZH2_R2 -i output/homer/hESC_WT_input_R1 > output/homer/hESC_WT_EZH2_outputPeaks.txt 
getDifferentialPeaksReplicates.pl -all -style factor -t output/homer/hESC_WT_QSER1_R1 output/homer/hESC_WT_QSER1_R2 -i output/homer/hESC_WT_input_R1 > output/homer/hESC_WT_QSER1_outputPeaks.txt
getDifferentialPeaksReplicates.pl -all -style factor -t output/homer/hESC_YAPKO_EZH2_R1 output/homer/hESC_YAPKO_EZH2_R2 -i output/homer/hESC_WT_input_R1 > output/homer/hESC_YAPKO_EZH2_outputPeaks.txt 
getDifferentialPeaksReplicates.pl -all -style factor -t output/homer/hESC_YAPKO_QSER1_R1 output/homer/hESC_YAPKO_QSER1_R2 -i output/homer/hESC_WT_input_R1 > output/homer/hESC_YAPKO_QSER1_outputPeaks.txt


## Call peaks with one bio rep (simplicate)
findPeaks output/homer/hESC_WT_EZH2_R1 -style factor -o auto -i output/homer/hESC_WT_input_R1
findPeaks output/homer/hESC_WT_EZH2_R2 -style factor -o auto -i output/homer/hESC_WT_input_R1
findPeaks output/homer/hESC_WT_QSER1_R1 -style factor -o auto -i output/homer/hESC_WT_input_R1
findPeaks output/homer/hESC_WT_QSER1_R2 -style factor -o auto -i output/homer/hESC_WT_input_R1
findPeaks output/homer/hESC_YAPKO_EZH2_R1 -style factor -o auto -i output/homer/hESC_WT_input_R1
findPeaks output/homer/hESC_YAPKO_EZH2_R2 -style factor -o auto -i output/homer/hESC_WT_input_R1
findPeaks output/homer/hESC_YAPKO_QSER1_R1 -style factor -o auto -i output/homer/hESC_WT_input_R1
findPeaks output/homer/hESC_YAPKO_QSER1_R2 -style factor -o auto -i output/homer/hESC_WT_input_R1
## 008*/003*
findPeaks ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1 -style factor -o auto -i ../003__ChIPseq_pluripotency/output/homer/hESC_WT_input_R1
findPeaks ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1 -style factor -o auto -i ../003__ChIPseq_pluripotency/output/homer/hESC_WT_input_R1

## Convert .txt to .bed
### simplicate
pos2bed.pl output/homer/hESC_WT_EZH2_R1/peaks.txt > output/homer/hESC_WT_EZH2_R1/peaks.bed
pos2bed.pl output/homer/hESC_WT_EZH2_R2/peaks.txt > output/homer/hESC_WT_EZH2_R2/peaks.bed
pos2bed.pl output/homer/hESC_WT_QSER1_R1/peaks.txt > output/homer/hESC_WT_QSER1_R1/peaks.bed
pos2bed.pl output/homer/hESC_WT_QSER1_R2/peaks.txt > output/homer/hESC_WT_QSER1_R2/peaks.bed
pos2bed.pl output/homer/hESC_YAPKO_EZH2_R1/peaks.txt > output/homer/hESC_YAPKO_EZH2_R1/peaks.bed
pos2bed.pl output/homer/hESC_YAPKO_EZH2_R2/peaks.txt > output/homer/hESC_YAPKO_EZH2_R2/peaks.bed
pos2bed.pl output/homer/hESC_YAPKO_QSER1_R1/peaks.txt > output/homer/hESC_YAPKO_QSER1_R1/peaks.bed
pos2bed.pl output/homer/hESC_YAPKO_QSER1_R2/peaks.txt > output/homer/hESC_YAPKO_QSER1_R2/peaks.bed
#### 008*/003*
pos2bed.pl ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks.txt > ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks.bed
pos2bed.pl ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks.txt > ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks.bed
### multiple replicate
pos2bed.pl output/homer/hESC_WT_EZH2_outputPeaks.txt > output/homer/hESC_WT_EZH2_outputPeaks.bed
pos2bed.pl output/homer/hESC_WT_QSER1_outputPeaks.txt > output/homer/hESC_WT_QSER1_outputPeaks.bed
pos2bed.pl output/homer/hESC_YAPKO_EZH2_outputPeaks.txt > output/homer/hESC_YAPKO_EZH2_outputPeaks.bed
pos2bed.pl output/homer/hESC_YAPKO_QSER1_outputPeaks.txt > output/homer/hESC_YAPKO_QSER1_outputPeaks.bed

 
 

```



--> All good






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



# Functional analysis with enrichR

- **Up and down reg genes**: THOR diff bound genes 
- **Unique list of genes**: QSER1, EZH2 (YAP1, TEAD4) bound genes in WT (macs2)
- **Up and down reg genes**: THOR diff bound genes for EZH2 and bound with YAP1 in WT
- **QSER1 and YAP1 co-bound genes _ homer**: 1192 genes co-bound by QSER1 and YAP1 (see `gastrulation paper/GastrulationPaper_peakOverlap_V4.pptx`)

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

output/homer/QSER1_YAP1_1192genes.txt

# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/none.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/homer/QSER1_YAP1_1192genes.txt", header=FALSE, stringsAsFactors=FALSE)
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

pdf("output/GO/enrichR_GO_Biological_Process_2023_QSER1_YAP1_1192genes.pdf", width=11, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  scale_fill_manual(name=" ",   # H3K27me3  H3K4me3
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
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_QSER1_YAP1_1192genes.txt", sep="\t", row.names=FALSE, quote=FALSE)





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
output/ChIPseeker/none.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1EnhancerANDORPromoterbinding.txt

output/ChIPseeker/EZH2_lost_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt
output/ChIPseeker/EZH2_pos_annot_promoterAnd5_geneSymbol_YAP1Promoterbinding.txt

output/homer/QSER1_YAP1_1192genes.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/none.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)
## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/homer/QSER1_YAP1_1192genes.txt", header=FALSE, stringsAsFactors=FALSE)
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

pdf("output/GO/enrichR_KEGG_2021_Human_QSER1_YAP1_1192genes.pdf", width=11, height=5)

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
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_QSER1_YAP1_1192genes.txt", sep="\t", row.names=FALSE, quote=FALSE)


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
- QSER1 comparison Conchi vs Dixon (`008002`) - MACS2
- EZH2 comparison Conchi vs Dixon (`008002`) - MACS2
- QSER1, EZH2, YAP1 
- gene with differential EZH2 binding (THORq4)
- QSER1 comparison Conchi vs Dixon (`008002`) - HOMER


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



sed 's/\r$//; s/.*/gene_name "&"/' ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_homer_hESC_WT_QSER1FLAG_promoterAnd5_geneSymbol.txt > ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_homer_hESC_WT_QSER1FLAG_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_homer_hESC_WT_QSER1FLAG_annot_noIntergenic_geneSymbol.txt > ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_homer_hESC_WT_QSER1FLAG_annot_noIntergenic_as_gtf_geneSymbol.txt


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


grep -Ff ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_homer_hESC_WT_QSER1FLAG_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_homer_hESC_WT_QSER1FLAG_promoterAnd5.gtf
grep -Ff ../002__ChIPseq_Dixon2021/output/ChIPseeker/annotation_homer_hESC_WT_QSER1FLAG_annot_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_homer_hESC_WT_QSER1FLAG_annot_noIntergenic.gtf




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




## homer QSER1 008001 vs 008002 with output/bigwig raw files
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008001noIntergenic_homer.sh # 40453678 ok
sbatch scripts/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008001promoterAnd5_homer.sh # 40453743 ok
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
sbatch scripts/matrix_TSS_10kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.sh # 31922267 ok



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


# FAIL WRONG TRAJ  Correlation EZH2 signal with Pseudotime 002*/003*

!! The wrong trajectory as been studied! It should start with Ectorderm, not end by it... !!

## FAIL WRONG TRAJ Test pseudotime_start_end_association() - FAIL


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





## FAIL WRONG TRAJ  Test Activation point with EZH2 - OK



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
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_logFClineageOver1__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol %>% 
    filter(
           logFClineage1 > 1) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_logFClineageOver05fdrdeg05__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_1kbTSS_geneSymbol %>% 
    filter(logFClineage1 > 0.5,
           fdr_DEG <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


## THIS plot is good looking, BTU it's biased because density considers the number of genes at each location in addition to the bc_max values.
#pdf("output/binBw/cov_pseudotime_traj2_peakSmooth_logFClineageOver1fdrdeg05__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=2)
pdf("output/binBw/cov_pseudotime_traj2_peakSmooth_fdrdeg05__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=2)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter( 
           fdr_DEG <0.05) %>%
    ggplot(aes(x = smooth_peak_pseudotime, weight = bc_max, fill = genotype, color = genotype)) +
    geom_density(alpha = 0.5, bw = 2, position = "identity") + 
    scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
    scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
    theme_bw()
dev.off()
########

pdf("output/binBw/histbin4_pseudotime_traj2_peakSmooth_logFClineageOver1fdrdeg05__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=2)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(logFClineage1 > 1, fdr_DEG < 0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = seq(0, max(smooth_peak_pseudotime), by =4), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Histogram of Median bc_max by Pseudotime Bin",
            x = "Pseudotime Bin",
            y = "Median bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()

#pdf("output/binBw/histBinCluster_pseudotime_traj2_peakSmooth_logFClineageOver05fdrdeg05__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=2)
pdf("output/binBw/histBinCluster_pseudotime_traj2_peakSmooth_logFClineageOver1__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=2)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(logFClineage1 > 1, fdr_DEG <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





#pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrStartEnd05logFClineage1Over05___WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_logFClineage1Over1___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(logFClineage1 > 1) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 300, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
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




## FAIL WRONG TRAJ  Test Activation point with H3K27me3 (from 001*/009*) - OK

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



pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_logFClineage1Over1fdrDEG05___WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    filter(logFClineage1 > 1, fdr_DEG <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max), color = "blue") +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 3000, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


pdf("output/binBw/histBinCluster_pseudotime_traj2_peakSmooth_logFClineageOver1__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=2)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    filter(logFClineage1 > 1) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
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


--> No correlation... Except when `logFClineage1 > 1`; even better associated with `fdr_DRG<0.05`, but only 92 genes



## FAIL WRONG TRAJ  Test TC-genotype diff with H3K27me3

In `002*/003*` at `# Heatmap clutering DEGs per traj _ REVISED METHOD` we identified TC-genotype DEG, l2fc2 is good treshold with 516 DEGs. We performed clustering and identified cluster induced early and lately upon DASA treatment.

--> Check the level of EZH2 in hESC in these cluster induced earlier and later. (We expect cluster induced earlier in DASA to have lower level of EZH2 in YAPKO, and cluster induced later in DASA to have higher level of EZH2 in YAPKO)


```R
# packages
library("tidyverse")
library("ggpubr")

# import files

l2fc2_clusterGeneSymbol <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/condRes_traj2_humangastruloid72hrs_l2fc2_clusterGeneSymbold.txt",
                      col_names = TRUE, trim_ws = TRUE)


# WT #########################
WT_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt")
WT_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt")
WT_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt")
l2fc2_clusterGeneSymbol_WT_EZH2_1kbTSS_geneSymbol = l2fc2_clusterGeneSymbol %>%
    left_join(WT_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol = l2fc2_clusterGeneSymbol %>%
    left_join(WT_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_WT_EZH2_250bpTSS_geneSymbol = l2fc2_clusterGeneSymbol %>%
    left_join(WT_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
# YAPKO #########################
YAPKO_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt")
YAPKO_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt")
YAPKO_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj2_peak_YAPKO_EZH2_1kbTSS_geneSymbol = l2fc2_clusterGeneSymbol %>%
    left_join(YAPKO_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol = l2fc2_clusterGeneSymbol %>%
    left_join(YAPKO_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj2_peak_YAPKO_EZH2_250bpTSS_geneSymbol = l2fc2_clusterGeneSymbol %>%
    left_join(YAPKO_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))


# plot WT and YPAKO




pdf("output/binBw/boxplot_condRes_traj2_humangastruloid72hrs_l2fc2_clusterGeneSymbol.pdf", width=10, height=8)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>%
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
ggplot(., aes(x = genotype, y = bc_max, fill = genotype)) +
    geom_boxplot() +
    facet_wrap(~ cluster, scales = "free_y") +
    scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
    theme_bw() +
    theme(strip.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(aes(group = genotype), method = "t.test", label = "p.format", label.y = -0.5)
dev.off()

```

--> All cluster tested show same level of H3K27me3 between WT and YAPKO







## FAIL WRONG TRAJ  Test Activation point condition specific with H3K27me3 (from 001*/009*) - OK

The level of H3K27me3 around TSS has been calculated in `001*/006*` at `# Quantify signal around TSS`. Only 1 Bio Rep for 1st slight test. (`output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt`)

--> If encouraging, could get another Rep from `001*/005*`

Here is the method Conchi proposed. Define pseudotime activation point separately for each condition (done at `### Condiments humangastru72hrs - pseudotime WT and DASA separated` in `002/003`); and check H3K27me3 level. Check whether DASA accelerate H3K27me3-target gene activation (eg. genes activated later in CONTROL will be activated earlier in DASA (higher postiive correlation)); because likely less H3K27me3 in YAPKO at hESC (*true for EZH2, when using logFCLineage > 1*)




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files



pseudotime_traj2_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidUNTREATED72hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj2_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidDASATINIB72hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")
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


pseudotime_traj2_peak_DEG_StartEnd = pseudotime_traj2_peak_UNTREATED %>%
    bind_rows(pseudotime_traj2_peak_DASATINIB) %>%
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


pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_logFClineageOVer1__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(logFClineage1 > 1 ) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 3000, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_logFClineageOVer1fdrDEG05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(logFClineage1 > 1 , fdr_DEG< 0.05) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 3000, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()


## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_100 = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 1000) %>%
  ungroup()
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))



pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_cellType_marker_1000__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_100 %>% 
    filter( ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 2000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()



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


pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 2000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()



### in peaks


pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_H3K27me3peaks__WT_H3K27me3_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
H3K27me3_peaks %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 3000, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()





```


--> No correlation... Except when `logFClineage1 > 1`; even better associated with `fdr_DRG<0.05`, but only 92 genes






## FAIL WRONG TRAJ  Test Activation point condition specific with EZH2 - OK


Here is the method Conchi proposed. Define pseudotime activation point separately for each condition (done at `### Condiments humangastru72hrs - pseudotime WT and DASA separated` in `002/003`); and check H3K27me3 level. Check whether DASA accelerate H3K27me3-target gene activation (eg. genes activated later in CONTROL will be activated earlier in DASA (higher postiive correlation)); because likely less H3K27me3 in YAPKO at hESC (*true for EZH2, when using logFCLineage > 1*)




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files



pseudotime_traj2_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidUNTREATED72hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj2_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidDASATINIB72hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")
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


pseudotime_traj2_peak_DEG_StartEnd = pseudotime_traj2_peak_UNTREATED %>%
    bind_rows(pseudotime_traj2_peak_DASATINIB) %>%
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



## signal H3K27me3 vs DEG timecourse


pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_logFClineageOVer1__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(logFClineage1 > 1 ) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 500, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_logFClineageOVer1fdrDEG05__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(logFClineage1 > 1 , fdr_DEG< 0.05) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 400, aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
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
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_100 = cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))



pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_100 %>% 
    filter( ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 2000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
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


pdf("output/binBw/corr_pseudotime_traj2_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05__WT_EZH2_500bpTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 2000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()


```


--> No correlation... Except when `logFClineage1 > 1`; even better associated with `fdr_DRG<0.05`, but only 92 genes



# V2 CORRECT TRAJ  Correlation EZH2 signal with Pseudotime 002*/003*




## V2 CORRECT TRAJ -Lineage3- Test Activation point with EZH2 - No correlation



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
pseudotime_traj3_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_noCondition_humangastruloid72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj3_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj3_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj3_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj3_noCondition_humangastruloid72hrs_V2.txt") %>% # NULL
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



pseudotime_traj3_peak_DEG_StartEnd = pseudotime_traj3_peak %>%
    left_join(pseudotime_traj3_DEG) %>%
    left_join(pseudotime_traj3_StartEnd)

# WT #########################
WT_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt")
WT_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt")
WT_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj3_peak_WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj3_peak_WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
# YAPKO #########################
YAPKO_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt")
YAPKO_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt")
YAPKO_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj3_peak_YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj3_peak_YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))

## signal EZH2 vs DEG timecourse

# pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdr0__WT_EZH2_1kbTSS_geneSymbol_bc_max.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_fdrDEG0__WT_EZH2_1kbTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj3_peak_WT_EZH2_1kbTSS_geneSymbol %>% 
    filter(
           fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


    filter(logFClineage1 > 0.5,
           fdr_DEG <0.05) 



# corr_pseudotime_traj3_peakSmooth_logFClineage1Over0fdrStartEnd05___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2
pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_fdrDEG0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()



pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_logFClineage1Over0fdrDEG05___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG<0.05, logFClineage1>0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()




pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_logFClineage1Over0fdrDEG05EZH2peaks___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4) 
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG<0.05, logFClineage1>0) %>%
    inner_join(EZH2_peaks) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()






# corr_pseudotime_traj3_peakSmooth_logFClineage1Over0fdrStartEnd05___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2
pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_fdrStartEnd05___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_StartEnd < 0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()




pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_EZH2peaks___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
EZH2_peaks %>%
    left_join(pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol) %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter() %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()





# histBinCluster_pseudotime_traj3_peakSmooth_fdrDEG05__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max
pdf("output/binBw/histBinCluster_pseudotime_traj3_peakSmooth_fdrStartEnd05__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_StartEnd < 0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 2.52, 3.31, 10.9, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()




## signal EZH2 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 200) %>%  # CHANGE HERE !!!!!!!!!!!!!!!
  ungroup()
pseudotime_traj3_peak_WTYAPKO_EZH2_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    add_column(genotype = "WT") %>%
    bind_rows(cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) )

#pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_cellType_marker_500__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_cellType_marker_200_logFClineage1OVer0fdrStartEnd05__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj3_peak_WTYAPKO_EZH2_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter(fdr_StartEnd<0.05, logFClineage1 > 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()





### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj3_peak_WTYAPKO_EZH2_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    add_column(genotype = "WT") %>%
    bind_rows(cellType_marker_100 %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) )


pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_cellType_marker_signif05_logFClineage1OVer0fdrStartEnd05__WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj3_peak_WTYAPKO_EZH2_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           fdr_StartEnd<0.05, logFClineage1 > 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
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




pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_fdrDEG0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = mean(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()






pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_fdrDEG05logFClineage1Over0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG < 0.05, logFClineage1>0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = mean(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()


# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_smoothOver0fdrDEG0__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) + 
    scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()



pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_smoothOver0fdrDEG0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj3_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0)  %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = mean(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()






```





## V2 CORRECT TRAJ -Lineage3-  Test Activation point with H3K27me3 (from 001*/009*) - correlation

The level of H3K27me3 around TSS has been calculated in `001*/006*` at `# Quantify signal around TSS`. Only 1 Bio Rep for 1st slight test. (`output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt`)

--> If encouraging, could get another Rep from `001*/005*`



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj3_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_noCondition_humangastruloid72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj3_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj3_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj3_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj3_noCondition_humangastruloid72hrs_V2.txt") %>%
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



pseudotime_traj3_peak_DEG_StartEnd = pseudotime_traj3_peak %>%
    left_join(pseudotime_traj3_DEG) %>%
    left_join(pseudotime_traj3_StartEnd)

# WT #########################
WT_H3K27me3_1kbTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_1kbTSS_geneSymbol.txt")
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_500bpTSS_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt")
pseudotime_traj3_peak_WT_H3K27me3_1kbTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj3_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))


## signal H3K27me3 vs DEG timecourse

# pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_logFClineageOVer1__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)


pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_logFClineageOVer0fdrDEG05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(logFClineage1 > 0, fdr_DEG <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 5000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG == 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 5000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()





pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_frdStartEnd05_H3K27me3peaks__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
H3K27me3_peaks %>% 
    left_join(pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol) %>%
    filter( fdr_StartEnd <0.05   ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 9000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



pdf("output/binBw/histBinCluster_pseudotime_traj3_peakSmooth_logFClineageOver0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    filter(logFClineage1 > 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 250) %>%
  ungroup()
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_cellType_marker_500_logFClineage1Over0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



pdf("output/binBw/histBinCluster_pseudotime_traj3_peakSmooth_cellType_marker_250_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()








### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.00000000000001)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_cellType_marker_signif00000000000001_fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           fdr_StartEnd <0.05 ) %>%
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







pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_fdrDEG0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_fdrDEG05logFClineage1Over0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG < 0.05, logFClineage1>0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_smoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=4)

pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 3000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()




pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_smoothOver0fdrDEG0___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)

pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    add_column(genotype = "WT") %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0)  %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = mean(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()




 

```



--> correlation... Especially with cell type marker, and together with `logFClineage1 > 0, fdr_StartEnd <0.05`
    --> Not much difference between treatment




## V2 CORRECT TRAJ -Lineage2-  Test Activation point with H3K27me3 (from 001*/009*) - correlation

The level of H3K27me3 around TSS has been calculated in `001*/006*` at `# Quantify signal around TSS`. Only 1 Bio Rep for 1st slight test. (`output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt`)

--> If encouraging, could get another Rep from `001*/005*`



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj2_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_noCondition_humangastruloid72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj2_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj2_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj2_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj2_noCondition_humangastruloid72hrs_V2.txt") %>%
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

# pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_logFClineageOVer1__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)


pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrDEG05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 5000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrStartEnd05logFClineage1OVer0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd <0.05, logFClineage1>0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 3000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()






pdf("output/binBw/histBinCluster_pseudotime_traj2_peakSmooth_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    filter(logFClineage1 > 0, fdr_StartEnd <0.05 ) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 2.62, 3.45, 10.7, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 500) %>%
  ungroup()
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_500_logFClineage1Over0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



pdf("output/binBw/histBinCluster_pseudotime_traj2_peakSmooth_cellType_marker_500_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 2.62, 3.45, 10.7, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()








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


pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_cellType_marker_signif05_logFClineage1Over0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           logFClineage1 > 0, fdr_StartEnd <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



pdf("output/binBw/histBinCluster_pseudotime_traj2_peakSmooth_cellType_marker_signif05_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 2.62, 3.45, 10.7, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()










pdf("output/binBw/histbin3_pseudotime_traj2_peakSmooth_fdrDEG0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



pdf("output/binBw/histbin3_pseudotime_traj2_peakSmooth_fdrDEG05logFClineage1Over0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG < 0.05, logFClineage1>0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()




# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_smoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=4)

pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 3000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()






pdf("output/binBw/histbin3_pseudotime_traj2_peakSmooth_smoothOver0fdrDEG0___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)

pseudotime_traj2_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    add_column(genotype = "WT") %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0)  %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = mean(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()






```



--> correlation... Especially with cell type marker, and together with `logFClineage1 > 0, fdr_StartEnd <0.05`





## V2 CORRECT TRAJ -Lineage2-  Test Activation point with EZH2 - xxx

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
pseudotime_traj2_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_noCondition_humangastruloid72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj2_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj2_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj2_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj2_noCondition_humangastruloid72hrs_V2.txt") %>% # NULL
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

# 

pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrDEG0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()



pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_fdrStartEnd05logFClineage1Over0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_StartEnd <0.05, logFClineage1>0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()








pdf("output/binBw/histbin3_pseudotime_traj2_peakSmooth_fdrDEG0median___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()






pdf("output/binBw/histbin3_pseudotime_traj2_peakSmooth_fdrDEG05logFClineage1Over0median___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG < 0.05, logFClineage1>0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()




# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj2_peakSmooth_smoothOver0fdrDEG0__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) + 
    scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()




pdf("output/binBw/histbin3_pseudotime_traj2_peakSmooth_smoothOver0fdrDEG0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj2_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj2_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0)  %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = mean(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





```







## V2 CORRECT TRAJ -Lineage1-  Test Activation point with EZH2 - xxx

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
pseudotime_traj1_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj1_noCondition_humangastruloid72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj1_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj1_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj1_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj1_noCondition_humangastruloid72hrs_V2.txt") %>% # NULL
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



pseudotime_traj1_peak_DEG_StartEnd = pseudotime_traj1_peak %>%
    left_join(pseudotime_traj1_DEG) %>%
    left_join(pseudotime_traj1_StartEnd)

# WT #########################
WT_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt")
WT_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt")
WT_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj1_peak_WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj1_peak_WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj1_peak_WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
# YAPKO #########################
YAPKO_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt")
YAPKO_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt")
YAPKO_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj1_peak_YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj1_peak_YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj1_peak_YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))

# 

pdf("output/binBw/corr_pseudotime_traj1_peakSmooth_fdrDEG0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj1_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj1_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()



pdf("output/binBw/corr_pseudotime_traj1_peakSmooth_logFClineage1Over0fdrDEG05___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj1_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj1_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG<0.05, logFClineage1>0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) +
  scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()








pdf("output/binBw/histbin3_pseudotime_traj1_peakSmooth_fdrDEG0median___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj1_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj1_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



pdf("output/binBw/histbin3_pseudotime_traj1_peakSmooth_fdrDEG05logFClineage1Over0median___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj1_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj1_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    filter(fdr_DEG < 0.05, logFClineage1>0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj1_peak_WT_EZH2_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj1_peakSmooth_smoothOver0fdrDEG0__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj1_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj1_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = genotype)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = genotype)) + 
    scale_color_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
  theme_bw()
dev.off()





pdf("output/binBw/histbin3_pseudotime_traj1_peakSmooth_smoothOver0fdrDEG0___WTYAPKO_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj1_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    bind_rows(pseudotime_traj1_peak_YAPKO_EZH2_500bpTSS_geneSymbol %>% add_column(genotype = "YAPKO")) %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0)  %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = mean(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue", "YAPKO" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()




```










## V2 CORRECT TRAJ -Lineage1-  Test Activation point with H3K27me3 - xxx


The level of H3K27me3 around TSS has been calculated in `001*/006*` at `# Quantify signal around TSS`. Only 1 Bio Rep for 1st slight test. (`output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt`)

--> If encouraging, could get another Rep from `001*/005*`




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj1_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj1_noCondition_humangastruloid72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj1_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj1_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj1_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj1_noCondition_humangastruloid72hrs_V2.txt") %>% # NULL
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



pseudotime_traj1_peak_DEG_StartEnd = pseudotime_traj1_peak %>%
    left_join(pseudotime_traj1_DEG) %>%
    left_join(pseudotime_traj1_StartEnd)




# WT #########################
WT_H3K27me3_1kbTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_1kbTSS_geneSymbol.txt")
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_500bpTSS_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt")
pseudotime_traj1_peak_WT_H3K27me3_1kbTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj1_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj1_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj1_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))




# 
pdf("output/binBw/corr_pseudotime_traj1_peakSmooth_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=4)
pseudotime_traj1_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG == 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 3000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()


# 
pdf("output/binBw/corr_pseudotime_traj1_peakSmooth_fdrStartEnd05logFClineage1Over0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=4)
pseudotime_traj1_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd <0.05, logFClineage1>0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 3000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()





pdf("output/binBw/histbin3_pseudotime_traj1_peakSmooth_fdrDEG0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj1_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



pdf("output/binBw/histbin3_pseudotime_traj1_peakSmooth_fdrDEG05logFClineage1Over0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)
pseudotime_traj1_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG < 0.05, logFClineage1>0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj1_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()




pdf("output/binBw/corr_pseudotime_traj1_peakSmooth_smoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=4)

pseudotime_traj1_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x = 0, label.y = 3000, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



pdf("output/binBw/histbin3_pseudotime_traj1_peakSmooth_smoothOver0fdrDEG0___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=3, height=2)

pseudotime_traj1_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    add_column(genotype = "WT") %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0)  %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = mean(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





```









## V2 CORRECT TRAJ -UNT_Lineage2-DAS_Lineage3-  Test Activation point condition specific with H3K27me3 (from 001*/009*) - correlation 

The level of H3K27me3 around TSS has been calculated in `001*/006*` at `# Quantify signal around TSS`. Only 1 Bio Rep for 1st slight test. (`output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt`)

--> If encouraging, could get another Rep from `001*/005*`

Here is the method Conchi proposed. Define pseudotime activation point separately for each condition (done at `### Condiments humangastru72hrs - pseudotime WT and DASA separated` in `002/003`); and check H3K27me3 level. Check whether DASA accelerate H3K27me3-target gene activation (eg. genes activated later in CONTROL will be activated earlier in DASA (higher postiive correlation)); because likely less H3K27me3 in YAPKO at hESC (*true for EZH2, when using logFCLineage > 1*)



**UNT_Traj2 and DAS_Traj3 = common traj3**


```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files



pseudotime_traj2_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidUNTREATED72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj3_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_humangastruloidDASATINIB72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")

pseudotime_traj3_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj3_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj3_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj3_noCondition_humangastruloid72hrs_V2.txt") %>% # NULL
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


pseudotime_traj23_peak_DEG_StartEnd = pseudotime_traj2_peak_UNTREATED %>%
    bind_rows(pseudotime_traj3_peak_DASATINIB) %>%
    left_join(pseudotime_traj3_DEG) %>%
    left_join(pseudotime_traj3_StartEnd)

# WT #########################
WT_H3K27me3_1kbTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_1kbTSS_geneSymbol.txt")
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_500bpTSS_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt")
pseudotime_traj23_peak_WT_H3K27me3_1kbTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj23_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))


## signal H3K27me3 vs DEG timecourse


pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_logFClineage1Over0fdrStartEnd0000000001__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(logFClineage1 > 0 , fdr_StartEnd<0.0000000001) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_fdrStartEnd05logFClineage1Over0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd < 0.05, logFClineage1>0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()



## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 500) %>%
  ungroup()
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))



pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_cellType_marker_100__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter(  logFClineage1 > 0 , fdr_StartEnd <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()





pdf("output/binBw/histBinCluster_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()








### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineage1Over0fdrStartEnd001__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           logFClineage1 > 0 , fdr_StartEnd <0.001 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()







pdf("output/binBw/histBinCluster_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





### in peaks


pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_H3K27me3peaks__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
H3K27me3_peaks %>%
    left_join(pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()



# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_smoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()



```


--> correlation... Especially with cell type marker, and together with `logFClineage1 > 0, fdr_StartEnd <0.05`
    --> Not much difference between treatment


## V2 CORRECT TRAJ -UNT_Lineage3-DAS_Lineage1-  Test Activation point condition specific with H3K27me3 (from 001*/009*) - correlation 

The level of H3K27me3 around TSS has been calculated in `001*/006*` at `# Quantify signal around TSS`. Only 1 Bio Rep for 1st slight test. (`output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt`)

--> If encouraging, could get another Rep from `001*/005*`

Here is the method Conchi proposed. Define pseudotime activation point separately for each condition (done at `### Condiments humangastru72hrs - pseudotime WT and DASA separated` in `002/003`); and check H3K27me3 level. Check whether DASA accelerate H3K27me3-target gene activation (eg. genes activated later in CONTROL will be activated earlier in DASA (higher postiive correlation)); because likely less H3K27me3 in YAPKO at hESC (*true for EZH2, when using logFCLineage > 1*)



**UNT_Traj3 and DAS_Traj1 = common traj2**


```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj3_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_humangastruloidUNTREATED72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj1_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj1_humangastruloidDASATINIB72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")

pseudotime_traj2_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj2_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj2_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj2_noCondition_humangastruloid72hrs_V2.txt") %>%   # NULL
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


pseudotime_traj31_peak_DEG_StartEnd = pseudotime_traj3_peak_UNTREATED %>%
    bind_rows(pseudotime_traj1_peak_DASATINIB) %>%
    left_join(pseudotime_traj2_DEG) %>%
    left_join(pseudotime_traj2_StartEnd)

# WT #########################
WT_H3K27me3_1kbTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_1kbTSS_geneSymbol.txt")
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_500bpTSS_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt")
pseudotime_traj31_peak_WT_H3K27me3_1kbTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj31_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))


## signal H3K27me3 vs DEG timecourse

#corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_logFClineage1Over0fdrStartEnd00001__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2

pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_logFClineage1Over0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_StartEnd < 0.05, logFClineage1 > 0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_DEG == 0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()



# Ectoderm/Mesoderm_4/Mesoderm_3/Mesoderm_2
pdf("output/binBw/histBinCluster_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_logFClineageOver0fdrStartEnd00001__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.00001) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 2.62,3.45,10.7, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 500) %>%
  ungroup()
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))



pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500_fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter(   fdr_DEG <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()





pdf("output/binBw/histBinCluster_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 2.62,3.45,10.7, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()








### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineage1Over0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           logFClineage1 > 0 , fdr_StartEnd <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()







pdf("output/binBw/histBinCluster_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





### in peaks


pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_H3K27me3peaks__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
H3K27me3_peaks %>%
    left_join(pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()




# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_smoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj31_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()





```


--> correlation... Especially with cell type marker, and together with `logFClineage1 > 0, fdr_StartEnd <0.05`


## V2 CORRECT TRAJ -UNT_Lineage1-DAS_Lineage2-  Test Activation point condition specific with H3K27me3 (from 001*/009*) - correlation 

The level of H3K27me3 around TSS has been calculated in `001*/006*` at `# Quantify signal around TSS`. Only 1 Bio Rep for 1st slight test. (`output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt`)

--> If encouraging, could get another Rep from `001*/005*`

Here is the method Conchi proposed. Define pseudotime activation point separately for each condition (done at `### Condiments humangastru72hrs - pseudotime WT and DASA separated` in `002/003`); and check H3K27me3 level. Check whether DASA accelerate H3K27me3-target gene activation (eg. genes activated later in CONTROL will be activated earlier in DASA (higher postiive correlation)); because likely less H3K27me3 in YAPKO at hESC (*true for EZH2, when using logFCLineage > 1*)



**UNT_Traj1 and DAS_Traj2 = common traj1**


```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj1_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj1_humangastruloidUNTREATED72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj2_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidDASATINIB72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")

pseudotime_traj1_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj1_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj1_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj1_noCondition_humangastruloid72hrs_V2.txt") %>% # NULL added !!!!!
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


pseudotime_traj12_peak_DEG_StartEnd = pseudotime_traj1_peak_UNTREATED %>%
    bind_rows(pseudotime_traj2_peak_DASATINIB) %>%
    left_join(pseudotime_traj1_DEG) %>%
    left_join(pseudotime_traj1_StartEnd)

# WT #########################
WT_H3K27me3_1kbTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_1kbTSS_geneSymbol.txt")
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_500bpTSS_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("../../001_EZH1_Project/006__CutRun_PSC_FA/output/binBw/WT_H3K27me3_250bpTSS_geneSymbol.txt")
pseudotime_traj12_peak_WT_H3K27me3_1kbTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj12_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))


## signal H3K27me3 vs DEG timecourse

#corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_logFClineage1Over0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2

pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_DEG==0 ) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_fdr_DEG05logFClineage1Over0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_DEG<0.05 , logFClineage1 > 0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()




pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_H3K27me3_peaks_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
H3K27me3_peaks %>%
  left_join(pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol) %>% 
    filter(fdr_DEG==0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()





# Ectoderm/Mesoderm_1/Mesoderm_2
pdf("output/binBw/histBin3_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    filter( fdr_DEG==0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = seq(0, max(smooth_peak_pseudotime), by =3), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_1", "Ectoderm")) %>%  # HERE MARKER CHANGED ACCORDINGLY
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 1000) %>%
  ungroup()
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))



pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_cellType_marker_1000_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter(  fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()





pdf("output/binBw/histBinCluster_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_n %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0,1.04, 11.1, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()








### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_1", "Ectoderm")) %>%  # HERE MARKER CHANGED ACCORDINGLY
  filter(max_pval<0.001)
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif001_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           fdr_DEG == 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()







pdf("output/binBw/histBinCluster_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0,1.04, 11.1, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()




# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_smoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj12_peak_WT_H3K27me3_500bpTSS_geneSymbol %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()




```


--> correlation... Especially with cell type marker, and together with `logFClineage1 > 0, fdr_StartEnd <0.05`







## V2 CORRECT TRAJ -UNT_Lineage2-DAS_Lineage3-  Test Activation point condition specific with EZH2 - No correlation


Here is the method Conchi proposed. Define pseudotime activation point separately for each condition (done at `### Condiments humangastru72hrs - pseudotime WT and DASA separated` in `002/003`); and check H3K27me3 level. Check whether DASA accelerate H3K27me3-target gene activation (eg. genes activated later in CONTROL will be activated earlier in DASA (higher postiive correlation)); because likely less H3K27me3 in YAPKO at hESC (*true for EZH2, when using logFCLineage > 1*)

**UNT_Lineage2-DAS_Lineage3 = common lineage3**


```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files

pseudotime_traj2_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidUNTREATED72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj3_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_humangastruloidDASATINIB72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")

pseudotime_traj3_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj3_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj3_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj3_noCondition_humangastruloid72hrs_V2.txt") %>%
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


pseudotime_traj23_peak_DEG_StartEnd = pseudotime_traj2_peak_UNTREATED %>%
    bind_rows(pseudotime_traj3_peak_DASATINIB) %>%
    left_join(pseudotime_traj3_DEG) %>%
    left_join(pseudotime_traj3_StartEnd)


# WT #########################
WT_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt")
WT_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt")
WT_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj23_peak_WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj23_peak_WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
# YAPKO #########################
YAPKO_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt")
YAPKO_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt")
YAPKO_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj23_peak_YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj23_peak_YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj23_peak_YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj23_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))




## signal H3K27me3 vs DEG timecourse


pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_fdrDEG0__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG == 0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_logFClineage1Over0fdrDEG05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(logFClineage1 > 0 , fdr_DEG<0.05) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()








## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 250) %>%
  ungroup()
pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))



pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_cellType_marker_250_logFClineageOver0fdrStartEnd05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter(  logFClineage1 > 0 , fdr_StartEnd <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()





pdf("output/binBw/histBinCluster_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_cellType_marker_250_logFClineageOver0fdrStartEnd05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()








### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineage1Over0fdrStartEnd05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           logFClineage1 > 0 , fdr_StartEnd <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()







pdf("output/binBw/histBinCluster_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol_cellType_marker_signif %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0, 5.05, 12.5, 13.2, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





### in peaks


pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_H3K27me3peaks__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
H3K27me3_peaks %>%
    left_join(pseudotime_traj23_peak_WT_H3K27me3_500bpTSS_geneSymbol) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()



# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj23_UNTREATEDDASATINIB_peakSmooth_smoothOver0fdrDEG0__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj23_peak_WT_EZH2_500bpTSS_geneSymbol %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()




```


--> No correlation... 






## V2 CORRECT TRAJ -UNT_Lineage3-DAS_Lineage1-  Test Activation point condition specific with EZH2 -  xxx



Here is the method Conchi proposed. Define pseudotime activation point separately for each condition (done at `### Condiments humangastru72hrs - pseudotime WT and DASA separated` in `002/003`); and check H3K27me3 level. Check whether DASA accelerate H3K27me3-target gene activation (eg. genes activated later in CONTROL will be activated earlier in DASA (higher postiive correlation)); because likely less H3K27me3 in YAPKO at hESC (*true for EZH2, when using logFCLineage > 1*)

**UNT_Lineage3-DAS_Lineage1 = common lineage2**


```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files

pseudotime_traj3_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_humangastruloidUNTREATED72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj1_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj1_humangastruloidDASATINIB72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")

pseudotime_traj2_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj2_noCondition_humangastruloid72hrs_V2.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj2_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj2_noCondition_humangastruloid72hrs_V2.txt") %>% # NULL
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


pseudotime_traj31_peak_DEG_StartEnd = pseudotime_traj3_peak_UNTREATED %>%
    bind_rows(pseudotime_traj1_peak_DASATINIB) %>%
    left_join(pseudotime_traj2_DEG) %>%
    left_join(pseudotime_traj2_StartEnd)


# WT #########################
WT_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt")
WT_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt")
WT_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj31_peak_WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj31_peak_WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
# YAPKO #########################
YAPKO_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt")
YAPKO_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt")
YAPKO_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj31_peak_YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj31_peak_YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj31_peak_YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj31_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))




## signal H3K27me3 vs DEG timecourse


pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_logFClineage1Over0fdrDEG05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG <0.05, logFClineage1>0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()





pdf("output/binBw/histBinCluster_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_logFClineageOver0fdrStartEnd00001__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.00001) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0,2.62,3.45,10.7, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()




## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 500) %>%
  ungroup()
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))



pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500_logFClineageOver0fdrStartEnd05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter(  logFClineage1 > 0 , fdr_StartEnd <0.05) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()





pdf("output/binBw/histBinCluster_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500_logFClineageOver0fdrStartEnd05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0,2.62,3.45,10.7, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()








### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_3", "Mesoderm_4", "Ectoderm")) %>%
  filter(max_pval<0.05)
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineage1Over0fdrStartEnd05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           logFClineage1 > 0 , fdr_StartEnd <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()







pdf("output/binBw/histBinCluster_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0,2.62,3.45,10.7, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()


# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj31_UNTREATEDDASATINIB_peakSmooth_smoothOver0fdrDEG0__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj31_peak_WT_EZH2_500bpTSS_geneSymbol %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()




```


--> Slight correlation ~0.2!











## V2 CORRECT TRAJ -UNT_Lineage1-DAS_Lineage2-  Test Activation point condition specific with EZH2 - Correlation


Here is the method Conchi proposed. Define pseudotime activation point separately for each condition (done at `### Condiments humangastru72hrs - pseudotime WT and DASA separated` in `002/003`); and check H3K27me3 level. Check whether DASA accelerate H3K27me3-target gene activation (eg. genes activated later in CONTROL will be activated earlier in DASA (higher postiive correlation)); because likely less H3K27me3 in YAPKO at hESC (*true for EZH2, when using logFCLineage > 1*)

**UNT_Lineage1-DAS_Lineage2 = common lineage1**


```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files

pseudotime_traj1_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj1_humangastruloidUNTREATED72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj2_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidDASATINIB72hrs_V2_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")

pseudotime_traj1_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj1_noCondition_humangastruloid72hrs_V2.txt") %>% 
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj1_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj1_noCondition_humangastruloid72hrs_V2.txt") %>%  # NULL
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


pseudotime_traj12_peak_DEG_StartEnd = pseudotime_traj1_peak_UNTREATED %>%
    bind_rows(pseudotime_traj2_peak_DASATINIB) %>%
    left_join(pseudotime_traj1_DEG) %>%
    left_join(pseudotime_traj1_StartEnd)


# WT #########################
WT_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_1kbTSS_geneSymbol.txt")
WT_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_500bpTSS_geneSymbol.txt")
WT_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj12_peak_WT_EZH2_1kbTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj12_peak_WT_EZH2_250bpTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(WT_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
# YAPKO #########################
YAPKO_EZH2_1kbTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_1kbTSS_geneSymbol.txt")
YAPKO_EZH2_500bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_500bpTSS_geneSymbol.txt")
YAPKO_EZH2_250bpTSS_geneSymbol <- read_tsv("output/binBw/YAPKO_EZH2_250bpTSS_geneSymbol.txt")
pseudotime_traj12_peak_YAPKO_EZH2_1kbTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_1kbTSS_geneSymbol) %>%
  filter(!is.na(bc_median))
pseudotime_traj12_peak_YAPKO_EZH2_500bpTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_median))
pseudotime_traj12_peak_YAPKO_EZH2_250bpTSS_geneSymbol = pseudotime_traj12_peak_DEG_StartEnd %>%
    left_join(YAPKO_EZH2_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_median))




## signal H3K27me3 vs DEG timecourse


pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_fdrDEG0__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG==0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()


pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_logFClineage1Over0fdrDEG05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol %>% 
    filter(logFClineage1>0, fdr_DEG< 0.05) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()





pdf("output/binBw/histBinCluster_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_logFClineageOver0fdrStartEnd00001__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.00001) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0,1.04, 11.1, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





## signal H3K27me3 vs marker genes
### Isolate top n marker genes for each cluster
#### 500bp TSS
cellType_marker_n = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_1", "Ectoderm")) %>%  # HERE MARKER CHANGED ACCORDINGLY
  group_by(cluster) %>%
  arrange(max_pval) %>%
  slice_head(n = 500) %>%
  ungroup()
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n = cellType_marker_n %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime))



pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500_logFClineageOver0fdrStartEnd00001__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n %>% 
    filter(  logFClineage1 > 0 , fdr_StartEnd <0.00001) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()





pdf("output/binBw/histBinCluster_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_cellType_marker_500_logFClineageOver0fdrStartEnd05__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_n %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0,1.04, 11.1, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()








### Select all signif marker genes
#### 500bp
cellType_marker_signif = cellType_marker %>%
  filter(cluster %in% c("Mesoderm_2", "Mesoderm_1", "Ectoderm")) %>%  # HERE MARKER CHANGED ACCORDINGLY
  filter(max_pval<0.05)
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif = cellType_marker_signif %>%
    dplyr::select(geneSymbol) %>%
    unique() %>%
    left_join(pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol)%>%
    filter(!is.na(smooth_peak_pseudotime)) 


pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineage1Over0fdrStartEnd00001__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif %>% 
    filter( 
           logFClineage1 > 0 , fdr_StartEnd <0.00001 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()







pdf("output/binBw/histBinCluster_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_cellType_marker_signif05_logFClineageOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=2)
pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol_cellType_marker_signif %>%
    filter( logFClineage1 > 0, fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(0,1.04, 11.1, max(smooth_peak_pseudotime)), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(median_bc_max = median(bc_max), .groups = 'drop') %>%
        ggplot(., aes(x = pseudotime_bin, y = median_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()





# remove smooth_peak 0 if present in one genotype
genesToRemove = pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol %>%
  filter(smooth_peak_pseudotime ==0) %>%
  dplyr::select(geneSymbol) %>%
  unique()



pdf("output/binBw/corr_pseudotime_traj12_UNTREATEDDASATINIB_peakSmooth_smoothOver0fdrDEG0__WT_EZH2_500bpTSS_geneSymbol_bc_max_V2.pdf", width=5, height=4)

pseudotime_traj12_peak_WT_EZH2_500bpTSS_geneSymbol %>%
    anti_join(genesToRemove) %>%
    filter( fdr_DEG == 0) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max,  color = condition)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
  theme_bw()
dev.off()







```


--> Slight correlation ~0.2!





# 2472hrs human gastruloid integration; correlation H3K27me3 with pseudotime 002*/003*

## ENCODE H3K27me3 file selection

Different H3K27me3 in H9 in ENCODE

- [ENCODE H9 H3K27me3](https://www.encodeproject.org/search/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=Histone+ChIP-seq&assay_title=Mint-ChIP-seq&status=released&target.label=H3K27me3&biosample_ontology.classification=cell+line&biosample_ontology.term_name=H9)
  - Mint-ChIP-seq from Bernstein lab; three different, sequenced with different machine:
    - [Illumina HiSeq 2500 with 3 Bio Rep](https://www.encodeproject.org/experiments/ENCSR725YWL/); done in 2021
    - [Illumina HiSeq 2500 with 2 Bio Rep](https://www.encodeproject.org/experiments/ENCSR016RFN/); done in 2021
    - [Illumina NextSeq 500 with 2 Bio Rep](https://www.encodeproject.org/experiments/ENCSR373TMA/); done in 2021
  - [Histone-ChIP-seq from Bernstein lab](https://www.encodeproject.org/experiments/ENCSR792GCH/); done in 2013


--> Let's prefer Mint-ChIP-seq. Let's see which bigwig file the best to use: *signal p-value* or *foldchange over control*?
  --> ENCODE bigwig files downloaded into `output/bigwig_ENCODE`

ENCFF043GTQ_signalpvalue123	
ENCFF201SJZ_fcovercontrol123

ENCFF850HXT_signalpvalue12
ENCFF130PLP_fcovercontrol12

ENCFF140RHC_signalpvalue12
ENCFF077DRH_fcovercontrol12


--> fcovercontrol look more homogeneous, signalpvalue show sharp very high peaks... Let's test fcovercontrol 1st but may be worth trying both...

--> Use: Mint-ChIP-seq fcovercontrol with 3 Bio Rep = `output/bigwig_ENCODE/ENCFF201SJZ_fcovercontrol123`


### Count signal around TSS ENCODE H3K27me3




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

# H3K27me3 250bp #####################################################

## Mint-ChIP-seq fcovercontrol with 3 Bio Rep = `output/bigwig_ENCODE/ENCFF201SJZ_fcovercontrol123`
bwFile <- "output/bigwig_ENCODE/ENCFF201SJZ.bigWig"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_250bpTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ")
#### Collect output and gene information
WT_H3K27me3_250bp_ENCFF201SJZ <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- unique(WT_H3K27me3_250bp_ENCFF201SJZ$gene)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)

WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol = WT_H3K27me3_250bp_ENCFF201SJZ %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    filter(geneSymbol != "") %>%
    mutate(bc_max = max(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_max) %>%
    unique()

write.table(WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol, file = "output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# H3K27me3 500bp #####################################################



## Mint-ChIP-seq fcovercontrol with 3 Bio Rep = `output/bigwig_ENCODE/ENCFF201SJZ_fcovercontrol123`
bwFile <- "output/bigwig_ENCODE/ENCFF201SJZ.bigWig"
regions <- read.table("../../001_EZH1_Project/003__CutRun/meta/ENCFF159KBI_gene_500bpTSS_sorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "gene", "strand")) 
counts <- bin.bw(bwFile, regions, outfile.prefix = "output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ")
#### Collect output and gene information
WT_H3K27me3_500bp_ENCFF201SJZ <- as_tibble(fread(cmd = "gunzip -c output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ.bgz") ) %>%
 left_join(regions) %>%
 separate(gene, into = c("gene", "version"), sep = "\\.") %>%
 dplyr::select(gene, bc) %>%
 unique()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- unique(WT_H3K27me3_500bp_ENCFF201SJZ$gene)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = gene_ids,
                   mart = ensembl)

WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol = WT_H3K27me3_500bp_ENCFF201SJZ %>%
    left_join(gene_info %>% dplyr::rename("gene" = "ensembl_gene_id",  "geneSymbol" = "hgnc_symbol")) %>%
    dplyr::select(geneSymbol, bc) %>%
    unique() %>%
    group_by(geneSymbol) %>%
    mutate(bc_max = max(bc)) %>%
    filter(geneSymbol != "NA") %>%
    dplyr::select(geneSymbol, bc_max) %>%
    unique()

write.table(WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol, file = "output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt", sep = "\t", row.names = FALSE, quote = FALSE)


```

--> Looks good, only bc_max generated




## Lineage3-  Test Activation point with H3K27me3 (from ENCODE) - correlation

The level of H3K27me3 around TSS has been calculated here at `### Count signal around TSS ENCODE H3K27me3`. 



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj3_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_noCondition_humangastruloid2472hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj3_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj3_noCondition_humangastruloid2472hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj3_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj3_noCondition_humangastruloid2472hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_StartEnd" = "fdr")  %>% 
    dplyr::select(geneSymbol, fdr_StartEnd, logFClineage1)





pseudotime_traj3_peak_DEG_StartEnd = pseudotime_traj3_peak %>%
    left_join(pseudotime_traj3_DEG) %>%
    left_join(pseudotime_traj3_StartEnd)

# H3K27me3 ENCODE #########################

WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt")

pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_max))
pseudotime_traj3_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_max))


## signal H3K27me3 vs DEG timecourse

# DEG
pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_fdrDEG05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_frdStartEnd0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd == 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()




# DEG and smoothpeak >0
pdf("output/binBw/corr_pseudotime_traj3_peakSmoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG == 0,
           smooth_peak_pseudotime > 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_peakSmoothOver0frdStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd <0.05 ,
           smooth_peak_pseudotime > 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



# bin
## DEG
pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_fdrDEG05median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()

## StartEnd
pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_fdrStartEnd0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_StartEnd == 0 ) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()


```


--> no correlation; not even when filtering smoothpeak>0


## Lineage5-  Test Activation point with H3K27me3 (from ENCODE) - correlation

The level of H3K27me3 around TSS has been calculated here at `### Count signal around TSS ENCODE H3K27me3`. 



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


# import files
pseudotime_traj5_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj5_noCondition_humangastruloid2472hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj5_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj5_noCondition_humangastruloid2472hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj5_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj5_noCondition_humangastruloid2472hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_StartEnd" = "fdr")  %>% 
    dplyr::select(geneSymbol, fdr_StartEnd, logFClineage1)





pseudotime_traj5_peak_DEG_StartEnd = pseudotime_traj5_peak %>%
    left_join(pseudotime_traj5_DEG) %>%
    left_join(pseudotime_traj5_StartEnd)

# H3K27me3 ENCODE #########################

WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt")

pseudotime_traj5_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj5_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_max))
pseudotime_traj5_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj5_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_max))


## signal H3K27me3 vs DEG timecourse
### DEG
pdf("output/binBw/corr_pseudotime_traj5_peakSmooth_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj5_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG == 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 500, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

### StartEnd
pdf("output/binBw/corr_pseudotime_traj5_peakSmooth_frdStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj5_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 400, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()




# DEG and smoothpeak >0
pdf("output/binBw/corr_pseudotime_traj5_peakSmoothOver0fdrDEG05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj5_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG <0.05,
           smooth_peak_pseudotime > 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj5_peakSmooth_peakSmoothOver0frdStartEnd0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj5_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd == 0 ,
           smooth_peak_pseudotime > 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



# bin
## DEG
pdf("output/binBw/histbin3_pseudotime_traj5_peakSmooth_fdrDEG0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=2)
pseudotime_traj5_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()

## StartEnd
pdf("output/binBw/histbin3_pseudotime_traj5_peakSmooth_fdrStartEnd05median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=2)
pseudotime_traj5_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_StartEnd <0.05 ) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()


```


--> no correlation; not even when filtering smoothpeak>0



## UNT_Lineage3-DAS_Lineage2 (Common3) -  Test Activation point with H3K27me3 (from ENCODE) - correlation

traj of interest COMMON=3 UNTREATED=traj3; DASATINIB=traj2



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")



# import files
pseudotime_traj3_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_humangastruloidUNTREATED2472hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj2_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj2_humangastruloidDASATINIB2472hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")

pseudotime_traj3_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj3_noCondition_humangastruloid2472hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj3_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj3_noCondition_humangastruloid2472hrs.txt") %>% # NULL added !!!!!
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_StartEnd" = "fdr")  %>% 
    dplyr::select(geneSymbol, fdr_StartEnd, logFClineage1)



pseudotime_traj32_peak_DEG_StartEnd = pseudotime_traj3_peak_UNTREATED %>%
    bind_rows(pseudotime_traj2_peak_DASATINIB) %>%
    left_join(pseudotime_traj3_DEG) %>%
    left_join(pseudotime_traj3_StartEnd)

# H3K27me3 ENCODE #########################

WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt")

pseudotime_traj32_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj32_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_max))
pseudotime_traj32_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj32_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_max))




## signal H3K27me3 vs DEG timecourse
### DEG

pdf("output/binBw/corr_pseudotime_traj32_UNTREATEDDASATINIB_peakSmooth_fdrDEG05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=5, height=4)
pseudotime_traj32_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_DEG <0.05 ) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()



### StartEnd

pdf("output/binBw/corr_pseudotime_traj32_UNTREATEDDASATINIB_peakSmooth_fdrStartEnd0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=5, height=4)
pseudotime_traj32_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_StartEnd ==0 ) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()






# DEG and smoothpeak >0

pdf("output/binBw/corr_pseudotime_traj32_UNTREATEDDASATINIB_peakSmooth_peakSmoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=5, height=4)
pseudotime_traj32_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_DEG == 0 ,
            smooth_peak_pseudotime > 0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()



### StartEnd and smoothpeak >0

pdf("output/binBw/corr_pseudotime_traj32_UNTREATEDDASATINIB_peakSmooth_peakSmoothOver0fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=5, height=4)
pseudotime_traj32_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_StartEnd <0.05 ,
            smooth_peak_pseudotime > 0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()






# bin
## DEG
pdf("output/binBw/histbin3_pseudotime_traj32_UNTREATEDDASATINIB_peakSmooth_fdrDEG05median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=6, height=4)
pseudotime_traj32_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG < 0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue","DASATINIB72hrs" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



## StartEnd
pdf("output/binBw/histbin3_pseudotime_traj32_UNTREATEDDASATINIB_peakSmooth_fdrStartEnd0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=6, height=4)
pseudotime_traj32_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue","DASATINIB72hrs" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



```


--> no correlation; not even when filtering smoothpeak>0





## UNT_Lineage5-DAS_Lineage5 (Common5) -  Test Activation point with H3K27me3 (from ENCODE) - correlation

traj of interest COMMON=5 UNTREATED=traj5; DASATINIB=traj5



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")



# import files
pseudotime_traj5_peak_UNTREATED <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj5_humangastruloidUNTREATED2472hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") %>%
    add_column(condition = "UNTREATED72hrs")
pseudotime_traj5_peak_DASATINIB <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj5_humangastruloidDASATINIB2472hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene")%>%
    add_column(condition = "DASATINIB72hrs")

pseudotime_traj5_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj5_noCondition_humangastruloid2472hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj5_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj5_noCondition_humangastruloid2472hrs.txt") %>% # NULL added !!!!!
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_StartEnd" = "fdr")  %>% 
    dplyr::select(geneSymbol, fdr_StartEnd, logFClineage1)



pseudotime_traj55_peak_DEG_StartEnd = pseudotime_traj5_peak_UNTREATED %>%
    bind_rows(pseudotime_traj5_peak_DASATINIB) %>%
    left_join(pseudotime_traj5_DEG) %>%
    left_join(pseudotime_traj5_StartEnd)

# H3K27me3 ENCODE #########################

WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt")

pseudotime_traj55_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj55_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_max))
pseudotime_traj55_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj55_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_max))




## signal H3K27me3 vs DEG timecourse
### DEG

pdf("output/binBw/corr_pseudotime_traj55_UNTREATEDDASATINIB_peakSmooth_fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=5, height=4)
pseudotime_traj55_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_DEG == 0 ) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()



### StartEnd

pdf("output/binBw/corr_pseudotime_traj55_UNTREATEDDASATINIB_peakSmooth_fdrStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=5, height=4)
pseudotime_traj55_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_StartEnd <0.05 ) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()






# DEG and smoothpeak >0

pdf("output/binBw/corr_pseudotime_traj55_UNTREATEDDASATINIB_peakSmooth_peakSmoothOver0fdrDEG05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=5, height=4)
pseudotime_traj55_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_DEG <0.05 ,
            smooth_peak_pseudotime > 0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()



### StartEnd and smoothpeak >0

pdf("output/binBw/corr_pseudotime_traj55_UNTREATEDDASATINIB_peakSmooth_peakSmoothOver0fdrStartEnd0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=5, height=4)
pseudotime_traj55_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter( fdr_StartEnd == 0 ,
            smooth_peak_pseudotime > 0) %>%
    ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max, color = condition)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", label.x.npc = c(0.1, 0.1), label.y.npc = c(0.9, 0.8), aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"), color = condition)) + 
    scale_color_manual(values = c("UNTREATED72hrs" = "blue", "DASATINIB72hrs" = "red")) +
    theme_bw()
dev.off()






# bin
## DEG
pdf("output/binBw/histbin3_pseudotime_traj55_UNTREATEDDASATINIB_peakSmooth_fdrDEG0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=6, height=4)
pseudotime_traj55_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG == 0) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue","DASATINIB72hrs" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



## StartEnd
pdf("output/binBw/histbin3_pseudotime_traj55_UNTREATEDDASATINIB_peakSmooth_fdrStartEnd05median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=6, height=4)
pseudotime_traj55_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(condition, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("UNTREATED72hrs" = "blue","DASATINIB72hrs" = "red")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()



```


--> no / weak correlation

# 2472hrs human gastruloid integration 3D paper; correlation H3K27me3 with pseudotime 002*/003*

Use H3K27me3 ENCODE files generated at  `## ENCODE H3K27me3 file selection`

## Lineage1 Epi/Blood -  Test Activation point with H3K27me3 (from ENCODE) - correlation - 3D paper

The level of H3K27me3 around TSS has been calculated here at `### Count signal around TSS ENCODE H3K27me3`. 



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("ggpubr")


XXXY HERE 

# import files
pseudotime_traj3_peak <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/traj3_noCondition_humangastruloid2472hrs_ActivationPoint.txt") %>%
    dplyr::rename("geneSymbol" = "gene") 
pseudotime_traj3_DEG = read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_association_traj3_noCondition_humangastruloid2472hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_DEG" = "fdr") %>% 
    dplyr::select(geneSymbol, meanLogFC, fdr_DEG)
pseudotime_traj3_StartEnd <- read_tsv("../../002_scRNAseq/003__YAP1/output/condiments/pseudotime_start_end_association_NULL_traj3_noCondition_humangastruloid2472hrs.txt") %>%
    dplyr::rename("geneSymbol" = "gene",
                  "fdr_StartEnd" = "fdr")  %>% 
    dplyr::select(geneSymbol, fdr_StartEnd, logFClineage1)





pseudotime_traj3_peak_DEG_StartEnd = pseudotime_traj3_peak %>%
    left_join(pseudotime_traj3_DEG) %>%
    left_join(pseudotime_traj3_StartEnd)

# H3K27me3 ENCODE #########################

WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt")

pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_500bpTSS_geneSymbol)  %>%
  filter(!is.na(bc_max))
pseudotime_traj3_peak_WT_H3K27me3_250bpTSS_geneSymbol = pseudotime_traj3_peak_DEG_StartEnd %>%
    left_join(WT_H3K27me3_250bpTSS_geneSymbol) %>%
  filter(!is.na(bc_max))


## signal H3K27me3 vs DEG timecourse

# DEG
pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_fdrDEG05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG <0.05 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_frdStartEnd0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd == 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()




# DEG and smoothpeak >0
pdf("output/binBw/corr_pseudotime_traj3_peakSmoothOver0fdrDEG0__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_DEG == 0,
           smooth_peak_pseudotime > 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()

pdf("output/binBw/corr_pseudotime_traj3_peakSmooth_peakSmoothOver0frdStartEnd05__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=4)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    filter(fdr_StartEnd <0.05 ,
           smooth_peak_pseudotime > 0 ) %>%
ggplot(., aes(x = smooth_peak_pseudotime, y = bc_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = 0, label.y = 750, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()
dev.off()



# bin
## DEG
pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_fdrDEG05median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_DEG <0.05) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()

## StartEnd
pdf("output/binBw/histbin3_pseudotime_traj3_peakSmooth_fdrStartEnd0median___WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs.pdf", width=3, height=2)
pseudotime_traj3_peak_WT_H3K27me3_500bpTSS_geneSymbol %>% 
    add_column(genotype = "WT") %>%
    filter(fdr_StartEnd == 0 ) %>%
    mutate(pseudotime_bin = cut(smooth_peak_pseudotime, breaks = c(seq(0, max(smooth_peak_pseudotime), by = 3), Inf), include.lowest = TRUE, right = FALSE)) %>%
    group_by(genotype, pseudotime_bin) %>%
    summarize(mean_bc_max = median(bc_max), se_bc_max = sd(bc_max) / sqrt(n()), .groups = 'drop') %>%
    ggplot(., aes(x = pseudotime_bin, y = mean_bc_max, fill = genotype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        geom_errorbar(aes(ymin = mean_bc_max - se_bc_max, ymax = mean_bc_max + se_bc_max),
                      position = position_dodge(0.9), width = 0.2) +
        scale_fill_manual(values = c("WT" = "blue")) +
        labs(title = "Barplot of Mean bc_max by Pseudotime Bin with Error Bars",
            x = "Pseudotime Bin",
            y = "Mean bc_max") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
dev.off()
```

# Level of H3K27me3 in cell type marker genes _ 3D paper _ Alternative method

## 2472hrs human gastruloid



```R
# packages
library("tidyverse")
library("ggpubr")

# DATA IMPORT #######################
## import cell type marker genes
markers_UNTREATED <- read.delim("../../002_scRNAseq/003__YAP1/output/seurat/srat_humangastruloid2472hrs_UNTREATED_3Dpaper_all_markers.txt", header = TRUE, row.names = 1) %>% as_tibble() %>%
dplyr::rename("geneSymbol" = "gene")
markers_DASATINIB <- read.delim("../../002_scRNAseq/003__YAP1/output/seurat/srat_humangastruloid2472hrs_DASATINIB_3Dpaper_all_markers.txt", header = TRUE, row.names = 1)%>% as_tibble()%>%
dplyr::rename("geneSymbol" = "gene")
## import H3K27me3 level in hESC
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt")


# Identify highly express genes in each cluster
## Group by cluster and arrange by p_val_adj (ascending) and avg_log2FC (descending)
marker_genes <- markers_DASATINIB %>% #  !!!!!!!!! CHANGE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 100) %>% #  !!!!!!!!! CHANGE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ungroup()
## Create a list to save each cluster's top 100 genes
cluster_top_genes <- list()
## Extract unique clusters
clusters <- unique(marker_genes$cluster)
## Loop through each cluster and save the data in the list
for (cl in clusters) {
  cluster_data <- marker_genes %>% filter(cluster == cl)
  cluster_top_genes[[cl]] <- cluster_data
}

cluster_top_genes[["CPC2"]]



# Plot all cluster for quick analysis
plot_data <- bind_rows(
  cluster_top_genes[["Epiblast"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Epiblast"),
  cluster_top_genes[["PrimitiveStreak"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "PrimitiveStreak"),
  cluster_top_genes[["ProliferatingCardiacMesoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "ProliferatingCardiacMesoderm"),
  cluster_top_genes[["CPC1"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CPC1"),
  cluster_top_genes[["CPC2"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CPC2"),   
  cluster_top_genes[["Cardiomyocyte"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Cardiomyocyte"),
  cluster_top_genes[["CadiacMesoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CadiacMesoderm"),
  cluster_top_genes[["Endoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Endoderm"),
  cluster_top_genes[["Mixed_Epiblast_Ectoderm_PrimitiveStreak"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Mixed_Epiblast_Ectoderm_PrimitiveStreak"),
  cluster_top_genes[["Unknown"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Unknown")
)

# Plot the data
plot_data$cluster <- factor(plot_data$cluster, levels = c("Epiblast", "PrimitiveStreak" ,"ProliferatingCardiacMesoderm", "CPC1", "CPC2", "Cardiomyocyte", "CadiacMesoderm", "Endoderm", "Mixed_Epiblast_Ectoderm_PrimitiveStreak", "Unknown")) # Reorder untreated 1st

comparisons <- list(
  c("Epiblast", "PrimitiveStreak"),
  c("Epiblast", "ProliferatingCardiacMesoderm"),
  c("Epiblast", "CPC1"),
  c("Epiblast", "CPC2"),
  c("Epiblast", "Cardiomyocyte"),
  c("Epiblast", "CadiacMesoderm"),
  c("Epiblast", "Endoderm"),
  c("Epiblast", "Mixed_Epiblast_Ectoderm_PrimitiveStreak"),
  c("Epiblast", "Unknown")
)

# Create the boxplot with statistical tests
pdf("output/binBw/barplot_all_top100_DASATINIB__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs_3Dpaper.pdf", width = 3, height = 3)

ggbarplot(
  plot_data, 
  x = "cluster", 
  y = "bc_max", 
  add = "mean_se"
) +
#  stat_compare_means(
#    comparisons = comparisons,
#    method = "wilcox.test", # Use Wilcoxon test for statistical comparisons
#    label = "p.format"
#  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  ) +
  labs(
    x = "Cluster",
    y = "H3K27me3 Level in hESC (Mean ± SE)"
  )
dev.off()






# Lineage2 _ Epi / PS /ProlCardMeso / CPC1 / CPC2 / Cardiomyocyte #############################
### 500bp TSS signal
plot_data <- bind_rows(
  cluster_top_genes[["Epiblast"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Epiblast"),
  cluster_top_genes[["PrimitiveStreak"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "PrimitiveStreak"),
  cluster_top_genes[["ProliferatingCardiacMesoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "ProliferatingCardiacMesoderm"),
  cluster_top_genes[["CPC1"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CPC1"),
  cluster_top_genes[["CPC2"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CPC2"),   
  cluster_top_genes[["Cardiomyocyte"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Cardiomyocyte")
)

# Plot the data
plot_data$cluster <- factor(plot_data$cluster, levels = c("Epiblast", "PrimitiveStreak" ,"ProliferatingCardiacMesoderm", "CPC1", "CPC2", "Cardiomyocyte")) # Reorder untreated 1st

comparisons <- list(
  c("Epiblast", "PrimitiveStreak"),
  c("Epiblast", "ProliferatingCardiacMesoderm"),
  c("Epiblast", "CPC1"),
  c("Epiblast", "CPC2"),
  c("Epiblast", "Cardiomyocyte")
)

# Create the boxplot with statistical tests
pdf("output/binBw/boxplot_traj2_top100_DASATINIB__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs_3Dpaper.pdf", width = 3, height = 3)
ggboxplot(
  plot_data, 
  x = "cluster", 
  y = "bc_max", 
  outlier.shape = NA # Remove outliers from the boxplot
) +
  stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test", # Use Wilcoxon test; you can change to "t.test" if needed
    label = "p.format",
    label.y = c(50, 75, 100, 125, 150)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  ) +
  theme(legend.position = "none") +
  labs(
    x = "Cluster",
    y = "H3K27me3 level\t 
    in hESC"
  ) +
  ylim (0,200)
dev.off()


pdf("output/binBw/barplot_traj2_top100_DASATINIB__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs_3Dpaper.pdf", width = 3, height = 3)
ggbarplot(
  plot_data, 
  x = "cluster", 
  y = "bc_max", 
  add = "mean_se"
) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test", # Use Wilcoxon test for statistical comparisons
    label = "p.format"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  ) +
  labs(
    x = "Cluster",
    y = "H3K27me3 Level in hESC (Mean ± SE)"
  )
dev.off()

### 250bp TSS signal
plot_data <- bind_rows(
  cluster_top_genes[["Epiblast"]] %>% 
    inner_join(WT_H3K27me3_250bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Epiblast"),
  cluster_top_genes[["PrimitiveStreak"]] %>% 
    inner_join(WT_H3K27me3_250bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "PrimitiveStreak"),
  cluster_top_genes[["ProliferatingCardiacMesoderm"]] %>% 
    inner_join(WT_H3K27me3_250bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "ProliferatingCardiacMesoderm"),
  cluster_top_genes[["CPC1"]] %>% 
    inner_join(WT_H3K27me3_250bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CPC1"),
  cluster_top_genes[["CPC2"]] %>% 
    inner_join(WT_H3K27me3_250bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CPC2"),   
  cluster_top_genes[["Cardiomyocyte"]] %>% 
    inner_join(WT_H3K27me3_250bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Cardiomyocyte")
)

# Plot the data
plot_data$cluster <- factor(plot_data$cluster, levels = c("Epiblast", "PrimitiveStreak" ,"ProliferatingCardiacMesoderm", "CPC1", "CPC2", "Cardiomyocyte")) # Reorder untreated 1st

comparisons <- list(
  c("Epiblast", "PrimitiveStreak"),
  c("Epiblast", "ProliferatingCardiacMesoderm"),
  c("Epiblast", "CPC1"),
  c("Epiblast", "CPC2"),
  c("Epiblast", "Cardiomyocyte")
)

# Create the boxplot with statistical tests
pdf("output/binBw/boxplot_traj2_top100_DASATINIB__WT_H3K27me3_250bpTSS_geneSymbol_bc_max_2472hrs_3Dpaper.pdf", width = 3, height = 3)
ggboxplot(
  plot_data, 
  x = "cluster", 
  y = "bc_max", 
  outlier.shape = NA # Remove outliers from the boxplot
) +
  stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test", # Use Wilcoxon test; you can change to "t.test" if needed
    label = "p.format",
    label.y = c(50, 75, 100, 125, 150)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  ) +
  theme(legend.position = "none") +
  labs(
    x = "Cluster",
    y = "H3K27me3 level\t 
    in hESC"
  ) +
  ylim (0,200)
dev.off()


# Lineage1 _ Epi / PS / Blood ############################# 

plot_data <- bind_rows(
  cluster_top_genes[["Epiblast"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Epiblast"),
  cluster_top_genes[["PrimitiveStreak"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "PrimitiveStreak"),
  cluster_top_genes[["CadiacMesoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CadiacMesoderm")
)

# Plot the data
plot_data$cluster <- factor(plot_data$cluster, levels = c("Epiblast", "PrimitiveStreak" ,"CadiacMesoderm")) # Reorder untreated 1st

comparisons <- list(
  c("Epiblast", "PrimitiveStreak"),
  c("Epiblast", "CadiacMesoderm")
)

# Create the boxplot with statistical tests
pdf("output/binBw/boxplot_traj1_top100_UNTREATED__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_2472hrs_3Dpaper.pdf", width = 3, height = 3)
ggboxplot(
  plot_data, 
  x = "cluster", 
  y = "bc_max", 
  outlier.shape = NA # Remove outliers from the boxplot
) +
  stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test", # Use Wilcoxon test; you can change to "t.test" if needed
    label = "p.format",
    label.y = c(50, 75)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  ) +
  theme(legend.position = "none") +
  labs(
    x = "Cluster",
    y = "H3K27me3 level\t 
    in hESC"
  ) +
  ylim (0,150)
dev.off()





```




## 72hrs human gastruloid _ OG clustering


```R
# packages
library("tidyverse")
library("ggpubr")

# DATA IMPORT #######################
## import cell type marker genes
markers_UNTREATED <- read.delim("../../002_scRNAseq/003__YAP1/output/seurat/srat_humangastruloid72hrs_UNTREATED_dim25kparam15res02_3Dpaper_all_markers.txt", header = TRUE, row.names = 1) %>% as_tibble() %>%
dplyr::rename("geneSymbol" = "gene")
markers_DASATINIB <- read.delim("../../002_scRNAseq/003__YAP1/output/seurat/srat_humangastruloid72hrs_DASATINIB_dim25kparam15res02_3Dpaper_all_markers.txt", header = TRUE, row.names = 1)%>% as_tibble()%>%
dplyr::rename("geneSymbol" = "gene")
## import H3K27me3 level in hESC
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt")


# Identify highly express genes in each cluster
## Group by cluster and arrange by p_val_adj (ascending) and avg_log2FC (descending)
marker_genes <- markers_DASATINIB %>% #  !!!!!!!!! CHANGE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 100) %>% #  !!!!!!!!! CHANGE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ungroup()
## Create a list to save each cluster's top 100 genes
cluster_top_genes <- list()
## Extract unique clusters
clusters <- unique(marker_genes$cluster)
## Loop through each cluster and save the data in the list
for (cl in clusters) {
  cluster_data <- marker_genes %>% filter(cluster == cl)
  cluster_top_genes[[cl]] <- cluster_data
}

cluster_top_genes[["CardiacMesoderm"]]

# Plot all cluster for quick analysis
plot_data <- bind_rows(
  cluster_top_genes[["Ectoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Ectoderm"),
  cluster_top_genes[["CardiacMesoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CardiacMesoderm"),
  cluster_top_genes[["CardiacProgenitors"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CardiacProgenitors"),
  cluster_top_genes[["NascentMesoderm1"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "NascentMesoderm1"),
  cluster_top_genes[["NascentMesoderm2"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "NascentMesoderm2"),   
  cluster_top_genes[["Endoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Endoderm")
)

# Plot the data
plot_data$cluster <- factor(plot_data$cluster, levels = c("Ectoderm", "CardiacMesoderm" ,"CardiacProgenitors", "NascentMesoderm1", "NascentMesoderm2", "Endoderm")) # Reorder untreated 1st

comparisons <- list(
  c("Ectoderm", "CardiacMesoderm"),
  c("Ectoderm", "CardiacProgenitors"),
  c("Ectoderm", "NascentMesoderm1"),
  c("Ectoderm", "NascentMesoderm2"),
  c("Ectoderm", "Endoderm")
)

# Create the boxplot with statistical tests
pdf("output/binBw/barplot_all_top100_DASATINIB__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_72hrsOGclustering_3Dpaper.pdf", width = 3, height = 3)
ggbarplot(
  plot_data, 
  x = "cluster", 
  y = "bc_max", 
  add = "mean_se"
) +
 # stat_compare_means(
  #  comparisons = comparisons,
  #  method = "wilcox.test", # Use Wilcoxon test for statistical comparisons
  #  label = "p.format"
  #) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  ) +
  labs(
    x = "Cluster",
    y = "H3K27me3 Level in hESC (Mean ± SE)"
  )
dev.off()

```

--> *NascentMesoderm12* are the population clearly showing low level of H3K27me3.




## 72hrs human gastruloid _ NEW clustering V1


```R
# packages
library("tidyverse")
library("ggpubr")

# DATA IMPORT #######################
## import cell type marker genes
markers_UNTREATED <- read.delim("../../002_scRNAseq/003__YAP1/output/seurat/srat_humangastruloid72hrs_UNTREATED_dim25kparam15res03_3Dpaper_all_markers.txt", header = TRUE, row.names = 1) %>% as_tibble() %>%
dplyr::rename("geneSymbol" = "gene")
markers_DASATINIB <- read.delim("../../002_scRNAseq/003__YAP1/output/seurat/srat_humangastruloid72hrs_DASATINIB_dim25kparam15res03_3Dpaper_all_markers.txt", header = TRUE, row.names = 1)%>% as_tibble()%>%
dplyr::rename("geneSymbol" = "gene")
## import H3K27me3 level in hESC
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt")


# Identify highly express genes in each cluster
## Group by cluster and arrange by p_val_adj (ascending) and avg_log2FC (descending)
marker_genes <- markers_DASATINIB %>% #  !!!!!!!!! CHANGE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 100) %>% #  !!!!!!!!! CHANGE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ungroup()
## Create a list to save each cluster's top 100 genes
cluster_top_genes <- list()
## Extract unique clusters
clusters <- unique(marker_genes$cluster)
## Loop through each cluster and save the data in the list
for (cl in clusters) {
  cluster_data <- marker_genes %>% filter(cluster == cl)
  cluster_top_genes[[cl]] <- cluster_data
}

cluster_top_genes[["CPC"]]

# Plot all cluster for quick analysis
plot_data <- bind_rows(
  cluster_top_genes[["Epiblast_Ectoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Epiblast_Ectoderm"),
  cluster_top_genes[["ProliferatingCardiacMesoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "ProliferatingCardiacMesoderm"),
  cluster_top_genes[["CardiacMesoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CardiacMesoderm"),
  cluster_top_genes[["CPC"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "CPC"),
  cluster_top_genes[["Cardiomyocyte"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Cardiomyocyte"),   
  cluster_top_genes[["Endoderm"]] %>% 
    inner_join(WT_H3K27me3_500bpTSS_geneSymbol, by = c("geneSymbol" = "geneSymbol")) %>% 
    mutate(cluster = "Endoderm")
)

# Plot the data
plot_data$cluster <- factor(plot_data$cluster, levels = c("Epiblast_Ectoderm", "ProliferatingCardiacMesoderm" ,"CardiacMesoderm", "CPC", "Cardiomyocyte", "Endoderm")) # Reorder untreated 1st

comparisons <- list(
  c("Ectoderm", "CardiacMesoderm"),
  c("Ectoderm", "CardiacProgenitors"),
  c("Ectoderm", "NascentMesoderm1"),
  c("Ectoderm", "NascentMesoderm2"),
  c("Ectoderm", "Endoderm")
)

# Create the boxplot with statistical tests
pdf("output/binBw/barplot_all_top100_DASATINIB__WT_H3K27me3_500bpTSS_geneSymbol_bc_max_72hrsNEWv1clustering_3Dpaper.pdf", width = 3, height = 3)
ggbarplot(
  plot_data, 
  x = "cluster", 
  y = "bc_max", 
  add = "mean_se"
) +
  #stat_compare_means(
  #  comparisons = comparisons,
  #  method = "wilcox.test", # Use Wilcoxon test for statistical comparisons
  #  label = "p.format"
  #) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  ) +
  labs(
    x = "Cluster",
    y = "H3K27me3 Level in hESC (Mean ± SE)"
  )
dev.off()

```







# Overlapping peaks EZH2, QSER1, YAP (`008003*/`), TEAD4 (`008003*/`)


Task from Conchi email 20240917:
- collect nb of peaks of EZH2, YAP, QSER1 and TEAD4 (Use ChIPseeker file so that we have gene information) *--> All optimal qval defined as 1.3*
- Identify peak that overlap (+/-200bp)
- Generate venn diagram of overlapping peaks

## macs2 overlapping peaks EZH2, QSER1, YAP (`008003*/`), TEAD4 (`008003*/`)

```bash
# Collect the nb of peaks qval 1.3
wc -l output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103.txt # 4284
wc -l output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.txt # 14166

wc -l ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103.txt # 1650
wc -l ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103.txt # 3351

# extend the peak of 200bp up and down
tail -n +2 output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103_extend200bp.txt # tail used top skip the header (1st line)
tail -n +2 output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103_extend200bp.txt # tail used top skip the header (1st line)
tail -n +2 ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103_extend200bp.txt # tail used top skip the header (1st line)
tail -n +2 ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103_extend200bp.txt # tail used top skip the header (1st line)



# extend the peak of 500bp up and down
tail -n +2 output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103_extend500bp.txt # tail used top skip the header (1st line)
tail -n +2 output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103_extend500bp.txt # tail used top skip the header (1st line)
tail -n +2 ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103_extend500bp.txt # tail used top skip the header (1st line)
tail -n +2 ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103_extend500bp.txt # tail used top skip the header (1st line)

```

--> Venn diagram of peak made using [Intervene webtool](https://www.bioinformatics.com.cn/plot_basic_genomic_regions_overlap_venn_diagram_026_en) From this [paper](https://www.bioinformatics.com.cn/static/papers/026_fig2A.pdf).
  --> See `002*/003*/gastrulation paper/GastrulationPaper_peakOverlap_V1.ppt` for detail  
  --> See `002*/003*/gastrulation paper/output/peakOverlap` for output files  

Let's use **R to add ChIPseeker output to VennDiagram output** (ie. to identify which gene/peak assocation); work in `output/peakOverlap`







```R
library("tidyverse")


#########################################################
# no Extension ##########################################
#########################################################
## import ChIPseeker output
ChIPseeker_QSER1 = read_tsv("output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "QSER1")
ChIPseeker_EZH2 = read_tsv("output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "EZH2")
ChIPseeker_YAP1 = read_tsv("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "YAP1")
ChIPseeker_TEAD4 = read_tsv("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "TEAD4")

ChIPseeker_noExtension = ChIPseeker_QSER1 %>%
  bind_rows(ChIPseeker_EZH2) %>%
  bind_rows(ChIPseeker_YAP1) %>%
  bind_rows(ChIPseeker_TEAD4) 

## import peak venn diagram output
X0001_YAP1 = read_tsv("output/peakOverlap/noExtension/sets/0001_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0001_YAP1")
X0010_TEAD4 = read_tsv("output/peakOverlap/noExtension/sets/0010_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0010_TEAD4")
X0011_TEAD4_YAP1 = read_tsv("output/peakOverlap/noExtension/sets/0011_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0011_TEAD4_YAP1")
X0100_EZH2 = read_tsv("output/peakOverlap/noExtension/sets/0100_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0100_EZH2")
X0110_EZH2_TEAD4 = read_tsv("output/peakOverlap/noExtension/sets/0110_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0110_EZH2_TEAD4")
X1000_QSER1 = read_tsv("output/peakOverlap/noExtension/sets/1000_QSER1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1000_QSER1")
X1001_QSER1_YAP1 = read_tsv("output/peakOverlap/noExtension/sets/1001_QSER1_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1001_QSER1_YAP1")
X1010_QSER1_TEAD4 = read_tsv("output/peakOverlap/noExtension/sets/1010_QSER1_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1010_QSER1_TEAD4")
X1011_QSER1_TEAD4_YAP1 = read_tsv("output/peakOverlap/noExtension/sets/1011_QSER1_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1011_QSER1_TEAD4_YAP1")
X1100_QSER1_EZH2 = read_tsv("output/peakOverlap/noExtension/sets/1100_QSER1_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1100_QSER1_EZH2")
X1110_QSER1_EZH2_TEAD4 = read_tsv("output/peakOverlap/noExtension/sets/1110_QSER1_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1110_QSER1_EZH2_TEAD4")
X1111_QSER1_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/noExtension/sets/1111_QSER1_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1111_QSER1_EZH2_TEAD4_YAP1")

VennDiagram_noExtension = X0001_YAP1 %>%
  bind_rows(X0010_TEAD4) %>%
  bind_rows(X0011_TEAD4_YAP1) %>%
  bind_rows(X0100_EZH2) %>%
  bind_rows(X0110_EZH2_TEAD4) %>%
  bind_rows(X1000_QSER1) %>%
  bind_rows(X1001_QSER1_YAP1) %>%
  bind_rows(X1010_QSER1_TEAD4) %>%
  bind_rows(X1011_QSER1_TEAD4_YAP1) %>%
  bind_rows(X1100_QSER1_EZH2) %>%
  bind_rows(X1110_QSER1_EZH2_TEAD4) %>%
  bind_rows(X1111_QSER1_EZH2_TEAD4_YAP1) 
  
# combine ChIPseeker and VennDiagram outputs
ChIPseeker_VennDiagram_noExtension = VennDiagram_noExtension %>%
  left_join(ChIPseeker_noExtension) %>%
  dplyr::select(VennDiagram, geneSymbol, gene, annotation, distanceToTSS, name, seqnames, start, end, ChIPseeker)

# Save output 
write.table(ChIPseeker_VennDiagram_noExtension, file = "output/peakOverlap/ChIPseeker_VennDiagram_noExtension.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#########################################################
# 200bp Extension ##########################################
#########################################################
## import ChIPseeker output
ChIPseeker_QSER1 = read_tsv("output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103_extend200bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X12", "distanceToTSS" = "X20", "geneSymbol" = "X21", "gene" = "X22") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "QSER1")
ChIPseeker_EZH2 = read_tsv("output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103_extend200bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X12", "distanceToTSS" = "X20", "geneSymbol" = "X21", "gene" = "X22") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "EZH2")
ChIPseeker_YAP1 = read_tsv("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103_extend200bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X12", "distanceToTSS" = "X20", "geneSymbol" = "X21", "gene" = "X22") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "YAP1")
ChIPseeker_TEAD4 = read_tsv("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103_extend200bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X12", "distanceToTSS" = "X20", "geneSymbol" = "X21", "gene" = "X22") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "TEAD4")

ChIPseeker_extend200bp = ChIPseeker_QSER1 %>%
  bind_rows(ChIPseeker_EZH2) %>%
  bind_rows(ChIPseeker_YAP1) %>%
  bind_rows(ChIPseeker_TEAD4) 

## import peak venn diagram output
X0001_YAP1 = read_tsv("output/peakOverlap/200npExtension/sets/0001_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0001_YAP1")
X0010_TEAD4 = read_tsv("output/peakOverlap/200npExtension/sets/0010_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0010_TEAD4")
X0011_TEAD4_YAP1 = read_tsv("output/peakOverlap/200npExtension/sets/0011_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0011_TEAD4_YAP1")
X0100_EZH2 = read_tsv("output/peakOverlap/200npExtension/sets/0100_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0100_EZH2")
X0110_EZH2_TEAD4 = read_tsv("output/peakOverlap/200npExtension/sets/0110_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0110_EZH2_TEAD4")
X1000_QSER1 = read_tsv("output/peakOverlap/200npExtension/sets/1000_QSER1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1000_QSER1")
X1001_QSER1_YAP1 = read_tsv("output/peakOverlap/200npExtension/sets/1001_QSER1_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1001_QSER1_YAP1")
X1010_QSER1_TEAD4 = read_tsv("output/peakOverlap/200npExtension/sets/1010_QSER1_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1010_QSER1_TEAD4")
X1011_QSER1_TEAD4_YAP1 = read_tsv("output/peakOverlap/200npExtension/sets/1011_QSER1_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1011_QSER1_TEAD4_YAP1")
X1100_QSER1_EZH2 = read_tsv("output/peakOverlap/200npExtension/sets/1100_QSER1_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1100_QSER1_EZH2")
X1110_QSER1_EZH2_TEAD4 = read_tsv("output/peakOverlap/200npExtension/sets/1110_QSER1_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1110_QSER1_EZH2_TEAD4")
X1111_QSER1_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/200npExtension/sets/1111_QSER1_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1111_QSER1_EZH2_TEAD4_YAP1")

VennDiagram_200npExtension = X0001_YAP1 %>%
  bind_rows(X0010_TEAD4) %>%
  bind_rows(X0011_TEAD4_YAP1) %>%
  bind_rows(X0100_EZH2) %>%
  bind_rows(X0110_EZH2_TEAD4) %>%
  bind_rows(X1000_QSER1) %>%
  bind_rows(X1001_QSER1_YAP1) %>%
  bind_rows(X1010_QSER1_TEAD4) %>%
  bind_rows(X1011_QSER1_TEAD4_YAP1) %>%
  bind_rows(X1100_QSER1_EZH2) %>%
  bind_rows(X1110_QSER1_EZH2_TEAD4) %>%
  bind_rows(X1111_QSER1_EZH2_TEAD4_YAP1) 
  
# combine ChIPseeker and VennDiagram outputs
ChIPseeker_VennDiagram_extend200bp = VennDiagram_200npExtension %>%
  left_join(ChIPseeker_extend200bp) %>%
  dplyr::select(VennDiagram, geneSymbol, gene, annotation, distanceToTSS, name, seqnames, start, end, ChIPseeker)

# Save output 
write.table(ChIPseeker_VennDiagram_extend200bp, file = "output/peakOverlap/ChIPseeker_VennDiagram_extend200bp.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)








#########################################################
# 500bp Extension ##########################################
#########################################################
## import ChIPseeker output
ChIPseeker_QSER1 = read_tsv("output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103_extend500bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X12", "distanceToTSS" = "X20", "geneSymbol" = "X21", "gene" = "X22") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "QSER1")
ChIPseeker_EZH2 = read_tsv("output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103_extend500bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X12", "distanceToTSS" = "X20", "geneSymbol" = "X21", "gene" = "X22") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "EZH2")
ChIPseeker_YAP1 = read_tsv("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103_extend500bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X12", "distanceToTSS" = "X20", "geneSymbol" = "X21", "gene" = "X22") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "YAP1")
ChIPseeker_TEAD4 = read_tsv("../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103_extend500bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X12", "distanceToTSS" = "X20", "geneSymbol" = "X21", "gene" = "X22") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "TEAD4")

ChIPseeker_extend500bp = ChIPseeker_QSER1 %>%
  bind_rows(ChIPseeker_EZH2) %>%
  bind_rows(ChIPseeker_YAP1) %>%
  bind_rows(ChIPseeker_TEAD4) 

## import peak venn diagram output
X0001_YAP1 = read_tsv("output/peakOverlap/500bpExtension/sets/0001_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0001_YAP1")
X0010_TEAD4 = read_tsv("output/peakOverlap/500bpExtension/sets/0010_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0010_TEAD4")
X0011_TEAD4_YAP1 = read_tsv("output/peakOverlap/500bpExtension/sets/0011_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0011_TEAD4_YAP1")
X0100_EZH2 = read_tsv("output/peakOverlap/500bpExtension/sets/0100_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0100_EZH2")
X0110_EZH2_TEAD4 = read_tsv("output/peakOverlap/500bpExtension/sets/0110_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0110_EZH2_TEAD4")
X0111_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/500bpExtension/sets/0111_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0111_EZH2_TEAD4_YAP1")
X1000_QSER1 = read_tsv("output/peakOverlap/500bpExtension/sets/1000_QSER1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1000_QSER1")
X1001_QSER1_YAP1 = read_tsv("output/peakOverlap/500bpExtension/sets/1001_QSER1_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1001_QSER1_YAP1")
X1010_QSER1_TEAD4 = read_tsv("output/peakOverlap/500bpExtension/sets/1010_QSER1_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1010_QSER1_TEAD4")
X1011_QSER1_TEAD4_YAP1 = read_tsv("output/peakOverlap/500bpExtension/sets/1011_QSER1_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1011_QSER1_TEAD4_YAP1")
X1100_QSER1_EZH2 = read_tsv("output/peakOverlap/500bpExtension/sets/1100_QSER1_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1100_QSER1_EZH2")
X1101_QSER1_EZH2_YAP1 = read_tsv("output/peakOverlap/500bpExtension/sets/1101_QSER1_EZH2_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1101_QSER1_EZH2_YAP1")
X1110_QSER1_EZH2_TEAD4 = read_tsv("output/peakOverlap/500bpExtension/sets/1110_QSER1_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1110_QSER1_EZH2_TEAD4")
X1111_QSER1_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/500bpExtension/sets/1111_QSER1_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1111_QSER1_EZH2_TEAD4_YAP1")

VennDiagram_500bpExtension = X0001_YAP1 %>%
  bind_rows(X0010_TEAD4) %>%
  bind_rows(X0011_TEAD4_YAP1) %>%
  bind_rows(X0100_EZH2) %>%
  bind_rows(X0110_EZH2_TEAD4) %>%
  bind_rows(X0111_EZH2_TEAD4_YAP1) %>%
  bind_rows(X1000_QSER1) %>%
  bind_rows(X1001_QSER1_YAP1) %>%
  bind_rows(X1010_QSER1_TEAD4) %>%
  bind_rows(X1011_QSER1_TEAD4_YAP1) %>%
  bind_rows(X1100_QSER1_EZH2) %>%
  bind_rows(X1101_QSER1_EZH2_YAP1) %>%
  bind_rows(X1110_QSER1_EZH2_TEAD4) %>%
  bind_rows(X1111_QSER1_EZH2_TEAD4_YAP1)

# combine ChIPseeker and VennDiagram outputs
ChIPseeker_VennDiagram_extend500bp = VennDiagram_500bpExtension %>%
  left_join(ChIPseeker_extend500bp) %>%
  dplyr::select(VennDiagram, geneSymbol, gene, annotation, distanceToTSS, name, seqnames, start, end, ChIPseeker)

# Save output 
write.table(ChIPseeker_VennDiagram_extend500bp, file = "output/peakOverlap/ChIPseeker_VennDiagram_extend500bp.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

```

--> File exported to Google Drive and converted to xcl `gastrulation paper/output/peakOverlap/ChIPseeker_VennDiagram.xlsx`




## homer overlapping peaks EZH2, QSER1, YAP (`008003*/`), TEAD4 (`008003*/`)

```bash
# Collect the nb of peaks homer
wc -l output/homer/hESC_WT_QSER1_outputPeaks.bed # 12462
wc -l output/homer/hESC_WT_EZH2_outputPeaks.bed # 2671
wc -l ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks.bed # 5597
wc -l ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks.bed # 3062



# extend the peak of 200bp up and down
tail -n +2 output/homer/hESC_WT_QSER1_outputPeaks.bed | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed # tail used top skip the header (1st line)
tail -n +2 output/homer/hESC_WT_EZH2_outputPeaks.bed | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  output/homer/hESC_WT_EZH2_outputPeaks_extend200bp.bed # tail used top skip the header (1st line)
tail -n +2 ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks.bed | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks_extend200bp.bed # tail used top skip the header (1st line)
tail -n +2 ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks.bed | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks_extend200bp.bed # tail used top skip the header (1st line)

tail -n +2 output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.txt # tail used top skip the header (1st line)
tail -n +2 output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_extend200bp.txt # tail used top skip the header (1st line)
tail -n +2 output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.txt # tail used top skip the header (1st line)
tail -n +2 output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot_extend200bp.txt # tail used top skip the header (1st line)



# extend the peak of 500bp up and down
tail -n +2 output/homer/hESC_WT_QSER1_outputPeaks.bed | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 >  output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed # tail used top skip the header (1st line)
tail -n +2 output/homer/hESC_WT_EZH2_outputPeaks.bed | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 >  output/homer/hESC_WT_EZH2_outputPeaks_extend500bp.bed # tail used top skip the header (1st line)
tail -n +2 ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks.bed | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 >  ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks_extend500bp.bed # tail used top skip the header (1st line)
tail -n +2 ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks.bed | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 >  ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks_extend500bp.bed # tail used top skip the header (1st line)

tail -n +2 output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 >  output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot_extend500bp.txt # tail used top skip the header (1st line)
tail -n +2 output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 >  output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_extend500bp.txt # tail used top skip the header (1st line)
tail -n +2 output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 >  output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot_extend500bp.txt # tail used top skip the header (1st line)
tail -n +2 output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot.txt | bedtools slop -i - -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 >  output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot_extend500bp.txt # tail used top skip the header (1st line)

```

--> Venn diagram of peak made using [Intervene webtool](https://www.bioinformatics.com.cn/plot_basic_genomic_regions_overlap_venn_diagram_026_en) From this [paper](https://www.bioinformatics.com.cn/static/papers/026_fig2A.pdf).
  --> See `002*/003*/gastrulation paper/GastrulationPaper_peakOverlap_V2.ppt` for detail  
  --> See `002*/003*/gastrulation paper/output/peakOverlap/homer_peakCoordinate` for output files  

Let's use **R to add ChIPseeker output to VennDiagram output** (ie. to identify which gene/peak assocation); work in `output/peakOverlap`





```R
library("tidyverse")


#########################################################
# no Extension ##########################################
#########################################################
## import ChIPseeker output
ChIPseeker_QSER1 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "QSER1")
ChIPseeker_EZH2 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "EZH2")
ChIPseeker_YAP1 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "YAP1")
ChIPseeker_TEAD4 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "TEAD4")

ChIPseeker_noExtension = ChIPseeker_QSER1 %>%
  bind_rows(ChIPseeker_EZH2) %>%
  bind_rows(ChIPseeker_YAP1) %>%
  bind_rows(ChIPseeker_TEAD4) 



## import peak venn diagram output
X0001_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/0001_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0001_YAP1")
X0010_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/0010_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0010_TEAD4")
X0011_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/0011_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0011_TEAD4_YAP1")
X0100_EZH2 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/0100_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0100_EZH2")
X0110_EZH2_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/0110_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0110_EZH2_TEAD4")
X0111_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/0111_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0111_EZH2_TEAD4_YAP1")
X1000_QSER1 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/1000_QSER1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1000_QSER1")
X1001_QSER1_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/1001_QSER1_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1001_QSER1_YAP1")
X1010_QSER1_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/1010_QSER1_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1010_QSER1_TEAD4")
X1011_QSER1_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/1011_QSER1_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1011_QSER1_TEAD4_YAP1")
X1100_QSER1_EZH2 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/1100_QSER1_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1100_QSER1_EZH2")
X1110_QSER1_EZH2_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/noExtension/1110_QSER1_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1110_QSER1_EZH2_TEAD4")
  

VennDiagram_noExtension = X0001_YAP1 %>%
  bind_rows(X0010_TEAD4) %>%
  bind_rows(X0011_TEAD4_YAP1) %>%
  bind_rows(X0100_EZH2) %>%
  bind_rows(X0110_EZH2_TEAD4) %>%
  bind_rows(X0111_EZH2_TEAD4_YAP1) %>%
  bind_rows(X1000_QSER1) %>%
  bind_rows(X1001_QSER1_YAP1) %>%
  bind_rows(X1010_QSER1_TEAD4) %>%
  bind_rows(X1011_QSER1_TEAD4_YAP1) %>%
  bind_rows(X1100_QSER1_EZH2) %>%
  bind_rows(X1110_QSER1_EZH2_TEAD4)
  
# combine ChIPseeker and VennDiagram outputs
ChIPseeker_VennDiagram_noExtension = VennDiagram_noExtension %>%
  left_join(ChIPseeker_noExtension) %>%
  dplyr::select(VennDiagram, geneSymbol, gene, annotation, distanceToTSS, name, seqnames, start, end, ChIPseeker)

# Save output 
write.table(ChIPseeker_VennDiagram_noExtension, file = "output/peakOverlap/homer_peakCoordinate/ChIPseeker_VennDiagram_homer_noExtension.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#########################################################
# 200bp Extension ##########################################
#########################################################
## import ChIPseeker output
ChIPseeker_QSER1 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X9", "distanceToTSS" = "X17", "geneSymbol" = "X18", "gene" = "X16") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "QSER1")
ChIPseeker_EZH2 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_extend200bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X9", "distanceToTSS" = "X17", "geneSymbol" = "X18", "gene" = "X16") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "EZH2")
ChIPseeker_YAP1 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X9", "distanceToTSS" = "X17", "geneSymbol" = "X18", "gene" = "X16") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "YAP1")
ChIPseeker_TEAD4 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot_extend200bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X9", "distanceToTSS" = "X17", "geneSymbol" = "X18", "gene" = "X16") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "TEAD4")

ChIPseeker_200bpExtension = ChIPseeker_QSER1 %>%
  bind_rows(ChIPseeker_EZH2) %>%
  bind_rows(ChIPseeker_YAP1) %>%
  bind_rows(ChIPseeker_TEAD4) 



## import peak venn diagram output
X0001_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/0001_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0001_YAP1")
X0010_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/0010_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0010_TEAD4")
X0011_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/0011_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0011_TEAD4_YAP1")
X0100_EZH2 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/0100_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0100_EZH2")
X0110_EZH2_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/0110_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0110_EZH2_TEAD4")
X0111_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/0111_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0111_EZH2_TEAD4_YAP1")
X1000_QSER1 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/1000_QSER1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1000_QSER1")
X1001_QSER1_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/1001_QSER1_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1001_QSER1_YAP1")
X1010_QSER1_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/1010_QSER1_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1010_QSER1_TEAD4")
X1011_QSER1_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/1011_QSER1_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1011_QSER1_TEAD4_YAP1")
X1100_QSER1_EZH2 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/1100_QSER1_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1100_QSER1_EZH2")
X1110_QSER1_EZH2_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/1110_QSER1_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1110_QSER1_EZH2_TEAD4")
X1111_QSER1_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/200bpExtension/1111_QSER1_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1111_QSER1_EZH2_TEAD4_YAP1")

VennDiagram_200bpExtension = X0001_YAP1 %>%
  bind_rows(X0010_TEAD4) %>%
  bind_rows(X0011_TEAD4_YAP1) %>%
  bind_rows(X0100_EZH2) %>%
  bind_rows(X0110_EZH2_TEAD4) %>%
  bind_rows(X0111_EZH2_TEAD4_YAP1) %>%
  bind_rows(X1000_QSER1) %>%
  bind_rows(X1001_QSER1_YAP1) %>%
  bind_rows(X1010_QSER1_TEAD4) %>%
  bind_rows(X1011_QSER1_TEAD4_YAP1) %>%
  bind_rows(X1100_QSER1_EZH2) %>%
  bind_rows(X1110_QSER1_EZH2_TEAD4) %>%
  bind_rows(X1111_QSER1_EZH2_TEAD4_YAP1)
  
# combine ChIPseeker and VennDiagram outputs
ChIPseeker_VennDiagram_200bpExtension = VennDiagram_200bpExtension %>%
  left_join(ChIPseeker_200bpExtension) %>%
  dplyr::select(VennDiagram, geneSymbol, gene, annotation, distanceToTSS, name, seqnames, start, end, ChIPseeker)

# Save output 
write.table(ChIPseeker_VennDiagram_200bpExtension, file = "output/peakOverlap/homer_peakCoordinate/ChIPseeker_VennDiagram_homer_200bpExtension.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)









#########################################################
# 500bp Extension ##########################################
#########################################################
## import ChIPseeker output
ChIPseeker_QSER1 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot_extend500bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X9", "distanceToTSS" = "X17", "geneSymbol" = "X18", "gene" = "X16") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "QSER1")
ChIPseeker_EZH2 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_extend500bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X9", "distanceToTSS" = "X17", "geneSymbol" = "X18", "gene" = "X16") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "EZH2")
ChIPseeker_YAP1 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot_extend500bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X9", "distanceToTSS" = "X17", "geneSymbol" = "X18", "gene" = "X16") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "YAP1")
ChIPseeker_TEAD4 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_TEAD4_R1_annot_extend500bp.txt", col_names = FALSE) %>% 
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3", "name" = "X6", "annotation" = "X9", "distanceToTSS" = "X17", "geneSymbol" = "X18", "gene" = "X16") %>%
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "TEAD4")

ChIPseeker_500bpExtension = ChIPseeker_QSER1 %>%
  bind_rows(ChIPseeker_EZH2) %>%
  bind_rows(ChIPseeker_YAP1) %>%
  bind_rows(ChIPseeker_TEAD4) 



## import peak venn diagram output
X0001_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/0001_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0001_YAP1")
X0010_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/0010_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0010_TEAD4")
X0011_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/0011_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0011_TEAD4_YAP1")
X0100_EZH2 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/0100_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0100_EZH2")
X0110_EZH2_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/0110_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0110_EZH2_TEAD4")
X0111_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/0111_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "0111_EZH2_TEAD4_YAP1")
X1000_QSER1 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/1000_QSER1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1000_QSER1")
X1001_QSER1_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/1001_QSER1_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1001_QSER1_YAP1")
X1010_QSER1_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/1010_QSER1_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1010_QSER1_TEAD4")
X1011_QSER1_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/1011_QSER1_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1011_QSER1_TEAD4_YAP1")
X1100_QSER1_EZH2 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/1100_QSER1_EZH2.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1100_QSER1_EZH2")
X1101_QSER1_EZH2_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/1101_QSER1_EZH2_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1101_QSER1_EZH2_YAP1")
X1110_QSER1_EZH2_TEAD4 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/1110_QSER1_EZH2_TEAD4.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1110_QSER1_EZH2_TEAD4")
X1111_QSER1_EZH2_TEAD4_YAP1 = read_tsv("output/peakOverlap/homer_peakCoordinate/500bpExtension/1111_QSER1_EZH2_TEAD4_YAP1.bed", col_names = FALSE) %>%
  dplyr::rename( "seqnames"= "X1",  "start" = "X2",  "end" = "X3") %>%
  add_column(VennDiagram = "1111_QSER1_EZH2_TEAD4_YAP1")

VennDiagram_500bpExtension = X0001_YAP1 %>%
  bind_rows(X0010_TEAD4) %>%
  bind_rows(X0011_TEAD4_YAP1) %>%
  bind_rows(X0100_EZH2) %>%
  bind_rows(X0110_EZH2_TEAD4) %>%
  bind_rows(X0111_EZH2_TEAD4_YAP1) %>%
  bind_rows(X1000_QSER1) %>%
  bind_rows(X1001_QSER1_YAP1) %>%
  bind_rows(X1010_QSER1_TEAD4) %>%
  bind_rows(X1011_QSER1_TEAD4_YAP1) %>%
  bind_rows(X1100_QSER1_EZH2) %>%
  bind_rows(X1101_QSER1_EZH2_YAP1) %>%
  bind_rows(X1110_QSER1_EZH2_TEAD4) %>%
  bind_rows(X1111_QSER1_EZH2_TEAD4_YAP1)
  
# combine ChIPseeker and VennDiagram outputs
ChIPseeker_VennDiagram_500bpExtension = VennDiagram_500bpExtension %>%
  left_join(ChIPseeker_500bpExtension) %>%
  dplyr::select(VennDiagram, geneSymbol, gene, annotation, distanceToTSS, name, seqnames, start, end, ChIPseeker)

# Save output 
write.table(ChIPseeker_VennDiagram_500bpExtension, file = "output/peakOverlap/homer_peakCoordinate/ChIPseeker_VennDiagram_homer_500bpExtension.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



```

--> File exported to Google Drive and converted to xcl `gastrulation paper/output/peakOverlap/homer_peakCoordinate/ChIPseeker_VennDiagram_homer.xlsx`




# homer YAP1:QSER1 1192 cobound genes with embryo E7 scRNAseq epiblast

-->  How many of the 1,192 YAP:QSER1 bound genes are transcriptionally regulated by YAP in the epiblast? Use the DEGs of the epiblast in the emryo e7 scRNAseq--How many of them are activated and how many are repressed in the YAP KO vs control?---This is to see if YAP:QSER1 binding correlates with YAP-activated or YAP-repressed  genes.


```R
# packages
library("tidyverse")
library("biomaRt")

# import input files
## import DEG Epiblast embryo E7
Epiblast_DEG = read.delim("../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-cYAPKO_response_E7_19dim_allGenes_V2.txt", header = TRUE, row.names = 1) %>% as_tibble(rownames = "gene")
## import the 1192 QSER1:YAP1 bound genes
QSER1_YAP1_genes= read_tsv("output/homer/QSER1_YAP1_1192genes.txt", col_names = FALSE) %>% 
  dplyr::rename( "gene"= "X1")



## convert  human gene id to mice gene id
### Initialize the biomaRt for human and mouse
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# if previous fail: #######
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror="asia")
mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror="asia")

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
############################

### Human gene list 
human_gene_list <- QSER1_YAP1_genes$gene #
### Corresponding mouse orthologs
human_to_mouse <- getLDS(attributes = c("hgnc_symbol"), 
                         filters = "hgnc_symbol", 
                         values = human_gene_list, 
                         mart = human,
                         attributesL = c("mgi_symbol"), 
                         martL = mouse,
                         uniqueRows = TRUE)


# Merge DEG with mice genes
QSER1_YAP1_genes_miceOrtholog = as_tibble(human_to_mouse) %>%
  dplyr::select("MGI.symbol") %>%
  dplyr::rename("gene" = "MGI.symbol" ) %>%
  left_join(Epiblast_DEG)


QSER1_YAP1_genes_miceOrtholog_signif = QSER1_YAP1_genes_miceOrtholog %>%
  filter(p_val_adj <0.05) # 27

sum(QSER1_YAP1_genes_miceOrtholog_signif$avg_log2FC < 0) # 24
sum(QSER1_YAP1_genes_miceOrtholog_signif$avg_log2FC > 0) # 3

Epiblast_DEG_signif = Epiblast_DEG %>%
    filter(p_val_adj <0.05) # 293


sum(Epiblast_DEG_signif$avg_log2FC < 0) # 142
sum(Epiblast_DEG_signif$avg_log2FC > 0) # 151
```


**Metrics**:
Among the 1,192 human genes bound by QSER1:YAP1:
- 910 have mice orthologs
Among the 32,285 mice genes:
- 293 are DEG in the Epiblast (padj <0.05)
 - 142 are downregulated (less express in YAPKO)
   - 151 are upregulated (more express in YAPKO)
Among the 910 mice orthologs:
- 27 are DEGs in the Epiblast (padj <0.05)
   - 24 are downregulated (less express in YAPKO)
   - 3 are upregulated (more express in YAPKO)





# Overlapping QSER1 peaks with CpG islands - macs2 and homer

--> *QSER1 peaks,  QSER1:YAP peaks, and QSER1:EZH2 peaks*, that lie on CpG islands (regardless of methylation status--just by genomic annotation of the regions)
  --> Collect the nb and reprent as histogram proportion CpG island overlap (x peaks, y prop)


CpG islands found here(https://www.biostars.org/p/236141/); in the **UCSC genome annotation database**; imported to `Master/meta`:

```bash
# download CpG islands
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz \
   | gunzip -c \
   | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 }' \
   | sort-bed - \
   > ../../Master/meta/cpgIslandExt.hg38.bed # 32038

```

```bash
# MACS2
## Transform .txt ChIPseeker file into bed file (remove header)
tail -n +2 output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.txt > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.bed
tail -n +2 output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103.txt > output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103.bed
tail -n +2 ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103.txt > ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103.bed
tail -n +2 ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103.txt > ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103.bed


## Overlap QSER1 with YAP1
bedtools intersect -wa -a output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.bed -b ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_YAP1_annot_qval1.30103.bed > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103-overlapYAP1_annot_qval1.30103.bed # 136
## Overlap QSER1 with EZH2
bedtools intersect -wa -a output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.bed -b output/ChIPseeker/annotation_macs2_hESC_WT_EZH2_annot_qval1.30103.bed > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103-overlapEZH2_annot_qval1.30103.bed # 2861
## Overlap QSER1 with TEAD4
bedtools intersect -wa -a output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.bed -b ../003__ChIPseq_pluripotency/output/ChIPseeker/annotation_macs2_hESC_WT_TEAD4_annot_qval1.30103.bed > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103-overlapTEAD4_annot_qval1.30103.bed # 335

## Overlap with CpG islands
### Overlap QSER1 (14165) with CpG islands
bedtools intersect -wa -a output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103-overlapcpgIslandExt.bed # 8506

### Overlap QSER1:YAP1 (136) with CpG islands
bedtools intersect -wa -a output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103-overlapYAP1_annot_qval1.30103.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103-overlapYAP1_annot_qval1.30103-overlapcpgIslandExt.bed # 31

### Overlap QSER1:EZH2 (2861) with CpG islands
bedtools intersect -wa -a output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103-overlapEZH2_annot_qval1.30103.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103-overlapEZH2_annot_qval1.30103-overlapcpgIslandExt.bed # 2728




# HOMER
## Transform .txt ChIPseeker file into bed file (remove header)
tail -n +2 output/homer/hESC_WT_QSER1_outputPeaks.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader.bed
tail -n +2 output/homer/hESC_WT_EZH2_outputPeaks.bed > output/homer/hESC_WT_EZH2_outputPeaks_noHeader.bed
tail -n +41 ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks.bed > ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks_noHeader.bed
tail -n +41 ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks.bed > ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks_noHeader.bed



## Overlap QSER1 with YAP1
bedtools intersect -wa -a output/homer/hESC_WT_QSER1_outputPeaks_noHeader.bed -b ../003__ChIPseq_pluripotency/output/homer/hESC_WT_YAP1_R1/peaks_noHeader.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapYAP1_peaks_noHeader.bed # 199
## Overlap QSER1 with EZH2
bedtools intersect -wa -a output/homer/hESC_WT_QSER1_outputPeaks_noHeader.bed -b output/homer/hESC_WT_EZH2_outputPeaks_noHeader.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapEZH2_peaks_noHeader.bed # 446
## Overlap QSER1 with TEAD4
bedtools intersect -wa -a output/homer/hESC_WT_QSER1_outputPeaks_noHeader.bed -b ../003__ChIPseq_pluripotency/output/homer/hESC_WT_TEAD4_R1/peaks_noHeader.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapTEAD4_peaks_noHeader.bed # 417

## Overlap with CpG islands
### Overlap QSER1 (12461) with CpG islands
bedtools intersect -wa -a output/homer/hESC_WT_QSER1_outputPeaks_noHeader.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapcpgIslandExt.bed # 6067

### Overlap QSER1:YAP1 (199) with CpG islands
bedtools intersect -wa -a output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapYAP1_peaks_noHeader.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapYAP1_peaks_noHeader-overlapcpgIslandExt.bed # 25

### Overlap QSER1:EZH2 (446) with CpG islands
bedtools intersect -wa -a output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapEZH2_peaks_noHeader.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapEZH2_peaks_noHeader-overlapcpgIslandExt.bed # 393
```

**macs2 peak caller**
- *QSER1*; 14,165 peaks: 8,506 overlapping with CpG islands (60.01%)
- *QSER1:EZH2*; 2,861 peaks: 2,728 overlapping with CpG islands (95.35%)
- *QSER1:YAP1*; 136 peaks: 31 overlapping with CpG islands (22.79%)

**homer peak caller**
- *QSER1*; 12,461 peaks: 6,067 overlapping with CpG islands (48.69%)
- *QSER1:EZH2*; 446 peaks: 393 overlapping with CpG islands (88.12%)
- *QSER1:YAP1*; 199 peaks: 25 overlapping with CpG islands (12.56%)

--> homer and macs2 comparable results (same tendency)


--> Hypothesis is that most QSER1 peaks lie on CpG islands regardless methylations status or association with EZH2, YAP1; let's increase QSER1 qvalue and check; and lets try using QSER1 region +/-200 and 500bp



```bash
# MACS2
## Overlap with CpG islands - increasing qvalue
### Overlap QSER1 qval 2.3 (9571) with CpG islands
bedtools intersect -wa -a output/macs2/broad/broad_blacklist_qval2.30103/hESC_WT_QSER1_pool_peaks.broadPeak -b ../../Master/meta/cpgIslandExt.hg38.bed > output/macs2/broad/broad_blacklist_qval2.30103/hESC_WT_QSER1_pool_peaks.broadPeak-overlapcpgIslandExt.bed # 5847
### Overlap QSER1 qval 3 (7297) with CpG islands
bedtools intersect -wa -a output/macs2/broad/broad_blacklist_qval3/hESC_WT_QSER1_pool_peaks.broadPeak -b ../../Master/meta/cpgIslandExt.hg38.bed > output/macs2/broad/broad_blacklist_qval3/hESC_WT_QSER1_pool_peaks.broadPeak-overlapcpgIslandExt.bed # 4614




## Overlap with CpG islands - using increased extended QSER1 region
### Overlap QSER1 qval 1.3 extended 200bp (14165) with CpG islands
bedtools intersect -wa -a output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103_extend200bp.txt -b ../../Master/meta/cpgIslandExt.hg38.bed > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103_extend200bp-overlapcpgIslandExt.bed # 8789
### Overlap QSER1 qval 1.3 extended 500bp (14165) with CpG islands
bedtools intersect -wa -a output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103_extend500bp.txt -b ../../Master/meta/cpgIslandExt.hg38.bed > output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103_extend500bp-overlapcpgIslandExt.bed # 9130


# HOMER
## Overlap with CpG islands - using increased extended QSER1 region
### Overlap QSER1 homer extended 200bp (12461) with CpG islands
bedtools intersect -wa -a output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/homer/hESC_WT_QSER1_outputPeaks_extend200bp-overlapcpgIslandExt.bed # 6597
### Overlap QSER1 homer extended 500bp (12461) with CpG islands
bedtools intersect -wa -a output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/homer/hESC_WT_QSER1_outputPeaks_extend500bp-overlapcpgIslandExt.bed # 7024


```

**macs2 peak caller**
- *+/-200bp QSER1*; 14,165 peaks: 8,789 overlapping with CpG islands (62.05%)
- *+/-500bp QSER1*; 14,165 peaks: 9,130 overlapping with CpG islands (64.45%)
- *qval 2.3*; QSER1; 9,571 peaks: 5,847 overlapping with CpG islands (61.09%)
- *qval3*; QSER1; 7,297 peaks: 4,614 overlapping with CpG islands (63.23%)

**homer peak caller**
- *+/-200bp QSER1*; 12,461 peaks: 6,597 overlapping with CpG islands (52.94%)
- *+/-500bp QSER1*;  12,461 peaks: 7,024 overlapping with CpG islands (56.37%)





--> QSER1 ~60% overlap with CpG islands seems correct (ie. also happened with changes param.)



**macs2 peak caller**
--> Let's confirm QSER1:YAP1 show very low CpG islands binding (~22.79%) by checking the 463 genes co-bound (not overlaping!) by YAP1 and QSER1.
  --> The genes have been copied to file `nano output/ChIPseeker/QSER1_YAP1_463genes.txt` from the Google drive file `gastrulation paper/output/ChIPseeker008001008003/VennDiagramGenes_QSER1EZH2TEAD4YAP1.xlsx`.


```R
# packages
library("tidyverse")

# import files
## import the 463 QSER1:YAP1 co-bound genes 
QSER1_YAP1 = read_tsv("output/ChIPseeker/QSER1_YAP1_463genes.txt", col_names = FALSE) %>% 
  dplyr::rename( "geneSymbol"= "X1")

## import the QSER1 ChIPseeker output to collect peak coordinates
ChIPseeker_QSER1 = read_tsv("output/ChIPseeker/annotation_macs2_hESC_WT_QSER1_annot_qval1.30103.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "QSER1")

# join files and save
ChIPseeker_QSER1 %>% inner_join(QSER1_YAP1) %>%
write.table(., file = "output/ChIPseeker/QSER1_YAP1_463genes_ChIPseeker.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

```


--> 463 genes bound by both YAP1 and QSER1 (not overlapping)
  --> 955 peaks from QSER1 (>463 as multiple peaks bound the same genes); proportion of these peaks overlapping with CpG islands




```bash
## Transform .txt ChIPseeker file into bed file (remove header)
tail -n +2 output/ChIPseeker/QSER1_YAP1_463genes_ChIPseeker.txt > output/ChIPseeker/QSER1_YAP1_463genes_ChIPseeker.bed


### Overlap (955) with CpG islands
bedtools intersect -wa -a output/ChIPseeker/QSER1_YAP1_463genes_ChIPseeker.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/ChIPseeker/QSER1_YAP1_463genes_ChIPseeker-overlapcpgIslandExt.bed # 501
```

- *QSER1 peaks assign to a gene also bound with YAP1 (peaks not overlapping)*; 955 peaks: 501 overlapping with CpG islands (52.46%)







**homer peak caller**
--> Let's confirm QSER1:YAP1 show very low CpG islands binding (~12.56%) by checking the 1192 genes co-bound (not overlaping!) by YAP1 and QSER1.
  --> The genes have been copied to file `nano output/homer/QSER1_YAP1_1192genes.txt` from the Google drive file `gastrulation paper/output/ChIPseeker008001008003/VennDiagramGenes_QSER1EZH2TEAD4YAP1_all_homer.xlsx`.


```R
# packages
library("tidyverse")

# import files
## import the 1192 QSER1:YAP1 co-bound genes 
QSER1_YAP1 = read_tsv("output/homer/QSER1_YAP1_1192genes.txt", col_names = FALSE) %>% 
  dplyr::rename( "geneSymbol"= "X1")

## import the QSER1 ChIPseeker output to collect peak coordinates
ChIPseeker_QSER1 = read_tsv("output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot.txt") %>% 
  dplyr::select(seqnames, start, end, name, annotation, distanceToTSS, geneSymbol, gene) %>%
  add_column(ChIPseeker = "QSER1")

# join files and save
ChIPseeker_QSER1 %>% inner_join(QSER1_YAP1) %>%
write.table(., file = "output/homer/QSER1_YAP1_1192genes_ChIPseeker.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

```


--> 1192 genes bound by both YAP1 and QSER1 (not overlapping)
  --> 2,073 peaks from QSER1 (>1192 as multiple peaks bound the same genes); proportion of these peaks overlapping with CpG islands




```bash
## Transform .txt ChIPseeker file into bed file (remove header)
tail -n +2 output/homer/QSER1_YAP1_1192genes_ChIPseeker.txt > output/homer/QSER1_YAP1_1192genes_ChIPseeker.bed


### Overlap (2,073) with CpG islands
bedtools intersect -wa -a output/homer/QSER1_YAP1_1192genes_ChIPseeker.bed -b ../../Master/meta/cpgIslandExt.hg38.bed > output/homer/QSER1_YAP1_1192genes_ChIPseeker-overlapcpgIslandExt.bed # 763
```

- *QSER1 peaks assign to a gene also bound with YAP1 (peaks not overlapping)*; 2,073 peaks: 763 overlapping with CpG islands (36.81%)
# GO humangastruloid2472hrs cluster on H3K27me3-enriched genes

- Collect top 50-500 marker genes (findAllMarker() done in` 002*/003*`), separatly for UNTREATED and DASATINIB
- Run enrichR and output the ENCODE and ChEA Consensus TFs from ChIP-X; SUZ12 CHEA
  - done in the [web tool app](https://maayanlab.cloud/Enrichr/enrich)
- Create a tidy xls with output result, column of interest: `Overlap, Adjusted P-value, treatment, cluster`

--> The 1st version here is with the top 100 genes; however some cluster does not have 100genes... So I re-run in slurm job the `allGenes` version to have 100 genes


```R
# packages
library("tidyverse")

# TOP 100
## Import output enrichR
top100 <- read_tsv("output/enrichR/PrelimVersion_Top100Genes.txt") %>%
  mutate(prop = (count/totalGene)*100)

top100$cluster <- factor(top100$cluster, levels = c("Endoderm", "Cardiac_Progenitors" ,"Neurogenic_Progenitors", "Muscle_Progenitors", "Progenitor_Cell_Undefined", "Progenitor_Cell_Mitotic", "Nascent_Mesoderm","Epiblast_ESC")) # Reorder untreated 1st

## dotplot
pdf("output/enrichR/dotplot_PrelimVersion_Top100Genes.pdf", width=5, height=2)
ggplot(top100, aes(x = -log10(adjPval), y = cluster, size = prop, fill = treatment)) +
  geom_point(shape = 21) +  # Use shape 21 for filled points
  scale_fill_manual(values = c("UNTREATED" = "blue", "DASATINIB" = "red")) +  # Fill color
  labs(x = "-log 10 Adjusted P-value", y = "Cluster", size = "Proportion") +  # Labels
  geom_vline(xintercept = 1.30103, linetype = "dotted", color = "black") +  # Add dotted black vertical line at x = 1.3
  theme_bw()  # Clean theme
dev.off()



# TOP 50
## Import output enrichR
top50 <- read_tsv("output/enrichR/PrelimVersion_Top50Genes.txt") %>%
  mutate(prop = (count/totalGene)*100)

top50$cluster <- factor(top50$cluster, levels = c("Endoderm", "Cardiac_Progenitors" ,"Neurogenic_Progenitors", "Muscle_Progenitors", "Progenitor_Cell_Undefined", "Progenitor_Cell_Mitotic", "Nascent_Mesoderm","Epiblast_ESC")) # Reorder untreated 1st

## dotplot
pdf("output/enrichR/dotplot_PrelimVersion_Top50Genes.pdf", width=5, height=2)
ggplot(top50, aes(x = -log10(adjPval), y = cluster, size = prop, fill = treatment)) +
  geom_point(shape = 21) +  # Use shape 21 for filled points
  scale_fill_manual(values = c("UNTREATED" = "blue", "DASATINIB" = "red")) +  # Fill color
  labs(x = "-log 10 Adjusted P-value", y = "Cluster", size = "Proportion") +  # Labels
  geom_vline(xintercept = 1.30103, linetype = "dotted", color = "black") +  # Add dotted black vertical line at x = 1.3
  theme_bw()  # Clean theme
dev.off()
```

--> Some sample does not include 100 genes, however I am showing proportion and not gene count, so that is not a problem






# GSEA humangastruloid2472hrs cluster on H3K27me3-enriched genes


- Take top 50-500 H3K27me3-enriched genes in hESC = GSEA input genes
- GSEA using method from `002*/003*` = fgsea() in R 
- Compare GSEA outputs cluster per cluster 




## Collect the top H3K27me3-enriched genes

--> Collect them from ENCODE (Mint-ChIP-seq fcovercontrol with 3 Bio Rep = `output/bigwig_ENCODE/ENCFF201SJZ_fcovercontrol123`)

```R
# packages
library("tidyverse")

# import H3K27me3 gene list
WT_H3K27me3_500bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol.txt") %>% filter(geneSymbol != "NA")
WT_H3K27me3_250bpTSS_geneSymbol <- read_tsv("output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol.txt") %>% filter(geneSymbol != "NA")

# Order and export top enriched genes

WT_H3K27me3_500bpTSS_geneSymbol %>%
  arrange(desc(bc_max)) %>%
  slice_head(n=200) %>%   ############ CHANGE HERE !!!!!!!!!!!!!!! #########
  dplyr::select(geneSymbol) %>%
write.table(., file = "output/binBw/WT_H3K27me3_500bp_ENCFF201SJZ_geneSymbol_topEnriched200.txt", sep = "\t", row.names = FALSE, quote = FALSE) ### CHANGE HERE !!!!!!!!!!!!!!! ####


WT_H3K27me3_250bpTSS_geneSymbol %>%
  arrange(desc(bc_max)) %>%
  slice_head(n=1000) %>%   ############ CHANGE HERE !!!!!!!!!!!!!!! #########
  dplyr::select(geneSymbol) %>%
write.table(., file = "output/binBw/WT_H3K27me3_250bp_ENCFF201SJZ_geneSymbol_topEnriched1000.txt", sep = "\t", row.names = FALSE, quote = FALSE) ### CHANGE HERE !!!!!!!!!!!!!!! ####

```
- *NOTE: 15-500 genes recommended for GSEA*

--> GSEA analysis done in `002*/003*` at xxx




# Upload files to GEO - gastrulation paper

Go [here](https://www.ncbi.nlm.nih.gov/geo/info/seq.html); and follow instructions in `Transfer Files`. Connect to my personal space (`uploads/thomasroule@orcid_A787EGG4`) and transfer files.

- Create a clean `GEO_gastrulationPaper_ChIP008001` folder with all cellranger files (barcode, features, matrix), and fq.gz files
  - I copy `input/hESC*QSER1*`
  - I copy `input/hESC*EZH2*`
  - I copy `input/hESC*input*`
  - I copy `output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-s1-rep*` = WT QSER1 bigwig
  - I copy `output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-s2-rep*` = YAPKO QSER1 bigwig
  - I copy `output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1-rep*` = WT EZH2 bigwig
  - I copy `output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s2-rep*` = YAPKO EZH2 bigwig
- Fill in the `seq_template_TR_gastrulationPaper_ChIPseq.xlsx` (`Metada` and `MD5` sheet notably)
- submit files

```bash
# do file integrity check with md5
md5sum * | awk '{print $2 "\t" $1}' > md5sums.txt

module load lftp

# connect to ftp
lftp -u geoftp,inAlwokhodAbnib5 ftp-private.ncbi.nlm.nih.gov # geoftp = username; inAlwokhodAbnib5 = pwd
cd uploads/thomasroule@orcid_A787EGG4

mirror -R GEO_gastrulationPaper_ChIP008001/
```

--> Done succesfully; release data: 2025-10-30

--> Files in `gastrulation paper` folder






# HOMER peak discovery

From Conchi email 1/8/2025(EMBOR and Cell Reports revisions), we need *DNA sequences/motifs* enriched in some of these QSER1 peaks shown in Figure S3D:
 
- QSER1 bound motifs on all ChIPseq peaks (~12K)
- QSER1 bound motifs on QSER1:YAP-cobound genes (993) --> From file Google Drive `002*/003*/gastrulation paper or QSER1 paper/output/ChIPseeker008001008003/VennDiagramGenes_QSER1EZH2TEAD4YAP1_all_homer.xlsx` there are 1192 genes QSER1 YAP1 bound; let's remove from that the 199 QSER1:YAP-perfectly aligned peak genes = 993 (= 1192 - 199). The 199 can be found at `002*/003*/gastrulation paper or QSER1 paper/output/peakOverlap/homer_peakCoordinate/noExtension/QSER1YAP1_199peaks.txt`. To generate bed file:
  - Collect QSER1 peaks within the 1192 QSER1 YAP1 bound genes. And remove the peak with direct overlap with YAP1.
- QSER1 bound motifs on QSER1:YAP-perfectly aligned peaks (199)
- QSER1 bound motifs on QSER1:EZH2 perfectly aligned peaks (445)

--> Follow [HOMER](http://homer.ucsd.edu/homer/motif/); recommend to use `findMotifs.pl, homer2` with FASTA file

```bash
# bed files of peaks
## Generate the QSER1:YAP-cobound genes (993)
nano output/ChIPseeker/QSER1_YAP1_1192coBoundGenes.txt # copy 002*/003*/gastrulation paper or QSER1 paper/output/ChIPseeker008001008003/VennDiagramGenes_QSER1EZH2TEAD4YAP1_all_homer.xlsx

### Extract the QSER1 peak coordinate of the QSER1_YAP1_1192coBoundGenes
output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot.txt
awk -F'\t' 'NR==FNR {genes[$1]; next} $18 in genes' output/ChIPseeker/QSER1_YAP1_1192coBoundGenes.txt output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot.txt > output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot__1192QSER1YAP1coBoundGenes.txt
#### Double check that it contain our 1192 genes
cut -f18 output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot__1192QSER1YAP1coBoundGenes.txt | sort | uniq | wc -l # 1192

### Remove the QSER1_YAP1 peak that directly overlap (199)


## Bed file locations
output/homer/hESC_WT_QSER1_outputPeaks.bed # all QSER1 peaks (~12k)
output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapYAP1_peaks_noHeader.bed  # QSER1:YAP-perfectly aligned peaks (199)
output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapEZH2_peaks_noHeader.bed # QSER1:EZH2 perfectly aligned peaks (445)



# Convert bed to FASTA
conda activate bowtie2 

bedtools getfasta -fi ../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot__1192QSER1YAP1coBoundGenes.txt > output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot__1192QSER1YAP1coBoundGenes.fa
bedtools getfasta -fi ../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed output/homer/hESC_WT_QSER1_outputPeaks.bed > output/homer/hESC_WT_QSER1_outputPeaks.fa
bedtools getfasta -fi ../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapYAP1_peaks_noHeader.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapYAP1_peaks_noHeader.fa
bedtools getfasta -fi ../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapEZH2_peaks_noHeader.bed > output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapEZH2_peaks_noHeader.fa


# Identify motif from BED
conda activate homer 

findMotifsGenome.pl output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot__1192QSER1YAP1coBoundGenes.txt hg38 output/motif -size 320
findMotifsGenome.pl output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapYAP1_peaks_noHeader.bed hg38 output/motif/hESC_WT_QSER1_outputPeaks_noHeader-overlapYAP1_peaks_noHeader -size 320

findMotifsGenome.pl output/homer/hESC_WT_QSER1_outputPeaks.bed hg38 output/motif/hESC_WT_QSER1_outputPeaks -size 320
findMotifsGenome.pl output/homer/hESC_WT_QSER1_outputPeaks_noHeader-overlapEZH2_peaks_noHeader.bed hg38 output/motif/hESC_WT_QSER1_outputPeaks_noHeader-overlapEZH2_peaks_noHeader -size 320


# Identify motif from FASTA

XXX not needed as bed worked XXX

conda activate homer 

findMotifs.pl output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot__1192QSER1YAP1coBoundGenes.fa fasta output/ChIPseeker/ -fasta  -find motif1.motif > outputfile.txt

```

- *NOTE: `awk -F'\t'` to indicate that file is tab separated!*




