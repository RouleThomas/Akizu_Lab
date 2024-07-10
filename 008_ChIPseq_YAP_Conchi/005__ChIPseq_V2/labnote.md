# Project

**SE**

Additional replicate for `001__ChIPseq_V1`:

- input WT and YAPKO
- IP for QSER1, TEAD4, YAP1


- *NOTE: Files weirdly small in size... ~50-500Mb. Weird Expected >Go*


**Objectives:**
- Add aditional replicate to `001__ChIPseq_V1` and integrate all into `006__ChIPseq_V1V2`


# Pipeline
- Download data (wget)
- Rename files
- FastQC (fastqc)
- Trimming (fastp)
- Histone content (R)
- Mapping (bowtie2)
- Spike-in scaling (DiffBind)
- Bigwig generation (deepTools)
- peak calling (MACS2)
- peak assignment to gene (ChIPseeker)

--> Detail of the overall pipeline in `Meeting_20230919_draft.xlsx` 



# Download / import data

- *NOTE: data downloaded in local from basespace/dropbox and then imported into the cluster*
--> Data transferred from EZH1-related ChIPseq (`001*/012*`)





# Rename file


I created a tab separated file with current (`sample_name.txt`) / new file names (keeping the .fq.gz sufix) and then:

**make sure to convert the `sample_name.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input_raw

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < sample_name.txt
```

--> All good 

--> Replicate named *R4/R5/R6* to facilitate integration with *R1/R2/R3* from `008*/001*`


# Fastp cleaning

```bash
sbatch scripts/fastp_1.sh # 21520625 ok
```


# FastQC


**raw**
```bash
sbatch scripts/fastqc_1_raw.sh # 21520628 ok
```


**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:21520625 scripts/fastqc_1_fastp.sh # 21520631 ok
```

--> Weird (below conclusion for IP and input):
- Low yield (<1M vs ~10M)
- Sequencing issue or sample problem
- High Ns in reads: possible degraded/low-quality DNA, low input
- Quality drop at read ends: potential degraded DNA
- Per tile sequence quality issues: possible flow cell overloading or dirt
- Consider re-sequencing or consulting Illumina? Or check bioanalyzer


Doc for (Interpretation of per tile sequence quality)[https://sequencing.qcfail.com/articles/position-specific-failures-of-flowcells/]


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)
--> *NOTE: I removed `--no-mixed --dovetail` and `-U` (instead of `-r`) for the fastq path as PE options* 


```bash
conda activate bowtie2

sbatch --dependency=afterany:21520625 scripts/bowtie2_1.sh # 21520634 ok
```


 
--> Looks ok; around 50-70% uniq aligned reads (70-90% total) 


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX BELO NOT XXXXXXXXXXX

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


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.

```bash
conda activate bowtie2

sbatch --dependency=afterany:21520634 scripts/samtools_unique_1.sh # 21520640 xxx
```

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX NOT MOD XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch scripts/samtools_MG1655_unique_1.sh # 9162457
sbatch scripts/samtools_MG1655_unique_2.sh # 9162461
sbatch scripts/samtools_MG1655_unique_3.sh # 9162467
```

--> More information on this step in the `005__CutRun` labnote



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



XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX





```bash
conda activate deeptools

# bigwig with extendReads 50 (Default ~250bp fragment length  150bp reads +100 bp)
sbatch --dependency=afterany:21520640 scripts/bamtobigwig_unique_extendReads100_1.sh # 21520651 xxx

# bigwig with extendReads from CHIPQC
```







XXXXXX BELOW NOT MOD






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


XXX copy from 008001

