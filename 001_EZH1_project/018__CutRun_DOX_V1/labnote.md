# Project

ESC

Three bio rep with DOX inducible system (KO cell line with DOX to create overexpression):
- WT (no dox)
- KO (no dox)
- OEKO (50ng/mL of dox) 

--> CutRun for: H3K27me3, EZH1, EZH2, IGG (WT only)

**Objectives:**
- Check if the DOX system work! 


# Pipeline
- Download data (wget)
- **combine samples (from multiple lanes)**
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


```bash
module load lftp

# Following email instructions
lftp -c 'set sftp:auto-confirm yes;set net:max-retries 20;open sftp://X202SC25079349-Z01-F001:g4kjwu8a@usftp23.novogene.com; mirror --verbose --use-pget-n=8 -c / ./'

# Copy all .fz.gz data from raw_data into input_raw_Novogene/ folder to be renames
find input_raw/01.RawData/ -type f -name "*.fq.gz" -exec cp {} input_raw_Novogene/ \;


```

--> All good, files fq.gz copied to `input/` to be renamed



# Combine files from multiple lanes

Concatenate fastq discuss [here](https://www.biostars.org/p/317385/): `cat string_L001_sampleID_R1.fastq.gz string_L002_sampleID_R1.fastq.gz  > string_sampleID_R1.fastq.gz`.

--> Several of our samples dispatched in two lanes (`L2` and `L3`; **the concatenated are names `L23`; all unique are `L2`**)
----> input sample: `input_raw_Novogene/` output to `input/`

```bash
# Concatenated samples are directly correctly name

sbatch scripts/concatenate.sh # 49713277 ok
```

--> All good, concatenated samples been rename to and output to `input/`





# Rename file

I created a tab separated file with current (`sample_name.txt`) / new file names (keeping the .fq.gz sufix) with `nano rename_map.txt`

--> **Rename all samples except the concatenated ones (already rename previously)**

```bash
cd input_raw_Novogene

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_map.txt
```

--> All good 



# Fastp cleaning

```bash
sbatch scripts/fastp.sh # 49800144 ok
```


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:49800144 scripts/bowtie2.sh # 49800169 ok
```

--> XXX Looks good; overall ~85% uniquely aligned reads


XXX Mapping on E coli  XXXXXXXXXXXXXXXXXXXXX
```bash
conda activate bowtie2
sbatch scripts/bowtie2_MG1655.sh # 29756394 xxx
```
--> Between 1 - 5% uniquely aligned reads (not a lot..; previously `005__CutRun` 10% (in `003__CutRun` was less than 1%) )
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-49800169.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_49800169.txt
```

Add these values to `/home/roulet/001_EZH1_project/018__CutRun_DOX_V1/samples_001018.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> XXX Overall >75% input reads as been uniquely mapped to the genome (90% non uniq) 



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.




```bash
conda activate bowtie2
sbatch --dependency=afterany:49800169 scripts/samtools_unique.sh # 49800224 ok
```


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
Let's do the same for E coli MG1655 spike in samples:
```bash
conda activate bowtie2
sbatch scripts/samtools_MG1655_unique.sh # 29772397 xxx
```
--> More information on this step in the `005__CutRun` labnote
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools
sbatch --dependency=afterany:49800224 scripts/bamtobigwig_unique.sh # 49800244 ok
```

- WT
*Pass*: H3K27me3 (R1 less good), EZH2 (R1 less good), EZH1 (noisy)
*Failed*: none
*NOTE:* H3K27me3 (R1 outlier less signal), EZH2 (R1 outlier less signal), EZH1 (R2 outlier more signal)
- KO
*Pass*: H3K27me3, EZH2, EZH1
*Failed*: none
*NOTE:* H3K27me3 (homogeneous), EZH2 (R1 outlier less signal), EZH1 (R3 outlier more signal)
- OEKO
*Pass*: H3K27me3, EZH2, EZH1
*Failed*: none
*NOTE:* H3K27me3 (homogeneous), EZH2 (R3 outlier less signal), EZH1 (R1 outlier more signal)

--> It works great, but a lot of disparities in between replicate for signal, many have heterogeneous signal: *normalization will be critical here!*







## Pearson correlation heatmap on bigwig signals

### Raw bigwig


```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_ESC.sh # 50086421 ok
sbatch scripts/multiBigwigSummary_WT.sh # 50086442 ok
sbatch scripts/multiBigwigSummary_KO.sh # 50086524 ok
sbatch scripts/multiBigwigSummary_OEKO.sh # 50086534 ok


############################################
# Plot ESC ###########
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_ESC.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_WT_EZH1_R1 ESC_WT_EZH1_R2 ESC_WT_EZH1_R3 ESC_WT_EZH2_R1 ESC_WT_EZH2_R2 ESC_WT_EZH2_R3 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 ESC_KO_EZH1_R1 ESC_KO_EZH1_R2 ESC_KO_EZH1_R3 ESC_KO_EZH2_R1 ESC_KO_EZH2_R2 ESC_KO_EZH2_R3 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_KO_H3K27me3_R3 ESC_OEKO_EZH1_R1 ESC_OEKO_EZH1_R2 ESC_OEKO_EZH1_R3 ESC_OEKO_EZH2_R1 ESC_OEKO_EZH2_R2 ESC_OEKO_EZH2_R3 ESC_OEKO_H3K27me3_R1 ESC_OEKO_H3K27me3_R2 ESC_OEKO_H3K27me3_R3 \
    --colors black black black black black black black black black black black black red red red red red red red red red blue blue blue blue blue blue blue blue blue \
    --markers 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' '<' '<' '<' 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' \
    -o output/bigwig/multiBigwigSummary_ESC_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_ESC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_WT_EZH1_R1 ESC_WT_EZH1_R2 ESC_WT_EZH1_R3 ESC_WT_EZH2_R1 ESC_WT_EZH2_R2 ESC_WT_EZH2_R3 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 ESC_KO_EZH1_R1 ESC_KO_EZH1_R2 ESC_KO_EZH1_R3 ESC_KO_EZH2_R1 ESC_KO_EZH2_R2 ESC_KO_EZH2_R3 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_KO_H3K27me3_R3 ESC_OEKO_EZH1_R1 ESC_OEKO_EZH1_R2 ESC_OEKO_EZH1_R3 ESC_OEKO_EZH2_R1 ESC_OEKO_EZH2_R2 ESC_OEKO_EZH2_R3 ESC_OEKO_H3K27me3_R1 ESC_OEKO_H3K27me3_R2 ESC_OEKO_H3K27me3_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_ESC_heatmap.pdf

#################################





############################################
# Plot WT ###########
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_WT.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_WT_EZH1_R1 ESC_WT_EZH1_R2 ESC_WT_EZH1_R3 ESC_WT_EZH2_R1 ESC_WT_EZH2_R2 ESC_WT_EZH2_R3 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --colors black black black black black black black black black black black black \
    --markers 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' '<' '<' '<' \
    -o output/bigwig/multiBigwigSummary_WT_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_WT.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_WT_EZH1_R1 ESC_WT_EZH1_R2 ESC_WT_EZH1_R3 ESC_WT_EZH2_R1 ESC_WT_EZH2_R2 ESC_WT_EZH2_R3 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_WT_heatmap.pdf

#################################




############################################
# Plot KO ###########
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_KO.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_KO_EZH1_R1 ESC_KO_EZH1_R2 ESC_KO_EZH1_R3 ESC_KO_EZH2_R1 ESC_KO_EZH2_R2 ESC_KO_EZH2_R3 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_KO_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --colors red red red red red red red red red black black black \
    --markers 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' '<' '<' '<' \
    -o output/bigwig/multiBigwigSummary_KO_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_KO.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_KO_EZH1_R1 ESC_KO_EZH1_R2 ESC_KO_EZH1_R3 ESC_KO_EZH2_R1 ESC_KO_EZH2_R2 ESC_KO_EZH2_R3 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_KO_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_KO_heatmap.pdf

#################################


############################################
# Plot OEKO ###########
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_OEKO.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_OEKO_EZH1_R1 ESC_OEKO_EZH1_R2 ESC_OEKO_EZH1_R3 ESC_OEKO_EZH2_R1 ESC_OEKO_EZH2_R2 ESC_OEKO_EZH2_R3 ESC_OEKO_H3K27me3_R1 ESC_OEKO_H3K27me3_R2 ESC_OEKO_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --colors blue blue blue blue blue blue blue blue blue black black black \
    --markers 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' '<' '<' '<' \
    -o output/bigwig/multiBigwigSummary_OEKO_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_OEKO.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_OEKO_EZH1_R1 ESC_OEKO_EZH1_R2 ESC_OEKO_EZH1_R3 ESC_OEKO_EZH2_R1 ESC_OEKO_EZH2_R2 ESC_OEKO_EZH2_R3 ESC_OEKO_H3K27me3_R1 ESC_OEKO_H3K27me3_R2 ESC_OEKO_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_OEKO_heatmap.pdf

#################################

```

--> look good



### Ferguson norm bigwig

--> Raw IGG used here.


```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_ESC_Ferguson.sh # 50198875 ok
sbatch scripts/multiBigwigSummary_WT_Ferguson.sh # 50199239 ok
sbatch scripts/multiBigwigSummary_KO_Ferguson.sh # 50199685 ok
sbatch scripts/multiBigwigSummary_OEKO_Ferguson.sh # 50200174 ok


############################################
# Plot ESC ###########
## PCA
plotPCA -in output/bigwig_Ferguson/multiBigwigSummary_ESC_Ferguson.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_WT_EZH1_R1 ESC_WT_EZH1_R2 ESC_WT_EZH1_R3 ESC_WT_EZH2_R1 ESC_WT_EZH2_R2 ESC_WT_EZH2_R3 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 ESC_KO_EZH1_R1 ESC_KO_EZH1_R2 ESC_KO_EZH1_R3 ESC_KO_EZH2_R1 ESC_KO_EZH2_R2 ESC_KO_EZH2_R3 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_KO_H3K27me3_R3 ESC_OEKO_EZH1_R1 ESC_OEKO_EZH1_R2 ESC_OEKO_EZH1_R3 ESC_OEKO_EZH2_R1 ESC_OEKO_EZH2_R2 ESC_OEKO_EZH2_R3 ESC_OEKO_H3K27me3_R1 ESC_OEKO_H3K27me3_R2 ESC_OEKO_H3K27me3_R3 \
    --colors black black black black black black black black black red red red red red red red red red blue blue blue blue blue blue blue blue blue \
    --markers 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' \
    -o output/bigwig_Ferguson/multiBigwigSummary_ESC_Ferguson_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig_Ferguson/multiBigwigSummary_ESC_Ferguson.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_WT_EZH1_R1 ESC_WT_EZH1_R2 ESC_WT_EZH1_R3 ESC_WT_EZH2_R1 ESC_WT_EZH2_R2 ESC_WT_EZH2_R3 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 ESC_KO_EZH1_R1 ESC_KO_EZH1_R2 ESC_KO_EZH1_R3 ESC_KO_EZH2_R1 ESC_KO_EZH2_R2 ESC_KO_EZH2_R3 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_KO_H3K27me3_R3 ESC_OEKO_EZH1_R1 ESC_OEKO_EZH1_R2 ESC_OEKO_EZH1_R3 ESC_OEKO_EZH2_R1 ESC_OEKO_EZH2_R2 ESC_OEKO_EZH2_R3 ESC_OEKO_H3K27me3_R1 ESC_OEKO_H3K27me3_R2 ESC_OEKO_H3K27me3_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_Ferguson/multiBigwigSummary_ESC_Ferguson_heatmap.pdf

#################################





############################################
# Plot WT ###########
## PCA
plotPCA -in output/bigwig_Ferguson/multiBigwigSummary_WT_Ferguson.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_WT_EZH1_R1 ESC_WT_EZH1_R2 ESC_WT_EZH1_R3 ESC_WT_EZH2_R1 ESC_WT_EZH2_R2 ESC_WT_EZH2_R3 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --colors black black black black black black black black black black black black \
    --markers 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' '<' '<' '<' \
    -o output/bigwig_Ferguson/multiBigwigSummary_WT_Ferguson_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig_Ferguson/multiBigwigSummary_WT_Ferguson.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_WT_EZH1_R1 ESC_WT_EZH1_R2 ESC_WT_EZH1_R3 ESC_WT_EZH2_R1 ESC_WT_EZH2_R2 ESC_WT_EZH2_R3 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_Ferguson/multiBigwigSummary_WT_Ferguson_heatmap.pdf

#################################




############################################
# Plot KO ###########
## PCA
plotPCA -in output/bigwig_Ferguson/multiBigwigSummary_KO_Ferguson.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_KO_EZH1_R1 ESC_KO_EZH1_R2 ESC_KO_EZH1_R3 ESC_KO_EZH2_R1 ESC_KO_EZH2_R2 ESC_KO_EZH2_R3 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_KO_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --colors red red red red red red red red red black black black \
    --markers 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' '<' '<' '<' \
    -o output/bigwig_Ferguson/multiBigwigSummary_KO_Ferguson_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig_Ferguson/multiBigwigSummary_KO_Ferguson.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_KO_EZH1_R1 ESC_KO_EZH1_R2 ESC_KO_EZH1_R3 ESC_KO_EZH2_R1 ESC_KO_EZH2_R2 ESC_KO_EZH2_R3 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_KO_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_Ferguson/multiBigwigSummary_KO_Ferguson_heatmap.pdf

#################################


############################################
# Plot OEKO ###########
## PCA
plotPCA -in output/bigwig_Ferguson/multiBigwigSummary_OEKO_Ferguson.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_OEKO_EZH1_R1 ESC_OEKO_EZH1_R2 ESC_OEKO_EZH1_R3 ESC_OEKO_EZH2_R1 ESC_OEKO_EZH2_R2 ESC_OEKO_EZH2_R3 ESC_OEKO_H3K27me3_R1 ESC_OEKO_H3K27me3_R2 ESC_OEKO_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --colors blue blue blue blue blue blue blue blue blue black black black \
    --markers 's' 's' 's' 'o' 'o' 'o' '>' '>' '>' '<' '<' '<' \
    -o output/bigwig_Ferguson/multiBigwigSummary_OEKO_Ferguson_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig_Ferguson/multiBigwigSummary_OEKO_Ferguson.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels ESC_OEKO_EZH1_R1 ESC_OEKO_EZH1_R2 ESC_OEKO_EZH1_R3 ESC_OEKO_EZH2_R1 ESC_OEKO_EZH2_R2 ESC_OEKO_EZH2_R3 ESC_OEKO_H3K27me3_R1 ESC_OEKO_H3K27me3_R2 ESC_OEKO_H3K27me3_R3 ESC_WT_IGG_R1 ESC_WT_IGG_R2 ESC_WT_IGG_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig_Ferguson/multiBigwigSummary_OEKO_Ferguson_heatmap.pdf

#################################

```








# MACS2 peak calling on bam unique


--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad.sh # 50090279 ok
sbatch scripts/macs2_broad_noIGG.sh # 50100707 xxx

# genotype per genotype
#XXX sbatch scripts/macs2_narrow.sh #  xxx

```


**broad** with / without IGG:
- *WT*:
    - ESC_WT_EZH1_R1; n(peaks)= 0 / 29
    - ESC_WT_EZH1_R2; n(peaks)= 90 / 27
    - ESC_WT_EZH1_R3; n(peaks)= 0 / 25
    - ESC_WT_EZH2_R1; n(peaks)= 264 / 23
    - ESC_WT_EZH2_R2; n(peaks)= 7223 / 4799
    - ESC_WT_EZH2_R3; n(peaks)= 6666 / 2257
    - ESC_WT_H3K27me3_R1; n(peaks)= 9025 / 6256
    - ESC_WT_H3K27me3_R2; n(peaks)= 26657 / 13328
    - ESC_WT_H3K27me3_R3; n(peaks)= 27812 / 13635
    - ESC_WT_IGG_R1; n(peaks)= NA / 29
    - ESC_WT_IGG_R2; n(peaks)= NA / 38
    - ESC_WT_IGG_R3; n(peaks)= NA / 30

- *KO*:
    - ESC_KO_EZH1_R1; n(peaks)= 0 / 28
    - ESC_KO_EZH1_R2; n(peaks)= 0 / 22
    - ESC_KO_EZH1_R3; n(peaks)= 0 / 27
    - ESC_KO_EZH2_R1; n(peaks)= 1371 / 128
    - ESC_KO_EZH2_R2; n(peaks)= 3929 / 1280
    - ESC_KO_EZH2_R3; n(peaks)= 4162 / 882
    - ESC_KO_H3K27me3_R1; n(peaks)= 12245 / 9265
    - ESC_KO_H3K27me3_R2; n(peaks)= 23307 / 14882
    - ESC_KO_H3K27me3_R3; n(peaks)= 27511 / 16760

- *OEKO*:
    - ESC_OEKO_EZH1_R1; n(peaks)= 2291 / 977
    - ESC_OEKO_EZH1_R2; n(peaks)= 608 / 20
    - ESC_OEKO_EZH1_R3; n(peaks)= 1068 / 25
    - ESC_OEKO_EZH2_R1; n(peaks)= 10324 / 8303
    - ESC_OEKO_EZH2_R2; n(peaks)= 8075 / 5677
    - ESC_OEKO_EZH2_R3; n(peaks)= 12764 / 6417
    - ESC_OEKO_H3K27me3_R1; n(peaks)= 15655 / 10875
    - ESC_OEKO_H3K27me3_R2; n(peaks)= 22711 / 13013
    - ESC_OEKO_H3K27me3_R3; n(peaks)= 25528 / 13393


--> Not sure IGG is better... Probably using IGG is better to  call peak here; less peaks recovered without using IGG!






XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX below not mod



--> OEF1aEZH1 in 50-day neurons: too noisy for EZH1, EZH2, and H3K27me3

--> NPC histone marks: OK for H3K4me3, H3K27ac, and H3K27me3; sharp and clear peaks.

--> NPC PRC2 components: too noisy...

*- NOTE: peak calling has been run 2 times adding the missing samples!*



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







# Normalization method from Ferguson et al 


Summary pipeline:
- bowtie2 unique
- convert bam to bigwig unique 1bp resolution (use unique reads)
- convert bigwig to bedGraph
- identify local maxima
- calculate 99th percentile of the signal


```bash
# Convert bigwig to bedgraph
conda activate BedToBigwig
## Unique bigwig (1bp resolution)
sbatch scripts/BedToBigwig_Ferguson_unique.sh # 50101629 xxx


# Remove blacklist regions
## Unique bigwig (1bp resolution)
sbatch --dependency=afterany:50101629 scripts/BedintersectBlacklist_Ferguson_unique.sh # 50101748 xxx


```

Use Python to identify local maxima, quantify the height for the 99th percentile peak

```bash
srun --mem=250g --pty bash -l


# Identify local maxima
## Unique bigwig (1bp resolution)
python scripts/LocalMaxima_Ferguson_unique.py

#  calculate the 99th percentile of the signal heights (score) in the local maxima files.
## Unique bigwig (1bp resolution)
python scripts/Percentile99_Ferguson_unique.py

# normalize AB per AB (using WT sample 1st replicate as reference)
## Unique bigwig (1bp resolution) APPLYING SF TO INITIAL BIGWIG
### 99th percentile
python scripts/norm_H3K27me3_Ferguson_Perc99_unique_initialBigwig.py
python scripts/norm_EZH2_Ferguson_Perc99_unique_initialBigwig.py
python scripts/norm_EZH1_Ferguson_Perc99_unique_initialBigwig.py

```

*H3K27me3*:
- WT_R1: SF= 1.0
- WT_R2: SF= 0.49166666666666664
- WT_R3: SF= 0.47580645161290325
- KO_R1: SF= 0.6082474226804123
- KO_R2: SF= 0.5042735042735043
- KO_R3: SF= 0.466403162055336
- OEKO_R1: SF= 0.6082474226804123
- OEKO_R2: SF= 0.6519337016574586
- OEKO_R3: SF= 0.5700483091787439
*EZH2*:
- WT_R1: SF= 1.0
- WT_R2: SF= 0.21052631578947367
- WT_R3: SF= 0.23529411764705882
- KO_R1: SF= 0.8
- KO_R2: SF= 0.3076923076923077
- KO_R3: SF= 0.36363636363636365
- OEKO_R1: SF= 0.0784313725490196
- OEKO_R2: SF= 0.1038961038961039
- OEKO_R3: SF= 0.16326530612244897
*EZH1*:
- WT_R1: SF= 1.0
- WT_R2: SF= 1.0
- WT_R3: SF= 1.0
- KO_R1: SF= 1.0
- KO_R2: SF= 1.0
- KO_R3: SF= 1.0
- OEKO_R1: SF= 0.21052631578947367
- OEKO_R2: SF= 1.0
- OEKO_R3: SF= 0.8888888888888888


Convert normalized bedGraph back to bigwig

```bash
conda activate BedToBigwig


# Unique bigwig (1bp resolution) APPLYING SF TO INITIAL BIGWIG
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_H3K27me3.sh # 50105248 ok
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_EZH2.sh # 50105274 ok
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_EZH1.sh # 50105291 ok
```


--> Looks great, replicate are homogeneous



# DIFFREPS for differential binding

According to `001*/016*` best paramters were `bin1000space100_gt_pval05_padj001` = *Bin 1000bp space 100bp, G test, pval 0.05 and padj 0.001: done without FC and with FC 1 treshold*
-> But lets re-test several ones to confirm


## WT vs KO - DIFFREPS - initialBigwig


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bedGraph 

output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bedGraph 


# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bed


## RUN NDIFFREPS ##
# 5000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt --window 5000 --step 100 --meth gt --pval 0.05

# 2000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt --window 2000 --step 100 --meth gt --pval 0.05

# 1000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05

# 500bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05

# 250bp every 50bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05
```




### Explore diffreps results in R




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin5000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin2000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 



# Replace Inf by min/max values
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == Inf] <- max(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == -Inf] <- min(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == Inf] <- max(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == -Inf] <- min(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == Inf] <- max(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == -Inf] <- min(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == Inf] <- max(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)
bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == -Inf] <- min(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)

bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == Inf] <- max(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)
bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == -Inf] <- min(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)



# List of dataset names
file_names <- c("bin5000space100_gt_pval05", "bin2000space100_gt_pval05", "bin1000space100_gt_pval05", "bin500space100_gt_pval05", "bin250space50_gt_pval05")

## Function to read and format each file
read_and_process <- function(file) {
  df <- get(file)  # Load dataset from environment
  df$dataset <- file  # Add dataset identifier
  return(df)
}

## Combine all datasets into one
combined_data <- bind_rows(lapply(file_names, read_and_process)) 



combined_data_counts <- combined_data %>% 
  filter(padj < 0.001, abs(log2FC) > 1, Control.avg > 100 | Treatment.avg > 100) %>%        # keep only |log2FC| > 1
  mutate(direction = if_else(log2FC < 0, "Negative", "Positive")) %>%
  group_by(dataset, direction) %>%
  summarise(count = n(), .groups = "drop")



    
  
## plot
pdf("output/diffreps/hist-WTvsKO-log2FC_distribution-padj001_gt_pval05_fc1_avg100_initialBigwig.pdf", width=8, height=2)
combined_data %>% 
  filter(padj<0.001, abs(log2FC) > 1, Control.avg > 100 | Treatment.avg > 100) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 4, face = "bold")) +
  geom_text(data = combined_data_counts, 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()



## Save output
write.table(combined_data, "output/diffreps/combined_data-WTvsKO-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(combined_data %>% 
  filter(
    padj < 0.001, 
    (log2FC > 1 | log2FC < -1), 
    dataset == "bin1000space100_gt_pval05",
    Control.avg > 100 | Treatment.avg > 100   # <- NEW FILTER
  ), "output/diffreps/combined_data-bin1000space100_gt_pval05_padj001_fc1_avg100-WTvsKO-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(combined_data %>% 
  filter(
    padj < 0.001, 
    (log2FC > 1 | log2FC < -1), 
    dataset == "bin2000space100_gt_pval05",
    Control.avg > 100 | Treatment.avg > 100   # <- NEW FILTER
  ), "output/diffreps/combined_data-bin2000space100_gt_pval05_padj001_fc1_avg100-WTvsKO-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

```


--> To obtain changes that looks real on IGV I had to only keep *.avg with at least 100; so that there is peak changes and not changes from noise to noise!
    --> Overall, this `output/diffreps/combined_data-bin1000space100_gt_pval05_padj001_fc1_avg100-WTvsKO-initialBigwig.txt` look like the best parameters





## WT vs OEKO - DIFFREPS - initialBigwig


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bedGraph 

output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bedGraph 


# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bed


## RUN NDIFFREPS ##
# 5000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt --window 5000 --step 100 --meth gt --pval 0.05

# 2000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt --window 2000 --step 100 --meth gt --pval 0.05

# 1000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05

# 500bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05

# 250bp every 50bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05
```




### Explore diffreps results in R



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin5000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin2000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 



# Replace Inf by min/max values
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == Inf] <- max(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin5000space100_gt_pval05$log2FC[bin5000space100_gt_pval05$log2FC == -Inf] <- min(bin5000space100_gt_pval05$log2FC[is.finite(bin5000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == Inf] <- max(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin2000space100_gt_pval05$log2FC[bin2000space100_gt_pval05$log2FC == -Inf] <- min(bin2000space100_gt_pval05$log2FC[is.finite(bin2000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == Inf] <- max(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)
bin1000space100_gt_pval05$log2FC[bin1000space100_gt_pval05$log2FC == -Inf] <- min(bin1000space100_gt_pval05$log2FC[is.finite(bin1000space100_gt_pval05$log2FC)], na.rm = TRUE)

bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == Inf] <- max(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)
bin500space100_gt_pval05$log2FC[bin500space100_gt_pval05$log2FC == -Inf] <- min(bin500space100_gt_pval05$log2FC[is.finite(bin500space100_gt_pval05$log2FC)], na.rm = TRUE)

bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == Inf] <- max(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)
bin250space50_gt_pval05$log2FC[bin250space50_gt_pval05$log2FC == -Inf] <- min(bin250space50_gt_pval05$log2FC[is.finite(bin250space50_gt_pval05$log2FC)], na.rm = TRUE)



# List of dataset names
file_names <- c("bin5000space100_gt_pval05", "bin2000space100_gt_pval05", "bin1000space100_gt_pval05", "bin500space100_gt_pval05", "bin250space50_gt_pval05")

## Function to read and format each file
read_and_process <- function(file) {
  df <- get(file)  # Load dataset from environment
  df$dataset <- file  # Add dataset identifier
  return(df)
}

## Combine all datasets into one
combined_data <- bind_rows(lapply(file_names, read_and_process)) 



combined_data_counts <- combined_data %>% 
  filter(padj < 0.001, abs(log2FC) > 1, Control.avg > 100 | Treatment.avg > 100) %>%        # keep only |log2FC| > 1
  mutate(direction = if_else(log2FC < 0, "Negative", "Positive")) %>%
  group_by(dataset, direction) %>%
  summarise(count = n(), .groups = "drop")



    
  
## plot
pdf("output/diffreps/hist-WTvsOEKO-log2FC_distribution-padj001_gt_pval05_fc1_avg100_initialBigwig.pdf", width=8, height=2)
combined_data %>% 
  filter(padj<0.001, abs(log2FC) > 1, Control.avg > 100 | Treatment.avg > 100) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
ggplot(., aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "black", color = "black", alpha = 0.7) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +  # Facet per dataset
  labs(title = "Log2FC Distribution Across Datasets",
       x = "Log2 Fold Change (log2FC)",
       y = "Frequency") +
  theme_bw() +
  theme(strip.text = element_text(size = 4, face = "bold")) +
  geom_text(data = combined_data_counts, 
            aes(x = ifelse(direction == "Negative", -6, 4),  # Fixed x positions
                y = Inf, 
                label = paste0(count)), 
            vjust = 1.5, 
            hjust = ifelse(combined_data_counts$direction == "Negative", 0, 1), 
            size = 3, fontface = "bold", color = "red")
dev.off()



## Save output
write.table(combined_data, "output/diffreps/combined_data-WTvsOEKO-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(combined_data %>% 
  filter(
    padj < 0.001, 
    (log2FC > 1 | log2FC < -1), 
    dataset == "bin1000space100_gt_pval05",
    Control.avg > 100 | Treatment.avg > 100   # <- NEW FILTER
  ), "output/diffreps/combined_data-bin1000space100_gt_pval05_padj001_fc1_avg100-WTvsOEKO-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(combined_data %>% 
  filter(
    padj < 0.001, 
    (log2FC > 1 | log2FC < -1), 
    dataset == "bin2000space100_gt_pval05",
    Control.avg > 100 | Treatment.avg > 100   # <- NEW FILTER
  ), "output/diffreps/combined_data-bin2000space100_gt_pval05_padj001_fc1_avg100-WTvsOEKO-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

```


--> Look very similar result to what happen in KO; both show mostly decrease of H3K27me3... And on IGV, look KO and OEKO are the same...








# THOR for differential binding
## Run THOR



```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge


# Default THOR TMM normalization - with input 
sbatch scripts/THOR-ESC_WTvsKO_H3K27me3-TMM.sh # 50510013 xxx
sbatch scripts/THOR-ESC_WTvsOEKO_H3K27me3-TMM.sh # 50510014 xxx

# Default THOR TMM normalization - without input 
sbatch scripts/THOR-ESC_WTvsKO_H3K27me3-TMMnoInput.sh # 50510036 xxx
sbatch scripts/THOR-ESC_WTvsOEKO_H3K27me3-TMMnoInput.sh # 50510042 xxx



# THOR genes normalization - HOX genes used - with input 
sbatch scripts/THOR-ESC_WTvsKO_H3K27me3-housekeepHOX.sh # 50510046 xxx
sbatch scripts/THOR-ESC_WTvsOEKO_H3K27me3-housekeepHOX.sh # 50510295 xxx

# THOR genes normalization - HOX genes used - without input 
sbatch scripts/THOR-ESC_WTvsKO_H3K27me3-housekeepHOXnoInput.sh # 50510049 xxx
sbatch scripts/THOR-ESC_WTvsOEKO_H3K27me3-housekeepHOXnoInput.sh # 50510463 xxx




```

**Conclusion for replicate similarity**:
- *Default THOR TMM normalization*: very bad, replicate very different (potential batch effect remaining, different signal noise ratio)
- *Housekeeping gene normalization with HK genes*: very bad, replicate very different (fail likely because HK genes lowly H3K27me3)
- *Housekeeping gene normalization with HOX genes*: very good, replicate cluster very well together
- *Housekeeping gene normalization with HK + HOX genes*: bad, replicate very different (fail likely because HK genes lowly H3K27me3 again, add noises)
- *Housekeeping gene normalization with HOX genes without input*: very good, replicate cluster very well together
- *Housekeeping gene normalization with HK + HOX genes without input*: bad, replicate very different
- *SpikeIn EpiCypher normalization*: very bad, replicate very different (potential over correction, SF are huge values)
- *SpikeIn EpiCypher DiffBind TMM normalization*: perform OK, less homogeneous than HOX normalization...
- *DiffBindTMM*: perform OK bad, less homogeneous than HOX normalization...
- *FergusonUniqueNorm99*: very good, replicate cluster very well together (*FergusonUniqueNorm99_noInput* perform similarly), let's prefer not using IGG, to be more in agreement with Ferguson method that do not use it.


--> *Error* `IndexError: cannot do a non-empty take from an empty axes.` on `scripts/THOR_PSC_WTvsKO_EZH1_TMM.sh`. Probably because WT vs KO, and KO EZH1 as no signal at all... Weird comparison! That is a control...

--> *Housekeeping genes* has been generated in `001*/002*` at `#### THOR with housekeeping genes normalization`. Collected from the [rgt-THOR tutorial](https://reg-gen.readthedocs.io/en/latest/thor/tool_usage.html)
    --> Let's also try **housekeeping gene normalization using the HOX genes**. Generate in `meta/`: Works great!!

























# ChIPseeker peak gene assignment


### On bin1000space100_gt_pval05_padj001_fc1_avg100 - H3K27me3

Let's assign peak to genes on the two best windowns/parameters as in `001*/009*`:
- Bin 1000bp space 100bp, G test, pval 0.05 and padj 0.001, fc tresh 1 and minimum avg 100 in at least one genotype



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


# Import diff peaks
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100 <- read.delim("output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) %>% 
  filter(padj < 0.001, abs(log2FC) > 1, Control.avg > 100 | Treatment.avg > 100)
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain = PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100 %>%
  filter(log2FC>1)
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost = PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100 %>%
  filter(log2FC<(-1))


PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100 <- read.delim("output/diffreps/ESC_WTvsOEKO_H3K27me3_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) %>% 
  filter(padj < 0.001, abs(log2FC) > 1, Control.avg > 100 | Treatment.avg > 100)
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain = PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100 %>%
  filter(log2FC>1)
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost = PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100 %>%
  filter(log2FC<(-1))


### SAVE Gain and Lost peaks
write.table(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100, file="output/diffreps/PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain, file="output/diffreps/PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost, file="output/diffreps/PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost.txt", sep="\t", quote=F, row.names=F) 

write.table(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100, file="output/diffreps/PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain, file="output/diffreps/PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost, file="output/diffreps/PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost.txt", sep="\t", quote=F, row.names=F) 
########

# Tidy peaks 
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_gr = makeGRangesFromDataFrame(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain,keep.extra.columns=TRUE)
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_gr = makeGRangesFromDataFrame(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost,keep.extra.columns=TRUE)
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_gr = makeGRangesFromDataFrame(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain,keep.extra.columns=TRUE)
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_gr = makeGRangesFromDataFrame(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost,keep.extra.columns=TRUE)

gr_list <- list(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain=PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_gr,PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost=PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_gr, PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain=PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_gr,PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost=PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_gr
)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100_initialBigwig.pdf", width = 16, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100_initialBigwig.pdf", width = 16, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot <- as.data.frame(peakAnnoList[["PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain"]]@anno)
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot <- as.data.frame(peakAnnoList[["PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost"]]@anno)
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot <- as.data.frame(peakAnnoList[["PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain"]]@anno)
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot <- as.data.frame(peakAnnoList[["PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost"]]@anno)

## Convert entrez gene IDs to gene symbols
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot, file="output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot, file="output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot, file="output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot, file="output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot.txt", sep="\t", quote=F, row.names=F)  


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5 = tibble(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5 = tibble(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5 = tibble(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5 = tibble(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol = PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol = PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol = PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol = PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```


Signal changes promoter and 5'  Gain / Lost - `bin1000space100_gt_pval05_padj001_fc1_avg100`:
- WT vs KO GAIN: gene (peak)= 72 (218)
- WT vs KO LOST: gene (peak)= 305 (1768)
- WT vs OEKO GAIN: gene (peak)= 69 (177)
- WT vs OEKO LOST: gene (peak)= 301 (1733)




