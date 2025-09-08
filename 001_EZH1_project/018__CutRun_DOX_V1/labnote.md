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
sbatch scripts/macs2_broad_noIGG.sh # 50100707 ok

sbatch scripts/macs2_broad_pool_1.sh # 50757176 ok
sbatch scripts/macs2_broad_pool_2.sh # 50757334 ok
sbatch scripts/macs2_broad_pool_3.sh # 50757381 ok


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





## MACS2 peak qvalue filtering

For **consensus peak** counting (ie Ferguson / local maxima method); I used the **pool peak** to generate the consensus peak file for the three genotype comparison

```bash
conda activate bowtie2 # for bedtools

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
- WT_H3K27me3: 2.3 or 3 (maybe more 3)
- KO_H3K27me3: 2.3 or 3 (maybe more 3)
- OEKO_H3K27me3: 2.3 or 3 (maybe more 3)

- WT_EZH2: 2.3 or 3 (maybe more 3)
- KO_EZH2: 2.3 or 3 (maybe more 3)
- OEKO_EZH2: 2.3 or 3 (maybe more 3)

- WT_EZH1: 2.3 
- KO_EZH1: 2.3
- OEKO_EZH1: 2.3




## Identify consensus peaks


Identify peak in WT, KO, OEKO, separately using MACS2, then merge overlapping peak = consensus peak. Then calculate signal in these regions


```bash
conda activate BedToBigwig

# concatenate and sort bed files
## Raw - non qvalue filtered ##############
cat output/macs2/broad/ESC_WT_H3K27me3_pool_peaks.broadPeak output/macs2/broad/ESC_KO_H3K27me3_pool_peaks.broadPeak output/macs2/broad/ESC_OEKO_H3K27me3_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak
cat output/macs2/broad/ESC_WT_EZH2_pool_peaks.broadPeak output/macs2/broad/ESC_KO_EZH2_pool_peaks.broadPeak output/macs2/broad/ESC_OEKO_EZH2_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak
cat output/macs2/broad/ESC_WT_EZH1_pool_peaks.broadPeak output/macs2/broad/ESC_KO_EZH1_pool_peaks.broadPeak output/macs2/broad/ESC_OEKO_EZH1_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak

## qvalue 2.3 ##############
### WT KO KOEF
cat output/macs2/broad/broad_blacklist_qval2.30103/ESC_WT_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/ESC_KO_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/ESC_OEKO_H3K27me3_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval2.30103/ESC_WT_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/ESC_KO_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/ESC_OEKO_EZH2_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval2.30103/ESC_WT_EZH1_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/ESC_KO_EZH1_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval2.30103/ESC_OEKO_EZH1_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak
## qvalue 3 ##############
### WT KO KOEF
cat output/macs2/broad/broad_blacklist_qval3/ESC_WT_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/ESC_KO_H3K27me3_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/ESC_OEKO_H3K27me3_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval3/ESC_WT_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/ESC_KO_EZH2_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/ESC_OEKO_EZH2_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak
cat output/macs2/broad/broad_blacklist_qval3/ESC_WT_EZH1_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/ESC_KO_EZH1_pool_peaks.broadPeak output/macs2/broad/broad_blacklist_qval3/ESC_OEKO_EZH1_pool_peaks.broadPeak | sort -k1,1 -k2,2n > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak



# merge = consensus peak identification
## Raw - non qvalue filtered ##############
### no merge extension
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge.bed
bedtools merge -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.merge.bed
### with 100bp peak merging
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge100bp.bed
bedtools merge -d 100 -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.merge100bp.bed
### with 500bp peak merging
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge500bp.bed
bedtools merge -d 500 -i output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.broadPeak > output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.merge500bp.bed

```

--> All good; consensus peak files are: `output/macs2/broad/ESC_WTKOOEKO_[ANTIBODY]_pool_peaks.sorted.merge[SIZE].bed`





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
python scripts/Percentile99_Ferguson_unique_AUTOSOMALCHR.py



# normalize AB per AB (using WT sample 1st replicate as reference)
## Unique bigwig (1bp resolution) APPLYING SF TO INITIAL BIGWIG
### 99th percentile
python scripts/norm_H3K27me3_Ferguson_Perc99_unique_initialBigwig.py
python scripts/norm_EZH2_Ferguson_Perc99_unique_initialBigwig.py
python scripts/norm_EZH1_Ferguson_Perc99_unique_initialBigwig.py

python scripts/norm_H3K27me3_Ferguson_Perc99_unique_initialBigwig_AUTOSOMALCHR.py
python scripts/norm_EZH2_Ferguson_Perc99_unique_initialBigwig_AUTOSOMALCHR.py
python scripts/norm_EZH1_Ferguson_Perc99_unique_initialBigwig_AUTOSOMALCHR.py

```

*H3K27me3*: (all CHR / autosomal chr only)
- WT_R1: SF= 1.0 / 1.0
- WT_R2: SF= 0.49166666666666664 / 0.48770491803278687
- WT_R3: SF= 0.47580645161290325 / 0.47222222222222222
- KO_R1: SF= 0.6082474226804123 / 0.6102564102564103
- KO_R2: SF= 0.5042735042735043 / 0.5063829787234042
- KO_R3: SF= 0.466403162055336 / 0.468503937007874
- OEKO_R1: SF= 0.6082474226804123 / 0.6071428571428571
- OEKO_R2: SF= 0.6519337016574586 / 0.6538461538461539
- OEKO_R3: SF= 0.5700483091787439 / 0.5721153846153846
*EZH2*:
- WT_R1: SF= 1.0 / 1.0
- WT_R2: SF= 0.21052631578947367 / 0.21052631578947367
- WT_R3: SF= 0.23529411764705882 / 0.23529411764705882
- KO_R1: SF= 0.8 / 0.8
- KO_R2: SF= 0.3076923076923077 / 0.3076923076923077
- KO_R3: SF= 0.36363636363636365 / 0.36363636363636365
- OEKO_R1: SF= 0.0784313725490196 / 0.07766990291262135
- OEKO_R2: SF= 0.1038961038961039 / 0.10256410256410256
- OEKO_R3: SF= 0.16326530612244897 / 0.16
*EZH1*:
- WT_R1: SF= 1.0 / 1.0 
- WT_R2: SF= 1.0 / 1.0
- WT_R3: SF= 1.0 / 1.0
- KO_R1: SF= 1.0 / 1.0
- KO_R2: SF= 1.0 / 1.0
- KO_R3: SF= 1.0 / 1.0
- OEKO_R1: SF= 0.21052631578947367 / 0.21052631578947367
- OEKO_R2: SF= 1.0 / 1.0
- OEKO_R3: SF= 0.8888888888888888 / 0.8888888888888888


--> Using autosomal chr only does not change a lot; probably because scaling based on maximum intensity signal, so sharp peak; not the X chromsome inactivation which is basal, spread, not so high
  --> **Let's use the version with all chr**



Convert normalized bedGraph back to bigwig

```bash
conda activate BedToBigwig


# Unique bigwig (1bp resolution) APPLYING SF TO INITIAL BIGWIG
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_H3K27me3.sh # 50105248 ok
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_EZH2.sh # 50105274 ok
sbatch scripts/BedToBigwig_Norm99_Ferguson_unique_initialBigwig_EZH1.sh # 50105291 ok


# Calculate median
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique_initialBigwig-H3K27me3.sh # 50654175 ok
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique_initialBigwig-EZH2.sh # 50654399 ok
sbatch scripts/bigwigmerge_Norm99_Ferguson_unique_initialBigwig-EZH1.sh # 50654529 ok


# Filter out noise from the bigwig (Consider value below 2 as 0)
sbatch scripts/bigwig_tresh2_Norm99_Ferguson_unique_initialBigwig-H3K27me3.sh # 51404091 ok
sbatch scripts/bigwig_tresh2_Norm99_Ferguson_unique_initialBigwig-EZH2.sh # 51404093 ok
sbatch scripts/bigwig_tresh2_Norm99_Ferguson_unique_initialBigwig-EZH1.sh # 51404094 ok



# Filter out noise from the bigwig (Consider value below 1 as 0)
sbatch scripts/bigwig_tresh1_Norm99_Ferguson_unique_initialBigwig-H3K27me3.sh # 51409597 ok
sbatch scripts/bigwig_tresh1_Norm99_Ferguson_unique_initialBigwig-EZH2.sh # 51409598 ok
sbatch scripts/bigwig_tresh1_Norm99_Ferguson_unique_initialBigwig-EZH1.sh # 51409599 ok

```


--> Looks great, replicate are homogeneous



# DIFFREPS for differential binding

According to `001*/016*` best paramters were `bin1000space100_gt_pval05_padj001` = *Bin 1000bp space 100bp, G test, pval 0.05 and padj 0.001: done without FC and with FC 1 treshold*
-> But lets re-test several ones to confirm


## WT vs KO - DIFFREPS - initialBigwig - H3K27me3


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




### Explore diffreps results in R  - H3K27me3




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







## WT vs KO - DIFFREPS - initialBigwig - EZH2


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bedGraph 

output/bigwig_Ferguson/ESC_KO_EZH2_R1_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_KO_EZH2_R2_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_KO_EZH2_R3_unique_norm99_initialBigwig.bedGraph 


# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_KO_EZH2_R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_KO_EZH2_R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_KO_EZH2_R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_KO_EZH2_R2_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_KO_EZH2_R3_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_KO_EZH2_R3_unique_norm99_initialBigwig.bed


## RUN NDIFFREPS ##
# 5000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt --window 5000 --step 100 --meth gt --pval 0.05

# 2000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt --window 2000 --step 100 --meth gt --pval 0.05

# 1000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05

# 500bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05

# 250bp every 50bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_KO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05
```




### Explore diffreps results in R  - EZH2




```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin5000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin2000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
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
  filter(padj < 0.001, abs(log2FC) > 1, Control.avg > 30 | Treatment.avg > 30) %>%        # keep only |log2FC| > 1
  mutate(direction = if_else(log2FC < 0, "Negative", "Positive")) %>%
  group_by(dataset, direction) %>%
  summarise(count = n(), .groups = "drop")



    
  
## plot
pdf("output/diffreps/hist-WTvsKO-EZH2-log2FC_distribution-padj001_gt_pval05_fc1_avg30_initialBigwig.pdf", width=8, height=2)
combined_data %>% 
  filter(padj<0.001, abs(log2FC) > 1, Control.avg > 30 | Treatment.avg > 30) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
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
write.table(combined_data, "output/diffreps/combined_data-WTvsKO-EZH2-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(combined_data %>% 
  filter(
    padj < 0.001, 
    (log2FC > 1 | log2FC < -1), 
    dataset == "bin1000space100_gt_pval05",
    Control.avg > 30 | Treatment.avg > 30   # <- NEW FILTER
  ), "output/diffreps/combined_data-bin1000space100_gt_pval05_padj001_fc1_avg30-WTvsKO-EZH2-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

```




--> To obtain changes that looks real on IGV I had to only keep *.avg with at least 30; so that there is peak changes and not changes from noise to noise!
    --> Overall, this `output/diffreps/combined_data-bin1000space100_gt_pval05_padj001_fc1_avg30-WTvsKO-EZH2-initialBigwig.txt` look like the best parameters







## WT vs OEKO - DIFFREPS - initialBigwig - H3K27me3


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




### Explore diffreps results in R - H3K27me3



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






## WT vs OEKO - DIFFREPS - initialBigwig - EZH2


```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bedGraph 

output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_unique_norm99_initialBigwig.bedGraph 
output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_unique_norm99_initialBigwig.bedGraph 


# Modify our bedGraph into bed (score in the 5th column); add dummy column 4
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed

awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_unique_norm99_initialBigwig.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "Row" NR, $4, "*"}' output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_unique_norm99_initialBigwig.bedGraph > output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_unique_norm99_initialBigwig.bed


## RUN NDIFFREPS ##
# 5000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt --window 5000 --step 100 --meth gt --pval 0.05

# 2000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt --window 2000 --step 100 --meth gt --pval 0.05

# 1000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05

# 500bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05

# 250bp every 50bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_unique_norm99_initialBigwig.bed -co output/bigwig_Ferguson/ESC_WT_EZH2_R1_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R2_unique_norm99_initialBigwig.bed output/bigwig_Ferguson/ESC_WT_EZH2_R3_unique_norm99_initialBigwig.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05
```




### Explore diffreps results in R - EZH2



```bash
conda activate deseq2
```

```R
# packages
library("tidyverse")
library("GenomicRanges")
set.seed(42)

# import files
bin5000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin5000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin2000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin2000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin1000space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin500space100_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin500space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) 
bin250space50_gt_pval05 <- read.delim("output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin250space50_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
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
  filter(padj < 0.001, abs(log2FC) > 1, Control.avg > 30 | Treatment.avg > 30) %>%        # keep only |log2FC| > 1
  mutate(direction = if_else(log2FC < 0, "Negative", "Positive")) %>%
  group_by(dataset, direction) %>%
  summarise(count = n(), .groups = "drop")



    
  
## plot
pdf("output/diffreps/hist-WTvsOEKO-EZH2-log2FC_distribution-padj001_gt_pval05_fc1_avg30_initialBigwig.pdf", width=8, height=2)
combined_data %>% 
  filter(padj<0.001, abs(log2FC) > 1, Control.avg > 30 | Treatment.avg > 30) %>%   ## !!!!!!!!!! CHANGE PVAL HERE !!!!!!!!!!!!!!!!!!!!!!
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
write.table(combined_data, "output/diffreps/combined_data-WTvsOEKO-EZH2-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(combined_data %>% 
  filter(
    padj < 0.001, 
    (log2FC > 1 | log2FC < -1), 
    dataset == "bin1000space100_gt_pval05",
    Control.avg > 30 | Treatment.avg > 30   # <- NEW FILTER
  ), "output/diffreps/combined_data-bin1000space100_gt_pval05_padj001_fc1_avg30-WTvsOEKO-EZH2-initialBigwig.txt", sep = "\t", quote = FALSE, row.names = FALSE)



```


--> Look very similar result to what happen in KO; both show mostly decrease of H3K27me3... And on IGV, look KO and OEKO are the same...




## WT vs KO - DIFFREPS without X chr - initialBigwig - H3K27me3


Let's directly remove all rows from X chr from the bed `output/bigwig_Ferguson/[SAMPLE ID]_unique_norm99_initialBigwig.bed`



```bash
conda activate ChIPseqSpikeInFree

## PREPARE BED FILE FOR QUANTIFICATION ##
output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed
output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed
output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed

output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bed
output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bed
output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bed

# Remove all rows with X chr
awk 'BEGIN{FS=OFS="\t"} $1!="chrX"' \
  output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bed \
  > output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed
awk 'BEGIN{FS=OFS="\t"} $1!="chrX"' \
  output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.bed \
  > output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed
awk 'BEGIN{FS=OFS="\t"} $1!="chrX"' \
  output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.bed \
  > output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed

awk 'BEGIN{FS=OFS="\t"} $1!="chrX"' \
  output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.bed \
  > output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed
awk 'BEGIN{FS=OFS="\t"} $1!="chrX"' \
  output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.bed \
  > output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed
awk 'BEGIN{FS=OFS="\t"} $1!="chrX"' \
  output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.bed \
  > output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed



## RUN NDIFFREPS ##
# 5000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.nochrX.bed-bin5000space100_gt_pval05-diff.nb.txt --window 5000 --step 100 --meth gt --pval 0.05

# 2000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.nochrX.bed-bin2000space100_gt_pval05-diff.nb.txt --window 2000 --step 100 --meth gt --pval 0.05

# 1000bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.nochrX.bed-bin1000space100_gt_pval05-diff.nb.txt --window 1000 --step 100 --meth gt --pval 0.05

# 500bp every 100bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.nochrX.bed-bin500space100_gt_pval05-diff.nb.txt --window 500 --step 100 --meth gt --pval 0.05

# 250bp every 50bp -  G test pval 0.05
diffReps.pl -tr output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed -co output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_unique_norm99_initialBigwig.nochrX.bed output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_unique_norm99_initialBigwig.nochrX.bed --chrlen ../../Master/meta/GRCh38_chrom_sizes_MAIN.tab -re output/diffreps/ESC_WTvsKO_H3K27me3_unique_norm99_initialBigwig.nochrX.bed-bin250space50_gt_pval05-diff.nb.txt --window 250 --step 50 --meth gt --pval 0.05
```

XXXXY HERE


### Explore diffreps results in R  - H3K27me3




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







# THOR for differential binding
## Run THOR



```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge


# Default THOR TMM normalization - with input 
sbatch scripts/THOR-ESC_WTvsKO_H3K27me3-TMM.sh # 50510013 ok
sbatch scripts/THOR-ESC_WTvsOEKO_H3K27me3-TMM.sh # 50510014 ok

# Default THOR TMM normalization - without input 
sbatch scripts/THOR-ESC_WTvsKO_H3K27me3-TMMnoInput.sh # 50510036 ok
sbatch scripts/THOR-ESC_WTvsOEKO_H3K27me3-TMMnoInput.sh # 50510042 ok



# THOR genes normalization - HOX genes used - with input 
sbatch scripts/THOR-ESC_WTvsKO_H3K27me3-housekeepHOX.sh # 50510046 ok
sbatch scripts/THOR-ESC_WTvsOEKO_H3K27me3-housekeepHOX.sh # 50510295 ok

# THOR genes normalization - HOX genes used - without input 
sbatch scripts/THOR-ESC_WTvsKO_H3K27me3-housekeepHOXnoInput.sh # 50510049 ok
sbatch scripts/THOR-ESC_WTvsOEKO_H3K27me3-housekeepHOXnoInput.sh # 50510463 ok




```

**Conclusion for replicate similarity**:
H3K27me3:
--> All method gave very good replicate homogeneity!
EZH2:
XXX
EZH1:
XXX

--> *Housekeeping genes* has been generated in `001*/002*` at `#### THOR with housekeeping genes normalization`. Collected from the [rgt-THOR tutorial](https://reg-gen.readthedocs.io/en/latest/thor/tool_usage.html)
    --> Let's also try **housekeeping gene normalization using the HOX genes**. Generate in `meta/`: Works great!!



### Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks


```R
# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# H3K27me3 WTvsKO TMM ##################
diffpeaks <- read_tsv("output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMM/ESCWTvsKOH3K27me3TMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMM/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMM/log2FC_qval30.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs KO_qval30") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 40) %>%
  write_tsv("output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMM/THOR_qval40.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 30) %>%
  group_by(X6) %>%
  summarise(n = n())
######################################################







# H3K27me3 WTvsKO TMM noInput ##################
diffpeaks <- read_tsv("output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMMnoInput/ESCWTvsKOH3K27me3TMMnoInput-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2", "count_KO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2+count_KO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMMnoInput/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMMnoInput/log2FC_qval50.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs KO_qval50") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMMnoInput/THOR_qval50.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())
######################################################






# H3K27me3 WTvsOEKO TMM ##################
diffpeaks <- read_tsv("output/THOR/THOR_ESC_WTvsOEKO_H3K27me3_TMM/ESCWTvsOEKOH3K27me3TMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_OEKO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2", "count_WT_3"), sep = ":", convert = TRUE) %>%
  separate(count_OEKO, into = c("count_OEKO_1","count_OEKO_2", "count_OEKO_3"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_OEKO_1+count_OEKO_2+count_OEKO_3) / (count_WT_1+count_WT_2+count_WT_3))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_ESC_WTvsOEKO_H3K27me3_TMM/log2FC.pdf", width=5, height=5)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs OEKO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_ESC_WTvsOEKO_H3K27me3_TMM/log2FC_qval30.pdf", width=5, height=5)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs OEKO_qval30") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_ESC_WTvsOEKO_H3K27me3_TMM/THOR_qval50.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())
######################################################



```


**Optimal qvalue**
- THOR_ESC_WTvsKO_H3K27me3_TMM: qval30
- THOR_ESC_WTvsKO_H3K27me3_TMMnoInput: qval50

--> It seems that the localization of the diff. peaks is always messed up, like in region with minimum to very low signal!
  --> So **I cannot trust THOR method as diff. binding always in WEIRD/ARTIFACTUAL regions...**












# Consensus peak DESEQ2 diff. binding

## Count signal in consensus peaks


I followed `001*/016*` for this part: 

Now calculate **signal in consensus peak**; calculate H3K27me3/EZH2 signal in H3K27me3/EZH2 consensus peaks 100bp merge

--> output files in `output/edgeR`



```bash
conda activate deeptools


# sample per sample (replicate per replicate)

########################################################
## merge no extension  ##############################
########################################################

## qvalue 2.3 ##############
#### WT
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx
#### KO
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx
#### OEKO
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx

## qvalue 3 ##############
#### WT
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx
#### KO
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx
#### OEKO
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_OEKO_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_OEKO_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx





### EZH2
## qvalue 2.3 ##############
#### WT
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R3-FergusonUniqueNorm99.sh #  xxx
#### KO
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R3-FergusonUniqueNorm99.sh #  xxx
#### OEKO
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R3-FergusonUniqueNorm99.sh #  xxx

## qvalue 3 ##############
#### WT
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_WT_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_WT_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_WT_EZH2_R3-FergusonUniqueNorm99.sh #  xxx
#### KO
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_KO_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_KO_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_KO_EZH2_R3-FergusonUniqueNorm99.sh #  xxx
#### OEKO
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_OEKO_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_OEKO_EZH2_R3-FergusonUniqueNorm99.sh #  xxx






########################################################
## merge 100bp extension  ##############################
########################################################

## qvalue 2.3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.sh # 51089151 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.sh # 51089153 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.sh # 51089157 ok
#### KO
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.sh # 51089158 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.sh # 51089161 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.sh # 51089163 ok
#### OEKO
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-FergusonUniqueNorm99.sh # 51089164 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R2-FergusonUniqueNorm99.sh # 51089168 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.sh # 51089172 ok

## qvalue 3 ##############
#### WT
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx
#### KO
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx
#### OEKO
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_OEKO_H3K27me3_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_OEKO_H3K27me3_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.sh #  xxx





### EZH2
## qvalue 2.3 ##############
#### WT
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R1-FergusonUniqueNorm99.sh # 51164275 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R2-FergusonUniqueNorm99.sh # 51164331 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R3-FergusonUniqueNorm99.sh # 51164332 ok
#### KO
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R1-FergusonUniqueNorm99.sh # 51164336 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R2-FergusonUniqueNorm99.sh # 51164391 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R3-FergusonUniqueNorm99.sh # 51164394 ok
#### OEKO
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99.sh # 51164396 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R2-FergusonUniqueNorm99.sh # 51164401 ok
sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R3-FergusonUniqueNorm99.sh # 51164466 ok

## qvalue 3 ##############
#### WT
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_WT_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_WT_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_WT_EZH2_R3-FergusonUniqueNorm99.sh #  xxx
#### KO
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_KO_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_KO_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_KO_EZH2_R3-FergusonUniqueNorm99.sh #  xxx
#### OEKO
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_OEKO_EZH2_R2-FergusonUniqueNorm99.sh #  xxx
#sbatch scripts/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval3merge100bp-ESC_OEKO_EZH2_R3-FergusonUniqueNorm99.sh #  xxx
```






## Perform DESEQ2 comparison


### H3K27me3 - merge100bp qval2.3 - R DESEQ2


```bash
conda activate deseq2
```


```R
library("tidyverse")
library("DESeq2")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot <- read.delim("output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)



# import SCORE 
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_WT_H3K27me3_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
  
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_KO_H3K27me3_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, score per row, coordinate and row

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R1")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R2")
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 = SCORE_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 %>%
  left_join(BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3 ) %>%
  left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R3")




# Tidy into a single tibble
SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp = SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R1 %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_WT_H3K27me3_R3) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R1) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_KO_H3K27me3_R3) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R1) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__ESC_OEKO_H3K27me3_R3)




######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsKO = SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score)) %>%
  filter(!startsWith(peakID, "chrX")) # 112,254 to 97200


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsKO %>%
  mutate(replicate = paste0(genotype, "_", replicate)) %>%  # Create unique column names
  select(-genotype) %>%  # Remove genotype column (since it's now part of replicate)
  pivot_wider(names_from = replicate, values_from = median_score, values_fill = 0)  
  


# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKO, -peakID), pull(countData_WTvsKO, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKO_raw <- SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsKO %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKO_raw, -sample), pull(colData_WTvsKO_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 100 # I TESTED 100, 250, 500: VERY COMPARABLE! So pick 100
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")



## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.58 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot)
# Export result 
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()




# Identify gain lost peak and genes
## 1) Keep only significant peaks, then keep only strong effects (|log2FC| >= threshold)
res_sig <- res_tibble %>%
  filter(padj < 0.05) %>%
  mutate(effect_dir = case_when(
    log2FoldChange >=  0.58 ~ "pos",
    log2FoldChange <= -0.58 ~ "neg",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(effect_dir))   # drop weak effects
## 2) Determine gene-level direction (Gain / Lost / Mix) from the remaining peaks
gene_direction <- res_sig %>%
  group_by(geneSymbol) %>%
  summarise(
    n_pos = sum(effect_dir == "pos"),
    n_neg = sum(effect_dir == "neg"),
    direction = case_when(
      n_pos > 0 & n_neg > 0 ~ "Mix",
      n_pos > 0             ~ "Gain",
      n_neg > 0             ~ "Lost",
      TRUE                  ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(direction))
### Here remove the 'Mix' genes
res_kept_peaks <- res_sig %>%
  inner_join(gene_direction %>% select(geneSymbol, direction), by = "geneSymbol")

# HERE WITHOUT MIX GENES:
upregulated_genes <- res_kept_peaks %>%
  dplyr::filter(direction == "Gain") %>%
  dplyr::select(geneSymbol) %>% unique()
downregulated_genes <- res_kept_peaks %>%
  dplyr::filter(direction == "Lost") %>%
  dplyr::select(geneSymbol) %>% unique()
write.table(upregulated_genes, file = "output/edgeR/upregulatedNoMix_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated_genes, file = "output/edgeR/downregulatedNoMix_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



# HERE KEEPING MIX GENES:
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.58)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.58)


# Save as coordinates for bed deeptools - including Mix peaks

write_tsv( upregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.bed", col_names = FALSE)
write_tsv( downregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.bed", col_names = FALSE)






######################################################
### WT vs OEKO ####################################
######################################################

SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsOEKO = SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp %>%
  filter(genotype %in% c("WT", "OEKO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score)) %>%
  filter(!startsWith(peakID, "chrX")) # 112,254 to 97200






# Convert to wide format
countData_WTvsOEKO <- SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsOEKO %>%
  mutate(replicate = paste0(genotype, "_", replicate)) %>%  # Create unique column names
  select(-genotype) %>%  # Remove genotype column (since it's now part of replicate)
  pivot_wider(names_from = replicate, values_from = median_score, values_fill = 0)  
  


# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsOEKO, -peakID), pull(countData_WTvsOEKO, peakID)) 


## Create colData file that describe all our samples
colData_WTvsOEKO_raw <- SCORE_BED_WTKOOEKO_H3K27me3_qval23merge100bp__WTvsOEKO %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsOEKO_raw, -sample), pull(colData_WTvsOEKO_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 100 # I TESTED 100, 250, 500: VERY COMPARABLE! So pick 100
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_OEKO_vs_WT", type="apeglm")


## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.58 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot)
# Export result 
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'OEKO vs WT, ESC, H3K27me3',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()



# Identify gain lost peak and genes
## 1) Keep only significant peaks, then keep only strong effects (|log2FC| >= threshold)
res_sig <- res_tibble %>%
  filter(padj < 0.05) %>%
  mutate(effect_dir = case_when(
    log2FoldChange >=  0.58 ~ "pos",
    log2FoldChange <= -0.58 ~ "neg",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(effect_dir))   # drop weak effects
## 2) Determine gene-level direction (Gain / Lost / Mix) from the remaining peaks
gene_direction <- res_sig %>%
  group_by(geneSymbol) %>%
  summarise(
    n_pos = sum(effect_dir == "pos"),
    n_neg = sum(effect_dir == "neg"),
    direction = case_when(
      n_pos > 0 & n_neg > 0 ~ "Mix",
      n_pos > 0             ~ "Gain",
      n_neg > 0             ~ "Lost",
      TRUE                  ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(direction))
### Here remove the 'Mix' genes
res_kept_peaks <- res_sig %>%
  inner_join(gene_direction %>% select(geneSymbol, direction), by = "geneSymbol")

# HERE WITHOUT MIX GENES:
upregulated_genes <- res_kept_peaks %>%
  dplyr::filter(direction == "Gain") %>%
  dplyr::select(geneSymbol) %>% unique()
downregulated_genes <- res_kept_peaks %>%
  dplyr::filter(direction == "Lost") %>%
  dplyr::select(geneSymbol) %>% unique()
write.table(upregulated_genes, file = "output/edgeR/upregulatedNoMix_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated_genes, file = "output/edgeR/downregulatedNoMix_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



# HERE KEEPING MIX GENES:
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)




res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.58)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.58)


# Save as coordinates for bed deeptools - including Mix peaks

write_tsv( upregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.bed", col_names = FALSE)
write_tsv( downregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.bed", col_names = FALSE)




```









### EZH2 - merge100bp qval2.3 - R DESEQ2


```bash
conda activate deseq2V2
```


```R
library("tidyverse")
library("DESeq2")
library("EnhancedVolcano")


set.seed(42)

# import bed reference to collect gene name
ESC_WTKOOEKO_EZH2_qval23merge100bp_annot <- read.delim("output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval23merge100bp_annot.txt", header=TRUE, sep="\t", skip=0) %>% 
  as_tibble() %>%
  dplyr::rename(chr = seqnames) %>%
  mutate(peakID = paste(chr, start, end, sep = "_")) %>%
  dplyr::select(chr, start, end, annotation, geneSymbol, gene, peakID)



# import SCORE 
SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())

SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R2-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())
SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R3-FergusonUniqueNorm99.txt", header=FALSE, sep="\t", skip=3) %>%
  as_tibble() %>%
  dplyr::rename(score = V1) %>%
  mutate(rowNumber = row_number())




# import BED position from matrix
BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_WT_EZH2_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
  
BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_KO_EZH2_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())

BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R1 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R2 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R2-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())
BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R3 <- read.delim("output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R3-FergusonUniqueNorm99.bed", header=TRUE, sep="\t", skip=0) %>%
  as_tibble() %>%
  dplyr::rename(chr = "X.chrom") %>%
  dplyr::select(chr, start, end) %>%
  mutate(rowNumber = row_number())



# Put together, gene name, score per row, coordinate and row

SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R1 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R1 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R1 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R1")
SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R2 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R2 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R2 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R2")
SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R3 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R3 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R3 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "WT", replicate = "R3")

SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R1 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R1 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R1 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R1")
SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R2 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R2 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R2 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R2")
SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R3 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R3 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R3 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "KO", replicate = "R3")

SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R1 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R1 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R1 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R1")
SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R2 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R2 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R2 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R2")
SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R3 = SCORE_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R3 %>%
  left_join(BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R3 ) %>%
  left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
  dplyr::select(peakID, score) %>%
  group_by(peakID) %>%  # Group by gene
  summarise(median_score = median(score, na.rm = TRUE)) %>%  # Compute median signal per gene
  unique() %>%
  add_column(genotype = "OEKO", replicate = "R3")




# Tidy into a single tibble
SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp = SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R1 %>%
  bind_rows(SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_WT_EZH2_R3) %>%
  bind_rows(SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R1) %>%
  bind_rows(SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_KO_EZH2_R3) %>%
  bind_rows(SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R1) %>%
  bind_rows(SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R2) %>%
  bind_rows(SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__ESC_OEKO_EZH2_R3)




######################################################
### WT vs KO ####################################
######################################################

SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__WTvsKO = SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp %>%
  filter(genotype %in% c("WT", "KO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score)) %>%
  filter(!startsWith(peakID, "chrX")) # 143,124 to 139,122


# Convert to wide format
countData_WTvsKO <- SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__WTvsKO %>%
  mutate(replicate = paste0(genotype, "_", replicate)) %>%  # Create unique column names
  select(-genotype) %>%  # Remove genotype column (since it's now part of replicate)
  pivot_wider(names_from = replicate, values_from = median_score, values_fill = 0)  
  


# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsKO, -peakID), pull(countData_WTvsKO, peakID)) 


## Create colData file that describe all our samples
colData_WTvsKO_raw <- SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__WTvsKO %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsKO_raw, -sample), pull(colData_WTvsKO_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 100 # I TESTED 100, 250, 500: VERY COMPARABLE! So pick 100
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_KO_vs_WT", type="apeglm")


## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.58 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot)
# Export result 
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'KO vs WT, ESC, EZH2',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.58)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.58)


# Save as coordinates for bed deeptools

write_tsv( upregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.bed", col_names = FALSE)
write_tsv( downregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.bed", col_names = FALSE)






######################################################
### WT vs OEKO ####################################
######################################################

SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__WTvsOEKO = SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp %>%
  filter(genotype %in% c("WT", "OEKO"),
         peakID != "NA") %>%
  mutate(median_score = round(median_score)) %>%
  filter(!startsWith(peakID, "chrX")) # 143,124 to 139,122






# Convert to wide format
countData_WTvsOEKO <- SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__WTvsOEKO %>%
  mutate(replicate = paste0(genotype, "_", replicate)) %>%  # Create unique column names
  select(-genotype) %>%  # Remove genotype column (since it's now part of replicate)
  pivot_wider(names_from = replicate, values_from = median_score, values_fill = 0)  
  


# Pre-requisetes for the DESeqDataSet
## Transform merged_data into a matrix
### Function to transform tibble into matrix
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}
### execute function
counts_all_matrix = make_matrix(dplyr::select(countData_WTvsOEKO, -peakID), pull(countData_WTvsOEKO, peakID)) 


## Create colData file that describe all our samples
colData_WTvsOEKO_raw <- SCORE_BED_WTKOOEKO_EZH2_qval23merge100bp__WTvsOEKO %>%
  distinct(replicate, genotype) %>%
  mutate(sample = paste(genotype, replicate, sep = "_"))
  
  
## transform df into matrix
coldata = make_matrix(dplyr::select(colData_WTvsOEKO_raw, -sample), pull(colData_WTvsOEKO_raw, sample))

## Check that row name of both matrix (counts and description) are the same
all(rownames(coldata) %in% colnames(counts_all_matrix)) # output TRUE is correct

## Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_all_matrix,
                              colData = coldata,
                              design= ~ genotype)

# DEGs
## Filter out gene with less than 5 reads
keep <- rowSums(counts(dds)) >= 100 # I TESTED 100, 250, 500: VERY COMPARABLE! So pick 100
dds <- dds[keep,]

## Specify the control sample
dds$genotype <- relevel(dds$genotype, ref = "WT")

## Differential expression analyses
dds <- DESeq(dds)
# res <- results(dds) # This is the classic version, but shrunk log FC is preferable
resultsNames(dds) # Here print value into coef below

res <- lfcShrink(dds, coef="genotype_OEKO_vs_WT", type="apeglm")


## Plot-volcano
# FILTER ON QVALUE 0.05 GOOD !!!! ###############################################
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 5e-2, 'Sky Blue',
    ifelse(res$log2FoldChange > 0.58 & res$padj < 5e-2, 'Orange',
      'grey'))



keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'Orange'] <- 'Up-regulated (q-val < 0.05; log2FC > 0.58)'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'
names(keyvals)[keyvals == 'Sky Blue'] <- 'Down-regulated (q-val < 0.05; log2FC < -0.58)'


res_tibble <- as_tibble(res, rownames = "peakID") %>% left_join(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot)
# Export result 
write.table(res_tibble, file="output/edgeR/DESEQ2-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.txt", sep="\t", row.names=FALSE, quote=FALSE)

pdf("output/edgeR/plotVolcano_res_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.pdf", width=3, height=4)    
EnhancedVolcano(res_tibble,
  lab = res_tibble$geneSymbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'OEKO vs WT, ESC, EZH2',
  pCutoff = 5e-2,         #
  FCcutoff = 0.58,
  pointSize = 1.0,
  labSize = 2,
  colCustom = keyvals,
  colAlpha = 1,
  legendPosition = 'none')  + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()


upregulated_genes <- sum(res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, na.rm = TRUE)
downregulated_genes <- sum(res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, na.rm = TRUE)

# Save as gene list for GO analysis:
upregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange > 0.58 & res_tibble$padj < 5e-2, ]
#### Filter for down-regulated genes
downregulated <- res_tibble[!is.na(res_tibble$log2FoldChange) & !is.na(res_tibble$padj) & res_tibble$log2FoldChange < -0.58 & res_tibble$padj < 5e-2, ]
#### Save
write.table(upregulated$geneSymbol, file = "output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(downregulated$geneSymbol, file = "output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange > 0.58)


res_tibble %>% dplyr::select(peakID, geneSymbol, log2FoldChange, padj) %>%
  filter(padj < 0.05, log2FoldChange < -0.58)




# Save as coordinates for bed deeptools

write_tsv( upregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.bed", col_names = FALSE)
write_tsv( downregulated %>%
  separate(peakID, into = c("chr", "start", "end"), sep = "_", remove = TRUE) %>%
  dplyr::select(chr, start, end), "output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.bed", col_names = FALSE)
```

--> For EZH2 I receive the error `Error in optimHess(par = init, fn = nbinomFn, gr = nbinomGr, x = x, y = y,  : non-finite value supplied by optim` when doing the `lfcShrink()`; error is discussed [here](https://www.biostars.org/p/9559740/) and is related to updated version of apeglm I need to install with `devtools::install_github("azhu513/apeglm")`
  --> For EZH2, let's not use `deseq2` env, but instead `deseq2V2` env!







## Complex heatmaps

Let's generate a heatmap that put together:
- Region gain/lost H3K27me3: only region associated with changes in a gene promoter
  - H3K27me3 signal for WT and KO
- RNAseq TPM value for the genes with H3K27me3 changes


### WT vs KO - Regions with H3K27me3 changes

```bash
conda activate deseq2

```

```R
# packages
library("tidyr")
library("ComplexHeatmap")
library("EnrichedHeatmap")
library("circlize")
library("GenomicRanges")
library("rtracklayer")


# ---- input paths  ----
bed_gain <- "output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.bed"     # Gain H3K27me3 (KO>WT)
bed_lost <- "output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.bed"  # Lost H3K27me3 (KO<WT)

bdg_WT <- "output/bigwig_Ferguson/ESC_WT_H3K27me3_unique_norm99_initialBigwig_median.sorted.bedGraph"
bdg_KO <- "output/bigwig_Ferguson/ESC_KO_H3K27me3_unique_norm99_initialBigwig_median.sorted.bedGraph"



# TPM table (GeneID, GeneSymbol, genotype, replicate, TPM, ENSG)
tpm_file <- "../019__RNAseq_ESC_V1/output/deseq2/txi_Kallisto-GeneLevel-TPM.txt"



# ---- read peak sets (BED -> GRanges) ----
read_bed3 <- function(path) {
  # 3-column BED (chr, start, end) without header
  df <- read_tsv(path, col_names = c("chr","start","end"), show_col_types = FALSE)
  GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
}

gr_gain <- read_bed3(bed_gain)
gr_lost <- read_bed3(bed_lost)
## Make one GRanges with a group label. Keep "Gain" first so it appears first in the heatmap.
regions_all <- c(gr_gain, gr_lost)
group <- factor(c(rep("Gain_H3K27me3", length(gr_gain)),
                  rep("Lost_H3K27me3", length(gr_lost))),
                levels = c("Gain_H3K27me3","Lost_H3K27me3"))

   
# ---- import bedGraphs ----
## bedGraph -> GRanges with a "score" column
cov_WT <- import(bdg_WT, format = "bedGraph")  # GRanges with score
cov_KO <- import(bdg_KO, format = "bedGraph")



# ---- make coverage functions for EnrichedHeatmap ----
## EnrichedHeatmap::normalizeToMatrix works with GRanges "signal" via cov2rle or providing cov as GRanges directly.
## We'll pass the GRanges bedGraph; EH will extract scores in windows.

target_center <- function(gr) {
  # center on peak middle (typical for histone marks)
  start <- floor((start(gr) + end(gr)) / 2)
  GRanges(seqnames = seqnames(gr), ranges = IRanges(start, width = 1))
}

centered_regions <- target_center(regions_all)

# parameters for the profile window
upstream   <- 2500   # bp
downstream <- 2500   # bp
bin_size   <- 50     # bp

mat_WT <- normalizeToMatrix(
  signal = cov_WT,
  target = centered_regions,
  extend = c(upstream, downstream),
  w = bin_size,
  mean_mode = "w0",   # window average (robust for bedGraph)
  empty_value = 0
)

mat_KO <- normalizeToMatrix(
  signal = cov_KO,
  target = centered_regions,
  extend = c(upstream, downstream),
  w = bin_size,
  mean_mode = "w0",
  empty_value = 0
)







# ---- RNA: KO vs WT log2FC per peak’s gene ----
## 1) summarize TPM per genotype (average across replicates)
tpm <- read_tsv(tpm_file, show_col_types = FALSE)
tpm_sum <- tpm %>%
  filter(genotype %in% c("WT","KO")) %>%
  group_by(GeneSymbol, genotype) %>%
  summarize(TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = genotype, values_from = TPM, values_fill = 0)

tpm_sum <- tpm_sum %>%
  mutate(RNA_log2FC_KOvsWT = log2((KO + 1e-6) / (WT + 1e-6)),
         WT_log2TPM = log2(WT + 1),
         KO_log2TPM = log2(KO + 1))

# 2) map peak -> geneSymbol from your DE result tibbles (upregulated / downregulated).
res_tibble <- read.delim("output/edgeR/DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt") %>%
  as_tibble() %>%
  mutate(peakID = as.character(peakID))

upregulated <- res_tibble %>%
  filter(!is.na(log2FoldChange), !is.na(padj),
         log2FoldChange > 0.58, padj < 5e-2)
downregulated <- res_tibble %>%
  filter(!is.na(log2FoldChange), !is.na(padj),
         log2FoldChange < -0.58, padj < 5e-2)

#    We need a vector of gene symbols aligned with 'regions_all' row order.
#    Build lookup: peakID -> geneSymbol
get_peak_gene <- function(df) {
  # df must have columns: peakID, geneSymbol
  df %>%
    select(peakID, geneSymbol) %>%
    distinct()
}

pk_gene_gain <- get_peak_gene(upregulated)
pk_gene_lost <- get_peak_gene(downregulated)


## Build character vector of peakIDs for gain & lost in the same order as BEDs
## We reconstruct peakIDs from BED 3-col by "chr_start_end" (your original convention)
mk_peak_id <- function(gr) paste0(as.character(seqnames(gr)), "_", start(gr), "_", end(gr))
peakID_gain <- mk_peak_id(gr_gain)
peakID_lost <- mk_peak_id(gr_lost)

# Join to fetch geneSymbol in the same order
gene_gain <- tibble(peakID = peakID_gain) %>%
  left_join(pk_gene_gain, by = "peakID") %>%
  pull(geneSymbol)

gene_lost <- tibble(peakID = peakID_lost) %>%
  left_join(pk_gene_lost, by = "peakID") %>%
  pull(geneSymbol)

genes_all <- c(gene_gain, gene_lost)


# 3) attach RNA metric to each peak row
rna_tbl <- tpm_sum %>% select(GeneSymbol, RNA_log2FC_KOvsWT, WT_log2TPM, KO_log2TPM)

rna_annot <- tibble(GeneSymbol = genes_all) %>%
  left_join(rna_tbl, by = "GeneSymbol")

RNA_log2FC <- rna_annot$RNA_log2FC_KOvsWT
WT_expr    <- rna_annot$WT_log2TPM
KO_expr    <- rna_annot$KO_log2TPM







## ---- colors ----
# signal color scale using actual range
sig_range <- range(c(as.numeric(mat_WT), as.numeric(mat_KO)), na.rm = TRUE)

# cap extreme max so color scale isn’t dominated by outliers
upper <- quantile(c(as.numeric(mat_WT), as.numeric(mat_KO)), 0.99, na.rm = TRUE)

sig_col <- circlize::colorRamp2(
  c(0, upper/2, upper),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

rna_col <- colorRamp2(c(-2, 0, 2), c("#2166AC","#F7F7F7","#B2182B"))

expr_col <- colorRamp2(quantile(c(WT_expr, KO_expr), c(0.05, 0.5, 0.95), na.rm = TRUE),
                       c("#2166AC","#F7F7F7","#B2182B"))


# ----- build robust color map from the actual matrix values -----
# use real dynamic range (capped at 99th pct)
sig_vals <- c(as.numeric(mat_WT), as.numeric(mat_KO))
sig_vals <- sig_vals[is.finite(sig_vals)]
upper <- quantile(sig_vals, 0.99, na.rm = TRUE)
mid   <- upper/2

sig_col <- colorRamp2(c(0, mid, upper),
                      c("#2166AC", "#F7F7F7", "#B2182B"))
legend_at <- c(0, round(mid, 1), round(upper, 1))

ht_WT <- EnrichedHeatmap(
  mat_WT,
  name = "WT_H3K27me3",
  col = sig_col,
  top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim = c(0, upper))),
  heatmap_legend_param = list(at = legend_at, labels = legend_at),
  row_split = group
)

ht_KO <- EnrichedHeatmap(
  mat_KO,
  name = "KO_H3K27me3",
  col = sig_col,
  top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim = c(0, upper))),
  heatmap_legend_param = list(at = legend_at, labels = legend_at)
)

## RNA as one or two side “heatmap columns”.
## Option A: one column with RNA KO vs WT log2FC
ht_RNAfc <- Heatmap(matrix(RNA_log2FC, ncol = 1),
                    name = "RNA log2FC\n(KO/WT)",
                    col = rna_col,
                    show_row_names = FALSE,
                    cluster_rows = FALSE,
                    width = unit(5, "mm"))

# Option B (uncomment to also show per-genotype expression columns)
# ht_RNA_WT <- Heatmap(matrix(WT_expr, ncol = 1), name = "WT log2(TPM+1)",
#                      col = expr_col, show_row_names = FALSE, cluster_rows = FALSE,
#                      width = unit(5, "mm"))
# ht_RNA_KO <- Heatmap(matrix(KO_expr, ncol = 1), name = "KO log2(TPM+1)",
#                      col = expr_col, show_row_names = FALSE, cluster_rows = FALSE,
#                      width = unit(5, "mm"))

# ---- draw ----
pdf("output/edgeR/EnrichedHeatmap-DESEQ2-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3-TPM.pdf", width = 8, height = 8)
draw(ht_WT + ht_KO + ht_RNAfc,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

```

--> BUG, it does not work...
  --> REquire troubleshooting


















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






### On bin1000space100_gt_pval05_padj001_fc1_avg30 - EZH2

Let's assign peak to genes on the two best windowns/parameters as in `001*/009*`:
- Bin 1000bp space 100bp, G test, pval 0.05 and padj 0.001, fc tresh 1 and minimum avg 30 in at least one genotype



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
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30 <- read.delim("output/diffreps/ESC_WTvsKO_EZH2_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) %>% 
  filter(padj < 0.001, abs(log2FC) > 1, Control.avg > 30 | Treatment.avg > 30)
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain = PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30 %>%
  filter(log2FC>1)
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost = PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30 %>%
  filter(log2FC<(-1))


PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30 <- read.delim("output/diffreps/ESC_WTvsOEKO_EZH2_unique_norm99_initialBigwig.bed-bin1000space100_gt_pval05-diff.nb.txt", sep = "\t", skip = 32, header = TRUE) %>%
  as_tibble() %>%
  dplyr::select(Chrom, Start, End, Length, Control.avg, Treatment.avg, log2FC, pval, padj) %>% 
  filter(padj < 0.001, abs(log2FC) > 1, Control.avg > 30 | Treatment.avg > 30)
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain = PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30 %>%
  filter(log2FC>1)
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost = PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30 %>%
  filter(log2FC<(-1))


### SAVE Gain and Lost peaks
write.table(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30, file="output/diffreps/PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain, file="output/diffreps/PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost, file="output/diffreps/PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost.txt", sep="\t", quote=F, row.names=F) 

write.table(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30, file="output/diffreps/PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain, file="output/diffreps/PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain.txt", sep="\t", quote=F, row.names=F) 
write.table(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost, file="output/diffreps/PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost.txt", sep="\t", quote=F, row.names=F) 
########

# Tidy peaks 
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_gr = makeGRangesFromDataFrame(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain,keep.extra.columns=TRUE)
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_gr = makeGRangesFromDataFrame(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost,keep.extra.columns=TRUE)
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_gr = makeGRangesFromDataFrame(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain,keep.extra.columns=TRUE)
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_gr = makeGRangesFromDataFrame(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost,keep.extra.columns=TRUE)

gr_list <- list(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain=PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_gr,PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost=PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_gr, PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain=PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_gr,PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost=PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_gr
)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_PSC_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30_initialBigwig.pdf", width = 16, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_PSC_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30_initialBigwig.pdf", width = 16, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot <- as.data.frame(peakAnnoList[["PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain"]]@anno)
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot <- as.data.frame(peakAnnoList[["PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost"]]@anno)
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot <- as.data.frame(peakAnnoList[["PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain"]]@anno)
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot <- as.data.frame(peakAnnoList[["PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost"]]@anno)

## Convert entrez gene IDs to gene symbols
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot$gene <- mapIds(org.Hs.eg.db, keys = PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot, file="output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot, file="output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot, file="output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot, file="output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot.txt", sep="\t", quote=F, row.names=F)  


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5 = tibble(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5 = tibble(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5 = tibble(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5 = tibble(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol = PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol = PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol = PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol = PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```


Signal changes promoter and 5'  Gain / Lost - EZH2: `bin1000space100_gt_pval05_padj001_fc1_avg30`:
- WT vs KO GAIN: gene (peak)= 6 (11)
- WT vs KO LOST: gene (peak)= 376 (543)
- WT vs OEKO GAIN: gene (peak)= 18 (27)
- WT vs OEKO LOST: gene (peak)= 530 (1295)




### On consensus peaks ESC_WTKOOEKO - H3K27me3 and EZH2


```bash
# files - consensus peaks H3K27me3 and EZH2 qval23 and 3 with merge100bp
output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed
output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge100bp.bed
output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed
output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge100bp.bed
```





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
ESC_WTKOOEKO_H3K27me3_qval23merge100bp <- read.delim("output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed", sep = "\t", header = FALSE) %>%
  as_tibble() %>%
  dplyr::rename("chr"= "V1", "start" = "V2", "end" = "V3")

ESC_WTKOOEKO_EZH2_qval23merge100bp <- read.delim("output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge100bp.bed", sep = "\t", header = FALSE) %>%
  as_tibble() %>%
  dplyr::rename("chr"= "V1", "start" = "V2", "end" = "V3")
  
ESC_WTKOOEKO_H3K27me3_qval3merge100bp <- read.delim("output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed", sep = "\t", header = FALSE) %>%
  as_tibble() %>%
  dplyr::rename("chr"= "V1", "start" = "V2", "end" = "V3")

ESC_WTKOOEKO_EZH2_qval3merge100bp <- read.delim("output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge100bp.bed", sep = "\t", header = FALSE) %>%
  as_tibble() %>%
  dplyr::rename("chr"= "V1", "start" = "V2", "end" = "V3")



# Tidy peaks 
ESC_WTKOOEKO_H3K27me3_qval23merge100bp_gr = makeGRangesFromDataFrame(ESC_WTKOOEKO_H3K27me3_qval23merge100bp,keep.extra.columns=TRUE)
ESC_WTKOOEKO_EZH2_qval23merge100bp_gr = makeGRangesFromDataFrame(ESC_WTKOOEKO_EZH2_qval23merge100bp,keep.extra.columns=TRUE)
ESC_WTKOOEKO_H3K27me3_qval3merge100bp_gr = makeGRangesFromDataFrame(ESC_WTKOOEKO_H3K27me3_qval3merge100bp,keep.extra.columns=TRUE)
ESC_WTKOOEKO_EZH2_qval3merge100bp_gr = makeGRangesFromDataFrame(ESC_WTKOOEKO_EZH2_qval3merge100bp,keep.extra.columns=TRUE)

gr_list <- list(ESC_WTKOOEKO_H3K27me3_qval23merge100bp=ESC_WTKOOEKO_H3K27me3_qval23merge100bp_gr,ESC_WTKOOEKO_EZH2_qval23merge100bp=ESC_WTKOOEKO_EZH2_qval23merge100bp_gr, ESC_WTKOOEKO_H3K27me3_qval3merge100bp=ESC_WTKOOEKO_H3K27me3_qval3merge100bp_gr,ESC_WTKOOEKO_EZH2_qval3merge100bp=ESC_WTKOOEKO_EZH2_qval3merge100bp_gr
)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_ESC_WTKOOEKO_H3K27me3EZH2_qval23qval3merge100bp.pdf", width = 16, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_ESC_WTKOOEKO_H3K27me3EZH2_qval23qval3merge100bp.pdf", width = 16, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot <- as.data.frame(peakAnnoList[["ESC_WTKOOEKO_H3K27me3_qval23merge100bp"]]@anno)
ESC_WTKOOEKO_EZH2_qval23merge100bp_annot <- as.data.frame(peakAnnoList[["ESC_WTKOOEKO_EZH2_qval23merge100bp"]]@anno)
ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot <- as.data.frame(peakAnnoList[["ESC_WTKOOEKO_H3K27me3_qval3merge100bp"]]@anno)
ESC_WTKOOEKO_EZH2_qval3merge100bp_annot <- as.data.frame(peakAnnoList[["ESC_WTKOOEKO_EZH2_qval3merge100bp"]]@anno)

## Convert entrez gene IDs to gene symbols
ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
ESC_WTKOOEKO_EZH2_qval23merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_EZH2_qval23merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
ESC_WTKOOEKO_EZH2_qval23merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_EZH2_qval23merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
ESC_WTKOOEKO_EZH2_qval3merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_EZH2_qval3merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
ESC_WTKOOEKO_EZH2_qval3merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_EZH2_qval3merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot, file="output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot, file="output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval23merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot, file="output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(ESC_WTKOOEKO_EZH2_qval3merge100bp_annot, file="output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval3merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5 = tibble(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5 = tibble(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot_promoterAnd5 = tibble(ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
ESC_WTKOOEKO_EZH2_qval3merge100bp_annot_promoterAnd5 = tibble(ESC_WTKOOEKO_EZH2_qval3merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5_geneSymbol = ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5_geneSymbol = ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot_promoterAnd5_geneSymbol = ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
ESC_WTKOOEKO_EZH2_qval3merge100bp_annot_promoterAnd5_geneSymbol = ESC_WTKOOEKO_EZH2_qval3merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval3merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(ESC_WTKOOEKO_EZH2_qval3merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval3merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```







### On EZH1 peaks; consensus and OEKO peaks - qval 2.3


```bash
# files - consensus peaks EZH1 qval2.3 merge100bp and OEKO peaks
output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.merge100bp.bed
output/macs2/broad/broad_blacklist_qval2.30103/ESC_OEKO_EZH1_pool_peaks.broadPeak
```





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
ESC_WTKOOEKO_EZH1_qval23merge100bp <- read.delim("output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH1_pool_peaks.sorted.merge100bp.bed", sep = "\t", header = FALSE) %>%
  as_tibble() %>%
  dplyr::rename("chr"= "V1", "start" = "V2", "end" = "V3")

ESC_OEKO_EZH1_qval23 <- read.delim("output/macs2/broad/broad_blacklist_qval2.30103/ESC_OEKO_EZH1_pool_peaks.broadPeak", sep = "\t", header = FALSE) %>%
  as_tibble() %>%
  dplyr::rename("chr"= "V1", "start" = "V2", "end" = "V3")
  

# Tidy peaks 
ESC_WTKOOEKO_EZH1_qval23merge100bp_gr = makeGRangesFromDataFrame(ESC_WTKOOEKO_EZH1_qval23merge100bp,keep.extra.columns=TRUE)
ESC_OEKO_EZH1_qval23_gr = makeGRangesFromDataFrame(ESC_OEKO_EZH1_qval23,keep.extra.columns=TRUE)

gr_list <- list(ESC_WTKOOEKO_EZH1_qval23merge100bp=ESC_WTKOOEKO_EZH1_qval23merge100bp_gr, ESC_OEKO_EZH1_qval23=ESC_OEKO_EZH1_qval23_gr
)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## plots
pdf("output/ChIPseeker/plotAnnoBar_ESC_WTKOOEKO_EZH1_qval23merge100bp-OEKO_EZH1_qval32.pdf", width = 16, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
pdf("output/ChIPseeker/plotDistToTSS_ESC_WTKOOEKO_EZH1_qval23merge100bp-OEKO_EZH1_qval32.pdf", width = 16, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()

## Get annotation data frame
ESC_WTKOOEKO_EZH1_qval23merge100bp_annot <- as.data.frame(peakAnnoList[["ESC_WTKOOEKO_EZH1_qval23merge100bp"]]@anno)
ESC_OEKO_EZH1_qval23_annot <- as.data.frame(peakAnnoList[["ESC_OEKO_EZH1_qval23"]]@anno)


## Convert entrez gene IDs to gene symbols
ESC_WTKOOEKO_EZH1_qval23merge100bp_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_EZH1_qval23merge100bp_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
ESC_WTKOOEKO_EZH1_qval23merge100bp_annot$gene <- mapIds(org.Hs.eg.db, keys = ESC_WTKOOEKO_EZH1_qval23merge100bp_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
ESC_OEKO_EZH1_qval23_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = ESC_OEKO_EZH1_qval23_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
ESC_OEKO_EZH1_qval23_annot$gene <- mapIds(org.Hs.eg.db, keys = ESC_OEKO_EZH1_qval23_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(ESC_WTKOOEKO_EZH1_qval23merge100bp_annot, file="output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH1_qval23merge100bp_annot.txt", sep="\t", quote=F, row.names=F)  
write.table(ESC_OEKO_EZH1_qval23_annot, file="output/ChIPseeker/annotation_ESC_OEKO_EZH1_qval23_annot.txt", sep="\t", quote=F, row.names=F)  


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5 = tibble(ESC_WTKOOEKO_EZH1_qval23merge100bp_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
ESC_OEKO_EZH1_qval23_annot_promoterAnd5 = tibble(ESC_OEKO_EZH1_qval23_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5_geneSymbol = ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
ESC_OEKO_EZH1_qval23_annot_promoterAnd5_geneSymbol = ESC_OEKO_EZH1_qval23_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(ESC_OEKO_EZH1_qval23_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annotation_ESC_OEKO_EZH1_qval23_annot_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```









# deepTools plots

## PEAKS 



```bash
conda activate deeptools






# Peak with DIFFREPS H3K27me3 changes bin1000space100_gt_pval05_padj001_fc1_avg100
## GAIN LOST PEAKS
### WT vs KO
output/diffreps/PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain.txt
output/diffreps/PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost.txt
output/diffreps/PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain.txt
output/diffreps/PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost.txt
### WT vs OEKO
output/diffreps/PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain.txt
output/diffreps/PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost.txt
output/diffreps/PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain.txt
output/diffreps/PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost.txt


# Peak with DESEQ2/MACS2 H3K27me3 changes
## GAIN LOST PEAKS
### WT vs KO
output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.bed
output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.bed
output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.bed
output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.bed
### WT vs OEKO
output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.bed
output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.bed
output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.bed
output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.bed


# Check signal WT vs KO regions with bigwig WT,KO,OEKO from FergusonUniqueNorm99
## DIFFREPS
sbatch scripts/matrix_PEAK_5kb-PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50662760 ok
sbatch scripts/matrix_PEAK_5kb-PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50742353 ok
### Check signal WT vs OEKO regions with bigwig WT,KO,OEKO from FergusonUniqueNorm99
sbatch scripts/matrix_PEAK_5kb-PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50662855 ok
sbatch scripts/matrix_PEAK_5kb-PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50742354 ok

## DESEQ2/MACS2 - all together
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51164829 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51164877 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51173339 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51166450 ok
## DESEQ2/MACS2 - IP per IP
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-H3K27me3.sh # 51174350 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2.sh # 51174488 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1.sh # 51174498 ok

sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-H3K27me3.sh # 51174613 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2.sh # 511746832 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1.sh # 511746839 ok


sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-H3K27me3.sh # 511746965 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2.sh # 51175086 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1.sh # 51175165 ok

sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-H3K27me3.sh # 51175284 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH2.sh # 51175395 ok
sbatch scripts/matrix_PEAK_5kb-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-q05fc058-WTKOOEKO-EZH1.sh # 51175525 ok




# check signal in MACS2 PEAKS
sbatch scripts/matrix_PEAK_5kb-macs2broad_WT_H3K27me3poolqval23-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50766577 ok
sbatch scripts/matrix_PEAK_5kb-macs2broad_WT_EZH2poolqval23-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50766992 ok
sbatch scripts/matrix_PEAK_5kb-macs2broad_OEKO_EZH1poolqval23-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50767482 ok
## consensus peaks
sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50895445 ok
sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-H3K27me3.sh # 51132562 ok
sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH2.sh # 51132659 ok
sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1.sh # 51132725 ok
sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_keepZero.sh # 51134009 ok
sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_minThreshold025.sh # 51398836 FAIL
sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_thresh2.sh # 51404482 ok
sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-EZH1_thresh1.sh # 51409782 xxx


sbatch scripts/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_EZH2poolqval23merge100bp-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50895514 ok


```






## GENES




```bash
# Generate gtf file from gene list:
## Gain Lost DIFFREPS H3K27me3
output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt
## Gain Lost DIFFREPS EZH2
output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt
## MACS2 consensus peaks qval2.3 merge100bp H3K27me3 EZH2
output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt
## DEGs from 001/019
output/deseq2/upregulated_q05fc058_ESC_KO_vs_ESC_WT.txt
output/deseq2/downregulated_q05fc058_ESC_KO_vs_ESC_WT.txt
output/deseq2/upregulated_q05fc058_ESC_OEKO_vs_ESC_WT.txt
output/deseq2/downregulated_q05fc058_ESC_OEKO_vs_ESC_WT.txt
## EZH1 peaks - consensus WTKOOEKO qval2.3 merge100bp - OEKO qval2.3 pool peaks
output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt
output/ChIPseeker/annotation_ESC_OEKO_EZH1_qval23_annot_promoterAnd5_geneSymbol.txt
# Peak with DESEQ2/MACS2 H3K27me3 changes
output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt
output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt
output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt
output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt
output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.txt
output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.txt
output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.txt
output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.txt


## put together Gain and Lost mix
cat output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt \
    output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt \
    | sort | uniq > output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5_geneSymbol.txt
cat output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt \
    output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt \
    | sort | uniq > output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5_geneSymbol.txt

cat output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol.txt \
    output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol.txt \
    | sort | uniq > output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5_geneSymbol.txt
cat output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol.txt \
    output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol.txt \
    | sort | uniq > output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5_geneSymbol.txt


### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt



sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt





sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5_as_gtf_geneSymbol.txt

sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt


sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5_as_gtf_geneSymbol.txt





sed 's/\r$//; s/.*/gene_name "&"/' ../019__RNAseq_ESC_V1/output/deseq2/upregulated_q05fc058_ESC_KO_vs_ESC_WT.txt > ../019__RNAseq_ESC_V1/output/deseq2/upregulated_q05fc058_ESC_KO_vs_ESC_WT_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../019__RNAseq_ESC_V1/output/deseq2/downregulated_q05fc058_ESC_KO_vs_ESC_WT.txt > ../019__RNAseq_ESC_V1/output/deseq2/downregulated_q05fc058_ESC_KO_vs_ESC_WT_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../019__RNAseq_ESC_V1/output/deseq2/upregulated_q05fc058_ESC_OEKO_vs_ESC_WT.txt > ../019__RNAseq_ESC_V1/output/deseq2/upregulated_q05fc058_ESC_OEKO_vs_ESC_WT_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../019__RNAseq_ESC_V1/output/deseq2/downregulated_q05fc058_ESC_OEKO_vs_ESC_WT.txt > ../019__RNAseq_ESC_V1/output/deseq2/downregulated_q05fc058_ESC_OEKO_vs_ESC_WT_as_gtf_geneSymbol.txt


sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/ChIPseeker/annotation_ESC_OEKO_EZH1_qval23_annot_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annotation_ESC_OEKO_EZH1_qval23_annot_promoterAnd5_as_gtf_geneSymbol.txt




sed 's/\r$//; s/.*/gene_name "&"/' output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt > output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.txt > output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt > output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.txt > output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.txt > output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.txt > output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.txt > output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.txt > output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2_as_gtf_geneSymbol.txt







## Filter the gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5.gtf

grep -Ff output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5.gtf



grep -Ff output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost_annot_promoterAnd5.gtf

grep -Ff output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5.gtf





grep -Ff output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5.gtf

grep -Ff output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5.gtf

grep -Ff output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost_annot_promoterAnd5.gtf

grep -Ff output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5.gtf




grep -Ff output/ChIPseeker/annotation_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC_WTKOOEKO_H3K27me3_qval23merge100bp_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC_WTKOOEKO_EZH2_qval23merge100bp_annot_promoterAnd5.gtf



grep -Ff ../019__RNAseq_ESC_V1/output/deseq2/upregulated_q05fc058_ESC_KO_vs_ESC_WT_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC001019_upregulated_q05fc058_ESC_KO_vs_ESC_WT.gtf
grep -Ff ../019__RNAseq_ESC_V1/output/deseq2/downregulated_q05fc058_ESC_KO_vs_ESC_WT_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC001019_downregulated_q05fc058_ESC_KO_vs_ESC_WT.gtf
grep -Ff ../019__RNAseq_ESC_V1/output/deseq2/upregulated_q05fc058_ESC_OEKO_vs_ESC_WT_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC001019_upregulated_q05fc058_ESC_OEKO_vs_ESC_WT.gtf
grep -Ff ../019__RNAseq_ESC_V1/output/deseq2/downregulated_q05fc058_ESC_OEKO_vs_ESC_WT_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC001019_downregulated_q05fc058_ESC_OEKO_vs_ESC_WT.gtf



grep -Ff output/ChIPseeker/annotation_ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC_WTKOOEKO_EZH1_qval23merge100bp_annot_promoterAnd5.gtf
grep -Ff output/ChIPseeker/annotation_ESC_OEKO_EZH1_qval23_annot_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC_OEKO_EZH1_qval23_annot_promoterAnd5.gtf





grep -Ff output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.gtf
grep -Ff output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.gtf
grep -Ff output/edgeR/upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.gtf
grep -Ff output/edgeR/downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.gtf
grep -Ff output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.gtf
grep -Ff output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.gtf
grep -Ff output/edgeR/upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.gtf
grep -Ff output/edgeR/downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.gtf









# DIFFREPS
## GAIN LOST GENES
### WT vs KO
meta/ENCFF159KBI_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5.gtf
meta/ENCFF159KBI_PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5.gtf
meta/ENCFF159KBI_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5.gtf
meta/ENCFF159KBI_PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5.gtf
### WT vs OEKO
meta/ENCFF159KBI_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Gain_annot_promoterAnd5.gtf
meta/ENCFF159KBI_PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__Lost_annot_promoterAnd5.gtf
meta/ENCFF159KBI_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Gain_annot_promoterAnd5.gtf
meta/ENCFF159KBI_PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__Lost_annot_promoterAnd5.gtf


# DESEQ2/MACS2
## GAIN LOST GENES
### WT vs KO
meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.gtf
meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-H3K27me3.gtf
meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.gtf
meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-EZH2.gtf
### WT vs OEKO
meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.gtf
meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-H3K27me3.gtf
meta/ENCFF159KBI_upregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.gtf
meta/ENCFF159KBI_downregulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-EZH2.gtf





# Check signal WT vs KO regions with bigwig WT,KO,OEKO from FergusonUniqueNorm99
## all genotypes and IP
sbatch scripts/matrix_GENETSS_5kb-PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50668786 ok
sbatch scripts/matrix_GENETSS_5kb-PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50752200 ok
## all genotypes but IP separated
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-H3K27me3.sh # 51084760 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-H3K27me3.sh # 51084806 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-EZH2.sh # 51084846 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-EZH2.sh # 51084927 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-EZH1.sh # 51084995 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-EZH1.sh # 51085074 xxx



# Check signal WT vs OEKO regions with bigwig WT,KO,OEKO from FergusonUniqueNorm99
## all genotypes and IP
sbatch scripts/matrix_GENETSS_5kb-PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50669023 ok
sbatch scripts/matrix_GENETSS_5kb-PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50752411 ok
## all genotypes but IP separated
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-H3K27me3.sh # 51085256 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-H3K27me3.sh # 51085328 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-EZH2.sh # 51085363 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-EZH2.sh # 51085392 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsOEKO_H3K27me3_bin1000space100_gt_pval05_padj001_fc1_avg100__GainLost-WTKOOEKO-EZH1.sh # 51085456 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-PSC_WTvsOEKO_EZH2_bin1000space100_gt_pval05_padj001_fc1_avg30__GainLost-WTKOOEKO-EZH1.sh # 51085472 xxx






# check signal in ALL GENES
sbatch scripts/matrix_GENETSS_5kb-ENCFF159KBI-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50768671 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-H3K27me3.sh # 50946785 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH2.sh # 50946787 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ENCFF159KBI-WTKOOEKO-EZH1.sh # 50946845 ok




# Check signal in MACS2 consensus peaks assign to genes - promoter and 5
## consensus H3K27me3/EZH2 all genotypes and IP
sbatch scripts/matrix_GENETSS_5kb-ESC_WTKOOEKO_H3K27me3_qval23merge100bp-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50900070 ok
sbatch scripts/matrix_GENETSS_5kb-ESC_WTKOOEKO_EZH2_qval23merge100bp-WTKOOEKO-H3K27me3EZH2EZH1.sh # 50900120 ok
## consensus H3K27me3/EZH2 all genotypes but IP separated
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_H3K27me3_qval23merge100bp-WTKOOEKO-H3K27me3.sh # 50947298 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_H3K27me3_qval23merge100bp-WTKOOEKO-EZH2.sh # 50947618 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_H3K27me3_qval23merge100bp-WTKOOEKO-EZH1.sh # 50947603 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH2_qval23merge100bp-WTKOOEKO-H3K27me3.sh # 50947622 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH2_qval23merge100bp-WTKOOEKO-EZH2.sh # 50947640 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH2_qval23merge100bp-WTKOOEKO-EZH1.sh # 50947650 ok


# Check signal in DEGs from 001/019
## all together
sbatch scripts/matrix_GENETSS_5kb-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51081723 ok
sbatch scripts/matrix_GENETSS_5kb-ESC001019_UpDownregulated_q05fc058_ESC_OEKO_vs_ESC_WT-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51081814 ok
## all genotypes but IP separated
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-H3K27me3.sh # 51081895 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.sh # 51081974 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_KO_vs_ESC_WT-WTKOOEKO-EZH1.sh # 51081979 ok

sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_OEKO_vs_ESC_WT-WTKOOEKO-H3K27me3.sh # 51082048 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_OEKO_vs_ESC_WT-WTKOOEKO-EZH2.sh # 51082118 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC001019_UpDownregulated_q05fc058_ESC_OEKO_vs_ESC_WT-WTKOOEKO-EZH1.sh # 51082139 ok


# Check signal in EZH1 peaks - consensus WTKOOEKO qval2.3 merge100bp - OEKO qval2.3 pool peaks - promoter and 5
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp-WTKOOEKO-H3K27me3.sh # 51129082 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp-WTKOOEKO-EZH2.sh # 51129343 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_WTKOOEKO_EZH1_qval23merge100bp-WTKOOEKO-EZH1.sh # 51129393 xxx
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_OEKO_EZH1_qval23-WTKOOEKO-H3K27me3.sh # 51129429 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_OEKO_EZH1_qval23-WTKOOEKO-EZH2.sh # 51129474 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-ESC_OEKO_EZH1_qval23-WTKOOEKO-EZH1.sh # 51129494 xxx








# DESEQ2/MACS2
## GAIN LOST GENES
### WT vs KO
## all genotypes and IP
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51178250 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51178361 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51178376 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-WTKOOEKO-H3K27me3EZH2EZH1.sh # 51178511 ok
## all genotypes but IP separated
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-H3K27me3.sh # 51179348 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.sh # 51179468 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH1.sh # 51179472 ok

sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-WTKOOEKO-H3K27me3.sh # 51179598 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-WTKOOEKO-EZH2.sh # 51179714 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_H3K27me3_qval23merge100bp-ESC_OEKO_vs_ESC_WT-WTKOOEKO-EZH1.sh # 51179723 ok

sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-H3K27me3.sh # 51182028 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH2.sh # 51179989 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_KO_vs_ESC_WT-WTKOOEKO-EZH1.sh # 51180108 ok

sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-WTKOOEKO-H3K27me3.sh # 51180740 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-WTKOOEKO-EZH2.sh # 51180889 ok
sbatch scripts/matrix_GENETSSTES_250bp100bp-regulated_q05fc058-WTKOOEKO_EZH2_qval23merge100bp-ESC_OEKO_vs_ESC_WT-WTKOOEKO-EZH1.sh # 51180904 ok





```

--> Changes of H3K27me3 is clear; but EZH2 does not clearly follow these changes.








