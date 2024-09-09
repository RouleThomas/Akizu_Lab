# Project

**H9 cell lines**

- ESC native; 
    - WT: SUZ12, EZH2, EZH1, H3K27me3
    - KOEF1aEZH1: SUZ12, EZH2, EZH1, H3K27me3
    - WTEF1aEZH1: SUZ12, EZH2, EZH1, H3K27me3
--> All in 2 Bio Rep

- ~34days Neurons;
    - WT: SUZ12, H3K27me3


**Objectives:**
- PSC experiment can be used to check the competition hypotheses between EZH1 EZH2 (see poster with Jasmine) + Spreading of H3K27me3
- NEU experiment to check whether CutRun in neurons work!



--> **Working samples can be added as aditional WT replicate!**




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
# Following email instructions
wget -r -b -c --user=X202SC24076724-Z01-F001 --password=jxc30ah4 ftp://usftp21.novogene.com:21/
```

--> All good, files created in `usftp21.novogene.com/`




# Rename file


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






# Fastp cleaning

```bash
sbatch scripts/fastp_1.sh # 24136634 ok
sbatch scripts/fastp_2.sh # 24136758 ok

```


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:24136634 scripts/bowtie2_1.sh # 24137045 ok
sbatch --dependency=afterany:24136758 scripts/bowtie2_2.sh # 24137877 partial fail
sbatch scripts/bowtie2_missing1.sh # 24306442 ok
```

--> Looks good; overall ~75% uniquely aligned reads 


**Mapping on E coli**

```bash
conda activate bowtie2

sbatch scripts/bowtie2_MG1655_1.sh # 24385988 xxx
sbatch scripts/bowtie2_MG1655_2.sh # 24386097 xxx

```

--> between 0.1 - 2% uniquely aligned reads (not a lot..; previously `005__CutRun` 10% (in `003__CutRun` was less than 1%) )





## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-24137045.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_24137045.txt

for file in slurm-24306442.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_24306442.txt

for file in slurm-24137877.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_24137877.txt
```

Add these values to `/home/roulet/001_EZH1_project/013__CutRun_PSC_native/samples_001013.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >75% input reads as been uniquely mapped to the genome (90% non uniq) 



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.




```bash
conda activate bowtie2

sbatch --dependency=afterany:24137045 scripts/samtools_unique_1.sh # 24137245 ok
sbatch --dependency=afterany:24137877 scripts/samtools_unique_2.sh # 24137923 partial fail

sbatch --dependency=afterany:24306442 scripts/samtools_missing1.sh # 24306518 ok

```



Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch --dependency=afterany:24385988 scripts/samtools_MG1655_unique_1.sh # 24387652 ok
sbatch --dependency=afterany:24386097 scripts/samtools_MG1655_unique_2.sh # 24387686 ok

# count the nb of uniq aligned reads
sbatch scripts/samtools_MG1655_unique_Count_1.sh # 26299778 ok
sbatch scripts/samtools_MG1655_unique_Count_2.sh # 26299851 ok




```

--> **IMPORTANT NOTE!!! I should use the MG1655 scaling factor by calculating uniq. aligned MG1655 / uniq. aligned human!! Before I was using non uniq. aligned reads for MG1655 scaling factor...**




# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:24137245 scripts/bamtobigwig_unique_1.sh # 24137373 ok
sbatch --dependency=afterany:24137923 scripts/bamtobigwig_unique_2.sh # 24138000 partial fail

sbatch --dependency=afterany:24306518 scripts/bamtobigwig_missing1.sh # 24306726 ok


```

- **NEU**
*Pass*: none
*Failed*: H3K27me3, SUZ12

- **ESC _ WT**
*Pass*: H3K27me3, SUZ12, EZH2
*Failed*: EZH1

- **ESC _ KOEF1aEZH1**
*Pass*: H3K27me3, SUZ12, EZH, EZH1
*Failed*: none

- **ESC _ WTEF1aEZH1**
*Pass*: H3K27me3, SUZ12, EZH, EZH1
*Failed*: none

--> Replicate 2 perform less well for ESC


## Pearson correlation heatmap on bigwig signals




```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all.sh # 24387938 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels NEU_WT_H3K27me3_R1 NEU_WT_IGG_R1 NEU_WT_SUZ12_R1 PSC_KO_EZH1_R1 PSC_KO_EZH1_R2 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KOEF1aEZH1_EZH1_R1 PSC_KOEF1aEZH1_EZH1_R2 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_KOEF1aEZH1_SUZ12_R1 PSC_KOEF1aEZH1_SUZ12_R2 PSC_WT_EZH1_R1 PSC_WT_EZH1_R2 PSC_WT_EZH2_R1 PSC_WT_EZH2_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WT_IGG_R1 PSC_WT_IGG_R2 PSC_WT_SUZ12_R1 PSC_WT_SUZ12_R2 PSC_WTEF1aEZH1_EZH1_R1 PSC_WTEF1aEZH1_EZH1_R2 PSC_WTEF1aEZH1_EZH2_R1 PSC_WTEF1aEZH1_EZH2_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_SUZ12_R1 PSC_WTEF1aEZH1_SUZ12_R2 \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels NEU_WT_H3K27me3_R1 NEU_WT_IGG_R1 NEU_WT_SUZ12_R1 PSC_KO_EZH1_R1 PSC_KO_EZH1_R2 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KOEF1aEZH1_EZH1_R1 PSC_KOEF1aEZH1_EZH1_R2 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_KOEF1aEZH1_SUZ12_R1 PSC_KOEF1aEZH1_SUZ12_R2 PSC_WT_EZH1_R1 PSC_WT_EZH1_R2 PSC_WT_EZH2_R1 PSC_WT_EZH2_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WT_IGG_R1 PSC_WT_IGG_R2 PSC_WT_SUZ12_R1 PSC_WT_SUZ12_R2 PSC_WTEF1aEZH1_EZH1_R1 PSC_WTEF1aEZH1_EZH1_R2 PSC_WTEF1aEZH1_EZH2_R1 PSC_WTEF1aEZH1_EZH2_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_SUZ12_R1 PSC_WTEF1aEZH1_SUZ12_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf





# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_PSC.sh # 24388068 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_PSC.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KO_EZH1_R1 PSC_KO_EZH1_R2 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KOEF1aEZH1_EZH1_R1 PSC_KOEF1aEZH1_EZH1_R2 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_KOEF1aEZH1_SUZ12_R1 PSC_KOEF1aEZH1_SUZ12_R2 PSC_WT_EZH1_R1 PSC_WT_EZH1_R2 PSC_WT_EZH2_R1 PSC_WT_EZH2_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WT_IGG_R1 PSC_WT_IGG_R2 PSC_WT_SUZ12_R1 PSC_WT_SUZ12_R2 PSC_WTEF1aEZH1_EZH1_R1 PSC_WTEF1aEZH1_EZH1_R2 PSC_WTEF1aEZH1_EZH2_R1 PSC_WTEF1aEZH1_EZH2_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_SUZ12_R1 PSC_WTEF1aEZH1_SUZ12_R2 \
    -o output/bigwig/multiBigwigSummary_PSC_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_PSC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KO_EZH1_R1 PSC_KO_EZH1_R2 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KOEF1aEZH1_EZH1_R1 PSC_KOEF1aEZH1_EZH1_R2 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_KOEF1aEZH1_SUZ12_R1 PSC_KOEF1aEZH1_SUZ12_R2 PSC_WT_EZH1_R1 PSC_WT_EZH1_R2 PSC_WT_EZH2_R1 PSC_WT_EZH2_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WT_IGG_R1 PSC_WT_IGG_R2 PSC_WT_SUZ12_R1 PSC_WT_SUZ12_R2 PSC_WTEF1aEZH1_EZH1_R1 PSC_WTEF1aEZH1_EZH1_R2 PSC_WTEF1aEZH1_EZH2_R1 PSC_WTEF1aEZH1_EZH2_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_SUZ12_R1 PSC_WTEF1aEZH1_SUZ12_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_PSC_heatmap.pdf




# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_PSC_WT.sh # 24388133 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_PSC_WT.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_WT_EZH1_R1 PSC_WT_EZH1_R2 PSC_WT_EZH2_R1 PSC_WT_EZH2_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WT_IGG_R1 PSC_WT_IGG_R2 PSC_WT_SUZ12_R1 PSC_WT_SUZ12_R2 \
    -o output/bigwig/multiBigwigSummary_PSC_WT_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_PSC_WT.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_WT_EZH1_R1 PSC_WT_EZH1_R2 PSC_WT_EZH2_R1 PSC_WT_EZH2_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WT_IGG_R1 PSC_WT_IGG_R2 PSC_WT_SUZ12_R1 PSC_WT_SUZ12_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_PSC_WT_heatmap.pdf





# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_PSC_KO.sh # 24388181 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_PSC_KO.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KO_EZH1_R1 PSC_KO_EZH1_R2 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 \
    -o output/bigwig/multiBigwigSummary_PSC_KO_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_PSC_KO.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KO_EZH1_R1 PSC_KO_EZH1_R2 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_PSC_KO_heatmap.pdf





# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_PSC_KOEF1aEZH1.sh # 24388236 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_PSC_KOEF1aEZH1.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KOEF1aEZH1_EZH1_R1 PSC_KOEF1aEZH1_EZH1_R2 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_KOEF1aEZH1_SUZ12_R1 PSC_KOEF1aEZH1_SUZ12_R2 \
    -o output/bigwig/multiBigwigSummary_PSC_KOEF1aEZH1_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_PSC_KOEF1aEZH1.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KOEF1aEZH1_EZH1_R1 PSC_KOEF1aEZH1_EZH1_R2 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_KOEF1aEZH1_SUZ12_R1 PSC_KOEF1aEZH1_SUZ12_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_PSC_KOEF1aEZH1_heatmap.pdf





# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_PSC_WTEF1aEZH1.sh # 24388291 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_PSC_WTEF1aEZH1.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_WTEF1aEZH1_EZH1_R1 PSC_WTEF1aEZH1_EZH1_R2 PSC_WTEF1aEZH1_EZH2_R1 PSC_WTEF1aEZH1_EZH2_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_SUZ12_R1 PSC_WTEF1aEZH1_SUZ12_R2 \
    -o output/bigwig/multiBigwigSummary_PSC_WTEF1aEZH1_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_PSC_WTEF1aEZH1.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_WTEF1aEZH1_EZH1_R1 PSC_WTEF1aEZH1_EZH1_R2 PSC_WTEF1aEZH1_EZH2_R1 PSC_WTEF1aEZH1_EZH2_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_SUZ12_R1 PSC_WTEF1aEZH1_SUZ12_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_PSC_WTEF1aEZH1_heatmap.pdf




# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_all.sh # 26305457 ok
# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_all.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KOEF1aEZH1_EZH1_R1 PSC_KOEF1aEZH1_EZH1_R2 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_KOEF1aEZH1_SUZ12_R1 PSC_KOEF1aEZH1_SUZ12_R2 PSC_WTEF1aEZH1_EZH1_R1 PSC_WTEF1aEZH1_EZH1_R2 PSC_WTEF1aEZH1_EZH2_R1 PSC_WTEF1aEZH1_EZH2_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_SUZ12_R1 PSC_WTEF1aEZH1_SUZ12_R2 \
    -o output/bigwig/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_all_plotPCA.pdf
## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KOEF1aEZH1_EZH1_R1 PSC_KOEF1aEZH1_EZH1_R2 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_KOEF1aEZH1_SUZ12_R1 PSC_KOEF1aEZH1_SUZ12_R2 PSC_WTEF1aEZH1_EZH1_R1 PSC_WTEF1aEZH1_EZH1_R2 PSC_WTEF1aEZH1_EZH2_R1 PSC_WTEF1aEZH1_EZH2_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_SUZ12_R1 PSC_WTEF1aEZH1_SUZ12_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_all_heatmap.pdf





# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_H3K27me3.sh # 26305500 ok
# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_H3K27me3.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 \
    -o output/bigwig/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_H3K27me3_plotPCA.pdf
## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_H3K27me3.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_KOEF1aEZH1_IGG_R1 PSC_KOEF1aEZH1_IGG_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 PSC_WTEF1aEZH1_IGG_R1 PSC_WTEF1aEZH1_IGG_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_WTEF1aEZH1vsKOEF1aEZH1_H3K27me3_heatmap.pdf






# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all_H3K27me3.sh # 26329533 ok


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all_H3K27me3.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 \
    --colors red red blue blue black black darkblue darkblue \
    --markers o o x x o o x x \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all_H3K27me3.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_heatmap.pdf






# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all_H3K27me3_DiffBindMG1655TMM.sh # 26396479 ok

# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all_H3K27me3_DiffBindMG1655TMM.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 \
    --colors red red blue blue black black darkblue darkblue \
    --markers o o x x o o x x \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_DiffBindMG1655TMM_plotPCA.pdf

plotPCA -in output/bigwig/multiBigwigSummary_all_H3K27me3_DiffBindMG1655TMM.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 \
    --colors red red blue blue black black darkblue darkblue \
    --markers o x o x o x o x \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_DiffBindMG1655TMM_plotPCA_2.pdf


## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all_H3K27me3_DiffBindMG1655TMM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_DiffBindMG1655TMM_heatmap.pdf






# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all_H3K27me3_TMM.sh # 26406742 ok

# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all_H3K27me3_TMM.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 \
    --colors red red blue blue black black darkblue darkblue \
    --markers o o x x o o x x \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_TMM_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all_H3K27me3_TMM.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KO_H3K27me3_R1 PSC_KO_H3K27me3_R2 PSC_KOEF1aEZH1_H3K27me3_R1 PSC_KOEF1aEZH1_H3K27me3_R2 PSC_WT_H3K27me3_R1 PSC_WT_H3K27me3_R2 PSC_WTEF1aEZH1_H3K27me3_R1 PSC_WTEF1aEZH1_H3K27me3_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_TMM_heatmap.pdf



```

--> Very noisy, not informative. 
    --> H3K27me3 cluster well together, then that is a mixed of samples/IP, even when checking per genotype

--> WTEF1aEZH1 and KOEF1aEZH1 look to cluster well together


# MACS2 peak calling on bam unique


--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad` and `narrow` **


```bash
conda activate macs2
# broad
sbatch scripts/macs2_broad_1.sh # 24400532 ok
sbatch scripts/macs2_broad_2.sh # 24400534 ok

# narrow
#--> NOT NEEDED


```






--> NEU no peak

--> PSC Rep1 work better than Rep2 (because less input material in Rep2)

**broad**: Check ppt for nb of peaks











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






# Spike in factor


## Calculate histone content

--> This histone content will be used to generate a scaling factor which will be used to histone-scaled our library size. The calcul/method to follow is from `003__CutRun/output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt`

**Pipeline:**
- Count the histone barcode on the clean reads
- Calculate SF (group by sample (replicate) and AB and calculate the total nb of reads. Then proportion of reads = nb read in sample / total reads. SF = min(proportion) / sample proportion)


## Count the histone barcode on the clean reads

Only for all H3K27me3 samples


```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_1.sh # 26196898 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_2.sh # 26197353 ok
```

--> It output the nb of reads found for each histone; then simply copy paste to the excell file `output/spikein/SpikeIn_QC_fastp_001013.xlsx` in GoogleDrive

- **PSC H3K27me3**: All enriched
- **NEU H3K27me3**: Not enriched (very lowly slightly for Reads2, not at all for Reads1)



## histone spike in factor


--> SF only calculating on WT and KO as KOEF1aEZH1 is NOT overexpressing..

```R
# package
library("tidyverse")
library("readxl")
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_001013.xlsx") 

## H3K27me3 with only WTEF1aEZH1 and KOEF1aEZH1
spikein_H3K27me3 = spikein %>%
    filter(Target == "H3K27me3",
    sample_ID %in% c("PSC_WTEF1aEZH1_H3K27me3_R1","PSC_WTEF1aEZH1_H3K27me3_R2", "PSC_KOEF1aEZH1_H3K27me3_R1", "PSC_KOEF1aEZH1_H3K27me3_R2")) %>%
    group_by(sample_ID, AB) %>%
    summarise(aligned=sum(counts))
# Total reads per IP
spikein_H3K27me3_total = spikein_H3K27me3 %>%
    ungroup() %>%
    group_by(AB) %>%
    mutate(total = sum(aligned)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_H3K27me3_read_prop = spikein_H3K27me3 %>%
    left_join(spikein_H3K27me3_total) %>%
    mutate(read_prop = aligned / total)
spikein_H3K27me3_read_prop_min = spikein_H3K27me3_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_H3K27me3_scaling_factor = spikein_H3K27me3_read_prop %>%
    left_join(spikein_H3K27me3_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_H3K27me3_scaling_factor, file="output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp-WTEF1aEZH1_KOEF1aEZH1.txt", sep="\t", quote=FALSE, row.names=FALSE)



## Collect total nb of spikein H3K27me3 reads
spikein_H3K27me3 = spikein %>%
    filter(Target == "H3K27me3") %>%
    group_by(sample_ID, AB) %>%
    summarise(aligned=sum(counts))


```

--> All good; higher SF for WTEF1aEZH2 as compared to KOEF1aEZH1




### Quality control plot

Then look at the xlsx file from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) to generate quality control plot. Use R cluster for vizualization (file is `spikein_QC.xlsx` in Google Drive), file in `output/spikein`.
```R
# package
library("tidyverse")
library("readxl")
# import df adn tidy to remove AB used in sample_ID
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_001013.xlsx") %>%
  separate(sample_ID, into = c("type", "condition", "tag"), sep = "_") %>%
  mutate(sample_ID = paste(type, condition, sep = "_")) %>%
  select(-type, -condition, -tag)


# data processing
spikein_sum_Barcode_read <- spikein %>%
	select(-Barcode) %>%
  	group_by(sample_ID, Target, AB) %>% 
	summarise(sum_read=sum(counts)) %>%
	unique() 

spikein_sum_Total_read <- spikein %>%
	select(-Barcode, -Target) %>%
  	group_by(sample_ID, AB) %>% 
	summarise(total_read=sum(counts)) %>%
	unique() 	

spikein_all <- spikein_sum_Barcode_read %>%
	left_join(spikein_sum_Total_read) %>%
	mutate(target_norm = (sum_read/total_read)*100)

## Histone scaling for H3K27me3
spikein_all_scale = spikein_all %>%
  group_by(sample_ID) %>%
  # Find the target_norm value when Target is H3K27me3 and AB is H3K27me3
  mutate(scaling_factor = ifelse(Target == "H3K27me3" & AB == "H3K27me3", target_norm, NA)) %>%
  # Fill the scaling_factor column with the appropriate value within each group
  fill(scaling_factor, .direction = "downup") %>%
  # Scale the target_norm values
  mutate(scaled_target_norm = target_norm / scaling_factor * 100) %>%
  # Remove the scaling_factor column
  select(-scaling_factor) %>%
  # Ungroup the data
  ungroup()
# Plot
pdf("output/spikein/QC_histone_spike_in_H3K27me3.pdf", width = 10, height = 3)
spikein_all_scale %>%
    filter(
           AB %in% c("H3K27me3", "IGG")) %>%
        ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
        geom_col(position = "dodge") +
        facet_wrap(~sample_ID, nrow=1) +
        geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


```


--> H3K27me3 is enriched for PSC. For NEU, very lowly...






# Ecoli scaling factor
## Mapping E coli
- Map our reads to the E. coli genome using same parameters as for human.
- Count the number of aligned reads to the spike-in control sequences for each sample `samtools view -S -F 4 -c sample_h3k27me3_rep1_spikein.sam > sample_h3k27me3_rep1_spikein_count.txt` (code in `bowtie2_MG1655.sh`)
- Do the math for scaling factor, same method as when using histone spike-in

--> **IMPORTANT NOTE: NEED to use the uniq. aligned reads for both MG1655 and hg38!!!**


```R
# package
library("tidyverse")
library("readxl")
library("ggpubr")

# all samples


# import df
spikein <- read_excel("output/spikein/SpikeIn_MG1655unique_001013.xlsx") 

## PSC
spikein <- read_excel("output/spikein/SpikeIn_MG1655unique_001013.xlsx") %>%
    filter(tissue == "PSC") %>%
    dplyr::select(-tissue)
# Total reads per IP
spikein_H3K27me3_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_H3K27me3_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655unique_PSC_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)



# WTEF1aEZH1 and KOEF1aEZH1 H3K27me3 only
spikein <- read_excel("output/spikein/SpikeIn_MG1655unique_001013.xlsx") %>%
    filter(sample_ID %in% c(
    "PSC_WTEF1aEZH1_H3K27me3_R1", "PSC_WTEF1aEZH1_H3K27me3_R2", "PSC_KOEF1aEZH1_H3K27me3_R1" ,"PSC_KOEF1aEZH1_H3K27me3_R2" ,"PSC_WTEF1aEZH1_IGG_R1",  "PSC_WTEF1aEZH1_IGG_R2", "PSC_KOEF1aEZH1_IGG_R1", "PSC_KOEF1aEZH1_IGG_R2"
    )) %>%
    dplyr::select(-tissue)
# Total reads per IP
spikein_H3K27me3_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_H3K27me3_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655unique_PSC_scaling_factor_WTEF1aEZH1andKOEF1aEZH1.txt", sep="\t", quote=FALSE, row.names=FALSE)

```

--> No need to calculate SF separately here! We obtain same values.



## Spike in Diffbind calcluation
### Histone spike in


**Using our scaling factor, let's estimate the 'new' library size** and provide it to `dba.normalize(library = c(1000, 12000))` = Like that our library size will be change taking into account our scaling factor! **Then we can normalize with library-size, RLE or TMM**... (issue discussed [here](https://support.bioconductor.org/p/9147040/)) 


### Adjust library size with histone scaling factor and apply normalization
Total number of reads is our library size (used samtools flagstat to double check) :

`samtools flagstat output/bowtie2/*unique.dupmark.sorted.bam` used to obtain library size (first value=library size)
--> Values save in GoogleDrive `007__*/sample_007.xlsx`. Histone-norm-library-size = library-size * SF. Using the non-reciprocal scaling factor, we increase the library-size; the more histone enriched, the more library size is increased, thus the more signal will decrease.

Now let's use these new histone-scaled library size and normalize with library-size,TMM or RLE. Let's use the **unique bam files** together with the **unique bam MACS2 raw files (xlsx, not the bed with pre-filtered qvalue)**

***Key points:***
- **Let's do 1 DiffBind per AB and tissue; otherwise the TMM normalization may take all, unrelated, samples into account!** --> Files are `meta_sample_macs2raw_unique*.txt`


```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 

# PSC WTEF1aEZH1 and KOEF1aEZH1 H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3_histoneSF.txt", header = TRUE, sep = "\t"))

### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)

## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3.RData")

### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()

pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist

sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)

### TMM 

sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(19509878, 21885102, 13598489,8632686), normalize = DBA_NORM_TMM) 

#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3.txt")


### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3_histoneSF.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3_histoneSF.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()





```



Now let's do the same method but using the **MG1655_library_scaled information and collect new MG1655_DiffBind_TMM_SF**:


```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 



load("output/DiffBind/sample_count_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3.RData")


### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist

sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibMG1655Scaled_TMM = dba.normalize(sample_count_blackgreylist, library = c( 13815063,18098857,12434028,15306766), normalize = DBA_NORM_TMM) 

#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibMG1655Scaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibMG1655Scaled_TMM, bRetrieve=TRUE)

console_output <- capture.output(print(sample_count_blackgreylist_LibMG1655Scaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3.txt")


### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3_MG1655SF.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibMG1655Scaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_WTEF1aEZH1KOEF1aEZH1_H3K27me3_MG1655SF.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibMG1655Scaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()






# WT and KO
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_WTKO_H3K27me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_WTKO_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_WTKO_H3K27me3.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_WTKO_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_WTKO_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_WTKO_H3K27me3_rep.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_TREATMENT, label=DBA_REPLICATE)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)

### TMM 
sample_count_blackgreylist_LibMG1655Scaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(13124625,14587838,13765544,9093591), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibMG1655Scaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibMG1655Scaled_TMM, bRetrieve=TRUE)

console_output <- capture.output(print(sample_count_blackgreylist_LibMG1655Scaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_PSC_WTKO_H3K27me3.txt")

### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_WTKO_H3K27me3_MG1655SF.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibMG1655Scaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_WTKO_H3K27me3_MG1655SF.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibMG1655Scaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_WTKO_H3K27me3_MG1655SF_rep.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibMG1655Scaled_TMM,DBA_TREATMENT, label=DBA_REPLICATE)
dev.off()








# all H3K27me3 samples (just to check whether how OEEZH1 genotypes cluster with the other samples)
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_all_H3K27me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_all_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_all_H3K27me3.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_all_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_all_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_all_H3K27me3_rep.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_TREATMENT, label=DBA_REPLICATE)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)

### TMM 
sample_count_blackgreylist_LibMG1655Scaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(13815063,18098857,12434028,15306766,13124625,14587838,13765544,9093591), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibMG1655Scaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibMG1655Scaled_TMM, bRetrieve=TRUE)

console_output <- capture.output(print(sample_count_blackgreylist_LibMG1655Scaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_PSC_all_H3K27me3.txt")

### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_all_H3K27me3_MG1655SF.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibMG1655Scaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_all_H3K27me3_MG1655SF.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibMG1655Scaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_all_H3K27me3_MG1655SF_rep.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibMG1655Scaled_TMM,DBA_TREATMENT, label=DBA_REPLICATE)
dev.off()

```



--> SF from histone and MG1655 are very comparable






# THOR diff peaks

Let's use THOR, notably to have IGG scaled bigwig...!

Comparison to do:
- H3K27me3 in WTEF1aEZH1 vs KOEF1aEZH1 --> Check whether they could be considered similar; if yes we will have 3 replicate for EZH1 overexpression (let's do without SF)


--> SF to use in THOR are the **reciprocal of MG1655_DiffBind_TMM** (histone are the same, so lets use the MG1655 ones)
--> Configs file created manually as `output/THOR/50dN_H3K27me3_WTvsKO.config` and `output/THOR/50dN_H3K27me3_WTvsKOEF1aEZH1.config` 



## Run THOR

*THOR is very buggy to make it work I need to temporaly change where to look for libraries lol.. So cannot use nano anymore for example...*

*Follow these parameters: `WTvsHET_unique_Keepdup` (perform best in previous CutRun)*

```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge

# AB per AB (TMM normalization from THOR)
sbatch scripts/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_TMM.sh # 26305339 ok
sbatch scripts/THOR_PSC_H3K27me3_WTvsKO_TMM.sh # 26331297 ok


# AB per AB (DiffBind SF TMM)
sbatch scripts/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_DiffBindMG1655TMM.sh # 26328774 ok
sbatch scripts/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_DiffBindhistoneTMM.sh # 26328776 ok
sbatch scripts/THOR_PSC_H3K27me3_WTvsKO_DiffBindMG1655TMM.sh # 26331379 ok

```


--> **WT vs KO; almost no diff peaks!** Could be true, as in KO only EZH2 is active, no compensation, EZH2 is already there working; EZH1 is NOT the main catalytic subunit in ESC, thus KO it do nothing






Generate median tracks:
```bash
conda activate BedToBigwig

sbatch scripts/bigwigmerge_THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_TMM.sh # 26308902 ok
```


## Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks


```R

# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# WTEF1aEZH1vsKOEF1aEZH1 TMM THOR
diffpeaks <- read_tsv("output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WTEF1aEZH1", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WTEF1aEZH1, into = c("count_WTEF1aEZH1_1","count_WTEF1aEZH1_2"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2) / (count_WTEF1aEZH1_1+count_WTEF1aEZH1_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WTEF1aEZH1 vs KOEF1aEZH1") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/log2FC_qval30.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WTEF1aEZH1 vs KOEF1aEZH1") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/THOR_qval50.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 30) %>%
  group_by(X6) %>%
  summarise(n = n())




# WTEF1aEZH1vsKOEF1aEZH1 DiffBindMG1655TMM 
diffpeaks <- read_tsv("output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_DiffBindMG1655TMM/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1DiffBindMG1655TMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WTEF1aEZH1", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WTEF1aEZH1, into = c("count_WTEF1aEZH1_1","count_WTEF1aEZH1_2"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2) / (count_WTEF1aEZH1_1+count_WTEF1aEZH1_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_DiffBindMG1655TMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WTEF1aEZH1 vs KOEF1aEZH1") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_DiffBindMG1655TMM/log2FC_qval30.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("PSC_WTEF1aEZH1 vs KOEF1aEZH1") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_DiffBindMG1655TMM/THOR_qval50.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 30) %>%
  group_by(X6) %>%
  summarise(n = n())



```

- *NOTE: FC positive = less in KO; negative = more in KO*

**Optimal qvalue:**
--> *PSC_WTEF1aEZH1 vs KOEF1aEZH1*; qval 20-30 look optimal for TMM and DiffBindMG1655TMM





--> In agreement with `003__CutRun`; in KO overall same number of gain and lost regions; and in KOEF1aEZH1 much more gain of H3K27me3 (act like the HET)

Isolate positive and negative THOR peaks to display deepTool plots

```bash
# positive negative peaks
## qval 10
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10_negative.bed

## qval 15
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval15.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval15_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval15.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval15_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval15.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval15_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval15.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval15_negative.bed

## qval 20
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval20.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval20_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval20.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval20_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval20.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval20_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval20.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval20_negative.bed

## qval 30
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval30.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval30_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval30.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval30_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval30.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval30_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval30.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval30_negative.bed

## qval 40
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval40.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval40_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval40.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval40_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval40.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval40_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval40.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval40_negative.bed

## qval 50
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval50.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval50_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval50.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval50_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval50.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval50_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval50.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval50_negative.bed

## qval 60
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval60.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval60_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval60.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval60_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval60.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval60_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval60.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval60_negative.bed

## qval 70
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval70.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval70_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval70.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval70_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval70.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval70_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval70.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval70_negative.bed

## qval 80
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval80.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval80_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval80.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval80_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval80.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval80_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval80.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval80_negative.bed

## qval 90
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval90.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval90_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval90.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval90_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval90.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval90_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval90.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval90_negative.bed

## qval 100
awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval100.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval100_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval100.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval100_negative.bed

awk -F'\t' '$16 > 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval100.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval100_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval100.bed > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval100_negative.bed
```













