# Project

2 conditions, 2 genotypes, plenty of AB
- NPC; WT and KO; AB: EZH1 (2 AB cs and pt), EZH2, H3K27me1/3, SUZ12 and IGG --> Objective, check whether SUZ12 overlap with EZH1 solely? Or also EZH2? Does EZH2 expression/binding increase with EZH1 KO?
- PSC (ESC); KOsynEZH1 (KO EZH1 with synapse-specific EZH1-HA tag=negative control) and KOEF1aEZH1 (KO EZH1 with ef1a strong promoter EZH1-HA); AB: EZH1 (2 AB cs and pt), EZH2, H3K27me3, SUZ12 and IGG


**Objectives:**
- Check EZH1, EZH2, H3K27me3 overlapping and notably upon EZH1 KO
- Validate the EZH1 AB specificity (compare EZH1-tag with EZH1 AB)


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

Go [there](http://data-deliver.novogene.com/login/X202SC23092052-Z01-F001) and enter credetnial: (check email Novogen)

I created a `nano url.txt` with all link and used `wget -i url.txt` to download them all (1 link per raw); then `mv input_raw_Novogene/*fq.gz input` .

# Rename file

I created a tab separated file with current / new file names (keeping the .fq.gz sufix) and then:

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_map.txt
```

--> All good

# Fastp cleaning

```bash
sbatch scripts/fastp_1.sh # 6677035 ok
sbatch scripts/fastp_2.sh # 6677038 ok
sbatch scripts/fastp_3.sh # 6677039 ok
```


# FastQC

**Raw:**
```bash
sbatch scripts/fastqc_1.sh # 6677168 ok
sbatch scripts/fastqc_2.sh # 6677169 ok
sbatch scripts/fastqc_3.sh # 6677171 ok
```

--> all good

**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:6677035 scripts/fastqc_fastp_1.sh # 6677362 ok
sbatch --dependency=afterany:6677038 scripts/fastqc_fastp_2.sh # 6677368 ok
sbatch --dependency=afterany:6677039 scripts/fastqc_fastp_3.sh # 6677369 ok
```

--> all good


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch scripts/bowtie2_1.sh # 6678343 ok
sbatch scripts/bowtie2_2.sh # 6678369 ok
sbatch scripts/bowtie2_3.sh # 6678383 ok
```

--> Looks good


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-6678343.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_6678343.txt

for file in slurm-6678369.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_6678369.txt

for file in slurm-6678383.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_6678383.txt
```

Add these values to `/home/roulet/001_EZH1_project/006__CutRun_PSC_FA/samples_006.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >60% input reads as been uniquely mapped to the genome (75% non uniq)



# Calculate histone content
**XX Not done; not necessary as only H3K27me3... XX**


# Samtools and read filtering

--> See `METHOD GOOD TO FOLLOW` in `003__CutRun` labnote

## Marking dupplicates
```bash
conda activate bowtie2

sbatch scripts/samtools_1.sh # 6712195 ok
sbatch scripts/samtools_2.sh # 6712214 ok
sbatch scripts/samtools_3.sh # 6712231 ok
```

## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.

```bash
conda activate bowtie2

sbatch --dependency=afterany:6712195 scripts/samtools_unique_1.sh # 6712296 ok
sbatch --dependency=afterany:6712214 scripts/samtools_unique_2.sh # 6712305 ok
sbatch --dependency=afterany:6712231 scripts/samtools_unique_3.sh # 6712321 ok
```

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch scripts/samtools_MG1655_unique_1.sh # xxx
sbatch scripts/samtools_MG1655_unique_2.sh # xxx
sbatch scripts/samtools_MG1655_unique_3.sh # xxx

```

--> More information on this step in the `005__CutRun` labnote

# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:6712296 scripts/bamtobigwig_unique_1.sh # 6712485 ok
sbatch --dependency=afterany:6712305 scripts/bamtobigwig_unique_2.sh # 6712486 ok
sbatch --dependency=afterany:6712321 scripts/bamtobigwig_unique_3.sh # 6712487 ok

sbatch scripts/bamtobigwig_unique_1_missing.sh # 6731623
```

--> XXX 




## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch --dependency=afterany:6731623 scripts/multiBigwigSummary_all.sh # 6731623



XXXXXX modify below :


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_PSC.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KOEF1aEZH1_EZH1cs PSC_KOsynEZH1_EZH1cs PSC_KOEF1aEZH1_EZH1pt PSC_KOsynEZH1_EZH1pt PSC_KOEF1aEZH1_H3K27me3 PSC_KOsynEZH1_H3K27me3 PSC_KOEF1aEZH1_HA PSC_KOsynEZH1_HA PSC_KOEF1aEZH1_SUZ12 PSC_KOsynEZH1_SUZ12 PSC_KOEF1aEZH1_IGG PSC_KOsynEZH1_IGG \
    -o output/bigwig/multiBigwigSummary_PSC_plotPCA.pdf
plotPCA -in output/bigwig/multiBigwigSummary_NPC.npz \
    --transpose \
    --ntop 0 \
    --labels NPC_WT_EZH1cs NPC_KO_EZH1cs NPC_WT_EZH1pt NPC_KO_EZH1pt NPC_WT_H3K27me1 NPC_KO_H3K27me1 NPC_WT_H3K27me3 NPC_KO_H3K27me3 NPC_WT_H3K4me3 NPC_KO_H3K4me3 NPC_WT_EZH2 NPC_KO_EZH2 NPC_WT_SUZ12 NPC_KO_SUZ12 NPC_WT_IGG NPC_KO_IGG \
    -o output/bigwig/multiBigwigSummary_NPC_plotPCA.pdf
plotPCA -in output/THOR/multiBigwigSummary_NPC_THOR.npz \
    --transpose \
    --ntop 0 \
    --labels WT_SUZ12 KO_SUZ12 WT_EZH2 KO_EZH2 WT_H3K27me3 KO_H3K27me3 WT_H3K4me3 KO_H3K4me3 \
    --colors blue blue green green red red gold gold \
    --markers o s o s o s o s \
    --plotWidth 7 \
    -o output/THOR/multiBigwigSummary_NPC_THOR_plotPCA.pdf


## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_PSC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KOEF1aEZH1_EZH1cs PSC_KOsynEZH1_EZH1cs PSC_KOEF1aEZH1_EZH1pt PSC_KOsynEZH1_EZH1pt PSC_KOEF1aEZH1_H3K27me3 PSC_KOsynEZH1_H3K27me3 PSC_KOEF1aEZH1_HA PSC_KOsynEZH1_HA PSC_KOEF1aEZH1_SUZ12 PSC_KOsynEZH1_SUZ12 PSC_KOEF1aEZH1_IGG PSC_KOsynEZH1_IGG \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_PSC_heatmap.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_PSC_subset_1.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KOEF1aEZH1_EZH1cs PSC_KOsynEZH1_EZH1cs   PSC_KOEF1aEZH1_H3K27me3 PSC_KOsynEZH1_H3K27me3 PSC_KOEF1aEZH1_SUZ12 PSC_KOsynEZH1_SUZ12 PSC_KOEF1aEZH1_IGG PSC_KOsynEZH1_IGG \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_PSC_subset_1_heatmap.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_NPC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels NPC_WT_EZH1cs NPC_KO_EZH1cs NPC_WT_EZH1pt NPC_KO_EZH1pt NPC_WT_H3K27me1 NPC_KO_H3K27me1 NPC_WT_H3K27me3 NPC_KO_H3K27me3 NPC_WT_H3K4me3 NPC_KO_H3K4me3 NPC_WT_EZH2 NPC_KO_EZH2 NPC_WT_SUZ12 NPC_KO_SUZ12 NPC_WT_IGG NPC_KO_IGG \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_NPC_heatmap.pdf
plotCorrelation \
    -in output/THOR/multiBigwigSummary_NPC_THOR.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_SUZ12 KO_SUZ12 WT_EZH2 KO_EZH2 WT_H3K27me3 KO_H3K27me3 WT_H3K4me3 KO_H3K4me3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/THOR/multiBigwigSummary_NPC_THOR_heatmap.pdf
```












