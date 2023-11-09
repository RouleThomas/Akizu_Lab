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

sbatch scripts/bamtobigwig_unique_1_missing.sh # 6731623 ok
```



- KOEF1aEZH1
*Pass*: SUZ12, EZH1cs, EZH2, H3K27me3
*Failed*: HA
- WT
*Pass*: SUZ12, EZH2, H3K27me3
*Failed*: HA, EZH1cs (better than earlier, more peaks, but still poorly)
- KO
*Pass*: H3K27me3
*Failed*: HA, SUZ12, EZH1cs, EZH2





## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch --dependency=afterany:6731623 scripts/multiBigwigSummary_all.sh # 6731623


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KOEF1aEZH1_SUZ12 PSC_KOEF1aEZH1_EZH2 PSC_KOEF1aEZH1_HA PSC_KOEF1aEZH1_EZH1cs PSC_KOEF1aEZH1_H3K27me3 PSC_KOEF1aEZH1_IGG PSC_WT_SUZ12 PSC_WT_EZH2 PSC_WT_HA PSC_WT_EZH1cs PSC_WT_H3K27me3 PSC_WT_IGG PSC_KO_SUZ12 PSC_KO_EZH2 PSC_KO_HA PSC_KO_EZH1cs PSC_KO_H3K27me3 PSC_KO_IGG \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf


## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KOEF1aEZH1_SUZ12 PSC_KOEF1aEZH1_EZH2 PSC_KOEF1aEZH1_HA PSC_KOEF1aEZH1_EZH1cs PSC_KOEF1aEZH1_H3K27me3 PSC_KOEF1aEZH1_IGG PSC_WT_SUZ12 PSC_WT_EZH2 PSC_WT_HA PSC_WT_EZH1cs PSC_WT_H3K27me3 PSC_WT_IGG PSC_KO_SUZ12 PSC_KO_EZH2 PSC_KO_HA PSC_KO_EZH1cs PSC_KO_H3K27me3 PSC_KO_IGG \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf


```






# MACS2 peak calling on bam unique

--> IGG samples used as control

--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_raw_1.sh # 6749578
sbatch scripts/macs2_raw_2.sh # 6749582
sbatch scripts/macs2_raw_3.sh # 6749583
```


XXX HBELOW TO MOD XXX


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
- NPC_H3K27me3 WT and KO; 2.30103
- NPC_SUZ12 WT and KO; 1.30103
- NPC_EZH2 WT and KO; 1.30103
- NPC_H3K4me3 WT and KO; 1.30103 
- PSC_EZH1cs KOEF1aEZH1; 1.30103 
- PSC_SUZ12 EF1aEZH1 and synEZH1; 1.30103 
- PSC_H3K27me3 EF1aEZH1 and synEZH1; 2.30103


--> No peak called for H3K27me1 samples







