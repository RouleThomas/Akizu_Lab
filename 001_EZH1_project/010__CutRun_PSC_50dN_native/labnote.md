# Project


- PSC and 50dN native (FA was not so good so back to native!):
    - WT: H3K27me3, H3K27me1, EZH2, IGG

--> All in simplicate


**Objectives:**
- Some issues with previous CutRun, here is more a test with few samples, only WT, but some tricky AB to check whether CutRun work again. **Working samples can be added as aditional WT replicate!**




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

Go [there](http://data-deliver.novogene.com/login/X202SC24031197-Z01-F001) and enter credetnial: (check email Novogen)

I created a `nano url.txt` with all link and used `wget -i url.txt` to download them all (1 link per raw); then `mv input_raw_Novogene/*fq.gz input` .



# Rename file

Renamed manually as only 8 samples



```bash
cp input_raw_Novogene/*.gz input/
```

--> All good 


# Fastp cleaning

```bash
sbatch scripts/fastp_1.sh # 16520720 xxx
sbatch scripts/fastp_2.sh # 16520721 xxx

```


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:16520720 scripts/bowtie2_1.sh # 16520751 xxx
sbatch --dependency=afterany:16520721 scripts/bowtie2_2.sh # 16520752 xxx
```

--> XXX Looks good; overall ~70% uniquely aligned reads XXX

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Mapping on E coli --> TO DO LATER! 

```bash
conda activate bowtie2

sbatch scripts/bowtie2_MG1655_1.sh # 13345349 ok
sbatch scripts/bowtie2_MG1655_2.sh # 13345352 ok
sbatch scripts/bowtie2_MG1655_3.sh # 13345353 ok
sbatch scripts/bowtie2_MG1655_missing.sh # 13345354 ok

```

--> between 0.5 - 2% uniquely aligned reads (not a lot..; previously `005__CutRun` 10% (in `003__CutRun` was less than 1%) )


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-16520751.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_16520751.txt

for file in slurm-16520752.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_16520752.txt

```

Add these values to `/home/roulet/001_EZH1_project/010__CutRun_PSC_50dN_native/samples_009.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> XXX Overall >75% input reads as been uniquely mapped to the genome (90% non uniq) XXX



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.




```bash
conda activate bowtie2

sbatch --dependency=afterany:16520751 scripts/samtools_unique_1.sh # 16520802 xxx
sbatch --dependency=afterany:16520752 scripts/samtools_unique_2.sh # 16520803 xxx

```

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch --dependency=afterany:13345349:13345352:13345353:13345354 scripts/samtools_MG1655_unique_1.sh # 13345712 xxx
sbatch --dependency=afterany:13345349:13345352:13345353:13345354 scripts/samtools_MG1655_unique_2.sh # 13345713 xxx
sbatch --dependency=afterany:13345349:13345352:13345353:13345354 scripts/samtools_MG1655_unique_3.sh # 13345715 xxx
sbatch --dependency=afterany:13345349:13345352:13345353:13345354 scripts/samtools_MG1655_unique_missing.sh # 13345737 xxx


```

--> More information on this step in the `005__CutRun` labnote

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:16520802 scripts/bamtobigwig_unique_1.sh # 16520886 xxx
sbatch --dependency=afterany:16520803 scripts/bamtobigwig_unique_2.sh # 16520889 xxx


```


- PSC native
*Pass*: 
*Failed*: 
- 50dN native
*Pass*: 
*Failed*:




--> Non unique (all raw reads!) vs unique bigwig (less signal or more noise?): Very similar. increase signal on the one that work but more bakground


## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch --dependency=afterany:16520886:16520889 scripts/multiBigwigSummary_all.sh # 16521108 xxx

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels 50dNFA_KOEF1aEZH1_EZH1cs 50dNnative_KOEF1aEZH1_EZH1cs 50dNFA_KOEF1aEZH1_EZH2 50dNnative_KOEF1aEZH1_EZH2 50dNFA_KOEF1aEZH1_H3K27me3 NPC_KO_EZH1cs NPC_KO_EZH2 NPC_KO_H3K27ac NPC_KO_IGG NPC_KOEF1aEZH1_EZH1cs NPC_KOEF1aEZH1_H3K27me3 NPC_KOEF1aEZH1_H3K4me3 NPC_KOEF1aEZH1_IGG NPC_KOEF1aEZH1_SUZ12 NPC_WT_EZH1cs NPC_WT_EZH2 NPC_WT_H3K27ac NPC_WT_H3K4me3 NPC_WT_SUZ12 \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels 50dNFA_KOEF1aEZH1_EZH1cs 50dNnative_KOEF1aEZH1_EZH1cs 50dNFA_KOEF1aEZH1_EZH2 50dNnative_KOEF1aEZH1_EZH2 50dNFA_KOEF1aEZH1_H3K27me3 NPC_KO_EZH1cs NPC_KO_EZH2 NPC_KO_H3K27ac NPC_KO_IGG NPC_KOEF1aEZH1_EZH1cs NPC_KOEF1aEZH1_H3K27me3 NPC_KOEF1aEZH1_H3K4me3 NPC_KOEF1aEZH1_IGG NPC_KOEF1aEZH1_SUZ12 NPC_WT_EZH1cs NPC_WT_EZH2 NPC_WT_H3K27ac NPC_WT_H3K4me3 NPC_WT_SUZ12 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


```

--> XXX Hard to conclude stuff XXX



# MACS2 peak calling on bam unique



--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch --dependency=afterany:16520802 scripts/macs2_broad_1.sh # 16521538 xxx
sbatch --dependency=afterany:16520803 scripts/macs2_broad_2.sh # 16521568 xxx

```



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