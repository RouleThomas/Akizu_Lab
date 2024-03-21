# Project

**H9 cell lines**

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
sbatch scripts/fastp_1.sh # 16520720 ok
sbatch scripts/fastp_2.sh # 16520721 ok

```


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:16520720 scripts/bowtie2_1.sh # 16520751 ok
sbatch --dependency=afterany:16520721 scripts/bowtie2_2.sh # 16520752 ok
```

--> Looks good; overall ~75% uniquely aligned reads

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

Add these values to `/home/roulet/001_EZH1_project/010__CutRun_PSC_50dN_native/samples_010.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >75% input reads as been uniquely mapped to the genome (90% non uniq) 



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.




```bash
conda activate bowtie2

sbatch --dependency=afterany:16520751 scripts/samtools_unique_1.sh # 16520802 ok
sbatch --dependency=afterany:16520752 scripts/samtools_unique_2.sh # 16520803 ok

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

sbatch --dependency=afterany:16520802 scripts/bamtobigwig_unique_1.sh # 16520886 ok
sbatch --dependency=afterany:16520803 scripts/bamtobigwig_unique_2.sh # 16520889 ok


```


- PSC native
*Pass*: PSC_H3K27me3, PSC_EZH2
*Failed*: PSC_H3K27me1
- 50dN native
*Pass*: NA
*Failed*: 50dN_H3K27me3, 50dN_EZH2, 50dN_H3K27me1





## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch --dependency=afterany:16520886:16520889 scripts/multiBigwigSummary_all.sh # 16521108 ok



# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_WT_EZH2 PSC_WT_H3K27me1 PSC_WT_H3K27me3 PSC_WT_IGG 50dN_WT_EZH2 50dN_WT_H3K27me1 50dN_WT_H3K27me3 50dN_WT_IGG \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_WT_EZH2 PSC_WT_H3K27me1 PSC_WT_H3K27me3 PSC_WT_IGG 50dN_WT_EZH2 50dN_WT_H3K27me1 50dN_WT_H3K27me3 50dN_WT_IGG \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf


```

--> 50dN form a group, with IGG PSC; confirming fail

--> H3K27me1 is apart, which show it may have work. Compare with ENCODE data, no ENCODE data.. Nor available data!! I only found [cancer cell lines HCT116 with H3K27me1 ChIP](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75217)





# MACS2 peak calling on bam unique



--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad_1.sh # 16581946 ok
sbatch scripts/macs2_broad_2.sh # 16581955 ok

sbatch scripts/macs2_narrow_1.sh # 16582595 ok

```



--> 50dN no to few peaks

--> PSC H3K27me1, no peaks

--> PSC EZH2 ~ 750 peaks broad and 543 narrow (3,580 in 006__CutRun...)





XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX below not mod





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

XXX


# deepTool plots


```bash
conda activate deeptools

# all genes 
sbatch scripts/matrix_TSS_10kb_allGenes.sh # 16677909 xxx
sbatch scripts/matrix_TSS_10kb_PSC_allGenes.sh # 16678016 xxx

# H3K27me3 peaks
sbatch scripts/matrix_TSS_10kb_PSC_H3K27me3Peaks.sh # 16678626 ok

```


--> PSC; H3K27me3 and EZH2 co-localize well


