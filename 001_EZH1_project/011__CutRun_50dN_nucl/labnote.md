# Project

**H9 cell lines**

- 50dN native (FA was not so good so back to native!):
    - WT: H3K27me3, H3K27me1 (2 AB tested), H3K2ac, EZH2, , EZH1, IGG

--> All in simplicate


**Objectives:**
- Some issues with previous CutRun, here is more a test with few samples, only WT.

This time, *no nuclear purification has been performed.*

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

Go [there](http://data-deliver.novogene.com/batchfiles/X202SC24035280-Z01-F001) and enter credetnial: (check email Novogen)

I created a `nano url.txt` with all link and used `wget -i url.txt` to download them all (1 link per raw); then `mv input_raw_Novogene/*fq.gz input` .



# Rename file

Renamed manually as only 8 samples



```bash
cp input_raw_Novogene/*.gz input/
```

--> All good 



# Fastp cleaning

```bash
sbatch scripts/fastp.sh # 17775901 xxx
```


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:17775901 scripts/bowtie2.sh # 17775970 xxx
```

--> XXX Looks good; overall ~75% uniquely aligned reads XXX

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Mapping on E coli --> TO DO LATER! 

```bash
conda activate bowtie2

sbatch scripts/bowtie2_MG1655_1.sh # 13345349 ok
```

--> between 0.5 - 2% uniquely aligned reads (not a lot..; previously `005__CutRun` 10% (in `003__CutRun` was less than 1%) )


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-17775901.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_17775901.txt
```

Add these values to `/home/roulet/001_EZH1_project/010__CutRun_PSC_50dN_native/samples_010.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >75% input reads as been uniquely mapped to the genome (90% non uniq) 



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.




```bash
conda activate bowtie2

sbatch --dependency=afterany:17775970 scripts/samtools_unique.sh # 17776169 xxx
```

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch --dependency=afterany:13345349:scripts/samtools_MG1655_unique_1.sh # 13345712 xxx
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

sbatch --dependency=afterany:17776169 scripts/bamtobigwig_unique.sh # 17776280 xxx


```

XXXXXXXXXXXXXXXXXXXXX below not mnod XXXXXXXXXXXXXXXXXXXXX



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

XXX