# Project

**H9 cell lines**

- PSC native:
    - KO: IGG (R1, R2, R3), SUZ12 (R1, R2, R3), H3K27me3 (R2), EZH2 (R1, R2), EZH1 (R2)
    - KOEF1aEZH1: IGG, EZH2
    - WT: IGG, EZH1cs, EZH2, SUZ12





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


```bash
module load lftp

cd input_raw/

# Following email instructions
lftp -c 'set sftp:auto-confirm yes;set net:max-retries 20;open sftp://X202SC25052722-Z01-F001:2mfm2na2@usftp23.novogene.com; mirror --verbose --use-pget-n=8 -c' 



XXXY HERE !!!!!!!!!!!!!!!!!!



# Copy all .fz.gz data into input_raw/ folder
rsync -av --include '*/' --include '*.fq.gz' --exclude '*' usftp21.novogene.com/01.RawData/ input_raw/ # copy from usftp21 folder to input_raw
find input_raw/ -mindepth 2 -type f -exec mv -t input_raw/ {} + # mv files from their folder to input_raw/ folder
find input_raw/ -type d -empty -delete # delete empty directory
```

--> All good, files created in `usftp21.novogene.com/`




# Rename file

I created a tab separated file with current (`sample_name.txt`) / new file names (keeping the .fq.gz sufix) with `nano rename_map.txt`


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
sbatch scripts/fastp.sh # 28205563 ok
```


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:28205563 scripts/bowtie2.sh # 28205600 ok
```

--> Looks good; overall ~75% uniquely aligned reads


Mapping on E coli 

```bash
conda activate bowtie2

sbatch scripts/bowtie2_MG1655.sh # 29756394 xxx
```

--> Between 1 - 5% uniquely aligned reads (not a lot..; previously `005__CutRun` 10% (in `003__CutRun` was less than 1%) )





## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-28205600.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_28205600.txt
```

Add these values to `/home/roulet/001_EZH1_project/014__CutRun_PSC_native/samples_001014.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >75% input reads as been uniquely mapped to the genome (90% non uniq) 



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.




```bash
conda activate bowtie2

sbatch --dependency=afterany:28205600 scripts/samtools_unique.sh # 28205653 ok
```



Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch scripts/samtools_MG1655_unique.sh # 29772397 xxx
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

sbatch --dependency=afterany:28205653 scripts/bamtobigwig_unique.sh # 28205735 ok


```

- PSC
*Pass*: PSC_WT_SUZ12_R1 PSC_WT_EZH2_R1 PSC_KOEF1aEZH1_EZH2_R1 PSC_KO_H3K27me3_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KO_SUZ12_R3 PSC_KO_EZH2_R1 PSC_KO_EZH1_R2 PSC_KO_EZH2_R2
*Failed*: PSC_WT_EZH1_R1 (maybe very low signal, unlikely)






## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all.sh # 28744326 ok
sbatch scripts/multiBigwigSummary_WT.sh # 28744549 ok
sbatch scripts/multiBigwigSummary_KO.sh # 28744586 ok


# Plot all
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_WT_SUZ12_R1 PSC_WT_EZH2_R1 PSC_WT_EZH1_R1 PSC_WT_IGG_R1 PSC_KO_H3K27me3_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KO_SUZ12_R3 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_EZH1_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_IGG_R3 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_IGG_R1 \
    --colors black black PSC_WT_EZH1_R1 black red red red red red red red red red red blue blue \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_WT_SUZ12_R1 PSC_WT_EZH2_R1 PSC_WT_EZH1_R1 PSC_WT_IGG_R1 PSC_KO_H3K27me3_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KO_SUZ12_R3 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_EZH1_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_IGG_R3 PSC_KOEF1aEZH1_EZH2_R1 PSC_KOEF1aEZH1_IGG_R1 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf



# Plot WT
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_WT.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_WT_SUZ12_R1 PSC_WT_EZH2_R1 PSC_WT_EZH1_R1 PSC_WT_IGG_R1 \
    -o output/bigwig/multiBigwigSummary_WT_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_WT.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_WT_SUZ12_R1 PSC_WT_EZH2_R1 PSC_WT_EZH1_R1 PSC_WT_IGG_R1 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_WT_heatmap.pdf



# Plot KO
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_KO.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_KO_H3K27me3_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KO_SUZ12_R3 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_EZH1_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_IGG_R3 \
    -o output/bigwig/multiBigwigSummary_KO_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_KO.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_KO_H3K27me3_R2 PSC_KO_SUZ12_R1 PSC_KO_SUZ12_R2 PSC_KO_SUZ12_R3 PSC_KO_EZH2_R1 PSC_KO_EZH2_R2 PSC_KO_EZH1_R2 PSC_KO_IGG_R1 PSC_KO_IGG_R2 PSC_KO_IGG_R3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_KO_heatmap.pdf



```

--> H3K27me3 which works, form a group appart.. All the other cluster together. 

--> As `010__CutRun_PSC_50dN_native` H3K27me1 form group apart... Which may indicate they barely kind of work but with a completely useless signal... Not sure what to conclude... But look a bit more different than a completely failed sample that is similar to IGG.




# MACS2 peak calling on bam unique


--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad` and `narrow` **


```bash
conda activate macs2
# genotype per genotype
sbatch --dependency=afterany:17776169 scripts/macs2_broad.sh # 17795289 ok

# genotype per genotype
sbatch scripts/macs2_narrow.sh # 17867455 ok

```

--> All fail, except *H3K27me3; barely with 5,324 peaks*

**broad**:
- 50dN_WT_EZH1; n(peaks)= 2
- 50dN_WT_EZH2; n(peaks)= 1
- 50dN_WT_H3K27ac; n(peaks)= 3
- 50dN_WT_H3K27me1AM; n(peaks)= 2
- 50dN_WT_H3K27me1OR; n(peaks)= 2
- 50dN_WT_H3K27me3; n(peaks)= 5,324
- 50dN_WT_SUZ12= n(peaks)= 0

**narrow**:
- 50dN_WT_EZH1; n(peaks)= 2
- 50dN_WT_EZH2; n(peaks)= 1
- 50dN_WT_H3K27ac; n(peaks)= 3
- 50dN_WT_H3K27me1AM; n(peaks)= 2
- 50dN_WT_H3K27me1OR; n(peaks)= 1
- 50dN_WT_H3K27me3; n(peaks)= 3,830
- 50dN_WT_SUZ12= n(peaks)= 0

--> *narrow* -mode does not help....





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




