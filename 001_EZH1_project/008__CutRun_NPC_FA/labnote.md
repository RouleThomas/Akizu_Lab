# Project
- NPC FA:
    - WT
    - KO
    - KOEF1aEZH1
- 50dN FA and native:
    - WT
    - KO

--> All in simplicate

**Objectives:**
- 50dN: check whether FA improve CutRun, notably for EZH1cs (CutRun 007 failed, technical issue)
- NPC FA: check whether we detect EZH1cs in NPC in WT or EZH1cs + test other AB




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

Go [there](http://data-deliver.novogene.com/login/X202SC24012448-Z01-F001) and enter credetnial: (check email Novogen)

I created a `nano url.txt` with all link and used `wget -i url.txt` to download them all (1 link per raw); then `mv input_raw_Novogene/*fq.gz input` .


# Combine files from multiple lanes

Concatenate fastq discuss [here](https://www.biostars.org/p/317385/): `cat string_L001_sampleID_R1.fastq.gz string_L002_sampleID_R1.fastq.gz  > string_sampleID_R1.fastq.gz`.

--> Several of our samples dispatched in two lanes (`L1` and `L4`; **the concatenated are names `L14`; all unique are `L4`**)
----> input sample: `input_raw_Novogene/` output to `input/`

```bash
sbatch scripts/concatenate_1.sh # 12228326 ok 
sbatch scripts/concatenate_2.sh # 12228327 ok 
sbatch scripts/concatenate_3.sh # 12228328 ok
```



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
sbatch scripts/fastp_1.sh # 12504999 ok
sbatch scripts/fastp_2.sh # 12505000 xxx
sbatch scripts/fastp_3.sh # 12505001 xxx
```

# FastQC

**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:12504999 scripts/fastqc_fastp_1.sh # 12505357 xxx
sbatch --dependency=afterany:12505000 scripts/fastqc_fastp_2.sh # 12505358 xxx
sbatch --dependency=afterany:12505001 scripts/fastqc_fastp_3.sh # 12505359 xxx
```

--> all good


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:12505357 scripts/bowtie2_1.sh # 12505521 ok
sbatch --dependency=afterany:12505358 scripts/bowtie2_2.sh # 12505546 ok
sbatch --dependency=afterany:12505359 scripts/bowtie2_3.sh # 12505547 ok

```

--> Looks good; overall ~70% uniquely aligned reads

Mapping on E coli --> TO DO LATER!  XXX

```bash
conda activate bowtie2

sbatch scripts/bowtie2_MG1655_1.sh # xxx
sbatch scripts/bowtie2_MG1655_2.sh # xxx
sbatch scripts/bowtie2_MG1655_3.sh # xxx

```

--> between 0.5 - 2% uniquely aligned reads (not a lot..; previously `005__CutRun` 10% (in `003__CutRun` was less than 1%) )


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-12505521.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_12505521.txt

for file in slurm-12505546.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_12505546.txt

for file in slurm-12505547.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_12505547.txt

```

Add these values to `/home/roulet/001_EZH1_project/008__CutRun_NPC_FA/samples_008.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >75% input reads as been uniquely mapped to the genome (90% non uniq)





## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.


```bash
conda activate bowtie2

sbatch --dependency=afterany:12505521 scripts/samtools_unique_1.sh # 12505895 ok
sbatch --dependency=afterany:12505546 scripts/samtools_unique_2.sh # 12505896 ok
sbatch --dependency=afterany:12505547 scripts/samtools_unique_3.sh # 12505897 ok

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

sbatch --dependency=afterany:12505895 scripts/bamtobigwig_unique_1.sh # 12514789 ok
sbatch --dependency=afterany:12505896 scripts/bamtobigwig_unique_2.sh # 12514790 ok
sbatch --dependency=afterany:12505897 scripts/bamtobigwig_unique_3.sh # 12514791 ok
```



- 50dN native and FA
*Pass*: 
*Failed*: H3K27me3, EZH1cs, EZH2
- NPC _ WT
*Pass*: H3K27ac, H3K4me3
*Failed*: EZH1cs, EZH2, SUZ12
- NPC _ KOEF1aEZH1
*Pass*: H3K27me3, H3K4me3, H
*Failed*: EZH1cs, SUZ12



## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all.sh # 12654374 ok

sbatch scripts/multiBigwigSummary_50dN.sh # 12654211 ok
sbatch scripts/multiBigwigSummary_NPC.sh # 12654309 ok

sbatch scripts/multiBigwigSummary_WT.sh # 12654318 ok
sbatch scripts/multiBigwigSummary_KOEF1aEZH1.sh # 12654324; 12654560 ok
sbatch scripts/multiBigwigSummary_KO.sh # 12654331; 12654558 ok



# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels 50dNFA_KOEF1aEZH1_EZH1cs 50dNnative_KOEF1aEZH1_EZH1cs 50dNFA_KOEF1aEZH1_EZH2 50dNnative_KOEF1aEZH1_EZH2 50dNFA_KOEF1aEZH1_H3K27me3 NPC_KO_EZH1cs NPC_KO_EZH2 NPC_KO_H3K27ac NPC_KO_IGG NPC_KOEF1aEZH1_EZH1cs NPC_KOEF1aEZH1_H3K27me3 NPC_KOEF1aEZH1_H3K4me3 NPC_KOEF1aEZH1_IGG NPC_KOEF1aEZH1_SUZ12 NPC_WT_EZH1cs NPC_WT_EZH2 NPC_WT_H3K27ac NPC_WT_H3K4me3 NPC_WT_SUZ12 \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

plotPCA -in output/bigwig/multiBigwigSummary_50dN.npz \
    --transpose \
    --ntop 0 \
    --labels 50dNFA_KOEF1aEZH1_EZH1cs 50dNnative_KOEF1aEZH1_EZH1cs 50dNFA_KOEF1aEZH1_EZH2 50dNnative_KOEF1aEZH1_EZH2 50dNFA_KOEF1aEZH1_H3K27me3 \
    -o output/bigwig/multiBigwigSummary_50dN_plotPCA.pdf

plotPCA -in output/bigwig/multiBigwigSummary_NPC.npz \
    --transpose \
    --ntop 0 \
    --labels NPC_KO_EZH1cs NPC_KO_EZH2 NPC_KO_H3K27ac NPC_KO_IGG NPC_KOEF1aEZH1_EZH1cs NPC_KOEF1aEZH1_H3K27me3 NPC_KOEF1aEZH1_H3K4me3 NPC_KOEF1aEZH1_IGG NPC_KOEF1aEZH1_SUZ12 NPC_WT_EZH1cs NPC_WT_EZH2 NPC_WT_H3K27ac NPC_WT_H3K4me3 NPC_WT_SUZ12 \
    -o output/bigwig/multiBigwigSummary_NPC_plotPCA.pdf


plotPCA -in output/bigwig/multiBigwigSummary_WT.npz \
    --transpose \
    --ntop 0 \
    --labels NPC_WT_EZH1cs NPC_WT_EZH2 NPC_WT_H3K27ac NPC_WT_H3K4me3 NPC_WT_SUZ12 \
    -o output/bigwig/multiBigwigSummary_WT_plotPCA.pdf



plotPCA -in output/bigwig/multiBigwigSummary_KOEF1aEZH1.npz \
    --transpose \
    --ntop 0 \
    --labels NPC_KOEF1aEZH1_EZH1cs NPC_KOEF1aEZH1_H3K27me3 NPC_KOEF1aEZH1_H3K4me3 NPC_KOEF1aEZH1_IGG NPC_KOEF1aEZH1_SUZ12 \
    -o output/bigwig/multiBigwigSummary_KOEF1aEZH1_plotPCA.pdf



plotPCA -in output/bigwig/multiBigwigSummary_KO.npz \
    --transpose \
    --ntop 0 \
    --labels NPC_KO_EZH1cs NPC_KO_EZH2 NPC_KO_H3K27ac NPC_KO_IGG \
    -o output/bigwig/multiBigwigSummary_KO_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels 50dNFA_KOEF1aEZH1_EZH1cs 50dNnative_KOEF1aEZH1_EZH1cs 50dNFA_KOEF1aEZH1_EZH2 50dNnative_KOEF1aEZH1_EZH2 50dNFA_KOEF1aEZH1_H3K27me3 NPC_KO_EZH1cs NPC_KO_EZH2 NPC_KO_H3K27ac NPC_KO_IGG NPC_KOEF1aEZH1_EZH1cs NPC_KOEF1aEZH1_H3K27me3 NPC_KOEF1aEZH1_H3K4me3 NPC_KOEF1aEZH1_IGG NPC_KOEF1aEZH1_SUZ12 NPC_WT_EZH1cs NPC_WT_EZH2 NPC_WT_H3K27ac NPC_WT_H3K4me3 NPC_WT_SUZ12 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf

plotCorrelation \
    -in output/bigwig/multiBigwigSummary_50dN.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels 50dNFA_KOEF1aEZH1_EZH1cs 50dNnative_KOEF1aEZH1_EZH1cs 50dNFA_KOEF1aEZH1_EZH2 50dNnative_KOEF1aEZH1_EZH2 50dNFA_KOEF1aEZH1_H3K27me3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_50dN_heatmap.pdf

plotCorrelation \
    -in output/bigwig/multiBigwigSummary_NPC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels NPC_KO_EZH1cs NPC_KO_EZH2 NPC_KO_H3K27ac NPC_KO_IGG NPC_KOEF1aEZH1_EZH1cs NPC_KOEF1aEZH1_H3K27me3 NPC_KOEF1aEZH1_H3K4me3 NPC_KOEF1aEZH1_IGG NPC_KOEF1aEZH1_SUZ12 NPC_WT_EZH1cs NPC_WT_EZH2 NPC_WT_H3K27ac NPC_WT_H3K4me3 NPC_WT_SUZ12 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_NPC_heatmap.pdf

plotCorrelation \
    -in output/bigwig/multiBigwigSummary_WT.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels NPC_WT_EZH1cs NPC_WT_EZH2 NPC_WT_H3K27ac NPC_WT_H3K4me3 NPC_WT_SUZ12 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_WT_heatmap.pdf

plotCorrelation \
    -in output/bigwig/multiBigwigSummary_KOEF1aEZH1.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels NPC_KOEF1aEZH1_EZH1cs NPC_KOEF1aEZH1_H3K27me3 NPC_KOEF1aEZH1_H3K4me3 NPC_KOEF1aEZH1_IGG NPC_KOEF1aEZH1_SUZ12 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_KOEF1aEZH1_heatmap.pdf


plotCorrelation \
    -in output/bigwig/multiBigwigSummary_KO.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels NPC_KO_EZH1cs NPC_KO_EZH2 NPC_KO_H3K27ac NPC_KO_IGG \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_KO_heatmap.pdf

```

--> Hard to conclude stuff



# MACS2 peak calling on bam unique



--> IGG samples used as control when available; **for NPC WT, IGG from KO has been used**
----> For 50dN, no IGG control used for peak calling...


--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad_50dN.sh # 12654885 ok
sbatch scripts/macs2_broad_KO_KOEF1aEZH1.sh # 12655091 ok
sbatch scripts/macs2_broad_WT.sh # 12655165 ok

sbatch scripts/macs2_narrow_50dN.sh # 12655462 ok
sbatch scripts/macs2_narrow_KO_KOEF1aEZH1.sh # 12655472 ok
sbatch scripts/macs2_narrow_WT.sh # 12655475 ok





```

--> OEF1aEZH1 in 50-day neurons: too noisy for EZH1, EZH2, and H3K27me3

--> NPC histone marks: OK for H3K4me3, H3K27ac, and H3K27me3; sharp and clear peaks.

--> NPC PRC2 components: too noisy...



XXXXXXXXXX below not mod





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



