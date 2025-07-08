# Project
- NPC FA:
    - WT
    - KO
    - KOEF1aEZH1
- 50dN FA and native:
    - WT
    - KO

--> All in simplicate

**!!! KOEF1aEZH1 do NOT overexpress EZH1 !!! So should be considered as KO !!!**

**Objectives:**
- 50dN: check whether FA improve CutRun, notably for EZH1cs (CutRun 007 failed, technical issue)
- NPC FA: check whether we detect EZH1cs in NPC in WT or EZH1cs + test other AB

# CONCLUSION samples QC

--> 50dN did NOT WORK!

--> NPC: H3K27me3 WT (1 Bio rep), H3K27me3 KO (1 Bio rep)


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


Rename the missing ones! The one that I did not concatenate!!:
`cp input_raw_Novogene/NPC_WT_IGG_CKDL240001568-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NPC_WT_IGG_CKDL240001568-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NPC_WT_H3K27me3_CKDL240001570-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NPC_WT_H3K27me3_CKDL240001570-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NPC_KO_H3K4me3_CKDL240001576-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NPC_KO_H3K4me3_CKDL240001576-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NPC_KO_H3K27me3_CKDL240001577-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NPC_KO_H3K27me3_CKDL240001577-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NPC_KO_SUZ12_CKDL240001581-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NPC_KO_SUZ12_CKDL240001581-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NPC_OE_KOH3K27ac_CKDL240001585-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NPC_OE_KOH3K27ac_CKDL240001585-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NPC_OE_KO_EZH2_CKDL240001587-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NPC_OE_KO_EZH2_CKDL240001587-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NEU_OE_KO_IGG_CL_CKDL240001589-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NEU_OE_KO_IGG_CL_CKDL240001589-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NEU_OE_KO_IGG_CKDL240001593-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NEU_OE_KO_IGG_CKDL240001593-1A_HWK3JDSX7_L4_2.fq.gz input_raw_Novogene/NEU_OE_KO_H3K27me3_CKDL240001594-1A_HWK3JDSX7_L4_1.fq.gz input_raw_Novogene/NEU_OE_KO_H3K27me3_CKDL240001594-1A_HWK3JDSX7_L4_2.fq.gz input/ `



```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_map_missing.txt
```

--> All good 


# Fastp cleaning

```bash
sbatch scripts/fastp_1.sh # 12504999 ok
sbatch scripts/fastp_2.sh # 12505000 ok
sbatch scripts/fastp_3.sh # 12505001 ok

sbatch scripts/fastp_missing.sh # 13300220 ok

```

# FastQC

**raw**
```bash
sbatch scripts/fastqc_raw.sh # 12696467 ok
```


**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:12504999 scripts/fastqc_fastp_1.sh # 12505357 ok
sbatch --dependency=afterany:12505000 scripts/fastqc_fastp_2.sh # 12505358 ok
sbatch --dependency=afterany:12505001 scripts/fastqc_fastp_3.sh # 12505359 ok
```

--> all good


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:12505357 scripts/bowtie2_1.sh # 12505521 ok
sbatch --dependency=afterany:12505358 scripts/bowtie2_2.sh # 12505546 ok
sbatch --dependency=afterany:12505359 scripts/bowtie2_3.sh # 12505547 ok

sbatch scripts/bowtie2_missing.sh # 13300398 ok
```

--> Looks good; overall ~70% uniquely aligned reads

Mapping on E coli --> TO DO LATER! 

```bash
conda activate bowtie2

sbatch scripts/bowtie2_MG1655_1.sh # 13345349 ok
sbatch scripts/bowtie2_MG1655_2.sh # 13345352 ok
sbatch scripts/bowtie2_MG1655_3.sh # 13345353 ok
sbatch scripts/bowtie2_MG1655_missing.sh # 13345354 ok

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

for file in slurm-13300398.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_13300398.txt


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

sbatch scripts/samtools_unique_missing.sh # 13343264 ok
```

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch --dependency=afterany:13345349:13345352:13345353:13345354 scripts/samtools_MG1655_unique_1.sh # 13345712 xxx
sbatch --dependency=afterany:13345349:13345352:13345353:13345354 scripts/samtools_MG1655_unique_2.sh # 13345713 xxx
sbatch --dependency=afterany:13345349:13345352:13345353:13345354 scripts/samtools_MG1655_unique_3.sh # 13345715 xxx
sbatch --dependency=afterany:13345349:13345352:13345353:13345354 scripts/samtools_MG1655_unique_missing.sh # 13345737 xxx


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

sbatch --dependency=afterany:13343264 scripts/bamtobigwig_unique_missing.sh # 13343408 xxx



# Non unique bigiwig
sbatch scripts/bamtobigwig_1.sh # 13164849 ok
sbatch scripts/bamtobigwig_2.sh # 13164850 ok
sbatch scripts/bamtobigwig_3.sh # 13164851 ok
```



- 50dN native and FA
*Failed*: H3K27me3, EZH1cs, EZH2
- NPC _ WT
*Pass*: H3K27ac, H3K4me3
*Failed*: EZH1cs, EZH2, SUZ12
- NPC _ KOEF1aEZH1
*Pass*: H3K27me3, H3K4me3, H
*Failed*: EZH1cs, SUZ12


--> Non unique (all raw reads!) vs unique bigwig (less signal or more noise?): Very similar. increase signal on the one that work but more bakground


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



--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch --dependency=afterany:13343264 scripts/macs2_broad_50dN.sh # 12654885 ok missed sample added; 13343914 ok
sbatch --dependency=afterany:13343264 scripts/macs2_broad_KO_KOEF1aEZH1.sh # 12655091 ok missed sample added; 13344031 ok
sbatch --dependency=afterany:13343264 scripts/macs2_broad_WT.sh # 12655165 ok missed sample added; 13344112 ok

sbatch --dependency=afterany:13343264 scripts/macs2_narrow_50dN.sh # 12655462 ok missed sample added; 13344314 ok
sbatch --dependency=afterany:13343264 scripts/macs2_narrow_KO_KOEF1aEZH1.sh # 12655472 ok missed sample added; 13344540 ok
sbatch --dependency=afterany:13343264 scripts/macs2_narrow_WT.sh # 12655475 ok missed sample added; 13344565 ok

```

--> OEF1aEZH1 in 50-day neurons: too noisy for EZH1, EZH2, and H3K27me3

--> NPC histone marks: OK for H3K4me3, H3K27ac, and H3K27me3; sharp and clear peaks.

--> NPC PRC2 components: too noisy...

*- NOTE: peak calling has been run 2 times adding the missing samples!*

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

Let's do the analysis for H3K27me3, H3K4me3 only; **compare WT vs KO only; as the KOEF1aEZH1 is NOT overexpressing!**

vs KOEF1a. Test 2 spikein normalization method (histone and Ecoli)


## Calculate histone content

--> This histone content will be used to generate a scaling factor which will be used to histone-scaled our library size. The calcul/method to follow is from `003__CutRun/output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt`

**Pipeline:**
- Count the histone barcode on the clean reads
- Calculate SF (group by sample (replicate) and AB and calculate the total nb of reads. Then proportion of reads = nb read in sample / total reads. SF = min(proportion) / sample proportion)


## Count the histone barcode on the clean reads



```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp.sh # 12922764 ok

sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_missing.sh # 13346147 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_missing_1.sh # 13354832 ok

```

--> It output the nb of reads found for each histone; then simply copy paste to the excell file `output/spikein/SpikeIn_QC_fastp_008.xlsx` in GoogleDrive

- `50dNFA_KOEF1aEZH1_H3K27me3`: enriched in H3K27me3, but much less nb of reads as compare to the NPC sample!
- `NPC_KO_H3K27ac`: not in histone control..
- `NPC_KOEF1aEZH1_H3K27me3`: enriched in H3K27me3
- `NPC_KOEF1aEZH1_H3K4me3`: enriched in H3K4me3
- `NPC_WT_H3K27ac`: not in histone control..
- `NPC_WT_H3K4me3`: enriched in H3K4me3



## histone spike in factor


--> SF only calculating on WT and KO as KOEF1aEZH1 is NOT overexpressing..

```R
# package
library("tidyverse")
library("readxl")
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_008.xlsx") 

## H3K27me3 with only WT and KO
spikein_H3K27me3 = spikein %>%
    filter(Target == "H3K27me3",
    sample_ID %in% c("NPC_WT_H3K27me3", "NPC_KO_H3K27me3")) %>%
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
write.table(spikein_H3K27me3_scaling_factor, file="output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt", sep="\t", quote=FALSE, row.names=FALSE)


## H3K4me3
spikein_H3K4me3 = spikein %>%
    filter(Target == "H3K4me3",
    sample_ID %in% c("NPC_WT_H3K4me3", "NPC_KO_H3K4me3")) %>%
    group_by(sample_ID, AB) %>%
    summarise(aligned=sum(counts))
# Total reads per IP
spikein_H3K4me3_total = spikein_H3K4me3 %>%
    ungroup() %>%
    group_by(AB) %>%
    mutate(total = sum(aligned)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_H3K4me3_read_prop = spikein_H3K4me3 %>%
    left_join(spikein_H3K4me3_total) %>%
    mutate(read_prop = aligned / total)
spikein_H3K4me3_read_prop_min = spikein_H3K4me3_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_H3K4me3_scaling_factor = spikein_H3K4me3_read_prop %>%
    left_join(spikein_H3K4me3_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_H3K4me3_scaling_factor, file="output/spikein/spikein_histone_H3K4me3_scaling_factor_fastp.txt", sep="\t", quote=FALSE, row.names=FALSE)


```

--> All good; in KO higher SF for H3K27me3




### Quality control plot

Then look at the xlsx file from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) to generate quality control plot. Use R cluster for vizualization (file is `spikein_QC.xlsx` in Google Drive), file in `output/spikein`.
```R
# package
library("tidyverse")
library("readxl")
# import df adn tidy to remove AB used in sample_ID
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_008.xlsx") %>%
  separate(sample_ID, into = c("type", "condition", "tag"), sep = "_") %>%
  mutate(sample_ID = paste(type, condition, sep = "_")) %>%
  select(-type, -condition, -tag, -tissue)


# NPC
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
pdf("output/spikein/QC_histone_spike_in_H3K27me3.pdf", width = 10, height = 4)
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


## Histone scaling for H3K4me3
spikein_all_scale = spikein_all %>%
  group_by(sample_ID) %>%
  # Find the target_norm value when Target is H3K4me3 and AB is H3K4me3
  mutate(scaling_factor = ifelse(Target == "H3K4me3" & AB == "H3K4me3", target_norm, NA)) %>%
  # Fill the scaling_factor column with the appropriate value within each group
  fill(scaling_factor, .direction = "downup") %>%
  # Scale the target_norm values
  mutate(scaled_target_norm = target_norm / scaling_factor * 100) %>%
  # Remove the scaling_factor column
  select(-scaling_factor) %>%
  # Ungroup the data
  ungroup()
# Plot
pdf("output/spikein/QC_histone_spike_in_H3K4me3.pdf", width = 10, height = 4)
spikein_all_scale %>%
    filter(
           AB %in% c("H3K4me3")) %>%
        ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
        geom_col(position = "dodge") +
        facet_wrap(~sample_ID, nrow=1) +
        geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


```


--> All good H3K27me3 and H3K4me3 enriched



# Ecoli scaling factor
## Mapping E coli

**Do it for H3K27me3, H3K4me3, H3K27ac for WT and KO**


- Map our reads to the E. coli genome using same parameters as for human.
- Count the number of aligned reads to the spike-in control sequences for each sample `samtools view -S -F 4 -c sample.sam > sample_spikein_count.txt`
- Do the math for scaling factor, same method as when using histone spike-in

```bash
# count nb of reads aligned to genome

samtools view -S -F 4 -c output/spikein/NPC_KO_H3K27ac_MG1655.sam > output/spikein/NPC_KO_H3K27ac-spikein_count.txt
samtools view -S -F 4 -c output/spikein/NPC_WT_H3K27ac_MG1655.sam > output/spikein/NPC_WT_H3K27ac-spikein_count.txt
samtools view -S -F 4 -c output/spikein/NPC_KO_H3K27me3_MG1655.sam > output/spikein/NPC_KO_H3K27me3-spikein_count.txt
samtools view -S -F 4 -c output/spikein/NPC_WT_H3K27me3_MG1655.sam > output/spikein/NPC_WT_H3K27me3-spikein_count.txt
samtools view -S -F 4 -c output/spikein/NPC_KO_H3K4me3_MG1655.sam > output/spikein/NPC_KO_H3K4me3-spikein_count.txt
samtools view -S -F 4 -c output/spikein/NPC_WT_H3K4me3_MG1655.sam > output/spikein/NPC_WT_H3K4me3-spikein_count.txt
```


--> There is some uniq mapped reads, around xxx% (in `003__CutRun` was less than 1%)

Now calculate SF in R, as for histone SF:


```R
# package
library("tidyverse")
library("readxl")
library("ggpubr")

# SF H3K27ac
spikein <- read_excel("output/spikein/SpikeIn_MG1655_008.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "H3K27ac")
# Total reads per IP
spikein_H3K27ac_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_H3K27ac_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_H3K27ac_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)



# SF H3K27me3
spikein <- read_excel("output/spikein/SpikeIn_MG1655_008.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "H3K27me3")
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
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_H3K27me3_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)


# SF H3K4me3
spikein <- read_excel("output/spikein/SpikeIn_MG1655_008.xlsx") %>%
    dplyr::select(-tissue) %>%
    filter(AB == "H3K4me3")
# Total reads per IP
spikein_H3K4me3_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_read_prop = spikein %>%
    left_join(spikein_H3K4me3_total) %>%
    mutate(read_prop = counts / total)
spikein_read_prop_min = spikein_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_scaling_factor = spikein_read_prop %>%
    left_join(spikein_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_H3K4me3_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)



```

--> histone vs MG1655 SF; same direction; histone a bit more *extreme*
----> GOOD!!




# Spike in scaling
## With MG1655 spike in for PTM CutRun


--> Let's use MG1655 as default method fopr spike in normalization! Seems more accurate and can be used for all AB!

**Using our scaling factor, let's estimate the 'new' library size** and provide it to `dba.normalize(library = c(1000, 12000))` = Like that our library size will be change taking into account our scaling factor! **Then we can normalize with library-size, RLE or TMM**... (issue discussed [here](https://support.bioconductor.org/p/9147040/)) 


### Adjust library size with MG1655 scaling factor and apply normalization
Total number of reads is our library size (used samtools flagstat to double check) :

`samtools flagstat output/bowtie2/*.dupmark.sorted.bam` used to obtain library size (first value=library size)
--> Values save in GoogleDrive `008__*/samples_008.xlsx`. Histone-norm-library-size = library-size * SF. Using the non-reciprocal scaling factor, we increase the library-size; the more histone enriched, the more library size is increased, thus the more signal will decrease.

Now let's use these new histone-scaled library size and normalize with library-size,TMM or RLE. Let's use the **unique bam files** together with the **unique bam MACS2 raw files (xlsx, not the bed with pre-filtered qvalue)**

***Key points:***
- **Let's do 1 DiffBind per AB (H3K27me3, H3K4me3,...) and tissue (PSC, NPC); otherwise the TMM normalization may take all, unrelated, samples into account!** --> Files are `meta_sample_macs2raw_unique*.txt`
- **For the non-histone CutRun, I will use the library size non histone scaled in DiffBind to collect TMM normalized SF**; I tested with and without specifying library size; and it does not change a lot the SF... Let's better use the one RiP method w/o providing the library size! Should provide BETTER correction


```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 

# ONE PER ONE
## NPC_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_H3K27me3.txt", header = TRUE, sep = "\t"))

### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)


## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3.RData")

### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()

pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist

sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)

### TMM 

sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(9181686,16861572), normalize = DBA_NORM_TMM) 

#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)


console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_H3K27me3.txt")


# NPC_H3K4me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_H3K4me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_H3K4me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K4me3.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(11433838,8912968), normalize = DBA_NORM_TMM) # 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_H3K4me3.txt")

# NPC_H3K27ac
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_H3K27ac.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27ac.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27ac.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(10279172,24996045), normalize = DBA_NORM_TMM) # 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_H3K27ac.txt")



# ALL TOGETHER FOR PCA/HEATMAP PLOT
## NPC
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_FACTOR, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_blackgreylist.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_blackgreylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_FACTOR, label=DBA_TREATMENT)
dev.off()
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(9181686,16861572,11433838,8912968,10279172,24996045), normalize = DBA_NORM_TMM) # 
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_NPC_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_NPC_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_FACTOR, label=DBA_TREATMENT)
dev.off()


```




# Generate Spike in scaled bigwig

--> Reciprocal from DiffBind_TMM is to be used when converting bam to bigwig!


### MG1655/E coli scaled bigwig


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_MG1655_DiffBind_TMM.sh # 13532257 fail; 13548881
```



Check some known target regulated in 2months neurons:
--> NEUROG2 seems less in KO which is good.
--> EFNA5 tiny decrease in KO (only in normalized data!)
--> GRIK3 tiny increase in KO

--> Something is WEIRD... When samples have MORE spike in, their signal should be reduced, as they overall have more DNA; but if I used the reciprocal from DiffBind_TMM; this is not respected (ie. sample with more spike in, we increased their signal...!)... That is true for both histone/MG1655-spike in DiffBind TMM norm...
----> What should be the BEST to use, is then the NON-reciprocal_DiffBind_TMM !!!
------> Let's try and compare with gene expression...! Maybe it is still good as we take into account the library size with the DiffBind_TMM method??
--------> YESSS in the end we correct the library size with the SF!!!!! So we 're good!!! reciprocal DiffBind_TMM IS TO BE USED!!


# THOR

Let's use THOR, notably to have IGG scaled bigwig...!

Comparison to do; NPC WT vs KO:
- H3K27me3
- H3K4me3
- H3K27ac


--> SF to use in THOR are the **reciprocal of MG1655_DiffBind_TMM**
--> Configs file created manually as `output/THOR/NPC_EZH2.config`

--> Lets also try to use the DiffBind spike in BAM method (similarly use the reciprocal from diffBind)




*THOR is very buggy to make it work I need to temporaly change where to look for libraries lol.. So cannot use nano anymore for example...*

*Follow these parameters: `WTvsHET_unique_Keepdup` (perform best in previous CutRun)*

```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge

# AB per AB
sbatch scripts/THOR_NPC_H3K27me3.sh # 13533943 fail; 13548195
sbatch scripts/THOR_NPC_H3K4me3.sh # 13534301 fail; 13535831 fail; 13548461
sbatch scripts/THOR_NPC_H3K27ac.sh # 13534612 fail; 13535918 fail; 13548462

```






## Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks

XXX Fuck that for now, if need to do it see `005`


# deepTools plot

On all genes, compare raw, DiffBind_TMM, THOR bigwigs


```bash
conda activate deeptools

# H3K27me3
sbatch scripts/matrix_TSS_10kb_H3K27me3_raw_allGenes.sh # 13550806 fail; 13562861 xxx
sbatch scripts/matrix_TSS_10kb_H3K27me3_DiffBindTMM_allGenes.sh # 13551680 fail; 13585180 xxx
sbatch scripts/matrix_TSS_10kb_H3K27me3_THOR_allGenes.sh # 13550735 fail; 13585432 xxx

# H3K4me3
sbatch scripts/matrix_TSS_10kb_H3K4me3_raw_allGenes.sh # 13550882 fail; 13585574 xxx
sbatch scripts/matrix_TSS_10kb_H3K4me3_DiffBindTMM_allGenes.sh # 13551681 fail; 13585759 xxx
sbatch scripts/matrix_TSS_10kb_H3K4me3_THOR_allGenes.sh # 13550762 fail; 13586013 xxx

# H3K27ac
sbatch scripts/matrix_TSS_10kb_H3K27ac_raw_allGenes.sh # 13550883 fail; 13586235 xxx
sbatch scripts/matrix_TSS_10kb_H3K27ac_DiffBindTMM_allGenes.sh # 13551682 fail; 13586308 xxx
sbatch scripts/matrix_TSS_10kb_H3K27ac_THOR_allGenes.sh # 13550775 fail; 13586565 xxx


```

--> H3K27me3 higher in KO (in agreement with 8wN data)

--> H3K4me3 similar in WT and KO

--> H3K27ac higher in KO (over -compensate the gain of H3K27me3?)





