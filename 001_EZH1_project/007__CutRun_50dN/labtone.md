# Project
- Neurons 50 days:
    - EF1a: EZH1, EZH2, H3K27me3, SUZ12, IGG
    - KO:  EZH1, EZH2, H3K27me3, SUZ12, IGG
    - Q731E (strong GOF):  EZH1, EZH2, H3K27me3, SUZ12, IGG
- PSC:
    - WT_ FA???: EZH1, H3K27me1
    - WT_ FA???: EZH1, H3K27me1


**Objectives:**
- Check whether Q731E indeed increase highly H3K27me3
- Check whether EZH1 AB work better in neurons
- test different FA condition in PSC for EZH1 AB 
- test different FA condition in PSC for H3K27me1 AB 




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

Go [there](http://data-deliver.novogene.com/batchfiles/X202SC23117109-Z01-F001) and enter credetnial: (check email Novogen)

I created a `nano url.txt` with all link and used `wget -i url.txt` to download them all (1 link per raw); then `mv input_raw_Novogene/*fq.gz input` .


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
sbatch scripts/fastp_1.sh # 9019614 ok
sbatch scripts/fastp_2.sh # 9019615 ok
sbatch scripts/fastp_3.sh # 9019616 ok
sbatch scripts/fastp_4.sh # 9019617 ok
```

# FastQC

**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:9019614 scripts/fastqc_fastp_1.sh # 9019677 ok
sbatch --dependency=afterany:9019615 scripts/fastqc_fastp_2.sh # 9019678 ok
sbatch --dependency=afterany:9019616 scripts/fastqc_fastp_3.sh # 9019680 ok
sbatch --dependency=afterany:9019617 scripts/fastqc_fastp_4.sh # 9019681 ok
```

--> all good


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:9019677 scripts/bowtie2_1.sh # 9020294 ok
sbatch --dependency=afterany:9019678 scripts/bowtie2_2.sh # 9020295 ok
sbatch --dependency=afterany:9019680 scripts/bowtie2_3.sh # 9020297 ok
sbatch --dependency=afterany:9019681 scripts/bowtie2_4.sh # 9020298 ok

```

--> Looks good

Mapping on E coli --> TO DO LATER! XXXXX

```bash
conda activate bowtie2

sbatch scripts/bowtie2_MG1655_1.sh # 9101160 ok
sbatch scripts/bowtie2_MG1655_2.sh # 9101602 ok
sbatch scripts/bowtie2_MG1655_3.sh # 9101606 ok

```

--> between 0.5 - 2% uniquely aligned reads (not a lot..; previously `005__CutRun` 10% (in `003__CutRun` was less than 1%) )


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-9020294.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_9020294.txt

for file in slurm-9020295.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_9020295.txt

for file in slurm-9020297.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_9020297.txt

for file in slurm-9020298.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_9020298.txt
```

Add these values to `/home/roulet/001_EZH1_project/007__CutRun_50dN/samples_007.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >72% input reads as been uniquely mapped to the genome (90% non uniq)





## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.

```bash
conda activate bowtie2

sbatch scripts/samtools_unique_1.sh # 9048648 ok
sbatch scripts/samtools_unique_2.sh # 9048649 ok
sbatch scripts/samtools_unique_3.sh # 9048650 ok
sbatch scripts/samtools_unique_4.sh # 9048651 ok

```

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch scripts/samtools_MG1655_unique_1.sh # 9162457
sbatch scripts/samtools_MG1655_unique_2.sh # 9162461
sbatch scripts/samtools_MG1655_unique_3.sh # 9162467
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

sbatch --dependency=afterany:9048648 scripts/bamtobigwig_unique_1.sh # 9048669 ok
sbatch --dependency=afterany:9048649 scripts/bamtobigwig_unique_2.sh # 9048672 ok
sbatch --dependency=afterany:9048650 scripts/bamtobigwig_unique_3.sh # 9048673 ok
sbatch --dependency=afterany:9048651 scripts/bamtobigwig_unique_4.sh # 9048674 ok
```



- KOEF1aEZH1
*Pass*: H3K27me3
*Failed*: EZH1cs, EZH2, SUZ12
- KO
*Pass*: H3K27me3
*Failed*: EZH1cs, EZH2, SUZ12
- WTQ731E
*Pass*: H3K27me3
*Failed*: EZH1cs, EZH2, SUZ12
- WT (PSC)
*Pass*: NA
*Failed*: EZH1cs and H3K27me1



## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_50dN.sh # 9064423 ok 
sbatch scripts/multiBigwigSummary_PSC.sh # 9064427 ok 
sbatch scripts/multiBigwigSummary_all.sh # 9064458 ok 


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels 50dN_KOEF1aEZH1_EZH1cs_R1 50dN_KOEF1aEZH1_EZH1cs_R2 50dN_KOEF1aEZH1_EZH2_R1 50dN_KOEF1aEZH1_EZH2_R2 50dN_KOEF1aEZH1_H3K27me3_R1 50dN_KOEF1aEZH1_H3K27me3_R2 50dN_KOEF1aEZH1_IGG_R1 50dN_KOEF1aEZH1_SUZ12_R1 50dN_KOEF1aEZH1_SUZ12_R2 50dN_KO_EZH1cs_R1 50dN_KO_EZH1cs_R2 50dN_KO_EZH2_R1 50dN_KO_EZH2_R2 50dN_KO_H3K27me3_R1 50dN_KO_H3K27me3_R2 50dN_KO_IGG_R1 50dN_KO_IGG_R2 50dN_KO_SUZ12_R1 50dN_KO_SUZ12_R2 50dN_WTQ731E_EZH1cs_R1 50dN_WTQ731E_EZH1cs_R2 50dN_WTQ731E_EZH2_R1 50dN_WTQ731E_EZH2_R2 50dN_WTQ731E_H3K27me3_R1 50dN_WTQ731E_H3K27me3_R2 50dN_WTQ731E_H3K27me3_R3 50dN_WTQ731E_IGG_R1 50dN_WTQ731E_IGG_R2 50dN_WTQ731E_SUZ12_R1 50dN_WTQ731E_SUZ12_R2 PSC_WT_EZH1cs_01FA PSC_WT_EZH1cs_1FA PSC_WT_H3K27me1_01FA PSC_WT_H3K27me1_1FA \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

plotPCA -in output/bigwig/multiBigwigSummary_PSC.npz \
    --transpose \
    --ntop 0 \
    --labels PSC_WT_EZH1cs_01FA PSC_WT_EZH1cs_1FA PSC_WT_H3K27me1_01FA PSC_WT_H3K27me1_1FA \
    -o output/bigwig/multiBigwigSummary_PSC_plotPCA.pdf

plotPCA -in output/bigwig/multiBigwigSummary_50dN.npz \
    --transpose \
    --ntop 0 \
    --labels 50dN_KOEF1aEZH1_EZH1cs_R1 50dN_KOEF1aEZH1_EZH1cs_R2 50dN_KOEF1aEZH1_EZH2_R1 50dN_KOEF1aEZH1_EZH2_R2 50dN_KOEF1aEZH1_H3K27me3_R1 50dN_KOEF1aEZH1_H3K27me3_R2 50dN_KOEF1aEZH1_IGG_R1 50dN_KOEF1aEZH1_SUZ12_R1 50dN_KOEF1aEZH1_SUZ12_R2 50dN_KO_EZH1cs_R1 50dN_KO_EZH1cs_R2 50dN_KO_EZH2_R1 50dN_KO_EZH2_R2 50dN_KO_H3K27me3_R1 50dN_KO_H3K27me3_R2 50dN_KO_IGG_R1 50dN_KO_IGG_R2 50dN_KO_SUZ12_R1 50dN_KO_SUZ12_R2 50dN_WTQ731E_EZH1cs_R1 50dN_WTQ731E_EZH1cs_R2 50dN_WTQ731E_EZH2_R1 50dN_WTQ731E_EZH2_R2 50dN_WTQ731E_H3K27me3_R1 50dN_WTQ731E_H3K27me3_R2 50dN_WTQ731E_H3K27me3_R3 50dN_WTQ731E_IGG_R1 50dN_WTQ731E_IGG_R2 50dN_WTQ731E_SUZ12_R1 50dN_WTQ731E_SUZ12_R2 \
    -o output/bigwig/multiBigwigSummary_50dN_plotPCA.pdf


## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels 50dN_KOEF1aEZH1_EZH1cs_R1 50dN_KOEF1aEZH1_EZH1cs_R2 50dN_KOEF1aEZH1_EZH2_R1 50dN_KOEF1aEZH1_EZH2_R2 50dN_KOEF1aEZH1_H3K27me3_R1 50dN_KOEF1aEZH1_H3K27me3_R2 50dN_KOEF1aEZH1_IGG_R1 50dN_KOEF1aEZH1_SUZ12_R1 50dN_KOEF1aEZH1_SUZ12_R2 50dN_KO_EZH1cs_R1 50dN_KO_EZH1cs_R2 50dN_KO_EZH2_R1 50dN_KO_EZH2_R2 50dN_KO_H3K27me3_R1 50dN_KO_H3K27me3_R2 50dN_KO_IGG_R1 50dN_KO_IGG_R2 50dN_KO_SUZ12_R1 50dN_KO_SUZ12_R2 50dN_WTQ731E_EZH1cs_R1 50dN_WTQ731E_EZH1cs_R2 50dN_WTQ731E_EZH2_R1 50dN_WTQ731E_EZH2_R2 50dN_WTQ731E_H3K27me3_R1 50dN_WTQ731E_H3K27me3_R2 50dN_WTQ731E_H3K27me3_R3 50dN_WTQ731E_IGG_R1 50dN_WTQ731E_IGG_R2 50dN_WTQ731E_SUZ12_R1 50dN_WTQ731E_SUZ12_R2 PSC_WT_EZH1cs_01FA PSC_WT_EZH1cs_1FA PSC_WT_H3K27me1_01FA PSC_WT_H3K27me1_1FA \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf

plotCorrelation \
    -in output/bigwig/multiBigwigSummary_PSC.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels PSC_WT_EZH1cs_01FA PSC_WT_EZH1cs_1FA PSC_WT_H3K27me1_01FA PSC_WT_H3K27me1_1FA \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_PSC_heatmap.pdf

plotCorrelation \
    -in output/bigwig/multiBigwigSummary_50dN.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels 50dN_KOEF1aEZH1_EZH1cs_R1 50dN_KOEF1aEZH1_EZH1cs_R2 50dN_KOEF1aEZH1_EZH2_R1 50dN_KOEF1aEZH1_EZH2_R2 50dN_KOEF1aEZH1_H3K27me3_R1 50dN_KOEF1aEZH1_H3K27me3_R2 50dN_KOEF1aEZH1_IGG_R1 50dN_KOEF1aEZH1_SUZ12_R1 50dN_KOEF1aEZH1_SUZ12_R2 50dN_KO_EZH1cs_R1 50dN_KO_EZH1cs_R2 50dN_KO_EZH2_R1 50dN_KO_EZH2_R2 50dN_KO_H3K27me3_R1 50dN_KO_H3K27me3_R2 50dN_KO_IGG_R1 50dN_KO_IGG_R2 50dN_KO_SUZ12_R1 50dN_KO_SUZ12_R2 50dN_WTQ731E_EZH1cs_R1 50dN_WTQ731E_EZH1cs_R2 50dN_WTQ731E_EZH2_R1 50dN_WTQ731E_EZH2_R2 50dN_WTQ731E_H3K27me3_R1 50dN_WTQ731E_H3K27me3_R2 50dN_WTQ731E_H3K27me3_R3 50dN_WTQ731E_IGG_R1 50dN_WTQ731E_IGG_R2 50dN_WTQ731E_SUZ12_R1 50dN_WTQ731E_SUZ12_R2 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_50dN_heatmap.pdf

```

--> two big groups: H3K27me3 IP versus the other
----> Seems only H3K27me3 IP has worked here




# MACS2 peak calling on bam unique

--> IGG samples used as control

--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad_1.sh # 9065217
sbatch scripts/macs2_broad_2.sh # 9065218
sbatch scripts/macs2_broad_3.sh # 9065220
sbatch scripts/macs2_broad_4.sh # xxx TO DO NO CNTROL !! XX

sbatch scripts/macs2_narrow_1.sh # 9065230
sbatch scripts/macs2_narrow_2.sh # 9065234
sbatch scripts/macs2_narrow_3.sh # 9065240
sbatch scripts/macs2_narrow_4.sh # xxx TO DO NO CNTROL !! XX
```

--> Very few peaks for all IP except H3K27me3... Technical issue...


```bash
conda activate bowtie2 # for bedtools
sbatch scripts/macs2_raw_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive

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

```bash
conda activate bowtie2 # for bedtools
sbatch scripts/macs2_raw_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive

# quick command to print median size of peak within a bed
awk '{print $3-$2}' your_bed_file.bed | sort -n | awk 'BEGIN {c=0; sum=0;} {a[c++]=$1; sum+=$1;} END {if (c%2) print a[int(c/2)]; else print (a[c/2-1]+a[c/2])/2;}'
```

**Optimal qvalue** according to IGV:
- 50dN_KOEF1aEZH1_H3K27me3: 1.30103 (2.3 more true peaks)
- 50dN_KO_H3K27me3: 1.30103 (2.3 more true peaks)
- 50dN_WTQ731E_H3K27me3: 1.30103 (2.3 more true peaks)



# Spike in factor

Let's do the analysis for H3K27me3 only; compare WT vs KO vs KOEF1a. Test 2 spikein normalization method (histone and Ecoli)


## Calculate histone content

--> This histone content will be used to generate a scaling factor which will be used to histone-scaled our library size. The calcul/method to follow is from `003__CutRun/output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt`

**Pipeline:**
- Count the histone barcode on the clean reads
- Calculate SF (group by sample (replicate) and AB and calculate the total nb of reads. Then proportion of reads = nb read in sample / total reads. SF = min(proportion) / sample proportion)


## Count the histone barcode on the clean reads



```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp.sh # 9128747 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_IGG_1.sh # 9156030 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_IGG_2.sh # 9156040 ok

sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_PSC.sh # 9155614 ok
```

--> It output the nb of reads found for each histone; then simply copy paste to the excell file `output/spikein/SpikeIn_QC_fastp_007.xlsx` in GoogleDrive

- *H3K27me3 50dN samples*: all samples well enriched for H3K27me3
- *H3K27me1 PSC samples*: all samples well enriched for H3K27me3


## histone spike in factor


```R
# package
library("tidyverse")
library("readxl")
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_007.xlsx") 

# 50dN
## H3K27me3
spikein_50dN_H3K27me3 = spikein %>%
    filter(tissue == "50dN",
           Target == "H3K27me3") %>%
    group_by(sample_ID, AB) %>%
    summarise(aligned=sum(counts))
# Total reads per IP
spikein_50dN_H3K27me3_total = spikein_50dN_H3K27me3 %>%
    ungroup() %>%
    group_by(AB) %>%
    mutate(total = sum(aligned)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_50dN_H3K27me3_read_prop = spikein_50dN_H3K27me3 %>%
    left_join(spikein_50dN_H3K27me3_total) %>%
    mutate(read_prop = aligned / total)
spikein_50dN_H3K27me3_read_prop_min = spikein_50dN_H3K27me3_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_50dN_H3K27me3_scaling_factor = spikein_50dN_H3K27me3_read_prop %>%
    left_join(spikein_50dN_H3K27me3_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_50dN_H3K27me3_scaling_factor, file="output/spikein/spikein_histone_50dN_H3K27me3_scaling_factor_fastp.txt", sep="\t", quote=FALSE, row.names=FALSE)


```

--> KO has a high SF!! High proportion of histone reads; so means that we will decrease it's signal! Which is good




### Quality control plot

Then look at the xlsx file from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) to generate quality control plot. Use R cluster for vizualization (file is `spikein_QC.xlsx` in Google Drive), file in `output/spikein`.
```R
# package
library("tidyverse")
library("readxl")
# import df adn tidy to remove AB used in sample_ID
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_007.xlsx") %>%
  separate(sample_ID, into = c("type", "condition", "tag"), sep = "_") %>%
  mutate(sample_ID = paste(type, condition, sep = "_")) %>%
  select(-type, -condition, -tag, -tissue)




# PSC
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
pdf("output/spikein/QC_histone_spike_in_H3K27me3.pdf", width = 14, height = 4)
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


## Histone scaling for H3K27me1
spikein_all_scale = spikein_all %>%
  group_by(sample_ID) %>%
  # Find the target_norm value when Target is H3K27me1 and AB is H3K27me1
  mutate(scaling_factor = ifelse(Target == "H3K27me1" & AB == "H3K27me1", target_norm, NA)) %>%
  # Fill the scaling_factor column with the appropriate value within each group
  fill(scaling_factor, .direction = "downup") %>%
  # Scale the target_norm values
  mutate(scaled_target_norm = target_norm / scaling_factor * 100) %>%
  # Remove the scaling_factor column
  select(-scaling_factor) %>%
  # Ungroup the data
  ungroup()
# Plot
pdf("output/spikein/QC_histone_spike_in_H3K27me1.pdf", width = 14, height = 4)
spikein_all_scale %>%
    filter(
           AB %in% c("H3K27me1", "IGG")) %>%
        ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
        geom_col(position = "dodge") +
        facet_wrap(~sample_ID, nrow=1) +
        geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


```


--> All good H3K27me3 IP H3K27me3 and H3K27me1 IP H3K27me1





## Ecoli MG1655 spike in factor


- Using the number of aligned reads to the spike-in control sequences for each sample `samtools view -S -F 4 -c sample_h3k27me3_rep1_spikein.sam > sample_h3k27me3_rep1_spikein_count.txt` --> *has been run in bowtie2 scripts* (**construct table on Google Drive**)
- Do the math for scaling factor, same method as when using histone spike-in

Now calculate SF in R, as for histone SF:

XXXXX


```R
# package
library("tidyverse")
library("readxl")
library("ggpubr")

# import df
spikein <- read_excel("output/spikein/SpikeIn_MG1655.xlsx") 

## 50dN
spikein <- read_table("output/spikein/SpikeIn_MG1655.txt") %>%
    filter(tissue == "NPC") %>%
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
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_NPC_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)


```







# depTools plot

deepTools plot to check H3K27me3 signal and compare the 3ng (R1, R2) and 6ng input used



```bash
conda activate deeptools

sbatch scripts/matrix_TSS_10kb_bigwig_unique_50dN_WT_H3K27me3.sh # 9760025 ok
sbatch scripts/matrix_TSS_10kb_bigwig_unique_50dN_WT_H3K27me3__ENCFF159KBI_WTpeaks_Promoter_5.sh # 9760026 ok

# Check median size of peaks

awk '{print $3-$2}' output/macs2/broad/broad_blacklist_qval1.30103/50dN_WTQ731E_H3K27me3_R1_peaks.broadPeak | sort -n | awk ' { a[i++]=$1; } END { x=int(i/2); if(i%2){print a[x];}else{print (a[x-1]+a[x])/2;} }' # 810
awk '{print $3-$2}' output/macs2/broad/broad_blacklist_qval1.30103/50dN_WTQ731E_H3K27me3_R2_peaks.broadPeak | sort -n | awk ' { a[i++]=$1; } END { x=int(i/2); if(i%2){print a[x];}else{print (a[x-1]+a[x])/2;} }' # 910
awk '{print $3-$2}' output/macs2/broad/broad_blacklist_qval1.30103/50dN_WTQ731E_H3K27me3_R3_peaks.broadPeak | sort -n | awk ' { a[i++]=$1; } END { x=int(i/2); if(i%2){print a[x];}else{print (a[x-1]+a[x])/2;} }' # 822

```


--> Overall the signal intensity is very comparable between samples, R1 a bit lower than R2 and R3...

--> The pearson correlation bigwig show a higher correlation between R1 R2 than with R3

--> peaks overall the same size between 3 and 6ng

--> Overall variability seems comparable to biological replicates




