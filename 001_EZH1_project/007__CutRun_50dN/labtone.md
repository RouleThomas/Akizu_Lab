# Project
- Neurons 50 days:
    - EF1a: EZH1, EZH2, H3K27me3, SUZ12, IGG
    - KO:  EZH1, EZH2, H3K27me3, SUZ12, IGG
    - WTQ731E:  EZH1, EZH2, H3K27me3, SUZ12, IGG
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


# H3K27me3 all 2 replicates together
sbatch scripts/macs2_broad_H3K27me3_pool.sh # 11037847 ok


```

--> Very few peaks for all IP except H3K27me3... Technical issue...


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




```R
# package
library("tidyverse")
library("readxl")
library("ggpubr")

# import df
spikein <- read_excel("output/spikein/SpikeIn_MG1655_007.xlsx")  %>%
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
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_50dN_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)


```

## Spike in Diffbind calcluation
### Histone spike in


**Using our scaling factor, let's estimate the 'new' library size** and provide it to `dba.normalize(library = c(1000, 12000))` = Like that our library size will be change taking into account our scaling factor! **Then we can normalize with library-size, RLE or TMM**... (issue discussed [here](https://support.bioconductor.org/p/9147040/)) 


### Adjust library size with histone scaling factor and apply normalization
Total number of reads is our library size (used samtools flagstat to double check) :

`samtools flagstat output/bowtie2/*unique.dupmark.sorted.bam` used to obtain library size (first value=library size)
--> Values save in GoogleDrive `007__*/sample_007.xlsx`. Histone-norm-library-size = library-size * SF. Using the non-reciprocal scaling factor, we increase the library-size; the more histone enriched, the more library size is increased, thus the more signal will decrease.

Now let's use these new histone-scaled library size and normalize with library-size,TMM or RLE. Let's use the **unique bam files** together with the **unique bam MACS2 raw files (xlsx, not the bed with pre-filtered qvalue)**

***Key points:***
- **Let's do 1 DiffBind per AB (H3K27me3) and tissue (50dN); otherwise the TMM normalization may take all, unrelated, samples into account!** --> Files are `meta_sample_macs2raw_unique*.txt`
- For **50dN_WTQ731E_H3K27me3_R3 I take 50dN_WTQ731E_H3K27me3_R1 IGG as control!**

```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 

# All 50dN H3K27me3
## 50dN_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_50dN_H3K27me3_histoneSF.txt", header = TRUE, sep = "\t"))

### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)


## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_50dN_H3K27me3_histoneSF.RData")
load("output/DiffBind/sample_count_macs2raw_unique_50dN_H3K27me3_histoneSF.RData")

### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_50dN_H3K27me3.pdf", width=14, height=20)  
plot(sample_count)
dev.off()

pdf("output/DiffBind/PCA_sample_macs2raw_unique_50dN_H3K27me3.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist

sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)

### TMM 

sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(8533246,15652994,47009725,51026263,26083225,27542538,30204595), normalize = DBA_NORM_TMM) 

#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)


console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_50dN_H3K27me3.txt")


### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_50dN_H3K27me3_histoneSF.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_50dN_H3K27me3_histoneSF.pdf", width=14, height=20) 
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



load("output/DiffBind/sample_count_macs2raw_unique_50dN_H3K27me3_histoneSF.RData")


### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist

sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(8533246,10079430,32469281,25669150,15181711,18083570,13174869), normalize = DBA_NORM_TMM) 

#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)

console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_50dN_H3K27me3.txt")


### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_50dN_H3K27me3_MG1655SF.pdf", width=14, height=20)  
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_50dN_H3K27me3_MG1655SF.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

```



--> SF from histone and MG1655 are very comparable

--> No striking difference in histone SF between H3K27me3 R1 and R3




## Generate Spike in scaled bigwig

--> Reciprocal from DiffBind_TMM is to be used when converting bam to bigwig!


### Histone scaled Bigwig (NON_Diffbind_TMM)


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_histone_1.sh # 9846985 ok
sbatch scripts/bamtobigwig_histone_2.sh # 9846998 ok
sbatch scripts/bamtobigwig_histone_3.sh # 9847001 ok
```

### MG1655 scaled Bigwig (NON_Diffbind_TMM)


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_MG1655_1.sh # 9847015 ok
sbatch scripts/bamtobigwig_MG1655_2.sh # 9847016 ok
sbatch scripts/bamtobigwig_MG1655_3.sh # 9847017 ok
```




### Histone scaled Bigwig DiffBind_TMM


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_histone_DiffBind_TMM_1.sh # 9846509 ok
sbatch scripts/bamtobigwig_histone_DiffBind_TMM_2.sh # 9846512 ok
sbatch scripts/bamtobigwig_histone_DiffBind_TMM_3.sh # 9846519 ok
```

Generate median tracks:
```bash
conda activate BedToBigwig

sbatch scripts/bigwigmerge_histone_DiffBind_TMM.sh # 10922925 ok
```

### MG1655/E coli scaled bigwig DiffBind_TMM


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_MG1655_DiffBind_TMM_1.sh # 9846532 ok
sbatch scripts/bamtobigwig_MG1655_DiffBind_TMM_2.sh # 9846537 ok
sbatch scripts/bamtobigwig_MG1655_DiffBind_TMM_3.sh # 9846538 ok
```


Generate median tracks:
```bash
conda activate BedToBigwig

sbatch scripts/bigwigmerge_MG1655_DiffBind_TMM.sh # 10922976 ok
```

--> Both bigwig norm method are very similar... 

--> **reciprocal DiffBind_TMM IS TO BE USED**!!




# depTools plot

## raw
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


## spike in with raw SF (NON- DiffBind TMM normalized)


```bash
conda activate deeptools

## all
sbatch --dependency=afterany:9846985:9846998:9847001 scripts/matrix_TSS_10kb_bigwig_unique_histone_50dN_H3K27me3.sh # 9847384 ok
sbatch --dependency=afterany:9847015:9847016:9847017 scripts/matrix_TSS_10kb_bigwig_unique_MG1655_50dN_H3K27me3.sh # 9847393 ok

## WT only
sbatch --dependency=afterany:9846985:9846998:9847001 scripts/matrix_TSS_10kb_bigwig_unique_histone_50dN_WT_H3K27me3.sh # 9847400 ok
sbatch --dependency=afterany:9847015:9847016:9847017 scripts/matrix_TSS_10kb_bigwig_unique_MG1655_50dN_WT_H3K27me3.sh # 9847410 ok
```


--> all samples are weird (notably WT). Profile are ugly and result changed when using histone or MG1655 -spike in SF...



## spike in (histone and MG1655) DiffBind Lib scaled (method 003__CutRun)


```bash
conda activate deeptools

## all
sbatch --dependency=afterany:9846509:9846512:9846519 scripts/matrix_TSS_10kb_bigwig_unique_histone_DiffBind_TMM_50dN_H3K27me3.sh # 9846649 ok
sbatch --dependency=afterany:9846532:9846537:9846538 scripts/matrix_TSS_10kb_bigwig_unique_MG1655_DiffBind_TMM_50dN_H3K27me3.sh # 9846669 ok

## WT only
sbatch --dependency=afterany:9846509:9846512:9846519 scripts/matrix_TSS_10kb_bigwig_unique_histone_DiffBind_TMM_50dN_WT_H3K27me3.sh # 9846767 ok
sbatch --dependency=afterany:9846532:9846537:9846538 scripts/matrix_TSS_10kb_bigwig_unique_MG1655_DiffBind_TMM_50dN_WT_H3K27me3.sh # 9846779 ok

## median H3K27me3
sbatch scripts/matrix_TSS_10kb_bigwig_unique_histone_DiffBind_TMM_50dN_H3K27me3_median.sh # 10929886 ok
sbatch scripts/matrix_TSS_10kb_bigwig_unique_MG1655_DiffBind_TMM_50dN_H3K27me3_median.sh # 10929206 ok



```


--> WT result are good! Samples are more closely related and R3 (6ng) seems to have a bit more signal than R1 and R2 (3ng). Overall very comparable histone and MG1655 -norm.

--> H3K27me3 level: KOEF1a > WT > KO (overall, as expected!!)

--> The profile is a bit flat, lowly enriched/higher around the TSS; maybe because non IGG norm?


## THOR - with SF spike in (MG1655) DiffBind TMM scaled (method 003__CutRun)


```bash
conda activate deeptools

## H3K27me3 all replicates
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3.sh # 10988809 ok

## H3K27me3 median
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median.sh # 11001463 ok

## H3K27me3 median with threshold 5 / 10- FAIL value remove!
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_threshold5.sh # 11037365 ok
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_threshold10.sh # 11037550 ok

## H3K27me3 median; genes with peak in WT only
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_WTH3K27me3peaks_median.sh # 11039919 ok

## H3K27me3 median; without skipZeros argument
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_noskipZeros.sh # 11041574 ok

## H3K27me3 median; with treshold 5 / 10 - bigwig modified (bedGraph)
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_bedGraph_threshold5.sh # 11043468 ok
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_bedGraph_threshold10.sh # 11043473 ok

## H3K27me3 median; raw THOR bigwig; region that gain vs lost H3K27me3
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_THORq10_positive_negative.sh # 11045278 ok
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_THORq15_positive_negative.sh # 11045279 ok
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_THORq20_positive_negative.sh # 11084288 ok

sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKOEF1aEZH1_THORq10_positive_negative.sh # 11045282 ok
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKOEF1aEZH1_THORq15_positive_negative.sh # 11045283 ok
sbatch scripts/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKOEF1aEZH1_THORq20_positive_negative.sh # 11084527 ok

```

--> The 2 replicates are comparable

--> The signal is NOT more centered into the TSS...
----> Let's try to add a threshold to the bigiwg, like not take signal under 5 and 10
------> Taking a threshold is not working, I got an oscilating profile WTF!!! and very few peaks. The issue is that the `--minThreshold` argument skip the value! I need to set them at 0...
---------> Maybe not use the argument `--skipZeros` ! Do not change anyhting!

Let's isolate the genes with peak in WT; create gtf and redo deepTool plots.
--> better, but still not a strong peak aound TSS

--> Using clean bigwig (removing low values) improve a bit, but we lose the diff between WT and KO... Overall using the **peaks in WT only is better (but still weird)...**




# clean bigwig file (remove low value)

The deeptool profile is weird probably due to high level of noise. Let's try to renmove the low value; assign value under 5 to 0:
- convert bigwig into bedgraph
- remove low value
- re-convert into bigwig


```bash
conda activate BedToBigwig

# convert bigwig to bedGrah
bigWigToBedGraph output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bw output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bedGraph
bigWigToBedGraph output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bw output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bedGraph
bigWigToBedGraph output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bw output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bedGraph

# put value less than 5 to 0
awk '{ if ($4 <= 5) $4 = 0; print }' output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bedGraph > output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median_threshold5.bedGraph
awk '{ if ($4 <= 5) $4 = 0; print }' output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bedGraph > output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median_threshold5.bedGraph
awk '{ if ($4 <= 5) $4 = 0; print }' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bedGraph > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median_threshold5.bedGraph

# put value less than 10 to 0
awk '{ if ($4 <= 10) $4 = 0; print }' output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bedGraph > output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median_threshold10.bedGraph
awk '{ if ($4 <= 10) $4 = 0; print }' output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bedGraph > output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median_threshold10.bedGraph
awk '{ if ($4 <= 10) $4 = 0; print }' output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bedGraph > output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median_threshold10.bedGraph


# re-convert bedGraph to bigwig
bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median_threshold5.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median_threshold5.bw
bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median_threshold5.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median_threshold5.bw
bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median_threshold5.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median_threshold5.bw

bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median_threshold10.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median_threshold10.bw
bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median_threshold10.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median_threshold10.bw
bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median_threshold10.bedGraph ../../Master/meta/GRCh38_chrom_sizes.tab output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median_threshold10.bw
```







# ChIPseeker peak gene assignment

## From optimal qval bed files peaks
Let's assign **peak to genes from MACS2 peak**:

**Optimal qvalue** according to IGV:
- 50dN_KOEF1aEZH1_H3K27me3: 1.30103 (2.3 more true peaks)
- 50dN_KO_H3K27me3: 1.30103 (2.3 more true peaks)
- 50dN_WTQ731E_H3K27me3: 1.30103 (2.3 more true peaks)



```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")



# Import macs2 peaks
## 50dN
WTQ731E = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/50dN_WTQ731E_H3K27me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
KO = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/50dN_KO_H3K27me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)    
KOEF1aEZH1 = as_tibble(read.table('output/macs2/broad/broad_blacklist_qval1.30103/50dN_KOEF1aEZH1_H3K27me3_pool_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
     

# Tidy peaks #-->> Re-Run from here with different qvalue!!
## 50dN
WTQ731E_gr = makeGRangesFromDataFrame(WTQ731E,keep.extra.columns=TRUE)
KO_gr = makeGRangesFromDataFrame(KO,keep.extra.columns=TRUE)
KOEF1aEZH1_gr = makeGRangesFromDataFrame(KOEF1aEZH1,keep.extra.columns=TRUE)
gr_list <- list(WTQ731E=WTQ731E_gr, KO=KO_gr,  KOEF1aEZH1=KOEF1aEZH1_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
WTQ731E_annot <- as.data.frame(peakAnnoList[["WTQ731E"]]@anno)
KO_annot <- as.data.frame(peakAnnoList[["KO"]]@anno)
KOEF1aEZH1_annot <- as.data.frame(peakAnnoList[["KOEF1aEZH1"]]@anno)


## Convert entrez gene IDs to gene symbols
WTQ731E_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WTQ731E_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
WTQ731E_annot$gene <- mapIds(org.Hs.eg.db, keys = WTQ731E_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KO_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KO_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KO_annot$gene <- mapIds(org.Hs.eg.db, keys = KO_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KOEF1aEZH1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KOEF1aEZH1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KOEF1aEZH1_annot$gene <- mapIds(org.Hs.eg.db, keys = KOEF1aEZH1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



## Save output table
write.table(WTQ731E_annot, file="output/ChIPseeker/annotation_macs2_WTQ731E_H3K27me3_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(KO_annot, file="output/ChIPseeker/annotation_macs2_KO_H3K27me3_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(KOEF1aEZH1_annot, file="output/ChIPseeker/annotation_macs2_KOEF1aEZH1_annot_H3K27me3_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
WTQ731E_annot_promoterAnd5 = tibble(WTQ731E_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
KO_annot_promoterAnd5 = tibble(KO_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
KOEF1aEZH1_annot_promoterAnd5 = tibble(KOEF1aEZH1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))


### Save output gene lists
WTQ731E_annot_promoterAnd5_geneSymbol = WTQ731E_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
KO_annot_promoterAnd5_geneSymbol = KO_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
KOEF1aEZH1_annot_promoterAnd5_geneSymbol = KOEF1aEZH1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()



write.table(WTQ731E_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_WTQ731E_H3K27me3_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(KO_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_KO_H3K27me3_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(KOEF1aEZH1_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_KOEF1aEZH1_H3K27me3_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


# Comparison peak position WT vs KO vs KOEF1aEZH1
## plots
pdf("output/ChIPseeker/plotAnnoBar_H3K27me3.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()


pdf("output/ChIPseeker/plotDistToTSS_H3K27me3.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()


```
**50dN**:
- Export gene list of genes with peak in WT
- Generate gtf 
- do deepTool plot



```bash
# Generate gtf from gene Symbol list
### Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annot_macs2_WTQ731E_H3K27me3_qval1.30103_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annot_macs2_WTQ731E_H3K27me3_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt

### Filter the gtf
grep -Ff output/ChIPseeker/annot_macs2_WTQ731E_H3K27me3_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt ../003__CutRun/meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_WT_H3K27me3peaks.gtf

```





## From THOR diff peaks
Let's assign **peak to genes from THOR positive (gain) and negative (lost) peaks** (for qval 10 and 15):

```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")



# Import THOR diff peaks
## qval 10
KO_gain = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
KO_lost = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval10_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)       

KOEF1aEZH1_gain = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
KOEF1aEZH1_lost = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)   

## qval 15
KO_gain = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval15_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
KO_lost = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval15_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)       

KOEF1aEZH1_gain = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval15_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
KOEF1aEZH1_lost = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval15_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)   

## qval 20
KO_gain = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval20_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
KO_lost = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval20_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)       

KOEF1aEZH1_gain = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval20_positive.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
KOEF1aEZH1_lost = as_tibble(read.table('output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval20_negative.bed') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)  

# Tidy peaks #-->> Re-Run from here with different qvalue!!
## 50dN
KO_gain_gr = makeGRangesFromDataFrame(KO_gain,keep.extra.columns=TRUE)
KO_lost_gr = makeGRangesFromDataFrame(KO_lost,keep.extra.columns=TRUE)
KOEF1aEZH1_gain_gr = makeGRangesFromDataFrame(KOEF1aEZH1_gain,keep.extra.columns=TRUE)
KOEF1aEZH1_lost_gr = makeGRangesFromDataFrame(KOEF1aEZH1_lost,keep.extra.columns=TRUE)

gr_list <- list(KO_gain=KO_gain_gr, KO_lost=KO_lost_gr,  KOEF1aEZH1_gain=KOEF1aEZH1_gain_gr,  KOEF1aEZH1_lost=KOEF1aEZH1_lost_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
KO_gain_annot <- as.data.frame(peakAnnoList[["KO_gain"]]@anno)
KO_lost_annot <- as.data.frame(peakAnnoList[["KO_lost"]]@anno)
KOEF1aEZH1_gain_annot <- as.data.frame(peakAnnoList[["KOEF1aEZH1_gain"]]@anno)
KOEF1aEZH1_lost_annot <- as.data.frame(peakAnnoList[["KOEF1aEZH1_lost"]]@anno)


## Convert entrez gene IDs to gene symbols
KO_gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KO_gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KO_gain_annot$gene <- mapIds(org.Hs.eg.db, keys = KO_gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KO_lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KO_lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KO_lost_annot$gene <- mapIds(org.Hs.eg.db, keys = KO_lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KOEF1aEZH1_gain_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KOEF1aEZH1_gain_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KOEF1aEZH1_gain_annot$gene <- mapIds(org.Hs.eg.db, keys = KOEF1aEZH1_gain_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KOEF1aEZH1_lost_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KOEF1aEZH1_lost_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KOEF1aEZH1_lost_annot$gene <- mapIds(org.Hs.eg.db, keys = KOEF1aEZH1_lost_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(KO_gain_annot, file="output/ChIPseeker/annotation_THOR_KO_gain_annot_qval20.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(KO_lost_annot, file="output/ChIPseeker/annotation_THOR_KO_lost_annot_qval20.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(KOEF1aEZH1_gain_annot, file="output/ChIPseeker/annotation_THOR_KOEF1aEZH1_gain_annot_qval20.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(KOEF1aEZH1_lost_annot, file="output/ChIPseeker/annotation_THOR_KOEF1aEZH1_lost_annot_qval20.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE


## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
KO_gain_annot_promoterAnd5 = tibble(KO_gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
KO_lost_annot_promoterAnd5 = tibble(KO_lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
KOEF1aEZH1_gain_annot_promoterAnd5 = tibble(KOEF1aEZH1_gain_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
KOEF1aEZH1_lost_annot_promoterAnd5 = tibble(KOEF1aEZH1_lost_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))

### Save output gene lists
KO_gain_annot_promoterAnd5_geneSymbol = KO_gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
KO_lost_annot_promoterAnd5_geneSymbol = KO_lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
KOEF1aEZH1_gain_annot_promoterAnd5_geneSymbol = KOEF1aEZH1_gain_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
KOEF1aEZH1_lost_annot_promoterAnd5_geneSymbol = KOEF1aEZH1_lost_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()


write.table(KO_gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_THOR_KO_gain_qval20_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(KO_lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_THOR_KO_lost_qval20_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(KOEF1aEZH1_gain_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_THOR_KOEF1aEZH1_gain_qval20_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(KOEF1aEZH1_lost_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_THOR_KOEF1aEZH1_lost_qva20_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

# Comparison peak position WT vs KO vs KOEF1aEZH1
## plots
pdf("output/ChIPseeker/plotAnnoBar_H3K27me3_THORq20_gain_lost.pdf", width = 8, height = 3)
plotAnnoBar(peakAnnoList)
dev.off()


pdf("output/ChIPseeker/plotDistToTSS_H3K27me3_THORq20_gain_lost.pdf", width = 8, height = 3)
plotDistToTSS(peakAnnoList, title="Distribution relative to TSS")
dev.off()
```





# THOR diff peaks

Let's use THOR, notably to have IGG scaled bigwig...!

Comparison to do; 50dN WT vs KO and WT vs KOEF1a:
- H3K27me3

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

# AB per AB (DiffBind TMM MG1655)
sbatch scripts/THOR_50dN_H3K27me3_WTvsKO.sh # 10983795 ok
sbatch scripts/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1.sh # 10984013 ok


```

Generate median tracks:
```bash
conda activate BedToBigwig

sbatch scripts/bigwigmerge_MG1655_THOR_DiffBind_TMM.sh # 11001294 ok
```



--> THOR bigwig tracks looks good! Prevous finding from `003__CutRun` (NEUROG2, GRIK3, GRIN1, EFNA5) found here too!!



## Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks


```R

# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# WTvsKO
diffpeaks <- read_tsv("output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_50dN_H3K27me3_WTvsKO/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("50dN_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_50dN_H3K27me3_WTvsKO/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO_qval25") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 25) %>%
  write_tsv("output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval25.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())




# WTvsKOEF1aEZH1
diffpeaks <- read_tsv("output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KOEF1aEZH1", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KOEF1aEZH1, into = c("count_KOEF1aEZH1_1","count_KOEF1aEZH1_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KOEF1aEZH1_1+count_KOEF1aEZH1_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("50dN_WT vs KOEF1aEZH1") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/log2FC_qval10.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 10) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO_qval10") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 10) %>%
  write_tsv("output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())
```

- *NOTE: FC positive = less in KO; negative = more in KO*

**Optimal qvalue:**
--> *H3K27me3*; qval 10 looks great!

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
```




# Functional analysis

## GO


### regular GO - for unique list of genes (not enrichR)

- H3K27me3 gain in KOEF1aEZH1 and lost in KO



```bash
conda activate deseq2
```

```R
# packages
library("clusterProfiler")
library("pathview")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("rtracklayer")
library("tidyverse")
library("biomaRt")

# import gene lists
list = read_csv("output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_KO_lost_qval10_promoterAnd5_geneSymbol_224.txt", col_names = FALSE)
list = read_csv("output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_KO_lost_qval15_promoterAnd5_geneSymbol_135.txt", col_names = FALSE)
list = read_csv("output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_KO_lost_qval20_promoterAnd5_geneSymbol_76.txt", col_names = FALSE)

## GO
ego <- enrichGO(gene = (list$X1), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)

pdf("output/GO/dotplot_BP_Venn_overlap_THOR_KOEF1aEZH1_gain_KO_lost_qval10_promoterAnd5_geneSymbol_224.pdf", width=5, height=12)
pdf("output/GO/dotplot_BP_Venn_overlap_THOR_KOEF1aEZH1_gain_KO_lost_qval20_promoterAnd5_geneSymbol_76.pdf", width=5, height=20)

dotplot(ego, showCategory = 50, title = "GO_Biological Process Enrichment Analysis")
dev.off()


GO_summary <- data.frame(ego)
write.csv(GO_summary, "output/GO/dotplot_BP_Venn_overlap_THOR_KOEF1aEZH1_gain_KO_lost_qval20_promoterAnd5_geneSymbol_76.csv")

## KEGG
### convert geneSymbol to entrezID
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#### Get the Entrez IDs corresponding to the gene symbols
genes_entrez <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                      filters = "hgnc_symbol",
                      values = list$X1, 
                      mart = mart)
#### Merge to retain the original order and to include all gene symbols
list_entrezid <- as_tibble(genes_entrez) 


ekegg <- enrichKEGG(gene = (list_entrezid$entrezgene_id), 
                keyType = "kegg",     # Use ENSEMBL if want to use ENSG000XXXX format
                organism = "hsa", 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05)

pdf("output/GO/dotplot_KEGG_Venn_overlap_THOR_KOEF1aEZH1_gain_KO_lost_qval15_promoterAnd5_geneSymbol_135.pdf", width=5, height=2)
dotplot(ekegg, showCategory = 15, title = "KEGG pathway Enrichment Analysis")
dev.off()



```


### enrichR GO - for up/down list of genes 


- specifically gain or down in KOEF1aEZH1 and KO


```R
# libr
library("tidyverse")
library("enrichR")

# Define databases for enrichment
dbs <- c("GO_Biological_Process_2023") # 


### GeneSymbol list of gain and lost
#### KO
output/ChIPseeker/Venn_overlap_THOR_KO_lost_specific_qval10_promoterAnd5_geneSymbol_1459.txt
output/ChIPseeker/Venn_overlap_THOR_KO_gain_specific_qval10_promoterAnd5_geneSymbol_2121.txt
output/ChIPseeker/Venn_overlap_THOR_KO_lost_specific_qval15_promoterAnd5_geneSymbol_1128.txt
output/ChIPseeker/Venn_overlap_THOR_KO_gain_specific_qval15_promoterAnd5_geneSymbol_1883.txt
output/ChIPseeker/Venn_overlap_THOR_KO_lost_specific_qval20_promoterAnd5_geneSymbol_901.txt
output/ChIPseeker/Venn_overlap_THOR_KO_gain_specific_qval20_promoterAnd5_geneSymbol_1631.txt
#### KOEF1aEZH1
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_lost_specific_qval10_promoterAnd5_geneSymbol_1459.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_specific_qval10_promoterAnd5_geneSymbol_4804.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_lost_specific_qval15_promoterAnd5_geneSymbol_1330.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_specific_qval15_promoterAnd5_geneSymbol_3265.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_lost_specific_qval20_promoterAnd5_geneSymbol_1087.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_specific_qval20_promoterAnd5_geneSymbol_2122.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_lost_specific_qval15_promoterAnd5_geneSymbol_1330.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_specific_qval15_promoterAnd5_geneSymbol_3265.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$GO_Biological_Process_2023
down <- edown$GO_Biological_Process_2023
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)
gos$Term <- gsub("\\(GO:[0-9]+\\)", "", gos$Term)  # Regular expression to match the GO pattern and replace it with an empty string

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)



# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_KO_specific_qval10_promoterAnd5_geneSymbol.pdf", width=10, height=10)
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_KO_specific_qval15_promoterAnd5_geneSymbol.pdf", width=10, height=10)
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_KO_specific_qval20_promoterAnd5_geneSymbol.pdf", width=10, height=10)

pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_KOEF1aEZH1_specific_qval10_promoterAnd5_geneSymbol.pdf", width=10, height=7)
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_KOEF1aEZH1_specific_qval15_promoterAnd5_geneSymbol.pdf", width=10, height=7)
pdf("output/GO/enrichR_GO_Biological_Process_2023_THOR_KOEF1aEZH1_specific_qval20_promoterAnd5_geneSymbol.pdf", width=10, height=7)


ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Peaks", 
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for GO BP") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_GO_Biological_Process_2023_THOR_KOEF1aEZH1_specific_qval20_promoterAnd5_geneSymbol.txt", sep="\t", row.names=FALSE, quote=FALSE)




# Define databases for enrichment
dbs <- c("KEGG_2021_Human") # 


### GeneSymbol list of gain and lost
#### KO
output/ChIPseeker/Venn_overlap_THOR_KO_lost_specific_qval10_promoterAnd5_geneSymbol_1459.txt
output/ChIPseeker/Venn_overlap_THOR_KO_gain_specific_qval10_promoterAnd5_geneSymbol_2121.txt
output/ChIPseeker/Venn_overlap_THOR_KO_lost_specific_qval15_promoterAnd5_geneSymbol_1128.txt
output/ChIPseeker/Venn_overlap_THOR_KO_gain_specific_qval15_promoterAnd5_geneSymbol_1883.txt
output/ChIPseeker/Venn_overlap_THOR_KO_lost_specific_qval20_promoterAnd5_geneSymbol_901.txt
output/ChIPseeker/Venn_overlap_THOR_KO_gain_specific_qval20_promoterAnd5_geneSymbol_1631.txt
#### KOEF1aEZH1
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_lost_specific_qval10_promoterAnd5_geneSymbol_1459.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_specific_qval10_promoterAnd5_geneSymbol_4804.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_lost_specific_qval15_promoterAnd5_geneSymbol_1330.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_specific_qval15_promoterAnd5_geneSymbol_3265.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_lost_specific_qval20_promoterAnd5_geneSymbol_1087.txt
output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_specific_qval20_promoterAnd5_geneSymbol_2122.txt


# IF starting with geneSymbol
## Read and preprocess data for downregulated genes
gene_names_down <- read.csv("output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_lost_specific_qval20_promoterAnd5_geneSymbol_1087.txt", header=FALSE, stringsAsFactors=FALSE)
list_down <- unique(as.character(gene_names_down$V1))
edown <- enrichr(list_down, dbs)

## Read and preprocess data for upregulated genes
gene_names_up <- read.csv("output/ChIPseeker/Venn_overlap_THOR_KOEF1aEZH1_gain_specific_qval20_promoterAnd5_geneSymbol_2122.txt", header=FALSE, stringsAsFactors=FALSE)
list_up <- unique(as.character(gene_names_up$V1))
eup <- enrichr(list_up, dbs)

# Extracting KEGG data and assigning types
up <- eup$KEGG_2021_Human
down <- edown$KEGG_2021_Human
up$type <- "up"
down$type <- "down"

# Get top enriched terms and sort by Combined.Score (Note: Adjust if you don't want the top 10)
up <- head(up[order(up$Combined.Score, decreasing = TRUE), ], 25)
down <- head(down[order(down$Combined.Score, decreasing = TRUE), ], 25)

# Convert adjusted p-values and differentiate direction for up and down
up$logAdjP <- -log10(up$Adjusted.P.value)
down$logAdjP <- -1 * -log10(down$Adjusted.P.value)

# Combine the two dataframes
gos <- rbind(down, up)
gos <- gos %>% arrange(logAdjP)

# Filter out rows where absolute logAdjP 1.3 = 0.05
gos <- gos %>% filter(abs(logAdjP) > 1.3)

# Create the order based on the approach given
up_pathways <- gos %>% filter(type == "up") %>% arrange(-logAdjP) %>% pull(Term)
down_pathways <- gos %>% filter(type == "down") %>% arrange(logAdjP) %>% pull(Term)
new_order <- c(down_pathways, up_pathways)
gos$Term <- factor(gos$Term, levels = new_order)


## FAIL as dupplicates:
up_pathways_suffixed <- paste0(up_pathways, "_up")
down_pathways_suffixed <- paste0(down_pathways, "_down")
new_order <- c(down_pathways_suffixed, up_pathways_suffixed)
gos$Term <- ifelse(gos$type == "up", paste0(gos$Term, "_up"), paste0(gos$Term, "_down"))
gos$Term <- factor(gos$Term, levels = new_order)

# Plotting with enhanced aesthetics
pdf("output/GO/enrichR_KEGG_2021_Human_THOR_KO_specific_qval10_promoterAnd5_geneSymbol.pdf", width=10, height=10)
pdf("output/GO/enrichR_KEGG_2021_Human_THOR_KO_specific_qval15_promoterAnd5_geneSymbol.pdf", width=10, height=10)
pdf("output/GO/enrichR_KEGG_2021_Human_THOR_KO_specific_qval20_promoterAnd5_geneSymbol.pdf", width=10, height=10)

pdf("output/GO/enrichR_KEGG_2021_Human_THOR_KOEF1aEZH1_specific_qval10_promoterAnd5_geneSymbol.pdf", width=10, height=5)
pdf("output/GO/enrichR_KEGG_2021_Human_THOR_KOEF1aEZH1_specific_qval15_promoterAnd5_geneSymbol.pdf", width=10, height=5)
pdf("output/GO/enrichR_KEGG_2021_Human_THOR_KOEF1aEZH1_specific_qval20_promoterAnd5_geneSymbol.pdf", width=10, height=5)

ggplot(gos, aes(x=Term, y=logAdjP, fill=type)) + 
  geom_bar(stat='identity', width=.7) +
  
  # Adjusted label position based on the type of gene (up/down) and increased separation
  geom_text(aes(label=Term, y=ifelse(type == "up", max(gos$logAdjP) + 2, min(gos$logAdjP) - 2)), hjust = ifelse(gos$type == "up", 1, 0), size = 7, color = "gray28") +
  
  geom_hline(yintercept = 0, linetype="solid", color = "black") +
  
  scale_fill_manual(name="Peaks", 
                    labels = c("Lost", "Gain"), 
                    values = c("down"="Sky Blue", "up"="Orange")) + 
  labs(title= "Diverging bars of -log10 Adjusted P-value for KEGG") + 
  coord_flip() + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15)
  )
dev.off()


## save output
write.table(gos, "output/GO/enrichR_KEGG_2021_Human_THOR_KOEF1aEZH1_specific_qval20_promoterAnd5_geneSymbol.txt", sep="\t", row.names=FALSE, quote=FALSE)



```


# Compare CutRun and RNAseq expression (from 2 months old neurons)


- import gene list that gain / lost H3K27me3 in KO: `output/ChIPseeker/annot_THOR_KO_*_qval*_promoterAnd5_geneSymbol.txt`
- join with expression data





```R
# packages
library("tidyverse")


XXXXXXXXXX below to mod




# Import RNAseq deseq2 output
## Raw FC ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_KO_vs_NPC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)
## Fitlered FC ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_KO_vs_NPC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)

# Merge files
H3K27me3_annot_gain_lost_RNA = H3K27me3_annot_gain_lost %>% 
    left_join(KO_vs_WT) %>%
    dplyr::select(gene, binding,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()
EZH2_annot_gain_lost_RNA = EZH2_annot_gain_lost %>% 
    left_join(KO_vs_WT) %>%
    dplyr::select(gene, binding,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()
SUZ12_annot_gain_lost_RNA = SUZ12_annot_gain_lost %>% 
    left_join(KO_vs_WT) %>%
    dplyr::select(gene, binding,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()
H3K4me3_annot_gain_lost_RNA = H3K4me3_annot_gain_lost %>% 
    left_join(KO_vs_WT) %>%
    dplyr::select(gene, binding,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()


# Volcano plot 
count_data <- H3K27me3_annot_gain_lost_RNA %>%     # CHANGE TITLE !!!!!!!
    group_by(binding, significance) %>%
    summarise(up = sum(log2FoldChange > 0),
              down = sum(log2FoldChange < 0),
              total = n()) %>%
    ungroup() %>%
    group_by(binding) %>%
    mutate(total_panel = sum(total)) %>%
    ungroup()


# pdf("output/ChIPseeker/THOR_qval20_H3K4me3_expression_promoterAnd5.pdf", width=7, height=4)  # CHANGE TITLE !!!!!!!
# pdf("output/ChIPseeker/THOR_qval20_H3K4me3_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE !!!!!!!

pdf("output/ChIPseeker/THOR_qval50_H3K27me3_expression_promoterAnd5_test.pdf", width=7, height=4)  # CHANGE TITLE !!!!!!!

H3K27me3_annot_gain_lost_RNA %>%         # CHANGE TITLE !!!!!!!
    ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
        geom_point(alpha = 0.8, size = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        labs(title = "KO vs WT",
             subtitle = "Expression level of diff. bound H3K27me3 genes",    # CHANGE TITLE !!!!!!!
             x = "Log2 Fold Change",
             y = "-log10(q-value)",
             color = "Significant (padj <= 0.05)") +
        facet_wrap(~binding) +
        theme_bw() +
        geom_text(data = count_data %>% filter(significance), 
                  aes(x = Inf, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
                  hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
        geom_text(data = count_data %>% distinct(binding, .keep_all = TRUE),
                  aes(x = Inf, y = -Inf, label = paste("Total:", total_panel, "genes")),
                  hjust = 1.1, vjust = -0.1, size = 3, color = "black")
dev.off()




```




