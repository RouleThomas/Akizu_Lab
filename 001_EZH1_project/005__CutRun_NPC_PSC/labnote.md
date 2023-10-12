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
sbatch scripts/fastp_1.sh # 5677308 ok
sbatch scripts/fastp_2.sh # 5677320 ok
sbatch scripts/fastp_3.sh # 5677323 ok
```


# FastQC

**Raw:**
```bash
sbatch scripts/fastqc_1.sh # 5677297 ok
sbatch scripts/fastqc_2.sh # 5677298 ok
sbatch scripts/fastqc_3.sh # 5677299 ok
```

--> some raw have over-represented sequence; otherwise mostly high quality >30

**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:5677308 scripts/fastqc_fastp_1.sh # 5678160 ok
sbatch --dependency=afterany:5677320 scripts/fastqc_fastp_2.sh # 5678206 ok
sbatch --dependency=afterany:5677323 scripts/fastqc_fastp_3.sh # 5678207 ok
```

Some weird GC content (+ overrepresented seq >1.5%):
- NPC_KO_H3K4me3 
- NPC_KO_H3K27me1 
- NPC_KO_IGG 
- NPC_WT_H3K4me3 
- NPC_WT_H3K27me3
- NPC_WT_IGG
- PSC_KOEF1aEZH1_H3K27me3
- PSC_KOEF1aEZH1_IGG
- PSC_KOsynEZH1_H3K27me3
- PSC_KOsynEZH1_HA
- PSC_KOsynEZH1_IGG

Maybe due to the CutRun spikein??

--> Overall, it is ok!!


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch scripts/bowtie2_1.sh # 5680487 ok
sbatch scripts/bowtie2_2.sh # 5680490 ok
sbatch scripts/bowtie2_3.sh # 5680491 ok
```


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-5680487.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_5680487.txt

for file in slurm-5680490.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_5680490.txt

for file in slurm-5680491.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_5680491.txt
```

Add these values to `/home/roulet/001_EZH1_project/005__CutRun_NPC_PSC/sample.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >75% input reads as been uniquely mapped to the genome


# Calculate histone content

--> This histone content will be used to generate a scaling factor which will be used to histone-scaled our library size. The calcul/method to follow is from `003__CutRun/output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt`

**Pipeline:**
- Count the histone barcode on the clean reads
- Calculate SF (group by sample (replicate) and AB and calculate the total nb of reads. Then proportion of reads = nb read in sample / total reads. SF = min(proportion) / sample proportion)


### Count the histone barcode on the clean reads


```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_1.sh # 5695005 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_2.sh # 5695008 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_3.sh # 5695009 ok
```

--> It output the nb of reads found for each histone; then simply copy paste to the excell file `output/spikein/SpikeIn_QC_fastp.xlsx` in GoogleDrive


### Calculate SF
Calculate the scaling factor as with the raw reads from `output/spikein/SpikeIn_QC_fastp.xlsx`

--> Separate the 2 experiments (tissue NPC and tissue PSC)

```R
# package
library("tidyverse")
library("readxl")
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp.xlsx") 

# PSC
## H3K27me3
spikein_PSC_H3K27me3 = spikein %>%
    filter(tissue == "PSC",
           Target == "H3K27me3") %>%
    group_by(sample_ID, AB) %>%
    summarise(aligned=sum(counts))
# Total reads per IP
spikein_PSC_H3K27me3_total = spikein_PSC_H3K27me3 %>%
    ungroup() %>%
    group_by(AB) %>%
    mutate(total = sum(aligned)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_PSC_H3K27me3_read_prop = spikein_PSC_H3K27me3 %>%
    left_join(spikein_PSC_H3K27me3_total) %>%
    mutate(read_prop = aligned / total)
spikein_PSC_H3K27me3_read_prop_min = spikein_PSC_H3K27me3_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_PSC_H3K27me3_scaling_factor = spikein_PSC_H3K27me3_read_prop %>%
    left_join(spikein_PSC_H3K27me3_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_PSC_H3K27me3_scaling_factor, file="output/spikein/spikein_histone_PSC_H3K27me3_scaling_factor_fastp.txt", sep="\t", quote=FALSE, row.names=FALSE)



# NPC
## H3K27me3
spikein_NPC_H3K27me3 = spikein %>%
    filter(tissue == "NPC",
           Target == "H3K27me3") %>%
    group_by(sample_ID, AB) %>%
    summarise(aligned=sum(counts))
# Total reads per IP
spikein_NPC_H3K27me3_total = spikein_NPC_H3K27me3 %>%
    ungroup() %>%
    group_by(AB) %>%
    mutate(total = sum(aligned)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)
# Read proportion
spikein_NPC_H3K27me3_read_prop = spikein_NPC_H3K27me3 %>%
    left_join(spikein_NPC_H3K27me3_total) %>%
    mutate(read_prop = aligned / total)
spikein_NPC_H3K27me3_read_prop_min = spikein_NPC_H3K27me3_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_NPC_H3K27me3_scaling_factor = spikein_NPC_H3K27me3_read_prop %>%
    left_join(spikein_NPC_H3K27me3_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_NPC_H3K27me3_scaling_factor, file="output/spikein/spikein_histone_NPC_H3K27me3_scaling_factor_fastp.txt", sep="\t", quote=FALSE, row.names=FALSE)

## H3K4me3
spikein_corr_H3K4me3_H3K27me1 <- read_excel("output/spikein/SpikeIn_QC_fastp_corr_H3K4me3_H3K27me1.xlsx") %>%
    drop_na()
spikein_NPC_H3K4me3 = spikein_corr_H3K4me3_H3K27me1 %>%
    filter(tissue == "NPC",
           Target == "H3K4me3") %>%
    group_by(new_sample_ID, new_AB) %>%
    summarise(aligned=sum(counts)) %>%
    drop_na() # na because of the other AB used
# Total reads per IP
spikein_NPC_H3K4me3_total = spikein_NPC_H3K4me3 %>%
    ungroup() %>%
    group_by(new_AB) %>%
    mutate(total = sum(aligned)) %>%
    ungroup() %>%
    distinct(new_AB, .keep_all = TRUE) %>%
    select(new_AB,total)
# Read proportion
spikein_NPC_H3K4me3_read_prop = spikein_NPC_H3K4me3 %>%
    left_join(spikein_NPC_H3K4me3_total) %>%
    mutate(read_prop = aligned / total)
spikein_NPC_H3K4me3_read_prop_min = spikein_NPC_H3K4me3_read_prop %>%
    group_by(new_AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_NPC_H3K4me3_scaling_factor = spikein_NPC_H3K4me3_read_prop %>%
    left_join(spikein_NPC_H3K4me3_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_NPC_H3K4me3_scaling_factor, file="output/spikein/spikein_histone_NPC_H3K4me3_scaling_factor_fastp.txt", sep="\t", quote=FALSE, row.names=FALSE)

## H3K27me1
spikein_corr_H3K4me3_H3K27me1 <- read_excel("output/spikein/SpikeIn_QC_fastp_corr_H3K4me3_H3K27me1.xlsx") %>%
    drop_na()

spikein_NPC_H3K27me1= spikein_corr_H3K4me3_H3K27me1 %>%
    filter(tissue == "NPC",
           Target == "H3K27me1") %>%
    group_by(new_sample_ID, new_AB) %>%
    summarise(aligned=sum(counts)) %>%
    drop_na() # na because of the other AB used
# Total reads per IP
spikein_NPC_H3K27me1_total = spikein_NPC_H3K27me1 %>%
    ungroup() %>%
    group_by(new_AB) %>%
    mutate(total = sum(aligned)) %>%
    ungroup() %>%
    distinct(new_AB, .keep_all = TRUE) %>%
    select(new_AB,total)
# Read proportion
spikein_NPC_H3K27me1_read_prop = spikein_NPC_H3K27me1 %>%
    left_join(spikein_NPC_H3K27me1_total) %>%
    mutate(read_prop = aligned / total)
spikein_NPC_H3K27me1_read_prop_min = spikein_NPC_H3K27me1_read_prop %>%
    group_by(new_AB) %>%
    summarise(min_prop=min(read_prop))
# Scaling factor
spikein_NPC_H3K27me1_scaling_factor = spikein_NPC_H3K27me1_read_prop %>%
    left_join(spikein_NPC_H3K27me1_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)
write.table(spikein_NPC_H3K27me1_scaling_factor, file="output/spikein/spikein_histone_NPC_H3K27me1_scaling_factor_fastp.txt", sep="\t", quote=FALSE, row.names=FALSE)



```
--> H3K27me1 seems it did not work; even though histone enrichment...

--> Analysis of histone content showed that there was a sample inversion. **NPC_KO_H3K27me1 = NPC_KO_H3K4me3**. This HAS BEEN TAKEN into account using a separate input file `output/spikein/SpikeIn_QC_fastp_corr_H3K4me3_H3K27me1.xlsx`
----> This file is GOOD and is the SF for NPC_H3K4me3; WT `1` and KO `1.09261114977754`


--> All scaling factors combined in `output/spikein/spikein_histone_scaling_factor_fastp.txt` 

***IMPORTANT NOTE***: The sample inversion is a mess. Let's instead, rename the bam file correctly; and add "nameCorr" at the end; and re-process EVERYTHING!!! to make sure I did not make weird stuff... I only make a new file for `NPC_KO_H3K4me3`


## Quality control plot

Then look at the xlsx file from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) to generate quality control plot. Use R cluster for vizualization (file is `spikein_QC.xlsx` in Google Drive), file in `output/spikein`.
```R
# package
library("tidyverse")
library("readxl")
# import df adn tidy to remove AB used in sample_ID
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp.xlsx") %>%
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
    filter(sample_ID %in% c("NPC_WT", "NPC_KO", "PSC_KOEF1aEZH1", "PSC_KOsynEZH1"),
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
  # Find the target_norm value when Target is H3K27me3 and AB is H3K27me3
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
pdf("output/spikein/QC_histone_spike_in_H3K27me1.pdf", width = 7, height = 4)
spikein_all_scale %>%
    filter(sample_ID %in% c("NPC_WT", "NPC_KO"),
           AB %in% c("H3K27me1", "IGG")) %>%
        ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
        geom_col(position = "dodge") +
        facet_wrap(~sample_ID, nrow=1, scale ='free') +
        geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()



## Histone scaling for H3K4me3
spikein_all_scale = spikein_all %>%
  group_by(sample_ID) %>%
  # Find the target_norm value when Target is H3K27me3 and AB is H3K27me3
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
pdf("output/spikein/QC_histone_spike_in_H3K4me3.pdf", width = 7, height = 4)
spikein_all_scale %>%
    filter(sample_ID %in% c("NPC_WT", "NPC_KO"),
           AB %in% c("H3K4me3", "IGG")) %>%
        ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
        geom_col(position = "dodge") +
        facet_wrap(~sample_ID, nrow=1, scale ='free') +
        geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


```



# Samtools and read filtering

--> See `METHOD GOOD TO FOLLOW` in `003__CutRun` labnote

## Marking dupplicates
```bash
conda activate bowtie2

sbatch scripts/samtools_1.sh # 5694290 ok
sbatch scripts/samtools_2.sh # 5694293 ok
sbatch scripts/samtools_3.sh # 5694297 ok
```

## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.

```bash
conda activate bowtie2

sbatch --dependency=afterany:5694290 scripts/samtools_unique_1.sh # 5695055 ok
sbatch --dependency=afterany:5694293 scripts/samtools_unique_2.sh # 5695056 ok
sbatch --dependency=afterany:5694297 scripts/samtools_unique_3.sh # 5695057 ok
```


# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch scripts/bamtobigwig_1.sh # 5696658 time limit; 
sbatch scripts/bamtobigwig_2.sh # 5696659 time limit; 
sbatch scripts/bamtobigwig_3.sh # 5696660 ok

sbatch scripts/bamtobigwig_missing_1.sh # 5712700 ok
sbatch scripts/bamtobigwig_missing_2.sh # 5712702 ok

```

--> All good.



## Histone spike in normalized bigiwg

That's for histone PTM samples only


XXX Need MACS2 output !!





## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_PSC.sh # 5780604 ok
sbatch scripts/multiBigwigSummary_PSC_subset_1.sh # 5780758
sbatch scripts/multiBigwigSummary_NPC.sh # 5780645 ok



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
```


# MACS2 peak calling on bam unique

--> IGG samples used as control

--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

***PEAK CALLING:***
--> **H3K27me3, SUZ12, EZH2, EZH1, H3K4me3 in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_raw_PSC.sh # 5784629 ok
sbatch scripts/macs2_raw_NPC.sh # 5784690 ok
sbatch scripts/macs2_raw_NPC_H3K4me3.sh # 5784847
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
- NPC_H3K27me3 WT and KO; 2.30103
- NPC_SUZ12 WT and KO; 1.30103
- NPC_EZH2 WT and KO; 1.30103
- NPC_H3K4me3 WT and KO; 1.30103 
- PSC_EZH1cs KOEF1aEZH1; 1.30103 
- PSC_SUZ12 EF1aEZH1 and synEZH1; 1.30103 
- PSC_H3K27me3 EF1aEZH1 and synEZH1; 2.30103


--> No peak called for H3K27me1 samples



# Ecoli scaling factor
## Mapping E coli
- Map our reads to the E. coli genome using same parameters as for human.
- Count the number of aligned reads to the spike-in control sequences for each sample `samtools view -S -F 4 -c sample_h3k27me3_rep1_spikein.sam > sample_h3k27me3_rep1_spikein_count.txt`
- Do the math for scaling factor, same method as when using histone spike-in



```bash
conda activate bowtie2

# Map and count all samples to the E Coli genome
sbatch scripts/bowtie2_spike_Ecoli_1.sh # 5796988 ok
sbatch scripts/bowtie2_spike_Ecoli_2.sh # 5796989 ok
sbatch scripts/bowtie2_spike_Ecoli_3.sh # 5796996 ok

```
--> There is some uniq mapped reads, around 10% (in `003__CutRun` was less than 1%)

Now calculate SF in R, as for histone SF:

--> SAMPLE INVERISON: *I put the 614496 NPC_KO_H3K27me1 value for NPC_KO_H3K4me3; and opposite; **so the output table is ALL GOOD no changes to make!!!**

```R
# package
library("tidyverse")
library("readxl")
library("ggpubr")

# import df
spikein <- read_excel("output/spikein/SpikeIn_MG1655.xlsx") 

## NPC
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


## PSC
spikein <- read_table("output/spikein/SpikeIn_MG1655.txt") %>%
    filter(tissue == "PSC") %>%
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
write.table(spikein_scaling_factor, file="output/spikein/spikein_MG1655_PSC_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)
```






## Ecoli/Exogeneous genome
```bash
sbatch scripts/samtools_MG1655_corr.sh # XXX
```



# Spike in scaling
## With histone spike in for PTM CutRun

**Using our scaling factor, let's estimate the 'new' library size** and provide it to `dba.normalize(library = c(1000, 12000))` = Like that our library size will be change taking into account our scaling factor! **Then we can normalize with library-size, RLE or TMM**... (issue discussed [here](https://support.bioconductor.org/p/9147040/)) 


### Adjust library size with histone scaling factor and apply normalization
Total number of reads is our library size (used samtools flagstat to double check) :

`samtools flagstat output/bowtie2/*.dupmark.sorted.bam` used to obtain library size (first value=library size)
--> Values save in GoogleDrive `005__*/sample.xlsx`. Histone-norm-library-size = library-size * SF. Using the non-reciprocal scaling factor, we increase the library-size; the more histone enriched, the more library size is increased, thus the more signal will decrease.

Now let's use these new histone-scaled library size and normalize with library-size,TMM or RLE. Let's use the **unique bam files** together with the **unique bam MACS2 raw files (xlsx, not the bed with pre-filtered qvalue)**

***Key points:***
- **Let's do 1 DiffBind per AB (H3K27me3, H3K4me3,...) and tissue (PSC, NPC); otherwise the TMM normalization may take all, unrelated, samples into account!** --> Files are `meta_sample_macs2raw_unique*.txt`
- **For the non-histone CutRun, I will use the library size non histone scaled in DiffBind to collect TMM normalized SF**; I tested with and without specifying library size; and it does not change a lot the SF... Let's better use the one RiP method w/o providing the library size! Should provide BETTER correction
- Not sure the **synEZH1 should be treated as an additional sample; it's a negative control like IGG; we should NOT use them to scale the normalization, notably for EZH1**; for SUZ12, H3K27me3 I think it is OK!

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

sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(8014384,10146334), normalize = DBA_NORM_TMM) 

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
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(5606820,9124745), normalize = DBA_NORM_TMM) # SAME HERE THE LIB for KO = the one from H3K27me1; 9124745  and not 8281344 or 12401403!!!
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_H3K4me3.txt")

# NPC_SUZ12
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_SUZ12.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_SUZ12.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_SUZ12.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(5681310,7250570), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_SUZ12.txt")


# NPC_EZH2
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_EZH2.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_EZH2.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_EZH2.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(6458862,7241732), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_NPC_EZH2.txt")



# PSC_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_H3K27me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_H3K27me3.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(9881249,8556476), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC_H3K27me3.txt")


# PSC_SUZ12
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_SUZ12.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_SUZ12.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_SUZ12.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(7595642,6895662), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC_SUZ12.txt")



# PSC_EZH1cs
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_EZH1cs.txt", header = TRUE, sep = "\t")) # count fail as no peak for synEZH1...
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_EZH1cs_EF1a.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba) 
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_EZH1cs.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_EZH1cs.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(5978472,7011234), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC_EZH1cs.txt")



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
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(8014384,10146334,5606820,9124745,5681310,7250570,6458862,7241732), normalize = DBA_NORM_TMM) # HERE sample inversion I make sure i use H3K27me1 for KO H3K4me3!
#### Here is to retrieve the scaling factor value
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


## PSC
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC.RData")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_FACTOR, label=DBA_TREATMENT)
dev.off()
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_blackgreylist.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_blackgreylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_FACTOR, label=DBA_TREATMENT)
dev.off()
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(9881249, 8556476,7595642,6895662, 5978472), normalize = DBA_NORM_TMM) # HERE sample inversion I make sure i use H3K27me1 for KO H3K4me3!
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC.txt")
### plot
pdf("output/DiffBind/clustering_sample_macs2raw_unique_PSC_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)  
plot(sample_count)
dev.off()
pdf("output/DiffBind/PCA_sample_macs2raw_unique_PSC_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_FACTOR, label=DBA_TREATMENT)
dev.off()
```

--> Cannot retrieve SF for PSC_EZH1cs because nothing to be compare with; as no peak detected in synEZH1, it cannot normalize the signal...

--> **sample inversion for NPC_H3K4me3** for this I modified accordingly the KO samples; which is NPC_KO_H3K4me3 = NPC_KO_H3K27me1. Modified in the meta file...
DETAIL MATH:
Then with this; calculate the histone-norm-library-size:
- library size: `samtools flagstat output/bowtie2/NPC_KO_H3K4me3_nameCorr.unique.dupmark.sorted.bam` = 8351320
- So correct histone-scale library size = 9124745
- Then correct SF:  0.803424 and 1.173694 for WT and KO


Now let's do the same method but using the **MG1655_library_scaled information and collect new MG1655_DiffBind_TMM_SF**:


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
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(12342151,9417384), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_NPC_H3K27me3.txt")


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
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(7064593,8351320), normalize = DBA_NORM_TMM) # SAME HERE THE LIB for KO = I USED THE GOOD ONE!
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_NPC_H3K4me3.txt")

# NPC_SUZ12
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_SUZ12.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_SUZ12.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_SUZ12.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(12726134,7250570), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_NPC_SUZ12.txt")


# NPC_EZH2
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_EZH2.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_EZH2.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_EZH2.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(8138166,7241732), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_NPC_EZH2.txt")



# PSC_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_H3K27me3.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_H3K27me3.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_H3K27me3.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(20925063,8556476), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LiMG1655Scaled_TMM_unique_SF_PSC_H3K27me3.txt")


# PSC_SUZ12
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_SUZ12.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_SUZ12.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_SUZ12.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(7595642,13308628), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibMG1655Scaled_TMM_unique_SF_PSC_SUZ12.txt")



# PSC_EZH1cs
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_EZH1cs.txt", header = TRUE, sep = "\t")) # count fail as no peak for synEZH1...
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_PSC_EZH1cs_EF1a.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba) 
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_PSC_EZH1cs.RData")
load("output/DiffBind/sample_count_macs2raw_unique_PSC_EZH1cs.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### TMM 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_TMM) 
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(5978472,7011234), normalize = DBA_NORM_TMM) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_unique_SF_PSC_EZH1cs.txt")

```

--> Cannot generate a MG1655_DiffBind_TMM SF for PSC_EZH1cs as no peak in the syn condition; so I put 1 for both...

--> The SF collected are NOT the same as the one when I used the histone-spike in methiod. However; I trust more the E coli method. As when there is more E coli DNA then the SF is more than 1 thus increasing the library size, thus decreasing the signal as there is just more DNA!
----> Interestingly, when doing the TMM normalization; it put very close togethe the SF generated from histone and from MG1655 method... Not sure, what to conclude....

## Generate Spike in scaled bigwig

--> Reciprocal from DiffBind_TMM is to be used when converting bam to bigwig!

### Histone scaled Bigwig


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_histone_DiffBind_TMM_1.sh # 5811580 ok
sbatch scripts/bamtobigwig_histone_DiffBind_TMM_2.sh # 5811582 ok
sbatch scripts/bamtobigwig_histone_DiffBind_TMM_3.sh # 5811584 ok
```



### MG1655/E coli scaled bigwig


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_MG1655_DiffBind_TMM_1.sh # 5812907 ok
sbatch scripts/bamtobigwig_MG1655_DiffBind_TMM_2.sh # 5812908 ok
sbatch scripts/bamtobigwig_MG1655_DiffBind_TMM_3.sh # 5812915 ok
```


--> Both bigwig norm method are very similar... 

Check some known target regulated in 2months neurons:
--> NEUROG2 seems less in KO which is good.
--> EFNA5 tiny decrease in KO (only in normalized data!)
--> GRIK3 tiny increase in KO

--> Something is WEIRD... When samples have MORE spike in, their signal should be reduced, as they overall have more DNA; but if I used the reciprocal from DiffBind_TMM; this is not respected (ie. sample with more spike in, we increased their signal...!)... That is true for both histone/MG1655-spike in DiffBind TMM norm...
----> What should be the BEST to use, is then the NON-reciprocal_DiffBind_TMM !!!
------> Let's try and compare with gene expression...! Maybe it is still good as we take into account the library size with the DiffBind_TMM method??
--------> YESSS in the end we correct the library size with the SF!!!!! So we 're good!!! reciprocal DiffBind_TMM IS TO BE USED!!



# Deeptools plot on DEGs


Let's use our RNAseq in NPC, and generate heatmap of H3K27me3 CutRun signal for up and down; to confirm the normalization works great. DEGs gtf already generated (`../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_KO_Up.gtf AND ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_KO_Down.gtf`)

```bash
conda activate deeptools

sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig.sh # 5827366 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig.sh # 5827369 ok

sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig_histone.sh # 5835369
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_histone.sh # 5835370

sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig_MG1655.sh # 5835363
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_MG1655.sh # 5835367


```

--> The **raw bigwig is not good**; non normalized file is NOT in agreement with gene expression; H3K27me3 signal always higher in the NPC_KO





# ChIPseeker peak gene assignment

## From optimal qval bed files peaks
Let's assign peak to genes from MACS2 peak:

**Optimal qvalue** according to IGV:
- NPC_H3K27me3 WT and KO; 2.30103
- NPC_SUZ12 WT and KO; 1.30103
- NPC_EZH2 WT and KO; 1.30103
- NPC_H3K4me3 WT and KO; 1.30103 
- PSC_EZH1cs KOEF1aEZH1; 1.30103 
- PSC_SUZ12 EF1aEZH1 and synEZH1; 1.30103 
- PSC_H3K27me3 EF1aEZH1 and synEZH1; 2.30103


XXX





