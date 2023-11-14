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

Let's do the same for E coli MG1655 spike in samples:

```bash
conda activate bowtie2

sbatch scripts/samtools_MG1655_unique_1.sh # 6544133 ok
sbatch scripts/samtools_MG1655_unique_2.sh # 6544173 ok
sbatch scripts/samtools_MG1655_unique_3.sh # 6544194 ok

# Create index for the genome to be load in IGV
samtools faidx ../003__CutRun/meta/MG1655_v2.fa

```

- *NOTE: to know what are the contig/chromosome names of MG1655; use `samtools idxstats output/spikein/*.sam`*
- *NOTE: The MG1655 fasta genome got an empty line at row 452 and many others... WTF, remove it and save the new fasta as 003__CutRun/meta/MG1655_v2.fa*

--> E coli MG1655 alignemnt seems to have work can see reads on IGV !!

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


sbatch scripts/bamtobigwig_unique_1.sh # 7022154
sbatch scripts/bamtobigwig_unique_2.sh # 7022155
sbatch scripts/bamtobigwig_unique_3.sh # 7022156
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
sbatch scripts/multiBigwigSummary_NPC_THOR.sh # 5931653 ok


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
plotPCA -in output/THOR/multiBigwigSummary_NPC_THOR.npz \
    --transpose \
    --ntop 0 \
    --labels WT_SUZ12 KO_SUZ12 WT_EZH2 KO_EZH2 WT_H3K27me3 KO_H3K27me3 WT_H3K4me3 KO_H3K4me3 \
    --colors blue blue green green red red gold gold \
    --markers o s o s o s o s \
    --plotWidth 7 \
    -o output/THOR/multiBigwigSummary_NPC_THOR_plotPCA.pdf


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
plotCorrelation \
    -in output/THOR/multiBigwigSummary_NPC_THOR.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_SUZ12 KO_SUZ12 WT_EZH2 KO_EZH2 WT_H3K27me3 KO_H3K27me3 WT_H3K4me3 KO_H3K4me3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/THOR/multiBigwigSummary_NPC_THOR_heatmap.pdf
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
sbatch scripts/samtools_MG1655_corr.sh # ok
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
----> Interestingly, when doing the TMM normalization; it put very close togethe the SF generated from histone and from MG1655 method... Not sure, what to conclude.... That looks weird!!


Lets try to use the **E coli bam spike in scaling method** (add MG1655 bam files in the meta file):
--> For NPC samples only
--> Vignete 7.6 from [DiffBind doc](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) to check how apply spikein

**spikein_LIB normalization**

```bash
conda activate DiffBind
```

```R
library("DiffBind") 


# ONE PER ONE
## NPC_H3K27me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_H3K27me3_MG1655bam.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3_MG1655bam.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3_MG1655bam.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### Spike in normalization 
sample_count_blackgreylist_LIB_spikein = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_LIB, spikein=TRUE) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LIB_spikein_SF = dba.normalize(sample_count_blackgreylist_LIB_spikein, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LIB_spikein_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LIB_spikein_SF_NPC_H3K27me3.txt")




# NPC_H3K4me3
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_H3K4me3_MG1655bam.txt", header = TRUE, sep = "\t")) # HERE I MADE SURE USING THE CORRECT SAMPLE FILE (sampel inversion!)
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_H3K4me3_MG1655bam.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K4me3_MG1655bam.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist 



sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### Spike in normalization 
sample_count_blackgreylist_LIB_spikein = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_LIB, spikein=TRUE) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LIB_spikein_SF = dba.normalize(sample_count_blackgreylist_LIB_spikein, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LIB_spikein_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LIB_spikein_SF_NPC_H3K4me3.txt")






# NPC_SUZ12
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_SUZ12_MG1655bam.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_SUZ12_MG1655bam.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_SUZ12_MG1655bam.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### Spike in normalization 
sample_count_blackgreylist_LIB_spikein = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_LIB, spikein=TRUE) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LIB_spikein_SF = dba.normalize(sample_count_blackgreylist_LIB_spikein, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LIB_spikein_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LIB_spikein_SF_NPC_SUZ12.txt")


# NPC_EZH2
### Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw_unique_NPC_EZH2_MG1655bam.txt", header = TRUE, sep = "\t"))
### Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw_unique_NPC_EZH2_MG1655bam.RData")
load("output/DiffBind/sample_count_macs2raw_unique_NPC_EZH2_MG1655bam.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### Spike in normalization 
sample_count_blackgreylist_LIB_spikein = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_LIB, spikein=TRUE) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_LIB_spikein_SF = dba.normalize(sample_count_blackgreylist_LIB_spikein, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_LIB_spikein_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LIB_spikein_SF_NPC_EZH2.txt")

```


**spikein_TMM and RLE normalization**

Some info [here](https://support.bioconductor.org/p/9148431/#9148434)

```R
library("DiffBind") 


# ONE PER ONE
## NPC_H3K27me3
load("output/DiffBind/sample_count_macs2raw_unique_NPC_H3K27me3_MG1655bam.RData")
### Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
### Spike in normalization        _ TMM
sample_count_blackgreylist_TMM_spikein = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_TMM, spikein=TRUE) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_TMM_spikein_SF = dba.normalize(sample_count_blackgreylist_TMM_spikein, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_TMM_spikein_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_TMM_spikein_SF_NPC_H3K27me3.txt")
### Spike in normalization        _ RLE       
sample_count_blackgreylist_RLE_spikein = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_RLE, spikein=TRUE) 
#### Here is to retrieve the scaling factor value
sample_count_blackgreylist_RLE_spikein_SF = dba.normalize(sample_count_blackgreylist_RLE_spikein, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_blackgreylist_RLE_spikein_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_RLE_spikein_SF_NPC_H3K27me3.txt")

```

--> SF are unchanged compared to the LIB method...





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

# NPC
## H3K27me3
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig.sh # 5827366 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig.sh # 5827369 ok

sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig_histone.sh # 5835369 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_histone.sh # 5835370 ok

sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig_MG1655.sh # 5835363 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_MG1655.sh # 5835367 ok

sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig_THOR.sh # 5861199 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_THOR.sh # 5861207 ok

sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig_LIB_spikein_THOR.sh # 6619767 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_LIB_spikein_THOR.sh # 6619963 ok

sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig_LIB_spikein_NotReciprocal_THOR.sh # 6664317 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K27me3_bigwig_LIB_spikein_NotReciprocal_THOR.sh # 6664367 ok


## SUZ12
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_SUZ12_bigwig_THOR.sh # 5932584 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_SUZ12_bigwig_THOR.sh # 5933018 ok
## EZH2
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_EZH2_bigwig_THOR.sh # 5933273 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_EZH2_bigwig_THOR.sh # 5933517 ok
## H3K4me3
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K4me3_bigwig_THOR.sh # 5933811 ok
sbatch scripts/matrix_TSS_10kb_DEGs_NPC_KO_Down_H3K4me3_bigwig_THOR.sh # 5933831 ok

# PSC
sbatch scripts/matrix_TSS_10kb_PSC_KOEF1aEZH1_bigwig_1.sh # 5850570 ok
sbatch scripts/matrix_TSS_10kb_PSC_KOEF1aEZH1_bigwig_2.sh # 5850574 ok
```
***NPC:***
--> The **raw bigwig is not good**; non normalized file is NOT in agreement with gene expression; H3K27me3 signal always higher in the NPC_KO
--> The **bigwig histone and MG1655 is OK**; slight differences in agreement with gene expression changes for H3K27me3. In CANNOT be perfect as the RNAseq do not match with this experiment; diff. protocol used not the same! **Let's prefer to use MG1655=same normalziation method for ALL samples!**
--> The **bigwig THOR is BETTER THAN MG1655**; clearly see increase of H3K27me3 for genes down; but no changes for gene Up.
--> The **matrix_TSS_10kb_DEGs_NPC_KO_Up_H3K27me3_bigwig_LIB_spikein_THOR** is SHIT. The DEGs H3K27me3 is NOT in agreement with profile; KO always higher...


***PSC:***
For AB comparions, Need to play with  `--zMin -2 -2 0 --zMax 2 2 3` in plotHeatmap to scale the heatmap individually! (`--perGroup` do not help; it just put them one under the other)

--> I used the raw bigwig, as different marks are compared.

--> EZH1cs and SUZ12 co-localize well

--> Order_2 (EZH1cs/SUZ12/H3K27me3) is better to assess the EZH1/H3K27me3 co-localization

Let's now compare the mark within and between genotypes; check signal genome-wide:

```bash
conda activate deeptools

# Within genotpes
sbatch scripts/matrix_TSS_10kb_NPC_WT_bigwig_MG1655.sh # 5849559 ok
sbatch scripts/matrix_TSS_10kb_NPC_KO_bigwig_MG1655.sh # 5849564 ok

# Between genotypes
sbatch scripts/matrix_TSS_10kb_NPC_H3K27me3_bigwig_MG1655.sh # 5849579 ok
sbatch scripts/matrix_TSS_10kb_NPC_EZH2_bigwig_MG1655.sh # 5849582 ok
sbatch scripts/matrix_TSS_10kb_NPC_SUZ12_bigwig_MG1655.sh # 5849585 ok
sbatch scripts/matrix_TSS_10kb_NPC_H3K4me3_bigwig_MG1655.sh # 5849589 ok

sbatch scripts/matrix_TSS_10kb_NPC_H3K27me3_bigwig_THOR.sh # 5861295 ok
sbatch scripts/matrix_TSS_10kb_NPC_EZH2_bigwig_THOR.sh # 5861354 ok
sbatch scripts/matrix_TSS_10kb_NPC_SUZ12_bigwig_THOR.sh # 5861385 ok
sbatch scripts/matrix_TSS_10kb_NPC_H3K4me3_bigwig_THOR.sh # 5861390 ok

sbatch scripts/matrix_TSS_10kb_NPC_EZH2_bigwig_raw.sh # 6472039 ok
sbatch scripts/matrix_TSS_10kb_NPC_SUZ12_bigwig_raw.sh # 6472070 ok
sbatch scripts/matrix_TSS_10kb_NPC_H3K27me3_bigwig_raw.sh # 6472080 ok
sbatch scripts/matrix_TSS_10kb_NPC_H3K4me3_bigwig_raw.sh # 6472160 ok

sbatch --dependency=afterany:6551938 scripts/matrix_TSS_10kb_NPC_EZH2_bigwig_LIB_spikein.sh # 6556104 ok
sbatch scripts/matrix_TSS_10kb_NPC_SUZ12_bigwig_LIB_spikein.sh # 6556141 ok
sbatch scripts/matrix_TSS_10kb_NPC_H3K27me3_bigwig_LIB_spikein.sh # 6556167 ok
sbatch scripts/matrix_TSS_10kb_NPC_H3K4me3_bigwig_LIB_spikein.sh # 6556182 ok

# Compare effect of crosslinking (native vs FA 005vs006)
sbatch --dependency=afterany:7022154:7022155:7022156 scripts/matrix_TSS_10kb_PSC_H3K27me3_bigwig_unique_005vs006.sh # 7022288 XXXX



```

--> For  **Within genotpes, I used the MG1655 version; but could have use the raw bigwig as I compare diff marks...**
- Nice co-localization centered around TSS of SUZ12, EZH2, H3K27me3

--> For **Between genotypes, important to use the MG1655 version however!**
- Genome-wide; signal is always higher for WT; except H3K27me3, same enrichment...
- EZH2 and SUZ12 are weird... The signal in KO is much lower than WT, we even see it on IGV! I suspect THOR normalization to account for IGG is NEEDED. Interestingly, much more IGG signal in KO than WT...
- **Bigwig THOR is PERFECT**; a btit higher H3K27me3 upon KO; H3K4me3 is not affected, SUZ12 and EZH2 decrease.


# ChIPseeker peak gene assignment

## From optimal qval bed files peaks
Let's assign **peak to genes from MACS2 peak**:

**Optimal qvalue** according to IGV:
- NPC_H3K27me3 WT and KO; 2.30103
- NPC_SUZ12 WT and KO; 1.30103
- NPC_EZH2 WT and KO; 1.30103
- NPC_H3K4me3 WT and KO; 1.30103 
- PSC_EZH1cs KOEF1aEZH1; 1.30103 
- PSC_SUZ12 EF1aEZH1 and synEZH1; 1.30103 
- PSC_H3K27me3 EF1aEZH1 and synEZH1; 2.30103

--> Assign peak to genes for NPC and PSC:
--> *NPC*: After assigning peak to genes; we will isolate putative EZH1 bound genes (SUZ12 without EZH2)

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
## NPC
H3K27me3 = as_tibble(read.table('output/macs2/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
EZH2 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/NPC_WT_EZH2_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)    
SUZ12 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)    
H3K4me3 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/NPC_WT_H3K4me3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)  


## PSC
EZH1 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_EZH1cs_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)  
SUZ12 = as_tibble(read.table('output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_SUZ12_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)  
H3K27me3 = as_tibble(read.table('output/macs2/broad_blacklist_qval2.30103/PSC_KOEF1aEZH1_H3K27me3_peaks.broadPeak') ) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)  

# Tidy peaks #-->> Re-Run from here with different qvalue!!
## NPC
H3K27me3_gr = makeGRangesFromDataFrame(H3K27me3,keep.extra.columns=TRUE)
EZH2_gr = makeGRangesFromDataFrame(EZH2,keep.extra.columns=TRUE)
SUZ12_gr = makeGRangesFromDataFrame(SUZ12,keep.extra.columns=TRUE)
H3K4me3_gr = makeGRangesFromDataFrame(H3K4me3,keep.extra.columns=TRUE)
gr_list <- list(H3K27me3=H3K27me3_gr, EZH2=EZH2_gr,  SUZ12=SUZ12_gr, H3K4me3=H3K4me3_gr)

## PSC
EZH1_gr = makeGRangesFromDataFrame(EZH1,keep.extra.columns=TRUE)
SUZ12_gr = makeGRangesFromDataFrame(SUZ12,keep.extra.columns=TRUE)
H3K27me3_gr = makeGRangesFromDataFrame(H3K27me3,keep.extra.columns=TRUE)
gr_list <- list(EZH1=EZH1_gr, SUZ12=SUZ12_gr,  H3K27me3=H3K27me3_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
H3K27me3_annot <- as.data.frame(peakAnnoList[["H3K27me3"]]@anno)
EZH2_annot <- as.data.frame(peakAnnoList[["EZH2"]]@anno)
SUZ12_annot <- as.data.frame(peakAnnoList[["SUZ12"]]@anno)
H3K4me3_annot <- as.data.frame(peakAnnoList[["H3K4me3"]]@anno)
EZH1_annot <- as.data.frame(peakAnnoList[["EZH1"]]@anno)

## Convert entrez gene IDs to gene symbols
H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
EZH2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH2_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
SUZ12_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = SUZ12_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
SUZ12_annot$gene <- mapIds(org.Hs.eg.db, keys = SUZ12_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
EZH1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH1_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(H3K27me3_annot, file="output/ChIPseeker/annotation_macs2_WT_H3K27me3_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(EZH2_annot, file="output/ChIPseeker/annotation_macs2_WT_EZH2_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(SUZ12_annot, file="output/ChIPseeker/annotation_macs2_WT_SUZ12_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_annot, file="output/ChIPseeker/annotation_macs2_WT_H3K4me3_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

write.table(H3K27me3_annot, file="output/ChIPseeker/annotation_macs2_KOEF1aEZH1_H3K27me3_qval2.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(SUZ12_annot, file="output/ChIPseeker/annotation_macs2_KOEF1aEZH1_SUZ12_qval1.30103.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K27me3_annot_promoterAnd5 = tibble(H3K27me3_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
EZH2_annot_promoterAnd5 = tibble(EZH2_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
SUZ12_annot_promoterAnd5 = tibble(SUZ12_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
H3K4me3_annot_promoterAnd5 = tibble(H3K4me3_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
EZH1_annot_promoterAnd5 = tibble(EZH1_annot) %>%
    filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
### Save output gene lists
H3K27me3_annot_promoterAnd5_geneSymbol = H3K27me3_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
EZH2_annot_promoterAnd5_geneSymbol = EZH2_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
SUZ12_annot_promoterAnd5_geneSymbol = SUZ12_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_annot_promoterAnd5_geneSymbol = H3K4me3_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()
EZH1_annot_promoterAnd5_geneSymbol = EZH1_annot_promoterAnd5 %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(H3K27me3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_WT_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(EZH2_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_WT_EZH2_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(SUZ12_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_WT_SUZ12_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_WT_H3K4me3_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

write.table(H3K27me3_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_KOEF1aEZH1_H3K27me3_qval2.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(SUZ12_annot_promoterAnd5_geneSymbol, file = "output/ChIPseeker/annot_macs2_KOEF1aEZH1_SUZ12_qval1.30103_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
```
**NPC**:
--> Export gene list to Online Venn diagram to isolate EZH1 specific (SUZ12 non EZH2); Then deepTools plot:
- Generate gtf of EZH1; EZH2 specifc gene sets --> `output/macs2/Venn_SUZ12-[no]EZH2_geneSymbol.txt`
- QC of our EZH1 vs EZH2 specific WT THOR bigwig with SUZ12, EZH2, H3K27me3 for bed EZH2-spe and bed EZH1-spe
----> Make sure no EZH2 signal in the putative EZH1 one!
- Compare H3K27me3 level for genes EZH1 and EZH2 -specific

**PSC**:
--> Export gene list to Online Venn diagram EZH1, SUZ12, H3K27me3
--> Then do GO/KEGG on EZH1 to see if steroid-related stuff


```bash
# Generate gtf from gene Symbol list
perl -p -i -e 's/\r$//' output/macs2/Venn_SUZ12-EZH2_geneSymbol.txt  # THIS TO CONVERT windowns to UNIX; as .txt from windows...
perl -p -i -e 's/\r$//' output/macs2/Venn_SUZ12-noEZH2_geneSymbol.txt 
### Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/macs2/Venn_SUZ12-EZH2_geneSymbol.txt > output/macs2/Venn_SUZ12-EZH2_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/macs2/Venn_SUZ12-noEZH2_geneSymbol.txt > output/macs2/Venn_SUZ12-noEZH2_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annot_macs2_WT_SUZ12_qval1.30103_promoterAnd5_geneSymbol.txt > output/ChIPseeker/annot_macs2_WT_SUZ12_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt
### Filter the gtf
grep -Ff output/macs2/Venn_SUZ12-EZH2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_SUZ12-EZH2.gtf
grep -Ff output/macs2/Venn_SUZ12-noEZH2_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_SUZ12-noEZH2.gtf
grep -Ff output/ChIPseeker/annot_macs2_WT_SUZ12_qval1.30103_promoterAnd5_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_macs2_WT_SUZ12_qval1.30103_promoterAnd5.gtf
```

Now deepTools

```bash
conda activate deeptools

# heatmap
## NPC
sbatch scripts/matrix_TSS_10kb_THOR_SUZ12-EZH2_WT.sh # 5904250 ok

sbatch scripts/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes.sh # 5904410
sbatch scripts/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak.sh # 5904491 ok
sbatch scripts/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_genes.sh # 5904530 ok

## PSC
sbatch scripts/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak.sh # 5952851 ok

```
--> Let's try k-means method with all genes and genes with a SUZ12 peaks (from macs2)

**conclusion to highlight PRC2-EZH1**:
- `matrix_TSS_10kb_THOR_SUZ12-EZH2_WT`: Fail at separating SUZ12 no EZH2; in the end SUZ12 binding is reduced here; as EZH2...
- `matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes`: XXX
- `matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_genes`: Fail at separating SUZ12 no EZH2. Kmeans lead to signal and no signal for both marks.
- `matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak` and `matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak`: GREAT; with 25 clusters we are able to isolate peaks that are SUZ12 without EZH2 !!! And with SUZ12 without EZH1 (surprisingly very few!)

Now from `matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak` and `matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak`; let's isolate the peak SUZ12-EZH2 (cluster 1-24) and the SUZ12-noEZH2 (cluster 25)

Also isolate cluster 22-25 as PRC2-EZH2 instead of only the cl25

```bash
# Isolate specific cluster and make bed
output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_heatmap_kmeans25.bed
## matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak cluster 1-24
awk '$13 != "cluster_25"' output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_heatmap_kmeans25.bed > output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_heatmap_kmeans25_cluster1-24.bed
awk '$13 != "cluster_25"' output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25.bed > output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster1-24.bed

awk '$13 != "cluster_25" && $13 != "cluster_24" && $13 != "cluster_23" && $13 != "cluster_22"' output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25.bed > output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster1-21.bed


### Re-format the bed from deepTools kmeans as it is buggy
conda activate bowtie2
bedtools intersect -wa -a output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_heatmap_kmeans25_cluster1-24.bed -b output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks.broadPeak | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' > output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_heatmap_kmeans25_cluster1-24_macs2Format.bed
bedtools intersect -wa -a output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster1-24.bed -b output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_SUZ12_peaks.broadPeak | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' > output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster1-24_macs2Format.bed

bedtools intersect -wa -a output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster1-21.bed -b output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_SUZ12_peaks.broadPeak | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' > output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster1-21_macs2Format.bed
bedtools intersect -wa -a output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster22-25.bed -b output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_SUZ12_peaks.broadPeak | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' > output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster22-25_macs2Format.bed


## matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak cluster 25
awk '$13 == "cluster_25"' output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_heatmap_kmeans25.bed > output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_heatmap_kmeans25_cluster25.bed
awk '$13 == "cluster_25"' output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25.bed > output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster25.bed

awk '$13 == "cluster_25" || $13 == "cluster_24" || $13 == "cluster_23" || $13 == "cluster_22"' output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25.bed > output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster22-25.bed

```

Now we have *SUZ12-EZH2; SUZ12-noEZH2; let's check the H3K27me3 level in WT; for the peaks (not the gene)*:
```bash
conda activate deeptools

# heatmap
## NPC
sbatch scripts/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_cl1-24_25.sh # 5908058 ok
sbatch scripts/matrix_TSS_10kb_THOR_SUZ12_EZH2_WTandKO_SUZ12_macs2_broadPeak_cl1-24_25.sh # 5908857 ok
sbatch scripts/matrix_TSS_10kb_THOR_SUZ12_EZH2_WTandKO_SUZ12_macs2_broadPeak_cl25.sh # 5911498 ok

## PSC
sbatch scripts/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_cl1-24_25.sh # 5958046 ok
sbatch scripts/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_cl1-21_22-25.sh # 6668812 ok

```

--> **WT control**: SUZ12-noEZH2 show a reduced level of H3K27me3 as compared to SUZ12-EZH2. Consistent as EZH2 is known as a better histone methyltransferase?

--> **KO effect**: The SUZ12 and EZH2 reduction upon KO is only visible for the EZH1-sepcific peaks!! H3K27me3 does NOT seems to be affected by EZH1 KO















# THOR 

Let's use THOR, notably to have IGG scaled bigwig...!

Comparison to do; NPC WT vs KO:
- EZH2
- SUZ12
- H3K27me3
- H3K4me3

--> SF to use in THOR are the **reciprocal of MG1655_DiffBind_TMM**
--> Configs file created manually as `output/THOR/NPC_EZH2.config`

--> Lets also try to use the DiffBind spike in BAM method (similarly use the reciprocal from diffBind)




## Run THOR

*THOR is very buggy to make it work I need to temporaly change where to look for libraries lol.. So cannot use nano anymore for example...*

*Follow these parameters: `WTvsHET_unique_Keepdup` (perform best in previous CutRun)*

```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge

# AB per AB
sbatch scripts/THOR_NPC_EZH2.sh # 5855649 ok
sbatch scripts/THOR_NPC_SUZ12.sh # 5856888 ok
sbatch scripts/THOR_NPC_H3K27me3.sh # 5856898 ok
sbatch scripts/THOR_NPC_H3K4me3.sh # 5856928 ok

# DiffBind MG1655 bam scaling factor
## doin the math (SF dibbfinb recipocal * lib size)
sbatch scripts/THOR_NPC_EZH2_LIB_spikein.sh # 6551938 ok
sbatch scripts/THOR_NPC_SUZ12_LIB_spikein.sh # 6552579 ok
sbatch scripts/THOR_NPC_H3K27me3_LIB_spikein.sh # 6552666 ok
sbatch scripts/THOR_NPC_H3K4me3_LIB_spikein.sh # 6552713 ok

sbatch scripts/THOR_NPC_H3K27me3_LIB_spikein_NotReciprocal.sh # 6625785 ok
```

--> (Using Recirpocl DiffBind TMM initial method) By eye we seems to still see the higher EZH2 enrichment in WT / KO... 
----> Using the DiffBind MG1655 bam scaling factor looks good, and value are VERY DIFFERENT thatn the 1st method. Lets' see peak gene assgihnent to gene and expression...





## Filter THOR peaks (qvalue)

Let's find the optimal qvalue for THOR diff peaks



```R

# load the file using the tidyverse
library("readr")
library("dplyr")
library("ggplot2")
library("tidyr")

# H3K27me3
diffpeaks <- read_tsv("output/THOR/THOR_NPC_H3K27me3/NPCH3K27me3-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE)  %>%
  mutate(FC = (count_KO) / (count_WT))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_H3K27me3/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_NPC_H3K27me3/log2FC_qval20.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO_qval20") +
  theme_bw()
dev.off()

## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 100) %>%
  write_tsv("output/THOR/THOR_NPC_H3K27me3/THOR_qval100.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())





# EZH2
diffpeaks <- read_tsv("output/THOR/THOR_NPC_EZH2/NPCEZH2-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE)  %>%
  mutate(FC = (count_KO) / (count_WT))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_EZH2/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_NPC_EZH2/log2FC_qval10.pdf", width=14, height=14)
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
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_NPC_EZH2/THOR_qval30.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())





# SUZ12
diffpeaks <- read_tsv("output/THOR/THOR_NPC_SUZ12/NPCSUZ12-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE)  %>%
  mutate(FC = (count_KO) / (count_WT))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_SUZ12/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_NPC_SUZ12/log2FC_qval25.pdf", width=14, height=14)
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
  filter(qval > 10) %>%
  write_tsv("output/THOR/THOR_NPC_SUZ12/THOR_qval10.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())





# H3K4me3
diffpeaks <- read_tsv("output/THOR/THOR_NPC_H3K4me3/NPCH3K4me3-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE)  %>%
  mutate(FC = (count_KO) / (count_WT))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_H3K4me3/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()

pdf("output/THOR/THOR_NPC_H3K4me3/log2FC_qval25.pdf", width=14, height=14)
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
  filter(qval > 10) %>%
  write_tsv("output/THOR/THOR_NPC_H3K4me3/THOR_qval10.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())
```

- *NOTE: FC positive = less in KO; negative = more in KO*

**Optimal qvalue:**
--> *H3K27me3*; qval 20/25 (may worth to check qval 10 too)
--> *EZH2*; qval 10
--> *SUZ12*; qval 10 
--> *H3K4me3*; qval 20/25 (may worth to check qval 10 too)


# Compare CutRun and RNAseq expression

Let's do:
- Assign diff peaks to genes (Assign with different qvalues)
- Save/output a list of genes in geneSymbol format (Gain Lost)
- Generate a GTF from this list of genes
--> Do heatmap and profile plot
- Check expression of this list of genes
--> Do volcano plot (ugly) Then (enhancedVolcano)



### Assign THOR-diff peaks to genes and check expression

Now let's compare RNAseq (expression) and CutRun for THOR qval 15 among others:
- Filter HETvsWT and KOvsWT diff bound genes into **gain and loss H3K27me3** (10 20 25 for all)
- **Keep only signal in Promoter, gene body and TES** (ie. filter out peak assigned to intergenic)
- **Merge with deseq2** log2FC data (tpm will not work as too variable; or log2tpm maybe?)
- Plot in x FC and y baseMean=deseq2-norm counts (+ color qvalue) with facet_wrap~gain or lost (ie. volcano plot gain/lost)

*--> FC < 1  More in WT/Less in KO | FC > 1 More in KO*

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


# Import diff. peaks



## qval10
H3K27me3 = read.table('output/THOR/THOR_NPC_H3K27me3/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
EZH2 = read.table('output/THOR/THOR_NPC_EZH2/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
SUZ12 = read.table('output/THOR/THOR_NPC_SUZ12/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
H3K4me3 = read.table('output/THOR/THOR_NPC_H3K4me3/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)




## qval20
H3K27me3 = read.table('output/THOR/THOR_NPC_H3K27me3/THOR_qval20.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
EZH2 = read.table('output/THOR/THOR_NPC_EZH2/THOR_qval20.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
SUZ12 = read.table('output/THOR/THOR_NPC_SUZ12/THOR_qval20.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
H3K4me3 = read.table('output/THOR/THOR_NPC_H3K4me3/THOR_qval20.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)


## qval25
H3K27me3 = read.table('output/THOR/THOR_NPC_H3K27me3/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
EZH2 = read.table('output/THOR/THOR_NPC_EZH2/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
SUZ12 = read.table('output/THOR/THOR_NPC_SUZ12/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
H3K4me3 = read.table('output/THOR/THOR_NPC_H3K4me3/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)




## qval30
H3K27me3 = read.table('output/THOR/THOR_NPC_H3K27me3/THOR_qval30.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
EZH2 = read.table('output/THOR/THOR_NPC_EZH2/THOR_qval30.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
SUZ12 = read.table('output/THOR/THOR_NPC_SUZ12/THOR_qval30.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
H3K4me3 = read.table('output/THOR/THOR_NPC_H3K4me3/THOR_qval30.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)


## qval40
H3K27me3 = read.table('output/THOR/THOR_NPC_H3K27me3/THOR_qval40.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)
## qval50
H3K27me3 = read.table('output/THOR/THOR_NPC_H3K27me3/THOR_qval50.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V13, FC=V14, count_WT = V11,count_KO =V12) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT,count_KO)


# Tidy peaks #-->> Re-Run from here with different qvalue!!
H3K27me3_gr = makeGRangesFromDataFrame(H3K27me3,keep.extra.columns=TRUE)
EZH2_gr = makeGRangesFromDataFrame(EZH2,keep.extra.columns=TRUE)
SUZ12_gr = makeGRangesFromDataFrame(SUZ12,keep.extra.columns=TRUE)
H3K4me3_gr = makeGRangesFromDataFrame(H3K4me3,keep.extra.columns=TRUE)

gr_list <- list(H3K27me3=H3K27me3_gr, EZH2=EZH2_gr,  SUZ12=SUZ12_gr, H3K4me3=H3K4me3_gr)

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
H3K27me3_annot <- as.data.frame(peakAnnoList[["H3K27me3"]]@anno)
EZH2_annot <- as.data.frame(peakAnnoList[["EZH2"]]@anno)
SUZ12_annot <- as.data.frame(peakAnnoList[["SUZ12"]]@anno)
H3K4me3_annot <- as.data.frame(peakAnnoList[["H3K4me3"]]@anno)

## Convert entrez gene IDs to gene symbols
H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
EZH2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
EZH2_annot$gene <- mapIds(org.Hs.eg.db, keys = EZH2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
SUZ12_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = SUZ12_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
SUZ12_annot$gene <- mapIds(org.Hs.eg.db, keys = SUZ12_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
H3K4me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = H3K4me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
H3K4me3_annot$gene <- mapIds(org.Hs.eg.db, keys = H3K4me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")

## Save output table
write.table(H3K27me3_annot, file="output/ChIPseeker/annotation_H3K27me3_qval50.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

write.table(EZH2_annot, file="output/ChIPseeker/annotation_EZH2_qval30.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(SUZ12_annot, file="output/ChIPseeker/annotation_SUZ12_qval30.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(H3K4me3_annot, file="output/ChIPseeker/annotation_H3K4me3_qval30.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE

# load annotation tables
# H3K27me3_annot <- read.table("output/ChIPseeker/annotation_H3K27me3_qval10.txt", sep="\t", header=TRUE)


# Filter Gain/Loss sites
## KEEP Distal Intergenic (keep ALL)   ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K27me3_annot_gain = tibble(H3K27me3_annot) %>%
    filter(FC > 1) %>%
    add_column(binding = "gain")
H3K27me3_annot_lost = tibble(H3K27me3_annot) %>%
    filter(FC < (1/1)) %>%
    add_column(binding = "lost")
H3K27me3_annot_gain_lost = H3K27me3_annot_gain %>% 
    bind_rows(H3K27me3_annot_lost) 
### Save output gene lists
H3K27me3_annot_gain_geneSymbol = H3K27me3_annot_gain %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_annot_lost_geneSymbol = H3K27me3_annot_lost %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(H3K27me3_annot_gain_geneSymbol, file = "output/ChIPseeker/H3K27me3_annot_gain_qval50_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_annot_lost_geneSymbol, file = "output/ChIPseeker/H3K27me3_annot_lost_qval50_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)

EZH2_annot_gain = tibble(EZH2_annot) %>%
    filter(FC > 1) %>%
    add_column(binding = "gain")
EZH2_annot_lost = tibble(EZH2_annot) %>%
    filter(FC < (1/1)) %>%
    add_column(binding = "lost")
EZH2_annot_gain_lost = EZH2_annot_gain %>% 
    bind_rows(EZH2_annot_lost) 
### Save output gene lists
EZH2_annot_gain_geneSymbol = EZH2_annot_gain %>%
    dplyr::select(geneSymbol) %>%
    unique()
EZH2_annot_lost_geneSymbol = EZH2_annot_lost %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(EZH2_annot_gain_geneSymbol, file = "output/ChIPseeker/EZH2_annot_gain_qval25_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(EZH2_annot_lost_geneSymbol, file = "output/ChIPseeker/EZH2_annot_lost_qval25_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


SUZ12_annot_gain = tibble(SUZ12_annot) %>%
    filter(FC > 1) %>%
    add_column(binding = "gain")
SUZ12_annot_lost = tibble(SUZ12_annot) %>%
    filter(FC < (1/1)) %>%
    add_column(binding = "lost")
SUZ12_annot_gain_lost = SUZ12_annot_gain %>% 
    bind_rows(SUZ12_annot_lost) 
### Save output gene lists
SUZ12_annot_gain_geneSymbol = SUZ12_annot_gain %>%
    dplyr::select(geneSymbol) %>%
    unique()
SUZ12_annot_lost_geneSymbol = SUZ12_annot_lost %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(SUZ12_annot_gain_geneSymbol, file = "output/ChIPseeker/SUZ12_annot_gain_qval25_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(SUZ12_annot_lost_geneSymbol, file = "output/ChIPseeker/SUZ12_annot_lost_qval25_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


H3K4me3_annot_gain = tibble(H3K4me3_annot) %>%
    filter(FC > 1) %>%
    add_column(binding = "gain")
H3K4me3_annot_lost = tibble(H3K4me3_annot) %>%
    filter(FC < (1/1)) %>%
    add_column(binding = "lost")
H3K4me3_annot_gain_lost = H3K4me3_annot_gain %>% 
    bind_rows(H3K4me3_annot_lost) 
### Save output gene lists
H3K4me3_annot_gain_geneSymbol = H3K4me3_annot_gain %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_annot_lost_geneSymbol = H3K4me3_annot_lost %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(H3K4me3_annot_gain_geneSymbol, file = "output/ChIPseeker/H3K4me3_annot_gain_qval25_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_annot_lost_geneSymbol, file = "output/ChIPseeker/H3K4me3_annot_lost_qval25_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)



## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
H3K27me3_annot_gain = tibble(H3K27me3_annot) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(binding = "gain")
H3K27me3_annot_lost = tibble(H3K27me3_annot) %>%
    filter(FC < (1/1), annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(binding = "lost")
H3K27me3_annot_gain_lost = H3K27me3_annot_gain %>% 
    bind_rows(H3K27me3_annot_lost) 
### Save output gene lists
H3K27me3_annot_gain_geneSymbol = H3K27me3_annot_gain %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K27me3_annot_lost_geneSymbol = H3K27me3_annot_lost %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(H3K27me3_annot_gain_geneSymbol, file = "output/ChIPseeker/H3K27me3_annot_gain_qval50_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K27me3_annot_lost_geneSymbol, file = "output/ChIPseeker/H3K27me3_annot_lost_qval50_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


EZH2_annot_gain = tibble(EZH2_annot) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(binding = "gain")
EZH2_annot_lost = tibble(EZH2_annot) %>%
    filter(FC < (1/1), annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(binding = "lost")
EZH2_annot_gain_lost = EZH2_annot_gain %>% 
    bind_rows(EZH2_annot_lost) 
### Save output gene lists
EZH2_annot_gain_geneSymbol = EZH2_annot_gain %>%
    dplyr::select(geneSymbol) %>%
    unique()
EZH2_annot_lost_geneSymbol = EZH2_annot_lost %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(EZH2_annot_gain_geneSymbol, file = "output/ChIPseeker/EZH2_annot_gain_qval25_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(EZH2_annot_lost_geneSymbol, file = "output/ChIPseeker/EZH2_annot_lost_qval25_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


SUZ12_annot_gain = tibble(SUZ12_annot) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(binding = "gain")
SUZ12_annot_lost = tibble(SUZ12_annot) %>%
    filter(FC < (1/1), annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(binding = "lost")
SUZ12_annot_gain_lost = SUZ12_annot_gain %>% 
    bind_rows(SUZ12_annot_lost) 
### Save output gene lists
SUZ12_annot_gain_geneSymbol = SUZ12_annot_gain %>%
    dplyr::select(geneSymbol) %>%
    unique()
SUZ12_annot_lost_geneSymbol = SUZ12_annot_lost %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(SUZ12_annot_gain_geneSymbol, file = "output/ChIPseeker/SUZ12_annot_gain_qval25_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(SUZ12_annot_lost_geneSymbol, file = "output/ChIPseeker/SUZ12_annot_lost_qval25_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


H3K4me3_annot_gain = tibble(H3K4me3_annot) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(binding = "gain")
H3K4me3_annot_lost = tibble(H3K4me3_annot) %>%
    filter(FC < (1/1), annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(binding = "lost")
H3K4me3_annot_gain_lost = H3K4me3_annot_gain %>% 
    bind_rows(H3K4me3_annot_lost) 
### Save output gene lists
H3K4me3_annot_gain_geneSymbol = H3K4me3_annot_gain %>%
    dplyr::select(geneSymbol) %>%
    unique()
H3K4me3_annot_lost_geneSymbol = H3K4me3_annot_lost %>%
    dplyr::select(geneSymbol) %>%
    unique()

write.table(H3K4me3_annot_gain_geneSymbol, file = "output/ChIPseeker/H3K4me3_annot_gain_qval25_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(H3K4me3_annot_lost_geneSymbol, file = "output/ChIPseeker/H3K4me3_annot_lost_qval25_promoterAnd5_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)






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

**Summary optimal qval** in agreement with expression:
- H3K27me3: qval50 (close to values we ve got at 8wN; around 25% not in agreement with gene expr; lower qvalue it's worst)
- EZH2: qval10 (higher huge drop in nb of genes)
- SUZ12: qval10 (higher huge drop in nb of genes)
- H3K4me3: qval10 (higher huge drop in nb of genes)


Let's make **clean enhanced volcano plot at the optimal qvalue**:

--> Done in `001__RNAseq` labnote at `### NPC KO vs WT`

## deepTools on THOR diff peaks

Let's do heatmap plot as this one from *003__CutRun*: `scripts/matrix_TSS_10kb_THOR_THORq15HETpeaks_positive_negative.sh`: All diff peaks; separating gain and lost

- Separate positive and negative binding and create 2 bed files for each IP
- Display as heatmap



```bash
conda activate deeptools

# Separate gain / lost positive/ negative peaks
awk -F'\t' '$14 > 1' output/THOR/THOR_NPC_H3K27me3/THOR_qval50.bed > output/THOR/THOR_NPC_H3K27me3/THOR_qval50_positive.bed
awk -F'\t' '$14 < 1' output/THOR/THOR_NPC_H3K27me3/THOR_qval50.bed > output/THOR/THOR_NPC_H3K27me3/THOR_qval50_negative.bed
awk -F'\t' '$14 > 1' output/THOR/THOR_NPC_EZH2/THOR_qval10.bed > output/THOR/THOR_NPC_EZH2/THOR_qval10_positive.bed
awk -F'\t' '$14 < 1' output/THOR/THOR_NPC_EZH2/THOR_qval10.bed > output/THOR/THOR_NPC_EZH2/THOR_qval10_negative.bed
awk -F'\t' '$14 > 1' output/THOR/THOR_NPC_SUZ12/THOR_qval10.bed > output/THOR/THOR_NPC_SUZ12/THOR_qval10_positive.bed
awk -F'\t' '$14 < 1' output/THOR/THOR_NPC_SUZ12/THOR_qval10.bed > output/THOR/THOR_NPC_SUZ12/THOR_qval10_negative.bed
awk -F'\t' '$14 > 1' output/THOR/THOR_NPC_H3K4me3/THOR_qval10.bed > output/THOR/THOR_NPC_H3K4me3/THOR_qval10_positive.bed
awk -F'\t' '$14 < 1' output/THOR/THOR_NPC_H3K4me3/THOR_qval10.bed > output/THOR/THOR_NPC_H3K4me3/THOR_qval10_negative.bed

# heatmap with differential peaks in KO H3

sbatch scripts/matrix_TSS_10kb_THOR_H3K27me3_q50_peaks_positive_negative.sh # 5900737
sbatch scripts/matrix_TSS_10kb_THOR_EZH2_q10_peaks_positive_negative.sh # 5901275
sbatch scripts/matrix_TSS_10kb_THOR_SUZ12_q10_peaks_positive_negative.sh # 5901363
sbatch scripts/matrix_TSS_10kb_THOR_H3K4me3_q10_peaks_positive_negative.sh # 5901385


```





# Functional analysis

Let's do functional analysis GO and KEGG on the genes bound with EZH1 (and SUZ12 (and H3K27me3)); to check whether anything related to steroids

--> Gene list has been collected from Venn diagram (**420 genes EZH1+SUZ12 / 333 EZH1+SUZ12+H3K27me3**)

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


# import genes list
EZH1_SUZ12 = read.table(file = c("output/GO/Venn_overlap_PSC_EZH1_SUZ12.txt"), header = FALSE)
EZH1_SUZ12_H3K27me3 = read.table(file = c("output/GO/Venn_overlap_PSC_EZH1_SUZ12_H3K27me3.txt"), header = FALSE)



ego <- enrichGO(gene = as.character(EZH1_SUZ12_H3K27me3$V1),   # CHANGE NAME HERE !!!
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)
                
pdf("output/GO/dotplot_GO_BP_overlap_PSC_EZH1_SUZ12.pdf", width=7, height=7)
pdf("output/GO/dotplot_GO_BP_overlap_PSC_EZH1_SUZ12_H3K27me3.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()


## convert SYMBOL to entrezID
EZH1_SUZ12_ENTREZID = mapIds(org.Hs.eg.db, keys = EZH1_SUZ12$V1,
                       column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
EZH1_SUZ12_H3K27me3_ENTREZID = mapIds(org.Hs.eg.db, keys = EZH1_SUZ12_H3K27me3$V1,
                       column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

enrichKEGG <- enrichKEGG(gene   = as.character(EZH1_SUZ12_H3K27me3_ENTREZID),        # CHANGE NAME HERE !!!
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/GO/dotplot_KEGG_overlap_PSC_EZH1_SUZ12.pdf", width=7, height=5)
pdf("output/GO/dotplot_KEGG_overlap_PSC_EZH1_SUZ12_H3K27me3.pdf", width=7, height=6)
dotplot(enrichKEGG, showCategory=20)
dev.off()



```




# Follow up tasks

## Meeting 20231025

Investigate why and whether true that loss EZH2 SUZ12 with KO EZH1.

1- Check peak overlap between gain/lost EZH2 SUZ12 peaks. With bedtools intersect
2- Take all WT SUZ12 peaks and check signal in WT and KO

**1- Check peak overlap between gain/lost EZH2 SUZ12 peaks. With bedtools intersect**


```bash
conda activate BedToBigwig

# peak location
output/THOR/THOR_NPC_SUZ12/THOR_qval10_positive.bed # 1,123 (median = 400bp)
output/THOR/THOR_NPC_SUZ12/THOR_qval10_negative.bed # 5,508 (median = 400bp)
output/THOR/THOR_NPC_EZH2/THOR_qval10_positive.bed # 1,267 (median = 450bp)
output/THOR/THOR_NPC_EZH2/THOR_qval10_negative.bed # 3,840 (median = 450bp)

# bedtools intersect exact peak
bedtools intersect -v -a output/THOR/THOR_NPC_SUZ12/THOR_qval10_positive.bed -b output/THOR/THOR_NPC_EZH2/THOR_qval10_positive.bed | wc -l # 745
bedtools intersect -v -a output/THOR/THOR_NPC_SUZ12/THOR_qval10_negative.bed -b output/THOR/THOR_NPC_EZH2/THOR_qval10_negative.bed | wc -l # 5,291

```

- nb of POSITIVE SUZ12 peaks that do NOT overlap with EZH2 peaks = 745 (1,123 total peaks) = 34 % of overlap...
- nb of NEGATIVE SUZ12 peaks that do NOT overlap with EZH2 peaks = 5,291 (5,508 total peaks) = 4 %


--> Very poor overlap but the peaks are very small; like median size around 400bp...


So let's assign these peak to genes and then check the targeted genes. For that use the already processed files:

```bash
output/ChIPseeker/SUZ12_annot_gain_qval10_promoterAnd5_geneSymbol.txt 
output/ChIPseeker/SUZ12_annot_lost_qval10_promoterAnd5_geneSymbol.txt 
output/ChIPseeker/EZH2_annot_gain_qval10_promoterAnd5_geneSymbol.txt 
output/ChIPseeker/EZH2_annot_lost_qval10_promoterAnd5_geneSymbol.txt 
```



**2- Take all WT SUZ12 peaks and check signal in WT and KO**


```bash
conda activate deeptools


sbatch scripts/matrix_TSS_10kb_THOR_NPC_SUZ12_macs2_broadPeak.sh # 6470946

```


















