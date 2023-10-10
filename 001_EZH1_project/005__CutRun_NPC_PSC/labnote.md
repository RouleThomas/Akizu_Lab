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

--> Analysis of histon content showed that there was a sample inversion. **NPC_KO_H3K27me1 = NPC_KO_H3K4me3**

--> H3K27me1 seems it did not work; even though histone enrichment...

--> All scaling factors combined in `output/spikein/spikein_histone_scaling_factor_fastp.txt`

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
















