# Project

**H9 cell lines**

- 50dN native (FA was not so good so back to native!):
    - WT: H3K27me3, H3K27me1 (2 AB tested), H3K2ac, EZH2, , EZH1, IGG

--> All in simplicate


**Objectives:**
- Some issues with previous CutRun, here is more a test with few samples, only WT.

This time, *no nuclear purification has been performed.*

--> **Working samples can be added as aditional WT replicate!**


# CONCLUSION samples QC

--> 50dN: H3K27me3 WT (1 Bio rep) works


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

sbatch --dependency=afterany:17775901 scripts/bowtie2.sh # 17775970 ok
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
for file in slurm-17775970.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_17775970.txt
```

Add these values to `/home/roulet/001_EZH1_project/011__CutRun_50dN_nucl/samples_011.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >75% input reads as been uniquely mapped to the genome (90% non uniq) 



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.




```bash
conda activate bowtie2

sbatch --dependency=afterany:17775970 scripts/samtools_unique.sh # 17776169 ok
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

sbatch --dependency=afterany:17776169 scripts/bamtobigwig_unique.sh # 17776280 ok


```

- 50dN
*Pass*: 50dN_WT_H3K27me3, 50dN_WT_IGG
*Failed*: 50dN_WT_EZH1, 50dN_WT_EZH2, 50dN_WT_H3K27ac, 50dN_WT_H3K27me1AM, 50dN_WT_H3K27me1OR, 50dN_WT_SUZ12




## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all.sh # 17863696 xxx


# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels 50dN_WT_EZH1 50dN_WT_EZH2 50dN_WT_H3K27ac 50dN_WT_H3K27me1AM 50dN_WT_H3K27me1OR 50dN_WT_H3K27me3 50dN_WT_IGG 50dN_WT_SUZ12 \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf

## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels 50dN_WT_EZH1 50dN_WT_EZH2 50dN_WT_H3K27ac 50dN_WT_H3K27me1AM 50dN_WT_H3K27me1OR 50dN_WT_H3K27me3 50dN_WT_IGG 50dN_WT_SUZ12 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf


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






# Spike in factor


## Calculate histone content

--> This histone content will be used to generate a scaling factor which will be used to histone-scaled our library size. The calcul/method to follow is from `003__CutRun/output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt`

**Pipeline:**
- Count the histone barcode on the clean reads
- Calculate SF (group by sample (replicate) and AB and calculate the total nb of reads. Then proportion of reads = nb read in sample / total reads. SF = min(proportion) / sample proportion)


## Count the histone barcode on the clean reads



```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp.sh # 17961302 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_IGG.sh # 18046054 ok


```


--> It output the nb of reads found for each histone; then simply copy paste to the excell file `output/spikein/SpikeIn_QC_fastp_011.xlsx` in GoogleDrive

- `50dN_WT_H3K27me1AM`: enriched in H3K27me1
- `50dN_WT_H3K27me1OR`: enriched in H3K27me1
- `50dN_WT_H3K27me3`:  enriched in H3K27me3





## histone spike in factor

XXXXXXXXX below not modified XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


### Quality control plot

Then look at the xlsx file from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) to generate quality control plot. Use R cluster for vizualization (file is `spikein_QC.xlsx` in Google Drive), file in `output/spikein`.
```R
# package
library("tidyverse")
library("readxl")
# import df adn tidy to remove AB used in sample_ID
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_011.xlsx") %>%
  separate(sample_ID, into = c("type", "condition", "tag"), sep = "_") %>%
  mutate(sample_ID = paste(type, condition, sep = "_")) %>%
  select(-type, -condition, -tag)


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
pdf("output/spikein/QC_histone_spike_in_H3K27me3.pdf", width = 6, height = 4)
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


## Histone scaling for H3K27me1AM
spikein_all_scale = spikein_all %>%
  group_by(sample_ID) %>%
  # Find the target_norm value when Target is H3K27me1 and AB is H3K27me1
  mutate(scaling_factor = ifelse(Target == "H3K27me1" & AB == "H3K27me1AM", target_norm, NA)) %>%
  # Fill the scaling_factor column with the appropriate value within each group
  fill(scaling_factor, .direction = "downup") %>%
  # Scale the target_norm values
  mutate(scaled_target_norm = target_norm / scaling_factor * 100) %>%
  # Remove the scaling_factor column
  select(-scaling_factor) %>%
  # Ungroup the data
  ungroup()
# Plot
pdf("output/spikein/QC_histone_spike_in_H3K27me1AM.pdf", width = 6, height = 4)
spikein_all_scale %>%
    filter(
           AB %in% c("H3K27me1AM", "IGG")) %>%
        ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
        geom_col(position = "dodge") +
        facet_wrap(~sample_ID, nrow=1) +
        geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        facet_wrap(~sample_ID)
dev.off()



## Histone scaling for H3K27me1OR
spikein_all_scale = spikein_all %>%
  group_by(sample_ID) %>%
  # Find the target_norm value when Target is H3K27me1 and AB is H3K27me1
  mutate(scaling_factor = ifelse(Target == "H3K27me1" & AB == "H3K27me1OR", target_norm, NA)) %>%
  # Fill the scaling_factor column with the appropriate value within each group
  fill(scaling_factor, .direction = "downup") %>%
  # Scale the target_norm values
  mutate(scaled_target_norm = target_norm / scaling_factor * 100) %>%
  # Remove the scaling_factor column
  select(-scaling_factor) %>%
  # Ungroup the data
  ungroup()
# Plot
pdf("output/spikein/QC_histone_spike_in_H3K27me1OR.pdf", width = 6, height = 4)
spikein_all_scale %>%
    filter(
           AB %in% c("H3K27me1OR", "IGG")) %>%
        ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
        geom_col(position = "dodge") +
        facet_wrap(~sample_ID, nrow=1) +
        geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        facet_wrap(~sample_ID)
dev.off()
```


--> H3K27me3 is enriched

--> H3K27me1 do not work, not enriched! New AB was testeed, that is a bad one, as previous one show enrichment (see `007__CutRun`)



