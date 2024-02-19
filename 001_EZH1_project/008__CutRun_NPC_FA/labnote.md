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
sbatch scripts/fastp_2.sh # 12505000 ok
sbatch scripts/fastp_3.sh # 12505001 ok
```

# FastQC

**raw**
```bash
sbatch scripts/fastqc_raw.sh # 12696467 xxx
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


# Non unique bigiwig
sbatch scripts/bamtobigwig_1.sh # 13164849 xxx
sbatch scripts/bamtobigwig_2.sh # 13164850 xxx
sbatch scripts/bamtobigwig_3.sh # 13164851 xxx
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


--> Non unique (all raw reads!) vs unique bigwig (less signal or more noise?): XXXXXXXXXXXXXXXXX


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

Let's do the analysis for H3K27me3 only; compare WT vs KO vs KOEF1a. Test 2 spikein normalization method (histone and Ecoli)


## Calculate histone content

--> This histone content will be used to generate a scaling factor which will be used to histone-scaled our library size. The calcul/method to follow is from `003__CutRun/output/spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt`

**Pipeline:**
- Count the histone barcode on the clean reads
- Calculate SF (group by sample (replicate) and AB and calculate the total nb of reads. Then proportion of reads = nb read in sample / total reads. SF = min(proportion) / sample proportion)


## Count the histone barcode on the clean reads



```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp.sh # 12922764 ok

```

--> It output the nb of reads found for each histone; then simply copy paste to the excell file `output/spikein/SpikeIn_QC_fastp_008.xlsx` in GoogleDrive

- `50dNFA_KOEF1aEZH1_H3K27me3`: enriched in H3K27me3, but much less nb of reads as compare to the NPC sample!
- `NPC_KO_H3K27ac`: not in histone control..
- `NPC_KOEF1aEZH1_H3K27me3`: enriched in H3K27me3
- `NPC_KOEF1aEZH1_H3K4me3`: enriched in H3K4me3
- `NPC_WT_H3K27ac`: not in histone control..
- `NPC_WT_H3K4me3`: enriched in H3K4me3



## histone spike in factor


```R
# package
library("tidyverse")
library("readxl")
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp_008.xlsx") 

## H3K27me3
spikein_H3K27me3 = spikein %>%
    filter(Target == "H3K27me3") %>%
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
    filter(Target == "H3K4me3") %>%
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

--> The H3K27me3 may not be taken into account as H3K27me3 FAIL for 50dN sample

--> The H3K4me3 SF looks good!




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


