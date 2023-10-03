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
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_1.sh # 5695005
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_2.sh # 5695008
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_3.sh # 5695009
```

--> It output the nb of reads found for each histone; then simply copy paste to the excell file `output/spikein/SpikeIn_QC_fastp.xlsx` in GoogleDrive

XXX

### Calculate SF
Calculate the scaling factor as with the raw reads from `output/spikein/SpikeIn_QC_fastp.xlsx`

```R
# package
library(tidyverse)
library(readxl)
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp.xlsx") # 

# Filter-in H3K27me3 counts for each samples
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

sbatch scripts/bamtobigwig_1.sh # 5696658
sbatch scripts/bamtobigwig_2.sh # 5696658
sbatch scripts/bamtobigwig_3.sh # 5696660
```



