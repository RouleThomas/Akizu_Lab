# pipeline used from previous analyses
Seems Shuo used same method as for the ChIP: 
- bowtie2: `-q --local --no-mixed --no-unal --dovetail --phred33 -x`
- samtools: `-bS -q 20 -f 0x2`

[Here](https://hbctraining.github.io/Intro-to-ChIPseq-flipped/CUT&RUN/Intro_to_CUT&RUN.html) different parameters are used. But whatever works is ok!

Compared to ChIP-seq experiments, CUT&RUN data has low background and low sequence depth, making peak calling more susceptible to false positives. Thus, peak calling from CUT&RUN data requires high specificity for true positive peaks instead of high sensitivity. Recommended not to use MACS2 but [SEACR](https://github.com/FredHutch/SEACR)
# Import files from Google drive to the cluster

Using `cp`; check integrity (volume) of files!! (`ls -lh`)

# File renaiming

**Move all fastq files within the input folder** (now [file] are in input/folder1/[file]) using a loop:
```bash
for file in input/*/;
    do mv "$file"* input/;
done

rmdir input/*
```
- `input/*/` this go 1 folder downstream; whatever their names
- we move file to `input/` directory 

**Rename files**
Detail about file renaiming can be found in `CutRun_infos.xlsx` in my Google drive folder.


```bash
# example for 1 file:
outdir="input"

x="R1_Het5_Igg_1"
raw_f="R1_Het5_Igg_CKDL220021572-1A_HF3THDSX5_L1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"

if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# Run slurm for all:
sbatch rename_raw.sh # 11452878 ok
```


# Spike-in quality control

Script from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) with quality control check has been modified to work with zipped fastq files and adapted to my specific nomenclature; also now specified samples analyzed...
```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript.sh # 11454380, for 8wN_HET_IGG_R1 8wN_HET_H3K27me3_R1 8wN_HET_IGG_R2 8wN_HET_H3K27me3_R2 8wN_KO_IGG_R1 only; launch other job for the rest of the samples:
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_1.sh # 11459027
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_2.sh # 11459026
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_3.sh # 11459025
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_4.sh # 11459024
```

Then look at the xlsx file from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) to generate quality control plot. Use R cluster for vizualization (file is `spikein_QC.xlsx` in Google Drive), file in `output/spikein`.
```R
# package
library(tidyverse)
library(readxl)
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC.xlsx")
spikein <- read_excel("output/spikein/SpikeIn_QC_update.xlsx") # slight mistake with missing values for 2 samples

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
pdf("output/spikein/QC_histone_spike_in.pdf", width = 10, height = 8)
spikein_all_scale %>%
  	ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
  	geom_col(position = "dodge") +
  	facet_wrap(~sample_ID) +
  	geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
  	theme_bw() +
  	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()
```

--> Spike in control analyses show that H3K27me3 is well enriched, all other Target and when using IGG showed less than 20% enrichment.


# Fastp trimming and fastqc
## Fastqc on raw reads
```bash
sbatch scripts/fastqc_WT.sh # 11454608
sbatch scripts/fastqc_HET.sh # 11454610
sbatch scripts/fastqc_KO.sh # 11454609
sbatch scripts/fastqc_patient.sh # 11454604
```

## Fastp to clean reads
```bash
sbatch scripts/fastp_WT.sh # 11454715
sbatch scripts/fastp_HET.sh # 11454712
sbatch scripts/fastp_KO.sh # 11454714
sbatch scripts/fastp_patient.sh # 11454699
```

## Fastqc on clean reads
```bash
sbatch scripts/fastqc_fastp_WT.sh # 11456187
sbatch scripts/fastqc_fastp_HET.sh # 11456189
sbatch scripts/fastqc_fastp_KO.sh # 11456188
sbatch scripts/fastqc_fastp_patient.sh # 11456186
```

--> fastqc showed clean and trimmed data after fastp


# Mapped clean reads
Let's test 3 different mapping parameters (same as the 3 tested for the ChIP; default, permissive-paired, permissive-unpaired) on `8wN_HET_H3K27me3_R1_1.fq.gz`

```bash
conda activate bowtie2

sbatch scripts/bowtie2_8wN_HET_H3K27me3_R1_default.sh # 11473163
sbatch scripts/bowtie2_8wN_HET_H3K27me3_R1_permissive_paired.sh # 11473165
sbatch scripts/bowtie2_8wN_HET_H3K27me3_R1_permissive_unpaired.sh # 11473166
```

- **permissive-paired** `--phred33 -q --local --no-mixed --no-unal --dovetail`: 
    - nb of uniquely mapped reads: 3286178 (55.04%)
    - >1 times: 2557545 (42.84%)
    - overall 98.5% alignment rate
- **default** `--phred33 -q --no-unal`: (bowtie2 default) nb of uniquely mapped reads:
    - nb of uniquely mapped reads: 4694225 (78.62%)
    - >1 times 1056952 (17.7%)
    - overall 97.73% alignement rate
- **permissive-unpaired** `--phred33 -q --local --no-unal --dovetail `: (can have uniquely map paired reads)  
    - nb of uniquely mapped reads: 3286178 (55.04%)
    - >1 times 2557545 (42.84%)
    - overall 98.8% alignement rate

--> Default better! Better concordant read alignment (with more uniquely mapped reads).

--> As for ChIP; let's try the following: `--phred33 -q --no-unal --no-mixed --dovetail` **endtoend**.

```bash
sbatch bowtie2_2dN_HET_H3K27me3_R1_endtoend_cutrun.sh # 11496424
```
- **endtoend** `--phred33 -q --no-unal --no-mixed --dovetail`: 
    - nb of uniquely mapped reads: 4697326 (78.66%)
    - >1 times: 1059589 (17.75%)
    - overall: 96.83%


Mapping for all samples with **endtoend** parameter:

```bash
sbatch scripts/bowtie2_HET.sh # 11498654 ok
sbatch scripts/bowtie2_KO.sh # 11498655 ok
sbatch scripts/bowtie2_WT.sh # 11498658 ok
sbatch scripts/bowtie2_patient.sh # 11826852 ok
```


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-11498654.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_11498654.txt

for file in slurm-11498655.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_11498655.txt

for file in slurm-11498658.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_11498658.txt

for file in slurm-11826852.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_11826852.txt
```

Add these values to `/home/roulet/001_EZH1_project/003__CutRun/mapping_QC.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >75% input reads as been uniquely mapped to the genome


# Spike-in for normalization
## Method1_Map with bowtie2 on histone spike in:
- Generate FASTA file for H3K27me3 barcode sequences (manually)
- Build bowtie2 index on it `bowtie2-build snap_cutana_h3k27me3.fasta snap_cutana_h3k27me3_index`
- Map all samples to these barcode with more permissive parameter to allow align paired-end reads to small index `--phred33 -q --no-unal`
- Count the number of aligned reads to the spike-in control sequences for each sample `samtools view -S -F 4 -c sample_h3k27me3_rep1_spikein.sam > sample_h3k27me3_rep1_spikein_count.txt`
```ruby
Sample 1 (Genotype A, Replicate 1): 10,000 aligned reads
Sample 2 (Genotype A, Replicate 2): 12,000 aligned reads
Sample 3 (Genotype B, Replicate 1): 8,000 aligned reads
Sample 4 (Genotype B, Replicate 2): 14,000 aligned reads
```
- Then do for each sample: nb of spike-in aligned reads / (sum of all spike-in aligned reads for all samples) --> This is relative proportion of spike-in for each samples 
```ruby
Sample 1: 10,000 / (10,000 + 12,000 + 8,000 + 14,000) = 0.25
Sample 2: 12,000 / (10,000 + 12,000 + 8,000 + 14,000) = 0.30
Sample 3: 8,000 / (10,000 + 12,000 + 8,000 + 14,000) = 0.20 # here smallest scaling factor
Sample 4: 14,000 / (10,000 + 12,000 + 8,000 + 14,000) = 0.35
```
- Divide each scaling factor with the smalest scaling factor among all our sample --> This are the scaling factor for each samples
```ruby
Sample 1: 0.25 / 0.20 = 1.25
Sample 2: 0.30 / 0.20 = 1.50
Sample 3: 0.20 / 0.20 = 1.00
Sample 4: 0.35 / 0.20 = 1.75
```
--> ChIP-seq peak intensities, read counts, or other metrics of interest by dividing the original values by the corresponding scaling factor for each sample.
- Apply the scaling correction after calling peaks --> normalize the resulting peak intensities using the scaling factors calculated from the spike-in controls

```bash
# Generate FASTA files with H3K27me3 barcodes
nano meta/spike_in_H3K27me3.fa # cp paste sequences in fasta format

# Index the histone sequences
bowtie2-build meta/spike_in_H3K27me3.fa meta/spike_in_H3K27me3_index

# Map and count all samples to the spike_in_H3K27me3_index
sbatch scripts/bowtie2_spike_in_WT.sh # 11783734
sbatch scripts/bowtie2_spike_in_HET.sh # 11783736
sbatch scripts/bowtie2_spike_in_KO.sh # 11783735
```
--> This method did not work, possible because the fasta generated is too small, thus bowtie2 cannot align our samples. 

## Method2_Used barcode count from the quality control

- Count the occurrences of the barcodes for the H3K27me3 spike-in controls in from the raw, untrimmed input files. = aligned reads
```ruby
# H3K27me3 IP
Sample 1 (Genotype A, Replicate 1): 10,000 aligned reads
Sample 2 (Genotype A, Replicate 2): 12,000 aligned reads
Sample 3 (Genotype B, Replicate 1): 8,000 aligned reads
Sample 4 (Genotype B, Replicate 2): 14,000 aligned reads

# IgG IP
Sample 5 (Genotype A, Replicate 1): 6,000 aligned reads
Sample 6 (Genotype A, Replicate 2): 7,000 aligned reads
Sample 7 (Genotype B, Replicate 1): 5,000 aligned reads
Sample 8 (Genotype B, Replicate 2): 9,000 aligned reads
```
- Calculate the proportion of reads corresponding to the H3K27me3 spike-in controls for each sample **per IP**
```ruby
# H3K27me3 IP
Sample 1: 10,000 / (10,000 + 12,000 + 8,000 + 14,000) = 0.25
Sample 2: 12,000 / (10,000 + 12,000 + 8,000 + 14,000) = 0.30
Sample 3: 8,000 / (10,000 + 12,000 + 8,000 + 14,000) = 0.20 # here smallest scaling factor
Sample 4: 14,000 / (10,000 + 12,000 + 8,000 + 14,000) = 0.35

# IgG IP
Sample 5: 6,000 / (6,000 + 7,000 + 5,000 + 9,000) = 0.24
Sample 6: 7,000 / (6,000 + 7,000 + 5,000 + 9,000) = 0.28
Sample 7: 5,000 / (6,000 + 7,000 + 5,000 + 9,000) = 0.20 # here smallest scaling factor
Sample 8: 9,000 / (6,000 + 7,000 + 5,000 + 9,000) = 0.36
```
- Use these proportions to calculate the scaling factors for normalization.
```ruby
Sample 1: 0.25 / 0.20 = 1.25
Sample 2: 0.30 / 0.20 = 1.50
Sample 3: 0.20 / 0.20 = 1.00
Sample 4: 0.35 / 0.20 = 1.75
```

--> Let's calculate scaling factor in R with the `output/spikein/SpikeIn_QC.xlsx` file:
```R
# package
library(tidyverse)
library(readxl)
# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_update.xlsx") # slight mistake with missing values for 2 samples

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


write.table(spikein_H3K27me3_scaling_factor, file="output/spikein/spikein_histone_H3K27me3_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)
```
--> Histone spike-in scaling factor can be found here: `output/spikein/spikein_histone_H3K27me3_scaling_factor.txt`

--> Probably better to use my trim fastp clean reads to calculate scaling factor

## Method2_Used barcode count from fastp clean reads
### Count the barcode on the clean reads

```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_1.sh # 11785968 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_2.sh # 11786001 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_3.sh # 11786009 ok
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript_fastp_4.sh # 11786010
```

### Calculate scaling factor
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

## Method4_Map with bowtie2 on Ecoli genome
Here we used the alternative method:
- Generate FASTA file for E. coli genome (K12,MG1655), download from [here](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3)
- Build bowtie2 index on it `bowtie2-build genome.fasta genome_index`
- Map all samples to these barcode with same parameter as I used `--phred33 -q --no-unal --no-mixed --dovetail` (Try more permissive parameter to allow align paired-end reads to small index if that fail `--phred33 -q --no-unal`)
- Count the number of aligned reads to the spike-in control sequences for each sample `samtools view -S -F 4 -c sample_h3k27me3_rep1_spikein.sam > sample_h3k27me3_rep1_spikein_count.txt`
- Do the math for scaling factor, same method as when using histone spike-in

```bash
conda activate bowtie2

# Index the E coli genome
bowtie2-build meta/MG1655.fa meta/MG1655

# Map and count all samples to the E Coli genome
sbatch scripts/bowtie2_spike_Ecoli_WT.sh # 11789096 ok
sbatch scripts/bowtie2_spike_Ecoli_HET.sh # 11789101 ok
sbatch scripts/bowtie2_spike_Ecoli_KO.sh # 11789095 ok
sbatch scripts/bowtie2_spike_Ecoli_patient.sh # 11789602 ok
```
--> There is some uniq mapped reads, around less than 1%

Collect all read counts from the `output/spikein/*MG1655_count.txt` files generated into `output/spikein/SpikeIn_MG1655.xlsx`. Then calculate scaling factor in R:

```R
# package
library(tidyverse)
library(readxl)
library(ggpubr)

# import df
spikein <- read_excel("output/spikein/SpikeIn_MG1655.xlsx") 


# Total reads per IP
spikein_H3K27me3_total = spikein %>%
    group_by(AB) %>%
    mutate(total = sum(counts)) %>%
    ungroup() %>%
    distinct(AB, .keep_all = TRUE) %>%
    select(AB,total)

# Read proportion
spikein_H3K27me3_read_prop = spikein %>%
    left_join(spikein_H3K27me3_total) %>%
    mutate(read_prop = counts / total)

spikein_H3K27me3_read_prop_min = spikein_H3K27me3_read_prop %>%
    group_by(AB) %>%
    summarise(min_prop=min(read_prop))

# Scaling factor
spikein_H3K27me3_scaling_factor = spikein_H3K27me3_read_prop %>%
    left_join(spikein_H3K27me3_read_prop_min) %>%
    mutate(scaling_factor = read_prop/min_prop)


write.table(spikein_H3K27me3_scaling_factor, file="output/spikein/spikein_MG1655_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)
```

Let's make a plot and compare the all 4 methods and their respective scaling_factor:
```R
# import all files and rename
spikein_H3K27me3_scaling_factor_fastp_tidy = spikein_H3K27me3_scaling_factor_fastp %>%
    select(sample_ID, AB, scaling_factor) %>%
    add_column(method = "histone_fastp")

spikein_H3K27me3_scaling_factor_raw_tidy = spikein_H3K27me3_scaling_factor_raw %>%
    select(sample_ID, AB, scaling_factor) %>%
    add_column(method = "histone_raw")

spikein_H3K27me3_scaling_factor_MG1655_tidy = spikein_H3K27me3_scaling_factor_MG1655 %>%
    select(sample_ID, AB, scaling_factor) %>%
    add_column(method = "Ecoli")

spikein_H3K27me3_scaling_factor_all = spikein_H3K27me3_scaling_factor_fastp_tidy %>%
    bind_rows(spikein_H3K27me3_scaling_factor_raw_tidy, spikein_H3K27me3_scaling_factor_MG1655_tidy)


# Plot
pdf("output/spikein/scaling_factor_all_methods.pdf", width = 14, height = 8)
pdf("output/spikein/scaling_factor_Ecolivshistone.pdf", width = 14, height = 8)
spikein_H3K27me3_scaling_factor_all %>%
    filter(method != "histone_raw") %>%
        ggplot(aes(x = sample_ID, y = scaling_factor, fill = method)) +
        geom_col(position = "dodge") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# Plot correlation
pdf("output/spikein/scaling_factor_Ecolivshistone_corr.pdf", width = 6, height = 6)
spikein_H3K27me3_scaling_factor_all %>%
    filter(method != "histone_raw") %>%
    select(sample_ID, method, scaling_factor) %>%
    spread(key = method, value = scaling_factor) %>%
        ggplot(.) +
            geom_point(aes(x = Ecoli, y = histone_fastp)) +
            geom_smooth(aes(x = Ecoli, y = histone_fastp), method = "lm", se = FALSE, linetype = "solid", color = "blue") +
            theme_bw() +
            xlab("Ecoli Scaling Factor") +
            ylab("Histone_fastp Scaling Factor") +
            stat_cor(aes(x = Ecoli, y = histone_fastp), label.x = Inf, label.y = Inf, hjust = 1, vjust = 1, size = 4)
dev.off()
```

--> Both histone and E coli normalization provide the overall similar scaling factor.

# Samtools and read filtering
## Reference/Endogeneous genome
```bash
sbatch scripts/samtools_HET.sh # 11578283 ok
sbatch scripts/samtools_KO.sh # 11578284, weirdly looong
sbatch scripts/samtools_WT.sh # 11578286 ok

sbatch scripts/samtools_patient.sh # 11850650 ok
```
--> `scripts/samtools_KO.sh` stop running for unknown reason. Or maybe just super-long, it is stuck at the `8wN_KO_IGG_R2` sample. Let's run again the whole script, with more memory and threads, and without sample 8wN_KO_IGG_R1 just to make sure the script is working. **Output in `output/tmp`**
```bash
sbatch scripts/samtools_KO_2.sh # ok
```
--> The long job has been cancel and samtools_KO_2.sh succesfully run.

--> New files transfered to the `output/bowtie2`. All is complete

## Ecoli/Exogeneous genome
```bash
sbatch scripts/samtools_MG1655.sh # 12046911; no I pick the histone spike in one; re-do with the correct MG1655 mapping
sbatch scripts/samtools_MG1655_corr.sh # 12098888 ok
```

# Generate wig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch scripts/bamtobigwig_WT.sh # 11827228 ok
sbatch scripts/bamtobigwig_HET.sh # 11827232 ok
sbatch scripts/bamtobigwig_KO.sh # 11827233 ok

sbatch scripts/bamtobigwig_patient.sh  # 11857843 ok
sbatch scripts/bamtobigwig_missing.sh # 12091320 ok
```

Generate median tracks:
```bash
conda activate BedToBigwig

sbatch scripts/bigwigmerge.sh # 39314 ok
```




Let's generate bigwig taking into account scaling factor
## Scaled bigwig
### E.coli-spike-in scaled-bigwig
```bash
conda activate deeptools

sbatch scripts/bamtobigwig_MG1655_scaled_WT.sh # 11829982 ok
sbatch scripts/bamtobigwig_MG1655_scaled_HET.sh # 11829983 ok
sbatch scripts/bamtobigwig_MG1655_scaled_KO.sh # 11829984 ok

sbatch scripts/bamtobigwig_MG1655_scaled_patient.sh # 11858276
```

### Histone-spike-in scaled-bigwig
```bash
conda activate deeptools

sbatch scripts/bamtobigwig_histone_scaled_WT.sh # 11830480 ok
sbatch scripts/bamtobigwig_histone_scaled_HET.sh # 11830482 ok
sbatch scripts/bamtobigwig_histone_scaled_KO.sh # 11830481 ok

sbatch scripts/bamtobigwig_histone_scaled_patient.sh # 11858339 ok
```
--> The bigwig spike-in normalized are very heterogeneous between replicates; smthing is wrong...

--> In the end, maybe they are not as heterogeneous when we correct with their respective IGG; they should be the good one as in such case KO show overall decrease of H3K27me3 which is expected!!; So let's make median and igg correction

```bash
conda activate BedToBigwig
sbatch scripts/bigwigmerge_histone_scaled.sh # 43449 ok
```

Now normalize with igg ratio and subtract:

**igg ratio + median**

```bash
conda activate deeptools
sbatch scripts/bigwig_histone_scaled_ratio.sh # 43476 ok

conda activate BedToBigwig
sbatch --dependency=afterany:43476 scripts/bigwigmerge_histone_scaled_ratio.sh # 43492 ok
```

Files looks good


**igg subtract + median**

```bash
conda activate deeptools
sbatch scripts/bigwig_histone_scaled_subtract.sh # 43501 ok

conda activate BedToBigwig
sbatch --dependency=afterany:43501 scripts/bigwigmerge_histone_scaled_subtract.sh # 43508 ok
```

Files looks good; subtract look seems more clean


### Histone-spike-in scaled-bigwig with RPGC

To mimick DiffBind let's do `--scaleFactor SF + RPGC normalization`. There is different norm method we can try (RPKM, CPM, BPM, RPGC, None); RPGC is 1x coverage (I saw many papers with it...)


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_histone_scaled_RPGC_WT.sh # 48330; re-run as overwrite 48706
sbatch scripts/bamtobigwig_histone_scaled_RPGC_HET.sh # 48332 ok
sbatch scripts/bamtobigwig_histone_scaled_RPGC_KO.sh # 48335 ok

sbatch scripts/bamtobigwig_histone_scaled_RPGC_patient.sh # 48336 ok
```
--> The replicates are very heterogeneous... Maybe will be corrected with the igg normalization but weird...

In `output/tmp`; let's try RPGC norm without applying `--scaleFactor`; to make sure scaleFactor is applying! `sbatch scripts/bamtobigwig_histone_scaled_RPGC_WT_test.sh # 48403`--> That's good, files are different (surprinsgly, much better in term of replicate homogeneity when applying RPGC solely lol). 


```bash
conda activate BedToBigwig
sbatch --dependency=afterany:48330:48332:48335:48336 scripts/bigwigmerge_histone_scaled_RPGC.sh # 48337 ok; re-run for WT as overwrite: ok
sbatch --dependency=afterany:48706 scripts/bigwigmerge_histone_scaled_RPGC_WT.sh # 48708 ok
```

--> Replicates are very heterogeneous; but the RPGC and normalization with scaleFactor are both taken into acount.


### Histone-spike-in and TMM-norm scaled-bigwig OR DiffBind-ScaleFactor

Here interesting discussion about [this](https://www.biostars.org/p/473442/)

To **account for both library size and spike-in** we can use the same **scaling factor** as we used with **DiffBind** (In DiffBind we corrected library size for spike-in and then apply TMM/LIB/RLE normalization for seq depth). 

I collected the different SF proposed by DiffBind, let's apply them on our files and see how they look (this should take into account the Igg!!!)

- LIB scaling factor:
 [1] 0.5996449 0.6214268 0.3329470 1.2442918 0.4196722 0.3844100 4.4879784
 [8] 0.3417159 0.7403840 1.0690635 0.7232678 1.0351977
- TMM scaling factor
 [1] 0.5753224 0.6668493 0.5153769 0.6837952 0.6657534 0.6721562 2.2223161
 [8] 0.5310329 0.7080880 0.8430827 0.7479317 0.7100108
- RLE scaling factor
 [1] 0.5839702 0.6672963 0.5075065 0.7000390 0.6664778 0.6647032 2.2382431
 [8] 0.5285152 0.7044995 0.8323578 0.7468121 0.7124129

--> HET; KO; WT; Rep 1 to 4

--> Reciprocal is to be used when converting bam to bigwig!


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_DiffBind_TMM_WT.sh # 48408 ok
sbatch scripts/bamtobigwig_DiffBind_TMM_HET.sh # 48409 ok
sbatch scripts/bamtobigwig_DiffBind_TMM_KO.sh # 48410 ok

# just for KO run non-reciprocal to see how it perform:
sbatch scripts/bamtobigwig_DiffBind_TMM_KO_NonReciprocal.sh # 48412 ok

# median
conda activate BedToBigwig
sbatch --dependency=afterany:48408:48409:48410 scripts/bigwigmerge_DiffBind_TMM.sh # 48710 ok
```
*NOTE: I used the same scaling factor for IP and IGG (not sure IGG should be used; maybe the SF already account for IGG)...*


--> The reciprocal are NOT good to use; so `output/bigwig_DiffBind_TMM` is good

--> replicates looks good!!!

--> The IGG looks good too; let's do the igg ratio/subtract



**igg ratio + median**

```bash
conda activate deeptools
sbatch scripts/bigwig_DiffBind_TMM_ratio.sh # 49452

conda activate BedToBigwig
sbatch --dependency=afterany:49452 scripts/bigwigmerge_DiffBind_TMM_ratio.sh # 49458
```

--> replicates looks good!!!

--> Not great for NEUROG2

**igg subtract + median**

```bash
conda activate deeptools
sbatch scripts/bigwig_DiffBind_TMM_subtract.sh # 49459

conda activate BedToBigwig
sbatch --dependency=afterany:49459 scripts/bigwigmerge_DiffBind_TMM_subtract.sh # 49460
```

--> replicates looks good!!!

--> Great for NEUROG2; still unsure whether IgG is already applied...




*NOTE: see note below; this below may be fail:*

What is wrong is that I should have use the reciprocal (1/n) and not n as scaling factor... Let' correct and save into `output/bigwig_histone_NotGenotypeGroup` (This is the true bigwig good to use, better than the *groupABgenotype* one)

```bash
conda activate deeptools
sbatch scripts/bamtobigwig_histone_scaled_WT_reciprocal.sh # 12370048 ok
sbatch scripts/bamtobigwig_histone_scaled_HET_reciprocal.sh # 12370046 ok
sbatch scripts/bamtobigwig_histone_scaled_KO_reciprocal.sh # 12370044 ok
sbatch scripts/bamtobigwig_histone_scaled_patient_reciprocal.sh # 12370043 ok
```

But that is not perfect yet as library size is not the same for all samples, for accurate comparison, we need the same final library size after spike in correction. 

**Not correct!! In the end, these versions (`output/bigwig_histone_NotGenotypeGroup`) are perfect!** Let's generate **median and Igg-norm** bigwig tracks:


```bash
conda activate BedToBigwig
sbatch scripts/bigwigmerge_histone_NotGenotypeGroup.sh # 19853 ok
```

**igg ratio + median**

```bash
conda activate deeptools
sbatch scripts/bigwig_histone_NotGenotypeGroup_ratio.sh # 19854 ok

conda activate BedToBigwig
sbatch --dependency=afterany:19854 scripts/bigwigmerge_histone_NotGenotypeGroup_ratio.sh # 19856 ok
```

Files looks good




**igg log2ratio + median**
```bash
conda activate deeptools
sbatch scripts/bigwig_histone_NotGenotypeGroup_log2ratio.sh # 19857 FAIL; 19865 ok

conda activate BedToBigwig
sbatch --dependency=afterany:19865 scripts/bigwigmerge_histone_NotGenotypeGroup_log2ratio.sh # 19858 FAIL; 19866 FAIL; 31427 ok
```

Files looks good



**igg substract + median**
```bash
conda activate deeptools
sbatch scripts/bigwig_histone_NotGenotypeGroup_subtract.sh # 19859 ok

conda activate BedToBigwig
sbatch --dependency=afterany:19859 scripts/bigwigmerge_histone_NotGenotypeGroup_subtract.sh # 19860 ok
```

Files looks great (overall better representation to assess the difference! Notably NEUROG2 region)



## Below is generation of the failed library depth scaled after histone scaling bigwig (igg_norm...)

```bash
conda activate BedToBigwig
sbatch scripts/bamtobigwig_histone_scaled_lib.sh # 12447710 ok
```

Save in `output/bigwig_histone_NotGenotypeGroup_lib`, the correct scaling factor to use are the `output/spikein/spikein_histone_groupABgenotype_scaling_factor.txt`.

Let's merge the bigwig into 1 file with wiggletools as median value (will do average of bigwig signal and not sum, many options see [github](https://github.com/Ensembl/WiggleTools)):

**Run wiggletools:**
```bash
conda activate BedToBigwig
sbatch scripts/bigwigmerge_histone_NotGenotypeGroup_lib.sh # 12450108 ok
```
*NOTE: bigwig are merge into 1 bedgraph which is then converted into 1 bigwig (wiggletools cannot output bigwig directly so need to pass by bedgraph or wiggle in between)*

--> No in the end the `output/bigwig_histone_NotGenotypeGroup` is better!

**IgG over H3K27me3 ratio + median**

Let's generate H3K27me3/IgG scaled histone-spike-in norm bigwig to be used with deepTools.

ratio and not log2ratio is better (= ratio of read counts per bin in the IP sample relative to the input sample). So run all samples; comand [here](https://deeptools.readthedocs.io/en/develop/content/tools/bigwigCompare.html):

```bash
conda activate deeptools
sbatch scripts/bigwig_histone_NotGenotypeGroup_lib_ratio.sh # 16250 ok

conda activate BedToBigwig
sbatch --dependency=afterany:16250 scripts/bigwigmerge_histone_NotGenotypeGroup_lib_ratio.sh # 16255 ok
```
Files in `output/bigwig_histone_NotGenotypeGroup_lib_IggNorm` and looks good


**IgG over H3K27me3 log2ratio + median**


```bash
conda activate deeptools
sbatch scripts/bigwig_histone_NotGenotypeGroup_lib_log2ratio.sh # 19359 ok

conda activate BedToBigwig
sbatch --dependency=afterany:19359 scripts/bigwigmerge_histone_NotGenotypeGroup_lib_log2ratio.sh # 19364 ok
```
Files in `output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio` and looks ok


**IgG over H3K27me3 subtract + median**


```bash
conda activate deeptools
sbatch scripts/bigwig_histone_NotGenotypeGroup_lib_subtract.sh # 19380 ok

conda activate BedToBigwig
sbatch --dependency=afterany:19380 scripts/bigwigmerge_histone_NotGenotypeGroup_lib_subtract.sh # 19382 ok
```
Files in `output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_subtract` and looks ok








# Peak calling_Version1 with failed normalized data

## SEACR peak calling
[Paper](https://doi.org/10.1186/s13072-019-0287-4) and [Github](https://github.com/FredHutch/SEACR)

Install required packages and prepare meta files
```bash
# Install bedtools within bowtie2 conda env
conda activate bowtie2
conda install -c bioconda bedtools
module load sam-bcf-tools/1.6
module load SAMtools/1.16.1* # New cluster

# Create a Tab delimited chr size file
samtools faidx ../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta # index genome
awk '{print $1 "\t" $2}' ../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.fai > ../../Master/meta/GRCh38_chrom_sizes.tab # extract chr name and size
```

**Preparing input bedgraph files** (Convert paired-end BAM to fragment bedgraph file):
```bash
# Example with 1 file (smaller in size)
bedtools bamtobed -bedpe -i output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam > output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bed > output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.clean.bed # Filter out >1000bp fragment
cut -f 1,2,6 output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.clean.bed | sort -k1,1 -k2,2n -k3,3n > output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.fragments.bed
bedtools genomecov -bg -i output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.fragments.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab > output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.fragments.bedgraph

# Run together
conda activate bowtie2

sbatch scripts/bamtobedgraph_WT_1.sh # 11826612 ok
sbatch scripts/bamtobedgraph_WT_2.sh # 11826613 ok
sbatch scripts/bamtobedgraph_HET_1.sh # 11826580 ok
sbatch scripts/bamtobedgraph_HET_2.sh # 11826581 ok
sbatch scripts/bamtobedgraph_KO_1.sh # 11826578 ok
sbatch scripts/bamtobedgraph_KO_2.sh # 11826579 ok

sbatch scripts/bamtobedgraph_patient.sh # 11858464 ok
```
*NOTE: At bedtools bamtobed warning that some reads have no mate. But that is because I filtered reads based on MAPQ and removed unmapped reads, secondary alignments, and reads failing quality check using samtools view command, so I might have removed one read from a pair while retaining the other*

**Run SCEA peak calling**:
```bash
# Example command
bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output

# Run together
## Stringeant
sbatch scripts/SEACR_WT.sh # 11826904; overwrite by mistake, relaunch: 11833019
sbatch scripts/SEACR_HET.sh # 11826911; overwrite by mistake, relaunch: 11833012
sbatch scripts/SEACR_KO.sh # 11826912; overwrite by mistake, relaunch: 11833018

sbatch scripts/SEACR_patient.sh # 11886323 ok

## Run all samples with relax (this is pretty fast)
sbatch scripts/SEACR_relax.sh # 11826920; overwrite by mistake, relaunch: 11833020

sbatch scripts/SEACR_relax_patient.sh # 11886498 ok
```
--> Very few peaks seems to have been called, so I did stringent and relax parameters

Maybe the very few peaks are due to the warnings about non-pair mate... Let's try to remove them and see if I have the same number of peaks. Let's do the test on the *8wN_WT_R1* samples

```bash
# Sort by query names
samtools sort -n -o output/bowtie2/8wN_WT_IGG_R1.dupmark.qname_sorted.bam \
    output/bowtie2/8wN_WT_IGG_R1.dupmark.bam

# Fix mate-pair information
samtools fixmate -m output/bowtie2/8wN_WT_IGG_R1.dupmark.qname_sorted.bam output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.bam

# Sort reads after fixing mate-pair information
samtools sort -o output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.bam \
    output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.bam

# Prepare input bedgraph
bedtools bamtobed -bedpe -i output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.bam > output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.bed > output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.clean.bed # Filter out >1000bp fragment
cut -f 1,2,6 output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.clean.bed | sort -k1,1 -k2,2n -k3,3n > output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.fragments.bed
bedtools genomecov -bg -i output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.fragments.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab > output/bowtie2/8wN_WT_IGG_R1.dupmark.fixmate.sorted.fragments.bedgraph
```
This does not change anything...

--> Need to investigate the wig files to see how my files looks

Issue may be because I used the non-spike-in normalized files. So let's transform the bedgraph files and multiply by the scaling factor, and run SEACR with the `norm` parameter.

```bash
conda activate bowtie2

# Scaled bedgraph
sbatch scripts/scaled_bedgraph_MG1655.sh # 11831789 ok
sbatch scripts/scaled_bedgraph_MG1655_patient.sh # 11887369 ok
sbatch scripts/scaled_bedgraph_histone.sh # 11832207 ok
sbatch scripts/scaled_bedgraph_histone_patient.sh # 11887755 ok


# Run together SEACR
## Stringeant
sbatch scripts/SEACR_MG1655_scaled.sh # 11833139 ok
sbatch scripts/SEACR_histone_scaled.sh # 11833144 ok

sbatch scripts/SEACR_MG1655_scaled_patient.sh # 11888295 ok
sbatch scripts/SEACR_histone_scaled_patient.sh # 11888924

## Run all samples with relax (this is pretty fast)
sbatch scripts/SEACR_MG1655_scaled_relax.sh # 11833148 ok
sbatch scripts/SEACR_histone_scaled_relax.sh # 11833149 ok

sbatch scripts/SEACR_MG1655_scaled_patient_relax.sh # 11889120 ok
sbatch scripts/SEACR_histone_scaled_patient_relax.sh # 11889642 ok
```

--> All looks good. But very few peaks... Cannot tweek the qvalue with SEACR also!



## SPIKER (MACS2) peak calling
Used blacklist from [ENCODE](https://github.com/Boyle-Lab/Blacklist). 

No need to downsamples samples here as I have my scaling factors.

Need to specify scaling factor using the macs2 command, Let's use spiker.sh for this. [Documentation](https://spiker.readthedocs.io/en/latest/usage.html) and [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8313745/).

**Install spiker**
```bash
conda create -n spiker 

conda activate spiker

cd ../../Master/software
pip3 install spiker
```
--> Failed at installing MACS2. So install MACS2 first, then, install spiker.
```bash
conda activate spiker
conda install -c bioconda macs2 # succeed
cd ../../Master/software
pip3 install spiker
```
--> Failed because zlib development files are missing, so install them within same conda env:
```bash 
conda install -c anaconda zlib
pip3 install spiker
```
--> Failed again because of failed macs2 and pysam dependencies, try install both together and retry:
```bash
conda install -c bioconda macs2 pysam # succeed
conda install -c bioconda pybigwig # succeed
pip3 install spiker
```
--> Failed again, looks like spiker try to install a more recent version of macs2, so let's try install it 1st
```bash
conda uninstall macs2
conda install -c bioconda macs2=2.2.7.1
pip3 install spiker
```
--> Cannot install it as my python version is too weak! Let's recreate an environment with python 3.11.2 which work with macs2 2.2.7.1:
```bash
conda env remove --name spiker
conda create -n spiker python=3.11.2

conda activate spiker
conda install -c bioconda macs2 # =2.2.7.1
```
--> No fail even worst cannot install macs2... Re-try looking at the [setup.py](https://github.com/liguowang/spiker/blob/main/setup.py) file from the github account to make sure we respect all dependencies
```bash
conda create --name spiker python=3.8
conda activate spiker
conda install -c bioconda numpy scipy pysam deeptools macs2 pyBigWig cython=0.29
pip3 install spiker
```
--> I also try installed with python3.11 and failed too. That tool is fuckin shit as hell to install; lets try downloading it from [here](https://pypi.org/project/spiker/#files) and mount it manually...
```bash
conda install -c anaconda zlib
conda install -c bioconda pysam 
python3 -m pip install spiker-1.0.5.tar.gz
```
--> Exact same fail. Let's try using a python envirnment
```bash
python3 -m venv spiker
source spiker/bin/activate
python3 -m pip install spiker-1.0.5.tar.gz
```
--> Fail again. Also tried within conda macs2 env + python spiker env... 

Let's try installing the package manually within a conda env, and remove from setup.py the requirement.
```bash
conda create --name spiker 
conda activate spiker
# Install the required dependencies using Conda (the one indicated in the setup.py from spiker Github)
conda install -c bioconda macs2
conda install -c bioconda scipy
conda install -c bioconda pysam
conda install -c bioconda deeptools  
conda install -c bioconda pyBigWig
# Clone the Spiker repository from GitHub
git clone https://github.com/liguowang/spiker.git
cd spiker
# Modify the setup.py to avoid it to run the install_require line by adding a "#"
## like this : install_requires=['macs2==2.2.7.1','pysam',
# install spiker
pip install .
# run it
python bin/spiker.py
```
Installation succesful, now try to run it. It failed, because it is not **using the python from my conda environment**. To correct this:
```bash
conda env list # To find path
conda activate test # test env where I create conda create --name test python=3.11
export PATH="/home/roulet/anaconda3/envs/test/bin:$PATH"
python --version # It now show 3.11 
```
So now here is how to run spiker:
```bash
conda activate spiker
export PATH="/home/roulet/anaconda3/envs/spiker/bin:$PATH"
spiker/bin/split_bam.py --help
```
It work, let's clean the installation so that I can use spiker to run it:
```bash
# Set spiker/bin/spiker.py executable (display in green)
chmod +x spiker/bin/spiker.py
```
It fail as spiker not installed, even though now pip install spiker work!! So here is the conclusion of this mess:

When creating a conda environment, make sure it use the python from the conda environment. Otherwise when you install python-dependent module, it will not work even though install... For that

**To reinstall it cleanly, do as follow:**

```bash
conda create --name spiker 
conda activate spiker
export PATH="/home/roulet/anaconda3/envs/spiker/bin:$PATH"
which python # To check the python from the environment is used
# Install the required dependencies using Conda (the one indicated in the setup.py from spiker Github)
conda install -c bioconda macs2
conda install -c bioconda scipy
conda install -c bioconda pysam
conda install -c bioconda deeptools  
conda install -c bioconda pyBigWig
# install spiker
pip install spiker
# run spiker
spiker.py --help
```

**Run spiker test**
Run [spiker](https://spiker.readthedocs.io/en/latest/usage.html#input-and-output) in broad as H3K27me3 peak; indicate scaling factor. 

Run as follow: `spiker.py --broad -t H3K27me3.sorted.bam -c IGG.sorted.bam --spikeIn --csf ControlScalingFactor --tsf treatmentScalingFactor -o H3K27ac`:
```bash
conda activate spiker
# Example/test for 1 file
spiker.py --broad -t output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam -c output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.bam --spikeIn --csf 2.0094484954298 --tsf 2.07929383602633 -o output/spiker/8wN_HET_H3K27me3_R1_MG1655scaled
```
--> It failed at step 5.1 because the write_bedGraph function in the MACS2 package is trying to write a string to a file that expects bytes-like objects. Let's use Python 3.6 instead of Python 3.7 which is more lenient when it comes to mixing bytes and strings.

**Re-install spiker**
```bash
conda env remove --name spiker 
conda create --name spiker python=3.6
conda activate spiker
which python # To check the python from the environment is used
export PATH="/home/roulet/anaconda3/envs/spiker/bin:$PATH"
which python # To check the python from the environment is used
conda install -c bioconda deeptoolsintervals pybigwig macs2
pip install pybit
conda install -c bioconda py2bit
pip install spiker
spiker.py --help
```
Still failed, so let's try a more recent version of macs2
```bash
conda uninstall macs2
pip3 install --upgrade pip
pip install MACS2
spiker.py --help
```
Fail again...
```bash
conda env remove --name spiker 
```
JUST GO TO HELL SPIKER!



# Peak calling on raw
## MACS2
**Create a macs2 environment**
```bash
conda create --name macs2
conda activate macs2
conda install -c bioconda macs2
```
The same parameter as for the ChIPseq are used; IGG as control instead of input.
```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_raw_WT.sh # 12075979 ok
sbatch scripts/macs2_raw_HET.sh # 12075980 ok
sbatch scripts/macs2_raw_KO.sh # 12075981 ok
sbatch scripts/macs2_raw_patient.sh # 12075982 ok
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
```

--> Overall the 0.005 qval is more accurate to call peak

## SEACR
Already tested previously and show very few peaks...


# Troubleshoot spike-in normalization
The previous normalization method used, seems shit. Bigwig are very heterogeneous between replicates, even more than the raw files. Let's try different approaches:
- Calculate **Histone spike-in factor from Cutana** (then generate scaled bigwig and scaled bam)
- **ChIPSeqSpike** (provide scaling factor and scaled bigwig; then need generate peaks)
- **[DiffBind]**(https://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) (provide scaling factor; then need generate bigwig and peaks)
- **SpikChIP** (provide scaled peak; need generate scaled bigwig)
- **Histone-EpiCypher guidelines**


## Histone spike-in factor from Cutana - Scaling factor
Calculate the scaling factor as with the raw reads from `output/spikein/SpikeIn_QC_fastp.xlsx`

Here we now group per genotype and AB:

```R
# package
library(tidyverse)
library(readxl)

# import df
spikein <- read_excel("output/spikein/SpikeIn_QC_fastp.xlsx") %>%
  mutate(genotype = str_extract(sample_ID, "HET|KO|WT|iPSCpatient")) %>% # Extract genotype from sample_ID
  filter(Target == "H3K27me3")


# Total reads per IP and genotype
spikein_H3K27me3_total = spikein %>%
  group_by(AB, genotype) %>%
  mutate(total = sum(counts)) %>%
  ungroup() %>%
  distinct(AB, genotype, .keep_all = TRUE) %>%
  select(AB, genotype, total)

# Read proportion
spikein_H3K27me3_read_prop = spikein %>%
  group_by(sample_ID) %>%
  mutate(counts = sum(counts)) %>%
  select(sample_ID,AB,genotype,counts) %>%
  unique() %>%
  left_join(spikein_H3K27me3_total) %>%
  mutate(read_prop = counts / total)

spikein_H3K27me3_read_prop_min = spikein_H3K27me3_read_prop %>%
  group_by(AB, genotype) %>%
  summarise(min_prop=min(read_prop))

# Scaling factor
spikein_H3K27me3_scaling_factor = spikein_H3K27me3_read_prop %>%
  left_join(spikein_H3K27me3_read_prop_min) %>%
  mutate(scaling_factor = read_prop/min_prop)

write.table(spikein_H3K27me3_scaling_factor, file="output/spikein/spikein_histone_groupABgenotype_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)

```

--> The scaling factor looks better (range between 1-3); let's generate bigwig and vizualize


### Histone spike-in factor from Cutana - Bigwig
```bash
sbatch scripts/bamtobigwig_histone_groupABgenotype_WT.sh # ok
```
--> Even though that is better that initially, the replicates are still very heterogeneous when loading on IGV.

--> Now let's divide instead of multiplying per the scaling factor...
```bash
conda activate deeptools
sbatch scripts/bamtobigwig_histone_groupABgenotype_WT_divide.sh # 12126291 ok
sbatch scripts/bamtobigwig_histone_groupABgenotype_allothers_divide_1.sh # 12343896
sbatch scripts/bamtobigwig_histone_groupABgenotype_allothers_divide_2.sh # 12343897
```

--> Also test normalization with the library size too as in ChIPseqSpikeInFree
```bash
conda activate deeptools
sbatch scripts/bamtobigwig_histone_groupABgenotype_WT_libscaled.sh # 12130826 ok
```

This is now MUCH BETTER!!! Both library-norm or not are almost the same. 



## Ecoli spike-in factor from Cutana - Scaling factor

Here we now group per genotype and AB:

```R
# package
library(tidyverse)
library(readxl)

# import df
spikein <- read_excel("output/spikein/SpikeIn_MG1655.xlsx") %>%
  mutate(genotype = str_extract(sample_ID, "HET|KO|WT|iPSCpatient"))


# Total reads per IP and genotype
spikein_H3K27me3_total = spikein %>%
  group_by(AB, genotype) %>%
  mutate(total = sum(counts)) %>%
  ungroup() %>%
  distinct(AB, genotype, .keep_all = TRUE) %>%
  select(AB, genotype, total)

# Read proportion
spikein_H3K27me3_read_prop = spikein %>%
  group_by(sample_ID) %>%
  mutate(counts = sum(counts)) %>%
  select(sample_ID,AB,genotype,counts) %>%
  unique() %>%
  left_join(spikein_H3K27me3_total) %>%
  mutate(read_prop = counts / total)

spikein_H3K27me3_read_prop_min = spikein_H3K27me3_read_prop %>%
  group_by(AB, genotype) %>%
  summarise(min_prop=min(read_prop))

# Scaling factor
spikein_H3K27me3_scaling_factor = spikein_H3K27me3_read_prop %>%
  left_join(spikein_H3K27me3_read_prop_min) %>%
  mutate(scaling_factor = read_prop/min_prop)

write.table(spikein_H3K27me3_scaling_factor, file="output/spikein/spikein_MG1655_groupABgenotype_scaling_factor.txt", sep="\t", quote=FALSE, row.names=FALSE)

```

--> The scaling factor looks ok



## Ecoli spike-in factor from Cutana - Bigwig

Do normalization taking into account library size:

```bash
conda activate deeptools
sbatch scripts/bamtobigwig_MG1655_groupABgenotype_libscaled.sh # 12202882 ok
```

--> looks not good; let's try the ChIPSeqSpike to obtain scaled bigwig directly




### ChIPSeqSpike
#### ChIPSeqSpike Installation
Need Bioconductor 3.10, doc [here](https://bioconductor.riken.jp/packages/3.10/bioc/html/ChIPSeqSpike.html)
```bash
srun --mem=50g --pty bash -l
conda create --name ChIPSeqSpike r-base=3.6.0 
conda install -c anaconda libxml2
conda install r-xml
conda install -c r r-lattice
conda activate ChIPSeqSpike
```
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPSeqSpike")
library("ChIPSeqSpike")
```
*NOTE: Many failed upon `BiocManager::install("ChIPSeqSpike")` in R; thus installed several R package through conda*

#### ChIPSeqSpike Run
Documentation can be found [here](https://rdrr.io/bioc/ChIPSeqSpike/f/inst/doc/ChIPSeqSpike.pdf)
```bash
srun --mem=50g --pty bash -l
conda activate ChIPSeqSpike
```
Generate sample metadata file (tab sep; created in .ods and copy paste to .txt) in `output/ChIPSeqSpike/` 

File renaiming
```bash
cd output/ChIPSeqSpike
for file in *.dupmark.sorted.bam; do
  new_filename="${file/.dupmark.sorted.bam/.bam}"
  mv "$file" "$new_filename"
done

for file in *.dupmark.sorted.bw; do
  new_filename="${file/.dupmark.sorted.bw/.bw}"
  mv "$file" "$new_filename"
done
```
*NOTE: ChIPSeqSpike required only one point per extension and that all files are in the same folder. To make it simple and avoid changing the name of my original files, dependent on other tool; I cp all bam/bigwig to `output/ChIPSeqSpike` and rename files*

Verification that all bam files contained aligned reads:
```bash
cd output/ChIPSeqSpike
current_folder=$(basename "$PWD")
output_file="aligned_read_counts.txt"

for file in *.bam; do
  echo -n "${current_folder}/${file}: " >> "$output_file"
  samtools view -c -F 4 "$file" >> "$output_file"
done
```
Let's make a quick plot here. I generated by hand a tidy version of the count aligned reads `output/ChIPSeqSpike/aligned_read_counts_tidy.xlsx`

```R
library("tidyverse")
library("readxl")

prop = read_excel("output/ChIPSeqSpike/aligned_read_counts_tidy.xlsx")
prop_percent = prop %>% mutate(prop=(count_MG1655/count_human) *100)

ggplot(prop_percent, aes(ID, prop, fill=genotype)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

```


```R
# Provide path to files
info_file_path <- "output/ChIPSeqSpike/meta_sample.txt"
bam_path <- "output/ChIPSeqSpike"
bigwig_path<- "output/ChIPSeqSpike"

# Create the ChIPSeqSpikeDatasetList
csds <- spikeDataset(info_file_path, bam_path, bigwig_path)
is(csds_test) # Should say "ChIPSeqSpikeDatasetList"

# Estimate scaling factor
csds_scale = estimateScalingFactors(csds, paired = TRUE)
## Save as R object
save(csds_scale, file = "output/ChIPSeqSpike/csds_scale.RData")


spikeSummary(csds_scale)
getRatio(csds_scale) # Give inf values... weird


# RPM scaling
csds_rpm <- scaling(csds)

# Input/IGG substraction
csds_input <- inputSubtraction(csds_rpm)

!!! IT FAIL HERE!!!

# Reverse rpm for exogenous scaling factor normalization
csds_reverse <- scaling(csds_input, reverse = TRUE)

# Application of the scaling factor
csds_exo <- scaling(csds_reverse, type = "exo")

```

--> The `csds_input <- inputSubtraction(csds_rpm)` has not been performed; need re-generate bigwig using same bin size between input/IGG. Let's re-generate the bigwig using another normalization (RPGC) and broader bin size (50 and not 1)

```bash
conda activate deeptools
sbatch scripts/bamtobigwig_ChIPSeqSpike.sh # 12287902 ok
```




### DiffBind
#### DiffBind Installation
Workshop [here](https://bioinformatics-core-shared-training.github.io/Quantitative-ChIPseq-Workshop/articles/Quantitative-ChIPseq-Workshop.html) and [manual](https://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf), [forum](https://support.bioconductor.org/p/9135565/) about spikein

```bash
conda create --name DiffBind -c bioconda bioconductor-diffbind
# FAIL
srun --mem=50g --pty bash -l
conda create --name DiffBind -c bioconda -c conda-forge bioconductor-diffbind libgcc-ng=9.3.0
conda install -c bioconda r-rlang
conda activate DiffBind
```
*troubleshoot: Upon `library("DiffBind")` it mention it needs r-rlang=>1.0.2 so I `conda install -c conda-forge r-rlang=1.0.2`; but it fail, as rlang part of tidyverse; I install tidyverse; but fail similarly with rlang. So update rlang and its dependencies in R `install.packages("rlang", dependencies = TRUE)`; try update it with conda `conda update rlang`; --> `conda install -c bioconda r-rlang` worked!*

```bash
srun --mem=100g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 
library("csaw") # For spikein norm

# Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample.txt", header = TRUE, sep = "\t"))

# Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count.RData")
load("output/DiffBind/sample_count.RData")

# plot
pdf("output/DiffBind/clustering_sample.pdf", width=14, height=20)  
plot(sample_count)
dev.off()

pdf("output/DiffBind/PCA_sample.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

# Blacklist is already applied, so let's generate GreyList (IGG)
## Greylist generation
sample_dba_greylist = dba.blacklist(sample_dba, blacklist=FALSE, greylist=TRUE)

# Now check how clustering look
sample_count_greylist = dba.count(sample_dba_greylist)
## This take time, here is checkpoint command to save/load:
save(sample_count_greylist, file = "output/DiffBind/sample_count_greylist.RData")
load("output/DiffBind/sample_count_greylist.RData")


# plot
pdf("output/DiffBind/clustering_greylist.pdf", width=14, height=20)
plot(sample_count_greylist)
dev.off()

pdf("output/DiffBind/PCA_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_greylist,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()


# Modeling and testing without spike in
## Define different contrast
sample_model <- dba.contrast(sample_count_greylist)
sample_model # Check the factor are correct (for us Treatment = genotype)
dba.show(sample_model, bContrast=T)

## Analyze the model
### Default setting
sample_model_analyze <- dba.analyze(sample_model, method=DBA_ALL_METHODS) # Here show DESEq2 and EDgeR method
dba.show(sample_model,bContrasts=TRUE) # HETvsKO = 88 (125); HETvsWT=25 (58); WTvsKO=132 (151): DESeq2 (edgeR)

#### Binding sites overlap with the 2 methods:
pdf("output/DiffBind/DBA_ALL_METHODS_Venn_contrast1.pdf", width=14, height=20) 
dba.plotVenn(sample_model,contrast=1,method=DBA_ALL_METHODS)
dev.off()
pdf("output/DiffBind/DBA_ALL_METHODS_Venn_contrast2.pdf", width=14, height=20) 
dba.plotVenn(sample_model,contrast=2,method=DBA_ALL_METHODS)
dev.off()
pdf("output/DiffBind/DBA_ALL_METHODS_Venn_contrast3.pdf", width=14, height=20) 
dba.plotVenn(sample_model,contrast=3,method=DBA_ALL_METHODS)
dev.off()


sample_model_analyze <- dba.analyze(sample_model, method=DBA_EDGER) # Let's pick the edgeR method

## Export the Diff Bind regions
### Convert to GR object
sample_model_report_contrast1 <- dba.report(sample_model_analyze,method=DBA_EDGER,contrast=1)
sample_model_report_contrast2 <- dba.report(sample_model_analyze,method=DBA_EDGER,contrast=2)
sample_model_report_contrast3 <- dba.report(sample_model_analyze,method=DBA_EDGER,contrast=3)
### Convert to bed and export
sample_model_report_contrast1_df <- data.frame(sample_model_report_contrast1)
sample_model_report_contrast2_df <- data.frame(sample_model_report_contrast2)
sample_model_report_contrast3_df <- data.frame(sample_model_report_contrast3)

colnames(sample_model_report_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_report_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_report_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_model_report_contrast1_df, file="output/DiffBind/sample_model_report_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_report_contrast2_df, file="output/DiffBind/sample_model_report_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_report_contrast3_df, file="output/DiffBind/sample_model_report_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_greylist_contrast1.pdf", width=14, height=20) 
pdf("output/DiffBind/plotMA_greylist_contrast2.pdf", width=14, height=20) 
pdf("output/DiffBind/plotMA_greylist_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_model_analyze,method=DBA_EDGER,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_greylist_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_model_analyze,method=DBA_EDGER, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_greylist_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_model_analyze, method=DBA_EDGER,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_greylist_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_model_analyze, method=DBA_EDGER,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_greylist_contrast1.pdf", width=14, height=20) 
plot(sample_model_analyze, method=DBA_EDGER, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_greylist_contrast2.pdf", width=14, height=20) 
plot(sample_model_analyze, method=DBA_EDGER, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_greylist_contrast3.pdf", width=14, height=20) 
plot(sample_model_analyze, method=DBA_EDGER, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_greylist_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_model_analyze, method=DBA_EDGER,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_greylist_contrast2.pdf", width=14, height=20) 
dba.plotPCA(sample_model_analyze, method=DBA_EDGER,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_greylist_contrast3.pdf", width=14, height=20) 
dba.plotPCA(sample_model_analyze, method=DBA_EDGER,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_greylist_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_model_analyze, method=DBA_EDGER,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_greylist_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_model_analyze, method=DBA_EDGER,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_greylist_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_model_analyze, method=DBA_EDGER,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()



# Now normalize with the spikein histone : THE CODE BELOW FAILED; let's try another approach...
## Collect the greylist counts and apply SF
SF <- c(
  `8wN_HET_H3K27me3_R1` = 1/1.50803188228081,
  `8wN_HET_H3K27me3_R2` = 1/1.79362354383814,
  `8wN_HET_H3K27me3_R3` = 1/1,
  `8wN_HET_H3K27me3_R4` = 1/2.71894543225015,
  `8wN_iPSCpatient_H3K27me3_R1` = 1/1.5370709981661,
  `8wN_KO_H3K27me3_R1` = 1/1.30657216494845,
  `8wN_KO_H3K27me3_R2` = 1/1,
  `8wN_KO_H3K27me3_R3` = 1/3.80953608247423,
  `8wN_KO_H3K27me3_R4` = 1/1.12899484536082,
  `8wN_WT_H3K27me3_R1` = 1/1,
  `8wN_WT_H3K27me3_R2` = 1/1.73020349058761,
  `8wN_WT_H3K27me3_R3` = 1/1.21332215532353,
  `8wN_WT_H3K27me3_R4` = 1/1.82966237329472
)



sample_count_greylist_spikein <- dba.normalize(sample_count_greylist, normalize = SF)




# plot
pdf("output/DiffBind/clustering_greylist_spikein.pdf", width=14, height=20)
plot(sample_count_greylist_spikein)
dev.off()

pdf("output/DiffBind/PCA_greylist_spikein.pdf", width=14, height=20) 
dba.plotPCA(sample_count_greylist_spikein,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()


# Modeling and testing without spike in
## Define different contrast
sample_model_spikein <- dba.contrast(sample_count_greylist_spikein)
sample_model_spikein # Check the factor are correct (for us Treatment = genotype)
dba.show(sample_model_spikein, bContrast=T)

## Analyze the model
### Default setting
sample_model_spikein_analyze <- dba.analyze(sample_model_spikein, method=DBA_ALL_METHODS) # Here show DESEq2 and EDgeR method
dba.show(sample_model_spikein_analyze,bContrasts=TRUE) # HETvsKO = 2 (125); HETvsWT=0 (58); WTvsKO=0 (151): DESeq2 (edgeR)


sample_model_spikein_analyze <- dba.analyze(sample_model_spikein, method=DBA_EDGER) # Let's pick the edgeR method

## Export the Diff Bind regions
### Convert to GR object
sample_model_spikein_report_contrast1 <- dba.report(sample_model_spikein_analyze,method=DBA_EDGER,contrast=1)
sample_model_spikein_report_contrast2 <- dba.report(sample_model_spikein_analyze,method=DBA_EDGER,contrast=2)
sample_model_spikein_report_contrast3 <- dba.report(sample_model_spikein_analyze,method=DBA_EDGER,contrast=3)
### Convert to bed and export
sample_model_spikein_report_contrast1_df <- data.frame(sample_model_spikein_report_contrast1)
sample_model_spikein_report_contrast2_df <- data.frame(sample_model_spikein_report_contrast2)
sample_model_spikein_report_contrast3_df <- data.frame(sample_model_spikein_report_contrast3)

colnames(sample_model_spikein_report_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_spikein_report_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_spikein_report_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_model_spikein_report_contrast1_df, file="output/DiffBind/sample_model_spikein_report_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_spikein_report_contrast2_df, file="output/DiffBind/sample_model_spikein_report_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_spikein_report_contrast3_df, file="output/DiffBind/sample_model_spikein_report_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_greylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_model_spikein_analyze,method=DBA_EDGER,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_greylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_model_spikein_analyze,method=DBA_EDGER,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_greylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_model_spikein_analyze,method=DBA_EDGER,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_greylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_model_spikein_analyze,method=DBA_EDGER, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_greylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_model_spikein_analyze, method=DBA_EDGER,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_greylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_model_spikein_analyze, method=DBA_EDGER,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_greylist_spikein_contrast1.pdf", width=14, height=20) 
plot(sample_model_spikein_analyze, method=DBA_EDGER, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_greylist_spikein_contrast2.pdf", width=14, height=20) 
plot(sample_model_spikein_analyze, method=DBA_EDGER, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_greylist_spikein_contrast3.pdf", width=14, height=20) 
plot(sample_model_spikein_analyze, method=DBA_EDGER, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_greylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_model_spikein_analyze, method=DBA_EDGER,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_greylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotPCA(sample_model_spikein_analyze, method=DBA_EDGER,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_greylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotPCA(sample_model_spikein_analyze, method=DBA_EDGER,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_greylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_model_spikein_analyze, method=DBA_EDGER,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_greylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_model_spikein_analyze, method=DBA_EDGER,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_greylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_model_spikein_analyze, method=DBA_EDGER,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
```
*NOTE: iPSC patient Rep2 has been removed*
*NOTE: default treshold for diffbind is FDR 0.05*
*NOTE: as iPSC has only 1 replicate I need to do it separately by hand. See doc here for [dba.contrast](https://www.rdocumentation.org/packages/DiffBind/versions/2.0.2/topics/dba.contrast)*
*NOTE: For the spike in, the samtools-clean bam did not work; seems there is no reads, so I use the raw bam file; but it does not work too*

**The step order to follow with DiffBind is:**
1. Count
2. Apply BlackList / GreyList
3. Count (with BlackList / GreyList applied)
4. Normalize (Spike in + TMM or RLE)
5. Identification of diff. bound regions


Let's Test **different normalization method (non spike in); using the bed 'TRUE peaks':**

Here is the generation of the Blacklist/Greylist count matrix:
```R
library("DiffBind") 
library("csaw") # For spikein norm

# Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample.txt", header = TRUE, sep = "\t"))

# Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)
## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count.RData")
load("output/DiffBind/sample_count.RData")

# plot
pdf("output/DiffBind/clustering_sample.pdf", width=14, height=20)  
plot(sample_count)
dev.off()

pdf("output/DiffBind/PCA_sample.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

# Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_dba, blacklist=FALSE, greylist=TRUE) # Here we apply blacklist and greylist

# Now check how clustering look
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
## This take time, here is checkpoint command to save/load:
save(sample_count_blackgreylist, file = "output/DiffBind/sample_count_blackgreylist.RData")
load("output/DiffBind/sample_count_blackgreylist.RData")


# plot
pdf("output/DiffBind/clustering_blackgreylist.pdf", width=14, height=20)
plot(sample_count_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

```

Now the spike-in normalization and test of TMM and lib depth (default) norm:

```R
library("DiffBind") 
library("csaw") # For spikein norm


# Load sample count Blacklist/Greylist applied
load("output/DiffBind/sample_count_blackgreylist.RData")

# Normalize 
sample_count_blackgreylist_spikein = dba.normalize(sample_count_blackgreylist, spikein=TRUE)

# plot
pdf("output/DiffBind/clustering_blackgreylist_spikein.pdf", width=14, height=20)
plot(sample_count_blackgreylist_spikein)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_spikein.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

# Set up contrast for comparison (all will be compare to the WT)
sample_count_blackgreylist_spikein_contrast = dba.contrast(sample_count_blackgreylist_spikein, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_spikein_contrast_analyze = dba.analyze(sample_count_blackgreylist_spikein_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)

# NO SPIKE IN FOR COMPARISON ##########
### Default lib-norm
sample_count_blackgreylist_contrast = dba.contrast(sample_count_blackgreylist, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_contrast_analyze = dba.analyze(sample_count_blackgreylist_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)

### TMM
sample_count_blackgreylist_TMM = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_TMM)

sample_count_blackgreylist_TMM_contrast = dba.contrast(sample_count_blackgreylist_TMM, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_TMM_contrast_analyze = dba.analyze(sample_count_blackgreylist_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)

################################


#### Binding sites overlap with the 2 methods:
pdf("output/DiffBind/DBA_ALL_METHODS_Venn_blackgreylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotVenn(sample_count_blackgreylist_spikein_contrast_analyze,contrast=1,method=DBA_ALL_METHODS)
dev.off()
pdf("output/DiffBind/DBA_ALL_METHODS_Venn_blackgreylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotVenn(sample_count_blackgreylist_spikein_contrast_analyze,contrast=2,method=DBA_ALL_METHODS)
dev.off()
pdf("output/DiffBind/DBA_ALL_METHODS_Venn_blackgreylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotVenn(sample_count_blackgreylist_spikein_contrast_analyze,contrast=3,method=DBA_ALL_METHODS)
dev.off()


sample_count_blackgreylist_spikein_contrast_DESEQ2 <- dba.analyze(sample_count_blackgreylist_spikein_contrast, method=DBA_DESEQ2) # Let's pick the edgeR method

## Export the Diff Bind regions
### Convert to GR object
sample_model_blackgreylist_spikein_report_contrast1 <- dba.report(sample_count_blackgreylist_spikein_contrast_DESEQ2,method=DBA_DESEQ2,contrast=1)
sample_model_blackgreylist_spikein_report_contrast2 <- dba.report(sample_count_blackgreylist_spikein_contrast_DESEQ2,method=DBA_DESEQ2,contrast=2)
sample_model_blackgreylist_spikein_report_contrast3 <- dba.report(sample_count_blackgreylist_spikein_contrast_DESEQ2,method=DBA_DESEQ2,contrast=3)
### Convert to bed and export
sample_model_blackgreylist_spikein_report_contrast1_df <- data.frame(sample_model_blackgreylist_spikein_report_contrast1)
sample_model_blackgreylist_spikein_report_contrast2_df <- data.frame(sample_model_blackgreylist_spikein_report_contrast2)
sample_model_blackgreylist_spikein_report_contrast3_df <- data.frame(sample_model_blackgreylist_spikein_report_contrast3)

colnames(sample_model_blackgreylist_spikein_report_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_blackgreylist_spikein_report_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_blackgreylist_spikein_report_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_model_blackgreylist_spikein_report_contrast1_df, file="output/DiffBind/sample_model_blackgreylist_spikein_report_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_blackgreylist_spikein_report_contrast2_df, file="output/DiffBind/sample_model_blackgreylist_spikein_report_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_blackgreylist_spikein_report_contrast3_df, file="output/DiffBind/sample_model_blackgreylist_spikein_report_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_blackgreylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_spikein_contrast_DESEQ2,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_spikein_contrast_DESEQ2,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_spikein_contrast_DESEQ2,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_blackgreylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_spikein_contrast_DESEQ2,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_blackgreylist_spikein_contrast1.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_spikein_contrast2.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_spikein_contrast3.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_blackgreylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_spikein_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_spikein_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_spikein_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_spikein_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3)
dev.off()



# Test different normalization method
## RLE normalization
sample_count_blackgreylist_spikein_RLE = dba.normalize(sample_count_blackgreylist, spikein=TRUE, normalize=DBA_NORM_RLE)

# plot
pdf("output/DiffBind/clustering_blackgreylist_spikein_RLE.pdf", width=14, height=20)
plot(sample_count_blackgreylist_spikein_RLE)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_spikein_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

# Set up contrast for comparison (all will be compare to the WT)
sample_count_blackgreylist_spikein_RLE_contrast = dba.contrast(sample_count_blackgreylist_spikein_RLE, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_spikein_RLE_contrast_analyze = dba.analyze(sample_count_blackgreylist_spikein_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)



## TMM normalization
sample_count_blackgreylist_spikein_TMM = dba.normalize(sample_count_blackgreylist, spikein=TRUE, normalize=DBA_NORM_TMM)

# plot
pdf("output/DiffBind/clustering_blackgreylist_spikein_TMM.pdf", width=14, height=20)
plot(sample_count_blackgreylist_spikein_TMM)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_spikein_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

# Set up contrast for comparison (all will be compare to the WT)
sample_count_blackgreylist_spikein_TMM_contrast = dba.contrast(sample_count_blackgreylist_spikein_TMM, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_spikein_TMM_contrast_analyze = dba.analyze(sample_count_blackgreylist_spikein_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)

# Generate plots


sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2 <- dba.analyze(sample_count_blackgreylist_spikein_TMM_contrast, method=DBA_DESEQ2) # Let's pick the edgeR method

## Export the Diff Bind regions
### Convert to GR object
sample_model_blackgreylist_spikein_report_contrast1 <- dba.report(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2,method=DBA_DESEQ2,contrast=1)
sample_model_blackgreylist_spikein_report_contrast2 <- dba.report(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2,method=DBA_DESEQ2,contrast=2)
sample_model_blackgreylist_spikein_report_contrast3 <- dba.report(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2,method=DBA_DESEQ2,contrast=3)
### Convert to bed and export
sample_model_blackgreylist_spikein_report_contrast1_df <- data.frame(sample_model_blackgreylist_spikein_report_contrast1)
sample_model_blackgreylist_spikein_report_contrast2_df <- data.frame(sample_model_blackgreylist_spikein_report_contrast2)
sample_model_blackgreylist_spikein_report_contrast3_df <- data.frame(sample_model_blackgreylist_spikein_report_contrast3)

colnames(sample_model_blackgreylist_spikein_report_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_blackgreylist_spikein_report_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_blackgreylist_spikein_report_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_model_blackgreylist_spikein_report_contrast1_df, file="output/DiffBind/sample_model_blackgreylist_spikein_TMM_report_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_blackgreylist_spikein_report_contrast2_df, file="output/DiffBind/sample_model_blackgreylist_spikein_TMM_report_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_blackgreylist_spikein_report_contrast3_df, file="output/DiffBind/sample_model_blackgreylist_spikein_TMM_report_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_blackgreylist_spikein_TMM_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_spikein_TMM_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_spikein_TMM_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_blackgreylist_spikein_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_spikein_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_spikein_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_blackgreylist_spikein_TMM_contrast1.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_spikein_TMM_contrast2.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_spikein_TMM_contrast3.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_spikein_TMM_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_spikein_TMM_contrast2.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_spikein_TMM_contrast3.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_blackgreylist_spikein_TMM_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_spikein_TMM_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_spikein_TMM_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_spikein_TMM_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()

```
--> Spike-in normalization help a lot for the identification of diff bound regions

--> Spike-in/Default(lib size)/DEseq2 and Spike-in/TMM/DEseq2 perform both great. Let's see how ChIPseq perform and use the same norm method for CutRun and ChIPseq

--> It appear that the E.coli/MG1655 spike in is very low (less than 1%), thus not appropriate for normalization


Let's troubleshoot Histone-EpiCypher guidelines for spike in normalization


### DiffBind - Standard method with histone spike-in normalization

Let's use the scaling factor I calculated previously from `spikein/spikein_histone_groupABgenotype_scaling_factor.txt`. Apply them in DiffBind and check on IGV that the scaling factor has been correctly applied (in the good way! Maybe test the value and its reciprocal)


Let's try generate a DBA object applying scaling factor as discussed [here](https://support.bioconductor.org/p/9147040/). Here I use **1/SF** as I want to downscale my samples.




```bash
srun --mem=100g --pty bash -l
conda activate DiffBind
```
```R
library("DiffBind") 
library("csaw") # For spikein norm

# Load the Blacklist/Greylist counts
load("output/DiffBind/sample_count.RData") # That is the raw (not greylist) one
load("output/DiffBind/sample_count_blackgreylist.RData")

# Apply scaling factor normalization (DiffBound n = edger 135/53/152 - deseq2 0/0/2)
sample_count_blackgreylist_histone_norm = dba.normalize(sample_count_blackgreylist, normalize = c(0.663115954,0.557530594,1,0.367789654,0.650588035,0.765361475,1,0.262499154,0.885743637,1,0.577966699,0.824183417,0.546548923))

# plot
pdf("output/DiffBind/clustering_blackgreylist_histone_norm.pdf", width=14, height=20)
plot(sample_count_blackgreylist_histone_norm)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_histone_norm.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()


sample_count_blackgreylist_histone_norm_contrast = dba.contrast(sample_count_blackgreylist_histone_norm, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_histone_norm_contrast_analyze = dba.analyze(sample_count_blackgreylist_histone_norm_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)




# Apply scaling factor normalization and lib size (DiffBound n = edger 135/53/152 - deseq2 88/26/124, backougrn TRUE deseq2= 88/26/124)
sample_count_blackgreylist_histone_norm = dba.normalize(sample_count_blackgreylist, normalize = c(0.663115954,0.557530594,1,0.367789654,0.650588035,0.765361475,1,0.262499154,0.885743637,1,0.577966699,0.824183417,0.546548923))

sample_count_blackgreylist_histone_norm_lib = dba.normalize(sample_count_blackgreylist_histone_norm, background = TRUE)


pdf("output/DiffBind/clustering_blackgreylist_histone_norm_lib.pdf", width=14, height=20)
plot(sample_count_blackgreylist_histone_norm_lib)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_histone_norm_lib.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_lib,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_histone_norm_lib_contrast = dba.contrast(sample_count_blackgreylist_histone_norm_lib, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_histone_norm_lib_contrast_analyze = dba.analyze(sample_count_blackgreylist_histone_norm_lib_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)




# Apply scaling factor normalization and RLE (DiffBound n = edger 135/53/152 - deseq2 92/11/62; background TRUE deseq2 = 75/26/105)
sample_count_blackgreylist_histone_norm = dba.normalize(sample_count_blackgreylist, normalize = c(0.663115954,0.557530594,1,0.367789654,0.650588035,0.765361475,1,0.262499154,0.885743637,1,0.577966699,0.824183417,0.546548923))

sample_count_blackgreylist_histone_norm_RLE = dba.normalize(sample_count_blackgreylist_histone_norm, normalize = DBA_NORM_RLE, background = TRUE)


pdf("output/DiffBind/clustering_blackgreylist_histone_norm_RLE.pdf", width=14, height=20)
plot(sample_count_blackgreylist_histone_norm_RLE)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_histone_norm_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_histone_norm_RLE_contrast = dba.contrast(sample_count_blackgreylist_histone_norm_RLE, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_histone_norm_RLE_contrast_analyze = dba.analyze(sample_count_blackgreylist_histone_norm_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)





# Apply scaling factor normalization and TMM (DiffBound n = edger 135/53/152 - deseq2 115/14/77; backougrn TRUE deseq2= 80/30/110)
sample_count_blackgreylist_histone_norm = dba.normalize(sample_count_blackgreylist, normalize = c(0.663115954,0.557530594,1,0.367789654,0.650588035,0.765361475,1,0.262499154,0.885743637,1,0.577966699,0.824183417,0.546548923))

sample_count_blackgreylist_histone_norm_TMM = dba.normalize(sample_count_blackgreylist_histone_norm, normalize = DBA_NORM_TMM, background =TRUE)


pdf("output/DiffBind/clustering_blackgreylist_histone_norm_TMM.pdf", width=14, height=20)
plot(sample_count_blackgreylist_histone_norm_TMM)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_histone_norm_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_histone_norm_TMM_contrast = dba.contrast(sample_count_blackgreylist_histone_norm_TMM, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_histone_norm_TMM_contrast_analyze = dba.analyze(sample_count_blackgreylist_histone_norm_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)



# I pick here the spike in normalization solely with EDGER

sample_count_blackgreylist_histone_norm = dba.normalize(sample_count_blackgreylist, normalize = c(0.663115954,0.557530594,1,0.367789654,0.650588035,0.765361475,1,0.262499154,0.885743637,1,0.577966699,0.824183417,0.546548923))


sample_count_blackgreylist_histone_norm_contrast = dba.contrast(sample_count_blackgreylist_histone_norm, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_histone_norm_contrast_EDGER <- dba.analyze(sample_count_blackgreylist_histone_norm_contrast, method=DBA_EDGER) # Let's pick the edgeR method



## Export the Diff Bind regions
### Convert to GR object
sample_model_blackgreylist_histone_norm_report_contrast1 <- dba.report(sample_count_blackgreylist_histone_norm_contrast_EDGER,method=DBA_EDGER,contrast=1)
sample_model_blackgreylist_histone_norm_report_contrast2 <- dba.report(sample_count_blackgreylist_histone_norm_contrast_EDGER,method=DBA_EDGER,contrast=2)
sample_model_blackgreylist_histone_norm_report_contrast3 <- dba.report(sample_count_blackgreylist_histone_norm_contrast_EDGER,method=DBA_EDGER,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_model_blackgreylist_histone_norm_report_contrast1_df <- data.frame(sample_model_blackgreylist_histone_norm_report_contrast1)
sample_model_blackgreylist_histone_norm_report_contrast2_df <- data.frame(sample_model_blackgreylist_histone_norm_report_contrast2)
sample_model_blackgreylist_histone_norm_report_contrast3_df <- data.frame(sample_model_blackgreylist_histone_norm_report_contrast3)

colnames(sample_model_blackgreylist_histone_norm_report_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_blackgreylist_histone_norm_report_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_model_blackgreylist_histone_norm_report_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_model_blackgreylist_histone_norm_report_contrast1_df, file="output/DiffBind/sample_model_blackgreylist_histone_norm_report_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_blackgreylist_histone_norm_report_contrast2_df, file="output/DiffBind/sample_model_blackgreylist_histone_norm_report_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_model_blackgreylist_histone_norm_report_contrast3_df, file="output/DiffBind/sample_model_blackgreylist_histone_norm_report_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)






# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_blackgreylist_histone_norm_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_histone_norm_contrast_EDGER,method=DBA_EDGER,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_histone_norm_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_histone_norm_contrast_EDGER,method=DBA_EDGER,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_histone_norm_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_histone_norm_contrast_EDGER,method=DBA_EDGER,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_blackgreylist_histone_norm_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_histone_norm_contrast_EDGER,method=DBA_EDGER, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_histone_norm_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_histone_norm_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_blackgreylist_histone_norm_contrast1.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_histone_norm_contrast2.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_histone_norm_contrast3.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_histone_norm_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_histone_norm_contrast2.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_histone_norm_contrast3.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_blackgreylist_histone_norm_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_histone_norm_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_histone_norm_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_histone_norm_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_histone_norm_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_histone_norm_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_histone_norm_contrast_EDGER, method=DBA_EDGER,contrast=3)
dev.off()



# Lets do plots for spike-in lib scaled

sample_count_blackgreylist_histone_norm = dba.normalize(sample_count_blackgreylist, normalize = c(0.663115954,0.557530594,1,0.367789654,0.650588035,0.765361475,1,0.262499154,0.885743637,1,0.577966699,0.824183417,0.546548923))

sample_count_blackgreylist_histone_norm_lib = dba.normalize(sample_count_blackgreylist_histone_norm)

sample_count_blackgreylist_histone_norm_lib_contrast = dba.contrast(sample_count_blackgreylist_histone_norm_lib, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))



sample_count_blackgreylist_histone_norm_contrast_DESEQ2 <- dba.analyze(sample_count_blackgreylist_histone_norm_lib_contrast, method=DBA_DESEQ2) # Let's pick 


# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_blackgreylist_histone_norm_lib_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_histone_norm_contrast_DESEQ2,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_histone_norm_lib_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_histone_norm_contrast_DESEQ2,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_histone_norm_lib_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_histone_norm_contrast_DESEQ2,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_blackgreylist_histone_norm_lib_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_histone_norm_contrast_DESEQ2,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_histone_norm_lib_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_histone_norm_lib_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_blackgreylist_histone_norm_lib_contrast1.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_histone_norm_lib_contrast2.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_histone_norm_lib_contrast3.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_histone_norm_lib_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_histone_norm_lib_contrast2.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_histone_norm_lib_contrast3.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_blackgreylist_histone_norm_lib_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_histone_norm_lib_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_histone_norm_lib_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_histone_norm_lib_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_histone_norm_lib_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_histone_norm_lib_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_histone_norm_contrast_DESEQ2, method=DBA_DESEQ2,contrast=3)
dev.off()
```
*NOTE: I played with the **library** parameter for normalization (`library = DBA_LIBSIZE_BACKGROUND, DBA_LIBSIZE_PEAKREADS`), but it does not change anythings. I also play with **background**: `background = TRUE`*

--> Whatever the method used, very few diff bound sites identified...

--> Also using spike-in norm but not library scaled, the result do not follow biological expectation (KO less H3K27me3 as compare to HET); here both are less than WT, but HET even less than KO. Same pattern occur when using the spike-in-lib scaled (or TMM or RLE) one.

Let's try to normalized my files without genotype grouping, in the end, it's weird to do this as it does not allow me to compare genotype!


```R
library("DiffBind") 
library("csaw") # For spikein norm

# Load the Blacklist/Greylist counts
load("output/DiffBind/sample_count.RData") # That is the raw (not greylist) one
load("output/DiffBind/sample_count_blackgreylist.RData")

# Apply scaling factor normalization (DiffBound n = edger 135/53/152 - deseq2 0/0/2)
sample_count_blackgreylist_histone_norm = dba.normalize(sample_count_blackgreylist, normalize = c(0.6309969100666774, 0.5305257400697344, 0.9515634580012263, 0.3499751950570516, 0.4408840406795075, 0.7653614754906817, 1, 0.2624991543197346, 0.8857436365711714, 0.5914183370169948, 0.3418201039555981, 0.4874371859296476, 0.3232390552755446))


# plot
pdf("output/DiffBind/clustering_blackgreylist_histone_norm_NOTgenotypegrouped.pdf", width=14, height=20)
plot(sample_count_blackgreylist_histone_norm)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_histone_norm_NOTgenotypegrouped.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()


sample_count_blackgreylist_histone_norm_contrast = dba.contrast(sample_count_blackgreylist_histone_norm, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_histone_norm_contrast_analyze = dba.analyze(sample_count_blackgreylist_histone_norm_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)


# Not lib size scaled lead to only 1 deg in deseq2... so let's normalize

sample_count_blackgreylist_histone_norm_lib = dba.normalize(sample_count_blackgreylist_histone_norm, normalize = DBA_NORM_RLE)


pdf("output/DiffBind/clustering_blackgreylist_histone_norm_RLE_NOTgenotypegrouped.pdf", width=14, height=20)
plot(sample_count_blackgreylist_histone_norm_lib)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_histone_norm_RLE_NOTgenotypegrouped.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_histone_norm_lib,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_histone_norm_lib_contrast = dba.contrast(sample_count_blackgreylist_histone_norm_lib, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_histone_norm_lib_contrast_analyze = dba.analyze(sample_count_blackgreylist_histone_norm_lib_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)


# Check whether not applying the histone normalization lead to same result:

sample_count_blackgreylist_RLE = dba.normalize(sample_count_blackgreylist, normalize = DBA_NORM_RLE, library = DBA_LIBSIZE_FULL)

normalize = c(0.6309969100666774, 0.5305257400697344, 0.9515634580012263, 0.3499751950570516, 0.4408840406795075, 0.7653614754906817, 1, 0.2624991543197346, 0.8857436365711714, 0.5914183370169948, 0.3418201039555981, 0.4874371859296476, 0.3232390552755446)



pdf("output/DiffBind/clustering_blackgreylist_RLE.pdf", width=14, height=20)
plot(sample_count_blackgreylist_RLE)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()



sample_count_blackgreylist_RLE_contrast = dba.contrast(sample_count_blackgreylist_RLE, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_RLE_contrast_analyze = dba.analyze(sample_count_blackgreylist_RLE, method=DBA_ALL_METHODS, bParallel = TRUE)

```
--> When normalizing with spike-in using `normalize = `, I cannot re-normalize the data with another normalize, it erase the normalization... 

--> Need to find a way to apply normalization SF AND then normalize per library-size

**Using our scaling factor, let's estimate the 'new' library size** and provide it to `dba.normalize(library = c(1000, 12000))` = Like that our library size will be change taking into account our scaling factor! **Then we can normalize with library-size, RLE or TMM**... (issue discussed [here](https://support.bioconductor.org/p/9147040/))


### Adjust library size with histone scaling factor and apply normalization
Total number of reads is our library size (used samtools flagstat to double check) :

`samtools flagstat output/bowtie2/*.dupmark.sorted.bam` used to obtain library size (first value=library size)

Calcul made in `/home/*R*/003_CutRun/Mapping_QC.xlsx`. Histone-norm-library-size = library-size * SF. Using the non-reciprocal scaling factor, we increase the library-size; the more histone enriched, the more library size is increased, thus the more signal will decrease.

Now let's use these new histone-scaled library size and normalize with library-size,TMM or RLE:

```R
library("DiffBind") 
library("csaw") # For spikein norm

# Load the Blacklist/Greylist counts
load("output/DiffBind/sample_count.RData") # That is the raw (not greylist) one
load("output/DiffBind/sample_count_blackgreylist.RData")

# LIB SIZE
sample_count_blackgreylist_LibHistoneScaled_LIB = dba.normalize(sample_count_blackgreylist, library = c(16435849, 17032874, 9125844, 34105171, 24753801, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047), normalize = DBA_NORM_LIB) # Default

pdf("output/DiffBind/clustering_blackgreylist_LibHistoneScaled_LIB.pdf", width=14, height=20)
plot(sample_count_blackgreylist_LibHistoneScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_LibHistoneScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_LIB,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_LibHistoneScaled_LIB_contrast = dba.contrast(sample_count_blackgreylist_LibHistoneScaled_LIB, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_LibHistoneScaled_LIB_contrast_analyze = dba.analyze(sample_count_blackgreylist_LibHistoneScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)


# RLE (DESEq2: 138, 19, 115)
sample_count_blackgreylist_LibHistoneScaled_RLE = dba.normalize(sample_count_blackgreylist, library = c(16435849, 17032874, 9125844, 34105171, 24753801, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047), normalize = DBA_NORM_RLE) # Default

pdf("output/DiffBind/clustering_blackgreylist_LibHistoneScaled_RLE.pdf", width=14, height=20)
plot(sample_count_blackgreylist_LibHistoneScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_LibHistoneScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_LibHistoneScaled_RLE_contrast = dba.contrast(sample_count_blackgreylist_LibHistoneScaled_RLE, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_LibHistoneScaled_RLE_contrast_analyze = dba.analyze(sample_count_blackgreylist_LibHistoneScaled_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)


# TMM (139/19/117, best one)

sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(16435849, 17032874, 9125844, 34105171, 24753801, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047), normalize = DBA_NORM_TMM) # Default

pdf("output/DiffBind/clustering_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_LibHistoneScaled_TMM_contrast = dba.contrast(sample_count_blackgreylist_LibHistoneScaled_TMM, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze = dba.analyze(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)




sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze = dba.analyze(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE)



# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()



## Export the Diff Bind regions
### Convert to GR object
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1 <- dba.report(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2 <- dba.report(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3 <- dba.report(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3)

colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



# test normalize with scaling library size only


sample_count_blackgreylist_LibHistoneScaledOnly = dba.normalize(sample_count_blackgreylist, library = c(16435849, 17032874, 9125844, 34105171, 24753801, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047)) # Default

pdf("output/DiffBind/clustering_blackgreylist_LibHistoneScaledOnly.pdf", width=14, height=20)
plot(sample_count_blackgreylist_LibHistoneScaledOnly)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_LibHistoneScaledOnly.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaledOnly,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()
```

--> Here the library are scaled taking into account histone spike in! It gave good clustering and OK results (few diff bound sites but I think these are correct)

--> TMM was overall the best normalization method (RLE perform well too): good clustering of samples and maximum nb of diff. sites

--> The diff. bound sites are very small, and also it seems there is very few detected but looking by eyes on the bigwig there should be more diff bound sites; notably, there should be more decrease in KO (it is clear we see it (on the `bigwig_histone` which in the end is the correct bigwig file to use (and not `bigwig_histone_NotGenotypeGroup`)))

So let's try to use the macs2 raw output to see whether we detect more differences; I removed patient as only 1 replicate cannot be used by DiffBind...

1st **re-generate count matrix:**
```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```


```R
library("DiffBind") 

# Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_macs2raw.txt", header = TRUE, sep = "\t"))

# Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)


## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_macs2raw.RData")
load("output/DiffBind/sample_count_macs2raw.RData")

# plot
pdf("output/DiffBind/clustering_sample_macs2raw.pdf", width=14, height=20)  
plot(sample_count)
dev.off()

pdf("output/DiffBind/PCA_sample_macs2raw.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

# Blacklist/Greylist generation
sample_dba_blackgreylist = dba.blacklist(sample_count, blacklist=TRUE, greylist=TRUE) # Here we apply blacklist and greylist

# Now check how clustering look
sample_count_blackgreylist = dba.count(sample_dba_blackgreylist)
## This take time, here is checkpoint command to save/load:
save(sample_count_blackgreylist, file = "output/DiffBind/sample_count_macs2raw_blackgreylist.RData")
load("output/DiffBind/sample_count_macs2raw_blackgreylist.RData")


# plot
pdf("output/DiffBind/clustering_macs2raw_blackgreylist.pdf", width=14, height=20)
plot(sample_count_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_macs2raw_blackgreylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()




# LIB SIZE (edgeR: 237/41/146; deseq2 10/1/16)
sample_count_blackgreylist_LibHistoneScaled_LIB = dba.normalize(sample_count_blackgreylist, library = c(16435849, 17032874, 9125844, 34105171, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047), normalize = DBA_NORM_LIB) # Default

## Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_LIB_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_LIB, bRetrieve=TRUE)


console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_LIB_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_LIB_SF.txt")
##


pdf("output/DiffBind/clustering_macs2raw_blackgreylist_LibHistoneScaled_LIB.pdf", width=14, height=20)
plot(sample_count_blackgreylist_LibHistoneScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_macs2raw_blackgreylist_LibHistoneScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_LIB,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_LibHistoneScaled_LIB_contrast = dba.contrast(sample_count_blackgreylist_LibHistoneScaled_LIB, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_LibHistoneScaled_LIB_contrast_analyze = dba.analyze(sample_count_blackgreylist_LibHistoneScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)


# RLE (edgeR: 237/41/146; deseq2 154/13/148)
sample_count_blackgreylist_LibHistoneScaled_RLE = dba.normalize(sample_count_blackgreylist, library = c(16435849, 17032874, 9125844, 34105171, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047), normalize = DBA_NORM_RLE)  # RLE

## Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_RLE_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_RLE, bRetrieve=TRUE)


console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_RLE_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_RLE_SF.txt")
##


pdf("output/DiffBind/clustering_macs2raw_blackgreylist_LibHistoneScaled_RLE.pdf", width=14, height=20)
plot(sample_count_blackgreylist_LibHistoneScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_macs2raw_blackgreylist_LibHistoneScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_LibHistoneScaled_RLE_contrast = dba.contrast(sample_count_blackgreylist_LibHistoneScaled_RLE, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_LibHistoneScaled_RLE_contrast_analyze = dba.analyze(sample_count_blackgreylist_LibHistoneScaled_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)


# TMM (edgeR: 237/41/146; deseq2 139/19/147)

sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(16435849, 17032874, 9125844, 34105171, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047), normalize = DBA_NORM_TMM) # TMM

## Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)


console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_blackgreylist_LibHistoneScaled_TMM_SF.txt")
##

pdf("output/DiffBind/clustering_macs2raw_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_macs2raw_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_blackgreylist_LibHistoneScaled_TMM_contrast = dba.contrast(sample_count_blackgreylist_LibHistoneScaled_TMM, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze = dba.analyze(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)


# TMM and RLE perform best; let's pick TMM
```
Still very few diff bound sites are detected; so let's tweak the qvalue! Discussion [here](https://www.biostars.org/p/444784/).


```R
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze)


## Export the Diff Bind regions
### Convert to GR object
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1 <- dba.report(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2 <- dba.report(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3 <- dba.report(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3)

colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)





# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
plot(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
```

--> with these parameters; neurog2 is significant (even FDR is significant... Weird). This looks great!!

--> Still issue is that the peak are small (small because it increase statistical power; longer the more background with mess the DEGs; the author even recommend keep 400bp lenght for broad mark as H3K27me3), it may be good to merge them in like a 1-5kb distance

Let's re-do this analysis but with the **MACS2 pre-filtered (qval 0.005; not the macs2raw) peaks and pvalue 0.05 DiffBind** (re-do a matrix; without the iPSC patient samples); meta is `meta_sample_q005.txt`:



```R
library("DiffBind") 


# Generate the sample metadata (in ods/copy paste to a .csv file)
sample_dba = dba(sampleSheet=read.table("output/DiffBind/meta_sample_q005.txt", header = TRUE, sep = "\t"))

# Batch effect investigation; heatmaps and PCA plots
sample_count = dba.count(sample_dba)



## This take time, here is checkpoint command to save/load:
save(sample_count, file = "output/DiffBind/sample_count_q005.RData")
load("output/DiffBind/sample_count_q005.RData")

# plot
pdf("output/DiffBind/clustering_sample_q005.pdf", width=14, height=20)  
plot(sample_count)
dev.off()

pdf("output/DiffBind/PCA_sample_q005.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()


# Greylist generation
sample_dba_greylist = dba.blacklist(sample_count, blacklist=FALSE, greylist=TRUE) # Here we apply greylist only

# Now check how clustering look
sample_count_greylist = dba.count(sample_dba_greylist)
## This take time, here is checkpoint command to save/load:
save(sample_count_greylist, file = "output/DiffBind/sample_count_q005_greylist.RData")
load("output/DiffBind/sample_count_q005_greylist.RData")


# plot
pdf("output/DiffBind/clustering_q005_greylist.pdf", width=14, height=20)
plot(sample_count_greylist)
dev.off()

pdf("output/DiffBind/PCA_q005_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_greylist,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()



# TMM norm
sample_count_greylist_LibHistoneScaled_TMM = dba.normalize(sample_count_greylist, library = c(16435849, 17032874, 9125844, 34105171, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047), normalize = DBA_NORM_TMM)

pdf("output/DiffBind/clustering_greylist_LibHistoneScaled_TMM.pdf", width=14, height=20)
plot(sample_count_greylist_LibHistoneScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_greylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_greylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_greylist_LibHistoneScaled_TMM_contrast = dba.contrast(sample_count_greylist_LibHistoneScaled_TMM, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze = dba.analyze(sample_count_greylist_LibHistoneScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE)



sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze)

# Examining results


## volcano plot with diff sites
pdf("output/DiffBind/volcano_greylist_LibHistoneScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_greylist_LibHistoneScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_greylist_LibHistoneScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()


## Export the Diff Bind regions
### Convert to GR object
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1 <- dba.report(sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2 <- dba.report(sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3 <- dba.report(sample_count_greylist_LibHistoneScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2)
sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df <- data.frame(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3)

colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_greylist_LibHistoneScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_greylist_LibHistoneScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_blackgreylist_LibHistoneScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_greylist_LibHistoneScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



# Collect SF:


sample_count_greylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_greylist_LibHistoneScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_greylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_greylist_LibHistoneScaled_TMM_SF.txt")
```
- **NOTE: The pre-filtered MACS2 analysis (qvalue 0.005) are labeled `*greylist*` or `*q005_greylist*`; the non-pre-filtered are label `*macs2raw*blacklistgreylist*` or `*blacklistgreylist*`**


--> The DiffBind scaling factor SF are Different than the macs2raw method!!!

--> Using qvalue 0.005 result in ~3-time less diff. bound sites identified ()



Now let's compare RNAseq (expression) and CutRun for macs2 qval 0.005:
- Filter HETvsWT and KOvsWT diff bound genes into **gain and loss H3K27me3**
- **Keep only signal in Promoter, gene body and TES** (ie. filter out peak assigned to intergenic)
- **Merge with deseq2** log2FC data (tpm will not work as too variable; or log2tpm maybe?)
- Plot in x FC and y baseMean=deseq2-norm counts (+ color qvalue) with facet_wrap~gain or lost (ie. volcano plot gain/lost)

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
library(VennDiagram)


# Import peaks
peaks_contrast1 = read.table('output/DiffBind/sample_count_greylist_LibHistoneScaled_TMM_contrast1_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast2 = read.table('output/DiffBind/sample_count_greylist_LibHistoneScaled_TMM_contrast2_df.bed')  %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast3 = read.table('output/DiffBind/sample_count_greylist_LibHistoneScaled_TMM_contrast3_df.bed')  %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11)  


# Tidy peaks
peaks_contrast1_gr = makeGRangesFromDataFrame(peaks_contrast1,keep.extra.columns=TRUE)
peaks_contrast2_gr = makeGRangesFromDataFrame(peaks_contrast2,keep.extra.columns=TRUE)
peaks_contrast3_gr = makeGRangesFromDataFrame(peaks_contrast3,keep.extra.columns=TRUE)

gr_list <- list(HETvsKO=peaks_contrast1_gr, HETvsWT=peaks_contrast2_gr, KOvsWT=peaks_contrast3_gr)


# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_DiffBind05_qval005.pdf", width=7, height=7)
vennplot(genes)
dev.off()

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
HETvsKO_annot <- as.data.frame(peakAnnoList[["HETvsKO"]]@anno)
HETvsWT_annot <- as.data.frame(peakAnnoList[["HETvsWT"]]@anno)
KOvsWT_annot <- as.data.frame(peakAnnoList[["KOvsWT"]]@anno)


## Convert entrez gene IDs to gene symbols
HETvsKO_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = HETvsKO_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
HETvsWT_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = HETvsWT_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KOvsWT_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KOvsWT_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")

HETvsKO_annot$gene <- mapIds(org.Hs.eg.db, keys = HETvsKO_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
HETvsWT_annot$gene <- mapIds(org.Hs.eg.db, keys = HETvsWT_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KOvsWT_annot$gene <- mapIds(org.Hs.eg.db, keys = KOvsWT_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(HETvsKO_annot, file="output/ChIPseeker/annotation_HETvsKO_qval005.txt", sep="\t", quote=F, row.names=F)
write.table(HETvsWT_annot, file="output/ChIPseeker/annotation_HETvsWT_qval005.txt", sep="\t", quote=F, row.names=F)
write.table(KOvsWT_annot, file="output/ChIPseeker/annotation_KOvsWT_qval005.txt", sep="\t", quote=F, row.names=F)


# load annotation tables
HETvsWT_annot <- read.table("output/ChIPseeker/annotation_HETvsWT_qval005.txt", sep="\t", header=TRUE)
KOvsWT_annot <- read.table("output/ChIPseeker/annotation_KOvsWT_qval005.txt", sep="\t", header=TRUE)

# Filter Gain/Loss sites
HETvsWT_annot_gain = tibble(HETvsWT_annot) %>%
    filter(FC > 0, annotation != "Distal Intergenic") %>%
    dplyr::select(-c(V6,V7,V8,strand,name)) %>%
    add_column(H3K27me3 = "gain")
HETvsWT_annot_lost = tibble(HETvsWT_annot) %>%
    filter(FC < 0, annotation != "Distal Intergenic") %>%
    dplyr::select(-c(V6,V7,V8,strand,name)) %>%
    add_column(H3K27me3 = "lost")
HETvsWT_annot_gain_lost = HETvsWT_annot_gain %>% 
    bind_rows(HETvsWT_annot_lost) %>%
    dplyr::select(seqnames,start,end,FC,pvalue,FDR,annotation,distanceToTSS,gene,geneSymbol,H3K27me3)

KOvsWT_annot_gain = tibble(KOvsWT_annot) %>%
    filter(FC > 0, annotation != "Distal Intergenic") %>%
    dplyr::select(-c(V6,V7,V8,strand,name)) %>%
    add_column(H3K27me3 = "gain")
KOvsWT_annot_lost = tibble(KOvsWT_annot) %>%
    filter(FC < 0, annotation != "Distal Intergenic") %>%
    dplyr::select(-c(V6,V7,V8,strand,name))  %>%
    add_column(H3K27me3 = "lost")
KOvsWT_annot_gain_lost = KOvsWT_annot_gain %>% 
    bind_rows(KOvsWT_annot_lost) %>%
    dplyr::select(seqnames,start,end,FC,pvalue,FDR,annotation,distanceToTSS,gene,geneSymbol,H3K27me3)

# Import RNAseq deseq2 output
HET_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_8wN_HET_vs_8wN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)
KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_8wN_KO_vs_8wN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

# Merge files
HETvsWT_annot_gain_lost_RNA = HETvsWT_annot_gain_lost %>% 
    left_join(HET_vs_WT) %>%
    dplyr::select(gene, H3K27me3,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()


KOvsWT_annot_gain_lost_RNA = KOvsWT_annot_gain_lost %>% 
    left_join(KO_vs_WT) %>%
    dplyr::select(gene, H3K27me3,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()


# Volcano plot
count_data <- HETvsWT_annot_gain_lost_RNA %>%
    group_by(H3K27me3, significance) %>%
    summarise(up = sum(log2FoldChange > 0),
              down = sum(log2FoldChange < 0),
              total = n()) %>%
    ungroup() %>%
    group_by(H3K27me3) %>%
    mutate(total_panel = sum(total)) %>%
    ungroup()

pdf("output/ChIPseeker/DiffBind05_qval005_HETvsWT_expression.pdf", width=7, height=4)
HETvsWT_annot_gain_lost_RNA %>%
    ggplot(aes(x = log2FoldChange, y = baseMean, color = significance)) +
        geom_point(alpha = 0.8, size = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        labs(title = "HET vs WT",
             subtitle = "Expression level of diff. bound H3K27me3 genes",
             x = "Log2 Fold Change",
             y = "Base Mean",
             color = "Significant (padj <= 0.05)") +
        facet_wrap(~H3K27me3) +
        theme_bw() +
        geom_text(data = count_data %>% filter(significance), 
                  aes(x = Inf, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
                  hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
        geom_text(data = count_data %>% distinct(H3K27me3, .keep_all = TRUE),
                  aes(x = Inf, y = -Inf, label = paste("Total:", total_panel, "genes")),
                  hjust = 1.1, vjust = -0.1, size = 3, color = "black")
dev.off()

count_data <- KOvsWT_annot_gain_lost_RNA %>%
    group_by(H3K27me3, significance) %>%
    summarise(up = sum(log2FoldChange > 0),
              down = sum(log2FoldChange < 0),
              total = n()) %>%
    ungroup() %>%
    group_by(H3K27me3) %>%
    mutate(total_panel = sum(total)) %>%
    ungroup()

pdf("output/ChIPseeker/DiffBind05_qval005_KOvsWT_expression.pdf", width=7, height=4)
KOvsWT_annot_gain_lost_RNA %>%
    ggplot(aes(x = log2FoldChange, y = baseMean, color = significance)) +
        geom_point(alpha = 0.8, size = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        labs(title = "KO vs WT",
             subtitle = "Expression level of diff. bound H3K27me3 genes",
             x = "Log2 Fold Change",
             y = "Base Mean",
             color = "Significant (padj <= 0.05)") +
        facet_wrap(~H3K27me3) +
        theme_bw() +
        geom_text(data = count_data %>% filter(significance), 
                  aes(x = Inf, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
                  hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
        geom_text(data = count_data %>% distinct(H3K27me3, .keep_all = TRUE),
                  aes(x = Inf, y = -Inf, label = paste("Total:", total_panel, "genes")),
                  hjust = 1.1, vjust = -0.1, size = 3, color = "black")
dev.off()

# Venn diagram to see whether gain H3K27me3 in HET are the same one in KO?
## GAIN
# Extract the unique gene lists as vectors
HET_gene_list <- unique(HETvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "gain") %>%
                          pull(gene))

KO_gene_list <- unique(KOvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "gain") %>%
                          pull(gene))


venn.plot <- venn.diagram(
  list(HET = HET_gene_list, KO = KO_gene_list),
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  category.names = c("HET Gain", "KO Gain"),
  cex = 2
)

# Plot the Venn diagram
pdf("output/ChIPseeker/Venn_DiffBind05_qval005_Gain_Het_KO.pdf", width=5, height=4)
grid.draw(venn.plot)
dev.off()

## LOST
# Extract the unique gene lists as vectors
HET_gene_list <- unique(HETvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "lost") %>%
                          pull(gene))

KO_gene_list <- unique(KOvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "lost") %>%
                          pull(gene))


venn.plot <- venn.diagram(
  list(HET = HET_gene_list, KO = KO_gene_list),
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  category.names = c("HET Lost", "KO Lost"),
  cex = 2
)

# Plot the Venn diagram
pdf("output/ChIPseeker/Venn_DiffBind05_qval005_Lost_Het_KO.pdf", width=5, height=4)
grid.draw(venn.plot)
dev.off()


## BOTH
# Extract the unique gene lists as vectors
HET_lost_gene_list <- unique(HETvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "lost") %>%
                          pull(gene))

KO_lost_gene_list <- unique(KOvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "lost") %>%
                          pull(gene))

HET_gain_gene_list <- unique(HETvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "gain") %>%
                          pull(gene))

KO_gain_gene_list <- unique(KOvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "gain") %>%
                          pull(gene))


venn.plot <- venn.diagram(
  list(HET_lost = HET_lost_gene_list, KO_lost = KO_lost_gene_list,
       HET_gain = HET_gain_gene_list, KO_gain = KO_gain_gene_list),
  filename = NULL,
  fill = c("blue", "red", "blue", "red"),
  alpha = 0.5,
  cex = 2
)

# Plot the Venn diagram
pdf("output/ChIPseeker/Venn_DiffBind05_qval005_LostGain_Het_KO.pdf", width=5, height=4)
grid.draw(venn.plot)
dev.off()
```

--> Using **macs2raw (qval0.05) is better**; more genes expression changes identifies; that goes in agreement with H3K27me3 changes







## Fine-tune diff bound peak calling (qvalue and peak merge)
Now that we have plenty of diff bound regions, we need to make sure what is detected is biologically signifcant. For that; we will verify that the pvalue 0.05 treshold is OK and increase it if necessary; then we will merge the peak in a 1-5kb distance.


--> The peak looks biologically relevant, even the one with small pvalue! Instead of merging the peaks, let's assign peak to genes with ChIPseeker and see whether diff bind sites follow gene expression changes




# Test other method for differential binding

Many different tools:
- use uniquely aligned reads only instead of MAPQ20 per default
- csaw modified [paper](https://f1000research.com/articles/4-1080) : looks shit as cannot deal with library size and all in R, looks painfull, and not famous...
- csaw original (look famous); [paper](https://academic.oup.com/nar/article/44/5/e45/2464481); not clear whether we can adjust lib size
- DiffChIPL [paper](https://academic.oup.com/bioinformatics/article/38/17/4062/6637512); method that claim to be the best (recent 2022); I think we can adjust lib size but not clear
- THOR [paper](https://academic.oup.com/nar/article/44/20/e153/2607977)
- ODIN [paper](https://academic.oup.com/bioinformatics/article/30/24/3467/2422257?login=false)


## use uniquely aligned reads only (instead of MAPQ20 per default)

- Re-generate alignment
- call peak with MACS2
- Proceed with DiffBind and compare the nb of diff bound sites identified (compare with )

YYY Not priority, test this ONLY if the ChIP show better results when using uniquely aligned reads...




## DiffChIPL

[paper](https://academic.oup.com/bioinformatics/article/38/17/4062/6637512) and [github](https://github.com/yancychy/DiffChIPL). Not a sliding window method, but claim working better than MACS2, and we can set library size!


### DiffChIPL-Installation
Create `DiffChIPL conda env`

```bash
conda create -n DiffChIPL -c conda-forge r-base=4.0.0
conda activate DiffChIPL
conda install -c conda-forge r-pkgdown # needed for devtools installation
conda install -c conda-forge r-devtools
```

```R
install.packages("usethis")
library("devtools")
install_github("yancychy/DiffChIPL")
```


Many fails... Dependencies fail...--> let's try copy `deseq2` conda env where many stuff installed: 

```bash
conda create --name DiffChIPL --clone deseq2
conda activate DiffChIPL
# was needed for install_github("yancychy/DiffChIPL") to work:
conda install -c bioconda bioconductor-limma 
conda install -c bioconda bioconductor-edger
conda install -c bioconda bioconductor-bamsignals
```

```R
install.packages("devtools")
library("devtools")

BiocManager::install("SGSeq")
install_github("OnofriAndreaPG/aomisc")
install_github("yancychy/DiffChIPL")

library("DiffChIPL")
```
*NOTE: these dependencies where problematic; so I installed them through conda or bioconductor or devtools_github: limma, edgeR, bamsignals, SGSeq, aomisc*

--> It works!!!

# DiffChIPL_Go 

Create the csv file as indicated in [github](https://github.com/yancychy/DiffChIPL/blob/main/example/simHist/sim_hist_1.csv) into `output/DiffChIPL/meta_sample_macs2raw.txt`; *DiffBind meta as been used!*


```bash
conda activate DiffChIPL
# convert tab sep file into csv
sed 's/\t/,/g' output/DiffChIPL/meta_sample_macs2raw.txt > output/DiffChIPL/meta_sample_macs2raw.csv
sed 's/\t/,/g' output/DiffChIPL/meta_sample_q005.txt > output/DiffChIPL/meta_sample_q005.csv
```
*NOTE: `spikein` column has been renmoved from the meta file*

```R
library("DiffChIPL")
library("statmod")

# meta file
flib="output/DiffChIPL/meta_sample_macs2raw.csv" # fail, cannot handle macs2 raw file
flib="output/DiffChIPL/meta_sample_q005.csv"



# WT vs HET 
flib="output/DiffChIPL/meta_sample_q005_WTvsHET.csv"
## Get read count from files
countL = getReadCount(inputF=flib)

save(countL, file = "output/DiffChIPL/countL_q005_WTvsHET.RData")

## Design matrix and normalization
# Define the variables
str1 = "8wN-H3K27me3-WT.vs.HET-Ridge"
group = c(1, 1, 1, 1, 0, 0, 0, 0)
ctrName = "WT"
treatName = "HET"
groupName = c(rep(treatName, 4), rep(ctrName, 4)) 

# Build the design matrix
design0 <- cbind(rep(1, 8), c(rep(1, 4), rep(0, 4)))
colnames(design0) <- c(ctrName, treatName)

# Display the design matrix
design0

# Update the count data (replace negative value with 0)
peakAll = cbind(as.character(countL$peakPos@seqnames), countL$peakPos@ranges@start,
                countL$peakPos@ranges@start+countL$peakPos@ranges@width-1)
rawid = paste0(countL$peakAll[,1],"_" ,peakAll[,2])
countAll = countL$countAll
rownames(countAll) = rawid


for(i in 1:ncol(countAll)){
  id = which(countAll[,i] < 1)
  countAll[id,i] = 0
}


# Normalize the count data
## cpmD = cpmNorm(countAll, libsize = fd$lsIP)
cpmD = cpmNorm(countAll, libsize = c(20293409, 29302285, 19824265, 28374047,16435849, 17032874, 9125844, 34105171)) # using spike-in corrected library-size; HET, KO, WT: c(16435849, 17032874, 9125844, 34105171, 24753801, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047)

cpmD = cpmNorm(countAll) # using spike-in corrected library-size; HET, KO, WT: c(16435849, 17032874, 9125844, 34105171, 24753801, 11502923, 10536410, 123012358, 9366195, 20293409, 29302285, 19824265, 28374047)


## differential analysis with DiffChIPL
pdf(file="output/DiffChIPL/resA.pdf")
resA = DiffChIPL(cpmD, design0, group0 = group)
dev.off()
fitRlimm3 = resA$fitDiffL
rtRlimm3 = resA$resDE

## Check the differential results
id_Rlimma_CPM = rownames(rtRlimm3[which(rtRlimm3$adj.P.Val < 0.05),])
 
rtRlimm3 = rtRlimm3[rownames(cpmD),]
aveE = rtRlimm3$AveExpr
logFC = rtRlimm3$logFC
padj = rtRlimm3$adj.P.Val

pdf(file="output/DiffChIPL/plotMAVoc2_WTvsHET.pdf")
plotMAVoc2(mean=aveE, logfc=logFC, adj.P.Val=padj, FC=1, padj=0.05, MA=TRUE,
          title=paste0("Rlimma-CPM \n", str1,"(padj<0.05)\n", 
                       length(id_Rlimma_CPM), " of ", nrow(rtRlimm3) ))
dev.off()

pdf(file="output/DiffChIPL/plotMAVoc2_WTvsHET_defaultLib.pdf")
plotMAVoc2(mean=aveE, logfc=logFC, adj.P.Val=padj, FC=1, padj=0.05, MA=TRUE,
          title=paste0("Rlimma-CPM \n", str1,"(padj<0.05)\n", 
                       length(id_Rlimma_CPM), " of ", nrow(rtRlimm3) ))
dev.off()


rtRlimm3_tibble = as_tibble(rtRlimm3)



# WT vs KO

```

- **NOTE: DiffChIPL by default remove dupplicates!!** So if perform better, maybe because it use file without duplicates? only unqiuely aligned reads?

--> Default parameter result in 20 differentially bound sites using scaled library size and 226 using default library-scaling 

--> There was multiple error and bug in this tool, shit as hell, do not recommend...



## THOR for diff. bound sites

THOR looks cool, used bam as input and we can provide scaling factor! in our case (spikein SF !!). Moreover we can use housekeeping genes for normalization! Looks great with many parameters that can be tweak. Here is the [paper](https://academic.oup.com/nar/article/44/20/e153/2607977) and some [tutorial](https://reg-gen.readthedocs.io/en/latest/thor/introduction.html) and [more](https://reg-gen.readthedocs.io/en/latest/thor/tool_usage.html).


### THOR installation

Install [RGT](https://reg-gen.readthedocs.io/en/latest/rgt/installation.html) in base env:

```bash
conda create --name RGT python # Always add python!!! Otherwise could use pyton from base!
conda activate RGT
pip install cython numpy scipy # will install within RGT conda env
pip install RGT --no-binary RGT
rgt-THOR --version

# error while running 1st time bigWigMerge cannot find libpng12.so.0

# Create a symbolic link to the missing libpng12.so.0; that direct toward ours libpng installed through conda
## First locate where the symbolic link should be created
ldd ~/anaconda3/envs/RGT/bin/bigWigMerge # this is to check if bigWigMerge found dependent libraries
## Second check which libpng should I use; is and where on my system
locate libpng16.so # to find where it is in my system
## Now create symbolic link that will use my libpng16; put the link where the software is looking for
ln -s /usr/lib64/libpng16.so ~/anaconda3/envs/RGT/lib/libpng12.so.0
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
ldd ~/anaconda3/envs/RGT/bin/bigWigMerge # now ldd found the libpng12.so.0 library !!!
bigWigMerge
```
--> Installation succesfful within `conda activate RGT`


### THOR diff. bound sites analysis

- THOR worked with one-by-one comparison, so let's start with *WT vs HET* then *WT vs KO*.
- We can provide SF for the bam... Which one to take...?? Maybe the histoneSpikeIn-DiffBind-TMM one as it is well normalized? Otherwise the histoneSpikeIn-NON-DiffBind-TMM normalized was weird (weird clustering)
**--> That way we use DiffBind to perform TMM and Spike in correction (and some data vizualization) and uses THOR for identification of diff. bound sites**

Here is histone-DiffBind-TMM SF for each bam (from which we scaled our bam to generate `output/bigwig_DiffBind_TMM`):
- output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam  1.412253844 
- output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam  1.186123259
- output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam  1.337020479 
- output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam  1.40842928 
- output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam  1.502057669
- output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam  1.487749425
- output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam  0.449980991 
- output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam  1.883122496 
- output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam  1.738155858
- output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam  1.499589188
- output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam  1.940327554
- output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam  1.462426177 


Here is histone-scaled reciprocal **NON-**DiffBind-TMM SF for each bam (from which we scaled our bam to generate `output/bigwig_histone_NotGenotypeGroup`); correct file to use, I double check for NEUROG2; and NOT `bigwig_histone`:
- output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam  0.5914183370169948
- output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam  0.3418201039555981
- output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam  0.4874371859296476
- output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam  0.3232390552755446
- output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam  0.7653614754906817
- output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam  1
- output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam  0.2624991543197346
- output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam  0.8857436365711714
- output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam  0.6309969100666774
- output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam  0.5305257400697344
- output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam  0.9515634580012263
- output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam  0.3499751950570516



**Create the CONFIG file** with nano `output/THOR/WTvsHET.config`:
```bash
#rep1
output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam
output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam
output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam
output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam
#rep2
output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam
output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam
output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam
output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam
#genome
../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
#chrom_sizes
../../Master/meta/GRCh38_chrom_sizes.tab
#inputs1
output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.bam
output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam
output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.bam
output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.bam
#inputs2
output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.bam
output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.bam
output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.bam
output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.bam
```
*NOTE: inputs are optional, it may worth trying without input too*

`output/THOR/WTvsKO.config`:
```bash
#rep1
output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam
output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam
output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam
output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam
#rep2
output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam
output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam
output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam
output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam
#genome
../../Master/meta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
#chrom_sizes
../../Master/meta/GRCh38_chrom_sizes.tab
#inputs1
output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.bam
output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam
output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.bam
output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.bam
#inputs2
output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.bam
output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.bam
output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.bam
output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.bam
```



**Run THOR**

*THOR is very buggy to make it work I need to temporaly change where to look for libraries lol.. So cannot use nano anymore for example...*

```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge

rgt-THOR --name WTvsHET --merge --output-dir output/THOR_WTvsHET --report --deadzones ../../Master/meta/hg38-blacklist.v2.bed --pvalue 0.1 --scaling-factors 1.412253844,1.186123259,1.337020479,1.40842928,1.738155858,1.499589188,1.940327554,1.462426177 output/THOR/WTvsHET.config

# seems to work so run within a slurm job with 200g mem
conda activate RGT
sbatch scripts/THOR_WTvsHET.sh # 1141742 bigwig FAIL; output path changed to output/THOR: 1254973 ok
sbatch scripts/THOR_WTvsHET_unique.sh # 1141743 bigwig FAIL
sbatch scripts/THOR_WTvsKO.sh # 1254977 ok

# Light test pvalue on WTvsHET
sbatch scripts/THOR_WTvsHETpval0.05.sh # 1254985 ok
sbatch scripts/THOR_WTvsHETpval0.01.sh # 1254987 ok
sbatch scripts/THOR_WTvsHETpval0.005.sh # 1256847 ok
sbatch scripts/THOR_WTvsHETpval0.001.sh # 1256848 ok
sbatch scripts/THOR_WTvsHETpval0.0005.sh # 1256849 ok
sbatch scripts/THOR_WTvsHETpval0.0001.sh # 1256850 ok
sbatch scripts/THOR_WTvsHETpval0.00001.sh # 1307914 ok
sbatch scripts/THOR_WTvsHETpval0.000001.sh # 1307915 ok
sbatch scripts/THOR_WTvsHETpval0.0000001.sh # 1307916 ok
sbatch scripts/THOR_WTvsHETpval0.00000001.sh # 1307917 ok


# Light test when using non-DiffBind TMM normalized scaling factors (bigwig_histone/ or bigwig_histone_):
sbatch scripts/THOR_WTvsHET_bigwig_histone_NotGenotypeGroupSF.sh # 1256830

# Light test for binsize and step with default pvalue 0.1
sbatch scripts/THOR_WTvsHETbinsize250.sh # 1256945 ok
sbatch scripts/THOR_WTvsHETbinsize500.sh # 1256946 ok
sbatch scripts/THOR_WTvsHETbinsize1000.sh # 1256947 ok
sbatch scripts/THOR_WTvsHETbinsize1500.sh # 1256948 ok


# Light test with poisson distribution (`--poisson`)
sbatch scripts/THOR_WTvsHETpval0.0001_poisson.sh  # 1307944 ok

# Poisson distribution default parameter
sbatch scripts/THOR_WTvsHET_poisson.sh # 1308030
sbatch scripts/THOR_WTvsKO_poisson.sh # 1308031

# Generate median bigwig files for default and poisson
conda activate BedToBigwig
sbatch scripts/bigwigmerge_THOR_WTvsHET.sh # 1308094
sbatch --dependency=afterany:1308030 scripts/bigwigmerge_THOR_WTvsHET_poisson.sh # 1308095
sbatch scripts/bigwigmerge_THOR_WTvsKO.sh # 1308096
sbatch --dependency=afterany:1308031 scripts/bigwigmerge_THOR_WTvsKO_poisson.sh # 1308097
```
- *NOTE: Options:  `-merge` option recommended for histone data. `–report` for HTML report (useless!), not super important; just to see how it look; `–deadzones` is blacklist; `-pvalue` 0.1 is default (can play with it);*
- *NOTE: do NOT put any "_" in the `--name` !! Or bug*
- *NOTE: unique is remove dupplicated reads with `--rmdup` option*
- *NOTE: got Ressource-related warnings at the end; after bigWigs processing; that is OK!*
- *NOTE: Let's re-run with the correct conda environment so that it can output the bigwig files (**output will now be within the` output/THOR` folder**)*
- *NOTE: I keep `binsize = 2*step` in the parameters; not sure optimal... But the default follow this.*
- *NOTE: the bigwig files are ALWAYS the same whatever the qvalue treshold used*
- *NOTE: bigwig files WT are very similar; nearly identical, when WT from WTvsHET than WTvsKO*

--> much LESS peaks are identified using uniquely aligned reads (`--rmdup` option)

--> Overall a LOT of diff peaks are identified, and looking at the IGV with `bigwig_DiffBind_TMM` it looks good !!! There may even be too many (114,886 diff peaks!!); so **we could filter qvalue**

--> pvalue from 0.1 to 0.0001 do not change a lot; then XXX. pvalue of XXX is optimal. In the end; **I could have just filter the qvalue in the bed output file**... (qvalue is in MACS2-like format)
----> 

--> increasing `--binsize` and `--step` more than default value result in MANY more false positive peaks; size overall not super longer in addition! --> **NOT GOOD; keep default**

--> using `output/bigwig_histone_NotGenotypeGroup` SF = non-DiffBind TMM normalized scaling factors show A LOT LESS diff bound regions (4.2Mo vs 900ko at same default qvalue) (in agreement with the deepTools profile showing lot of disparities between replicates when using these normalization factors) --> **NOT GOOD**

--> The bigwig files generated looks good (**when using DiffBind_TMM SF**); replicate are comparable

--> playing with binsize and step (default is 100, 50) is XXX

Let's instead of using different pvalue to call for peak, put the raw bed output in R and filter p-value here + generate FC (found [here](http://ginolhac.github.io/chip-seq/peak/)):

XXXX

```R
# load the file using the tidyverse
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
diffpeaks <- read_tsv("TC1-I-A-D0vsD3-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
# split the last field into three
diffpeaks %>%
  separate(X11, into = c("count1", "count2", "third"), sep = ";", convert = TRUE) %>%
  mutate(FC = count2 / count1) -> thor_splitted
# plot the histogram of the fold-change computed above, count second condition / count 1st condition
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1))

# create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(log2(FC) > 0.5) %>%
  select(X1, X2, X3) %>%
  mutate(X1 = paste0("chr", X1)) %>%
  write_tsv("THOR_logFC0.5.bed", col_names = FALSE)

```


**Check on IGV how it look with different qvalue; FC treshold**

--> Optimal qvalue and FC treshold seems to be XXX

--> Poisson seems to perform XXX


## Histone-EpiCypher guidelines for scaling normalization

Here is the guidelines from the [website](https://www.epicypher.com/content/documents/SNAP-CUTANA_K-MetStat_user_guide.pdf) *For spike-in normalization, a scale factor was calculated for each sample by dividing the percent of total reads aligned to human genome by the percent of total reads aligned to the spike-in barcodes (Scale Factor = % Human Reads / % Spike-in Reads) and applying this factor to adjust the total sequencing reads of each respective sample*.

- In `spikein/spikein_histone_groupABgenotype_scaling_factor.txt` or `spikein/spikein_histone_H3K27me3_scaling_factor_fastp.txt` use align_reads or total = nb of reads aligned to spike in
- In `/home/roulet/001_EZH1_project/003__CutRun/mapping_QC.xlsx` = nb of reads aligned to human genome and total number

Let's make an excell file to collect all the counts in `output/Epicypher/read_counts_human_histone.xlsx`

Now let's apply the scaling factor in Bam to bigwig and output to `output/bigwig_Epicypher`

```bash
sbatch scripts/bamtobigwig_Epicypher.sh # 12336358; output was same re-run 12341847 ok
sbatch scripts/bamtobigwig_Epicypher_unique.sh # 12336357; output was same re-run 12341850 ok
```

--> It seems too exagerated, like when there is low histone; signal is super enriched.






# PCA on Bigwig files
Let's do PCA with [multiBigwigSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html) and [PCAplot](https://deeptools.readthedocs.io/en/2.4.1/content/tools/plotPCA.html)/[plotCorrelation](https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html) from deeptools:

- Use bin mode (score of the compile bigwig is calculated at a 10kb bins (default), on the entire genome)
- Let's do it all samples together and per genotype for better vizualization

## PCA on raw files

```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch --dependency=afterany:12092810 scripts/multiBigwigSummary_all.sh # 12092852 ok
sbatch --dependency=afterany:12092810 scripts/multiBigwigSummary_H3K27me3.sh # 12092853 ok

# Plot
## All genotypes all points
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels 8wN_WT_IGG_R1 8wN_WT_IGG_R2 8wN_WT_IGG_R3 8wN_WT_IGG_R4 8wN_WT_H3K27me3_R1 8wN_WT_H3K27me3_R2 8wN_WT_H3K27me3_R3 8wN_WT_H3K27me3_R4 8wN_KO_IGG_R1 8wN_KO_IGG_R2 8wN_KO_IGG_R3 8wN_KO_IGG_R4 8wN_KO_H3K27me3_R1 8wN_KO_H3K27me3_R2 8wN_KO_H3K27me3_R3 8wN_KO_H3K27me3_R4 8wN_HET_IGG_R1 8wN_HET_IGG_R2 8wN_HET_IGG_R3 8wN_HET_IGG_R4 8wN_HET_H3K27me3_R1 8wN_HET_H3K27me3_R2 8wN_HET_H3K27me3_R3 8wN_HET_H3K27me3_R4 8wN_iPSCpatient_IGG_R1 8wN_iPSCpatient_IGG_R2 8wN_iPSCpatient_H3K27me3_R1 8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf
plotPCA -in output/bigwig/multiBigwigSummary_H3K27me3.npz \
    --transpose \
    --ntop 0 \
    --labels 8wN_WT_H3K27me3_R1 8wN_WT_H3K27me3_R2 8wN_WT_H3K27me3_R3 8wN_WT_H3K27me3_R4 8wN_KO_H3K27me3_R1 8wN_KO_H3K27me3_R2 8wN_KO_H3K27me3_R3 8wN_KO_H3K27me3_R4 8wN_HET_H3K27me3_R1 8wN_HET_H3K27me3_R2 8wN_HET_H3K27me3_R3 8wN_HET_H3K27me3_R4 8wN_iPSCpatient_H3K27me3_R1 8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    -o output/bigwig/multiBigwigSummary_H3K27me3_plotPCA.pdf

plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_H3K27me3.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_H3K27me3_heatmap.pdf
```

## PCA on normalized files
### Histone
The code below has not been run; I think useless to do it as it will not show any differences with the raw one...
```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_histone_all.sh # 
sbatch scripts/multiBigwigSummary_histone_H3K27me3.sh # 

# Plot
## All genotypes all points
plotPCA -in output/bigwig/multiBigwigSummary_histone_all.npz \
    --transpose \
    --ntop 0 \
    --labels 8wN_WT_IGG_R1 8wN_WT_IGG_R2 8wN_WT_IGG_R3 8wN_WT_IGG_R4 8wN_WT_H3K27me3_R1 8wN_WT_H3K27me3_R2 8wN_WT_H3K27me3_R3 8wN_WT_H3K27me3_R4 8wN_KO_IGG_R1 8wN_KO_IGG_R2 8wN_KO_IGG_R3 8wN_KO_IGG_R4 8wN_KO_H3K27me3_R1 8wN_KO_H3K27me3_R2 8wN_KO_H3K27me3_R3 8wN_KO_H3K27me3_R4 8wN_HET_IGG_R1 8wN_HET_IGG_R2 8wN_HET_IGG_R3 8wN_HET_IGG_R4 8wN_HET_H3K27me3_R1 8wN_HET_H3K27me3_R2 8wN_HET_H3K27me3_R3 8wN_HET_H3K27me3_R4 8wN_iPSCpatient_IGG_R1 8wN_iPSCpatient_IGG_R2 8wN_iPSCpatient_H3K27me3_R1 8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    -o output/bigwig/multiBigwigSummary_histone_all_plotPCA.pdf
plotPCA -in output/bigwig/multiBigwigSummary_histone_H3K27me3.npz \
    --transpose \
    --ntop 0 \
    --labels 8wN_WT_H3K27me3_R1 8wN_WT_H3K27me3_R2 8wN_WT_H3K27me3_R3 8wN_WT_H3K27me3_R4 8wN_KO_H3K27me3_R1 8wN_KO_H3K27me3_R2 8wN_KO_H3K27me3_R3 8wN_KO_H3K27me3_R4 8wN_HET_H3K27me3_R1 8wN_HET_H3K27me3_R2 8wN_HET_H3K27me3_R3 8wN_HET_H3K27me3_R4 8wN_iPSCpatient_H3K27me3_R1 8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    -o output/bigwig/multiBigwigSummary_histone_H3K27me3_plotPCA.pdf

plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_histone_all_heatmap.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_H3K27me3.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_histone_H3K27me3_heatmap.pdf
```


Let's do **plotCorrelation clustering** with our **bigwig_histone** scaled and **bigwig_DiffBind_TMM** scaled files (H3K27me3 only and no patient):

```bash
conda activate deeptools

# Generate compile bigwig (.npz) files and plots
sbatch scripts/multiBigwigSummary_bigwig_histone.sh # 61437 ok
sbatch scripts/multiBigwigSummary_bigwig_DiffBind_TMM.sh # 65646 ok
```

--> both files are exactly the same...




# ChIPseqSpikeInFree on CutRun

Let's see what is the prediction from ChIPseqSpikeInFree on our data

```bash
conda activate ChIPseqSpikeInFree
sbatch scripts/ChIPseqSpikeInFree.sh # 12379092
```

It does not work super well ("ok" correlated: r=0.4, p0.1). Possibly because CutRun is too clean, very few low proportion reads as compare to ChIPseq so may fail...



# ChIPseeker for binding profiles
[Tutorial](http://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html) and [documentation]()


## Pool the peak into 1 file
Let's use IDR to pool our replicates and identify 'high-confidence' peak (used in the ENCODE pipeline)

[Tutorial](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html) and [github](https://github.com/nboley/idr) for IDR.

**Installation:**
```bash
conda create -n idr -c bioconda idr # command can be launch from anywhere (directory and node)
conda activate idr
```
**Generate peak files:**
qvalue 2.3 (0.005) was optimal for CutRun:

```bash
conda activate idr

idr --samples output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_R1_peaks.broadPeak \
    output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_R2_peaks.broadPeak \
    output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_R3_peaks.broadPeak \
    output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_R4_peaks.broadPeak \
    --input-file-type broadPeak \
    --output-file output/idr/8wN_WT_H3K27me3_idr \
    --plot \
    --log-output-file output/idr/8wN_WT_H3K27me3_idr.log
```

IDR can only handle 2 replicates... We can either do 2 per 2 comparison and keep merge (see [here](https://github.com/nboley/idr/issues/35)). Or use another method...:


Re-run MACS2 but in pool to have 1 file per condition (keep blacklist and qval2.30103 (0.005) filtering):
```bash
sbatch scripts/macs2_pool.sh # 12447366 ok
```
*NOTE: for patient, I only take Rep1 as Rep2 failed, no enrichment*

Now let's filter out blacklist and qvalue:
```bash
sbatch scripts/macs2_pool_peak_signif.sh # ok

```
*NOTE: I relaunch the script and changed the qvalue; 2.30103 (q0.005) is good as observed looking at individual replicates*





## Run ChIPseeker
ChIPseeker docu [here for plotAvgProf and heatmap around TSS](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/chipseeker_visualization.html) and [here more information plot and functional analyses](https://bioconductor.statistik.tu-dortmund.de/packages/3.7/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html) and [here for peak assignment](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html)

For **ChIPseeker**, use **conda base and R 4.2.2 module**
```bash
conda deactivate
module load R/4.2.2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")


# Import peaks
peaks_WT =  read.table('output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_pool_peaks.broadPeak') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9) 
peaks_KO =  read.table('output/macs2/broad_blacklist_qval2.30103/8wN_KO_H3K27me3_pool_peaks.broadPeak') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9) 
peaks_HET =  read.table('output/macs2/broad_blacklist_qval2.30103/8wN_HET_H3K27me3_pool_peaks.broadPeak') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9) 
peaks_patient =  read.table('output/macs2/broad_blacklist_qval2.30103/8wN_iPSCpatient_H3K27me3_pool_peaks.broadPeak') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9) 


# Tidy peaks
WT_gr = makeGRangesFromDataFrame(peaks_WT,keep.extra.columns=TRUE)
KO_gr = makeGRangesFromDataFrame(peaks_KO,keep.extra.columns=TRUE)
HET_gr = makeGRangesFromDataFrame(peaks_HET,keep.extra.columns=TRUE)
patient_gr = makeGRangesFromDataFrame(peaks_patient,keep.extra.columns=TRUE)

gr_list <- list(WT=WT_gr, KO=KO_gr, HET=HET_gr, patient=patient_gr)

# Isolate TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=10000, downstream=10000) # region around TSS

# Calculate the tag matrix
tagMatrixList <- lapply(gr_list, getTagMatrix, windows=promoter)


# Vizualization
## Profile plots around TSS
### Profile multiple plot
pdf("output/ChIPseeker/profile_TSS_10kb.pdf", width=14, height=20)
plotAvgProf(tagMatrixList, xlim=c(-10000, 10000), conf=0.95,resample=200, facet="row") # last 10-15min
dev.off()
### Profile one plot
pdf("output/ChIPseeker/profile_oneplot_TSS_10kb.pdf", width=14, height=14)
plotAvgProf(tagMatrixList, xlim=c(-10000, 10000))
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=50000, downstream=50000) # region around TSS
tagMatrixList <- lapply(gr_list, getTagMatrix, windows=promoter)
pdf("output/ChIPseeker/profile_oneplot_TSS_50kb.pdf", width=14, height=14)
plotAvgProf(tagMatrixList, xlim=c(-50000, 50000))
dev.off()
### Heatmap
pdf("output/ChIPseeker/heatmap_TSS_10kb.pdf", width=14, height=20)
tagHeatmap(tagMatrixList, xlim=c(-10000, 10000), color=NULL)
dev.off()
pdf("output/ChIPseeker/heatmap_TSS_50kb.pdf", width=14, height=20) # last 30min + failed with 150g memory; work with 400g 
tagHeatmap(tagMatrixList, xlim=c(-50000, 50000), color=NULL)
dev.off()


## Profile plots around gene body
### Profile
pdf("output/ChIPseeker/profile_gene_body.pdf", width=14, height=20)
plotPeakProf2(peak = gr_list, upstream = rel(0.2), downstream = rel(0.2),    # last 10-15min
             conf = 0.95, by = "gene", type = "body", nbin = 100,
             TxDb = txdb, ignore_strand = F) # nbin = filter out peak smaller than 100
dev.off()

## Genomic Annotation ALL TOGETHER
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

### Barplot
pdf("output/ChIPseeker/annotation_barplot.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()
### Barplot
plotDistToTSS()

## Peak distance to TSS
pdf("output/ChIPseeker/DistToTSS.pdf", width=14, height=5)
plotDistToTSS(peakAnnoList)
dev.off()

## Genomic Annotation ONE BY ONE
peakAnno_WT <- annotatePeak(WT_gr, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_KO <- annotatePeak(KO_gr, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_HET <- annotatePeak(HET_gr, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_patient <- annotatePeak(patient_gr, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
### Pie
pdf("output/ChIPseeker/annotation_pie_WT.pdf", width=14, height=20)
plotAnnoPie(peakAnno_WT)
dev.off()
pdf("output/ChIPseeker/annotation_pie_KO.pdf", width=14, height=20)
plotAnnoPie(peakAnno_KO)
dev.off()
pdf("output/ChIPseeker/annotation_pie_HET.pdf", width=14, height=20)
plotAnnoPie(peakAnno_HET)
dev.off()
pdf("output/ChIPseeker/annotation_pie_patient.pdf", width=14, height=20)
plotAnnoPie(peakAnno_patient)
dev.off()

### Barplot
pdf("output/ChIPseeker/annotation_barplot_WT.pdf", width=14, height=5)
plotAnnoBar(peakAnno_WT)
dev.off()
pdf("output/ChIPseeker/annotation_barplot_KO.pdf", width=14, height=5)
plotAnnoBar(peakAnno_KO)
dev.off()
pdf("output/ChIPseeker/annotation_barplot_HET.pdf", width=14, height=5)
plotAnnoBar(peakAnno_HET)
dev.off()
pdf("output/ChIPseeker/annotation_barplot_patient.pdf", width=14, height=5)
plotAnnoBar(peakAnno_patient)
dev.off()

### Vennpie
pdf("output/ChIPseeker/annotation_vennpie_WT.pdf", width=14, height=20)
vennpie(peakAnno_WT)
dev.off()
pdf("output/ChIPseeker/annotation_vennpie_KO.pdf", width=14, height=20)
vennpie(peakAnno_KO)
dev.off()
pdf("output/ChIPseeker/annotation_vennpie_HET.pdf", width=14, height=20)
vennpie(peakAnno_HET)
dev.off()
pdf("output/ChIPseeker/annotation_vennpie_patient.pdf", width=14, height=20)
vennpie(peakAnno_patient)
dev.off()

### Upsetplot # Need upset plot package, not dramatic weird looking...
pdf("output/ChIPseeker/annotation_upsetplot_WT.pdf", width=14, height=20)
upsetplot(peakAnno_WT, vennpie=TRUE)
dev.off()


## Peak distance to TSS
pdf("output/ChIPseeker/DistToTSS_WT.pdf", width=14, height=20)
plotDistToTSS(peakAnno_WT,
              title="Distribution of H3K27me3\nrelative to TSS")
dev.off()
pdf("output/ChIPseeker/DistToTSS_KO.pdf", width=14, height=20)
plotDistToTSS(peakAnno_KO,
              title="Distribution of H3K27me3\nrelative to TSS")
dev.off()
pdf("output/ChIPseeker/DistToTSS_HET.pdf", width=14, height=20)
plotDistToTSS(peakAnno_HET,
              title="Distribution of H3K27me3\nrelative to TSS")
dev.off()
pdf("output/ChIPseeker/DistToTSS_patient.pdf", width=14, height=20)
plotDistToTSS(peakAnno_patient,
              title="Distribution of H3K27me3\nrelative to TSS")
dev.off()


# Functional profiles
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))

## KEGG
compKEGG <- compareCluster(geneCluster   = genes,
                         fun           = "enrichKEGG",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/ChIPseeker/functional_KEGG.pdf", width=7, height=15)
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()

## GO
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_BP.pdf", width=7, height=15)
dotplot(compGO, showCategory = 15, title = "GO_Biological Process Enrichment Analysis")
dev.off()
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "CC") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_CC.pdf", width=7, height=15)
dotplot(compGO, showCategory = 15, title = "GO_Cellular Component Enrichment Analysis")
dev.off()
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "MF") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_MF.pdf", width=7, height=15)
dotplot(compGO, showCategory = 15, title = "GO_Molecular Function Enrichment Analysis")
dev.off()


## Disease
compDO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichDO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
pdf("output/ChIPseeker/functional_DO.pdf", width=7, height=15)
dotplot(compDO, showCategory = 15, title = "Disease Ontology Enrichment Analysis")
dev.off()

## ReactomePA
compReactome <- compareCluster(geneCluster   = genes,
                               fun           = "enrichPathway",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH",
                               organism      = "human")
pdf("output/ChIPseeker/functional_Reactome.pdf", width=7, height=15)
dotplot(compReactome, showCategory = 15, title = "Reactome Pathway Enrichment Analysis")
dev.off()

# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes.pdf", width=7, height=7)
vennplot(genes)
dev.off()

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
WT_annot <- as.data.frame(peakAnnoList[["WT"]]@anno)
KO_annot <- as.data.frame(peakAnnoList[["KO"]]@anno)
HET_annot <- as.data.frame(peakAnnoList[["HET"]]@anno)
patient_annot <- as.data.frame(peakAnnoList[["patient"]]@anno)

## Convert entrez gene IDs to gene symbols
WT_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = WT_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KO_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KO_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
HET_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = HET_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
patient_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = patient_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")

WT_annot$gene <- mapIds(org.Hs.eg.db, keys = WT_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KO_annot$gene <- mapIds(org.Hs.eg.db, keys = KO_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
HET_annot$gene <- mapIds(org.Hs.eg.db, keys = HET_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
patient_annot$gene <- mapIds(org.Hs.eg.db, keys = patient_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")




## Save output table
write.table(WT_annot, file="output/ChIPseeker/annotation_WT.txt", sep="\t", quote=F, row.names=F)
write.table(KO_annot, file="output/ChIPseeker/annotation_KO.txt", sep="\t", quote=F, row.names=F)
write.table(HET_annot, file="output/ChIPseeker/annotation_HET.txt", sep="\t", quote=F, row.names=F)
write.table(patient_annot, file="output/ChIPseeker/annotation_patient.txt", sep="\t", quote=F, row.names=F)
```
*NOTE: the `resample` in `plotAvgProf()` is how many time you bootstrap; nb need to be a divisor of the total nb of datapoint (eg. 100/200 for 20k datapoints (=-1000/1000) is ok, but 500 is not)* \
*NOTE: the gene peak assignment use the nearest TSS method, then other order of priority, all genes are assigned except if they do not fall into Promoter, 5’ UTR, 3’ UTR, Exon, Intron, Downstream. In such case peak assign as 'intergenic'*


Let's assign the **differentially bound peak to genes**:
```bash
srun --mem=300g --pty bash -l
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


# Import peaks
peaks_contrast1 = read.table('output/DiffBind/sample_count_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast1_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast2 = read.table('output/DiffBind/sample_count_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast2_df.bed')  %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast3 = read.table('output/DiffBind/sample_count_macs2raw_blackgreylist_LibHistoneScaled_TMM_contrast3_df.bed')  %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11)  


# Tidy peaks
peaks_contrast1_gr = makeGRangesFromDataFrame(peaks_contrast1,keep.extra.columns=TRUE)
peaks_contrast2_gr = makeGRangesFromDataFrame(peaks_contrast2,keep.extra.columns=TRUE)
peaks_contrast3_gr = makeGRangesFromDataFrame(peaks_contrast3,keep.extra.columns=TRUE)

gr_list <- list(HETvsKO=peaks_contrast1_gr, HETvsWT=peaks_contrast2_gr, KOvsWT=peaks_contrast3_gr)



# Functional profiles
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # 

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))

## KEGG
compKEGG <- compareCluster(geneCluster   = genes,
                         fun           = "enrichKEGG",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/ChIPseeker/functional_KEGG_DiffBind05.pdf", width=7, height=15)
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()

## GO
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_BP_DiffBind05.pdf", width=7, height=15)
dotplot(compGO, showCategory = 15, title = "GO_Biological Process Enrichment Analysis")
dev.off()
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "CC") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_CC_DiffBind05.pdf", width=7, height=15)
dotplot(compGO, showCategory = 15, title = "GO_Cellular Component Enrichment Analysis")
dev.off()
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "MF") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_MF_DiffBind05.pdf", width=7, height=15)
dotplot(compGO, showCategory = 15, title = "GO_Molecular Function Enrichment Analysis")
dev.off()


## Disease
compDO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichDO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
pdf("output/ChIPseeker/functional_DO_DiffBind05.pdf", width=7, height=15)
dotplot(compDO, showCategory = 15, title = "Disease Ontology Enrichment Analysis")
dev.off()

## ReactomePA
compReactome <- compareCluster(geneCluster   = genes,
                               fun           = "enrichPathway",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH",
                               organism      = "human")
pdf("output/ChIPseeker/functional_Reactome_DiffBind05.pdf", width=7, height=15)
dotplot(compReactome, showCategory = 15, title = "Reactome Pathway Enrichment Analysis")
dev.off()

# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_DiffBind05.pdf", width=7, height=7)
vennplot(genes)
dev.off()

# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
HETvsKO_annot <- as.data.frame(peakAnnoList[["HETvsKO"]]@anno)
HETvsWT_annot <- as.data.frame(peakAnnoList[["HETvsWT"]]@anno)
KOvsWT_annot <- as.data.frame(peakAnnoList[["KOvsWT"]]@anno)


## Convert entrez gene IDs to gene symbols
HETvsKO_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = HETvsKO_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
HETvsWT_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = HETvsWT_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
KOvsWT_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = KOvsWT_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")

HETvsKO_annot$gene <- mapIds(org.Hs.eg.db, keys = HETvsKO_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
HETvsWT_annot$gene <- mapIds(org.Hs.eg.db, keys = HETvsWT_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
KOvsWT_annot$gene <- mapIds(org.Hs.eg.db, keys = KOvsWT_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(HETvsKO_annot, file="output/ChIPseeker/annotation_HETvsKO.txt", sep="\t", quote=F, row.names=F)
write.table(HETvsWT_annot, file="output/ChIPseeker/annotation_HETvsWT.txt", sep="\t", quote=F, row.names=F)
write.table(KOvsWT_annot, file="output/ChIPseeker/annotation_KOvsWT.txt", sep="\t", quote=F, row.names=F)
```

Now let's compare RNAseq and CutRun for MACS2raw qval 0.05:
- Filter HETvsWT and KOvsWT diff bound genes into **gain and loss H3K27me3**
- **Keep only signal in Promoter, gene body and TES** (ie. filter out peak assigned to intergenic)
- **Merge with deseq2** log2FC data (tpm will not work as too variable; or log2tpm maybe?)
- Plot in x FC and y baseMean=deseq2-norm counts (+ color qvalue) with facet_wrap~gain or lost (ie. volcano plot gain/lost)

Plot functional analysis for the genes that gain HET and decreased in expression (NaiaraPlot)

```R
# lib
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")

# load CutRun gainLoss table
HETvsKO_annot <- read.table("../003__CutRun/output/ChIPseeker/annotation_HETvsKO.txt", sep="\t", header=TRUE)
HETvsWT_annot <- read.table("../003__CutRun/output/ChIPseeker/annotation_HETvsWT.txt", sep="\t", header=TRUE)
KOvsWT_annot <- read.table("../003__CutRun/output/ChIPseeker/annotation_KOvsWT.txt", sep="\t", header=TRUE)

# Filter Gain/Loss sites
HETvsWT_annot_gain = tibble(HETvsWT_annot) %>%
    filter(FC > 0, annotation != "Distal Intergenic") %>%
    dplyr::select(-c(V6,V7,V8,strand,name)) %>%
    add_column(H3K27me3 = "gain")
HETvsWT_annot_lost = tibble(HETvsWT_annot) %>%
    filter(FC < 0, annotation != "Distal Intergenic") %>%
    dplyr::select(-c(V6,V7,V8,strand,name)) %>%
    add_column(H3K27me3 = "lost")
HETvsWT_annot_gain_lost = HETvsWT_annot_gain %>% 
    bind_rows(HETvsWT_annot_lost) %>%
    dplyr::select(seqnames,start,end,FC,pvalue,FDR,annotation,distanceToTSS,gene,geneSymbol,H3K27me3)

KOvsWT_annot_gain = tibble(KOvsWT_annot) %>%
    filter(FC > 0, annotation != "Distal Intergenic") %>%
    dplyr::select(-c(V6,V7,V8,strand,name)) %>%
    add_column(H3K27me3 = "gain")
KOvsWT_annot_lost = tibble(KOvsWT_annot) %>%
    filter(FC < 0, annotation != "Distal Intergenic") %>%
    dplyr::select(-c(V6,V7,V8,strand,name))  %>%
    add_column(H3K27me3 = "lost")
KOvsWT_annot_gain_lost = KOvsWT_annot_gain %>% 
    bind_rows(KOvsWT_annot_lost) %>%
    dplyr::select(seqnames,start,end,FC,pvalue,FDR,annotation,distanceToTSS,gene,geneSymbol,H3K27me3)

# Import RNAseq deseq2 output
HET_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_8wN_HET_vs_8wN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)
KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_8wN_KO_vs_8wN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

# Merge files
HETvsWT_annot_gain_lost_RNA = HETvsWT_annot_gain_lost %>% 
    left_join(HET_vs_WT) %>%
    dplyr::select(gene, H3K27me3,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()


KOvsWT_annot_gain_lost_RNA = KOvsWT_annot_gain_lost %>% 
    left_join(KO_vs_WT) %>%
    dplyr::select(gene, H3K27me3,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()


# Volcano plot
count_data <- HETvsWT_annot_gain_lost_RNA %>%
    group_by(H3K27me3, significance) %>%
    summarise(up = sum(log2FoldChange > 0),
              down = sum(log2FoldChange < 0),
              total = n()) %>%
    ungroup() %>%
    group_by(H3K27me3) %>%
    mutate(total_panel = sum(total)) %>%
    ungroup()

pdf("output/ChIPseeker/DiffBind05_HETvsWT_expression.pdf", width=7, height=4)
HETvsWT_annot_gain_lost_RNA %>%
    ggplot(aes(x = log2FoldChange, y = baseMean, color = significance)) +
        geom_point(alpha = 0.8, size = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        labs(title = "HET vs WT",
             subtitle = "Expression level of diff. bound H3K27me3 genes",
             x = "Log2 Fold Change",
             y = "Base Mean",
             color = "Significant (padj <= 0.05)") +
        facet_wrap(~H3K27me3) +
        theme_bw() +
        geom_text(data = count_data %>% filter(significance), 
                  aes(x = Inf, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
                  hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
        geom_text(data = count_data %>% distinct(H3K27me3, .keep_all = TRUE),
                  aes(x = Inf, y = -Inf, label = paste("Total:", total_panel, "genes")),
                  hjust = 1.1, vjust = -0.1, size = 3, color = "black")
dev.off()


pdf("output/ChIPseeker/DiffBind05_HETvsWT_expression_qvalue.pdf", width=7, height=4)
HETvsWT_annot_gain_lost_RNA %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
        geom_point(alpha = 0.8, size = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        labs(title = "HET vs WT",
             subtitle = "Expression level of diff. bound H3K27me3 genes",
             x = "Log2 Fold Change",
             y = "-log10(q-value)",
             color = "Significant (padj <= 0.05)") +
        facet_wrap(~H3K27me3) +
        theme_bw() +
        geom_text(data = count_data %>% filter(significance), 
                  aes(x = Inf, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
                  hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
        geom_text(data = count_data %>% distinct(H3K27me3, .keep_all = TRUE),
                  aes(x = Inf, y = -Inf, label = paste("Total:", total_panel, "genes")),
                  hjust = 1.1, vjust = -0.1, size = 3, color = "black")
dev.off()



count_data <- KOvsWT_annot_gain_lost_RNA %>%
    group_by(H3K27me3, significance) %>%
    summarise(up = sum(log2FoldChange > 0),
              down = sum(log2FoldChange < 0),
              total = n()) %>%
    ungroup() %>%
    group_by(H3K27me3) %>%
    mutate(total_panel = sum(total)) %>%
    ungroup()

pdf("output/ChIPseeker/DiffBind05_KOvsWT_expression.pdf", width=7, height=4)
KOvsWT_annot_gain_lost_RNA %>%
    ggplot(aes(x = log2FoldChange, y = baseMean, color = significance)) +
        geom_point(alpha = 0.8, size = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        labs(title = "KO vs WT",
             subtitle = "Expression level of diff. bound H3K27me3 genes",
             x = "Log2 Fold Change",
             y = "Base Mean",
             color = "Significant (padj <= 0.05)") +
        facet_wrap(~H3K27me3) +
        theme_bw() +
        geom_text(data = count_data %>% filter(significance), 
                  aes(x = Inf, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
                  hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
        geom_text(data = count_data %>% distinct(H3K27me3, .keep_all = TRUE),
                  aes(x = Inf, y = -Inf, label = paste("Total:", total_panel, "genes")),
                  hjust = 1.1, vjust = -0.1, size = 3, color = "black")
dev.off()

# Venn diagram to see whether gain H3K27me3 in HET are the same one in KO?
## GAIN
# Extract the unique gene lists as vectors
HET_gene_list <- unique(HETvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "gain") %>%
                          pull(gene))

KO_gene_list <- unique(KOvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "gain") %>%
                          pull(gene))


venn.plot <- venn.diagram(
  list(HET = HET_gene_list, KO = KO_gene_list),
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  category.names = c("HET Gain", "KO Gain"),
  cex = 2
)

# Plot the Venn diagram
pdf("output/ChIPseeker/Venn_DiffBind05_Gain_Het_KO.pdf", width=5, height=4)
grid.draw(venn.plot)
dev.off()

## LOST
# Extract the unique gene lists as vectors
HET_gene_list <- unique(HETvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "lost") %>%
                          pull(gene))

KO_gene_list <- unique(KOvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "lost") %>%
                          pull(gene))


venn.plot <- venn.diagram(
  list(HET = HET_gene_list, KO = KO_gene_list),
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  category.names = c("HET Lost", "KO Lost"),
  cex = 2
)

# Plot the Venn diagram
pdf("output/ChIPseeker/Venn_DiffBind05_Lost_Het_KO.pdf", width=5, height=4)
grid.draw(venn.plot)
dev.off()


## BOTH
# Extract the unique gene lists as vectors
HET_lost_gene_list <- unique(HETvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "lost") %>%
                          pull(gene))

KO_lost_gene_list <- unique(KOvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "lost") %>%
                          pull(gene))

HET_gain_gene_list <- unique(HETvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "gain") %>%
                          pull(gene))

KO_gain_gene_list <- unique(KOvsWT_annot_gain_lost_RNA %>%
                          filter(H3K27me3 == "gain") %>%
                          pull(gene))


venn.plot <- venn.diagram(
  list(HET_lost = HET_lost_gene_list, KO_lost = KO_lost_gene_list,
       HET_gain = HET_gain_gene_list, KO_gain = KO_gain_gene_list),
  filename = NULL,
  fill = c("blue", "red", "blue", "red"),
  alpha = 0.5,
  cex = 2
)

# Plot the Venn diagram
pdf("output/ChIPseeker/Venn_DiffBind05_LostGain_Het_KO.pdf", width=5, height=4)
grid.draw(venn.plot)
dev.off()


# Functional analyses
## For the WT vs HET only (NaiaraPlot)

HETvsWT_annot_gain_RNA = HETvsWT_annot_gain_lost_RNA %>%
    filter(H3K27me3 == "gain", 
           log2FoldChange < 0,
           significance == "TRUE") 

HETvsWT_annot_gain_RNA = HETvsWT_annot_gain_lost_RNA %>%
    filter(H3K27me3 == "gain", 
           log2FoldChange < 0)

## Read GTF file
gtf_file <- "../../Master/meta/ENCFF159KBI.gtf"
gtf_data <- import(gtf_file)

## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name

## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble() %>%
  separate(gene_id, into = c("gene_id", "trash"), sep = "\\.", remove = TRUE) %>%
  dplyr::select(-trash)

## Import genes_cluster list and background list
HETvsWT_annot_gain_RNA_gene = HETvsWT_annot_gain_RNA %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  left_join(gene_id_name) %>%
  dplyr::select(gene_name) 



background = read_csv("../001__RNAseq/output/deseq2_hg38/raw_2dN_HET_vs_2dN_WT.txt") %>%
  dplyr::select(gene) %>%
  rename(gene_id = gene) %>%
  inner_join(gene_id_name) %>%
  dplyr::select(gene_name)


## Run GO enrichment analysis 
ego <- enrichGO(gene = as.character(HETvsWT_annot_gain_RNA_gene$gene_name), 
                universe = as.character(background$gene_name),
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "BP",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   # can do FDR
                pvalueCutoff = 0.05, 
                readable = TRUE)

ego <- enrichGO(gene   = as.character(HETvsWT_annot_gain_RNA_gene$gene_name),
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         keyType = "SYMBOL",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP")

## Save GO analyses
GO_summary <- data.frame(ego)

write.csv(GO_summary, "output/GO_hg38/cluster4_BP.csv")
write.csv(GO_summary, "output/GO_hg38/cluster4_MF.csv")
write.csv(GO_summary, "output/GO_hg38/cluster4_CC.csv")



## Vizualization

pdf("output/ChIPseeker/functional_GO_dotplot_BP_HETvsWT_annot_gain_RNA.pdf", width=7, height=5)
dotplot(ego, showCategory=50)
dev.off()
pdf("output/ChIPseeker/functional_GO_emmaplot_BP_HETvsWT_annot_gain_RNA.pdf", width=8, height=10)
emapplot(pairwise_termsim(ego), showCategory = 50)
dev.off()



## Create a list of vector:
gene_list <- list(HET_lost = HET_lost_gene_list,
                  KO_lost = KO_lost_gene_list,
                  HET_gain = HET_gain_gene_list,
                  KO_gain = KO_gain_gene_list)

## Convert ENSEMBL geneID to Entrez gene ID:
convert_ensembl_to_entrez <- function(gene_list) {
  result <- select(org.Hs.eg.db, keys = gene_list, keytype = "ENSEMBL", columns = "ENTREZID")
  return(result$ENTREZID)
}

converted_gene_list <- lapply(gene_list, convert_ensembl_to_entrez)
names(converted_gene_list) <- names(gene_list)

 
 
## KEGG

compKEGG <- compareCluster(geneCluster   = converted_gene_list,
                         fun           = "enrichKEGG",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/ChIPseeker/functional_KEGG_DiffBind05_HET_KO.pdf", width=7, height=7)
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()


## GO
compGO <- compareCluster(geneCluster   = converted_gene_list,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_BP_DiffBind05_HET_KO.pdf", width=7, height=15)
dotplot(compGO, showCategory = 15, title = "GO_Biological Process Enrichment Analysis")
dev.off()

compGO <- compareCluster(geneCluster   = converted_gene_list,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "CC") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_CC_DiffBind05_HET_KO.pdf", width=7, height=10)
dotplot(compGO, showCategory = 15, title = "GO_Cellular Component Enrichment Analysis")
dev.off()

compGO <- compareCluster(geneCluster   = converted_gene_list,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "MF") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_MF_DiffBind05_HET_KO.pdf", width=7, height=10)
dotplot(compGO, showCategory = 15, title = "GO_Molecular Function Enrichment Analysis")
dev.off()


## Disease
compDO <- compareCluster(geneCluster   = converted_gene_list,
                         fun           = "enrichDO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
pdf("output/ChIPseeker/functional_DO_DiffBind05_HET_KO.pdf", width=7, height=5)
dotplot(compDO, showCategory = 15, title = "Disease Ontology Enrichment Analysis")
dev.off()

## ReactomePA
compReactome <- compareCluster(geneCluster   = converted_gene_list,
                               fun           = "enrichPathway",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH",
                               organism      = "human")
pdf("output/ChIPseeker/functional_Reactome_DiffBind05_HET_KO.pdf", width=7, height=10)
dotplot(compReactome, showCategory = 15, title = "Reactome Pathway Enrichment Analysis")
dev.off()


# Output the gene list
## Convert ENSEMBL geneID to SYMBOL gene ID:
convert_ensembl_to_entrez <- function(gene_list) {
  result <- select(org.Hs.eg.db, keys = gene_list, keytype = "ENSEMBL", columns = "SYMBOL")
  return(result$SYMBOL)
}

converted_gene_list <- lapply(gene_list, convert_ensembl_to_entrez)
names(converted_gene_list) <- names(gene_list)

## Function to create a tibble with two columns: gene ID and type
create_gene_tibble <- function(gene_vector, list_name) {
  tibble(gene = gene_vector, type = list_name)
}

## Create a list of tibbles for each gene list
tibble_list <- mapply(create_gene_tibble, converted_gene_list, names(converted_gene_list), SIMPLIFY = FALSE)

## Combine the tibbles into a single tibble
resulting_tibble <- bind_rows(tibble_list)

## SAVE
write.table(resulting_tibble, file = "output/ChIPseeker/gene_list_DiffBind05_HET_KO.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
```

--> Gain/Lost in HET and KO; 590/432 genes and 765/749 (Promoter, gene body and TES);

removing gene dupplicates we got: \
--> Gain/Lost in HET and KO; 462/352 genes and 576/581 (Promoter, gene body and TES);

--> Gain/Lost in H3K27me3 correlates with reduced/increased expression for HET and KO.

--> Few genes are shared like gain in HET and KO; we got specificity of the mutation regarding changing of H3K27me3!

--> NEUROG2 is in HET Lost...

--> NaiaraPlot, HET vs WT; no GO hit for the 77 genes that gain H3K27me3 and decrease in expr; However I found hits when I do not use **`universe`, `background` list of genes!**


**IMPORTANT NOTE: When doing GO, do NOT set a universe (background list of genes) it perform better!**

## DiffBind pvalue05 adjustment


Let's use our pvalue0.05-filtered files where diff. bound sites have been assigned to genes. However, re-do plots of expression and function; using adjusted qvalue 0.05 FDR.
 


```bash
srun --mem=500g --pty bash -l
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

# Import Diff. bind sites pvalue 0.05
HETvsKO_annot <- read_tsv("output/ChIPseeker/annotation_HETvsKO.txt") %>%
    dplyr::select(-c(strand, name, V6, V7, V8))
HETvsWT_annot <- read_tsv("output/ChIPseeker/annotation_HETvsWT.txt") %>%
    dplyr::select(-c(strand, name, V6, V7, V8))
KOvsWT_annot <- read_tsv("output/ChIPseeker/annotation_KOvsWT.txt") %>%
    dplyr::select(-c(strand, name, V6, V7, V8))


# Adjust the qvalue with FDR 0.05 and keep only these one
HETvsKO_annot_FDR05 = HETvsKO_annot %>%
    mutate(FDR05 = p.adjust(pvalue, method = "bonferroni")) %>% # change here fdr or bonferroni
    filter(FDR05 <= 0.05) # change qvalue treshold
HETvsWT_annot_FDR05 = HETvsWT_annot %>%
    mutate(FDR05 = p.adjust(pvalue, method = "bonferroni")) %>%
    filter(FDR05 <= 0.05)
cannot_FDR05 = KOvsWT_annot %>%
    mutate(FDR05 = p.adjust(pvalue, method = "bonferroni")) %>%
    filter(FDR05 <= 0.05)

```

--> Adjusting the pvalue with FDR do NOT change anything, all the pvalue less than 0.05 remains less than 0.05 once FDR adjusted

However; here is the nb of diff detected when playing with the qvalue and test:
- HETvsKO / HETvsWT / KOvsWT
- pvalue 0.05: 2,200 / 1,641 / 2,263
- FDR 0.05: 2,200 / 1,641 / 2,263
- FDR 0.01: 234 / 17 / 221
- BH 0.05: 44 / 11 / 37


Important note about the size of the peak called 400bp per default [here](https://support.bioconductor.org/p/9150459/). For the author, identifying Diff. bound sites in small region is more powerfull for statistic, and likely to identify real different region; as less background (indeed the more you increase >400bp the more you can have background signal...); regions are surrounded around peak summit; then I can enlarge them manually.














# deepTools CutRun vizualization with not the good bigwig

Let's try to use deepTools to explore the data: tutorial [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html) and [here](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html#reference-point).


DeepTools can used bed or bigwig to estimate signal (heatmap or profile) around a point of interest (eg. TSS). Let's use our **median-histone-scaled bigwig**! Generate a matrix for WT,KO,HET and WT,KO,HET,patient for 10 and 50kb around the TSS.


**Profile around the TSS:**
Matrix:

```bash
conda activate deeptools

# example for 1 file 10kb up down TSS:
computeMatrix reference-point --referencePoint TSS \
-b 10000 -a 10000 \
-R [GTF] \
-S [all bigwig] \
--skipZeros \
--blackListFileName [BED] \
-o ~/chipseq/results/visualization/matrixNanog_TSS_chr12.gz \
-p 6 \
--outFileSortedRegions ~/chipseq/results/visualization/regions_TSS_chr12.bed

# Run the different matrix using 200g mem each (last < 48 hrs)
## 10kb
sbatch scripts/matrix_TSS_10kb_all.sh # include the patient # 15041 ok
sbatch scripts/matrix_TSS_10kb.sh # 15042 ok

## 50kb
sbatch scripts/matrix_TSS_50kb_all.sh # 15043 ; FAILED
sbatch scripts/matrix_TSS_50kb.sh # 15044 ; FAILED

## 10kb replicates
sbatch scripts/matrix_WT_Rep_10kb.sh # 15045 ok
sbatch scripts/matrix_KO_Rep_10kb.sh # 15046 ok
sbatch scripts/matrix_HET_Rep_10kb.sh # 15047 ok
```
*NOTE: The succesfull one, show a lot of genes with 'Skipping GENE due to being absent in the compute matrix. Discussion on the issue [here](https://github.com/deeptools/deepTools/issues/447), seems to be caused by dupplicates*



Data vizualization:

```bash
# Check the replicates
## profile; last long so run sbatch script
### example for 1 file;
plotProfile -m output/deeptools/matrix_WT_Rep_10kb.gz \
-out output/deeptools/matrix_WT_Rep_10kb.png \
--perGroup \
--plotTitle "" --samplesLabel "Rep1" "Rep2" "Rep3" "Rep4" \
--refPointLabel "TSS" \
-T "WT read density" \
-z ""
### Run 
sbatch scripts/matrix_WT_Rep_10kb_profile.sh # interactive
sbatch scripts/matrix_HET_Rep_10kb_profile.sh # 15465
sbatch scripts/matrix_KO_Rep_10kb_profile.sh # 15466

## heatmap; last long so run sbatch script
### example for 1 file;
plotHeatmap -m output/deeptools/matrix_WT_Rep_10kb.gz \  
-out output/deeptools/matrix_WT_Rep_10kb_heatmap.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar' \
--zMin -2 --zMax 2  
### Run 
sbatch scripts/matrix_WT_Rep_10kb_heatmap.sh # 15467
sbatch scripts/matrix_HET_Rep_10kb_heatmap.sh # 15468
sbatch scripts/matrix_KO_Rep_10kb_heatmap.sh # 15469

# Check the genotypes
## profile; last long so run sbatch script
sbatch scripts/matrix_TSS_10kb_all_profile.sh # 15525
sbatch scripts/matrix_TSS_10kb_profile.sh # 15527

## profile with clustering
plotProfile -m output/deeptools/matrix_TSS_10kb.gz \
    -out output/deeptools/matrix_TSS_10kb_profile_cluster4.png \
    --kmeans 4 \
    --numPlotsPerRow 4 \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""
### Run 
```

The heatmaps are completely black, because I should have set `--missingDataAsZero` in `computeMatrix`. In addition the clustering rise an error message related to Nan values and then failed (probably related). 

Let's do it again and re-generate profile and heatmap without Nas and zero and see whether that something changed:


```bash
conda activate deeptools

# Replicates
sbatch scripts/matrix_TSS_10kb_WT_missingDataAsZero.sh # 16258 ok
sbatch --dependency=afterany:16258 scripts/matrix_TSS_10kb_WT_missingDataAsZero_profile.sh # 16259 ok
sbatch --dependency=afterany:16258 scripts/matrix_TSS_10kb_WT_missingDataAsZero_heatmap.sh # 16260 ok

sbatch scripts/matrix_TSS_10kb_HET_missingDataAsZero.sh # 17627
sbatch --dependency=afterany:17627 scripts/matrix_TSS_10kb_HET_missingDataAsZero_profile.sh # 17628
sbatch --dependency=afterany:17627 scripts/matrix_TSS_10kb_HET_missingDataAsZero_heatmap.sh # 17630

sbatch scripts/matrix_TSS_10kb_KO_missingDataAsZero.sh # 17631
sbatch --dependency=afterany:17631 scripts/matrix_TSS_10kb_KO_missingDataAsZero_profile.sh # 17632
sbatch --dependency=afterany:17631 scripts/matrix_TSS_10kb_KO_missingDataAsZero_heatmap.sh # 17633

# Genotype TSS (10kb)
sbatch scripts/matrix_TSS_10kb_missingDataAsZero.sh # 15543 ok
sbatch --dependency=afterany:15543 scripts/matrix_TSS_10kb_missingDataAsZero_profile.sh # 15551 ok
sbatch --dependency=afterany:15543 scripts/matrix_TSS_10kb_missingDataAsZero_heatmap.sh # 15552 ok
sbatch --dependency=afterany:15543 scripts/matrix_TSS_10kb_missingDataAsZero_heatmap_cluster4.sh # 15553 FAIL


# Genotype gene body (-5 / +5 kb - TSS / TES)
sbatch scripts/matrix_gene_5kb_missingDataAsZero.sh # 18533
sbatch --dependency=afterany:18533 scripts/matrix_gene_5kb_missingDataAsZero_profile.sh # 18534
sbatch --dependency=afterany:18533 scripts/matrix_gene_5kb_missingDataAsZero_heatmap.sh # 18535
```
*NOTE: still the message Skipping GENE due to being absent in the replicate/genotype compute matrix*

--> Profile are different 

--> Heatmap is working 

--> So let's use these parameters instead (`--missingDataAsZero --skipZeros`)

--> The heatmap color should be fine-tune (maybe do a -2 / 2 limit as in the tutorial)


**With the ratio IgG/H3K27me3 bigwig**
Files in `output/bigwig_histone_NotGenotypeGroup_lib_IggNorm`

```bash
conda activate deeptools
# Replicates
sbatch --dependency=afterany:16250 scripts/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm.sh # 16261 ok
sbatch --dependency=afterany:16261 scripts/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm_profile.sh # 16262 ok 
sbatch --dependency=afterany:16261 scripts/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm_heatmap.sh # 16263 ok

sbatch scripts/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm.sh # 17635
sbatch --dependency=afterany:17635 scripts/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_profile.sh # 17638
sbatch --dependency=afterany:17635 scripts/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_heatmap.sh # 17639

sbatch scripts/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm.sh # 17636
sbatch --dependency=afterany:17636 scripts/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_profile.sh # 17640
sbatch --dependency=afterany:17636 scripts/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_heatmap.sh # 17641


# Genotype TSS (10kb)
sbatch --dependency=afterany:16255 scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm.sh # 16264 ok
sbatch --dependency=afterany:16264 scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm_profile.sh # 16265 ok
sbatch --dependency=afterany:16264 scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm_heatmap.sh # 16266 ok
sbatch --dependency=afterany:16264 scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm_profile_cluster4.sh # 16268 FAIL



# Genotype gene body (-5 / +5 kb - TSS / TES)
sbatch scripts/matrix_gene_5kb_missingDataAsZero_IggNorm.sh # 18738
sbatch --dependency=afterany:18738 scripts/matrix_gene_5kb_missingDataAsZero_IggNorm_profile.sh # 18739
sbatch --dependency=afterany:18738 scripts/matrix_gene_5kb_missingDataAsZero_IggNorm_heatmap.sh # 18740
```
*NOTE: still the message Skipping GENE due to being absent in the replicate/genotype compute matrix*

--> The heatmap color should be fine-tune (maybe do a -2 / 2 limit as in the tutorial)

--> Also maybe better to use GTF with gene only instead of all transcripts here! (`awk '$3 == "gene"' input.gtf > genes.gtf`)

--> Replicate 3 KO is weird





**With the log2ratio IgG/H3K27me3 bigwig**
Files in `output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_log2ratio`


```bash
conda activate deeptools
# Replicates
sbatch --dependency=afterany:19364 scripts/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm_log2ratio.sh # 19392 FAIL; 19826
sbatch --dependency=afterany:19826 scripts/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm_log2ratio_profile.sh # 19393 FAIL; 19827
sbatch --dependency=afterany:19826 scripts/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm_log2ratio_heatmap.sh # 19394 FAIL; fuck heatmap need finetuning

sbatch --dependency=afterany:19364 scripts/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_log2ratio.sh # 19395 FAIL; 19828
sbatch --dependency=afterany:19828 scripts/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_log2ratio_profile.sh # 19399 FAIL; 19829
sbatch --dependency=afterany:19828 scripts/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_log2ratio_heatmap.sh # 19400 FAIL; fuck heatmap need finetuning 

sbatch --dependency=afterany:19364 scripts/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_log2ratio.sh # 19404 FAIL; 19830
sbatch --dependency=afterany:19830 scripts/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_log2ratio_profile.sh # 19406 FAIL; 19831
sbatch --dependency=afterany:19404 scripts/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_log2ratio_heatmap.sh # 19407 FAIL; fuck heatmap need finetuning 


# Genotype TSS (10kb)
sbatch --dependency=afterany:19364 scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm_log2ratio.sh # 19422 ok
sbatch --dependency=afterany:19422 scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm_log2ratio_profile.sh # 19425 ok
sbatch --dependency=afterany:19422 scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm_log2ratio_heatmap.sh # 19427 ok




# Genotype gene body (-5 / +5 kb - TSS / TES)
sbatch --dependency=afterany:19364 scripts/matrix_gene_5kb_missingDataAsZero_IggNorm_log2ratio.sh # 19439 ok
sbatch --dependency=afterany:19439 scripts/matrix_gene_5kb_missingDataAsZero_IggNorm_log2ratio_profile.sh # 19442 ok
sbatch --dependency=afterany:19439 scripts/matrix_gene_5kb_missingDataAsZero_IggNorm_log2ratio_heatmap.sh # 19444 ok
```
*NOTE: Final files for median are called `*_ratio_median.bw`; even the log2ratio and subtract; but they are good, in respectie correct folders*

--> FAIL at matrix because path to bigwig files was not good; corrected and relaunched.





**With the subtract IgG/H3K27me3 bigwig**
Files in `output/bigwig_histone_NotGenotypeGroup_lib_IggNorm_subtract`

```bash
conda activate deeptools
# Replicates
sbatch --dependency=afterany:19382 scripts/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm_subtract.sh # 19460 FAIL; 19832
sbatch --dependency=afterany:19832 scripts/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm_subtract_profile.sh # 19833




sbatch --dependency=afterany:19382 scripts/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_subtract.sh # 19463 FAIL; 19834
sbatch --dependency=afterany:19834 scripts/matrix_TSS_10kb_HET_missingDataAsZero_IggNorm_subtract_profile.sh # 19835


sbatch --dependency=afterany:19382 scripts/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_subtract.sh # 19465 FAIL; 19836
sbatch --dependency=afterany:19836 scripts/matrix_TSS_10kb_KO_missingDataAsZero_IggNorm_subtract_profile.sh # 19837



# Genotype TSS (10kb)
sbatch --dependency=afterany:19382 scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm_subtract.sh # 19468 ok
sbatch scripts/matrix_TSS_10kb_missingDataAsZero_IggNorm_subtract_profile.sh # 19838




# Genotype gene body (-5 / +5 kb - TSS / TES)
sbatch --dependency=afterany:19382 scripts/matrix_gene_5kb_missingDataAsZero_IggNorm_subtract.sh # 19471 ok
sbatch scripts/matrix_gene_5kb_missingDataAsZero_IggNorm_subtract_profile.sh # 19839
```
*NOTE: Final files for median are called `*_ratio_median.bw`; even the log2ratio and subtract; but they are good, in respectie correct folders*


The correct bigwig to take are the `output/bigwig_histone_NotGenotypeGroup` and NOT the `output/bigwig_histone_NotGenotypeGroup_lib`... --> Let's repeat using the correct ones... 


# deepTools CutRun vizualization with the bigwig_histone_NotGenotypeGroup (not bigwig_histone_NotGenotypeGroup_lib)
`output/bigwig_histone_NotGenotypeGroup` and NOT the `output/bigwig_histone_NotGenotypeGroup_lib`

Let's do 5kb around TSS and 1kb around gene body; should look better.

## Non-igg normalized



```bash
conda activate deeptools
# Replicates
sbatch scripts/matrix_TSS_5kb_WT_corr.sh # 19847 ok
sbatch --dependency=afterany:19847 scripts/matrix_TSS_5kb_WT_corr_profile.sh # 19848 ok

sbatch scripts/matrix_TSS_5kb_HET_corr.sh # 19849 ok
sbatch --dependency=afterany:19849 scripts/matrix_TSS_5kb_HET_corr_profile.sh # 19850 ok

sbatch scripts/matrix_TSS_5kb_KO_corr.sh # 19851 ok
sbatch --dependency=afterany:19851 scripts/matrix_TSS_5kb_KO_corr_profile.sh # 19852 ok


# Genotype TSS (10kb) 
sbatch scripts/matrix_TSS_5kb_corr.sh # 19867 ok
sbatch --dependency=afterany:19867 scripts/matrix_TSS_5kb_corr_profile.sh # 19868 ok


# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch scripts/matrix_gene_1kb_corr.sh # 19869 ok
sbatch --dependency=afterany:19869 scripts/matrix_gene_1kb_corr_profile.sh # 19870 ok
```

Looks good, but KO has more reads do not know why...

Let's filter out the very low reads (>2 looks appropriate looking at IGV; `--minThreshold 2` added to computeMatrix); 

```bash
conda activate deeptools
# Replicates (WT only 1st; do HET and KO if needed)
sbatch scripts/matrix_TSS_5kb_WT_corr_min2.sh # 33029 FAIL because treshold of 2 is too high: no value remains
sbatch --dependency=afterany:33029 scripts/matrix_TSS_5kb_WT_corr_min2_profile.sh # 33042 FAIL
sbatch --dependency=afterany:33029 scripts/matrix_TSS_5kb_WT_corr_min2_profileNoPerGroup.sh # 35297 FAIL


# Genotype TSS (10kb) 
sbatch scripts/matrix_TSS_5kb_corr_min2.sh # 33058 ok
sbatch --dependency=afterany:33058 scripts/matrix_TSS_5kb_corr_min2_profile.sh # 33063 ok
sbatch --dependency=afterany:33058 scripts/matrix_TSS_5kb_corr_min2_profileNoPerGroup.sh # 35306 ok


# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch scripts/matrix_gene_1kb_corr_min2.sh # 33067 ok
sbatch --dependency=afterany:33067 scripts/matrix_gene_1kb_corr_min2_profile.sh # 33072 ok
```

--> min2 is faaar too high; this is not value of 2 like in IGV; that is the read counts... Looking at the plots >0.2 should be better.





## igg ratio

```bash
conda activate deeptools
# Replicates
sbatch --dependency=afterany:19856 scripts/matrix_TSS_5kb_WT_corr_IggNorm.sh # 19892 ok
sbatch --dependency=afterany:19892 scripts/matrix_TSS_5kb_WT_corr_IggNorm_profile.sh # 19893 ok

sbatch --dependency=afterany:19856 scripts/matrix_TSS_5kb_HET_corr_IggNorm.sh # 19894 ok
sbatch --dependency=afterany:19894 scripts/matrix_TSS_5kb_HET_corr_IggNorm_profile.sh # 19895 ok

sbatch --dependency=afterany:19856 scripts/matrix_TSS_5kb_KO_corr_IggNorm.sh # 19897 ok
sbatch --dependency=afterany:19897 scripts/matrix_TSS_5kb_KO_corr_IggNorm_profile.sh # 19899 ok


# Genotype TSS (10kb) 
sbatch --dependency=afterany:19856 scripts/matrix_TSS_5kb_corr_IggNorm.sh # 19922 ok
sbatch --dependency=afterany:19922 scripts/matrix_TSS_5kb_corr_IggNorm_profile.sh # 19944 ok


# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch --dependency=afterany:19856 scripts/matrix_gene_1kb_corr_IggNorm.sh # 19946 ok
sbatch --dependency=afterany:19946 scripts/matrix_gene_1kb_corr_IggNorm_profile.sh # 19947 ok
```

Looks good, but KO is kind of weird; Rep2 this time (drop down before TSS)


Let's filter out the very low reads (>2 looks appropriate looking at IGV; `--minThreshold 2` added to computeMatrix); 

```bash
conda activate deeptools
# Replicates (WT only 1st; do HET and KO if needed)
sbatch scripts/matrix_TSS_5kb_WT_corr_IggNorm_min2.sh # 33110 ok
sbatch --dependency=afterany:33110 scripts/matrix_TSS_5kb_WT_corr_IggNorm_min2_profile.sh # 33116 ok
sbatch --dependency=afterany:33110 scripts/matrix_TSS_5kb_WT_corr_IggNorm_min2_profileNoPerGroup.sh # 35139 ok

# Genotype TSS (10kb) 
sbatch scripts/matrix_TSS_5kb_corr_IggNorm_min2.sh # 33130 ok
sbatch --dependency=afterany:33130 scripts/matrix_TSS_5kb_corr_IggNorm_min2_profile.sh # 33131 ok
sbatch --dependency=afterany:33130 scripts/matrix_TSS_5kb_corr_IggNorm_min2_profileNoPerGroup.sh # 35203 ok

# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch scripts/matrix_gene_1kb_corr_IggNorm_min2.sh # 33132 ok
sbatch --dependency=afterany:33132 scripts/matrix_gene_1kb_corr_IggNorm_min2_profile.sh # 33145 ok
```

Removing the `--perGroup` is unknown still as the min2 parameter is failed... 

Lets rerun with min 0.5 instead. No let's use other bigwig


## igg log2ratio

log2ratio is weird looking (liek even the bigwig), fuck that one; would be great if drastic strong differences but that is not the case.

## igg subtract 


```bash
conda activate deeptools
# Replicates
sbatch scripts/matrix_TSS_5kb_WT_corr_IggNorm_subtract.sh # 32940 ok
sbatch --dependency=afterany:32940 scripts/matrix_TSS_5kb_WT_corr_IggNorm_subtract_profile.sh # 32946 ok

sbatch scripts/matrix_TSS_5kb_HET_corr_IggNorm_subtract.sh # 32960 ok
sbatch --dependency=afterany:32960 scripts/matrix_TSS_5kb_HET_corr_IggNorm_subtract_profile.sh # 32961 ok

sbatch scripts/matrix_TSS_5kb_KO_corr_IggNorm_subtract.sh # 32964 ok
sbatch --dependency=afterany:32964 scripts/matrix_TSS_5kb_KO_corr_IggNorm_subtract_profile.sh # 32967 ok


# Genotype TSS (10kb) 
sbatch scripts/matrix_TSS_5kb_corr_IggNorm_subtract.sh # 32971 ok
sbatch --dependency=afterany:32971 scripts/matrix_TSS_5kb_corr_IggNorm_subtract_profile.sh # 32981 ok


# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch scripts/matrix_gene_1kb_corr_IggNorm_subtract.sh # 32996 ok
sbatch --dependency=afterany:32996 scripts/matrix_gene_1kb_corr_IggNorm_subtract_profile.sh # 33000 ok
```



Weird! many negative values! 


Let's filter out the very low reads (>2 looks appropriate looking at IGV; `--minThreshold 2` added to computeMatrix); 

```bash
conda activate deeptools
# Replicates (WT only 1st; do HET and KO if needed)
sbatch scripts/matrix_TSS_5kb_WT_corr_IggNorm_subtract_min2.sh # 33839 ok
sbatch --dependency=afterany:33839 scripts/matrix_TSS_5kb_WT_corr_IggNorm_subtract_min2_profile.sh # 33841 ok
sbatch --dependency=afterany:33839 scripts/matrix_TSS_5kb_WT_corr_IggNorm_subtract_min2_profileNoPerGroup.sh # 35256 ok

# Genotype TSS (10kb) 
sbatch scripts/matrix_TSS_5kb_corr_IggNorm_subtract_min2.sh # 33847 ok
sbatch --dependency=afterany:33847 scripts/matrix_TSS_5kb_corr_IggNorm_subtract_min2_profile.sh # 33850 ok
sbatch --dependency=afterany:33847 scripts/matrix_TSS_5kb_corr_IggNorm_subtract_min2_profileNoPerGroup.sh # 35296 ok


# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch scripts/matrix_gene_1kb_corr_IggNorm_subtract_min2.sh # 33852 ok
sbatch --dependency=afterany:33852 scripts/matrix_gene_1kb_corr_IggNorm_subtract_min2_profile.sh # 33854 ok
```


The filtering is far too high, I even have many negative value...

--> In the end, these bigwigs are not the good one lol!!!

# deepTools CutRun vizualization with the bigwig_histone


## raw bigwig_histone

```bash
conda activate deeptools
# Replicates
sbatch scripts/matrix_TSS_5kb_WT_histone.sh # 47836 ok
sbatch --dependency=afterany:47836 scripts/matrix_TSS_5kb_WT_histone_profile.sh # 47837 ok
sbatch --dependency=afterany:47836 scripts/matrix_TSS_5kb_WT_histone_profile_noPerGroup.sh # 47838 ok
sbatch scripts/matrix_TSS_5kb_WT_histone_ZeroKept.sh # 47840 ok
sbatch --dependency=afterany:47840 scripts/matrix_TSS_5kb_WT_histone_profile_ZeroKept.sh # 47841 ok

sbatch scripts/matrix_TSS_5kb_HET_histone.sh # 47842 ok
sbatch --dependency=afterany:47842 scripts/matrix_TSS_5kb_HET_histone_profile.sh # 47843 ok

sbatch scripts/matrix_TSS_5kb_KO_histone.sh # 47844 ok
sbatch --dependency=afterany:47844 scripts/matrix_TSS_5kb_KO_histone_profile.sh # 47845 ok

# Genotype TSS (5kb) 
sbatch scripts/matrix_TSS_5kb_histone.sh # 47846 ok
sbatch --dependency=afterany:47846 scripts/matrix_TSS_5kb_histone_profile.sh # 47847 ok

# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch scripts/matrix_gene_1kb_histone.sh # 47848 ok
sbatch --dependency=afterany:47848 scripts/matrix_gene_1kb_histone_profile.sh # 47849 ok
```

--> Keepin the zero with removing `--skipZeros` do not change a shit.

--> Removing the `--perGroup` in profile is making facet_wrap per sample!

Still unclear why replicate are different...



## log igg

```bash
conda activate deeptools
# Replicates
sbatch scripts/matrix_TSS_5kb_WT_histone_ratio.sh # 47851 ok
sbatch --dependency=afterany:47851 scripts/matrix_TSS_5kb_WT_histone_ratio_profile.sh # 47852 ok
sbatch --dependency=afterany:47851 scripts/matrix_TSS_5kb_WT_histone_ratio_profile_noPerGroup.sh # 47853 ok
sbatch scripts/matrix_TSS_5kb_WT_histone_ratio_ZeroKept.sh # 47854 ok
sbatch --dependency=afterany:47854 scripts/matrix_TSS_5kb_WT_histone_ratio_ZeroKept_profile.sh # 47856 ok

sbatch scripts/matrix_TSS_5kb_HET_histone_ratio.sh # 47857 ok
sbatch --dependency=afterany:47857 scripts/matrix_TSS_5kb_HET_histone_ratio_profile.sh # 47858 ok

sbatch scripts/matrix_TSS_5kb_KO_histone_ratio.sh # 47859 ok
sbatch --dependency=afterany:47859 scripts/matrix_TSS_5kb_KO_histone_ratio_profile.sh # 47860 ok

# Genotype TSS (5kb) 
sbatch scripts/matrix_TSS_5kb_histone_ratio.sh # 47861 ok
sbatch --dependency=afterany:47861 scripts/matrix_TSS_5kb_histone_ratio_profile.sh # 47877 ok

# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch scripts/matrix_gene_1kb_histone_ratio.sh # 47879 ok
sbatch --dependency=afterany:47879 scripts/matrix_gene_1kb_histone_ratio_profile.sh # 47883 ok
```




## subtract igg

```bash
conda activate deeptools
# Replicates
sbatch scripts/matrix_TSS_5kb_WT_histone_subtract.sh # 47889 ok
sbatch --dependency=afterany:47889 scripts/matrix_TSS_5kb_WT_histone_subtract_profile.sh # 47890 ok
sbatch --dependency=afterany:47889 scripts/matrix_TSS_5kb_WT_histone_subtract_profile_noPerGroup.sh # 47892 ok
sbatch scripts/matrix_TSS_5kb_WT_histone_subtract_ZeroKept.sh # 47893 ok
sbatch --dependency=afterany:47893 scripts/matrix_TSS_5kb_WT_histone_subtract_ZeroKept_profile.sh # 47894 ok

sbatch scripts/matrix_TSS_5kb_HET_histone_subtract.sh # 47895 ok
sbatch --dependency=afterany:47895 scripts/matrix_TSS_5kb_HET_histone_subtract_profile.sh # 47896 ok

sbatch scripts/matrix_TSS_5kb_KO_histone_subtract.sh # 47897 ok
sbatch --dependency=afterany:47897 scripts/matrix_TSS_5kb_KO_histone_subtract_profile.sh # 47898 ok

# Genotype TSS (5kb) 
sbatch scripts/matrix_TSS_5kb_histone_subtract.sh # 47899 ok
sbatch --dependency=afterany:47899 scripts/matrix_TSS_5kb_histone_subtract_profile.sh # 47900 ok

# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch scripts/matrix_gene_1kb_histone_subtract.sh # 47901 ok
sbatch --dependency=afterany:47901 scripts/matrix_gene_1kb_histone_subtract_profile.sh # 47902 ok
```

--> KO and HET go down at TSS lol WTF!



# deepTools CutRun vizualization with the bigwig_DiffBind_TMM


## bigwig_DiffBind_TMM

```bash
conda activate deeptools
# All genes
## Replicates
sbatch --dependency=afterany:48710 scripts/matrix_TSS_5kb_WT_DiffBind_TMM.sh # 48724 ok
sbatch --dependency=afterany:48724 scripts/matrix_TSS_5kb_WT_DiffBind_TMM_profile.sh # 48725 ok

sbatch --dependency=afterany:48710 scripts/matrix_TSS_5kb_HET_DiffBind_TMM.sh # 48726 ok
sbatch --dependency=afterany:48726 scripts/matrix_TSS_5kb_HET_DiffBind_TMM_profile.sh # 48728 ok

sbatch --dependency=afterany:48710 scripts/matrix_TSS_5kb_KO_DiffBind_TMM.sh # 48729 ok
sbatch --dependency=afterany:48729 scripts/matrix_TSS_5kb_KO_DiffBind_TMM_profile.sh # 48730 ok

## Genotype TSS (5kb) 
sbatch --dependency=afterany:48710 scripts/matrix_TSS_5kb_DiffBind_TMM.sh # 48731 ok
sbatch --dependency=afterany:48731 scripts/matrix_TSS_5kb_DiffBind_TMM_profile.sh # 48732 ok

## Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch --dependency=afterany:48710 scripts/matrix_gene_1kb_DiffBind_TMM.sh # 48733 ok
sbatch --dependency=afterany:48733 scripts/matrix_gene_1kb_DiffBind_TMM_profile.sh # 48734 ok

## clustering
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_profile_kmeans.sh # 49710; 49731 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans.sh # 49734 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_profile_hclust4.sh # 49707  TOO LONG; fuck it ok

### Compute how many genes per cluster

grep -c "cluster_1" output/deeptools/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans6.txt # n=1,755
grep -c "cluster_2" output/deeptools/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans6.txt # n=2,999
grep -c "cluster_3" output/deeptools/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans6.txt # n=9,201
grep -c "cluster_4" output/deeptools/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans6.txt # n=32,626
grep -c "cluster_5" output/deeptools/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans6.txt # n=50,003
grep -c "cluster_6" output/deeptools/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans6.txt # n=104,009
```
*NOTE: for the clustering to work the `--plotTitle` and other `title stuff` arguments need to be removed*

--> The replicates are much more better!!! Very comparable!

--> The genotype comparison do not show striking difference (almost nothing)

--> hclust method take forever, fuck it, let;s use kmeans


Now let's display only the DEGs, or diff. bound genes or both (show 6 clusters):

Need generate gtf files that contain (keep gtf in `output/deseq2_hg38`):
**only the DEGs**:
```R
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# Collect all DEGs (both in WT vs KO and WT vs HET and HET vs KO) as a gene list
## Import deseq2 output and filter qvalue 0.05
HETvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_8wN_HET_vs_8wN_WT.txt')) %>%
    filter(padj <= 0.05) %>%
    dplyr::select(-"...1") %>%
    add_column(contrast = "HETvsWT")
KOvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_8wN_KO_vs_8wN_WT.txt')) %>%
    filter(padj <= 0.05) %>%
    dplyr::select(-"...1") %>%
    add_column(contrast = "KOvsWT")
KOvsHET = as_tibble(read_csv('output/deseq2_hg38/raw_8wN_KO_vs_8wN_HET.txt')) %>%
    filter(padj <= 0.05) %>%
    dplyr::select(-"...1") %>%
    add_column(contrast = "KOvsHET")

DEGs = HETvsWT %>%
    bind_rows(KOvsWT, KOvsHET)

## Import the GTF file
gtf <- import('../../Master/meta/ENCFF159KBI.gtf')

## Filter in the GTF file with only the DEGs
gtf_DEGs <- gtf %>% 
  as_tibble() %>% 
  filter(gene_id %in% DEGs$gene)

## Save the new GTF
export(gtf_DEGs, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN.gtf")
```

**Only the diffbound sites DiffBind05**:
```R
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
# Collect all diff. bound sites (HET and KO changes vs the WT) as a gene list
HET_KO = read_tsv('../003__CutRun/output/ChIPseeker/gene_list_DiffBind05_HET_KO.txt')

## Import the GTF file
gtf <- import('../../Master/meta/ENCFF159KBI.gtf')

## Filter in the GTF file with only the DEGs
gtf_HET_KO <- gtf %>% 
  as_tibble() %>% 
  filter(gene_name %in% HET_KO$gene)

## Save the new GTF
export(gtf_HET_KO, con = "output/deseq2_hg38/ENCFF159KBI_DiffBind05_8wN.gtf")
```

**only the DEGs that are also DiffBound05**:
Let's fusion our two gtfs and keep the genes that overlap:

```R
# import gtf
gtf_HET_KO = import("output/deseq2_hg38/ENCFF159KBI_DiffBind05_8wN.gtf")
gtf_DEGs = import("output/deseq2_hg38/ENCFF159KBI_DEGs_8wN.gtf")

# Find overlap genes
overlapping_genes <- findOverlaps(gtf_DEGs, gtf_HET_KO)

# Extract overlapping genes and combine them
gtf_DEGs_overlapping <- gtf_DEGs[queryHits(overlapping_genes)]
gtf_HET_KO_overlapping <- gtf_HET_KO[subjectHits(overlapping_genes)]
combined_gtf <- c(gtf_DEGs_overlapping, gtf_HET_KO_overlapping)

# Export gtf
export(combined_gtf, con = "output/deseq2_hg38/ENCFF159KBI_DiffBind05_DEGs_8wN.gtf")
```
*NOTE: dupplicated rows in the `output/deseq2_hg38/ENCFF159KBI_DiffBind05_DEGs_8wN.gtf`; deepTools seems to bug with that*

Let's use instead bedtools intersect:

```bash
conda activate BedToBigwig

bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DiffBind05_8wN.gtf -b ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN.gtf > ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DiffBind05_DEGs_8wN_bedtools.gtf
```
*NOTE: -u to display unique hit between a and b*


--> Files looks all good (no more dupplicated rows and only includes DEGs/DiffBind genes)



Generate **clustering matrix with the DEGs or/and DiffBound05**:

```bash
conda activate deeptools
# DEGs
## clustering
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_DEGs.sh # 157192 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_DEGs_profile.sh # 158363 ok
sbatch --dependency=afterany:157192 scripts/matrix_gene_1kb_DiffBind_TMM_DEGs_heatmap_kmeans.sh # 157287 ok


# Diff bound genes
## clustering
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_DiffBind05.sh # 157193 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_DiffBind05_profile.sh # 158064 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_DiffBind05_heatmap_kmeans.sh # 157296  ok

# DEGs and diff bound genes
## clustering
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_DiffBind05_DEGs.sh # 158062 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_DiffBind05_DEGs_profile.sh # 158065 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_DiffBind05_DEGs_heatmap_kmeans.sh #  ok
```
- *NOTE: dupplicated rows in the `output/deseq2_hg38/ENCFF159KBI_DiffBind05_DEGs_8wN.gtf` but deepTools do NOT take them into account* 
- *NOTE: To display different nb of clusters; kmeans clustering number has been changed manually using `nano`*
- *NOTE: If color shit can play with --colorNumber and --colorList (like to chose which color and when changing it)*

**Understand deepTools profiling**:

Plot NEUROG2 and other genes to understand how deepTools work:



```bash
conda activate deeptools
# Generate the NEUROG2 gtf file
grep "NEUROG2" meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_NEUROG2.gtf

# Generate matrix/plot
## NEUROG2 TSS
computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_NEUROG2.gtf \
    -S output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_DiffBind_TMM_NEUROG2.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_DiffBind_TMM_NEUROG2.bed
plotProfile -m output/deeptools/matrix_TSS_5kb_DiffBind_TMM_NEUROG2.gz \
    -out output/deeptools/matrix_TSS_5kb_DiffBind_TMM_NEUROG2_profile.png \
    --perGroup \
    --colors black blue red \
    --plotTitle "" --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


## NEUROG2 gene body
computeMatrix scale-regions \
    -b 1000 -a 1000 \
    -R meta/ENCFF159KBI_NEUROG2.gtf \
    -S output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb_DiffBind_TMM_NEUROG2.gz \
    -p max/2 \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_DiffBind_TMM_NEUROG2.bed
plotProfile -m output/deeptools/matrix_gene_1kb_DiffBind_TMM_NEUROG2.gz \
    -out output/deeptools/matrix_gene_1kb_DiffBind_TMM_NEUROG2_profile.png \
    --perGroup \
    --colors black blue red \
    --plotTitle "" --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""
```
--> deepTools uses the exact value from the bigwig; it really plot the bigwig, literaly!

--> So I realize NEUROG2 as value of around 40, but plot I generated for all genes as value of 2 so lot of transcripts that do not contain peaks are plotted!

Let's **filter out transcripts that do not contain any peaks**; whatever their genotypes.


```bash
# For each genotype collect geneSymbol (gene name) list
## Print 20th column in each rows; sort; remov dupplicates
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_WT.txt | sort | uniq > output/ChIPseeker/annotation_WT_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_HET.txt | sort | uniq > output/ChIPseeker/annotation_HET_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_KO.txt | sort | uniq > output/ChIPseeker/annotation_KO_geneSymbol.txt

## Concatenate all gene into 1 file
cat output/ChIPseeker/annotation_WT_geneSymbol.txt output/ChIPseeker/annotation_HET_geneSymbol.txt output/ChIPseeker/annotation_KO_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_geneSymbol.txt
```
*NOTE: awk process file line per line; `-F'\t'` = tab-separated; `(NR==1 || FNR>1)` tell that header line and all lines after are processed; `{print $20}` will only print the 20th column=gene name `sort` sort output because `uniq` only remove adjacent dupplicated rows*


Now generate a **new gtf file that contains these genes**:

```bash
# Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_geneSymbol.txt > output/ChIPseeker/annotation_as_gtf_geneSymbol.txt
# Filter the gtf
grep -Ff output/ChIPseeker/annotation_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_peak.gtf
```
- *NOTE: `sed` here add `gene name "` before and `"` after*
- *NOTE: `-Ff` tel grep to read the pattern in each row of the file*


Generate **clustering matrix with the peak-containing genes**:

```bash
conda activate deeptools
# all genotypes
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_peaks.sh # 159551 ok
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_peaks.sh # 159550 ok
```
*NOTE: each command contain plotProfile and plotHeatmap with 6 clusters*

Mean value is increased, but still median profile is around 5-6 which is very low/small peaks. Let's filter for peak > 5 value (looks real)


```bash
conda activate deeptools
# all genotypes
## min 5
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_peaks_min5.sh # 122660 ok
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_peaks_min5.sh # 122661 fail; not enough value

## min 1
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_peaks_min1.sh # 164081 ok

## keep value of 0
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_peaks_keepZero.sh # 164102 ok

```
--> Seems that the filtering is now too strong; indeed we remove all value smaller than 5 and even 1; but when there is a peak; the side around TSS/TES is smaller than 5 so also removed...

--> Using not filtering-out 0 (keepZero) is almost same as removing them

Generate **clustering matrix with the peak (>5 score)-containing genes**:

For each genotype; use scaled bigwig (Bigwig_DiffBind_TMM) with pool MACS2 peaks and do the following:
- Convert scaled bigwig to bedGraph --> Not needed; already got bedGraph! `output/bigwig_DiffBind_TMM/`
- With `bedtools map` compute the maximum score found in each peak regions (use significant peaks: `output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_pool_peaks.broadPeak`)
- Filter the bed to keep only the peak that have a max value of 5

*NOTE: install `bigWigToBedGraph` with `conda install -c bioconda ucsc-bigwigtobedgraph` in BedToBigwig conda env.*

```bash
conda activate BedToBigwig

# Collect bedGraph
output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted.bedGraph 
output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted.bedGraph 
output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted.bedGraph
# Collect bed peaks
output/macs2/broad_blacklist_qval2.30103/8wN_HET_H3K27me3_pool_peaks.broadPeak
output/macs2/broad_blacklist_qval2.30103/8wN_KO_H3K27me3_pool_peaks.broadPeak
output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_pool_peaks.broadPeak

# Compute maximum score with bedtools map (last columns = max value)
bedtools map -a output/macs2/broad_blacklist_qval2.30103/8wN_HET_H3K27me3_pool_peaks.broadPeak -b output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted.bedGraph -c 4 -o max > output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted_maxScorePoolpeaks.bed
bedtools map -a output/macs2/broad_blacklist_qval2.30103/8wN_KO_H3K27me3_pool_peaks.broadPeak -b output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted.bedGraph -c 4 -o max > output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks.bed
bedtools map -a output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_pool_peaks.broadPeak -b output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted.bedGraph -c 4 -o max > output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks.bed

# Filter-out peaks with max score below 1/2/5/10/20
## HET
awk '$10 >= 1' output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted_maxScorePoolpeaks_min1.bed
## KO
awk '$10 >= 1' output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min1.bed
awk '$10 >= 2' output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min2.bed
awk '$10 >= 5' output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min5.bed
awk '$10 >= 10' output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min10.bed
awk '$10 >= 20' output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min20.bed
## WT
awk '$10 >= 1' output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min1.bed
awk '$10 >= 2' output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min2.bed
awk '$10 >= 5' output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min5.bed
awk '$10 >= 10' output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min10.bed
awk '$10 >= 20' output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks.bed > output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min20.bed

# Filter-in the peak with score >1/2/5/10/20 in the ChIPseeker annotation files
## HET
bedtools intersect -wa -a output/ChIPseeker/annotation_HET.bed -b output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted_maxScorePoolpeaks_min1.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min1.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_HET.bed -b output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted_maxScorePoolpeaks_min2.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min2.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_HET.bed -b output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted_maxScorePoolpeaks_min5.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min5.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_HET.bed -b output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted_maxScorePoolpeaks_min10.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min10.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_HET.bed -b output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.sorted_maxScorePoolpeaks_min20.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min20.bed
## KO
bedtools intersect -wa -a output/ChIPseeker/annotation_KO.bed -b output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min1.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min1.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_KO.bed -b output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min2.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min2.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_KO.bed -b output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min5.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min5.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_KO.bed -b output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min10.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min10.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_KO.bed -b output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.sorted_maxScorePoolpeaks_min20.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min20.bed
## WT
bedtools intersect -wa -a output/ChIPseeker/annotation_WT.bed -b output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min1.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_WT.bed -b output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min2.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min2.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_WT.bed -b output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min5.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min5.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_WT.bed -b output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min10.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min10.bed
bedtools intersect -wa -a output/ChIPseeker/annotation_WT.bed -b output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.sorted_maxScorePoolpeaks_min20.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min20.bed

# Filter out all rows that contain "Intergenic"
## HET
grep -v "Intergenic" output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min1.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min1_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min2.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min2_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min5.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min5_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min10.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min10_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min20.bed > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min20_noIntergenic.bed
## KO
grep -v "Intergenic" output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min1.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min1_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min2.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min2_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min5.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min5_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min10.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min10_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min20.bed > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min20_noIntergenic.bed
## WT
grep -v "Intergenic" output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min2.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min2_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min5.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min5_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min10.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min10_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min20.bed > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min20_noIntergenic.bed

# For each genotype collect geneSymbol (gene name) list
## Print 20th column in each rows; sort; remov dupplicates
### HET
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min1_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min2_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min2_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min5_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min5_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min10_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min10_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min20_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min20_noIntergenic_geneSymbol.txt
### KO
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min1_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min2_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min2_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min5_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min5_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min10_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min10_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min20_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min20_noIntergenic_geneSymbol.txt



### WT
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min2_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min2_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min5_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min5_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min10_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min10_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min20_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min20_noIntergenic_geneSymbol.txt

## Concatenate all gene into 1 file
### min1
cat output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt
### min2
cat output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min2_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min2_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min2_noIntergenic_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_maxScorePoolpeaks_min2_noIntergenic_geneSymbol.txt

### min5
cat output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min5_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min5_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min5_noIntergenic_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_maxScorePoolpeaks_min5_noIntergenic_geneSymbol.txt

### min10
cat output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min10_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min10_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min10_noIntergenic_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_maxScorePoolpeaks_min10_noIntergenic_geneSymbol.txt

### min20
cat output/ChIPseeker/annotation_HET_maxScorePoolpeaks_min20_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_KO_maxScorePoolpeaks_min20_noIntergenic_geneSymbol.txt output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min20_noIntergenic_geneSymbol.txt | sort | uniq > output/ChIPseeker/annotation_maxScorePoolpeaks_min20_noIntergenic_geneSymbol.txt

# Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_maxScorePoolpeaks_min1_noIntergenic_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_maxScorePoolpeaks_min2_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_maxScorePoolpeaks_min2_noIntergenic_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_maxScorePoolpeaks_min5_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_maxScorePoolpeaks_min5_noIntergenic_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_maxScorePoolpeaks_min10_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_maxScorePoolpeaks_min10_noIntergenic_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_maxScorePoolpeaks_min20_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_maxScorePoolpeaks_min20_noIntergenic_as_gtf_geneSymbol.txt

# Filter the gtf
grep -Ff output/ChIPseeker/annotation_maxScorePoolpeaks_min1_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_maxScorePoolpeaks_min1_noIntergenic.gtf
grep -Ff output/ChIPseeker/annotation_maxScorePoolpeaks_min2_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_maxScorePoolpeaks_min2_noIntergenic.gtf
grep -Ff output/ChIPseeker/annotation_maxScorePoolpeaks_min5_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_maxScorePoolpeaks_min5_noIntergenic.gtf
grep -Ff output/ChIPseeker/annotation_maxScorePoolpeaks_min10_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_maxScorePoolpeaks_min10_noIntergenic.gtf
grep -Ff output/ChIPseeker/annotation_maxScorePoolpeaks_min20_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_maxScorePoolpeaks_min20_noIntergenic.gtf
```
- *NOTE: `bedtools map -c` = the column from the bedGraph to use for the operation; -o = max; collect maximum score*
- *NOTE: `grep -v "something" input > output` --> Remove all rows containing "something"*
- *NOTE: I generated the `output/ChIPseeker/annotation_*.bed` by changing `.txt` to `.bed` and removing the 1st row*
- *NOTE: using `wc -l [FILE]` ; I checked the nb of rows in each file, the more I filter, the less the transcripts count; so looks good*

Let's now generate the matrix and deepTools plots:

```bash
conda activate deeptools
# all genotypes
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic.sh # 222843 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min2_noIntergenic.sh # 222847 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min5_noIntergenic.sh # 222851 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min10_noIntergenic.sh # 222854 ok
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min20_noIntergenic.sh # 222855 ok

sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic.sh # 222861 ok
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_maxScorePoolpeaks_min2_noIntergenic.sh # 222865 ok
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_maxScorePoolpeaks_min5_noIntergenic.sh # 222869 ok
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_maxScorePoolpeaks_min10_noIntergenic.sh # 222870 ok
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_maxScorePoolpeaks_min20_noIntergenic.sh # 222874 ok
```
*NOTE: each command contain plotProfile and plotHeatmap with 6 clusters*

--> min20 peak TSS is around 9-10 (got KO decrease upstream TSS; HET increase downstrem TSS)
--> min2/5/10 peak TSS is around 7 (got KO decrease upstream TSS; HET increase downstrem TSS; more marked than min20)
--> min1 is around 6-7 (best to choose as keep as much genes)

nb of unique genes in each category: *20= 3,894; 10= 7,624; 1= 8,625* `awk '$3 ~ /gene/' file.txt | wc -l`

Matrix and deeptools for **DEGs and diffbind and non-intergenic at least 1 genotype peaks**:


```bash
# Generate gtf DEGs diffbind peak no intergenic
conda activate BedToBigwig

bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DiffBind05_DEGs_8wN_bedtools.gtf -b meta/ENCFF159KBI_maxScorePoolpeaks_min1_noIntergenic.gtf > meta/ENCFF159KBI_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs_8wN_bedtools.gtf

# deepTools plot
conda activate deeptools

sbatch scripts/matrix_gene_1kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs.sh # 255664
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_maxScorePoolpeaks_min1_noIntergenic_DiffBind05_DEGs.sh # 255776
```
*NOTE: -u to display unique hit between a and b*


Matrix and deeptools for **Non-intergenic peak in WT (whether or not found in HET/KO)**:


```bash
conda activate deeptools
# File with genes from WT with min score of 1; no intergenic
output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt 

# Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1_noIntergenic_as_gtf_geneSymbol.txt

# Filter the gtf
grep -Ff output/ChIPseeker/annotation_WT_maxScorePoolpeaks_min1_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_WT_maxScorePoolpeaks_min1_noIntergenic.gtf

# deepTools plot
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_WT_maxScorePoolpeaks_min1_noIntergenic.sh # 589600 
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_WT_maxScorePoolpeaks_min1_noIntergenic.sh # 589874
```

--> KO less spread around TSS; HET more signal within gene body



Matrix and deeptools for genes where **H3K27me3 goes up in HET and H3K27me3 goes down in KO**

```bash
conda activate BedToBigwig
# File where the diffbound sites has been assigned to genes
output/ChIPseeker/annotation_HETvsWT.txt
output/ChIPseeker/annotation_KOvsWT.txt

# Filter the peaks
## Filter HET Diff bound sites that goes Up (contrast2 positif=Up in HET)
awk '$10 > 0' output/ChIPseeker/annotation_HETvsWT.txt > output/ChIPseeker/annotation_HETvsWT_UpHET.txt
## Filter HET Diff bound sites that goes Down (contrast2 negatif= Down in HET)
awk '$10 < 0' output/ChIPseeker/annotation_HETvsWT.txt > output/ChIPseeker/annotation_HETvsWT_DownHET.txt

## Filter KO Diff bound sites that goes Up (contrast2 positif=Up in KO)
awk '$10 > 0' output/ChIPseeker/annotation_KOvsWT.txt > output/ChIPseeker/annotation_KOvsWT_UpKO.txt
## Filter KO Diff bound sites that goes Down (contrast2 negatif= Down in KO)
awk '$10 < 0' output/ChIPseeker/annotation_KOvsWT.txt > output/ChIPseeker/annotation_KOvsWT_DownKO.txt

# filter out the intergenic and convert to bed (remove header manually)
grep -v "Intergenic" output/ChIPseeker/annotation_HETvsWT_UpHET.txt > output/ChIPseeker/annotation_HETvsWT_UpHET_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_HETvsWT_DownHET.txt > output/ChIPseeker/annotation_HETvsWT_DownHET_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_KOvsWT_UpKO.txt > output/ChIPseeker/annotation_KOvsWT_UpKO_noIntergenic.bed
grep -v "Intergenic" output/ChIPseeker/annotation_KOvsWT_DownKO.txt > output/ChIPseeker/annotation_KOvsWT_DownKO_noIntergenic.bed


# Collect gene ID
awk -F'\t' '(NR==1 || FNR>1) {print $22}' output/ChIPseeker/annotation_HETvsWT_UpHET_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_HETvsWT_UpHET_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $22}' output/ChIPseeker/annotation_HETvsWT_DownHET_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_HETvsWT_DownHET_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $22}' output/ChIPseeker/annotation_KOvsWT_UpKO_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_KOvsWT_UpKO_noIntergenic_geneSymbol.txt
awk -F'\t' '(NR==1 || FNR>1) {print $22}' output/ChIPseeker/annotation_KOvsWT_DownKO_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_KOvsWT_DownKO_noIntergenic_geneSymbol.txt

# Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_HETvsWT_UpHET_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_HETvsWT_UpHET_noIntergenic_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_HETvsWT_DownHET_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_HETvsWT_DownHET_noIntergenic_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_KOvsWT_UpKO_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_KOvsWT_UpKO_noIntergenic_as_gtf_geneSymbol.txt
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_KOvsWT_DownKO_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_KOvsWT_DownKO_noIntergenic_as_gtf_geneSymbol.txt

# Filter the gtf
grep -Ff output/ChIPseeker/annotation_HETvsWT_UpHET_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_HETvsWT_UpHET_noIntergenic.gtf
grep -Ff output/ChIPseeker/annotation_HETvsWT_DownHET_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_HETvsWT_DownHET_noIntergenic.gtf
grep -Ff output/ChIPseeker/annotation_KOvsWT_UpKO_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_KOvsWT_UpKO_noIntergenic.gtf
grep -Ff output/ChIPseeker/annotation_KOvsWT_DownKO_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_KOvsWT_DownKO_noIntergenic.gtf

# Combine GTF
## Up HET and Down in KO
bedtools intersect -wa -u -a meta/ENCFF159KBI_HETvsWT_UpHET_noIntergenic.gtf -b meta/ENCFF159KBI_KOvsWT_DownKO_noIntergenic.gtf > meta/ENCFF159KBI_UpHET_DownKO_noIntergenic.gtf
## Down in HET and up in KO
bedtools intersect -wa -u -a meta/ENCFF159KBI_HETvsWT_DownHET_noIntergenic.gtf -b meta/ENCFF159KBI_KOvsWT_UpKO_noIntergenic.gtf > meta/ENCFF159KBI_DownHET_UpKO_noIntergenic.gtf
## Down HET and Down KO
bedtools intersect -wa -u -a meta/ENCFF159KBI_HETvsWT_DownHET_noIntergenic.gtf -b meta/ENCFF159KBI_KOvsWT_DownKO_noIntergenic.gtf > meta/ENCFF159KBI_DownHET_DownKO_noIntergenic.gtf
## Up HET and Up KO
bedtools intersect -wa -u -a meta/ENCFF159KBI_HETvsWT_UpHET_noIntergenic.gtf -b meta/ENCFF159KBI_KOvsWT_UpKO_noIntergenic.gtf > meta/ENCFF159KBI_UpHET_UpKO_noIntergenic.gtf


# deepTools plot
conda activate deeptools
## Up HET and Down in KO
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_UpHET_DownKO_noIntergenic.sh # 591041
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_UpHET_DownKO_noIntergenic.sh # 591040
```

--> KO less spread around TSS; HET more signal within gene body; but less genes/transcripts so plot more variable...


Matrix and deeptools for genes where **H3K27me3 goes up in HET AND expression goes down; and H3K27me3 goes down in KO AND expression goes Up**

Filter the gene Up/Down in each mutant in R and generate GTF:
```bash
conda activate deseq2
cd /scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq
```


```R
# library
library("rtracklayer")
library("tidyverse")

## Import deseq2 output and filter qvalue 0.05
HETvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_8wN_HET_vs_8wN_WT.txt')) %>%
    filter(padj <= 0.001) %>%  # HERE CHANGE qvalue if NEEDED
    dplyr::select(-"...1") %>%
    add_column(contrast = "HETvsWT")
KOvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_8wN_KO_vs_8wN_WT.txt')) %>%
    filter(padj <= 0.001) %>% # HERE CHANGE qvalue if NEEDED
    dplyr::select(-"...1") %>%
    add_column(contrast = "KOvsWT")

## Filter Up/Down
HET_Down = HETvsWT %>% 
    filter(log2FoldChange < 0) 
HET_Up = HETvsWT %>% 
    filter(log2FoldChange > 0) 
KO_Down = KOvsWT %>% 
    filter(log2FoldChange < 0) 
KO_Up = KOvsWT %>% 
    filter(log2FoldChange > 0) 


## Import the GTF file
gtf <- import('../../Master/meta/ENCFF159KBI.gtf')

## Filter in the GTF file with only the DEGs
gtf_DEGs_HET_Down <- gtf %>% 
  as_tibble() %>% 
  filter(gene_id %in% HET_Down$gene)
gtf_DEGs_HET_Up <- gtf %>% 
  as_tibble() %>% 
  filter(gene_id %in% HET_Up$gene)
gtf_DEGs_KO_Down <- gtf %>% 
  as_tibble() %>% 
  filter(gene_id %in% KO_Down$gene)
gtf_DEGs_KO_Up <- gtf %>% 
  as_tibble() %>% 
  filter(gene_id %in% KO_Up$gene)

## Save the new GTF
export(gtf_DEGs_HET_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Down.gtf")
export(gtf_DEGs_HET_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Up.gtf")
export(gtf_DEGs_KO_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Down.gtf")
export(gtf_DEGs_KO_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Up.gtf")

# For other qvalue:
export(gtf_DEGs_HET_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Down_qval001.gtf")
export(gtf_DEGs_HET_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Up_qval001.gtf")
export(gtf_DEGs_KO_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Down_qval001.gtf")
export(gtf_DEGs_KO_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Up_qval001.gtf")
```

Then, `bedools intersect` the **DEGs GTF (HET down and KO up)** and generate the deepTools plots:

```bash
conda activate BedToBigwig

## expression Down HET and Up in KO DEGs only
cat ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Down.gtf ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Up.gtf > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_unsort.gtf
sort meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_unsort.gtf | uniq > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_unsort_unique.gtf
bedtools sort -i meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_unsort_unique.gtf > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_sort.gtf
### qval 0.01
cat ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Down_qval01.gtf ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Up_qval01.gtf > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval01_unsort.gtf
sort meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval01_unsort.gtf | uniq > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval01_unsort_unique.gtf
bedtools sort -i meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval01_unsort_unique.gtf > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval01_sort.gtf
### qval 0.001
cat ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Down_qval001.gtf ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Up_qval001.gtf > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval001_unsort.gtf
sort meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval001_unsort.gtf | uniq > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval001_unsort_unique.gtf
bedtools sort -i meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval001_unsort_unique.gtf > meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval001_sort.gtf


## expression Up HET and Down in KO DEGs only
cat ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Up.gtf ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Down.gtf > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_unsort.gtf
sort meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_unsort.gtf | uniq > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_unsort_unique.gtf
bedtools sort -i meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_unsort_unique.gtf > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_sort.gtf
### qval 0.01
cat ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Up_qval01.gtf ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Down_qval01.gtf > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval01_unsort.gtf
sort meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval01_unsort.gtf | uniq > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval01_unsort_unique.gtf
bedtools sort -i meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval01_unsort_unique.gtf > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval01_sort.gtf
### qval 0.001
cat ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_HET_Up_qval001.gtf ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_8wN_KO_Down_qval001.gtf > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval001_unsort.gtf
sort meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval001_unsort.gtf | uniq > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval001_unsort_unique.gtf
bedtools sort -i meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval001_unsort_unique.gtf > meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval001_sort.gtf

```
*NOTE: `bedtools intersect -wa -u` write `a` if overlap with `b`; only once*

--> `meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_sort.gtf` is file with DEGs Down HET and Up KO but still contain intergenic peak gene

--> The gtf files looks good in term of expression changes (checked on IGV)

Now let's generate **gtf file that contains non-intergenic peaks; in at least 1 genotype**:

```bash
# collect all annotated genes for each genotype
cat output/ChIPseeker/annotation_WT.txt output/ChIPseeker/annotation_HET.txt output/ChIPseeker/annotation_KO.txt | sort | uniq > output/ChIPseeker/annotation.txt

# Filter out Intergenic
grep -v "Intergenic" output/ChIPseeker/annotation.txt > output/ChIPseeker/annotation_noIntergenic.txt


# Collect geneSymbol (gene name) list
## Print 20th column in each rows; sort; remov dupplicates
awk -F'\t' '(NR==1 || FNR>1) {print $20}' output/ChIPseeker/annotation_noIntergenic.txt | sort | uniq > output/ChIPseeker/annotation_noIntergenic_geneSymbol.txt

# Modify the .txt file that list all genes so that it match gtf structure
sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_noIntergenic_as_gtf_geneSymbol.txt
# Filter the gtf
grep -Ff output/ChIPseeker/annotation_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_peak_noIntergenic.gtf
```
- *NOTE: awk process file line per line; `-F'\t'` = tab-separated; `(NR==1 || FNR>1)` tell that header line and all lines after are processed; `{print $20}` will only print the 20th column=gene name `sort` sort output because `uniq` only remove adjacent dupplicated rows*
- *NOTE: `sed` here add `gene name "` before and `"` after*
- *NOTE: `-Ff` tel grep to read the pattern in each row of the file*


Let's combine the **Up/down DEGs genes with the peak non intergenic GTF** and deepTools:


```bash
## Combine DEGs and peak non intergenic
bedtools intersect -wa -u -a meta/ENCFF159KBI_peak_noIntergenic.gtf -b meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_sort.gtf > meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_sort.gtf

bedtools intersect -wa -u -a meta/ENCFF159KBI_peak_noIntergenic.gtf -b meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_sort.gtf > meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_sort.gtf

bedtools intersect -wa -u -a meta/ENCFF159KBI_peak_noIntergenic.gtf -b meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval01_sort.gtf > meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_qval01_sort.gtf

bedtools intersect -wa -u -a meta/ENCFF159KBI_peak_noIntergenic.gtf -b meta/ENCFF159KBI_DEGs_HET_Down_KO_Up_qval001_sort.gtf > meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_qval001_sort.gtf

bedtools intersect -wa -u -a meta/ENCFF159KBI_peak_noIntergenic.gtf -b meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval01_sort.gtf > meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_qval01_sort.gtf

bedtools intersect -wa -u -a meta/ENCFF159KBI_peak_noIntergenic.gtf -b meta/ENCFF159KBI_DEGs_HET_Up_KO_Down_qval001_sort.gtf > meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_qval001_sort.gtf


# deepTools plot
conda activate deeptools
## DEGs down HET and up in KO and peak non intergenic
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_noIntergenic_DEGs_HET_Down_KO_Up.sh # 681175
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_noIntergenic_DEGs_HET_Up_KO_Down.sh # 1308054

sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_noIntergenic_DEGs_HET_Down_KO_Up.sh # 681177

```

--> looks cool; seems enough genes and we clearly see up HET and down KO


Then, `bedools intersect` the **DEGs GTF with the DiffBind GTF** and generate the deepTools plots:

```bash
conda activate BedToBigwig
## DEG down HET and up in KO and diffbind up HET and down in KO
bedtools intersect -wa -u -a meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_sort.gtf -b meta/ENCFF159KBI_UpHET_DownKO_noIntergenic.gtf > meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.gtf

# deepTools plot
conda activate deeptools
## Up HET and Down in KO
sbatch scripts/matrix_gene_1kb_DiffBind_TMM_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.sh # 681256
sbatch scripts/matrix_TSS_5kb_DiffBind_TMM_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.sh # 681257
```

--> looks cool; maybe not enough genes and we clearly see up HET and down KO


Let's **count the nb of unique genes in various our gtf**:
```bash
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_sort.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_UpHET_DownKO_noIntergenic.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_WT_maxScorePoolpeaks_min1_noIntergenic.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_sort.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l

awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_qval01_sort.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_qval001_sort.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
awk -F'\t' '{split($9,a,";"); for(i in a) if(a[i] ~ /gene_id/) print a[i]}' meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_qval001_sort.gtf | tr -d ' ' | tr -d '\"' | sort | uniq | wc -l
```
nb of unique genes:
- meta/ENCFF159KBI_WT_maxScorePoolpeaks_min1_noIntergenic.gtf: 4,665 (peak in WT)
- meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_sort.gtf: 2,100 (DEGs)
- meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_qval01_sort.gtf: 1,517 (DEGs)
- meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_qval001_sort.gtf: 1,009 (DEGs)
- meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_sort.gtf: 1,818 (DEGs)
- meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_qval01_sort.gtf: 1,393 (DEGs)
- meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Up_KO_Down_qval001_sort.gtf: 1,047 (DEGs)
- meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.gtf: 28 (DEGs + Diff. bound)
- meta/ENCFF159KBI_UpHET_DownKO_noIntergenic.gtf: 55 (Diff. bound)

### Functional analyses for the expected gene list

Let's see what are the genes where H3K27me3 goes Up in HET and Down in KO; in agreement with expression changes

```bash
conda activate deseq2
```

```R
# library
library("rtracklayer")
library("tidyverse")
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("DOSE")
library("pathview")
library("enrichplot")


# Collect the gene list from the gtf
## Import the GTF file
gtf <- import("meta/ENCFF159KBI_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.gtf")
## Extract gene names and convert to Entrez ID
gene_symbols <- unique(elementMetadata(gtf)$gene_name) # = elementMetadata(gtf)$gene_name %>% unique() 

genes <- mapIds(org.Hs.eg.db, 
                   keys = gene_symbols, 
                   column = "ENTREZID", 
                   keytype = "SYMBOL", 
                   multiVals = "first")



# Functional profiles
## KEGG
enrichKEGG <- enrichKEGG(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/ChIPseeker/functional_KEGG_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.pdf", width=7, height=5)
dotplot(enrichKEGG, showCategory = 13, title = "KEGG Pathway Enrichment Analysis")
dev.off()

pdf("output/ChIPseeker/emapplot_functional_KEGG_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.pdf", width=10, height=8)
emapplot(pairwise_termsim(enrichKEGG), showCategory = 13)
dev.off()

## GO
enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_BP_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.pdf", width=7, height=15)
dotplot(enrichGO, showCategory = 1, title = "GO_Biological Process Enrichment Analysis")
dev.off()
### NO ENRICHMENT


enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "CC") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_CC_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.pdf", width=7, height=15)
dotplot(enrichGO, showCategory = 1, title = "GO_Cellular Component Enrichment Analysis")
dev.off()
### 1 enrichment for Potassium Channel Complex

enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "MF") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)
pdf("output/ChIPseeker/functional_GO_MF_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.pdf", width=7, height=15)
dotplot(enrichGO, showCategory = 15, title = "GO_Molecular Function Enrichment Analysis")
dev.off()
### NO ENRICHMENT

## Disease
enrichDO <- enrichDO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
pdf("output/ChIPseeker/functional_DO_peak_noIntergenic_DEGs_HET_Down_KO_Up_DiffBind_UpHET_DownKO.pdf", width=7, height=15)
dotplot(enrichDO, showCategory = 15, title = "Disease Ontology Enrichment Analysis")
dev.off()
### NO ENRICHMENT



```

--> Not interesting in term of GO... Some hits for KEGG; and GO Cellular Component Potassium Channel Complex (GRIK2, DPP6, KCNQ1)




## bigwig_DiffBind_TMM_ratio


```bash
conda activate deeptools
# Replicates
sbatch --dependency=afterany:49458 scripts/matrix_TSS_5kb_WT_DiffBind_TMM_ratio.sh # 49467 ok
sbatch --dependency=afterany:49467 scripts/matrix_TSS_5kb_WT_DiffBind_TMM_ratio_profile.sh # 49468 ok

sbatch --dependency=afterany:49458 scripts/matrix_TSS_5kb_HET_DiffBind_TMM_ratio.sh # 49469 ok
sbatch --dependency=afterany:49469 scripts/matrix_TSS_5kb_HET_DiffBind_TMM_ratio_profile.sh # 49470 ok

sbatch --dependency=afterany:49458 scripts/matrix_TSS_5kb_KO_DiffBind_TMM_ratio.sh # 49471 ok
sbatch --dependency=afterany:49471 scripts/matrix_TSS_5kb_KO_DiffBind_TMM_ratio_profile.sh # 49472 ok

# Genotype TSS (5kb) 
sbatch --dependency=afterany:49458 scripts/matrix_TSS_5kb_DiffBind_TMM_ratio.sh # 49473 ok
sbatch --dependency=afterany:49473 scripts/matrix_TSS_5kb_DiffBind_TMM_ratio_profile.sh # 49475 ok

# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch --dependency=afterany:49458 scripts/matrix_gene_1kb_DiffBind_TMM_ratio.sh # 49477 ok
sbatch --dependency=afterany:49477 scripts/matrix_gene_1kb_DiffBind_TMM_ratio_profile.sh # 49478 ok
```

--> Replicate are quite OK but let's forget it as the NEUROG2 region looks shit

## bigwig_DiffBind_TMM_subtract


```bash
conda activate deeptools
# Replicates
sbatch --dependency=afterany:49460 scripts/matrix_TSS_5kb_WT_DiffBind_TMM_subtract.sh # 49479 ok
sbatch --dependency=afterany:49479 scripts/matrix_TSS_5kb_WT_DiffBind_TMM_subtract_profile.sh # 49480 ok

sbatch --dependency=afterany:49460 scripts/matrix_TSS_5kb_HET_DiffBind_TMM_subtract.sh # 49482 ok
sbatch --dependency=afterany:49482 scripts/matrix_TSS_5kb_HET_DiffBind_TMM_subtract_profile.sh # 49483 ok

sbatch --dependency=afterany:49460 scripts/matrix_TSS_5kb_KO_DiffBind_TMM_subtract.sh # 49485 ok
sbatch --dependency=afterany:49485 scripts/matrix_TSS_5kb_KO_DiffBind_TMM_subtract_profile.sh # 49486 ok

# Genotype TSS (5kb) 
sbatch --dependency=afterany:49460 scripts/matrix_TSS_5kb_DiffBind_TMM_subtract.sh # 49487 ok
sbatch --dependency=afterany:49487 scripts/matrix_TSS_5kb_DiffBind_TMM_subtract_profile.sh # 49488 ok

# Genotype gene body (-1 / +1 kb - TSS / TES)
sbatch --dependency=afterany:49460 scripts/matrix_gene_1kb_DiffBind_TMM_subtract.sh # 49489 ok
sbatch --dependency=afterany:49489 scripts/matrix_gene_1kb_DiffBind_TMM_subtract_profile.sh # 49490 ok
```

--> Very bad! Some values are negative!!!

--> Seems that IgG is already taken into account; or I failed the IgG scaling factor to be used in the bigwig

--> DiffBind_TMM non IgG corr seems the best to use





