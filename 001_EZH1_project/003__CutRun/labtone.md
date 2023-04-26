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

Let's generate bigwig taking into account scaling factor:
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


What is wrong is that I should have use the reciprocal (1/n) and not n as scaling factor... Let' correct and save into output/bigwig_histone_NotGenotypeGroup (This is the true bigwig good to use, better than the *groupABgenotype* one)

```bash
conda activate deeptools
sbatch scripts/bamtobigwig_histone_scaled_WT_reciprocal.sh # 12370048 ok
sbatch scripts/bamtobigwig_histone_scaled_HET_reciprocal.sh # 12370046 ok
sbatch scripts/bamtobigwig_histone_scaled_KO_reciprocal.sh # 12370044 ok
sbatch scripts/bamtobigwig_histone_scaled_patient_reciprocal.sh # 12370043 ok
```




# Peak calling_Version1 with failed normalized data

## SEACR peak calling
[Paper](https://doi.org/10.1186/s13072-019-0287-4) and [Github](https://github.com/FredHutch/SEACR)

Install required packages and prepare meta files
```bash
# Install bedtools within bowtie2 conda env
conda activate bowtie2
conda install -c bioconda bedtools
module load sam-bcf-tools/1.6

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
sbatch scripts/bamtobigwig_ChIPSeqSpike.sh # 12287902 XXX
```



Re run ChIPSeqSpike:

XXX


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


### Histone lib depth


XXX



# ChIPseqSpikeInFree on CutRun

Let's see what is the prediction from ChIPseqSpikeInFree on our data

```bash
conda activate ChIPseqSpikeInFree
sbatch scripts/ChIPseqSpikeInFree.sh # 12379092
```

It does not work super well ("ok" correlated: r=0.4, p0.1). Possibly because CutRun is too clean, very few low proportion reads as compare to ChIPseq so may fail...



# ChIPseeker for binding profiles
[Tutorial](http://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html) and [documentation]()

For **ChIPseeker**, use **conda base and R 4.2.2 module**
```bash
conda deactivate
module load R/4.2.2
```

```R
library("ChIPseeker")

XXX

```


