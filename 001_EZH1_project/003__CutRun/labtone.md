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
sbatch scripts/bowtie2_patient.sh # 11826852 XXX
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


```bash
sbatch scripts/samtools_HET.sh # 11578283 ok
sbatch scripts/samtools_KO.sh # 11578284, weirdly looong
sbatch scripts/samtools_WT.sh # 11578286 ok

sbatch scripts/samtools_patient.sh # XXX
```
--> `scripts/samtools_KO.sh` stop running for unknown reason. Or maybe just super-long, it is stuck at the `8wN_KO_IGG_R2` sample. Let's run again the whole script, with more memory and threads, and without sample 8wN_KO_IGG_R1 just to make sure the script is working. **Output in `output/tmp`**
```bash
sbatch scripts/samtools_KO_2.sh # ok
```
--> The long job has been cancel and samtools_KO_2.sh succesfully run.

--> New files transfered to the `output/bowtie2`. All is complete

# Generate wig coverage files
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch scripts/bamtobigwig_WT.sh # 11827228
sbatch scripts/bamtobigwig_HET.sh # 11827232
sbatch scripts/bamtobigwig_KO.sh # 11827233

sbatch scripts/bamtobigwig_patient.sh  # XXX
```


# Peak calling

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
sbatch scripts/bamtobedgraph_WT_1.sh # 11826612 ok
sbatch scripts/bamtobedgraph_WT_2.sh # 11826613 ok
sbatch scripts/bamtobedgraph_HET_1.sh # 11826580 ok
sbatch scripts/bamtobedgraph_HET_2.sh # 11826581 ok
sbatch scripts/bamtobedgraph_KO_1.sh # 11826578 ok
sbatch scripts/bamtobedgraph_KO_2.sh # 11826579 ok

sbatch scripts/bamtobedgraph_patient.sh # XXX
```
*NOTE: At bedtools bamtobed warning that some reads have no mate. But that is because I filtered reads based on MAPQ and removed unmapped reads, secondary alignments, and reads failing quality check using samtools view command, so I might have removed one read from a pair while retaining the other*

**Run SCEA peak calling**:
```bash
# Example command
bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output

# Run together
## Stringeant
sbatch scripts/SEACR_WT.sh # 11826904
sbatch scripts/SEACR_HET.sh # 11826911
sbatch scripts/SEACR_KO.sh # 11826912

sbatch scripts/SEACR_patient.sh # XXX

## Run all samples with relax (this is pretty fast)
sbatch scripts/SEACR_relax.sh # 11826920

sbatch scripts/SEACR_relax_patient.sh # XXX


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


## MACS2 peak calling







