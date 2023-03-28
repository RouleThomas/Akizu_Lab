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
	mutate(target_norm = (sum_read/total_read)*100) %>%


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
spikein_all_scale %>%
  	ggplot(aes(x = Target, y = scaled_target_norm, fill = AB)) +
  	geom_col(position = "dodge") +
  	facet_wrap(~sample_ID) +
  	geom_hline(yintercept = 20, color = "red", linetype = "longdash") +
  	theme_bw() +
  	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
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

--> fastqc are XXX


# Mapped clean reads

XXX
