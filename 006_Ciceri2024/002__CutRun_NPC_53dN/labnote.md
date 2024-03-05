# Project and goals 

Re-analysis of CutRun dataset from (Ciceri et al)[https://www.nature.com/articles/s41586-023-06984-8].

- Identify histone -mod bound genes in WT and compare with our data


# Download data


- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject PRJNA803355 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 14679607 ok

```



## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename_002.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_002.txt
```

--> All good 





# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 14681036 ok
```



# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:14681036 scripts/bowtie2_NPC.sh # 14681230 ok
sbatch --dependency=afterany:14681036 scripts/bowtie2_53dN.sh # 14681246 ok
```

--> Looks good; overall ~30-80% uniquely aligned reads
----> Seems less uniquel mapped reads than us but they sequence FAR more depth (~20m reads vs 5 for us)


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-14681230.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_14681230.txt

for file in slurm-14681246.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_14681246.txt


```

Add these values to `/home/roulet/006_Ciceri2024/002__CutRun_NPC_53dN/samples_002.xlsx`\
Then in R; see `/home/roulet/006_Ciceri2024/006_Ciceri2024.R`.

--> Overall >60% input reads as been uniquely mapped to the genome (90% non uniq)



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.


```bash
conda activate bowtie2

sbatch --dependency=afterany:14681230 scripts/samtools_unique_NPC.sh # 14681286 NODE failure
sbatch --dependency=afterany:14681246 scripts/samtools_unique_53dN.sh # 14681396 NODE failure


sbatch scripts/samtools_unique_NPC_1.sh # 15101660 fail; 15116487 ok
sbatch scripts/samtools_unique_NPC_2.sh # 15275393 ok
sbatch scripts/samtools_unique_53dN_1.sh # 15101878 ok

```


# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:14681286 scripts/bamtobigwig_unique_NPC.sh # 14681462 node failure
sbatch --dependency=afterany:14681396 scripts/bamtobigwig_unique_53dN.sh # 14681480 node failure; 

sbatch --dependency=afterany:15101878 scripts/bamtobigwig_unique_53dN_1.sh # 15102289 ok
sbatch scripts/bamtobigwig_unique_NPC_1.sh # 15102387 fail; 15116725 fail; 15137560 ok
sbatch --dependency=afterany:15275393 scripts/bamtobigwig_unique_NPC_2.sh # 15275543 ok



```



- NPC
PASS: H3K4m3 (rep very diff.), H3K9me3 (a bit noisy), H3K27me3
FAIL: H3K27ac (very low signal and noisy, seems R2 work better)
- 53dN
PASS: H3K4me3, H3K9me3 (a bit noisy), H3K27ac, H3K27me3
FAIL: *H3K9me3* could be there


--> The failed one, are also failed in the bigwig Ciceri files...

--> Compare with H3K27ac from (ENCODE)[https://www.encodeproject.org/experiments/ENCSR799SRL/]


```bash
# gunzip the file

# Convert wig to bigwig
srun --mem=500g --pty bash -l

## install wigtobigwig
conda activate BedToBigwig
conda install bioconda::ucsc-wigtobigwig # fail
conda install bioconda/label/cf201901::ucsc-wigtobigwig  # fail
#### --> fail create a new conda env
conda create -n wigtobigwig -c bioconda ucsc-wigtobigwig
conda activate wigtobigwig

## convert wig to bigwig
wigToBigWig output/bigwig_hg19/GSM767343_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.SK504.wig ../../Master/meta/hg19.chrom.sizes output/bigwig_hg19/GSM767343_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.SK504.bw
wigToBigWig output/bigwig_hg19/GSM818031_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.AK220.wig ../../Master/meta/hg19.chrom.sizes output/bigwig_hg19/GSM818031_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.AK220.bw
wigToBigWig output/bigwig_hg19/GSM896162_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.AK319.wig ../../Master/meta/hg19.chrom.sizes output/bigwig_hg19/GSM896162_UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K27ac.AK319.bw


```
NOTE: hg19 chrom size copy from [ucsc](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/)



## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools

# Generate compile bigwig (.npz) files _ hg38 Akizu analysis
sbatch scripts/multiBigwigSummary_NPC.sh # 15280856 ok
sbatch scripts/multiBigwigSummary_53dN.sh # 15280864 ok

# Generate compile bigwig (.npz) files _ hg19 Ciceri analysis
sbatch scripts/multiBigwigSummary_NPC_Ciceri.sh # 15280901 ok
sbatch scripts/multiBigwigSummary_53dN_Ciceri.sh # 15281092 ok


# Akizu with Ciceri _ NPC
sbatch scripts/multiBigwigSummary_NPC_CutRun001008_Ciceri.sh # 15290030 ok


```

**Good to use**:
- *NPC*: H3K27me3, H3K4me3, H3K9me3, IGG
- *53dN*: H3K27me3, H3K9me3

**Bad to use**:
- *NPC*: H3K27ac: no signal! Like IGG
- *53dN*: AB mix between IGG, H3K4me3, H3K27ac

--> Same observation between Akizu and Ciceri analysis(hg38 hg19)

--> The *good to use* ones nicely correlate with our data in NPC WT `CutRun__005008` (`CutRun__009`)



# MACS2 peak calling on bam unique



--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch scripts/macs2_broad_53dN.sh # 15401244 ok
sbatch scripts/macs2_broad_NPC.sh # 15401308 ok

sbatch scripts/macs2_broad_53dN_noIGG.sh # 15401429 ok
```

--> H3K27ac in NPC show 3 peak in R1! And a ~7k peaks in R2. R2 is better, but still very ugly and noisy!

--> 53dN IGG vs not using IGG: almost the same, so let's better use IGG






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



# deepTool plots


On all genes


```bash
conda activate deeptools

# All genes all histone marks
## NPC
sbatch scripts/matrix_TSS_10kb_NPC_raw_allGenes.sh # 15402417 ok
sbatch scripts/matrix_TSS_5kb_NPC_H3K27ac_raw_allGenes.sh # 15402511 ok
sbatch scripts/matrix_TSS_10kb_NPC_H3K27me3_H3K4me3_raw_allGenes.sh # 15402602 ok


## 53dN
sbatch scripts/matrix_TSS_10kb_53dN_raw_allGenes.sh # 15402437 ok
sbatch scripts/matrix_TSS_5kb_53dN_H3K27ac_raw_allGenes.sh # 15402531 ok
sbatch scripts/matrix_TSS_10kb_53dN_H3K27me3_H3K4me3_raw_allGenes.sh # 15402648 ok


# Akizu and Ciceri H3K27me3, H3K4me3, IGG
sbatch scripts/matrix_TSS_10kb_H3K27me3_H3K4me3_CutRun001008_Ciceri_raw_allGenes.sh # 15432941 ok


```


--> signal is very poor for H3K27ac in NPC as compared to 53dN, notably for R1, almost like IGG...
















