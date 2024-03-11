# Project and goals 

H9 (WA09) at d28
Inhibitor: DOT1Linh, EHMTinh, EZH2inh
IP: H3K27me3, H3K4me3

Re-analysis of CutRun dataset from (Ciceri et al)[https://www.nature.com/articles/s41586-023-06984-8].

EZH2 inhibitors should also target EZH1. So could be considered as an EZH1 and EZH2 inhibitors.
- Identify genes that lose H3K27me3 with inhibitor treatment = putative EZH2 target
- Check whether putative EZH2 target are the same genes that gain H3K27me3 in our EZH1 KO (if yes, would suggest our genes that gain H3K27me3 indeed gain it because of increase EZH2 activity)



# Download data




- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject PRJNA803355 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 16188011 ok

```



## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename_003.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_003.txt
```

--> All good 





# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 16197773 xxx
```



# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:16197773 scripts/bowtie2_1.sh # 16198162 xxx
sbatch --dependency=afterany:16197773 scripts/bowtie2_2.sh # 16198215 xxx
```

-->  XXX Looks good; overall ~30-80% uniquely aligned reads XXX
----> Seems less uniquel mapped reads than us but they sequence FAR more depth (~20m reads vs 5 for us)


## Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-16198162.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_16198162.txt

for file in slurm-16198215.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2/alignment_counts_16198215.txt


```

Add these values to `/home/roulet/006_Ciceri2024/003__CutRun_EpigeneticInhibitors/samples_003.xlsx`\
Then in R; see `/home/roulet/006_Ciceri2024/006_Ciceri2024.R`.

--> XXX Overall >60% input reads as been uniquely mapped to the genome (90% non uniq) XXX



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.


```bash
conda activate bowtie2

sbatch --dependency=afterany:16198162 scripts/samtools_unique_1.sh # 16198634 xxx
sbatch --dependency=afterany:16198215 scripts/samtools_unique_2.sh # 16198636 xxx


```


# Generate bigwig coverage files
## Raw bigwig
Paramaters:
- `--binSize 1` for good resolution
- `--scaleFactor 0.5` to obtain the exact number of reads respective to the bam, otherwise it count two instead of 1
- `--extendReads` Reads extented taking into account mean fragment size of all mated reads.

```bash
conda activate deeptools

sbatch --dependency=afterany:16198634 scripts/bamtobigwig_unique_1.sh # 16198742 xxx
sbatch --dependency=afterany:16198636 scripts/bamtobigwig_unique_2.sh # 16198743 xxx

```



- EZH2inh
PASS: 
FAIL: 
- EHMTinh
PASS: 
FAIL: 
- DOT1Linh
PASS: 
FAIL: 
- DMSO
PASS: 
FAIL: 





## Pearson correlation heatmap on bigwig signals



```bash
conda activate deeptools

# Generate compile bigwig (.npz) files _ hg38 Akizu analysis
sbatch --dependency=afterany:16198742:16198743 scripts/multiBigwigSummary_all.sh # 16200148 xxx



```

XXXXXXXXXXXXX CHUI AL below not mod XXXXXXXXXXXXX

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







