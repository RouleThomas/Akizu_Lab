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
sbatch scripts/fastp_raw.sh # 14681036 xxx
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

--> Overall >xxx% input reads as been uniquely mapped to the genome (90% non uniq)



## Removing dupplicates (only uniquely aligned reads)
This is prefered for THOR bam input.


```bash
conda activate bowtie2

sbatch --dependency=afterany:14681230 scripts/samtools_unique_NPC.sh # 14681286 NODE failure
sbatch --dependency=afterany:14681246 scripts/samtools_unique_53dN.sh # 14681396 NODE failure


sbatch scripts/samtools_unique_NPC_1.sh # 15101660 fail; 15116487 ok
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

sbatch --dependency=afterany:15101878 scripts/bamtobigwig_unique_53dN_1.sh # 15102289 xxx
sbatch scripts/bamtobigwig_unique_NPC_1.sh # 15102387 fail; 15116725 fail; 15137560 xxx

```


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX BELOW NOT MD XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



- 50dN native and FA
*Pass*: 
*Failed*: H3K27me3, EZH1cs, EZH2
- NPC _ WT
*Pass*: H3K27ac, H3K4me3
*Failed*: EZH1cs, EZH2, SUZ12
- NPC _ KOEF1aEZH1
*Pass*: H3K27me3, H3K4me3, H
*Failed*: EZH1cs, SUZ12


--> Non unique (all raw reads!) vs unique bigwig (less signal or more noise?): Very similar. increase signal on the one that work but more bakground


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

```

--> Hard to conclude stuff



# MACS2 peak calling on bam unique



--> The **peaks are called on the uniquely aligned reads** (it performed better on our previous CutRun)

**PEAK CALLING  in `broad`**


```bash
conda activate macs2
# genotype per genotype
sbatch --dependency=afterany:13343264 scripts/macs2_broad_50dN.sh # 12654885 ok missed sample added; 13343914 ok
sbatch --dependency=afterany:13343264 scripts/macs2_broad_KO_KOEF1aEZH1.sh # 12655091 ok missed sample added; 13344031 ok
sbatch --dependency=afterany:13343264 scripts/macs2_broad_WT.sh # 12655165 ok missed sample added; 13344112 ok

sbatch --dependency=afterany:13343264 scripts/macs2_narrow_50dN.sh # 12655462 ok missed sample added; 13344314 ok
sbatch --dependency=afterany:13343264 scripts/macs2_narrow_KO_KOEF1aEZH1.sh # 12655472 ok missed sample added; 13344540 ok
sbatch --dependency=afterany:13343264 scripts/macs2_narrow_WT.sh # 12655475 ok missed sample added; 13344565 ok

```

--> OEF1aEZH1 in 50-day neurons: too noisy for EZH1, EZH2, and H3K27me3

--> NPC histone marks: OK for H3K4me3, H3K27ac, and H3K27me3; sharp and clear peaks.

--> NPC PRC2 components: too noisy...

*- NOTE: peak calling has been run 2 times adding the missing samples!*

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













