# Project
- Neurons 50 days:
    - EF1a: EZH1, EZH2, H3K27me3, SUZ12, IGG
    - KO:  EZH1, EZH2, H3K27me3, SUZ12, IGG
    - Q731E (strong GOF):  EZH1, EZH2, H3K27me3, SUZ12, IGG
- PSC:
    - WT_ FA???: EZH1, H3K27me1
    - WT_ FA???: EZH1, H3K27me1


**Objectives:**
- Check whether Q731E indeed increase highly H3K27me3
- Check whether EZH1 AB work better in neurons
- test different FA condition in PSC for EZH1 AB 
- test different FA condition in PSC for H3K27me1 AB 




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

Go [there](http://data-deliver.novogene.com/batchfiles/X202SC23117109-Z01-F001) and enter credetnial: (check email Novogen)

I created a `nano url.txt` with all link and used `wget -i url.txt` to download them all (1 link per raw); then `mv input_raw_Novogene/*fq.gz input` .


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
sbatch scripts/fastp_1.sh # 9019614
sbatch scripts/fastp_2.sh # 9019615
sbatch scripts/fastp_3.sh # 9019616
sbatch scripts/fastp_4.sh # 9019617
```

# FastQC

**Fastp-cleaned:**
```bash
sbatch --dependency=afterany:9019614 scripts/fastqc_fastp_1.sh # 9019677
sbatch --dependency=afterany:9019615 scripts/fastqc_fastp_2.sh # 9019678
sbatch --dependency=afterany:9019616 scripts/fastqc_fastp_3.sh # 9019680
sbatch --dependency=afterany:9019617 scripts/fastqc_fastp_4.sh # 9019681
```

--> all good


# Mapping

Let's map with endtoend parameter as for `003__CutRun` (`--phred33 -q --no-unal --no-mixed --dovetail`)

```bash
conda activate bowtie2

sbatch --dependency=afterany:9019677 scripts/bowtie2_1.sh # 9020294
sbatch --dependency=afterany:9019678 scripts/bowtie2_2.sh # 9020295
sbatch --dependency=afterany:9019680 scripts/bowtie2_3.sh # 9020297
sbatch --dependency=afterany:9019681 scripts/bowtie2_4.sh # 9020298

```

--> Looks good


XXX




