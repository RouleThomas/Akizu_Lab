# Project

**H9 cell lines**
- ESC native:
    - WT: 3 Bio Rep (A1-3)
    - KO: 3 Bio Rep (A4-6)
    - KOEF1aEZH1: 3 Bio Rep (A7-9)


--> Directional mRNA library preparation (poly A enrichment), NovaSeq X Plus Series (PE150)




**Objectives:**
- Put together with CutRun 018


Novogene Input	Sample Name	
A1	R1 EZH1 WT	ESC_WT_R1
A2	R1 EZH1 KO	ESC_KO_R1
A3	R1 EZH1 KO (50 ng dox)	R1 EZH1 OEKO    ESC_OEKO_R1
A4	R2 EZH1 WT	    ESC_WT_R2
A5	R2 EZH1 KO	    ESC_KO_R2
A6	R2 EZH1 KO (50 ng dox)	R2 EZH1 OEKO    ESC_OEKO_R2
A7	R3 EZH1 WT	    ESC_WT_R3
A8	R3 EZH1 KO	    ESC_KO_R3
A9	R3 EZH1 KO (50 ng dox)	R3 EZH1 OEKO    ESC_OEKO_R3





# Pipeline
- Download data (wget)
- Rename files
- FastQC (fastqc)
- Trimming (fastp)
- Count with Kallisto (better than Salmon for switch variation analysis, and better than featureCounts; best to do transcriptome level coutning and then to gene
)

--> Detail of the overall pipeline in `Meeting_20230919_draft.xlsx` 

# Download / import data


```bash
# Following email instructions
module load lftp
lftp -c 'set sftp:auto-confirm yes;set net:max-retries 20;open sftp://X202SC25078875-Z01-F001:wfcxcay4@usftp23.novogene.com; mirror --verbose --use-pget-n=8 -c'

# Copy all .fz.gz data into input/ folder
rsync -av --include '*/' --include '*.fq.gz' --exclude '*' usftp23.novogene.com/ input/ # copy from usftp23 folder to input
find input/ -mindepth 2 -type f -exec mv -t input/ {} + # mv files from their folder to input/ folder
find input/ -type d -empty -delete # delete empty directory

```

--> All good, files created in `usftp23.novogene.com/`




# Rename file

Renamed manually as only 8 samples

--> All good 



# Fastp cleaning

```bash
sbatch scripts/fastp.sh # 50212025 xxx
```



## mapping fastp trim

```bash
sbatch --dependency=afterany:50212025 scripts/STAR_mapping_fastp.sh # 50212217 xxx
```

-->  xxx






# Count with Kallisto

Follow instruction [here](https://pachterlab.github.io/kallisto/download.html)

## count with Kallisto

```bash
conda activate kallisto


## run in sbatch
sbatch --dependency=afterany:50212217 scripts/kallisto_count_gtf.sh # 50212654 xxx
```








