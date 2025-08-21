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
sbatch scripts/fastp.sh # 50212025 ok
```

## mapping fastp trim

```bash
sbatch --dependency=afterany:50212025 scripts/STAR_mapping_fastp.sh # 50212217 ok
```

-->  xxx


# Count with Kallisto

Follow instruction [here](https://pachterlab.github.io/kallisto/download.html)

## count with Kallisto

```bash
conda activate kallisto


## run in sbatch
sbatch scripts/kallisto_count_gtf.sh # 50212654 FAIL SHOULD USE STRANDED!; 50302968 xxx

# Convert pseudoalignment to bigwig
sbatch --dependency=afterany:50302968 scripts/TPM_bw.sh # 50306142 xxx
```

- *NOTE: Added `--rf-stranded --genomebam` options for strandness and pseudobam alignemt generation*



## Gene quantification with txImport


XXXY HERE RUN BELOW!



Follow [this](https://nbisweden.github.io/workshop-RNAseq/2011/lab_kallisto.html#2_Quantification)


```bash
# Extract all transcriptnames (1st) and genenames (4th) from  GTF and write to a file.   
awk 'BEGIN{OFS=","; print "TXNAME,GENEID"}
     $3=="transcript"{
       match($0,/transcript_id "([^"]+)"/,t);
       match($0,/gene_id "([^"]+)"/,g);
       if(t[1]!="" && g[1]!=""){
         tx=t[1]; gn=g[1];
         sub(/\.[0-9]+$/,"",tx);
         sub(/\.[0-9]+$/,"",gn);
         print tx, gn;
       }
     }' ../../Master/meta/gencode.v47.annotation.gtf \
| sort -u > ../../Master/meta/gencode.v47.annotation.tx2gene.csv


awk 'BEGIN{OFS=","; print "TXNAME,GENEID"}
     $3=="transcript"{
       match($0,/transcript_id "([^"]+)"/,t);
       match($0,/gene_id "([^"]+)"/,g);
       if(t[1]!="" && g[1]!="") print t[1], g[1];
     }' ../../Master/meta/gencode.v47.annotation.gtf \
| sort -u > ../../Master/meta/gencode.v47.annotation.tx2gene.csv



conda activate deseq2
```

Go in R to create metadata file; and convert transcript to gene count

Metadata file format as: SampleName, SampleID, No, Model, Day, Group, Replicate


```R
# packages
library("tidyverse")
library("dplyr") # data wrangling
library("ggplot2") # plotting
library("DESeq2") # rna-seq
library("tximport") # importing kalisto transcript counts to geneLevels
library("readr") # Fast readr of files.
library("rhdf5") # read/convert kalisto output files.  


# Create metadata
samples <- c(
  "ESC_WT_R1","ESC_KO_R1","ESC_OEKO_R1",
  "ESC_WT_R2","ESC_KO_R2","ESC_OEKO_R2",
  "ESC_WT_R3","ESC_KO_R3","ESC_OEKO_R3"
)
mr <- data.frame(
  SampleID   = samples,
  No         = 1:length(samples),
  Model      = "ESC",
  Day        = "ESC",
  Group      = sub("ESC_","", sub("_R[0-9]","", samples)),
  Replicate  = sub(".*_R","", samples),
  row.names  = samples         # <-- set SampleName as rownames
)
mr


# List all abundance.tsv files
files <- list.files(
  path = "output/kallisto",
  pattern = "abundance.tsv",
  recursive = TRUE,
  full.names = TRUE
)
# Name the files with your sample names
names(files) <- rownames(mr)
files



# Convert transcript to gene ID
tx2gene <- read_csv("../../Master/meta/gencode.v47.annotation.tx2gene.csv")
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

```


























