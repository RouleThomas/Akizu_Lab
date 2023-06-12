# Pipeline used from previous analyses

- ChIP-seq reads were aligned using Bowtie2 with parameters: -q --local --no-mixed --no-unal --dovetail69
- Uniquely aligned reads with a mapping quality score >= 20 and concordant alignments were kept for downstream analysis using samtools -q 20 -f 0x2
- H3K27me3 peaks were called using MACS2 with parameters: --broad --keep-dup all -p 1e-5 -- broad-cutoff 1e-570
- ChIP-seq signal was normalized to sequencing depth using deepTools: bamCoverage --normalizeUsing CPM71
- ChIP-seq signal around peaks was computed using deepTools computeMatrix
- ChIP-seq signal around peaks was visualized using deepTools plotHeatmap function


# Import files from Google drive to the cluster


I cannot use a bash script to import files as if disconnection the transfer will fail. So cp the old-fashion way, let computer running o/n.
```bash
cp -r /home/roulet/tsclient/roule/Google\ Drive\ Streaming/Shared\ drives/akizulab/Primary\ Data/ChIPseqs/Fastq\ files/ESC\ NPC\ and\ Neurons/ Dec/ 2020/LopezN_FASTQ-216179970/* \
/scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq/input
``` 

# Concatenate fastq files
Here, each sample sequenced from 2 Illumina lanes, need concatenate them into one.

**Move all fastq files within the input folder** (now [file] are in input/folder1/folder2/[file]) using a loop:
```bash
for file in input/*/*/;
    do mv "$file"* input/;
done
```
- `input/*/*/` this go 2 folder downstream; whatever their names
- we move file to `input/` directory 

## Light test of concatenation with 1 sample to make sure we keep file integrity
```bash
# copy files as backup
## Read1 (lane1 lane2)
cp input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz input/backup/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz
cp input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz input/backup/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz
## Read2 (lane1 lane2)
cp input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz input/backup/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz
cp input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz input/backup/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz

# concatenate lane 
## Read1 (lane1 lane2)
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz > input/ESC_KO_input_R2_1.fq.gz
## Read2
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz > input/ESC_KO_input_R2_2.fq.gz
# fastqc
## Read1
fastqc -o output/fastqc input/ESC_KO_input_R2_1.fq.gz
fastqc -o output/fastqc input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz
## Read2
fastqc -o output/fastqc input/ESC_KO_input_R2_2.fq.gz
fastqc -o output/fastqc input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz
```
--> The number of reads is similar in the 2 (paired)-read; and reduced/2 when looking at only 1 of the 2 lanes; so concatenation keep data integrity and combine all reads. Just in case, backup has been performed from all files: **To be deleted after mapping**


## Concatenate the 2 lanes of all samples
```bash
# example for 1 file:
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R1_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R1_001.fastq.gz > input/ESC_KO_input_R2_1.fq.gz
cat input/10-ESCs-EZH1-KO-2-INPUT_S8_L001_R2_001.fastq.gz input/10-ESCs-EZH1-KO-2-INPUT_S8_L002_R2_001.fastq.gz > input/ESC_KO_input_R2_2.fq.gz

# time per time:
sbatch scripts/concat_1.sh # 11256971; fail with sample17, repeat it manually
sbatch scripts/concat_2.sh # 11256972 ok
sbatch scripts/concat_3.sh # 11256973 ok
```
--> `2dN_KO_input_R2_2.fq.gz` is missing... So I look for it `grep -rn '2dN_KO_input_R2' scripts/`, mistake corrected in concat_2, and re-run manually for both 2dN_KO_input_R1 and 2dN_KO_input_R2. Will double check all is good at the fastqc step



# Trimming and Fastqc
## Fastqc on raw reads
```bash
# example for 1 file:
fastqc -o output/fastqc input/2dN_HET_H3K27me3_R1_1.fq.gz

# grouped:
sbatch scripts/fastqc_1.sh # 11262560 ok 
sbatch scripts/fastqc_2.sh # 11262458 2dN_KO_input_R2_2 fastqc fail (file was empty)
sbatch scripts/fastqc_3.sh # 11262515 ok
```

--> Re-generate the 2dN_KO_input_R2:
```bash
# Copy backup
cp input/backup/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L001_R2_001.fastq.gz input/
cp input/backup/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L002_R2_001.fastq.gz input/

# Concatenate backup
cat input/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L001_R2_001.fastq.gz input/24-INPUT-Neurons--DID2--EZH1-KO-2_S34_L002_R2_001.fastq.gz > input/2dN_KO_input_R2_2.fq.gz

# fastqc
fastqc -o output/fastqc input/2dN_KO_input_R2_2.fq.gz
```
--> 2dN_KO_input_R1 and 2dN_KO_input_R2 are not the same, so we're good

ESC_KO_H3K27me3_R2_2 and ESC_KO_H3K27me3_R2_1 are the exact same one according to fastqc... I did `grep -rn 'ESC_KO_H3K27me3_R2' scripts` and found the cat command was wrong for ESC_KO_H3K27me3_R2_2 file. Let's regenerate ESC_KO_H3K27me3_R2_2:
```bash
# Copy backup
cp input/backup/4-ESCs-EZH1-KO-2-H3k27me3_S3_L001_R2_001.fastq.gz input/
cp input/backup/4-ESCs-EZH1-KO-2-H3k27me3_S3_L002_R2_001.fastq.gz input/

# Concatenate backup
cat input/4-ESCs-EZH1-KO-2-H3k27me3_S3_L001_R2_001.fastq.gz input/4-ESCs-EZH1-KO-2-H3k27me3_S3_L002_R2_001.fastq.gz > input/ESC_KO_H3K27me3_R2_2.fq.gz


# fastqc
fastqc -o output/fastqc input/ESC_KO_H3K27me3_R2_2.fq.gz
```
--> All good now



## Fastqc on fastp-trimmed-reads
### Trimming with Fastp
```bash
# example for 1 file:
fastp -i input/2dN_HET_H3K27me3_R1_1.fq.gz -I input/2dN_HET_H3K27me3_R1_2.fq.gz \
      -o output/fastp/2dN_HET_H3K27me3_R1_1.fq.gz -O output/fastp/2dN_HET_H3K27me3_R1_2.fq.gz \
	  -h output/fastp/2dN_HET_H3K27me3_R1 -j output/fastp/2dN_HET_H3K27me3_R1

# grouped:
sbatch fastp_ESC.sh # 11262881 ok
sbatch fastp_NPC.sh # 11262882 ok 
sbatch fastp_2dN.sh # 11262874, 2dN_KO_input_R2 fail
```
Let's repeat fastp-triming for ESC_KO_H3K27me3_R2_2:
```bash
fastp -i input/ESC_KO_H3K27me3_R2_1.fq.gz -I input/ESC_KO_H3K27me3_R2_2.fq.gz \
    -o output/fastp/ESC_KO_H3K27me3_R2_1.fq.gz -O output/fastp/ESC_KO_H3K27me3_R2_2.fq.gz \
	-h output/fastp/ESC_KO_H3K27me3_R2 -j output/fastp/ESC_KO_H3K27me3_R2
```


### Fastqc trimmed reads
Scripts adapted so that they can only start when there respective `fastp.sh` job is finish:
```bash
# grouped:
sbatch --dependency=afterany:11262881 scripts/fastqc_fastp_ESC.sh # 11264158 ok
sbatch --dependency=afterany:11262882 scripts/fastqc_fastp_NPC.sh # 11264163 ok
sbatch --dependency=afterany:11262874 scripts/fastqc_fastp_2dN.sh # 11264174, 2dN_KO_input_R2
```
--> Fail, I took raw fastq input for fastqc, not the fastp-trim reads, script corrected and re-launch:
```bash
sbatch scripts/fastqc_fastp_ESC.sh # 11313011 ok
sbatch scripts/fastqc_fastp_NPC.sh # 11313006 ok
sbatch scripts/fastqc_fastp_2dN.sh # 11313014 ok
```

--> Re-generate the 2dN_KO_input_R2 fastp and fastqc: 
```bash
sbatch fastp_fastqc_2dN_KO_input_R2.sh # 11312794
```
Let's repeat fastqc for fastp-triming ESC_KO_H3K27me3_R2_2:
```bash
fastqc -o output/fastqc/fastp output/fastp/ESC_KO_H3K27me3_R2_2.fq.gz # ok
fastqc -o output/fastqc/fastp output/fastp/ESC_KO_H3K27me3_R2_1.fq.gz # ok
```
--> Same nb of files in fastqc/raw and fastqc/fastp; read trim so size not homogeneous.


# Mapping GRCh38
## Install and setup pre-requisets
**Download and install Bowtie2**
Create a **Bowtie2; conda environment** 
```bash
conda create -n bowtie2 -c bioconda bowtie2 # command can be launch from anywhere (directory and node)
```

**Bowtie2 genome indexation**
```bash
conda activate bowtie2
cd /scr1/users/roulet/Akizu_Lab/Master/meta

# comand
bowtie2-build GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta bowtie2_genome_dir/GRCh38

# run into
sbatch bowtie2_index.sh # 11351555 ok (~3 hours)
```

**Picard and samtools**
```bash
module load picard/2.26.10-Java-15
module load sam-bcf-tools/1.6
module load SAMtools/1.16.1* # New cluster
```

## Mapping
### Test on 1 sample
Let's map 2dN_HET_H3K27me3_R1, let's run command per command in interactive for light test of the script
```bash
bowtie2 --phred33 -q \
	-x ../../Master/meta/bowtie2_genome_dir/GRCh38 \
        -S output/bowtie2/2dN_HET_H3K27me3_R1.sam \
        -1 output/fastp/2dN_HET_H3K27me3_R1_1.fq.gz  \
        -2 output/fastp/2dN_HET_H3K27me3_R1_2.fq.gz | \ # here we stream the sam output directly to samtools to avoid creating the file
samtools view -Sb - | \
samtools sort -o output/bowtie2/2dN_HET_H3K27me3_R1.bam
```
--> It seems to work but it last forever so I am not sure, double check log files!

--> Only sam file work, then it failed. 

### Mapping on all samples with Shuo parameter
Bowtie2 parameters:
Run all mapping (sam > sort sam and generate bam > removing of dupplicates with picard > bam indexation). Parameters (I keep the same as Shuo):
- `--local` does not require that the entire read align from one end to the other (I was not using this)
- `--no-mixed` Do NOT align read if the read pair do not align (I was not using this)
- `--no-unal` To save space: Suppress SAM records for reads that failed to align. (I was not using this)
- `--dovetail` To allow reads to align even if they overlap (I was not using this)

Samtools parameters:
- `-f 0x2` only includes properly paired reads (I was using `-F 772` to exclude unmapped reads, secondary alignments, and reads failing quality checks)
- `-q 20` (I was using 30)

```bash
sbatch bowtie2_map_ESC.sh # 11452938 (~10hrs per sample!) ok
sbatch bowtie2_map_NPC.sh # 11452939 (~10hrs per sample!) ok
sbatch bowtie2_map_2dN.sh # 11452937 (~10hrs per sample!) ok
```

--> The sam files have been well generated but then it failed. Re-run script from samtools:

--> Uniquely mapped reads is > 50% (which look ok according to data from [encode](https://academic.oup.com/bib/article/18/2/279/2453282))

### Compare mapping efficacy using different parameters

The more uniquely mapped reads, the better. Let's compare with the `2dN_HET_H3K27me3_R1.sam` sample
```bash
sbatch scripts/bowtie2_2dN_HET_H3K27me3_R1_param1.sh # bowtie2 default parameter # 11472979 ok
sbatch scripts/bowtie2_2dN_HET_H3K27me3_R1_param2.sh # parameter fine-tuned from Shuo; less stringeant # 11472978 ok
```

- Shuo parameter/**permissive-paired** `--phred33 -q --local --no-mixed --no-unal --dovetail`: (cannot have uniquely map paired reads)
    - nb of uniquely mapped reads: 23071868 (51.18%)
    - 154229070 (34.22%) >1 times  
    - 6580562 (14.6%) 0 times
    - overall 89.18%
- param1/**default** `--phred33 -q --no-unal`: (bowtie2 default) nb of uniquely mapped reads:
    - nb of uniquely mapped reads: 31938471 (70.85%)
    - 5743226 (12.74%) >1 times
    - 6580562 (38.23%) 0 times
    - overall 92.66%
- param2/**permissive-unpaired** `--phred33 -q --local --no-unal --dovetail`: (can have uniquely map paired reads)  
    - nb of uniquely mapped reads: 23071868 (51.18%)
    - 5743226 (34.22%) >1 times
    - 15429070 (12.74%) 0 times
    - overall 93.43%


--> Default is better! Better overall alignment rate, better concordant read alignment (with more uniquely mapped reads).

--> Probably because end-to-end method is used and not --local, aligning the entire read, not cliped may help the alignment. And we allow mix alignment

--> Let's try the following: `--phred33 -q --no-unal --no-mixed --dovetail` **endtoend**. Here it is more clean, we do not allow mix alignemnt and keep dovetail option possible.

```bash
sbatch bowtie2_2dN_HET_H3K27me3_R1_endtoend.sh # 11496376 ok
```
- param2/**endtoend** `--phred33 -q --no-unal --no-mixed --dovetail`:  
    - nb of uniquely mapped reads: 31927165 (70.82%)
    - >1 times 5761022 (12.78%)
    - overall 85.28%

*NOTE: All previous mapping have been moved to `/output/tmp/` folder. Can be deleted.*



## Mapping and read filtering with endtoend parameter
As I use `--no-mixed` the unpaired-mapped reads are already remove, so does not make sense to remove them again using -f 0x2 in samtools! Instead, let's use `-F 772` to exclude unmapped reads, secondary alignments (read that may map elsewhere), and reads failing quality checks.

```bash
# example for 1 file
sbatch scripts/samtools_2dN_HET_H3K27me3_R1.sh # 11506999 ok

# run job per time
sbatch scripts/bowtie2_samtools_2dN_1.sh # 11537670 ok
sbatch scripts/bowtie2_samtools_2dN_2.sh # 11537904 ok
sbatch scripts/bowtie2_samtools_ESC_1.sh # 11538307 cancel time limit; ESC_KO_input_R1
sbatch scripts/bowtie2_samtools_ESC_2.sh # 11538448 cancel time limit; ESC_WT_input_R2, ESC_WT_input_R3
sbatch scripts/bowtie2_samtools_NPC_1.sh # 11538814 ok
sbatch scripts/bowtie2_samtools_NPC_2.sh # 11538966 ok
```
--> The test example has worked well. Can be run for all samples once mapping is done.

Some samples stop due to time limit, rerun them:
```bash
sbatch scripts/bowtie2_samtools_canceled.sh # 11828352 ok
```

### Quality control metrics
Quality control plot (total read before trimming/ total read after trimming/ uniquely aligned reads)

Collect nb of reads from the slurm bowtie2 jobs:
```bash
for file in slurm-11537670.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2_endtoend/alignment_counts_11537670.txt

for file in slurm-11537904.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2_endtoend/alignment_counts_11537904.txt

for file in slurm-11538307.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2_endtoend/alignment_counts_11538307.txt

for file in slurm-11538448.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2_endtoend/alignment_counts_11538448.txt

for file in slurm-11538814.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2_endtoend/alignment_counts_11538814.txt

for file in slurm-11538966.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2_endtoend/alignment_counts_11538966.txt

for file in slurm-11828352.out; do
    total_reads=$(grep "reads; of these" $file | awk '{print $1}')
    aligned_exactly_1_time=$(grep "aligned concordantly exactly 1 time" $file | awk '{print $1}')
    aligned_more_than_1_time=$(grep "aligned concordantly >1 times" $file | awk '{print $1}')
    echo -e "$total_reads\t$aligned_exactly_1_time\t$aligned_more_than_1_time"
done > output/bowtie2_endtoend/alignment_counts_11828352.txt
```

Add these values to `/home/roulet/001_EZH1_project/002__CutRun/mapping_QC.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >80% input reads as been uniquely mapped to the genome


## Mapping and read filtering with endtoend parameter and higher samtools stringency (high quality mapping)

For ChIPSeqSpikeInFree to work best, we need very high quality reads; let's upgrade the MAPQ treshold from 20 to 30.


```bash
sbatch scripts/samtools_highquality_1.sh # 12345522 ok
sbatch scripts/samtools_highquality_2.sh # 12345523 ok
```

This work but I should also have filter the dupplicates. I did it later in the part **ChIPseqSpikeInFree - High quality reads or uniq mapp**





# Read depth normalization
Download jvarkit [here](https://github.com/lindenb/jvarkit). Transfer to `Master/software/` and `unzip JVARKIT.zip`; to use it, simply `java -jar jvarkit.jar`.

Use our `downsampleBAM.sh`; simply adapt the path to JVARKIT and sample names
```bash
sbatch scripts/downsampleBAM.sh # 11934723 ok
```

# ChIPseqSpikeInFree normalization
## ChIPseqSpikeInFree installation
See [github](https://github.com/stjude/ChIPseqSpikeInFree) and [paper](https://doi.org/10.1093/bioinformatics/btz720)

Create specific conda env for the installation
```bash
conda create --name ChIPseqSpikeInFree r-base=3.6.1 # 3.6.1 as in the R sessioninfo from github
conda activate ChIPseqSpikeInFree 
```
Follow Github-prerequisted for R packages installation.
```R
library(Rsamtools)
```
It failed in R when using `BiocManager::install`, lets try:
```bash
srun --mem=50g --pty bash -l
conda create --name ChIPseqSpikeInFree -c conda-forge -c bioconda r-base=3.6.1 bioconductor-rsamtools bioconductor-genomicranges bioconductor-genomicalignments
conda activate ChIPseqSpikeInFree 
conda install -c conda-forge r-usethis # Was needed to install devtools in R
conda install -c conda-forge r-gert # Was needed to install devtools in R

```
troubleshoot: I failed installing devtools in R; because of multiple missing packages, among them r-glue is problematic, need version =>1.6.1; I tried this `conda install -c conda-forge r-glue=1.6.1 # Was needed to install devtools in R` but fail and also try installing devtools though conda but failed `conda install -c conda-forge r-devtools` --> In R `install.packages("glue")` worked!

Installation within R:
```R
# Load dependencies
library("Rsamtools")
library("GenomicAlignments")

# Install ChIPseqSpikeInFree
install.packages("glue")
install.packages("devtools")
library("devtools")
install_github("stjude/ChIPseqSpikeInFree")
packageVersion('ChIPseqSpikeInFree') # 1.2.4
library("ChIPseqSpikeInFree")
```

Now run ChIPseqSpikeInFree: `conda activate ChIPseqSpikeInFree`:
```R
# Load packages
library("Rsamtools")
library("GenomicAlignments")
library("ChIPseqSpikeInFree")

# Create sample_meta.txt; tab delimited format `output/ChIPseqSpikeInFree/sample_meta.txt
metaFile <- "output/ChIPseqSpikeInFree/sample_meta.txt"
bams <- c("output/bowtie2_endtoend/2dN_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_HET_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_HET_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_HET_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_H3K27me3_R3.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_input_R3.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_input_R2.dupmark.sorted.bam")


# Run ChIPSpikeInFree
ChIPseqSpikeInFree(bamFiles = bams, chromFile = "hg38", metaFile = metaFile, prefix = "all_sample")
```
Job get terminated; let's run it within a script
```bash
conda activate ChIPseqSpikeInFree
sbatch scripts/ChIPseqSpikeInFree_all.sh # 12100044 FAIL; 12125487
```
*FAIL: `Oooops: you need to change metadata file **\n") stop("Please check whether you IDs in metaFile match with colnames(bam filenames) in parsedMatrix.`; I modifiy the `output/ChIPseqSpikeInFree/sample_meta.txt` and include only file name for ID instead of full path.*

Worked!

Let's re-run without the input and ESC_WT_H3K27me3_R3 that failed, so that reference is re-calculated maybe:

```bash
conda activate ChIPseqSpikeInFree
sbatch scripts/ChIPseqSpikeInFree_H3K27me3.sh # 12332255 ok
```

--> The SF are different

## ChIPseqSpikeInFree - Bigwig generation

Now let's produce **bigwig files as recommended by the [paper](https://github.com/stjude/ChIPseqSpikeInFree):**
1. Convert bam to bedfiles
2. Scaled the bedfile
3. Generate bigwig
**Original code:**
```bash
libSize=`cat sample1.bed|wc -l`
scale=15000000/($libSize*$SF)
genomeCoverageBed -bg -scale $scale -i sample1.bed  -g mm9.chromSizes > sample1.bedGraph
bedGraphToBigWig sample1.bedGraph mm9.chromSizes sample1.bw
```
*NOTE: Let's create a new environment to scale bed and generate bigwig with the ChIPseqSpikeInFree method:*
```bash
conda create -n BedToBigwig
conda install -c bioconda bedtools
conda install -c bioconda ucsc-bedgraphtobigwig
conda install -n ucsc openssl=1.0 # needed for bedgraphtobigwig to work
```

```bash
conda activate BedToBigwig
# Convert bam to bedfiles
sbatch scripts/bamToBed_1.sh # 12330666 ok
sbatch scripts/bamToBed_2.sh # 12330667 ok

# Apply scaling factor and transform to bigwig
sbatch scripts/BedScaledToBigwig_1.sh # 12342448; bedgraph ok bigwig transformation failed
sbatch scripts/BedScaledToBigwig_2.sh # 12342447; bedgraph ok bigwig transformation failed

# Sort the bedgraph and transform to bigwig
sbatch --dependency=afterany:12342448 scripts/SortBedToBigwig_1.sh # 12342750 ok
sbatch scripts/SortBedToBigwig_2.sh # 12342718 ok
```
*NOTE: I forget to sort the bedgraph before generate bigwig; to do next time*

--> The bigwig looks OK, except for ESC_WT_R2 that have a very different profiles.




### ChIPseqSpikeInFree - High quality reads or uniq mapp


```bash
conda activate ChIPseqSpikeInFree
sbatch scripts/ChIPseqSpikeInFree_highquality.sh # 12373851 ok
```

--> We still have the same SF tendancy. Let's now try to keep MAPQ>20 but keep only uniquely mapped reads:

Remove dupplicates and re-run ChIPseqSpikeInFree
```bash
sbatch scripts/samtools_unique_1.sh # 12378180 ok
sbatch scripts/samtools_unique_2.sh # 12378181 ok

sbatch scripts/ChIPseqSpikeInFree_unique.sh #  12385784 ok
```

--> It is much better!! So for using **ChIPseqSpikeInFree, better to use uniquely aligned reads only (remove dupplicates with PICARD)**

Generate the bigwig using these scaling factor; but generate bigwig from the `*.dupmark.sorted.bam` and NOT the `*.unique.dupmark.sorted.bam` (take as input our bam converted to bedgraph in `output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig` and output new bigwig in `output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF`)

```bash
conda activate BedToBigwig
# Apply scaling factor and transform to bigwig
sbatch scripts/BedScaledToBigwig_uniqueSF_1.sh # 12390183 ok
sbatch scripts/BedScaledToBigwig_uniqueSF_2.sh # 12390184 ok
```

--> The bigwig looks great! We observed increase signal from ESC to 2dN that we cannot observe when using the raw bigwig! Also we see decrease H3K27me3 signal in ESC vs mutants

Let's merge the bigwig into 1 file with wiggletools (will do average of bigwig signal and not sum, many options see [github](https://github.com/Ensembl/WiggleTools)):

**Installation wiggletools:**
```bash
conda activate BedToBigwig
conda install -c bioconda wiggletools
```
**Run wiggletools:**
```bash
conda activate BedToBigwig
sbatch scripts/bigwigmerge_uniqueSF.sh # 12450081 ok
sbatch scripts/bigwigmerge_uniqueSF_input.sh # FUCK THE INPUT, we do not care about them as not ChIPseqSpikeInFree norm...
```
*NOTE: bigwig are merge into 1 bedgraph which is then converted into 1 bigwig (wiggletools cannot output bigwig directly so need to pass by bedgraph or wiggle in between)*


# Coverage bigwig file
## Raw coverage bigwig

```bash
conda activate deeptools

sbatch scripts/bamtobigwig_ESC.sh # 11974332  time limit ; Rerun on scripts/bamtobigwig_canceled.sh
sbatch scripts/bamtobigwig_NPC.sh # 11974335 time limit ; Rerun on scripts/bamtobigwig_canceled.sh
sbatch scripts/bamtobigwig_2dN.sh # 11974336 ok
sbatch scripts/bamtobigwig_canceled.sh # 12019087 ok
```


## depth-norm/downsample coverage bigwig
```bash
conda activate deeptools

sbatch scripts/bamtobigwig_downsample_ESC.sh # 11980007 time limit ; 12003265 ok
sbatch scripts/bamtobigwig_downsample_ESC_input_WT.sh # 12003264
sbatch scripts/bamtobigwig_downsample_NPC.sh # 11980008 time limit ; 12003264 ok
sbatch scripts/bamtobigwig_downsample_NPC_input_WT.sh # 12003265
sbatch scripts/bamtobigwig_downsample_2dN.sh # 11980009 ok
```

## input-normalize coverage bigwig
bamCompare normalize per sequence depth and then perform calculation

Light test with 1 file to see how it perform
```bash
conda activate deeptools
bamCompare -b1 output/bowtie2_endtoend/ESC_WT_H3K27me3_R1.dupmark.sorted.bam -b2 output/bowtie2_endtoend/ESC_WT_input_R1.dupmark.sorted.bam -o output/bigwig_inputNorm/ESC_WT_R1_log2ratio.bw
```
Works great! ratio and not log2ratio is better (= ratio of read counts per bin in the IP sample relative to the input sample). So run all samples:

```bash
conda activate deeptools
sbatch scripts/bamtobigwig_inputNorm_ESC.sh # 12141577 ok
sbatch scripts/bamtobigwig_inputNorm_NPC.sh # 12141576 ok
sbatch scripts/bamtobigwig_inputNorm_2dN.sh # 12141578 ok
```



## ChIPseqSpikeInFree coverage bigwig
Follow recommendation from [github](https://github.com/stjude/ChIPseqSpikeInFree):
```bash
libSize=`cat sample1.bed|wc -l`
scale=15000000/($libSize*$SF)
genomeCoverageBed -bg -scale $scale -i sample1.bed  -g mm9.chromSizes > sample1.bedGraph
bedGraphToBigWig sample1.bedGraph mm9.chromSizes sample1.bw
```
We can adapt this to work with bam and bamCoverage to make it more straightforward:
```bash
conda activate deeptools

# Command example:
libSize=`samtools view -c -F 260 sample1.bam`
SF=<your_scaling_factor>
scale=$(echo "15000000/($libSize*$SF)" | bc -l)
bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor $scale

# All sample together
sbatch scripts/bamtobigwig_ChIPseqSpikeInFree.sh # 12141584; cancel time limit (only 2dN/ESC have been processed...)
sbatch scripts/bamtobigwig_ChIPseqSpikeInFree_canceled.sh # 12146384 ok

```
*NOTE: 15000000 is a reference to normalize the read counts in the ChIP-seq data. It represents a target library size to which the actual library size will be scaled. Could have choose any number, but 15m is commonly used*





## input-normalize ChIPseqSpikeInFree coverage bigwig

XXX

Can apply ChIPseqSpikeInFree scaling factor and then normalize


# Peak calling

## MACS2 peak calling raw
```bash
conda activate macs2
# example for 1 file
macs2 callpeak -t output/bowtie2_endtoend/ESC_HET_H3K27me3_R1.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/ESC_HET_input_R1.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir output/macs2 -n ESC_HET_H3K27me3_R1 --broad

# run per time
sbatch scripts/macs2_ESC.sh # 11979368 ok
sbatch scripts/macs2_NPC.sh # 11979372 ok
sbatch scripts/macs2_2dN.sh # 11979373 ok
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
sbatch scripts/macs2_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive
```

--> Overall the 0.005/0.001 qval is more accurate to call peak

## MACS2 peak calling depth-norm/downsample
```bash
conda activate macs2

# run per time
sbatch scripts/macs2_downsample_ESC.sh # 11980099 ok
sbatch scripts/macs2_downsample_NPC.sh # 11980109 ok
sbatch scripts/macs2_downsample_2dN.sh # 11980110 ok
```

```bash
conda activate bowtie2 # for bedtools
sbatch scripts/macs2_downsample_peak_signif.sh # 1.30103/2/2.30103/3/4/5 # Run in interactive
```

--> Overall the non downsample method show more peaks at the same qvalue 



# PCA on Bigwig files

Let's do PCA with [multiBigwigSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html) and [PCAplot](https://deeptools.readthedocs.io/en/2.4.1/content/tools/plotPCA.html)/[plotCorrelation](https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html) from deeptools:

- Use bin mode (score of the compile bigwig is calculated at a 10kb bins (default), on the entire genome)
- Let's do it all samples together and per genotype for better vizualization

## PCA on raw files

```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all_H3K27me3.sh # 12003437 ok
sbatch scripts/multiBigwigSummary_all.sh # 12003478 FAIL; some input files missing (ESC_WT_input_R1.dupmark.sorted.bw ESC_WT_input_R2.dupmark.sorted.bw ESC_WT_input_R3.dupmark.sorted.bw NPC_WT_input_R2.dupmark.sorted.bw) --> Re-generated and job re-launch 12024171 ok
## Time per time (for genotype effect)
sbatch scripts/multiBigwigSummary_ESC.sh # 12037427 ok
sbatch scripts/multiBigwigSummary_NPC.sh # 12037490 ok
sbatch scripts/multiBigwigSummary_2dN.sh # 12037421 ok

sbatch scripts/multiBigwigSummary_ESC_noinput.sh # 12046629 ok
sbatch scripts/multiBigwigSummary_NPC_noinput.sh # 12046630 ok
sbatch scripts/multiBigwigSummary_2dN_noinput.sh # 12046584 ok

## Genotype per genotype (for time effect)
sbatch scripts/multiBigwigSummary_WT.sh # 12038201 ok
sbatch scripts/multiBigwigSummary_HET.sh # 12038200 ok
sbatch scripts/multiBigwigSummary_KO.sh # 12038199 ok

sbatch scripts/multiBigwigSummary_WT_noinput.sh # 12046734 ok
sbatch scripts/multiBigwigSummary_HET_noinput.sh # 12046735 ok
sbatch scripts/multiBigwigSummary_KO_noinput.sh # 12046733 ok


# Plot
## All genotypes all points
plotPCA -in output/bigwig/multiBigwigSummary_all.npz \
    --transpose \
    --ntop 0 \
    --labels 2dN_HET_input_R1 2dN_HET_input_R2 2dN_KO_input_R1 2dN_KO_input_R2 2dN_WT_input_R1 2dN_WT_input_R2 ESC_HET_input_R1 ESC_HET_input_R2 ESC_KO_input_R1 ESC_KO_input_R2 ESC_WT_input_R1 ESC_WT_input_R2 ESC_WT_input_R3 NPC_HET_input_R1 NPC_HET_input_R2 NPC_KO_input_R1 NPC_KO_input_R2 NPC_WT_input_R1 NPC_WT_input_R2 2dN_HET_H3K27me3_R1 2dN_HET_H3K27me3_R2 2dN_KO_H3K27me3_R1 2dN_KO_H3K27me3_R2 2dN_WT_H3K27me3_R1 2dN_WT_H3K27me3_R2 ESC_HET_H3K27me3_R1 ESC_HET_H3K27me3_R2 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 NPC_HET_H3K27me3_R1 NPC_HET_H3K27me3_R2 NPC_KO_H3K27me3_R1 NPC_KO_H3K27me3_R2 NPC_WT_H3K27me3_R1 NPC_WT_H3K27me3_R2 \
    -o output/bigwig/multiBigwigSummary_all_plotPCA.pdf
plotPCA -in output/bigwig/multiBigwigSummary_all_H3K27me3.npz \
    --transpose \
    --ntop 0 \
    --labels 2dN_HET_H3K27me3_R1 2dN_HET_H3K27me3_R2 2dN_KO_H3K27me3_R1 2dN_KO_H3K27me3_R2 2dN_WT_H3K27me3_R1 2dN_WT_H3K27me3_R2 ESC_HET_H3K27me3_R1 ESC_HET_H3K27me3_R2 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 NPC_HET_H3K27me3_R1 NPC_HET_H3K27me3_R2 NPC_KO_H3K27me3_R1 NPC_KO_H3K27me3_R2 NPC_WT_H3K27me3_R1 NPC_WT_H3K27me3_R2 \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_plotPCA.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_heatmap.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_all_H3K27me3.npz \
    --corMethod pearson \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_all_H3K27me3_heatmap.pdf

## Time per time (for genotype effect)
plotPCA -in output/bigwig/multiBigwigSummary_ESC_noinput.npz \
    --transpose \
    --ntop 0 \
    --labels ESC_HET_H3K27me3_R1 ESC_HET_H3K27me3_R2 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 \
    -o output/bigwig/multiBigwigSummary_ESC_noinput_plotPCA.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_ESC_noinput.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_ESC_noinput_heatmap.pdf

plotPCA -in output/bigwig/multiBigwigSummary_NPC_noinput.npz \
    --transpose \
    --ntop 0 \
    --labels NPC_HET_H3K27me3_R1 NPC_HET_H3K27me3_R2 NPC_KO_H3K27me3_R1 NPC_KO_H3K27me3_R2 NPC_WT_H3K27me3_R1 NPC_WT_H3K27me3_R2 \
    -o output/bigwig/multiBigwigSummary_NPC_noinput_plotPCA.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_NPC_noinput.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_NPC_noinput_heatmap.pdf

plotPCA -in output/bigwig/multiBigwigSummary_2dN_noinput.npz \
    --transpose \
    --ntop 0 \
    --labels 2dN_HET_H3K27me3_R1 2dN_HET_H3K27me3_R2 2dN_KO_H3K27me3_R1 2dN_KO_H3K27me3_R2 2dN_WT_H3K27me3_R1 2dN_WT_H3K27me3_R2 \
    -o output/bigwig/multiBigwigSummary_2dN_noinput_plotPCA.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_2dN_noinput.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_2dN_noinput_heatmap.pdf
## Genotype per genotype (for time effect)
plotPCA -in output/bigwig/multiBigwigSummary_WT_noinput.npz \
    --transpose \
    --ntop 0 \
    --labels 2dN_WT_H3K27me3_R1 2dN_WT_H3K27me3_R2 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 NPC_WT_H3K27me3_R1 NPC_WT_H3K27me3_R2 \
    -o output/bigwig/multiBigwigSummary_WT_noinput_plotPCA.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_WT_noinput.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_WT_noinput_heatmap.pdf

plotPCA -in output/bigwig/multiBigwigSummary_HET_noinput.npz \
    --transpose \
    --ntop 0 \
    --labels 2dN_HET_H3K27me3_R1 2dN_HET_H3K27me3_R2 ESC_HET_H3K27me3_R1 ESC_HET_H3K27me3_R2 NPC_HET_H3K27me3_R1 NPC_HET_H3K27me3_R2 \
    -o output/bigwig/multiBigwigSummary_HET_noinput_plotPCA.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_HET_noinput.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_HET_noinput_heatmap.pdf

plotPCA -in output/bigwig/multiBigwigSummary_KO_noinput.npz \
    --transpose \
    --ntop 0 \
    --labels 2dN_KO_H3K27me3_R1 2dN_KO_H3K27me3_R2 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 NPC_KO_H3K27me3_R1 NPC_KO_H3K27me3_R2 \
    -o output/bigwig/multiBigwigSummary_KO_noinput_plotPCA.pdf
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_KO_noinput.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig/multiBigwigSummary_KO_noinput_heatmap.pdf
```


## PCA on downsample files
```bash
conda activate deeptools
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_all_H3K27me3_downsample.sh # 12019080 ok
sbatch scripts/multiBigwigSummary_all_downsample.sh # 12019081 ok
# Time per time (for genotype effect)

# Genotype per genotype (for time effect)

# Plot
## All genotypes all points
plotPCA -in output/bigwig_downsample/multiBigwigSummary_all_downsample.npz \
    --transpose \
    --ntop 0 \
    --labels 2dN_HET_input_R1 2dN_HET_input_R2 2dN_KO_input_R1 2dN_KO_input_R2 2dN_WT_input_R1 2dN_WT_input_R2 ESC_HET_input_R1 ESC_HET_input_R2 ESC_KO_input_R1 ESC_KO_input_R2 ESC_WT_input_R1 ESC_WT_input_R2 ESC_WT_input_R3 NPC_HET_input_R1 NPC_HET_input_R2 NPC_KO_input_R1 NPC_KO_input_R2 NPC_WT_input_R1 NPC_WT_input_R2 2dN_HET_H3K27me3_R1 2dN_HET_H3K27me3_R2 2dN_KO_H3K27me3_R1 2dN_KO_H3K27me3_R2 2dN_WT_H3K27me3_R1 2dN_WT_H3K27me3_R2 ESC_HET_H3K27me3_R1 ESC_HET_H3K27me3_R2 ESC_KO_H3K27me3_R1 ESC_KO_H3K27me3_R2 ESC_WT_H3K27me3_R1 ESC_WT_H3K27me3_R2 ESC_WT_H3K27me3_R3 NPC_HET_H3K27me3_R1 NPC_HET_H3K27me3_R2 NPC_KO_H3K27me3_R1 NPC_KO_H3K27me3_R2 NPC_WT_H3K27me3_R1 NPC_WT_H3K27me3_R2 \
    -o output/bigwig_downsample/multiBigwigSummary_all_downsample_plotPCA.pdf
plotCorrelation \
    -in output/bigwig_downsample/multiBigwigSummary_all_downsample.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig_downsample/multiBigwigSummary_all_downsample_heatmap.pdf
plotCorrelation \
    -in output/bigwig_downsample/multibigwigSummary_all_H3K27me3_downsample.npz \
    --corMethod pearson \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o output/bigwig_downsample/multiBigwigSummary_all_H3K27me3_downsample_heatmap.pdf
## Time per time (for genotype effect)

## Genotype per genotype (for time effect)
```
*NOTE: I tested Spearman and Pearson correlation. Pearson perform better (more accurate); I also remove outlier as Pearson are more sensitive; could come from blacklist regions.*

--> Let's **use the non-downsample/raw** as it give more peaks; and the replicate looks more 'similar' on IGV. Moreover MACS2 already account for these differences while calling peak.




# DiffBind

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
sample_count = dba.count(sample_dba, bParallel=TRUE) # bParallel=TRUE is to count using multiple processors

```
The job is terminated, probably because there is too many samples! 
Let's try to run it in a Rscript with shit tons of memory (500G)


```bash
conda activate DiffBind
sbatch scripts/DiffBind_ChIP.sh # 122003317 ok
```
It worked!

```R
# Load the raw counts
load("output/DiffBind/sample_count.RData")


# plot
pdf("output/DiffBind/clustering_sample.pdf", width=14, height=20)  
plot(sample_count)
dev.off()

pdf("output/DiffBind/PCA_sample.pdf", width=14, height=20) 
dba.plotPCA(sample_count,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()

# Blacklist is already applied, so let's generate GreyList (IGG)
## Greylist generation
sample_dba_greylist = dba.blacklist(sample_count, blacklist=FALSE, greylist=TRUE, cores=1)
```
*NOTE: I added cores=1 to avoid a weird error: `sample_dba_greylist = dba.blacklist(sample_count, blacklist=FALSE, greylist=TRUE) Genome detected: Hsapiens.UCSC.hg38 Counting control reads for greylist... Error in value[[3L]](cond) :    GreyListChIP error: Error in result[[njob]] <- value: attempt to select less than one element in OneIndex In addition: Warning message: In parallel::mccollect(wait = FALSE, timeout = 1) :   1 parallel job did not deliver a result > sample_dba_greylist Error: object 'sample_dba_greylist' not found`*

The greylist has been killed again... Let's run it within a script:

```bash
conda activate DiffBind
sbatch scripts/DiffBind_ChIP_greylist.sh # 12270854 ok
```

It worked!


```R
# Load the greylist filtered counts
load("output/DiffBind/sample_dba_greylist.RData")

# Now check how clustering look
sample_count_blackgreylist = dba.count(sample_dba_greylist)


# plot
pdf("output/DiffBind/clustering_blackgreylist.pdf", width=14, height=20)
plot(sample_count_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()



# Test different normalization and pick the one that identify most diff. bound sites
## default lib-depth normalization
sample_count_blackgreylist_lib = dba.normalize(sample_count_blackgreylist)

# plot
pdf("output/DiffBind/clustering_blackgreylist_lib.pdf", width=14, height=20)
plot(sample_count_blackgreylist_lib)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_lib.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_lib,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


##  RLE-depth normalization
sample_count_blackgreylist_RLE = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_RLE)

# plot
pdf("output/DiffBind/clustering_blackgreylist_RLE.pdf", width=14, height=20)
plot(sample_count_blackgreylist_RLE)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_RLE,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


##  TMM-depth normalization
sample_count_blackgreylist_TMM = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_TMM)

# plot
pdf("output/DiffBind/clustering_blackgreylist_TMM.pdf", width=14, height=20)
plot(sample_count_blackgreylist_TMM)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_TMM,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


##  TMM-depth normalization; library PEARKREADS
sample_count_blackgreylist_TMM_PEAKREADS = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_TMM, library=DBA_LIBSIZE_PEAKREADS)

# plot
pdf("output/DiffBind/clustering_blackgreylist_TMM_PEAKREADS.pdf", width=14, height=20)
plot(sample_count_blackgreylist_TMM_PEAKREADS)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_TMM_PEAKREADS.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_TMM_PEAKREADS,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()

##  TMM-depth normalization; library PEARKREADS
sample_count_blackgreylist_TMM_FULL = dba.normalize(sample_count_blackgreylist, normalize=DBA_NORM_TMM, library=DBA_LIBSIZE_FULL)

# plot
pdf("output/DiffBind/clustering_blackgreylist_TMM_FULL.pdf", width=14, height=20)
plot(sample_count_blackgreylist_TMM_FULL)
dev.off()

pdf("output/DiffBind/PCA_blackgreylist_TMM_FULL.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_TMM_FULL,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


# Let's pick TMM for the next (as also work for Cut Run) and the Greylist one:
### Greylist matrix contrast
sample_count_blackgreylist
# Set up contrast for comparison (all will be compare to the WT)
sample_count_blackgreylist_contrast = dba.contrast(sample_count_blackgreylist, reorderMeta = list(Treatment="WT", Condition ="ESC"), categories=c(DBA_TREATMENT,DBA_CONDITION))

sample_count_blackgreylist_contrast_analyze = dba.analyze(sample_count_blackgreylist_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)


### TMM matrix contrast
sample_count_blackgreylist_TMM
# Set up contrast for comparison (all will be compare to the WT)
sample_count_blackgreylist_TMM_contrast = dba.contrast(sample_count_blackgreylist_TMM, reorderMeta = list(Treatment="WT", Condition ="ESC"), categories=c(DBA_TREATMENT,DBA_CONDITION))

sample_count_blackgreylist_TMM_contrast_analyze = dba.analyze(sample_count_blackgreylist_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE)

### --> Pick TMM DESEQ2
# Set up contrast for comparison (all will be compare to the WT)
sample_count_blackgreylist_TMM_contrast = dba.contrast(sample_count_blackgreylist_TMM, reorderMeta = list(Treatment="WT", Condition ="ESC"), categories=c(DBA_TREATMENT,DBA_CONDITION), minMembers=2)

sample_count_blackgreylist_TMM_contrast_analyze = dba.analyze(sample_count_blackgreylist_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE)

dba.contrast(sample_count_blackgreylist_TMM, minMembers=2)


```
*NOTE: I tested playing with the library parameter for normalization (Full,DBA_LIBSIZE_PEAKREADS,etc) ant it do not change anythings*

--> RLE and TMM normalization perform well (Lib depth poorly separate our data)

The multi-factor analyses is not clear, let's generate count matrix at each time point, and from this identify the diff bound sites (*For ESC; replicate 3 has not been taken*)



Let's apply the ChIPseqSpikeInFree scaling factor to normalize data (normalize library size, then normalize per seq. depth; as for the CutRun). ChIPseqSpikeInFree-norm-library-size = library-size * SF. `samtools flagstat output/bowtie2_endtoend/*.dupmark.sorted.bam` used to obtain library size (first value=library size):
- sample / library size * SF = scaled library size
- 2dN_HET_H3K27me3_R1 / 70887858 1.97 = 139649080
- 2dN_HET_H3K27me3_R2 / 67346082 1.75 = 117855644
- 2dN_KO_H3K27me3_R1 / 58363794 1.46 = 85211139
- 2dN_KO_H3K27me3_R2 / 71964648 1 = 71964648
- 2dN_WT_H3K27me3_R1 / 71682964 1.29 = 92471024
- 2dN_WT_H3K27me3_R2 / 60792560 1.69 = 102739426
- ESC_HET_H3K27me3_R1 / 85933694 10.51 = 903163124
- ESC_HET_H3K27me3_R2 / 66583922 23.35 = 1554734579
- ESC_KO_H3K27me3_R1 / 86014942 10.06 = 865310316
- ESC_KO_H3K27me3_R2 / 57643920 15.78 = 909621058
- ESC_WT_H3K27me3_R1 / 90968406 7 = 636778842
- ESC_WT_H3K27me3_R2 / 79650052 4.31 = 343291724
- NPC_HET_H3K27me3_R1 / 40423510 1.13 = 45678566
- NPC_HET_H3K27me3_R2 / 82710030 1.43 = 118275343
- NPC_KO_H3K27me3_R1 / 67904376 1.55 = 105251783
- NPC_KO_H3K27me3_R2 / 84619268 2.64 = 223394867
- NPC_WT_H3K27me3_R1 / 77698696 1.45 = 112663109
- NPC_WT_H3K27me3_R2 / 72126718 1.51 = 108911344

*NOTE: I run samtools flagstat manually and copy/paste, I run a sbatch job to keep track also as:*
 ```bash
 sbatch scripts/libsize_dupmark.sh # 12397755 ok 
 ```
	
Let's generate different meta_sample files:
- meta_sample_all = All samples (except ESC WT R3)
- meta_sample_ESC = All ESC samples (WT, HET, KO)
- meta_sample_NPC
- meta_sample_2dN
- meta_sample_WT
- meta_sample_HET
- meta_sample_KO

Generate/save **count matrix for each class**:
```bash
conda activate DiffBind

sbatch scripts/DiffBind_all.sh # 12398837 ok
sbatch scripts/DiffBind_ESC.sh # 12398844 ok
sbatch scripts/DiffBind_NPC.sh # 12398846 ok
sbatch scripts/DiffBind_2dN.sh # 12398848 ok
sbatch scripts/DiffBind_WT.sh # 12398852 ok
sbatch scripts/DiffBind_HET.sh # 12398854 ok
sbatch scripts/DiffBind_KO.sh # 12398856 ok
```

Apply/save **GreyList count matrix for each class**:
```bash
conda activate DiffBind

sbatch scripts/DiffBind_greylist.sh # 12436659 ok
```
*NOTE: For all analyses less than 12 peaks have been removed with the greylist filtering*

Now let's make plot (greylist and ChIPseqSpikeInFree (CSIF) corrected) and diff analyses for:
- **genotype effect** per time-point:
    - ESC: WT vs KO, WT vs HET,…
    - NPC: WT vs KO, WT vs HET,…
    - 2dN: WT vs KO, WT vs HET,…
- **time course effect** per genotype:
    - WT: ESC vs NPC, ESC, 2dN,…
    - HET: ESC vs NPC, ESC, 2dN,…
    - KO: ESC vs NPC, ESC, 2dN,…
- **genotype effect over the time course**



```R
library("DiffBind")

# Genotype effect ~ ESC state
## Load the greylist
load("output/DiffBind/sample_count_ESC_greylist.RData")

sample_count_ESC_greylist

## plot greylist
pdf("output/DiffBind/clustering_ESC_greylist.pdf", width=14, height=20)
plot(sample_count_ESC_greylist)
dev.off()

pdf("output/DiffBind/PCA_ESC_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_greylist,DBA_REPLICATE, label=c(DBA_TREATMENT))
dev.off()

## Lib-size ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; edgeR: 35/399/362; DESEQ2: 45/4956/4480; non-lib norm= 2/2/2 and 39/440/417)

sample_count_ESC_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_ESC_greylist, library = c(903163124, 1554734579, 865310316, 909621058, 636778842, 343291724), normalize = DBA_NORM_LIB) # Default


pdf("output/DiffBind/clustering_ESC_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20)
plot(sample_count_ESC_greylist_LibCSIFScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_ESC_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_greylist_LibCSIFScaled_LIB,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_ESC_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


## TMM ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; edgeR: 35/399/362 ; DESEQ2: 82/737/787)

sample_count_ESC_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_ESC_greylist, library = c(903163124, 1554734579, 865310316, 909621058, 636778842, 343291724), normalize = DBA_NORM_TMM) # Default


pdf("output/DiffBind/clustering_ESC_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_ESC_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_ESC_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_greylist_LibCSIFScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_ESC_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


## RLE ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; edgeR: 35/399/362 ; DESEQ2: 84/719/770)

sample_count_ESC_greylist_LibCSIFScaled_RLE = dba.normalize(sample_count_ESC_greylist, library = c(903163124, 1554734579, 865310316, 909621058, 636778842, 343291724), normalize = DBA_NORM_RLE) # Default


pdf("output/DiffBind/clustering_ESC_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20)
plot(sample_count_ESC_greylist_LibCSIFScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_ESC_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_greylist_LibCSIFScaled_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_ESC_greylist_LibCSIFScaled_RLE_contrast = dba.contrast(sample_count_ESC_greylist_LibCSIFScaled_RLE, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_ESC_greylist_LibCSIFScaled_RLE_contrast_analyze = dba.analyze(sample_count_ESC_greylist_LibCSIFScaled_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)

### Diff analyses with LIB/DESEQ2:


sample_count_ESC_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_ESC_greylist, library = c(903163124, 1554734579, 865310316, 909621058, 636778842, 343291724), normalize = DBA_NORM_LIB) # Default

sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_ESC_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_ESC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_ESC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_ESC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_ESC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_ESC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_ESC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_ESC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
plot(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_ESC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
plot(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_ESC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
plot(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_ESC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_ESC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_ESC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_ESC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_ESC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_ESC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_ESC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_ESC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_ESC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()



## Export the Diff Bind regions
### Convert to GR object
sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast1 <- dba.report(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast2 <- dba.report(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast3 <- dba.report(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast1_df <- data.frame(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast1)
sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast2_df <- data.frame(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast2)
sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast3_df <- data.frame(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast3)

colnames(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast1_df, file="output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast2_df, file="output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast3_df, file="output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_LIB_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)




# Genotype effect ~ NPC state
## Load the greylist
load("output/DiffBind/sample_count_NPC_greylist.RData")

sample_count_NPC_greylist

## plot greylist
pdf("output/DiffBind/clustering_NPC_greylist.pdf", width=14, height=20)
plot(sample_count_NPC_greylist)
dev.off()

pdf("output/DiffBind/PCA_NPC_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_greylist,DBA_REPLICATE, label=c(DBA_TREATMENT))
dev.off()

## Lib-size ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; edgeR: 25/9/75; DESEQ2: 20/18/109; non-lib norm= 25/9/75 and 46/19/94)

sample_count_NPC_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_NPC_greylist, library = c(45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_LIB) # Default


pdf("output/DiffBind/clustering_NPC_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20)
plot(sample_count_NPC_greylist_LibCSIFScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_NPC_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_greylist_LibCSIFScaled_LIB,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_NPC_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


## TMM ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; edgeR: 25/9/75; DESEQ2: 27/31/148)

sample_count_NPC_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_NPC_greylist, library = c(45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_TMM) # Default


pdf("output/DiffBind/clustering_NPC_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_NPC_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_NPC_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_greylist_LibCSIFScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_NPC_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



## RLE ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; edgeR: 25/9/75; DESEQ2: 25/31/149)

sample_count_NPC_greylist_LibCSIFScaled_RLE = dba.normalize(sample_count_NPC_greylist, library = c(45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_RLE) # Default


pdf("output/DiffBind/clustering_NPC_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20)
plot(sample_count_NPC_greylist_LibCSIFScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_NPC_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_greylist_LibCSIFScaled_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_NPC_greylist_LibCSIFScaled_RLE_contrast = dba.contrast(sample_count_NPC_greylist_LibCSIFScaled_RLE, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_NPC_greylist_LibCSIFScaled_RLE_contrast_analyze = dba.analyze(sample_count_NPC_greylist_LibCSIFScaled_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



### Diff analyses with LIB/DESEQ2:


sample_count_NPC_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_NPC_greylist, library = c(45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_LIB) # Default

sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_NPC_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_NPC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_NPC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_NPC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_NPC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_NPC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_NPC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_NPC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
plot(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_NPC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
plot(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_NPC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
plot(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_NPC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_NPC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_NPC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_NPC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_NPC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_NPC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_NPC_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_NPC_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_NPC_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()



## Export the Diff Bind regions
### Convert to GR object
sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast1 <- dba.report(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast2 <- dba.report(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast3 <- dba.report(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast1_df <- data.frame(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast1)
sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast2_df <- data.frame(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast2)
sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast3_df <- data.frame(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast3)

colnames(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast1_df, file="output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast2_df, file="output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast3_df, file="output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_LIB_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



# Genotype effect ~ 2dN state
## Load the greylist
load("output/DiffBind/sample_count_2dN_greylist.RData")

sample_count_2dN_greylist

## plot greylist
pdf("output/DiffBind/clustering_2dN_greylist.pdf", width=14, height=20)
plot(sample_count_2dN_greylist)
dev.off()

pdf("output/DiffBind/PCA_2dN_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_greylist,DBA_REPLICATE, label=c(DBA_TREATMENT))
dev.off()

## Lib-size ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; edgeR: 15/20/1; DESEQ2: 2/14/15; non-lib norm= 25/9/75 and 117/106/0)

sample_count_2dN_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_2dN_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426), normalize = DBA_NORM_LIB) # Default


pdf("output/DiffBind/clustering_2dN_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20)
plot(sample_count_2dN_greylist_LibCSIFScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_2dN_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_greylist_LibCSIFScaled_LIB,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_2dN_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



## TMM ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; DESEQ2: 14/65/2)

sample_count_2dN_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_2dN_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426), normalize = DBA_NORM_TMM) # Default


pdf("output/DiffBind/clustering_2dN_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_2dN_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_2dN_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_greylist_LibCSIFScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_2dN_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



## RLE ChIPseqSpikeInFree norm (HETvsKO/HETvsWT/KOvsWT; DESEQ2: 11/56/1)

sample_count_2dN_greylist_LibCSIFScaled_RLE = dba.normalize(sample_count_2dN_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426), normalize = DBA_NORM_RLE) # Default


pdf("output/DiffBind/clustering_2dN_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20)
plot(sample_count_2dN_greylist_LibCSIFScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_2dN_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_greylist_LibCSIFScaled_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_2dN_greylist_LibCSIFScaled_RLE_contrast = dba.contrast(sample_count_2dN_greylist_LibCSIFScaled_RLE, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_2dN_greylist_LibCSIFScaled_RLE_contrast_analyze = dba.analyze(sample_count_2dN_greylist_LibCSIFScaled_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



### Diff analyses with LIB/DESEQ2:


sample_count_2dN_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_2dN_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426), normalize = DBA_NORM_LIB) # Default

sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_2dN_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_2dN_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_2dN_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_2dN_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_2dN_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_2dN_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_2dN_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_2dN_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
plot(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_2dN_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
plot(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_2dN_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
plot(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_2dN_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_2dN_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_2dN_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_2dN_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_2dN_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_2dN_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_2dN_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_2dN_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_2dN_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()



## Export the Diff Bind regions
### Convert to GR object
sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast1 <- dba.report(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast2 <- dba.report(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast3 <- dba.report(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast1_df <- data.frame(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast1)
sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast2_df <- data.frame(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast2)
sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast3_df <- data.frame(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast3)

colnames(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast1_df, file="output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast2_df, file="output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast3_df, file="output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_LIB_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

**genotype effect** per time-point:
    - ESC: WT vs KO, WT vs HET,… --> LIB/DESEQ2 norm perform best (DESEQ2/LIB n = 45/4956/4480)
    - NPC: WT vs KO, WT vs HET,… --> Very few diff sites! (TMM perform best, but no drastic diff. with the other methods...; so for uniformity I used LIB norm that perform very well for ESC) (DESEQ2/LIB n = 20/18/109)
    - 2dN: WT vs KO, WT vs HET,… --> Non CSIF give more diff sites, but 0 for KOvsWT, weird, let's keep CSIF-norm LIB-DESEQ2 (DESEQ2/LIB n = 2/14/5)

--> Only differences between genotypes at ESC (do not follow biological expectation; ie. KO and HET decrease H3K27me3 vs WT)

--> At NPC and 2dN very few diff bound sites; in agreement with the CutRun at 8wN

Now lets look **Time effect per genotype:**

```R

# Time effect ~ WT genotype
## Load the greylist
load("output/DiffBind/sample_count_WT_greylist.RData")

sample_count_WT_greylist

## plot greylist
pdf("output/DiffBind/clustering_WT_greylist.pdf", width=14, height=20)
plot(sample_count_WT_greylist)
dev.off()

pdf("output/DiffBind/PCA_WT_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_greylist,DBA_REPLICATE, label=c(DBA_CONDITION))
dev.off()

## Lib-size ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; edgeR: 1888/2/1690; DESEQ2: 10016/0/10030; non-lib norm= 1717/0/1594)

sample_count_WT_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_WT_greylist, library = c(92471024,102739426,636778842,343291724,112663109,108911344), normalize = DBA_NORM_LIB) # Default


pdf("output/DiffBind/clustering_WT_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20)
plot(sample_count_WT_greylist_LibCSIFScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_WT_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_greylist_LibCSIFScaled_LIB,DBA_REPLICATE, label=DBA_CONDITION)
dev.off()

sample_count_WT_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_WT_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)




## TMM ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; edgeR: 25/9/75; DESEQ2: 1702/1/1549)

sample_count_WT_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_WT_greylist, library = c(92471024,102739426,636778842,343291724,112663109,108911344), normalize = DBA_NORM_TMM) # Default


pdf("output/DiffBind/clustering_WT_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_WT_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_WT_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_greylist_LibCSIFScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_WT_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_WT_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



## RLE ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; DESEQ2: 1688/1/1532)

sample_count_WT_greylist_LibCSIFScaled_RLE = dba.normalize(sample_count_WT_greylist, library = c(92471024,102739426,636778842,343291724,112663109,108911344), normalize = DBA_NORM_RLE) # Default


pdf("output/DiffBind/clustering_WT_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20)
plot(sample_count_WT_greylist_LibCSIFScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_WT_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_greylist_LibCSIFScaled_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_WT_greylist_LibCSIFScaled_RLE_contrast = dba.contrast(sample_count_WT_greylist_LibCSIFScaled_RLE, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_WT_greylist_LibCSIFScaled_RLE_contrast_analyze = dba.analyze(sample_count_WT_greylist_LibCSIFScaled_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



### Diff analyses with LIB/DESEQ2:


sample_count_WT_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_WT_greylist, library = c(92471024,102739426,636778842,343291724,112663109,108911344), normalize = DBA_NORM_LIB) # Default

sample_count_WT_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_WT_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_WT_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_WT_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_WT_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_WT_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_WT_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_WT_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_WT_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
plot(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_WT_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
plot(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_WT_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
plot(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_WT_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_WT_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_WT_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_WT_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_WT_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_WT_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_WT_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_WT_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_WT_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()



## Export the Diff Bind regions
### Convert to GR object
sample_count_WT_greylist_LibCSIFScaled_LIB_contrast1 <- dba.report(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_WT_greylist_LibCSIFScaled_LIB_contrast2 <- dba.report(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_WT_greylist_LibCSIFScaled_LIB_contrast3 <- dba.report(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_WT_greylist_LibCSIFScaled_LIB_contrast1_df <- data.frame(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast1)
sample_count_WT_greylist_LibCSIFScaled_LIB_contrast2_df <- data.frame(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast2)
sample_count_WT_greylist_LibCSIFScaled_LIB_contrast3_df <- data.frame(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast3)

colnames(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast1_df, file="output/DiffBind/sample_count_WT_greylist_LibCSIFScaled_LIB_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast2_df, file="output/DiffBind/sample_count_WT_greylist_LibCSIFScaled_LIB_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast3_df, file="output/DiffBind/sample_count_WT_greylist_LibCSIFScaled_LIB_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)





##### some non-norm plots:

### Diff analyses with LIB/DESEQ2:


sample_count_WT_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_WT_greylist, normalize = DBA_NORM_LIB) # Default

sample_count_WT_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_WT_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



# Examining results
## volcano plot with diff sites
pdf("output/DiffBind/volcano_WT_greylist_LIB_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()



## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_WT_greylist_LIB_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_WT_greylist_LIB_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_WT_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()



# Time effect ~ HET genotype
## Load the greylist
load("output/DiffBind/sample_count_HET_greylist.RData")

sample_count_HET_greylist

## plot greylist
pdf("output/DiffBind/clustering_HET_greylist.pdf", width=14, height=20)
plot(sample_count_HET_greylist)
dev.off()

pdf("output/DiffBind/PCA_HET_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_greylist,DBA_REPLICATE, label=c(DBA_CONDITION))
dev.off()

## Lib-size ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; edgeR: 2711/3/2948; DESEQ2: 6354/2/6474; non-lib norm= 2875/3/3015)

sample_count_HET_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_HET_greylist, library = c(139649080,117855644,903163124,1554734579,45678566,118275343), normalize = DBA_NORM_LIB) # Default


pdf("output/DiffBind/clustering_HET_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20)
plot(sample_count_HET_greylist_LibCSIFScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_HET_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_greylist_LibCSIFScaled_LIB,DBA_REPLICATE, label=DBA_CONDITION)
dev.off()


sample_count_HET_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_HET_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)




## TMM ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; DESEQ2: 3846/16/3716)

sample_count_HET_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_HET_greylist, library = c(139649080,117855644,903163124,1554734579,45678566,118275343), normalize = DBA_NORM_TMM) # Default


pdf("output/DiffBind/clustering_HET_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_HET_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_HET_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_greylist_LibCSIFScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_HET_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_HET_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



## RLE ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; DESEQ2: 3848/27/3720)

sample_count_HET_greylist_LibCSIFScaled_RLE = dba.normalize(sample_count_HET_greylist, library = c(139649080,117855644,903163124,1554734579,45678566,118275343), normalize = DBA_NORM_RLE) # Default


pdf("output/DiffBind/clustering_HET_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20)
plot(sample_count_HET_greylist_LibCSIFScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_HET_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_greylist_LibCSIFScaled_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_HET_greylist_LibCSIFScaled_RLE_contrast = dba.contrast(sample_count_HET_greylist_LibCSIFScaled_RLE, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_HET_greylist_LibCSIFScaled_RLE_contrast_analyze = dba.analyze(sample_count_HET_greylist_LibCSIFScaled_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



### Diff analyses with LIB/DESEQ2:


sample_count_HET_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_HET_greylist, library = c(139649080,117855644,903163124,1554734579,45678566,118275343), normalize = DBA_NORM_LIB) # Default

sample_count_HET_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_HET_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_HET_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_HET_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_HET_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_HET_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_HET_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_HET_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_HET_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
plot(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_HET_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
plot(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_HET_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
plot(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_HET_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_HET_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_HET_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_HET_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_HET_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_HET_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_HET_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_HET_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_HET_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()



## Export the Diff Bind regions
### Convert to GR object
sample_count_HET_greylist_LibCSIFScaled_LIB_contrast1 <- dba.report(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_HET_greylist_LibCSIFScaled_LIB_contrast2 <- dba.report(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_HET_greylist_LibCSIFScaled_LIB_contrast3 <- dba.report(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_HET_greylist_LibCSIFScaled_LIB_contrast1_df <- data.frame(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast1)
sample_count_HET_greylist_LibCSIFScaled_LIB_contrast2_df <- data.frame(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast2)
sample_count_HET_greylist_LibCSIFScaled_LIB_contrast3_df <- data.frame(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast3)

colnames(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast1_df, file="output/DiffBind/sample_count_HET_greylist_LibCSIFScaled_LIB_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast2_df, file="output/DiffBind/sample_count_HET_greylist_LibCSIFScaled_LIB_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_HET_greylist_LibCSIFScaled_LIB_contrast3_df, file="output/DiffBind/sample_count_HET_greylist_LibCSIFScaled_LIB_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)





# Time effect ~ KO genotype
## Load the greylist
load("output/DiffBind/sample_count_KO_greylist.RData")

sample_count_KO_greylist

## plot greylist
pdf("output/DiffBind/clustering_KO_greylist.pdf", width=14, height=20)
plot(sample_count_KO_greylist)
dev.off()

pdf("output/DiffBind/PCA_KO_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_greylist,DBA_REPLICATE, label=c(DBA_CONDITION))
dev.off()

## Lib-size ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; edgeR: 2457/3/2139; DESEQ2: 5916/5/5586; non-lib norm= 2581/0/2304)

sample_count_KO_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_KO_greylist, library = c(85211139,71964648,865310316,909621058,105251783,223394867), normalize = DBA_NORM_LIB) # Default


pdf("output/DiffBind/clustering_KO_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20)
plot(sample_count_KO_greylist_LibCSIFScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_KO_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_greylist_LibCSIFScaled_LIB,DBA_REPLICATE, label=DBA_CONDITION)
dev.off()

sample_count_KO_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_KO_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)




## TMM ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; DESEQ2: 3231/23/3063 )

sample_count_KO_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_KO_greylist, library = c(85211139,71964648,865310316,909621058,105251783,223394867), normalize = DBA_NORM_TMM) # Default


pdf("output/DiffBind/clustering_KO_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_KO_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_KO_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_greylist_LibCSIFScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_KO_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_KO_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



## RLE ChIPseqSpikeInFree norm (2dNvsESC/2dNvsNPC/NPCvsESC; DESEQ2: 3229/18/3035)

sample_count_KO_greylist_LibCSIFScaled_RLE = dba.normalize(sample_count_KO_greylist, library = c(85211139,71964648,865310316,909621058,105251783,223394867), normalize = DBA_NORM_RLE) # Default


pdf("output/DiffBind/clustering_KO_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20)
plot(sample_count_KO_greylist_LibCSIFScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_KO_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_greylist_LibCSIFScaled_RLE,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

sample_count_KO_greylist_LibCSIFScaled_RLE_contrast = dba.contrast(sample_count_KO_greylist_LibCSIFScaled_RLE, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_KO_greylist_LibCSIFScaled_RLE_contrast_analyze = dba.analyze(sample_count_KO_greylist_LibCSIFScaled_RLE_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



### Diff analyses with LIB/DESEQ2:


sample_count_KO_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_KO_greylist, library = c(85211139,71964648,865310316,909621058,105251783,223394867), normalize = DBA_NORM_LIB) # Default

sample_count_KO_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_KO_greylist_LibCSIFScaled_LIB, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)



# Examining results
## MA plot with diff sites
pdf("output/DiffBind/plotMA_KO_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotMA(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/plotMA_KO_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotMA(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/plotMA_KO_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotMA(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
dev.off()


## volcano plot with diff sites
pdf("output/DiffBind/volcano_KO_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_KO_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_KO_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()
## PCA/heatmat only on the diff sites
pdf("output/DiffBind/clustering_KO_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
plot(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/clustering_KO_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
plot(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=2)
dev.off()
pdf("output/DiffBind/clustering_KO_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
plot(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2, contrast=3)
dev.off()

pdf("output/DiffBind/PCA_KO_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_KO_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_KO_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, label=DBA_TREATMENT)
dev.off()


## Read concentration heatmap
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

pdf("output/DiffBind/clustering_reads_KO_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_KO_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()
pdf("output/DiffBind/clustering_reads_KO_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotHeatmap(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3, correlations=FALSE,
                              scale="row", colScheme = hmap)
dev.off()


## Read distribution within diff bound sites
pdf("output/DiffBind/BoxPlotBindAffinity_reads_KO_greylist_LibCSIFScaled_LIB_contrast1.pdf", width=14, height=20) 
dba.plotBox(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=1)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_KO_greylist_LibCSIFScaled_LIB_contrast2.pdf", width=14, height=20) 
dba.plotBox(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/BoxPlotBindAffinity_reads_KO_greylist_LibCSIFScaled_LIB_contrast3.pdf", width=14, height=20) 
dba.plotBox(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()



## Export the Diff Bind regions
### Convert to GR object
sample_count_KO_greylist_LibCSIFScaled_LIB_contrast1 <- dba.report(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_KO_greylist_LibCSIFScaled_LIB_contrast2 <- dba.report(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_KO_greylist_LibCSIFScaled_LIB_contrast3 <- dba.report(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_KO_greylist_LibCSIFScaled_LIB_contrast1_df <- data.frame(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast1)
sample_count_KO_greylist_LibCSIFScaled_LIB_contrast2_df <- data.frame(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast2)
sample_count_KO_greylist_LibCSIFScaled_LIB_contrast3_df <- data.frame(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast3)

colnames(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast1_df, file="output/DiffBind/sample_count_KO_greylist_LibCSIFScaled_LIB_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast2_df, file="output/DiffBind/sample_count_KO_greylist_LibCSIFScaled_LIB_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_KO_greylist_LibCSIFScaled_LIB_contrast3_df, file="output/DiffBind/sample_count_KO_greylist_LibCSIFScaled_LIB_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

**time course effect** per genotype:
    - WT: ESC vs NPC, ESC, 2dN,… --> (far less diff bound when using non-norm); LIB work far better; 2times more diff bound sites (DESEQ2/LIB n = 10016/0/10030; 2dNvsESC/2dNvsNPC/NPCvsESC)
    - HET: ESC vs NPC, ESC, 2dN,… --> (far less diff bound when using non-norm); LIB work far better; 2times more diff bound sites (DESEQ2/LIB n = 6354/2/6474; 2dNvsESC/2dNvsNPC/NPCvsESC)
    - KO: ESC vs NPC, ESC, 2dN,… --> (far less diff bound when using non-norm); LIB work far better; 2times more diff bound sites (DESEQ2/LIB n = 5916/5/5586)


--> Overall for all genotype almost no difference for NPC and 2dN but a lot for ESC vs NPC/2dN. Most of the sites are increased in H3K27me3 upon differentation (in agreement with the biology); a bit more mark for WT


**genotype effect over the time course**
Let's try to use the 'block' parameter (see [here](https://support.bioconductor.org/p/61354/)) to identify diff bound sites between genotypes over the time-course
```R
## Load the greylist
load("output/DiffBind/sample_count_all_greylist.RData")

sample_count_all_greylist

## plot greylist
pdf("output/DiffBind/clustering_all_greylist.pdf", width=14, height=20)
plot(sample_count_all_greylist)
dev.off()

pdf("output/DiffBind/PCA_all_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_all_greylist,DBA_TREATMENT, label=c(DBA_CONDITION))
dev.off()



## Lib-size ChIPseqSpikeInFree norm (

sample_count_all_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_all_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426,903163124,1554734579,865310316,909621058,636778842,343291724,45678566,118275343,105251783,223394867,112663109,108911344), normalize = DBA_NORM_LIB) # Default


pdf("output/DiffBind/clustering_all_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20)
plot(sample_count_all_greylist_LibCSIFScaled_LIB)
dev.off()

pdf("output/DiffBind/PCA_all_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_all_greylist_LibCSIFScaled_LIB,DBA_TREATMENT, label=DBA_CONDITION)
dev.off()


## TMM ChIPseqSpikeInFree norm 

sample_count_all_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_all_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426,903163124,1554734579,865310316,909621058,636778842,343291724,45678566,118275343,105251783,223394867,112663109,108911344), normalize = DBA_NORM_TMM) # Default


pdf("output/DiffBind/clustering_all_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_all_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_all_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_all_greylist_LibCSIFScaled_TMM,DBA_TREATMENT, label=DBA_CONDITION)
dev.off()


## RLE ChIPseqSpikeInFree norm 

sample_count_all_greylist_LibCSIFScaled_RLE = dba.normalize(sample_count_all_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426,903163124,1554734579,865310316,909621058,636778842,343291724,45678566,118275343,105251783,223394867,112663109,108911344), normalize = DBA_NORM_RLE) # Default


pdf("output/DiffBind/clustering_all_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20)
plot(sample_count_all_greylist_LibCSIFScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_all_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_all_greylist_LibCSIFScaled_RLE,DBA_TREATMENT, label=DBA_CONDITION)
dev.off()


# Time-course analyses


sample_count_all_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_all_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426,903163124,1554734579,865310316,909621058,636778842,343291724,45678566,118275343,105251783,223394867,112663109,108911344), normalize = DBA_NORM_LIB) # Default


sample_count_all_greylist_LibCSIFScaled_LIB_contrast = dba.contrast(sample_count_all_greylist_LibCSIFScaled_LIB, minMembers=2, , design = "~Condition + Treatment", reorderMeta = list(Condition="ESC",Treatment="WT"))

sample_count_all_greylist_LibCSIFScaled_LIB_contrast_analyze = dba.analyze(sample_count_all_greylist_LibCSIFScaled_LIB_contrast, method=DBA_ALL_METHODS, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


## Retrive the diff-sites for the genotypes differences (contrast 5 and 6 = HETvsWT and KOvsWT)

## Export the Diff Bind regions
### Convert to GR object
sample_count_all_greylist_LibCSIFScaled_LIB_contrast5 <- dba.report(sample_count_all_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=5)
sample_count_all_greylist_LibCSIFScaled_LIB_contrast6 <- dba.report(sample_count_all_greylist_LibCSIFScaled_LIB_contrast_analyze,method=DBA_DESEQ2,contrast=6)

### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_all_greylist_LibCSIFScaled_LIB_contrast5_df <- data.frame(sample_count_all_greylist_LibCSIFScaled_LIB_contrast5)
sample_count_all_greylist_LibCSIFScaled_LIB_contrast6_df <- data.frame(sample_count_all_greylist_LibCSIFScaled_LIB_contrast6)


colnames(sample_count_all_greylist_LibCSIFScaled_LIB_contrast5_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_all_greylist_LibCSIFScaled_LIB_contrast6_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")


write.table(sample_count_all_greylist_LibCSIFScaled_LIB_contrast5_df, file="output/DiffBind/sample_count_all_greylist_LibCSIFScaled_LIB_contrast5_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_all_greylist_LibCSIFScaled_LIB_contrast6_df, file="output/DiffBind/sample_count_all_greylist_LibCSIFScaled_LIB_contrast6_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

--> Clustering is not amazing, ESC WT clearly separated; then ESC mutants then the other samples.

--> 2dNvsESC/2dNvsNPC/ESCvsNPC/HETvsKO/HETvsWT/KOvsWT : n DESEQ2/LIB = 35256/0/35210/0/99/862

--> On IGV the diff look real

# DiffBind using macs2 raw files (As for CutRun)

Let's use the **exact same analysis method as CutRun** for convenience (macs2 raw file as input, ChIpseqSpikeInFree library scaling, TMM normalization, pvalue for Diff sites)

## DiffBind macs2 raw_Time effect

Generate the raw count matrix:


```bash
conda activate DiffBind

sbatch scripts/DiffBind_WT_macs2raw.sh # 163174 ok
sbatch scripts/DiffBind_HET_macs2raw.sh # 163176 ok
sbatch scripts/DiffBind_KO_macs2raw.sh # 163177 ok
```

Then let's make plots PCA; clustering; collect diff sites; **with pvalue 0.05**:

```R
# Genotype WT
## Load count matrix
load("output/DiffBind/sample_count_WT_macs2raw_greylist.Rdata")
# plots
pdf("output/DiffBind/clustering_WT_macs2raw_greylist.pdf", width=14, height=20)  
plot(sample_count_WT_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_WT_macs2raw_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_blackgreylist,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


# lib size scaling
sample_count_WT_blackgreylist_LibCSIFScaled_TMM = dba.normalize(sample_count_WT_blackgreylist, library = c(92471024,102739426,636778842,343291724,112663109,108911344), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_WT_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_WT_blackgreylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_WT_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_blackgreylist_LibCSIFScaled_TMM,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_WT_blackgreylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_WT_blackgreylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_WT_blackgreylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_WT_blackgreylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Genotype HET
## Load count matrix
load("output/DiffBind/sample_count_HET_macs2raw_greylist.Rdata")
# plots
pdf("output/DiffBind/clustering_HET_macs2raw_greylist.pdf", width=14, height=20)  
plot(sample_count_HET_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_HET_macs2raw_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_blackgreylist,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


# lib size scaling
sample_count_HET_blackgreylist_LibCSIFScaled_TMM = dba.normalize(sample_count_HET_blackgreylist, library = c(139649080,117855644,903163124,1554734579,45678566,118275343), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_HET_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_HET_blackgreylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_HET_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_blackgreylist_LibCSIFScaled_TMM,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_HET_blackgreylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_HET_blackgreylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_HET_blackgreylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_HET_blackgreylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Genotype KO
## Load count matrix
load("output/DiffBind/sample_count_KO_macs2raw_greylist.Rdata")
# plots
pdf("output/DiffBind/clustering_KO_macs2raw_greylist.pdf", width=14, height=20)  
plot(sample_count_KO_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_KO_macs2raw_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_blackgreylist,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


# lib size scaling
sample_count_KO_blackgreylist_LibCSIFScaled_TMM = dba.normalize(sample_count_KO_blackgreylist, library = c(85211139,71964648,865310316,909621058,105251783,223394867), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_KO_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_KO_blackgreylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_KO_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_blackgreylist_LibCSIFScaled_TMM,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_KO_blackgreylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_KO_blackgreylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_KO_blackgreylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_KO_blackgreylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

--> The mutants are more clean (clustering), thus showing more diff. bound sites over the time-course


## DiffBind macs2 raw_Genotype effect

Generate the raw count matrix:

```bash
conda activate DiffBind

sbatch scripts/DiffBind_ESC_macs2raw.sh # 163181 ok
sbatch scripts/DiffBind_NPC_macs2raw.sh # 163182 ok
sbatch scripts/DiffBind_2dN_macs2raw.sh # 163183 ok
```

```R
# Time ESC
## Load count matrix
load("output/DiffBind/sample_count_ESC_macs2raw_greylist.Rdata")
# plots
pdf("output/DiffBind/clustering_ESC_macs2raw_greylist.pdf", width=14, height=20)  
plot(sample_count_ESC_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_ESC_macs2raw_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_blackgreylist,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


# lib size scaling
sample_count_ESC_blackgreylist_LibCSIFScaled_TMM = dba.normalize(sample_count_ESC_blackgreylist, library = c(903163124, 1554734579, 865310316, 909621058, 636778842, 343291724), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_ESC_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_ESC_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM,DBA_TREATMENT, label=DBA_CONDITION)
dev.off()

# Diff. bounds

sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_ESC_blackgreylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_ESC_blackgreylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_ESC_blackgreylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Time NPC
## Load count matrix
load("output/DiffBind/sample_count_NPC_macs2raw_greylist.Rdata")
# plots
pdf("output/DiffBind/clustering_NPC_macs2raw_greylist.pdf", width=14, height=20)  
plot(sample_count_NPC_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_NPC_macs2raw_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_blackgreylist,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


# lib size scaling
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM = dba.normalize(sample_count_NPC_blackgreylist, library = c(45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_NPC_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_NPC_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_NPC", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_NPC", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_NPC", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Time 2dN
## Load count matrix
load("output/DiffBind/sample_count_2dN_macs2raw_greylist.Rdata")
# plots
pdf("output/DiffBind/clustering_2dN_macs2raw_greylist.pdf", width=14, height=20)  
plot(sample_count_2dN_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_2dN_macs2raw_greylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_blackgreylist,DBA_CONDITION, label=c(DBA_TREATMENT))
dev.off()


# lib size scaling
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM = dba.normalize(sample_count_2dN_blackgreylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_2dN_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_2dN_blackgreylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_2dN", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_2dN", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_2dN", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

--> NPC/2dN show a very bad clustering between genotypes; also high pvalue show lot of false-positive looking at some random on IGV (the low pvalue looks good otherwise)



Let's make a **quality check on NPC/2dN sample**; notably the clustering using macs2 pre-filtered peak is much better; let's check if also more diff. bound sites. If true --> Use MACS2 pre-filtered peaks for all analyses... Including for Cut Run....

Repeat analysis for NPC/2dN sample using pre-filtered macs2 peak and pvalue 0.05; compare the nb of Diff. bound sites identified and whether looks true or not.


**Genotype effect with macs2 pre-filtered and pvalue 0.05**

```R
# Time ESC
## Load count matrix
load("output/DiffBind/sample_count_ESC_greylist.Rdata")


# plots already processed (PCA_ESC_greylist.pdf and clustering_ESC_greylist.pdf)

# lib size scaling
sample_count_ESC_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_ESC_greylist, library = c(903163124, 1554734579, 865310316, 909621058, 636778842, 343291724), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_ESC_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_ESC_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_ESC_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_ESC_greylist_LibCSIFScaled_TMM,DBA_TREATMENT, label=DBA_CONDITION)
dev.off()



# Diff. bounds

sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_ESC_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_ESC_greylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_ESC_greylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_ESC_greylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast1)
sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast2)
sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Time NPC
## Load count matrix
load("output/DiffBind/sample_count_NPC_greylist.Rdata")

# plots already processed (PCA_NPC_greylist.pdf and clustering_NPC_greylist.pdf)



# lib size scaling
sample_count_NPC_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_NPC_greylist, library = c(45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_NPC_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_NPC_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_NPC_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_NPC_greylist_LibCSIFScaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_NPC_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_NPC_greylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_NPC_greylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_NPC_greylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_NPC", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_NPC", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_NPC", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Time 2dN
## Load count matrix
load("output/DiffBind/sample_count_2dN_greylist.Rdata")

# plots already processed (PCA_2dN_greylist.pdf and clustering2dN_greylist.pdf)



# lib size scaling
sample_count_2dN_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_2dN_greylist, library = c(139649080,117855644,85211139,71964648,92471024,102739426), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_2dN_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_2dN_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_2dN_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_2dN_greylist_LibCSIFScaled_TMM,DBA_TREATMENT, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_2dN_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_TREATMENT, reorderMeta = list(Treatment="WT"))

sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_2dN_greylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_2dN_greylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_2dN_greylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_2dN", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_2dN", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_2dN", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```


**Time effect with macs2 pre-filtered and pvalue 0.05**


```R
# Genotype WT
## Load count matrix
load("output/DiffBind/sample_count_WT_greylist.Rdata")



# plots already processed (PCA_WT_greylist.pdf and clustering_WT_greylist.pdf)


# lib size scaling
sample_count_WT_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_WT_greylist, library = c(92471024,102739426,636778842,343291724,112663109,108911344), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_WT_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_WT_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_WT_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_WT_greylist_LibCSIFScaled_TMM,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_WT_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_WT_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_WT_greylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_WT_greylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_WT_greylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_WT_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_WT_greylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_WT_greylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_WT_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_WT_greylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Genotype HET
## Load count matrix
load("output/DiffBind/sample_count_HET_greylist.Rdata")

# plots already processed (PCA_HET_greylist.pdf and clustering_HET_greylist.pdf)



# lib size scaling
sample_count_HET_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_HET_greylist, library = c(139649080,117855644,903163124,1554734579,45678566,118275343), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_HET_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_HET_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_HET_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_HET_greylist_LibCSIFScaled_TMM,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_HET_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_HET_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_HET_greylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_HET_greylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_HET_greylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_HET_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_HET", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_HET_greylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_HET_greylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_HET_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_HET_greylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


# Genotype KO
## Load count matrix
load("output/DiffBind/sample_count_KO_greylist.Rdata")

# plots already processed (PCA_KO_greylist.pdf and clustering_KO_greylist.pdf)



# lib size scaling
sample_count_KO_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_KO_greylist, library = c(85211139,71964648,865310316,909621058,105251783,223394867), normalize = DBA_NORM_TMM)


# plots
pdf("output/DiffBind/clustering_KO_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_KO_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_KO_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_KO_greylist_LibCSIFScaled_TMM,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()

# Diff. bounds

sample_count_KO_greylist_LibCSIFScaled_TMM_contrast = dba.contrast(sample_count_KO_greylist_LibCSIFScaled_TMM, minMembers=2, categories = DBA_CONDITION, reorderMeta = list(Condition="ESC"))

sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze = dba.analyze(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast, method=DBA_DESEQ2, bParallel = TRUE, bBlacklist=FALSE, bGreylist=FALSE)


sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze$config$bUsePval <- TRUE

dba.report(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze)


## volcano plot with diff sites
pdf("output/DiffBind/volcano_KO_greylist_LibCSIFScaled_TMM_contrast1.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2, contrast=1)
dev.off()
pdf("output/DiffBind/volcano_KO_greylist_LibCSIFScaled_TMM_contrast2.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=2)
dev.off()
pdf("output/DiffBind/volcano_KO_greylist_LibCSIFScaled_TMM_contrast3.pdf", width=14, height=20) 
dba.plotVolcano(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze, method=DBA_DESEQ2,contrast=3)
dev.off()

# Print nb of diff sites
dba.report(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=1, bGain=TRUE,bLoss=TRUE) 
dba.report(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=2, bGain=TRUE,bLoss=TRUE)
dba.report(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze, contrast=3, bGain=TRUE,bLoss=TRUE)


## Export the Diff Bind regions
### Convert to GR object
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1 <- dba.report(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=1)
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2 <- dba.report(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=2)
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3 <- dba.report(sample_count_KO_greylist_LibCSIFScaled_TMM_contrast_analyze,method=DBA_DESEQ2,contrast=3)
### Convert to bed and exportsample_model__blackgreylist_histone_norm_report_contrast1
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1_df <- data.frame(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1)
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2_df <- data.frame(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2)
sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3_df <- data.frame(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3)

colnames(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")
colnames(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3_df) <- c("seqnames", "start", "end", "width", "strand", "Conc", "Conc_KO", "Conc_KO", "Fold", "p.value", "FDR")

write.table(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast1_df, file="output/DiffBind/sample_count_KO_greylist_LibCSIFScaled_TMM_contrast1_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast2_df, file="output/DiffBind/sample_count_KO_greylist_LibCSIFScaled_TMM_contrast2_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sample_count_KO_blackgreylist_LibCSIFScaled_TMM_contrast3_df, file="output/DiffBind/sample_count_KO_greylist_LibCSIFScaled_TMM_contrast3_df.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



# Collect all scaling factor genotype/time effect:

sample_count_WT_greylist_LibCSIFScaled_TMM_SF = dba.normalize(sample_count_WT_greylist_LibCSIFScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_WT_greylist_LibCSIFScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_WT_greylist_LibCSIFScaled_TMM_SF.txt")

sample_count_HET_greylist_LibCSIFScaled_TMM_SF = dba.normalize(sample_count_HET_greylist_LibCSIFScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_HET_greylist_LibCSIFScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_HET_greylist_LibCSIFScaled_TMM_SF.txt")

sample_count_KO_greylist_LibCSIFScaled_TMM_SF = dba.normalize(sample_count_KO_greylist_LibCSIFScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_KO_greylist_LibCSIFScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_KO_greylist_LibCSIFScaled_TMM_SF.txt")

sample_count_ESC_greylist_LibCSIFScaled_TMM_SF = dba.normalize(sample_count_ESC_greylist_LibCSIFScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_ESC_greylist_LibCSIFScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_TMM_SF.txt")

sample_count_NPC_greylist_LibCSIFScaled_TMM_SF = dba.normalize(sample_count_NPC_greylist_LibCSIFScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_NPC_greylist_LibCSIFScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_TMM_SF.txt")

sample_count_2dN_greylist_LibCSIFScaled_TMM_SF = dba.normalize(sample_count_2dN_greylist_LibCSIFScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_2dN_greylist_LibCSIFScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_TMM_SF.txt")
```
- **NOTE: The pre-filtered MACS2 analysis (qvalue 0.005) are labeled `*greylist*`; the non-pre-filtered are label `*macs2raw*blacklistgreylist*` or `*blacklistgreylist*`**
- *NOTE: DiffBind scaling factor SF can be found in `output/DiffBind/*SF.txt`*



Now let's compare RNAseq (expression) and CutRun for macs2 filtered qval 0.005:
- Filter HETvsWT and KOvsWT diff bound genes into **gain and loss H3K27me3** at **each time-point**
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


# ESC



# Import peaks
peaks_contrast1 = read.table('output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast1_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast2 = read.table('output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast2_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast3 = read.table('output/DiffBind/sample_count_ESC_greylist_LibCSIFScaled_TMM_contrast3_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 


# Tidy peaks
peaks_contrast1_gr = makeGRangesFromDataFrame(peaks_contrast1,keep.extra.columns=TRUE)
peaks_contrast2_gr = makeGRangesFromDataFrame(peaks_contrast2,keep.extra.columns=TRUE)
peaks_contrast3_gr = makeGRangesFromDataFrame(peaks_contrast3,keep.extra.columns=TRUE)

gr_list <- list(HETvsKO=peaks_contrast1_gr, HETvsWT=peaks_contrast2_gr, KOvsWT=peaks_contrast3_gr)


# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_DiffBind05_qval005_ESC.pdf", width=7, height=7)
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
write.table(HETvsKO_annot, file="output/ChIPseeker/annotation_HETvsKO_qval005_ESC.txt", sep="\t", quote=F, row.names=F)
write.table(HETvsWT_annot, file="output/ChIPseeker/annotation_HETvsWT_qval005_ESC.txt", sep="\t", quote=F, row.names=F)
write.table(KOvsWT_annot, file="output/ChIPseeker/annotation_KOvsWT_qval005_ESC.txt", sep="\t", quote=F, row.names=F)


# load annotation tables
HETvsWT_annot <- read.table("output/ChIPseeker/annotation_HETvsWT_qval005_ESC.txt", sep="\t", header=TRUE)
KOvsWT_annot <- read.table("output/ChIPseeker/annotation_KOvsWT_qval005_ESC.txt", sep="\t", header=TRUE)

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
HET_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_ESC_HET_vs_ESC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_ESC_KO_vs_ESC_WT.txt')) %>%
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

pdf("output/ChIPseeker/DiffBind05_qval005_ESC_HETvsWT_expression.pdf", width=7, height=4)
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

pdf("output/ChIPseeker/DiffBind05_qval005_ESC_KOvsWT_expression.pdf", width=7, height=4)
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


# NPC
# Import peaks
peaks_contrast1 = read.table('output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast1_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast2 = read.table('output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast2_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast3 = read.table('output/DiffBind/sample_count_NPC_greylist_LibCSIFScaled_TMM_contrast3_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 


# Tidy peaks
peaks_contrast1_gr = makeGRangesFromDataFrame(peaks_contrast1,keep.extra.columns=TRUE)
peaks_contrast2_gr = makeGRangesFromDataFrame(peaks_contrast2,keep.extra.columns=TRUE)
peaks_contrast3_gr = makeGRangesFromDataFrame(peaks_contrast3,keep.extra.columns=TRUE)

gr_list <- list(HETvsKO=peaks_contrast1_gr, HETvsWT=peaks_contrast2_gr, KOvsWT=peaks_contrast3_gr)


# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_DiffBind05_qval005_NPC.pdf", width=7, height=7)
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
write.table(HETvsKO_annot, file="output/ChIPseeker/annotation_HETvsKO_qval005_NPC.txt", sep="\t", quote=F, row.names=F)
write.table(HETvsWT_annot, file="output/ChIPseeker/annotation_HETvsWT_qval005_NPC.txt", sep="\t", quote=F, row.names=F)
write.table(KOvsWT_annot, file="output/ChIPseeker/annotation_KOvsWT_qval005_NPC.txt", sep="\t", quote=F, row.names=F)


# load annotation tables
HETvsWT_annot <- read.table("output/ChIPseeker/annotation_HETvsWT_qval005_NPC.txt", sep="\t", header=TRUE)
KOvsWT_annot <- read.table("output/ChIPseeker/annotation_KOvsWT_qval005_NPC.txt", sep="\t", header=TRUE)

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
HET_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_HET_vs_NPC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_KO_vs_NPC_WT.txt')) %>%
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

pdf("output/ChIPseeker/DiffBind05_qval005_NPC_HETvsWT_expression.pdf", width=7, height=4)
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

pdf("output/ChIPseeker/DiffBind05_qval005_NPC_KOvsWT_expression.pdf", width=7, height=4)
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


# 2dN 
# Import peaks
peaks_contrast1 = read.table('output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast1_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast2 = read.table('output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast2_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast3 = read.table('output/DiffBind/sample_count_2dN_greylist_LibCSIFScaled_TMM_contrast3_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 


# Tidy peaks
peaks_contrast1_gr = makeGRangesFromDataFrame(peaks_contrast1,keep.extra.columns=TRUE)
peaks_contrast2_gr = makeGRangesFromDataFrame(peaks_contrast2,keep.extra.columns=TRUE)
peaks_contrast3_gr = makeGRangesFromDataFrame(peaks_contrast3,keep.extra.columns=TRUE)

gr_list <- list(HETvsKO=peaks_contrast1_gr, HETvsWT=peaks_contrast2_gr, KOvsWT=peaks_contrast3_gr)


# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_DiffBind05_qval005_2dN.pdf", width=7, height=7)
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
write.table(HETvsKO_annot, file="output/ChIPseeker/annotation_HETvsKO_qval005_2dN.txt", sep="\t", quote=F, row.names=F)
write.table(HETvsWT_annot, file="output/ChIPseeker/annotation_HETvsWT_qval005_2dN.txt", sep="\t", quote=F, row.names=F)
write.table(KOvsWT_annot, file="output/ChIPseeker/annotation_KOvsWT_qval005_2dN.txt", sep="\t", quote=F, row.names=F)


# load annotation tables
HETvsWT_annot <- read.table("output/ChIPseeker/annotation_HETvsWT_qval005_2dN.txt", sep="\t", header=TRUE)
KOvsWT_annot <- read.table("output/ChIPseeker/annotation_KOvsWT_qval005_2dN.txt", sep="\t", header=TRUE)

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
HET_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_2dN_HET_vs_2dN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_2dN_KO_vs_2dN_WT.txt')) %>%
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

pdf("output/ChIPseeker/DiffBind05_qval005_2dN_HETvsWT_expression.pdf", width=7, height=4)
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

pdf("output/ChIPseeker/DiffBind05_qval005_2dN_KOvsWT_expression.pdf", width=7, height=4)
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
```

Now let's compare RNAseq (expression) and CutRun for macs2raw (qval 0.05):
- Filter HETvsWT and KOvsWT diff bound genes into **gain and loss H3K27me3** at **each time-point**
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


# ESC
# Import peaks
peaks_contrast1 = read.table('output/DiffBind/sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast2 = read.table('output/DiffBind/sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast3 = read.table('output/DiffBind/sample_count_ESC_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 


# Tidy peaks
peaks_contrast1_gr = makeGRangesFromDataFrame(peaks_contrast1,keep.extra.columns=TRUE)
peaks_contrast2_gr = makeGRangesFromDataFrame(peaks_contrast2,keep.extra.columns=TRUE)
peaks_contrast3_gr = makeGRangesFromDataFrame(peaks_contrast3,keep.extra.columns=TRUE)

gr_list <- list(HETvsKO=peaks_contrast1_gr, HETvsWT=peaks_contrast2_gr, KOvsWT=peaks_contrast3_gr)


# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_DiffBind05_qval05_ESC.pdf", width=7, height=7)
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
write.table(HETvsKO_annot, file="output/ChIPseeker/annotation_HETvsKO_qval05_ESC.txt", sep="\t", quote=F, row.names=F)
write.table(HETvsWT_annot, file="output/ChIPseeker/annotation_HETvsWT_qval05_ESC.txt", sep="\t", quote=F, row.names=F)
write.table(KOvsWT_annot, file="output/ChIPseeker/annotation_KOvsWT_qval05_ESC.txt", sep="\t", quote=F, row.names=F)


# load annotation tables
HETvsWT_annot <- read.table("output/ChIPseeker/annotation_HETvsWT_qval005_ESC.txt", sep="\t", header=TRUE)
KOvsWT_annot <- read.table("output/ChIPseeker/annotation_KOvsWT_qval005_ESC.txt", sep="\t", header=TRUE)

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
HET_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_ESC_HET_vs_ESC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_ESC_KO_vs_ESC_WT.txt')) %>%
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

pdf("output/ChIPseeker/DiffBind05_qval05_ESC_HETvsWT_expression.pdf", width=7, height=4)
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

pdf("output/ChIPseeker/DiffBind05_qval05_ESC_KOvsWT_expression.pdf", width=7, height=4)
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



# NPC
# Import peaks
peaks_contrast1 = read.table('output/DiffBind/sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast2 = read.table('output/DiffBind/sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast3 = read.table('output/DiffBind/sample_count_NPC_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 


# Tidy peaks
peaks_contrast1_gr = makeGRangesFromDataFrame(peaks_contrast1,keep.extra.columns=TRUE)
peaks_contrast2_gr = makeGRangesFromDataFrame(peaks_contrast2,keep.extra.columns=TRUE)
peaks_contrast3_gr = makeGRangesFromDataFrame(peaks_contrast3,keep.extra.columns=TRUE)

gr_list <- list(HETvsKO=peaks_contrast1_gr, HETvsWT=peaks_contrast2_gr, KOvsWT=peaks_contrast3_gr)


# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_DiffBind05_qval05_NPC.pdf", width=7, height=7)
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
write.table(HETvsKO_annot, file="output/ChIPseeker/annotation_HETvsKO_qval05_NPC.txt", sep="\t", quote=F, row.names=F)
write.table(HETvsWT_annot, file="output/ChIPseeker/annotation_HETvsWT_qval05_NPC.txt", sep="\t", quote=F, row.names=F)
write.table(KOvsWT_annot, file="output/ChIPseeker/annotation_KOvsWT_qval05_NPC.txt", sep="\t", quote=F, row.names=F)


# load annotation tables
HETvsWT_annot <- read.table("output/ChIPseeker/annotation_HETvsWT_qval05_NPC.txt", sep="\t", header=TRUE)
KOvsWT_annot <- read.table("output/ChIPseeker/annotation_KOvsWT_qval05_NPC.txt", sep="\t", header=TRUE)

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
HET_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_HET_vs_NPC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_KO_vs_NPC_WT.txt')) %>%
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

pdf("output/ChIPseeker/DiffBind05_qval05_NPC_HETvsWT_expression.pdf", width=7, height=4)
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

pdf("output/ChIPseeker/DiffBind05_qval05_NPC_KOvsWT_expression.pdf", width=7, height=4)
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


# 2dN 
# Import peaks
peaks_contrast1 = read.table('output/DiffBind/sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast1_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast2 = read.table('output/DiffBind/sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast2_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 
peaks_contrast3 = read.table('output/DiffBind/sample_count_2dN_blackgreylist_LibCSIFScaled_TMM_contrast3_df.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V5, V6=V6, V7=V7, V8=V8, FC=V9, pvalue=V10, FDR=V11) 


# Tidy peaks
peaks_contrast1_gr = makeGRangesFromDataFrame(peaks_contrast1,keep.extra.columns=TRUE)
peaks_contrast2_gr = makeGRangesFromDataFrame(peaks_contrast2,keep.extra.columns=TRUE)
peaks_contrast3_gr = makeGRangesFromDataFrame(peaks_contrast3,keep.extra.columns=TRUE)

gr_list <- list(HETvsKO=peaks_contrast1_gr, HETvsWT=peaks_contrast2_gr, KOvsWT=peaks_contrast3_gr)


# Overlap assigned genes btwn my gr_list
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Fine-tune here gene peak assignemnt

genes= lapply(peakAnnoList, function(i) unique(as.data.frame(i)$geneId))

pdf("output/ChIPseeker/overlap_genes_DiffBind05_qval05_2dN.pdf", width=7, height=7)
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
write.table(HETvsKO_annot, file="output/ChIPseeker/annotation_HETvsKO_qval05_2dN.txt", sep="\t", quote=F, row.names=F)
write.table(HETvsWT_annot, file="output/ChIPseeker/annotation_HETvsWT_qval05_2dN.txt", sep="\t", quote=F, row.names=F)
write.table(KOvsWT_annot, file="output/ChIPseeker/annotation_KOvsWT_qval05_2dN.txt", sep="\t", quote=F, row.names=F)


# load annotation tables
HETvsWT_annot <- read.table("output/ChIPseeker/annotation_HETvsWT_qval05_2dN.txt", sep="\t", header=TRUE)
KOvsWT_annot <- read.table("output/ChIPseeker/annotation_KOvsWT_qval05_2dN.txt", sep="\t", header=TRUE)

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
HET_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_2dN_HET_vs_2dN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

KO_vs_WT = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_2dN_KO_vs_2dN_WT.txt')) %>%
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

pdf("output/ChIPseeker/DiffBind05_qval05_2dN_HETvsWT_expression.pdf", width=7, height=4)
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

pdf("output/ChIPseeker/DiffBind05_qval05_2dN_KOvsWT_expression.pdf", width=7, height=4)
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
```

--> The macs2raw (qval 0.05) show much more DEGs and Diff. bound sites; however many goes in the wrong direction (eg. gain H3K27me3 up expre)

--> We could use macs2raw (qval 0.05) and follow on the one that goes in the right direction?

--> To identify diff. bound sites and associated genes; let's use the macs2raw qval 0.05 and focus on the gene where expression is changed in agreement

--> to show clustering let's use the qval 0.005 pre-filtered

### DiffBind_Collect scaling factor (SF) to generate clean DiffBind TMM bigwig

Collect all samples; 2 replicates per samples (not include the 3 ESC WT Rep). 

--> Greylist filtered counts (using qval 0.005) has already been generated as `save(sample_count_all_greylist, file = "output/DiffBind/sample_count_all_greylist.RData")`; so let's import it; normalized with library size/TMM; do PCA/clustering and collect SF

```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```

```R
# Load package
library("DiffBind")

# Load greylist count
load("output/DiffBind/sample_count_all_greylist.RData")

# Apply library-size correction and TMM normalization
sample_count_all_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_all_greylist, library = c(139649080, 117855644, 85211139, 71964648, 92471024, 102739426, 903163124, 1554734579, 865310316, 909621058, 636778842, 343291724, 45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_TMM) # TMM norm
sample_count_all_greylist_LibCSIFScaled_LIB = dba.normalize(sample_count_all_greylist, library = c(139649080, 117855644, 85211139, 71964648, 92471024, 102739426, 903163124, 1554734579, 865310316, 909621058, 636778842, 343291724, 45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_LIB)
sample_count_all_greylist_LibCSIFScaled_RLE = dba.normalize(sample_count_all_greylist, library = c(139649080, 117855644, 85211139, 71964648, 92471024, 102739426, 903163124, 1554734579, 865310316, 909621058, 636778842, 343291724, 45678566, 118275343, 105251783, 223394867, 112663109, 108911344), normalize = DBA_NORM_RLE)

# Plot PCA / clustering
pdf("output/DiffBind/clustering_all_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20)
plot(sample_count_all_greylist_LibCSIFScaled_TMM)
dev.off()
pdf("output/DiffBind/clustering_all_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20)
plot(sample_count_all_greylist_LibCSIFScaled_LIB)
dev.off()
pdf("output/DiffBind/clustering_all_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20)
plot(sample_count_all_greylist_LibCSIFScaled_RLE)
dev.off()

pdf("output/DiffBind/PCA_all_greylist_LibCSIFScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_all_greylist_LibCSIFScaled_TMM,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_all_greylist_LibCSIFScaled_LIB.pdf", width=14, height=20) 
dba.plotPCA(sample_count_all_greylist_LibCSIFScaled_LIB,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()
pdf("output/DiffBind/PCA_all_greylist_LibCSIFScaled_RLE.pdf", width=14, height=20) 
dba.plotPCA(sample_count_all_greylist_LibCSIFScaled_RLE,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()

# Collect scaling factor
sample_count_all_greylist_LibCSIFScaled_TMM_SF = dba.normalize(sample_count_all_greylist_LibCSIFScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_all_greylist_LibCSIFScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_all_greylist_LibCSIFScaled_TMM_SF.txt")

sample_count_all_greylist_LibCSIFScaled_LIB_SF = dba.normalize(sample_count_all_greylist_LibCSIFScaled_LIB, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_all_greylist_LibCSIFScaled_LIB_SF))
writeLines(console_output, "output/DiffBind/sample_count_all_greylist_LibCSIFScaled_LIB_SF.txt")

sample_count_all_greylist_LibCSIFScaled_RLE_SF = dba.normalize(sample_count_all_greylist_LibCSIFScaled_RLE, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_all_greylist_LibCSIFScaled_RLE_SF))
writeLines(console_output, "output/DiffBind/sample_count_all_greylist_LibCSIFScaled_RLE_SF.txt")
```
--> TMM and RLE SF are comparable. LIB is very different

--> The clustering is OK, not amazing; at least ESC clustered together, WT clsutered together then HET/KO/NPC/2dN a bit mixed; LIB/RLE/TMM perform similarly; so stay with TMM

- sample / library size * SF = scaled library size --> DiffBind_TMM_SF ; reciprocal
- 2dN_HET_H3K27me3_R1 / 70887858 1.97 = 139649080 --> 0.6596462 ; 1.515964164
- 2dN_HET_H3K27me3_R2 / 67346082 1.75 = 117855644 --> 0.6285088 ; 1.591067619
- 2dN_KO_H3K27me3_R1 / 58363794 1.46 = 85211139 --> 0.5460104 ; 1.831466946
- 2dN_KO_H3K27me3_R2 / 71964648 1 = 71964648 --> 0.5574423 ; 1.793907639
- 2dN_WT_H3K27me3_R1 / 71682964 1.29 = 92471024 --> 0.6458845 ; 1.548264434
- 2dN_WT_H3K27me3_R2 / 60792560 1.69 = 102739426 --> 0.5910808 ; 1.691816077
- ESC_HET_H3K27me3_R1 / 85933694 10.51 = 903163124 --> 0.4633394 ; 2.158245122
- ESC_HET_H3K27me3_R2 / 66583922 23.35 = 1554734579 --> 0.3576010 ; 2.796412762
- ESC_KO_H3K27me3_R1 / 86014942 10.06 = 865310316 --> 0.4468665 ; 2.237804803
- ESC_KO_H3K27me3_R2 / 57643920 15.78 = 909621058 --> 0.3132557 ; 3.1922803
- ESC_WT_H3K27me3_R1 / 90968406 7 = 636778842 --> 0.7957969 ; 1.25660203
- ESC_WT_H3K27me3_R2 / 79650052 4.31 = 343291724 --> 0.7333544 ; 1.363597191
- NPC_HET_H3K27me3_R1 / 40423510 1.13 = 45678566 --> 0.3083840 ; 3.242710387
- NPC_HET_H3K27me3_R2 / 82710030 1.43 = 118275343 --> 0.7630974 ; 1.310448705
- NPC_KO_H3K27me3_R1 / 67904376 1.55 = 105251783 --> 0.5464982 ; 1.829832193
- NPC_KO_H3K27me3_R2 / 84619268 2.64 = 223394867 --> 0.7914168 ; 1.2635567
- NPC_WT_H3K27me3_R1 / 77698696 1.45 = 112663109 --> 0.7014795 ; 1.425558409
- NPC_WT_H3K27me3_R2 / 72126718 1.51 = 108911344 --> 0.6423245 ; 1.556845489






Generate **Bigwig DiffBind TMM scaled** (LibCSIFScaled_TMM):

Use reciprocal of SF !!

```bash
conda activate deeptools

sbatch scripts/bamtobigwig_DiffBind_TMM_ESC.sh # 894655
sbatch scripts/bamtobigwig_DiffBind_TMM_NPC.sh # 894656
sbatch scripts/bamtobigwig_DiffBind_TMM_2dN.sh # 894657
sbatch scripts/bamtobigwig_DiffBind_TMM_input.sh # 894659

# median
conda activate BedToBigwig
sbatch --dependency=afterany:894655:894656:894657:894659 scripts/bigwigmerge_DiffBind_TMM.sh # 894697
```
- *NOTE: I also generated input DiffBind_LibCSIFScaled_TMM corrected*

--> The replicates looks OK

--> We do NOT see the increase of H3K27me3 in ESC versus diff samples; something is wrong! (maybe we cannot apply both ChIPseqSpikeInFree and TMM; I tested TMM only next)







# ChIPseeker
[Tutorial](http://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html) and [documentation]()
## ChIPseeker installation
In conda base; Install within R 4.2.2 module
```R
BiocManager::Install("ChIPseeker")
library("ChIPseeker")
```

**Re-install on the new cluster RES-RHEL-RH9HPC:**
Installed together with deseq2; use **deseq2 conda env for ChIPseeker now**


## Run ChIPseeker

XXX


# Pool the peak into 1 file
We could use IDR to pool our replicates and identify 'high-confidence' peak (used in the ENCODE pipeline)

[Tutorial](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html) and [github](https://github.com/nboley/idr) for IDR.

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
But for homogeneity as I cannot use IDR for CutRun as more than 4 replicates let's use instead macs2 directly.


Re-run MACS2 but in pool to have 1 file per condition (keep blacklist and qval2.30103 (0.005) filtering):
```bash
sbatch scripts/macs2_pool.sh # 15059 ok
```
*NOTE: for ESC WT I took Rep1 and Rep2 only (as Rep3 is another WT clone)*


Now let's filter out blacklist and qvalue:
```bash
sbatch scripts/macs2_pool_peak_signif.sh # ok
```
*NOTE: I relaunch the script and changed the qvalue; 2.30103 (q0.005) is good as observed looking at individual replicates*



## Explore with DeepTools


Let's try to use deepTools to explore the data: tutorial [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html) and [here](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html#reference-point).



Let's compare deepTools for `bigwig_DiffBind_TMM` versus `bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF` (check replicates)


```bash
conda activate deeptools
# deepTools plot
## bigwig_DiffBind_TMM
sbatch --dependency=afterany:894697 scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_ESC_noIntergenic_Rep.sh # 894880 ok
sbatch --dependency=afterany:894697 scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_NPC_noIntergenic_Rep.sh # 894900 ok
sbatch --dependency=afterany:894697 scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_2dN_noIntergenic_Rep.sh # 894901 ok


## bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF
sbatch scripts/matrix_gene_1kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_ESC_noIntergenic_Rep.sh # 894902 ok
sbatch scripts/matrix_gene_1kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_NPC_noIntergenic_Rep.sh # 894904 ok
sbatch scripts/matrix_gene_1kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_2dN_noIntergenic_Rep.sh # 894905 ok
```
*NOTE: for comparison here I used the `ENCFF159KBI_peak_noIntergenic.gtf` (peak in at least 1 genotype non intergenic)*

--> `bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF` is ok for replicates (not perfect but overall seems ok); except 
NPC KO very bad. Also I feel ESC WT is VERY different as compare to ESC KO and HET; the difference is HUGE with WT much more H3K27me3 (maybe true but that's strong)


--> `bigwig_DiffBind_TMM` is overall better for replicates clustering; at all time-point. However; the conclusion differ for ESC:
- `bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF` = WT is MORE H3K27me3 / KO-HET
- `bigwig_DiffBind_TMM` = WT is LESS H3K27me3 / KO-HET


--> `bigwig_DiffBind_TMM` show that ESC samples are MORE H3K27me3 than NPC and 2dN. That is not true, so something wrong with the `bigwig_DiffBind_TMM` method...


Maybe we cannot do both ChIPseqInFree and DiffBind TMM correction; may result in a weird double correction... Let's try using DiffBind TMM without scaling the library and see if the SF we obtain increase H3K27me3 upon differentiation


```bash
srun --mem=500g --pty bash -l
conda activate DiffBind
```

```R
# Load package
library("DiffBind")

# Load greylist count
load("output/DiffBind/sample_count_all_greylist.RData")

# Apply ONLY TMM normalization
sample_count_all_greylist_LibCSIFScaled_TMM = dba.normalize(sample_count_all_greylist, normalize = DBA_NORM_TMM) # TMM norm


# Plot PCA / clustering
pdf("output/DiffBind/clustering_all_greylist_TMM.pdf", width=14, height=20)
plot(sample_count_all_greylist_LibCSIFScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_all_greylist_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_all_greylist_LibCSIFScaled_TMM,DBA_CONDITION, label=DBA_TREATMENT)
dev.off()


# Collect scaling factor
sample_count_all_greylist_LibCSIFScaled_TMM_SF = dba.normalize(sample_count_all_greylist_LibCSIFScaled_TMM, bRetrieve=TRUE)
console_output <- capture.output(print(sample_count_all_greylist_LibCSIFScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_all_greylist_TMM_SF.txt")
```


- sample / library size * SF = scaled library size --> DiffBind_TMM_SF ; reciprocal | DiffBind_TMM_ONLY_SF ; reciprocal \ CSsIF_DiffBind_LIB-rec
- 2dN_HET_H3K27me3_R1 / 70887858 1.97 = 139649080 --> 0.6596462 ; 1.515964164 | 1.1452498 ; 0.8731719490367953 \ 0.3845337-2.600552305298599
- 2dN_HET_H3K27me3_R2 / 67346082 1.75 = 117855644 --> 0.6285088 ; 1.591067619 | 1.0914982 ; 0.9161719185611117 \ 0.3245239-3.081437145307326
- 2dN_KO_H3K27me3_R1 / 58363794 1.46 = 85211139 --> 0.5460104 ; 1.831466946 | 0.9471807 ; 1.055764755341827 \ 0.2346349-4.261940572353047
- 2dN_KO_H3K27me3_R2 / 71964648 1 = 71964648 --> 0.5574423 ; 1.793907639 | 0.9744307 ; 1.026240244688514 \ 0.1981598-5.046432222882744
- 2dN_WT_H3K27me3_R1 / 71682964 1.29 = 92471024 --> 0.6458845 ; 1.548264434 | 1.1265964 ; 0.887629323154237 \ 0.2546255-3.927336421528873
- 2dN_WT_H3K27me3_R2 / 60792560 1.69 = 102739426 --> 0.5910808 ; 1.691816077 | 1.0306665 ; 0.9702459524977284 \ 0.2829003-3.5348142083978
- ESC_HET_H3K27me3_R1 / 85933694 10.51 = 903163124 --> 0.4633394 ; 2.158245122 | 0.8051032 ; 1.24207679214292 \ 2.4869239-0.4021031765386951
- ESC_HET_H3K27me3_R2 / 66583922 23.35 = 1554734579 --> 0.3576010 ; 2.796412762 |  0.6234852 ; 1.603887309594518 \ 4.2810723-0.2335863377032899
- ESC_KO_H3K27me3_R1 / 86014942 10.06 = 865310316 --> 0.4468665 ; 2.237804803 | 0.7811012; 1.280243840362811 \ 2.3826935-0.4196930910333201
- ESC_KO_H3K27me3_R2 / 57643920 15.78 = 909621058 --> 0.3132557 ; 3.1922803 | 0.5454599 ; 1.833315336287782 \ 2.5047063-0.3992484068890632
- ESC_WT_H3K27me3_R1 / 90968406 7 = 636778842 --> 0.7957969 ; 1.25660203 |  1.3835209 ; 0.7227935624246804 \ 1.7534159-0.5703153484578302
- ESC_WT_H3K27me3_R2 / 79650052 4.31 = 343291724 --> 0.7333544 ; 1.363597191 | 1.2753760 ; 0.7840824980241121 \ 0.9452782-1.057889624451299
- NPC_HET_H3K27me3_R1 / 40423510 1.13 = 45678566 --> 0.3083840 ; 3.242710387 |  0.5362265 ; 1.864883589304147 \ 0.1257792-7.950440136365949
- NPC_HET_H3K27me3_R2 / 82710030 1.43 = 118275343 --> 0.7630974 ; 1.310448705 |  1.3247137 ; 0.754880092204074 \ 0.3256796-3.070502420170008
- NPC_KO_H3K27me3_R1 / 67904376 1.55 = 105251783 --> 0.5464982 ; 1.829832193 |  0.9526521 ; 1.04970114483556 \ 0.2898183-3.450437739783858
- NPC_KO_H3K27me3_R2 / 84619268 2.64 = 223394867 --> 0.7914168 ; 1.2635567 |  1.3738308 ; 0.7278916734142225 \ 0.6151337-1.625662843703735
- NPC_WT_H3K27me3_R1 / 77698696 1.45 = 112663109 --> 0.7014795 ; 1.425558409 |  1.2176666 ; 0.8212428590880295 \ 0.3102259-3.223457486947415
- NPC_WT_H3K27me3_R2 / 72126718 1.51 = 108911344 --> 0.6423245 ; 1.556845489 | 1.1169502 ; 0.8952950632893033 \ 0.2998951-3.334499296587373

--> Median value for DiffBind_TMM_ONLY_SF; ESC: 1.26116031625287, NPC: 0.858268961188666, 2dN: 0.9462061048648195. Which is good as it will decrease signal in ESC


Now generate **bigwig DiffBind_TMM_ONLY (non ChIPseqSpikeInFree + TMM) and deepTools plots**:

*NOTE: Use reciprocal of SF (from the SF given by DiffBind)!!*

```bash
conda activate deeptools

sbatch scripts/bamtobigwig_DiffBind_TMM_ONLY_ESC.sh # 1011491 ok
sbatch scripts/bamtobigwig_DiffBind_TMM_ONLY_NPC.sh # 1011508 ok
sbatch scripts/bamtobigwig_DiffBind_TMM_ONLY_2dN.sh # 1011514 ok
sbatch scripts/bamtobigwig_DiffBind_TMM_ONLY_input.sh # 1011585 ok

# median
conda activate BedToBigwig
sbatch --dependency=afterany:1011491:1011508:1011514:1011585 scripts/bigwigmerge_DiffBind_TMM_ONLY.sh # 1011632 ok

conda activate deeptools
# deepTools plot
## bigwig_DiffBind_TMM
sbatch --dependency=afterany:1011632 scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_ONLY_ESC_noIntergenic_Rep.sh # 1011757 ok
sbatch --dependency=afterany:1011632 scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_ONLY_NPC_noIntergenic_Rep.sh # 1011830 ok
sbatch --dependency=afterany:1011632 scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_ONLY_2dN_noIntergenic_Rep.sh # 1011903 ok
```
*NOTE: I also generated input DiffBind_LibCSIFScaled_TMM corrected*

--> The replicates looks good

--> We do NOT see the increase of H3K27me3 in ESC versus diff samples...

The **ChIPseqSpikeInFree+DiffBind_LIB normalization** seems to provide correct SF; let's try it out:


```bash
conda activate deeptools

sbatch scripts/bamtobigwig_DiffBind_LIB_ESC.sh # 1063646 ok
sbatch scripts/bamtobigwig_DiffBind_LIB_NPC.sh # 1063648 ok
sbatch scripts/bamtobigwig_DiffBind_LIB_2dN.sh # 1063662 ok
sbatch scripts/bamtobigwig_DiffBind_LIB_input.sh # 1063684 ok

# median
conda activate BedToBigwig
sbatch --dependency=afterany:1063646:1063648:1063662:1063684 scripts/bigwigmerge_DiffBind_LIB.sh # 1063692 ok

conda activate deeptools
# deepTools plot
## bigwig_DiffBind_LIB
sbatch --dependency=afterany:1063692 scripts/matrix_gene_1kb_bigwig_DiffBind_LIB_ESC_noIntergenic_Rep.sh # 1063693 ok
sbatch --dependency=afterany:1063692 scripts/matrix_gene_1kb_bigwig_DiffBind_LIB_NPC_noIntergenic_Rep.sh # 1063694 ok
sbatch --dependency=afterany:1063692 scripts/matrix_gene_1kb_bigwig_DiffBind_LIB_2dN_noIntergenic_Rep.sh # 1063696 ok
```
*NOTE: I also generated input DiffBind_LibCSIFScaled_LIB corrected*


--> The replicates are not amazing, looks as in the ChIPseqSpikeInFree bigwigs

--> We do see the increase of H3K27me3 in ESC versus diff samples, which is good

--> At ESC, mutants are very low in H3K27me3 vs WT

## QC ChIPSeq and RNAseq

Let's do **QC for the ESC state**, let' compare RNAseq and ChIPseq profile at ESC; the genes up-regulated should show decrease of H3K27me3, the genes down-regulated should show increase of H3K27me3:

- Generate single gtf files that contains only:
    - up-regulated genes in HET at ESC
    - up-regulated genes in KO at ESC
    - down-regulated genes in HET at ESC
    - down-regulated genes in KO at ESC
- Generate gtf file that only contains gene assign with H3K27me3 in WT (and another file for in at least 1 genotype)
- Put together gtf expression and gtf peak assign in WT (or any other genotype)
- deepTool profile on these gtf files using `Bigwig_LibCSIFScaled_LIB`



### single gtf files for expression at ESC,NPC,2dN and 4wN

Filter the gene Up/Down in each mutant in R and generate GTF:
```bash
conda activate deseq2
cd /scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq
```


```R
# library
library("rtracklayer")
library("tidyverse")

# ESC time DEGs mutants 
## Import deseq2 output and filter qvalue 0.05
HETvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_ESC_HET_vs_ESC_WT.txt')) %>%
    filter(padj <= 0.05) %>%
    dplyr::select(-"...1") %>%
    add_column(contrast = "HETvsWT")
KOvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_ESC_KO_vs_ESC_WT.txt')) %>%
    filter(padj <= 0.05) %>%
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
export(gtf_DEGs_HET_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_HET_Down.gtf")
export(gtf_DEGs_HET_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_HET_Up.gtf")
export(gtf_DEGs_KO_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_KO_Down.gtf")
export(gtf_DEGs_KO_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_KO_Up.gtf")

 
# NPC time DEGs mutants 
## Import deseq2 output and filter qvalue 0.05
HETvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_NPC_HET_vs_NPC_WT.txt')) %>%
    filter(padj <= 0.05) %>%
    dplyr::select(-"...1") %>%
    add_column(contrast = "HETvsWT")
KOvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_NPC_KO_vs_NPC_WT.txt')) %>%
    filter(padj <= 0.05) %>%
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
export(gtf_DEGs_HET_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_HET_Down.gtf")
export(gtf_DEGs_HET_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_HET_Up.gtf")
export(gtf_DEGs_KO_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_KO_Down.gtf")
export(gtf_DEGs_KO_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_KO_Up.gtf")


# 2dN time DEGs mutants 
## Import deseq2 output and filter qvalue 0.05
HETvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_2dN_HET_vs_2dN_WT.txt')) %>%
    filter(padj <= 0.05) %>%
    dplyr::select(-"...1") %>%
    add_column(contrast = "HETvsWT")
KOvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_2dN_KO_vs_2dN_WT.txt')) %>%
    filter(padj <= 0.05) %>%
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
export(gtf_DEGs_HET_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_2dN_HET_Down.gtf")
export(gtf_DEGs_HET_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_2dN_HET_Up.gtf")
export(gtf_DEGs_KO_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_2dN_KO_Down.gtf")
export(gtf_DEGs_KO_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_2dN_KO_Up.gtf")


# 4wN time DEGs mutants 
## Import deseq2 output and filter qvalue 0.05
HETvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_4wN_HET_vs_4wN_WT.txt')) %>%
    filter(padj <= 0.05) %>%
    dplyr::select(-"...1") %>%
    add_column(contrast = "HETvsWT")
KOvsWT = as_tibble(read_csv('output/deseq2_hg38/raw_4wN_KO_vs_4wN_WT.txt')) %>%
    filter(padj <= 0.05) %>%
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
export(gtf_DEGs_HET_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_4wN_HET_Down.gtf")
export(gtf_DEGs_HET_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_4wN_HET_Up.gtf")
export(gtf_DEGs_KO_Down, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_4wN_KO_Down.gtf")
export(gtf_DEGs_KO_Up, con = "output/deseq2_hg38/ENCFF159KBI_DEGs_4wN_KO_Up.gtf")
```


deepTools profile for these DEGs only (**if result weird; let's filter out the genes that contain peak in at least 1 genotype or in WT only at ESC**)


```bash
# deepTools plot
conda activate deeptools
## ESC DEGs
## bigwig_DiffBind_LIB keeping replicates
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_LIB_DEGs_ESC_HET_Down_Rep.sh # 1072657
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_LIB_DEGs_ESC_HET_up_Rep.sh # 1072658
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_LIB_DEGs_ESC_KO_Down_Rep.sh # 1072664
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_LIB_DEGs_ESC_KO_up_Rep.sh # 1072666 FAIl gtf; 

## bigwig_DiffBind_TMM keeping replicates
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_DEGs_ESC_HET_Down_Rep.sh # 1074398
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_DEGs_ESC_HET_up_Rep.sh # 1074399
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_DEGs_ESC_KO_Down_Rep.sh # 1074400
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_DEGs_ESC_KO_up_Rep.sh # 1074402
```
**bigwig_DiffBind_LIB:**
--> Looking at HET or KO samples only, we indeed see a changes in aggrement with gene expression (eg. more express in HET show less H3K27me3 signal)

--> When looking WT vs mutant we still have the WT MUCH higher than the mutants... Something is wrong.

--> Western blot do not show any drastic differences at ESC WT vs mutants so its unlikely to have STRONG differences...


XXX Maybe we should look at the ratio with input?

XXX Let's filter-out to keep only genes enriched in WT tu rule-out this possibility


**bigwig_DiffBind_TMM:**

--> Seems to be working: when more express in mutant (DEGs Up), H3K27me3 is higher in WT; When less express in mutants (DEGs Less) H3K27me3 is comparable WT-mutants

--> Mutants profile a bit heterogeneous and weird upstream TSS

Let's filter out and keep only the genes that contain H3K27me3 peak in WT or at least 1 genotype


### ChIPseeker - Assign peak to genes each time-point

Let's:
- assign peak to genes at each time-point
- generate gtf file with these coordinates genes (WT only, and all genotypes)
- Overlap gtf peak assign with gtf DEGs
- Profile deepTools H3K27me3



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

# ESC
# Define a list of file paths
files <- c('output/macs2/broad_blacklist_qval2.30103/ESC_WT_H3K27me3_pool_peaks.broadPeak', 
           'output/macs2/broad_blacklist_qval2.30103/ESC_HET_H3K27me3_pool_peaks.broadPeak', 
           'output/macs2/broad_blacklist_qval2.30103/ESC_KO_H3K27me3_pool_peaks.broadPeak')

names <- c("WT", "HET", "KO")

# Function to read, process and save annotation
process_file <- function(file, name){
  
  # Read file
  peaks <- read.table(file) %>% 
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
  
  # Make GRanges object
  peaks_gr = makeGRangesFromDataFrame(peaks,keep.extra.columns=TRUE)
  
  # Annotate peaks
  peaks_annot = annotatePeak(peaks_gr, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
  
  # Get annotation data frame
  annot_df <- as.data.frame(peaks_annot@anno)
  
  # Convert entrez gene IDs to gene symbols
  annot_df$geneSymbol <- mapIds(org.Hs.eg.db, keys = annot_df$geneId, column = "SYMBOL", keytype = "ENTREZID")
  annot_df$gene <- mapIds(org.Hs.eg.db, keys = annot_df$geneId, column = "ENSEMBL", keytype = "ENTREZID")
  
  # Save output table
  write.table(annot_df, file=paste0("output/ChIPseeker/annotation_ESC_", name, ".txt"), sep="\t", quote=F, row.names=F)
}

# Apply the function to all files
mapply(process_file, file = files, name = names)


# NPC
# Define a list of file paths
files <- c('output/macs2/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_pool_peaks.broadPeak', 
           'output/macs2/broad_blacklist_qval2.30103/NPC_HET_H3K27me3_pool_peaks.broadPeak', 
           'output/macs2/broad_blacklist_qval2.30103/NPC_KO_H3K27me3_pool_peaks.broadPeak')

names <- c("WT", "HET", "KO")

# Function to read, process and save annotation
process_file <- function(file, name){
  
  # Read file
  peaks <- read.table(file) %>% 
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
  
  # Make GRanges object
  peaks_gr = makeGRangesFromDataFrame(peaks,keep.extra.columns=TRUE)
  
  # Annotate peaks
  peaks_annot = annotatePeak(peaks_gr, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
  
  # Get annotation data frame
  annot_df <- as.data.frame(peaks_annot@anno)
  
  # Convert entrez gene IDs to gene symbols
  annot_df$geneSymbol <- mapIds(org.Hs.eg.db, keys = annot_df$geneId, column = "SYMBOL", keytype = "ENTREZID")
  annot_df$gene <- mapIds(org.Hs.eg.db, keys = annot_df$geneId, column = "ENSEMBL", keytype = "ENTREZID")
  
  # Save output table
  write.table(annot_df, file=paste0("output/ChIPseeker/annotation_NPC_", name, ".txt"), sep="\t", quote=F, row.names=F)
}

# Apply the function to all files
mapply(process_file, file = files, name = names)



# 2dN
# Define a list of file paths
files <- c('output/macs2/broad_blacklist_qval2.30103/2dN_WT_H3K27me3_pool_peaks.broadPeak', 
           'output/macs2/broad_blacklist_qval2.30103/2dN_HET_H3K27me3_pool_peaks.broadPeak', 
           'output/macs2/broad_blacklist_qval2.30103/2dN_KO_H3K27me3_pool_peaks.broadPeak')

names <- c("WT", "HET", "KO")

# Function to read, process and save annotation
process_file <- function(file, name){
  
  # Read file
  peaks <- read.table(file) %>% 
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4)
  
  # Make GRanges object
  peaks_gr = makeGRangesFromDataFrame(peaks,keep.extra.columns=TRUE)
  
  # Annotate peaks
  peaks_annot = annotatePeak(peaks_gr, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
  
  # Get annotation data frame
  annot_df <- as.data.frame(peaks_annot@anno)
  
  # Convert entrez gene IDs to gene symbols
  annot_df$geneSymbol <- mapIds(org.Hs.eg.db, keys = annot_df$geneId, column = "SYMBOL", keytype = "ENTREZID")
  annot_df$gene <- mapIds(org.Hs.eg.db, keys = annot_df$geneId, column = "ENSEMBL", keytype = "ENTREZID")
  
  # Save output table
  write.table(annot_df, file=paste0("output/ChIPseeker/annotation_2dN_", name, ".txt"), sep="\t", quote=F, row.names=F)
}

# Apply the function to all files
mapply(process_file, file = files, name = names)
```
--> I used a function to assign peak to genes and save

**Filter-out intergenic-peak, convert assnotation file to gtfs**



```bash
conda activate BedToBigwig
# File where the diffbound sites has been assigned to genes
output/ChIPseeker/annotation_ESC_WT.txt
output/ChIPseeker/annotation_ESC_HET.txt
output/ChIPseeker/annotation_ESC_KO.txt
output/ChIPseeker/annotation_NPC_WT.txt
output/ChIPseeker/annotation_NPC_HET.txt
output/ChIPseeker/annotation_NPC_KO.txt
output/ChIPseeker/annotation_2dN_WT.txt
output/ChIPseeker/annotation_2dN_HET.txt
output/ChIPseeker/annotation_2dN_KO.txt

# Define array of sample names
samples=("ESC_WT" "ESC_HET" "ESC_KO" "NPC_WT" "NPC_HET" "NPC_KO" "2dN_WT" "2dN_HET" "2dN_KO")

# Iterate over each sample
for sample in ${samples[@]}
do
    # Filter out intergenic regions
    grep -v "Intergenic" output/ChIPseeker/annotation_${sample}.txt > output/ChIPseeker/annotation_${sample}_noIntergenic.bed

    # Collect gene ID
    awk -F'\t' '(NR==1 || FNR>1) {print $21}' output/ChIPseeker/annotation_${sample}_noIntergenic.bed | sort | uniq > output/ChIPseeker/annotation_${sample}_noIntergenic_geneSymbol.txt

    # Modify the txt file
    sed 's/^/gene_name "/; s/$/"/' output/ChIPseeker/annotation_${sample}_noIntergenic_geneSymbol.txt > output/ChIPseeker/annotation_${sample}_noIntergenic_as_gtf_geneSymbol.txt

    # Filter the gtf
    grep -Ff output/ChIPseeker/annotation_${sample}_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_${sample}_noIntergenic.gtf
done

# Combine GTF for at least in 1 genotype for each time-point
## ESC
### Combine the three gtfs and output only gene names
cat meta/ENCFF159KBI_ESC_WT_noIntergenic.gtf meta/ENCFF159KBI_ESC_HET_noIntergenic.gtf meta/ENCFF159KBI_ESC_KO_noIntergenic.gtf | sort | uniq | cut -f9 | sed -n 's/.*gene_name "\([^"]*\)".*/\1/p' | sort | uniq > meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic_geneSymbol.txt
### filter-in the gtf the gene names
###  Modify the txt file
sed 's/^/gene_name "/; s/$/"/' meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic_geneSymbol.txt > meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic_as_gtf_geneSymbol.txt
###  Filter the gtf
grep -Ff meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic.gtf
## NPC
### Combine the three gtfs and output only gene names
cat meta/ENCFF159KBI_NPC_WT_noIntergenic.gtf meta/ENCFF159KBI_NPC_HET_noIntergenic.gtf meta/ENCFF159KBI_NPC_KO_noIntergenic.gtf | sort | uniq | cut -f9 | sed -n 's/.*gene_name "\([^"]*\)".*/\1/p' | sort | uniq > meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic_geneSymbol.txt
### filter-in the gtf the gene names
###  Modify the txt file
sed 's/^/gene_name "/; s/$/"/' meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic_geneSymbol.txt > meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic_as_gtf_geneSymbol.txt
###  Filter the gtf
grep -Ff meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic.gtf
## 2dN
### Combine the three gtfs and output only gene names
cat meta/ENCFF159KBI_2dN_WT_noIntergenic.gtf meta/ENCFF159KBI_2dN_HET_noIntergenic.gtf meta/ENCFF159KBI_2dN_KO_noIntergenic.gtf | sort | uniq | cut -f9 | sed -n 's/.*gene_name "\([^"]*\)".*/\1/p' | sort | uniq > meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic_geneSymbol.txt
### filter-in the gtf the gene names
###  Modify the txt file
sed 's/^/gene_name "/; s/$/"/' meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic_geneSymbol.txt > meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic_as_gtf_geneSymbol.txt
###  Filter the gtf
grep -Ff meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic.gtf





# Combine GTF for at least in 1 genotype at any time-point
### Combine the three gtfs and output only gene names
cat meta/ENCFF159KBI_ESC_WT_noIntergenic.gtf meta/ENCFF159KBI_ESC_HET_noIntergenic.gtf meta/ENCFF159KBI_ESC_KO_noIntergenic.gtf meta/ENCFF159KBI_NPC_WT_noIntergenic.gtf meta/ENCFF159KBI_NPC_HET_noIntergenic.gtf meta/ENCFF159KBI_NPC_KO_noIntergenic.gtf meta/ENCFF159KBI_2dN_WT_noIntergenic.gtf meta/ENCFF159KBI_2dN_HET_noIntergenic.gtf meta/ENCFF159KBI_2dN_KO_noIntergenic.gtf | sort | uniq | cut -f9 | sed -n 's/.*gene_name "\([^"]*\)".*/\1/p' | sort | uniq > meta/ENCFF159KBI_ESC_NPC_2dN_WT_HET_KO_noIntergenic_geneSymbol.txt
### filter-in the gtf the gene names
###  Modify the txt file
sed 's/^/gene_name "/; s/$/"/' meta/ENCFF159KBI_ESC_NPC_2dN_WT_HET_KO_noIntergenic_geneSymbol.txt > meta/ENCFF159KBI_ESC_NPC_2dN_WT_HET_KO_noIntergenic_as_gtf_geneSymbol.txt
###  Filter the gtf
grep -Ff meta/ENCFF159KBI_ESC_NPC_2dN_WT_HET_KO_noIntergenic_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_ESC_NPC_2dN_WT_HET_KO_noIntergenic.gtf


# Combine GTF for in WT for each time-point
XXX


# Combine GTF for in WT at any time-point
XXX


# Combine DEGs with gene-peak assign
## per time-point
### ESC_ Up/down express and gene-peak assigned (at least 1 genotype)
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_HET_Down.gtf -b meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic_DEGs_ESC_HET_Down.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_HET_Up.gtf -b meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic_DEGs_ESC_HET_Up.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_KO_Down.gtf -b meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic_DEGs_ESC_KO_Down.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_ESC_KO_Up.gtf -b meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_ESC_WT_HET_KO_noIntergenic_DEGs_ESC_KO_Up.gtf
### NPC_ Up/down express and gene-peak assigned (at least 1 genotype)
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_HET_Down.gtf -b meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic_DEGs_NPC_HET_Down.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_HET_Up.gtf -b meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic_DEGs_NPC_HET_Up.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_KO_Down.gtf -b meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic_DEGs_NPC_KO_Down.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_NPC_KO_Up.gtf -b meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_NPC_WT_HET_KO_noIntergenic_DEGs_NPC_KO_Up.gtf
### 2dN_ Up/down express and gene-peak assigned (at least 1 genotype)
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_2dN_HET_Down.gtf -b meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic_DEGs_2dN_HET_Down.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_2dN_HET_Up.gtf -b meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic_DEGs_2dN_HET_Up.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_2dN_KO_Down.gtf -b meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic_DEGs_2dN_KO_Down.gtf
bedtools intersect -wa -u -a ../001__RNAseq/output/deseq2_hg38/ENCFF159KBI_DEGs_2dN_KO_Up.gtf -b meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_2dN_WT_HET_KO_noIntergenic_DEGs_2dN_KO_Up.gtf
```



Generate deepTools plots (**DEGs with gene peak assigned**):


```bash
# deepTools plot
conda activate deeptools
## ESC
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_ESC_WT_HET_KO_noIntergenic_DEGs_ESC_HET_Down.sh # 1076547 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_ESC_WT_HET_KO_noIntergenic_DEGs_ESC_HET_Up.sh # 1076548 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_ESC_WT_HET_KO_noIntergenic_DEGs_ESC_KO_Down.sh # 1076549 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_ESC_WT_HET_KO_noIntergenic_DEGs_ESC_KO_Up.sh # 1076550 FAIL; 1076584 ok
## NPC
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_NPC_WT_HET_KO_noIntergenic_DEGs_NPC_HET_Down.sh # 1076554 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_NPC_WT_HET_KO_noIntergenic_DEGs_NPC_HET_Up.sh # 1076556 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_NPC_WT_HET_KO_noIntergenic_DEGs_NPC_KO_Down.sh # 1076557 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_NPC_WT_HET_KO_noIntergenic_DEGs_NPC_KO_Up.sh # 1076558 ok
## 2dN
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_2dN_WT_HET_KO_noIntergenic_DEGs_2dN_HET_Down.sh # 1076559 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_2dN_WT_HET_KO_noIntergenic_DEGs_2dN_HET_Up.sh # 1076560 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_2dN_WT_HET_KO_noIntergenic_DEGs_2dN_KO_Down.sh # 1076561 ok
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_2dN_WT_HET_KO_noIntergenic_DEGs_2dN_KO_Up.sh # 1076562 ok
```

--> ESC; DEGs Down looks good, we see higher H3K27me3 signal in mutants. DEGs Up is not super clear; WT vs mutants show very different profile

--> NPC; HET is always higher than WT (even when genes are downregulated...). KO is working great (KO DEGs Up show less H3K27me3 and DEGs Down show more H3K27me3)

--> 2dN; HET is always higher than WT (even when genes are downregulated...). KO is working great (KO DEGs Up show less H3K27me3 and DEGs Down show more H3K27me3)

--> Overall seems to be working pretty great; as expected. 

XXX could try to look at TSS only instead of gene body to avoid the very different H3K27me3 distribution profile WT vs mutants at ESC





### Generate deepTools plot for time-course clustered genes

As quality control let's check cluster 13 and 18:

- Collect gene list from each cluster
- Combine with the gtf (genes assigned to a peak at any time point in any genotypes: `meta/ENCFF159KBI_ESC_NPC_2dN_WT_HET_KO_noIntergenic.gtf`)
- deepTools plot


```bash
conda activate BedToBigwig
# Generate list of gene in ENSG format from teh gene cluster file
## Extract gene from cluster 13 and 18 as a list of genes
awk -F, '$3==13 {gsub(/"/, "", $1); print $1}' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_25cl.txt > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_25cl_genesFromCl13.txt
awk -F, '$3==18 {gsub(/"/, "", $1); print $1}' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_25cl.txt > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_25cl_genesFromCl18.txt

## Generate the gtf
### filter-in the gtf the gene names
###  Modify the txt file
sed 's/^/gene_id "/; s/$/"/' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_25cl_genesFromCl13.txt > meta/cluster_gene_rlog_25cl_genesFromCl13_as_gtf_geneId.txt
sed 's/^/gene_id "/; s/$/"/' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_25cl_genesFromCl18.txt > meta/cluster_gene_rlog_25cl_genesFromCl18_as_gtf_geneId.txt
###  Filter the gtf
grep -Ff meta/cluster_gene_rlog_25cl_genesFromCl13_as_gtf_geneId.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_cluster_gene_rlog_25cl_genesFromCl13.gtf
grep -Ff meta/cluster_gene_rlog_25cl_genesFromCl18_as_gtf_geneId.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_cluster_gene_rlog_25cl_genesFromCl18.gtf
### Only keep the genes that have a H3K27me3 peak in at least 1 gentoype in at least 1 time-point
bedtools intersect -wa -u -a meta/ENCFF159KBI_cluster_gene_rlog_25cl_genesFromCl13.gtf -b meta/ENCFF159KBI_ESC_NPC_2dN_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_cluster_gene_rlog_25cl_genesFromCl13_ESC_NPC_2dN_WT_HET_KO_noIntergenic.gtf
bedtools intersect -wa -u -a meta/ENCFF159KBI_cluster_gene_rlog_25cl_genesFromCl18.gtf -b meta/ENCFF159KBI_ESC_NPC_2dN_WT_HET_KO_noIntergenic.gtf > meta/ENCFF159KBI_cluster_gene_rlog_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic.gtf


# deepTools plot
conda activate deeptools
## Per genotype over the time course
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl13_ESC_NPC_2dN_WT_HET_KO_noIntergenic_WT.sh # 1076713
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl13_ESC_NPC_2dN_WT_HET_KO_noIntergenic_HET.sh # 1076714
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl13_ESC_NPC_2dN_WT_HET_KO_noIntergenic_KO.sh # 1076715

sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic_WT.sh # 1076717
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic_HET.sh # 1076718
sbatch scripts/matrix_gene_1kb_bigwig_DiffBind_TMM_25cl_genesFromCl18_ESC_NPC_2dN_WT_HET_KO_noIntergenic_KO.sh # 1076719


```

--> The ESC is weird again; always higher in H3K27me3; even for cluster 13 where expression decrease..

--> Otherwise; NPC to 2dN is logical as expected for each cluster. 

--> Still, overall profile is weird, look very noisy. Need to find a way to clean our data!


Try to use the **uniquely mapped reads only** (reads I used to generate the ChIpseqSpikeInFree scaling factor)

- Collect DiffBind_TMM scaling factor but generated with the bam unique aligned reads only (not all MAPQ20 per default)
- Apply DiffBind TMM scaling factor on the unique aligned read files and generate bigwig files
- Check their profile with deepTools


```bash
conda activate DiffBind

XXX

```
