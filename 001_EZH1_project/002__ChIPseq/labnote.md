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
--> The number of reads is similar in the 2 (paired)-read; and reduced/2 when looking at only 1 of the 2 lanes; so concatenation keep data integrity and combine all reads. Just in case, backup has been performed from all files: **To be deleted after mapping (XXX)**


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
- param2/**permissive-unpaired** `--phred33 -q --local --no-unal --dovetail`: (can have uniquely map paired reads) XXX 
    - nb of uniquely mapped reads: 23071868 (51.18%)
    - 5743226 (34.22%) >1 times
    - 15429070 (12.74%) 0 times
    - overall 93.43%


--> Default is better! Better overall alignment rate, better concordant read alignment (with more uniquely mapped reads).

--> Probably because end-to-end method is used and not --local, aligning the entire read, not cliped may help the alignment. And we allow mix alignment

--> Let's try the following: `--phred33 -q --no-unal --no-mixed --dovetail` **endtoend**. Here it is more clean, we do not allow mix alignemnt and keep dovetail option possible.

```bash
sbatch bowtie2_2dN_HET_H3K27me3_R1_endtoend.sh # 11496376
```
- param2/**endtoend** `--phred33 -q --no-unal --no-mixed --dovetail`: XXX 
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
library("tidyverse")
# Create sample_meta.txt; tab delimited format `output/ChIPseqSpikeInFree/sample_meta.txt
metaFile <- "output/ChIPseqSpikeInFree/sample_meta.txt"
bams <- c("output/bowtie2_endtoend/2dN_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_HET_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_HET_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_HET_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_H3K27me3_R3.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_input_R3.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_input_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_input_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_input_R2.dupmark.sorted.bam")


# Run ChIPSpikeInFree
ChIPseqSpikeInFree(bamFiles = bams, chromFile = "hg38", metaFile = metaFile, prefix = "all_sample")
```
Job get terminated; let's run it within a script
```bash
conda activate ChIPseqSpikeInFree
sbatch scripts/ChIPseqSpikeInFree_all.sh # 12100044
```
XXX


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



# ChIPseeker
[Tutorial](http://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html) and [documentation]()
## ChIPseeker installation
In conda base; Install within R 4.2.2 module
```R
BiocManager::Install("ChIPseeker")
library("ChIPseeker")
```

## Run ChIPseeker







