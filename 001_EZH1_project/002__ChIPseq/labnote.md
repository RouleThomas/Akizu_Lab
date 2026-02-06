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

# Re-launch pgbam with uniquely aligned reads for ESC and NPC WT and KO only (test 013/001)
sbatch scripts/samtools_highquality_ESCNPCWTKO.sh # 65660504 ok; but fail naming at samtools sort 
sbatch scripts/samtools_highquality_ESCNPCWTKO_corr.sh # 65792238 xxx

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
*NOTE: for comparison here I used the `ENCFF159KBI_peak_noIntergenic.gtf` (peak in at least 1 genotype non intergenic; **from CutRun**)*

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


YYY Maybe we should look at the ratio with input?

YYY Let's filter-out to keep only genes enriched in WT tu rule-out this possibility


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

YYY


# Combine GTF for in WT at any time-point
YYY


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

YYY could try to look at TSS only instead of gene body to avoid the very different H3K27me3 distribution profile WT vs mutants at ESC





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


# Use the uniquely mapped reads only (reads I used to generate the ChIpseqSpikeInFree scaling factor) for downstream analysis
- Generate macs2 analysis from uniquely mapped reads bam file
- Collect DiffBind_TMM scaling factor but generated with the bam unique aligned reads only (not all MAPQ20 per default)
- Apply DiffBind TMM scaling factor on the unique aligned read files and generate bigwig files
- Check their profile with deepTools
- Try using THOR with these scaling factors if looks good

## MACS2 peak calling

in `output/macs2_unique`

```bash
conda activate macs2

# call peak per replicate
sbatch scripts/macs2_ESC_unique.sh # 1077363 ok
sbatch scripts/macs2_NPC_unique.sh # 1077364 ok
sbatch scripts/macs2_2dN_unique.sh # 1077365 ok
## stat
conda activate BedToBigwig
sbatch scripts/macs2_peak_signif_unique.sh # Re-run for different qval

# call peak per sample
sbatch scripts/macs2_pool_unique_unique.sh # 1077367 ok
## stat
sbatch scripts/macs2_pool_peak_signif_unique.sh # Re-run for different qval
```

--> All good, many LESS peaks identified...



## DiffBind_TMM (uniquely mapped reads)

Generate new meta file in `output/DiffBind/meta_sample_all_macs2raw_unique.txt`

```bash
conda activate DiffBind
# scripts for counting/greylist/blacklist
sbatch scripts/DiffBind_all_macs2raw_unique.sh # 1123762 ok
```

Let's apply the ChIPseqSpikeInFree scaling factor to normalize data (normalize library size, then normalize per seq. depth; as for the CutRun). **THIS TIME WE NORMALIZE THE UNIQUELY ALIGNED DATA**. ChIPseqSpikeInFree-norm-library-size = library-size * SF. `samtools flagstat output/bowtie2_endtoend/*unique.dupmark.sorted.bam` used to obtain library size (first value=library size):
- sample / library_size_UNIQUE * SF = scaled library size
- 2dN_HET_H3K27me3_R1 / 40818892 1.97 = 80413217
- 2dN_HET_H3K27me3_R2 / 38757156 1.75 = 67825023
- 2dN_KO_H3K27me3_R1 / 36478530 1.46 = 53258654
- 2dN_KO_H3K27me3_R2 / 46363646 1 = 46363646
- 2dN_WT_H3K27me3_R1 / 43088276 1.29 = 55583876
- 2dN_WT_H3K27me3_R2 / 43478292 1.69 = 73478313
- ESC_HET_H3K27me3_R1 / 44064716 10.51 = 463120165
- ESC_HET_H3K27me3_R2 / 29954726 23.35 = 699442852
- ESC_KO_H3K27me3_R1 / 27538348 10.06 = 277035781
- ESC_KO_H3K27me3_R2 / 23084678 15.78 = 364276219
- ESC_WT_H3K27me3_R1 / 57347074 7 = 401429518
- ESC_WT_H3K27me3_R2 / 58089688 4.31 = 250366555
- NPC_HET_H3K27me3_R1 / 28967884 1.13 = 32733709
- NPC_HET_H3K27me3_R2 / 53067608 1.43 = 75886679
- NPC_KO_H3K27me3_R1 / 42573738 1.55 = 65989294
- NPC_KO_H3K27me3_R2 / 48730202 2.64 = 128647733
- NPC_WT_H3K27me3_R1 / 52824596 1.45 = 76595664
- NPC_WT_H3K27me3_R2 / 51005778 1.51 = 77018725



```R
library("DiffBind") 
# Load counting files
load("output/DiffBind/count_all_macs2raw_unique_blackgreylist.RData")


# plot
pdf("output/DiffBind/clustering_all_macs2raw_unique_blackgreylist.pdf", width=14, height=20)
plot(sample_count_blackgreylist)
dev.off()

pdf("output/DiffBind/PCA_all_macs2raw_unique_blackgreylist.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()

# Apply TMM normalization 
 
# TMM norm
sample_count_blackgreylist_LibHistoneScaled_TMM = dba.normalize(sample_count_blackgreylist, library = c(80413217, 67825023,53258654, 46363646,55583876,73478313, 463120165, 699442852,277035781, 364276219, 401429518, 250366555, 32733709, 75886679, 65989294, 128647733, 76595664, 77018725), normalize = DBA_NORM_TMM)

pdf("output/DiffBind/clustering_all_macs2raw_unique_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20)
plot(sample_count_blackgreylist_LibHistoneScaled_TMM)
dev.off()

pdf("output/DiffBind/PCA_all_macs2raw_unique_blackgreylist_LibHistoneScaled_TMM.pdf", width=14, height=20) 
dba.plotPCA(sample_count_blackgreylist_LibHistoneScaled_TMM,DBA_REPLICATE, label=DBA_TREATMENT)
dev.off()


## Here is to retrieve the scaling factor value
sample_count_blackgreylist_LibHistoneScaled_TMM_SF = dba.normalize(sample_count_blackgreylist_LibHistoneScaled_TMM, bRetrieve=TRUE)


console_output <- capture.output(print(sample_count_blackgreylist_LibHistoneScaled_TMM_SF))
writeLines(console_output, "output/DiffBind/sample_count_all_macs2raw_unique_blackgreylist_LibHistoneScaled_TMM_SF.txt")
```
--> The PCA/clustering is OK; genotype well clustered for ESC, then 2dN and NPC are mixed; WT tends to be together then the other genotypes are mixed...


Here is the SF and updated table for **uniquely aligned reads**:

- sample / library_size_UNIQUE * SF = scaled library size = DiffBind_TMM_SF / Reciprocal_DiffBind_TMM_SF
- 2dN_HET_H3K27me3_R1 / 40818892 1.97 = 80413217 = 0.7527256 / 1.328505368
- 2dN_HET_H3K27me3_R2 / 38757156 1.75 = 67825023 = 0.7213495 / 1.386290557
- 2dN_KO_H3K27me3_R1 / 36478530 1.46 = 53258654 = 0.6986599 / 1.431311572
- 2dN_KO_H3K27me3_R2 / 46363646 1 = 46363646 = 0.7226374 / 1.38381988
- 2dN_WT_H3K27me3_R1 / 43088276 1.29 = 55583876 = 0.7803635 / 1.281454092
- 2dN_WT_H3K27me3_R2 / 43478292 1.69 = 73478313 = 0.9313169 / 1.073748366
- ESC_HET_H3K27me3_R1 / 44064716 10.51 = 463120165 = 0.4348701 / 2.299537264
- ESC_HET_H3K27me3_R2 / 29954726 23.35 = 699442852 = 0.2897601 / 3.45113078
- ESC_KO_H3K27me3_R1 / 27538348 10.06 = 277035781 = 0.2699819 / 3.703952006
- ESC_KO_H3K27me3_R2 / 23084678 15.78 = 364276219 = 0.2244825 / 4.454690232
- ESC_WT_H3K27me3_R1 / 57347074 7 = 401429518 = 0.9900251 / 1.010075401
- ESC_WT_H3K27me3_R2 / 58089688 4.31 = 250366555 = 1.1046768 / 0.905242149
- NPC_HET_H3K27me3_R1 / 28967884 1.13 = 32733709 = 0.4501949 / 2.22126017
- NPC_HET_H3K27me3_R2 / 53067608 1.43 = 75886679 = 1.0141581 / 0.986039553
- NPC_KO_H3K27me3_R1 / 42573738 1.55 = 65989294 = 0.6830601 / 1.46400002
- NPC_KO_H3K27me3_R2 / 48730202 2.64 = 128647733 = 0.9034272 / 1.106896051
- NPC_WT_H3K27me3_R1 / 52824596 1.45 = 76595664 = 0.9716131 / 1.029216259
- NPC_WT_H3K27me3_R2 / 51005778 1.51 = 77018725 = 0.9549168 / 1.047211652

```
$norm.factors
 [1] 0.7527256 0.7213495 0.6986599 0.7226374 0.7803635 0.9313169 0.4348701
 [8] 0.2897601 0.2699819 0.2244825 0.9900251 1.1046768 0.4501949 1.0141581
[15] 0.6830601 0.9034272 0.9716131 0.9549168
```

### Bigwig_ ChIPseqInFree-DiffBind-TMM-uniqueReads

Apply these ChIPseqInFree-DiffBind-TMM-uniqueReads new SF and check whether deepTools plot are weird or not (**apply the reciprocal**)...

In `output/bigwig_UniqueBamUniqueSF_DiffBind_TMM`

```bash
conda activate deeptools
sbatch scripts/bamtobigwig_UniqueBamUniqueSF_DiffBind_TMM.sh # 1652289 cancel time limit; 
sbatch scripts/bamtobigwig_UniqueBamUniqueSF_DiffBind_TMM_1.sh # 1654848
sbatch scripts/bamtobigwig_UniqueBamUniqueSF_DiffBind_TMM_2.sh # 1654849
sbatch scripts/bamtobigwig_UniqueBamUniqueSF_DiffBind_TMM_3.sh # 1654850
# deepTools plot
## For ESC only, all genotypes
### Less express in HET
sbatch --dependency=afterany:1654848:1654849:1654850 scripts/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_DEGs_ESC_HET_Down_Rep.sh # 1652290 FAIL time limit with 1652289; 1654864
### More express in HET
sbatch --dependency=afterany:1654848:1654849:1654850 scripts/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_DEGs_ESC_HET_Up_Rep.sh # 1652314; 1654920


## per genotype over time-course; all genes
sbatch --dependency=afterany:1654848:1654849:1654850 scripts/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_WT_Rep.sh # 1652366; 1654922
sbatch --dependency=afterany:1654848:1654849:1654850 scripts/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_HET_Rep.sh # 1652367; 1654927
sbatch --dependency=afterany:1654848:1654849:1654850 scripts/matrix_gene_1kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_KO_Rep.sh # 1652370; 1654934




```
*NOTE: for comparison here I used the `ENCFF159KBI_peak_noIntergenic.gtf` (peak in at least 1 genotype non intergenic; **from CutRun**)*

--> Mutants profile (HET/KO) looks like SHIT; same weird profile, WT is good

--> At ESC DEGs HET; almost no changes of H3K27me3 WT vs HET (ie. WT always more H3K27me3). It's SHIT!!!

--> We observe decrease H3K27me3 from ESC to NPC/2dN for all 3 genotypes. BAD!!

Seems that using the **`bigwig_UniqueBamUniqueSF_DiffBind_TMM` do NOT improve the analysis** (HET and KO profile are weird for all time-points + we do NOT observe increase H3K27me3 upon diff.)

# THOR diff. bound sites
Let's try to use **THOR with TMM-scaled scaling factor**, calculated from uniquely aligned reads bam with DiffBind_TMM; or from non-uniquely aligned reads bam with DiffBind_TMM. Both **scaling factor are initially from ChIPseqSpikeInFree** and re-scaled the library before applying the TMM nornmalization.

*NOTE: In THOR; use scaling factor value used to convert Bam to Bigwig (reciprocal value of DiffBind output)*

DiffBind_TMM SF from **NON-uniquely aligned bam** (used to generate `output/bigwig_DiffBind_TMM`)
- 2dN_HET_H3K27me3_R1 = 1.515964164
- 2dN_HET_H3K27me3_R2 = 1.591067619
- 2dN_KO_H3K27me3_R1 = 1.831466946
- 2dN_KO_H3K27me3_R2 = 1.793907639
- 2dN_WT_H3K27me3_R1 = 1.548264434
- 2dN_WT_H3K27me3_R2 = 1.691816077
- ESC_HET_H3K27me3_R1 = 2.158245122
- ESC_HET_H3K27me3_R2 = 2.796412762
- ESC_KO_H3K27me3_R1 = 2.237804803
- ESC_KO_H3K27me3_R2 = 3.1922803
- ESC_WT_H3K27me3_R1 = 1.25660203
- ESC_WT_H3K27me3_R2 = 1.363597191
- NPC_HET_H3K27me3_R1 = 3.242710387
- NPC_HET_H3K27me3_R2 = 1.310448705
- NPC_KO_H3K27me3_R1 = 1.829832193
- NPC_KO_H3K27me3_R2 = 1.2635567
- NPC_WT_H3K27me3_R1 = 1.425558409
- NPC_WT_H3K27me3_R2 = 1.556845489

DiffBind_TMM SF from **uniquely aligned bam** (used to generate `output/bamtobigwig_UniqueBamUniqueSF_DiffBind_TMM`)
- 2dN_HET_H3K27me3_R1 = 1.328505368
- 2dN_HET_H3K27me3_R2 = 1.386290557
- 2dN_KO_H3K27me3_R1 = 1.431311572
- 2dN_KO_H3K27me3_R2 = 1.38381988
- 2dN_WT_H3K27me3_R1 = 1.281454092
- 2dN_WT_H3K27me3_R2 = 1.073748366
- ESC_HET_H3K27me3_R1 = 2.299537264
- ESC_HET_H3K27me3_R2 = 3.45113078
- ESC_KO_H3K27me3_R1 = 3.703952006
- ESC_KO_H3K27me3_R2 = 4.454690232
- ESC_WT_H3K27me3_R1 = 1.010075401
- ESC_WT_H3K27me3_R2 = 0.905242149
- NPC_HET_H3K27me3_R1 = 2.22126017
- NPC_HET_H3K27me3_R2 = 0.986039553
- NPC_KO_H3K27me3_R1 = 1.46400002
- NPC_KO_H3K27me3_R2 = 1.106896051
- NPC_WT_H3K27me3_R1 = 1.029216259
- NPC_WT_H3K27me3_R2 = 1.047211652

**ChIPseqSpikeInFree SF_RAW** so not pass-by DiffBind, no TMM-normalized (used to generate `output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF`)
- 2dN_HET_H3K27me3_R1 = 1.97
- 2dN_HET_H3K27me3_R2 = 1.75
- 2dN_KO_H3K27me3_R1 = 1.46
- 2dN_KO_H3K27me3_R2 = 1
- 2dN_WT_H3K27me3_R1 = 1.29
- 2dN_WT_H3K27me3_R2 = 1.69
- ESC_HET_H3K27me3_R1 = 10.51
- ESC_HET_H3K27me3_R2 = 23.35
- ESC_KO_H3K27me3_R1 = 10.06
- ESC_KO_H3K27me3_R2 = 15.78
- ESC_WT_H3K27me3_R1 = 7
- ESC_WT_H3K27me3_R2 = 4.31
- NPC_HET_H3K27me3_R1 = 1.13
- NPC_HET_H3K27me3_R2 = 1.43
- NPC_KO_H3K27me3_R1 = 1.55
- NPC_KO_H3K27me3_R2 = 2.64
- NPC_WT_H3K27me3_R1 = 1.45
- NPC_WT_H3K27me3_R2 = 1.51

**ChIPseqSpikeInFree SF_RAW** But Reciprocal!
- 2dN_HET_H3K27me3_R1 = 1.97 = 0.5076142131979695
- 2dN_HET_H3K27me3_R2 = 1.75 = 0.5714285714285714
- 2dN_KO_H3K27me3_R1 = 1.46 = 0.6849315068493151
- 2dN_KO_H3K27me3_R2 = 1 = 1
- 2dN_WT_H3K27me3_R1 = 1.29 = 0.7751937984496124
- 2dN_WT_H3K27me3_R2 = 1.69 = 0.5917159763313609
- ESC_HET_H3K27me3_R1 = 10.51 = 0.0951474785918173
- ESC_HET_H3K27me3_R2 = 23.35 = 0.0428265524625268
- ESC_KO_H3K27me3_R1 = 10.06 = 0.099403578528827
- ESC_KO_H3K27me3_R2 = 15.78 = 0.0633713561470215
- ESC_WT_H3K27me3_R1 = 7 = 0.1428571428571429
- ESC_WT_H3K27me3_R2 = 4.31 = 0.2320185614849188
- NPC_HET_H3K27me3_R1 = 1.13 = 0.8849557522123894
- NPC_HET_H3K27me3_R2 = 1.43 = 0.6993006993006993
- NPC_KO_H3K27me3_R1 = 1.55 = 0.6451612903225806
- NPC_KO_H3K27me3_R2 = 2.64 = 0.3787878787878788
- NPC_WT_H3K27me3_R1 = 1.45 = 0.6896551724137931
- NPC_WT_H3K27me3_R2 = 1.51 = 0.6622516556291391

**ChIPseqSpikeInFree SF_FromUniqueBAM** But Reciprocal!
- 2dN_HET_H3K27me3_R1  5.360881149 = 0.186536499
- 2dN_HET_H3K27me3_R2  4.5216682	= 0.221157315345
- 2dN_KO_H3K27me3_R1  3.55057692 = 0.28164437
- 2dN_KO_H3K27me3_R2  3.090909733 = 0.323529345
- 2dN_WT_H3K27me3_R1  3.705591736 = 0.269862433
- 2dN_WT_H3K27me3_R2  4.898554232 = 0.204141866
- ESC_HET_H3K27me3_R1  30.87467768 = 0.032389002
- ESC_HET_H3K27me3_R2  46.62952347 = 0.021445641
- ESC_KO_H3K27me3_R1  18.46905206 = 0.054144631
- ESC_KO_H3K27me3_R2  24.28508126 = 0.041177544
- ESC_WT_H3K27me3_R1  26.76196787 = 0.03736646
- ESC_WT_H3K27me3_R2  16.69110369 = 0.059912156
- NPC_HET_H3K27me3_R1  2.182247261 = 0.458243215
- NPC_HET_H3K27me3_R2  5.059111963 = 0.197663149
- NPC_KO_H3K27me3_R1  4.39928626 = 0.2273096
- NPC_KO_H3K27me3_R2  8.576515552 = 0.116597468
- NPC_WT_H3K27me3_R1  5.106377613 = 0.195833539
- NPC_WT_H3K27me3_R2  5.134581652 = 0.194757834




**Create the CONFIG file** with nano; *(I put them in meta in VSC/Github but are in output/THOR in HPC)*:
- `output/THOR/ESC_WTvsHET_UniqueBamDiffBindTMM.config` # Light test
- `output/THOR/ESC_WTvsHET_DiffBindTMM.config` # Light test
- `output/THOR/ESC_WTvsKO_UniqueBamDiffBindTMM.config` # Light test
- `output/THOR/ESC_WTvsKO_DiffBindTMM.config` # Light test
- `output/THOR/NPC_WTvsHET.config`
- `output/THOR/NPC_WTvsKO.config`
- `output/THOR/2dN_WTvsHET.config`
- `output/THOR/2dN_WTvsKO.config`
- `output/THOR/WT_ESCvsNPC_UniqueBamDiffBindTMM.config` # Light test
- `output/THOR/WT_ESCvsNPC_DiffBindTMM.config` # Light test

- `output/THOR/WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree.config`; not created as = `output/THOR/WT_ESCvsNPC_UniqueBamDiffBindTMM.config` # Light test
- `output/THOR/WT_ESCvsNPC_ChIPseqSpikeInFree.config` ; not created as = `output/THOR/WT_ESCvsNPC_DiffBindTMM.config` # Light test

- `output/THOR/WT_ESCvs2dN_UniqueBamDiffBindTMM.config` # Light test
- `output/THOR/WT_ESCvs2dN_DiffBindTMM.config` # Light test
- `output/THOR/HET_ESCvsNPC.config` 
- `output/THOR/HET_ESCvs2dN.config`
- `output/THOR/KO_ESCvsNPC.config`
- `output/THOR/KO_ESCvs2dN.config`

--> Do the *# Light test* first and check wehter it is good; test both SF and corresponding bam files for each *# Light test* sample




## THOR with different SF (DiffBindTMM from uniquely/NON-uniquely aligned bam AND ChIPseqSpikeInFree ony)
- try with DiffBind_TMM where ChIPseqSpikeInFree was applied on library size
- try with ChIPseqSpikeInFree SF; from `bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF`; they are NOT TMM-normalized though... But as discussed earlier maybe we cannot use both ChIPseqSpikeInFree and TMM normalization...
- try without scaling factor (TMM default normalziation)
- try with housekeeping genes normalization

*THOR is very buggy to make it work I need to temporaly change where to look for libraries lol.. So cannot use nano anymore for example...*

```bash
# Needed step to change where THOR look for libraries
conda activate RGT
export LD_LIBRARY_PATH=~/anaconda3/envs/RGT/lib:$LD_LIBRARY_PATH
bigWigMerge # for testing

# run 200g mem
conda activate RGT
## Genotype comparison at ESC
sbatch scripts/THOR_ESC_WTvsHET_UniqueBamDiffBindTMM.sh # 1655774 ok
sbatch scripts/THOR_ESC_WTvsHET_DiffBindTMM.sh # 1655779 ok
sbatch scripts/THOR_ESC_WTvsKO_UniqueBamDiffBindTMM.sh # 1655812 ok
sbatch scripts/THOR_ESC_WTvsKO_DiffBindTMM.sh # 1655834 ok
## Genotype comparison at ESC, Default TMM-normalization (NO SF)
sbatch scripts/THOR_ESC_WTvsHET_UniqueBamTMM.sh # 1673898
sbatch scripts/THOR_ESC_WTvsHET_TMM.sh # 1673902
sbatch scripts/THOR_ESC_WTvsKO_UniqueBamTMM.sh # 1673903
sbatch scripts/THOR_ESC_WTvsKO_TMM.sh # 1673934
## Genotype comparison at ESC ChIPseqSpikeInFree
sbatch scripts/THOR_ESC_WTvsHET_UniqueBamChIPseqSpikeInFree_Corr.sh # 1661296
sbatch scripts/THOR_ESC_WTvsHET_ChIPseqSpikeInFree_Corr.sh # 1661299
sbatch scripts/THOR_ESC_WTvsKO_UniqueBamChIPseqSpikeInFree_Corr.sh # 1661300
sbatch scripts/THOR_ESC_WTvsKO_ChIPseqSpikeInFree_Corr.sh # 1661301
## Time-effect for WT
sbatch scripts/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM.sh # 1655862 ok
sbatch scripts/THOR_WT_ESCvsNPC_DiffBindTMM.sh # 1655868 ok
sbatch scripts/THOR_WT_ESCvs2dN_UniqueBamDiffBindTMM.sh # 1655880
sbatch scripts/THOR_WT_ESCvs2dN_DiffBindTMM.sh # 1655883
## Time-effect for WT ChIPseqSpikeInFree
sbatch scripts/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree.sh # 1656686 ok
sbatch scripts/THOR_WT_ESCvsNPC_ChIPseqSpikeInFree.sh # 1656714 ok
sbatch scripts/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Corr.sh # 1661293 ok
sbatch scripts/THOR_WT_ESCvsNPC_ChIPseqSpikeInFree_Corr.sh # 1661294 ok
## Time-effect for WT, Default TMM-normalization (NO SF)
sbatch scripts/THOR_WT_ESCvsNPC_UniqueBamTMM.sh # 1673882 ok
sbatch scripts/THOR_WT_ESCvsNPC_TMM.sh # 1673885 ok
## Time-effect for WT, Default TMM-normalization (NO SF) with rmdup
sbatch scripts/THOR_WT_ESCvsNPC_rmdup_TMM.sh # 1674683 ok

### Optimal parameters have been found with TMM-norm; now here are the missing samples:
#### Time-effect for KO and HET, Default TMM-normalization (NO SF)
sbatch scripts/THOR_HET_ESCvsNPC_UniqueBamTMM.sh # 1681090 ok
sbatch scripts/THOR_KO_ESCvsNPC_UniqueBamTMM.sh # 1681115 ok
sbatch scripts/THOR_WT_NPCvs2dN_UniqueBamTMM.sh # 1700276 ok
sbatch scripts/THOR_HET_NPCvs2dN_UniqueBamTMM.sh # 1700332 ok
sbatch scripts/THOR_KO_NPCvs2dN_UniqueBamTMM.sh # 1700357 ok
#### Time-effect for KO and HET, UniqueBamDiffBindTMM norm
sbatch scripts/THOR_HET_ESCvsNPC_UniqueBamDiffBindTMM.sh # 3833376
sbatch scripts/THOR_KO_ESCvsNPC_UniqueBamDiffBindTMM.sh # 3833447
### Time-effect for WT, HET and KO; UniqueBam with reciprocal ChIPseqSpikeInFree (From uniqueBAM) SF
sbatch scripts/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Reciprocal.sh # 3868853
sbatch scripts/THOR_HET_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Reciprocal.sh # 3868879 FAIL; 3869149 FAIL; 3870062
sbatch scripts/THOR_KO_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Reciprocal.sh #  WAIT IF WT GOOD





#### Genotype comparison at NPC and 2dN, Default TMM-normalization (NO SF)
sbatch scripts/THOR_NPC_WTvsHET_UniqueBamTMM.sh # 1681129 ok
sbatch scripts/THOR_NPC_WTvsKO_UniqueBamTMM.sh # 1681138 ok
sbatch scripts/THOR_2dN_WTvsHET_UniqueBamTMM.sh # 1681142 ok
sbatch scripts/THOR_2dN_WTvsKO_UniqueBamTMM.sh # 1681198 ok






## housekeeping genes-norm
#### time-effect
sbatch scripts/THOR_WT_ESCvsNPC_housekeep.sh # 1915845 ok
sbatch scripts/THOR_WT_ESCvsNPC_uniqueBAMhousekeep.sh # 1915944 ok
sbatch scripts/THOR_HET_ESCvsNPC_housekeep.sh # 1916367 ok
sbatch scripts/THOR_HET_ESCvsNPC_uniqueBAMhousekeep.sh # 1916383 ok
sbatch scripts/THOR_KO_ESCvsNPC_housekeep.sh # 1916443 ok
sbatch scripts/THOR_KO_ESCvsNPC_uniqueBAMhousekeep.sh # 1916509 ok

sbatch scripts/THOR_WT_NPCvs2dN_housekeep.sh # 1916999 ok
sbatch scripts/THOR_WT_NPCvs2dN_uniqueBAMhousekeep.sh # 1917008 ok
sbatch scripts/THOR_HET_NPCvs2dN_housekeep.sh # 1917186 ok
sbatch scripts/THOR_HET_NPCvs2dN_uniqueBAMhousekeep.sh # 1917268 ok
sbatch scripts/THOR_KO_NPCvs2dN_housekeep.sh # 1917343 ok
sbatch scripts/THOR_KO_NPCvs2dN_uniqueBAMhousekeep.sh # 1917401 ok
#### genotype-effect
sbatch scripts/THOR_ESC_WTvsHET_housekeep.sh # 1916031 ok
sbatch scripts/THOR_ESC_WTvsHET_uniqueBAMhousekeep.sh # 1916074 ok
sbatch scripts/THOR_ESC_WTvsKO_housekeep.sh # 1917516
sbatch scripts/THOR_ESC_WTvsKO_uniqueBAMhousekeep.sh # 2020601

sbatch scripts/THOR_NPC_WTvsHET_housekeep.sh # 1917668 ok 
sbatch scripts/THOR_NPC_WTvsHET_uniqueBAMhousekeep.sh # 1917721 ok 
sbatch scripts/THOR_NPC_WTvsKO_housekeep.sh # 1917810 ok 
sbatch scripts/THOR_NPC_WTvsKO_uniqueBAMhousekeep.sh # 1917820 ok 

sbatch scripts/THOR_2dN_WTvsHET_housekeep.sh #  1918597 ok 
sbatch scripts/THOR_2dN_WTvsHET_uniqueBAMhousekeep.sh # 1918605 ok 
sbatch scripts/THOR_2dN_WTvsKO_housekeep.sh # 1918654 ok 
sbatch scripts/THOR_2dN_WTvsKO_uniqueBAMhousekeep.sh # 1918827 ok 
```

Go in R to explore the data real quick within `conda activate deseq2`:

```R
# load the file using the tidyverse
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# WT_ESCvsNPC_UniqueBamDiffBindTMM
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/log2FC_qval10.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 10) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()

# WT_ESCvsNPC_DiffBindTMM
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_DiffBindTMM/WTESCvsNPCDiffBindTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_DiffBindTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()


# WT_ESCvsNPC_ChIPseqSpikeInFree
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_ChIPseqSpikeInFree/WTESCvsNPCChIPseqSpikeInFree-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_ChIPseqSpikeInFree/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()

# WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree/WTESCvsNPCUniqueBamChIPseqSpikeInFree-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()


# WT_ESCvsNPC_ChIPseqSpikeInFree_Corr
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_ChIPseqSpikeInFree_Corr/WTESCvsNPCChIPseqSpikeInFreeCorr-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_ChIPseqSpikeInFree_Corr/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()

# WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Corr
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Corr/WTESCvsNPCUniqueBamChIPseqSpikeInFreeCorr-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Corr/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()



# WT_ESCvsNPC_TMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_TMM/WTESCvsNPCTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_TMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_WT_ESCvsNPC_TMM/log2FC_qval10.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 15) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 20) %>%
  write_tsv("output/THOR/THOR_WT_ESCvsNPC_TMM/THOR_qval20.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())

# WT_ESCvsNPC_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/log2FC_qval20.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 20) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 20) %>%
  write_tsv("output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval20.bed", col_names = FALSE)
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval30.bed", col_names = FALSE)

## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())



# HET_ESCvsNPC_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_HET_ESCvsNPC_UniqueBamTMM/HETESCvsNPCUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_HET_ESCvsNPC_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("HET_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_HET_ESCvsNPC_UniqueBamTMM/log2FC_qval30.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("HET_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_HET_ESCvsNPC_UniqueBamTMM/THOR_qval30.bed", col_names = FALSE)


## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())



# KO_ESCvsNPC_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_KO_ESCvsNPC_UniqueBamTMM/KOESCvsNPCUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_KO_ESCvsNPC_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("KO_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_KO_ESCvsNPC_UniqueBamTMM/log2FC_qval30.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 30) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("KO_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_KO_ESCvsNPC_UniqueBamTMM/THOR_qval30.bed", col_names = FALSE)
  

## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())



# ESC_WTvsHET_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_ESC_WTvsHET_UniqueBamTMM/ESCWTvsHETUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_HET", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_HET, into = c("count_HET_1","count_HET_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_HET_1+count_HET_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_ESC_WTvsHET_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs HET") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_ESC_WTvsHET_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs HET") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 5) %>%
  write_tsv("output/THOR/THOR_ESC_WTvsHET_UniqueBamTMM/THOR_qval5.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())


# ESC_WTvsKO_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_ESC_WTvsKO_UniqueBamTMM/ESCWTvsKOUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_ESC_WTvsKO_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_ESC_WTvsKO_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs KO") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 5) %>%
  write_tsv("output/THOR/THOR_ESC_WTvsKO_UniqueBamTMM/THOR_qval5.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())



# NPC_WTvsHET_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_NPC_WTvsHET_UniqueBamTMM/NPCWTvsHETUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_HET", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_HET, into = c("count_HET_1","count_HET_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_HET_1+count_HET_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_WTvsHET_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs HET") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_NPC_WTvsHET_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 5) %>%
  write_tsv("output/THOR/THOR_NPC_WTvsHET_UniqueBamTMM/THOR_qval5.bed", col_names = FALSE)
  ## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())



# NPC_WTvsKO_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_NPC_WTvsKO_UniqueBamTMM/NPCWTvsKOUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_WTvsKO_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_NPC_WTvsKO_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 5) %>%
  write_tsv("output/THOR/THOR_NPC_WTvsKO_UniqueBamTMM/THOR_qval5.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())



# 2dN_WTvsHET_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_2dN_WTvsHET_UniqueBamTMM/2dNWTvsHETUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_HET", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_HET, into = c("count_HET_1","count_HET_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_HET_1+count_HET_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_2dN_WTvsHET_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs HET") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_2dN_WTvsHET_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs HET") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 5) %>%
  write_tsv("output/THOR/THOR_2dN_WTvsHET_UniqueBamTMM/THOR_qval5.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())



# 2dN_WTvsKO_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_2dN_WTvsKO_UniqueBamTMM/2dNWTvsKOUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_2dN_WTvsKO_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_2dN_WTvsKO_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs KO") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 5) %>%
  write_tsv("output/THOR/THOR_2dN_WTvsKO_UniqueBamTMM/THOR_qval5.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())

# WT_NPCvs2dN_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_WT_NPCvs2dN_UniqueBamTMM/WTNPCvs2dNUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_NPC", "count_2dN", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  separate(count_2dN, into = c("count_2dN_1","count_2dN_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_2dN_1+count_2dN_2) / (count_NPC_1+count_NPC_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_NPCvs2dN_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_NPC vs 2dN") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_WT_NPCvs2dN_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_NPC vs 2dN") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 20) %>%
  write_tsv("output/THOR/THOR_WT_NPCvs2dN_UniqueBamTMM/THOR_qval20.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())

# 2dN_WTvsKO_UniqueBamTMM
diffpeaks <- read_tsv("output/THOR/THOR_2dN_WTvsKO_UniqueBamTMM/2dNWTvsKOUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
  
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_2dN_WTvsKO_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_2dN_WTvsKO_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs KO") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 5) %>%
  write_tsv("output/THOR/THOR_2dN_WTvsKO_UniqueBamTMM/THOR_qval5.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())



# HET_NPCvs2dN_UniqueBamTMM (Defult TMM norm)
diffpeaks <- read_tsv("output/THOR/THOR_HET_NPCvs2dN_UniqueBamTMM/HETNPCvs2dNUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_NPC", "count_2dN", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  separate(count_2dN, into = c("count_2dN_1","count_2dN_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_2dN_1+count_2dN_2) / (count_NPC_1+count_NPC_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_HET_NPCvs2dN_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("HET_NPC vs 2dN") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_HET_NPCvs2dN_UniqueBamTMM/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("HET_NPC vs 2dN") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_HET_NPCvs2dN_UniqueBamTMM/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 25) %>%
  group_by(X6) %>%
  summarise(n = n())


# KO_NPCvs2dN_UniqueBamTMM
diffpeaks <- read_tsv("output/THOR/THOR_KO_NPCvs2dN_UniqueBamTMM/KONPCvs2dNUniqueBamTMM-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_NPC", "count_2dN", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  separate(count_2dN, into = c("count_2dN_1","count_2dN_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_2dN_1+count_2dN_2) / (count_NPC_1+count_NPC_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_KO_NPCvs2dN_UniqueBamTMM/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("KO_NPC vs 2dN") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_KO_NPCvs2dN_UniqueBamTMM/log2FC_qval15.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 15) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("KO_NPC vs 2dN") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 30) %>%
  write_tsv("output/THOR/THOR_KO_NPCvs2dN_UniqueBamTMM/THOR_qval30.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 10) %>%
  group_by(X6) %>%
  summarise(n = n())


# WT_ESCvsNPC_housekeep
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_housekeep/WTESCvsNPChousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_housekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_WT_ESCvsNPC_housekeep/log2FC_qval40.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 40) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 60) %>%
  write_tsv("output/THOR/THOR_WT_ESCvsNPC_housekeep/THOR_qval60.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 40) %>%
  group_by(X6) %>%
  summarise(n = n())

# WT_ESCvsNPC_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_WT_ESCvsNPC_uniqueBAMhousekeep/WTESCvsNPCuniqueBAMhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_ESCvsNPC_uniqueBAMhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_WT_ESCvsNPC_uniqueBAMhousekeep/log2FC_qval60.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 60) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_WT_ESCvsNPC_uniqueBAMhousekeep/THOR_qval60.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())


# HET_ESCvsNPC_housekeep
diffpeaks <- read_tsv("output/THOR/THOR_HET_ESCvsNPC_housekeep/HETESCvsNPChousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_HET_ESCvsNPC_housekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_HET_ESCvsNPC_housekeep/log2FC_qval50.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("HET_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 100) %>%
  write_tsv("output/THOR/THOR_HET_ESCvsNPC_housekeep/THOR_qval100.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())


# HET_ESCvsNPC_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_HET_ESCvsNPC_UniqueBamhousekeep/HETESCvsNPCUniqueBamhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_HET_ESCvsNPC_UniqueBamhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_HET_ESCvsNPC_UniqueBamhousekeep/log2FC_qval50.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("HET_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 40) %>%
  write_tsv("output/THOR/THOR_HET_ESCvsNPC_UniqueBamhousekeep/THOR_qval40.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())


# KO_ESCvsNPC_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_KO_ESCvsNPC_UniqueBamhousekeep/KOESCvsNPCUniqueBamhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_ESC", "count_NPC", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_ESC, into = c("count_ESC_1","count_ESC_2"), sep = ":", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_NPC_1+count_NPC_2) / (count_ESC_1+count_ESC_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_KO_ESCvsNPC_UniqueBamhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("KO_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_KO_ESCvsNPC_UniqueBamhousekeep/log2FC_qval50.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("KO_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 60) %>%
  write_tsv("output/THOR/THOR_KO_ESCvsNPC_UniqueBamhousekeep/THOR_qval60.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())


# WT_NPCvs2dN_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_WT_NPCvs2dN_UniqueBamhousekeep/WTNPCvs2dNUniqueBamhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_NPC", "count_2dN", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_NPC, into = c("count_NPC_1","count_NPC_2"), sep = ":", convert = TRUE) %>%
  separate(count_2dN, into = c("count_2dN_1","count_2dN_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_2dN_1+count_2dN_2) / (count_NPC_1+count_NPC_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_WT_NPCvs2dN_UniqueBamhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_WT_NPCvs2dN_UniqueBamhousekeep/log2FC_qval60.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 60) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 60) %>%
  write_tsv("output/THOR/THOR_WT_NPCvs2dN_UniqueBamhousekeep/THOR_qval60.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())


# ESC_WTvsHET_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_ESC_WTvsHET_uniqueBAMhousekeep/ESCWTvsHETuniqueBAMhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_HET", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_HET, into = c("count_HET_1","count_HET_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_HET_1+count_HET_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_ESC_WTvsHET_uniqueBAMhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs HET") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_ESC_WTvsHET_uniqueBAMhousekeep/log2FC_qval60.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 60) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 40) %>%
  write_tsv("output/THOR/THOR_ESC_WTvsHET_uniqueBAMhousekeep/THOR_qval40.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())


# ESC_WTvsKO_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_ESC_WTvsKO_UniqueBamhousekeep/ESCWTvsKOUniqueBamhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_HET", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_HET, into = c("count_HET_1","count_HET_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_HET_1+count_HET_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_ESC_WTvsHET_uniqueBAMhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("ESC_WT vs HET") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_ESC_WTvsHET_uniqueBAMhousekeep/log2FC_qval60.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 60) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("WT_ESC vs NPC") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 40) %>%
  write_tsv("output/THOR/THOR_ESC_WTvsHET_uniqueBAMhousekeep/THOR_qval40.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())


# NPC_WTvsHET_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_NPC_WTvsHET_UniqueBamhousekeep/NPCWTvsHETUniqueBamhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_HET", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_HET, into = c("count_HET_1","count_HET_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_HET_1+count_HET_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_WTvsHET_UniqueBamhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs HET") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_NPC_WTvsHET_UniqueBamhousekeep/log2FC_qval50.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs HET") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 25) %>%
  write_tsv("output/THOR/THOR_NPC_WTvsHET_UniqueBamhousekeep/THOR_qval25.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())



# NPC_WTvsKO_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_NPC_WTvsKO_UniqueBamhousekeep/NPCWTvsKOUniqueBamhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_NPC_WTvsKO_UniqueBamhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_NPC_WTvsKO_UniqueBamhousekeep/log2FC_qval50.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("NPC_WT vs HET") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_NPC_WTvsKO_UniqueBamhousekeep/THOR_qval50.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())




# 2dN_WTvsHET_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_2dN_WTvsHET_UniqueBamhousekeep/2dNWTvsHETUniqueBamhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_HET", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_HET, into = c("count_HET_1","count_HET_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_HET_1+count_HET_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_2dN_WTvsHET_UniqueBamhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs HET") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_2dN_WTvsHET_UniqueBamhousekeep/log2FC_qval50.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 50) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs HET") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_2dN_WTvsHET_UniqueBamhousekeep/THOR_qval50.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())




# 2dN_WTvsKO_uniqueBAMhousekeep
diffpeaks <- read_tsv("output/THOR/THOR_2dN_WTvsKO_UniqueBamhousekeep/2dNWTvsKOUniqueBamhousekeep-diffpeaks.bed",
                      col_names = FALSE, trim_ws = TRUE, col_types = cols(X1 = col_character()))
## split the last field and calculate FC
thor_splitted = diffpeaks %>%
  separate(X11, into = c("count_WT", "count_KO", "qval"), sep = ";", convert = TRUE) %>%
  separate(count_WT, into = c("count_WT_1","count_WT_2"), sep = ":", convert = TRUE) %>%
  separate(count_KO, into = c("count_KO_1","count_KO_2"), sep = ":", convert = TRUE) %>%
  mutate(FC = (count_KO_1+count_KO_2) / (count_WT_1+count_WT_2))
## plot the histogram of the fold-change computed above, count second condition / count 1st condition
pdf("output/THOR/THOR_2dN_WTvsKO_UniqueBamhousekeep/log2FC.pdf", width=14, height=14)
thor_splitted %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs KO") +
  theme_bw()
dev.off()
pdf("output/THOR/THOR_2dN_WTvsKO_UniqueBamhousekeep/log2FC_qval25.pdf", width=14, height=14)
thor_splitted %>%
  filter(qval > 25) %>%
  ggplot(aes(x = log2(FC))) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(-5, 3, 1)) +
  ggtitle("2dN_WT vs KO") +
  theme_bw()
dev.off()
## create a bed file, append chr to chromosome names and write down the file
thor_splitted %>%
  filter(qval > 50) %>%
  write_tsv("output/THOR/THOR_2dN_WTvsKO_UniqueBamhousekeep/THOR_qval50.bed", col_names = FALSE)
## how many minus / plus
thor_splitted %>%
  filter(qval > 50) %>%
  group_by(X6) %>%
  summarise(n = n())


```

--> `THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM` and `THOR_WT_ESCvsNPC_DiffBindTMM` = bad; much more Decrease H3K27me3 upon differentation...
----> Genotype comparison I do not even look as differentiation control increase H3K27me3 fail...

--> `THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree` and `THOR_WT_ESCvsNPC_ChIPseqSpikeInFree` = bad; show only decrease H3K27me3 upon differentiation... THIS IS WEIRD... As scaling factor used should show increase!! 
----> Mistake is that I directly put the output SF value from ChIPseqSpikeInFree, but should have normalize it for library-size before; with:   `1 / (15000000/($libSize*$SF))`

--> Analysis with **UniqueBamTMM** (no SF applied):
----> **time**-effect; WT similar amount gain/lost from ESC to NPC; however for HET and KO, most are lost H3K27me3 upon diff
----> **genotype**-effect; at ESC/NPC/2dN HET and KO lead to overall decrease of H3K27me3 (this was NOT observed with the CutRun data...)
------> WT vs KO at 2dN failed; almost no diff. sites

**Housekeeping gene-norm**:

--> More diff bound sites identified with non-uniqueBAM; BUT **replicates are MORE clean with uniqueBAM**

--> The 2 replicates of **HET ESC are very heterogeneous** (uniqueBam and non-UniqueBam): Something needs to be corected; check the housekeeping genes for this sampple : XXX

--> The 2 replicates of **KO ESC and NPC are a bit heterogeneous**

--> The 2 replicates of **WT 2dN are a bit heterogeneous**

----> shows decrease H3K27me3 from ESC to NPC at NEUROG2.

----> Much more diff. bound sites when using the non-uniqueBAM one!

----> Let's pick optimal qvalue for WT and apply fotr the rest:
- time-effect: qval 50 (good for unique and non-uniqueBAM)
- genotype-effect: qval XXX



Let's collect the correct SF from ChIPseqSpikeInFree:

**ChIPseqSpikeInFree SF_CORRECT** so not pass-by DiffBind, no TMM-normalized but manually LIB-size normalized (used to generate `output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF`)
- sample / libSize_UNIQUE @ libSize_NotUnique | SF(From Unique ChIPseqSpikeInFree) | 1/(15000000/($libSize*$SF))_UniqueLibsize ? 1/(15000000/($libSize*$SF))_NotUniqueLibsize
- 2dN_HET_H3K27me3_R1 / 40818892 @ 70887858 | 1.97 | 5.360881149 ? 9.309938684
- 2dN_HET_H3K27me3_R2 / 38757156 @ 67346082 | 1.75 | 4.5216682	? 7.8570429
- 2dN_KO_H3K27me3_R1 / 36478530 @ 58363794 | 1.46 | 3.55057692 ? 5.680742616
- 2dN_KO_H3K27me3_R2 / 46363646 @ 71964648 | 1 | 3.090909733 ? 4.7976432
- 2dN_WT_H3K27me3_R1 / 43088276 @ 71682964 | 1.29 | 3.705591736 ? 6.164734904
- 2dN_WT_H3K27me3_R2 / 43478292 @ 60792560 | 1.69 | 4.898554232 ? 6.849295093
- ESC_HET_H3K27me3_R1 / 44064716 @ 85933694 | 10.51 | 30.87467768 ? 60.21087493
- ESC_HET_H3K27me3_R2 / 29954726 @ 66583922 | 23.35 | 46.62952347 ? 103.6489719
- ESC_KO_H3K27me3_R1 / 27538348 @ 86014942 | 10.06 | 18.46905206 ? 57.68735443
- ESC_KO_H3K27me3_R2 / 23084678 @ 57643920 | 15.78 | 24.28508126 ?	60.64140384
- ESC_WT_H3K27me3_R1 / 57347074 @ 90968406 | 7 | 26.76196787 ? 42.4519228
- ESC_WT_H3K27me3_R2 / 58089688 @ 79650052 | 4.31 | 16.69110369 ? 22.88611494
- NPC_HET_H3K27me3_R1 / 28967884 @ 40423510 | 1.13 | 2.182247261 ?	3.045237753
- NPC_HET_H3K27me3_R2 / 53067608 @ 82710030 | 1.43 | 5.059111963 ? 7.88502286
- NPC_KO_H3K27me3_R1 / 42573738 @ 67904376 | 1.55 | 4.39928626 ? 7.01678552 
- NPC_KO_H3K27me3_R2 / 48730202 @ 84619268 | 2.64 | 8.576515552 ? 14.89299117
- NPC_WT_H3K27me3_R1 / 52824596 @ 77698696 | 1.45 | 5.106377613 ? 7.510873947
- NPC_WT_H3K27me3_R2 / 51005778 @ 72126718 | 1.51 | 5.134581652 ? 7.260756279


***NOTE: ChIPseqSpikeInFree bigwig files have NOT been generated using BamCoverage, but with  `genomeCoverageBed -bg -scale $scale ` where the `$scale` is a multiplying factor; thus I did: `1/15000000/($libSize*$SF)` --> Otherwise the non-reciprocal SF return will increase signal for ESC...*** --> File for calculation in GoogleDrive under `002__ChIPseq/ChIPseqSpikeInFactor_ScalingFactor.xlsx` and new THOR-related files as `*_Corr`

--> Using LIB-normalize ChIPseqSpikeInFree scaling factors is similarly weird, we still have only less H3K27me3 in NPC than ESC...

--> Using TMM-Default we got comparable number of Gain/Lost from ESC to NPC in WT

--> qval 10-15 for WT_ESCvsNPC is good for both uniqueBam and not uniqueBam. However uniqueBam shows more H3K27me3 in ESC and non-uniqueBam show more H3K27me3 in NPC...

--> uniqueBam (using pre-filtered bam files that only contain uniquely mapepd reads) versus using `--rmdup` argument in THOR is better; as it work LOL... `--rmdup` fail and did not give any diff bound peaks



### Assign THOR-diff peaks to genes and check expression

Now let's compare RNAseq (expression) and CutRun for THOR qval 15 among others:
- Filter HETvsWT and KOvsWT diff bound genes into **gain and loss H3K27me3**
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


# Import diff. peaks
## TIME-EFFECT ESC_NPC
## qval10_WT_ESCvsNPC_TMM
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_TMM/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval10_WT_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval15_WT_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval15.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval20_WT_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval20.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval25_WT_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval30_WT_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval30.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval20_WT_ESCvsNPC_TMM
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_TMM/THOR_qval20.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval40_WT_ESCvsNPC_housekeep
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_housekeep/THOR_qval40.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval60_WT_ESCvsNPC_uniqueBAMhousekeep
ESCvsNPC = read.table('output/THOR/THOR_WT_ESCvsNPC_uniqueBAMhousekeep/THOR_qval60.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval50_HET_ESCvsNPC_uniqueBAMhousekeep
ESCvsNPC = read.table('output/THOR/THOR_HET_ESCvsNPC_UniqueBamhousekeep/THOR_qval50.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval50_KO_ESCvsNPC_uniqueBAMhousekeep
ESCvsNPC = read.table('output/THOR/THOR_KO_ESCvsNPC_UniqueBamhousekeep/THOR_qval50.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)



## qval10_HET_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_HET_ESCvsNPC_UniqueBamTMM/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval25_HET_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_HET_ESCvsNPC_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval10_KO_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_KO_ESCvsNPC_UniqueBamTMM/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval25_KO_ESCvsNPC_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_KO_ESCvsNPC_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)

## TIME-EFFECT NPC_2dN
## qval25_WT_NPCvs2dN_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_WT_NPCvs2dN_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)
## qval25_HET_NPCvs2dN_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_HET_NPCvs2dN_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_ESC_1= V11, count_ESC_2=V12, count_NPC_1=V13, count_NPC_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_ESC_1,count_ESC_2,count_NPC_1,count_NPC_2)




## GENOTYPE-EFFECT : !!! TO MAKE IT SIMPLE, I KEPT the 'ESCvsNPC' title !!!!
## qval5_ESC_WTvsHET_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_ESC_WTvsHET_UniqueBamTMM/THOR_qval5.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval25_ESC_WTvsHET_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_ESC_WTvsHET_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval25_ESC_WTvsKO_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_ESC_WTvsKO_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval10_ESC_WTvsKO_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_ESC_WTvsKO_UniqueBamTMM/THOR_qval10.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval5_NPC_WTvsHET_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_NPC_WTvsHET_UniqueBamTMM/THOR_qval5.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval25_NPC_WTvsHET_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_NPC_WTvsHET_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval25_NPC_WTvsKO_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_NPC_WTvsKO_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval25_2dN_WTvsHET_UniqueBamTMM
ESCvsNPC = read.table('output/THOR/THOR_2dN_WTvsHET_UniqueBamTMM/THOR_qval25.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval50_NPC_WTvsHET_UniqueBamhousekeep
ESCvsNPC = read.table('output/THOR/THOR_NPC_WTvsHET_UniqueBamhousekeep/THOR_qval50.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval50_NPC_WTvsKO_UniqueBamhousekeep
ESCvsNPC = read.table('output/THOR/THOR_NPC_WTvsKO_UniqueBamhousekeep/THOR_qval50.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval50_2dN_WTvsHET_UniqueBamhousekeep
ESCvsNPC = read.table('output/THOR/THOR_2dN_WTvsHET_UniqueBamhousekeep/THOR_qval50.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)
## qval50_2dN_WTvsKO_UniqueBamhousekeep
ESCvsNPC = read.table('output/THOR/THOR_2dN_WTvsKO_UniqueBamhousekeep/THOR_qval50.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, V7=V7, V8=V8, qvalue=V15, FC=V16, count_WT_1= V11, count_WT_2=V12, count_HET_1=V13, count_HET_2=V14) %>% dplyr::select(Chr, start,end,qvalue,FC,count_WT_1,count_WT_2,count_HET_1,count_HET_2)


# Tidy peaks #-->> Re-Run from here with different qvalue!!
ESCvsNPC_gr = makeGRangesFromDataFrame(ESCvsNPC,keep.extra.columns=TRUE)
gr_list <- list(ESCvsNPC=ESCvsNPC_gr)



# Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here
## Get annotation data frame
ESCvsNPC_annot <- as.data.frame(peakAnnoList[["ESCvsNPC"]]@anno)



## Convert entrez gene IDs to gene symbols
ESCvsNPC_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = ESCvsNPC_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")

ESCvsNPC_annot$gene <- mapIds(org.Hs.eg.db, keys = ESCvsNPC_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
### TIME-EFFECT
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_WT_ESCvsNPC_qval10_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_HET_ESCvsNPC_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F)  # CHANGE TITLE
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_KO_ESCvsNPC_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F) 
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_WT_NPCvs2dN_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F) 
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_HET_NPCvs2dN_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F) 
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_WT_ESCvsNPC_qval40_uniqueBAMhousekeep.txt", sep="\t", quote=F, row.names=F) 

### GENOTYPE-EFFECT
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_ESC_WTvsHET_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F) 
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_ESC_WTvsKO_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F) 
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_NPC_WTvsHET_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F) 
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_NPC_WTvsKO_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F) 
write.table(ESCvsNPC_annot, file="output/ChIPseeker/annotation_2dN_WTvsHET_qval25_UniqueBamTMM.txt", sep="\t", quote=F, row.names=F) 



# Filter Gain/Loss sites
## KEEP Distal Intergenic (keep ALL)   ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
WTvsHET_annot_gain = tibble(WTvsHET_annot) %>%
    filter(FC > 1.5) %>%
    add_column(H3K27me3 = "gain")
WTvsHET_annot_lost = tibble(WTvsHET_annot) %>%
    filter(FC < (1/1.5)) %>%
    add_column(H3K27me3 = "lost")
WTvsHET_annot_gain_lost = WTvsHET_annot_gain %>% 
    bind_rows(WTvsHET_annot_lost) 
    


## Remove Distal Intergenic   ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
ESCvsNPC_annot_gain = tibble(ESCvsNPC_annot) %>%
    filter(FC > 1, annotation != "Distal Intergenic") %>%
    add_column(H3K27me3 = "gain")
ESCvsNPC_annot_lost = tibble(ESCvsNPC_annot) %>%
    filter(FC < (1/1), annotation != "Distal Intergenic") %>%
    add_column(H3K27me3 = "lost")
ESCvsNPC_annot_gain_lost = ESCvsNPC_annot_gain %>% 
    bind_rows(ESCvsNPC_annot_lost) 
    

## Keep only signals in promoter of 5'UTR ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
ESCvsNPC_annot_gain = tibble(ESCvsNPC_annot) %>%
    filter(FC > 1, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "gain")
ESCvsNPC_annot_lost = tibble(ESCvsNPC_annot) %>%
    filter(FC < (1/1), annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "lost")
ESCvsNPC_annot_gain_lost = ESCvsNPC_annot_gain %>% 
    bind_rows(ESCvsNPC_annot_lost) 

ESCvsNPC_annot_gain = tibble(ESCvsNPC_annot) %>%
    filter(FC > 8, annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "gain")
ESCvsNPC_annot_lost = tibble(ESCvsNPC_annot) %>%
    filter(FC < (1/8), annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR")) %>%
    add_column(H3K27me3 = "lost")
ESCvsNPC_annot_gain_lost = ESCvsNPC_annot_gain %>% 
    bind_rows(ESCvsNPC_annot_lost) 


# Import RNAseq deseq2 output
## Raw FC ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_WT_ESC_vs_WT_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_HET_ESC_vs_HET_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)

RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_KO_ESC_vs_KO_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj)
## Fitlered FC ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
#### TIME-EFFECT
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_WT_ESC_vs_WT_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 1 | log2FoldChange <= -1)

RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_WT_ESC_vs_WT_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_HET_ESC_vs_HET_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_KO_ESC_vs_KO_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_WT_NPC_vs_WT_2dN.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)   
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_HET_NPC_vs_HET_2dN.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)   
#### GENOTYPE-EFFECT
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_ESC_HET_vs_ESC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)    
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_ESC_KO_vs_ESC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_HET_vs_NPC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_KO_vs_NPC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_2dN_HET_vs_2dN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)
RNA_expression = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_2dN_KO_vs_2dN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene, baseMean,log2FoldChange,padj) %>%
    filter(log2FoldChange >= 0.5 | log2FoldChange <= -0.5)


# Merge files
ESCvsNPC_annot_gain_lost_RNA = ESCvsNPC_annot_gain_lost %>% 
    left_join(RNA_expression) %>%
    dplyr::select(gene, H3K27me3,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()

ESCvsNPC_annot_gain_lost_RNA = ESCvsNPC_annot_gain_lost %>% 
    left_join(RNA_expression) %>%
    dplyr::select(gene, H3K27me3,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.001) %>%  # add signif TRUE if 0.05
    unique()
# Volcano plot
count_data <- ESCvsNPC_annot_gain_lost_RNA %>%
    group_by(H3K27me3, significance) %>%
    summarise(up = sum(log2FoldChange > 0),
              down = sum(log2FoldChange < 0),
              total = n()) %>%
    ungroup() %>%
    group_by(H3K27me3) %>%
    mutate(total_panel = sum(total)) %>%
    ungroup()


## TIME-EFFECT
pdf("output/ChIPseeker/THOR_qval10_WT_NPCvsESC_TMM_expression.pdf", width=7, height=4)  # CHANGE TITLE
pdf("output/ChIPseeker/THOR_qval10_WT_NPCvsESC_UniqueBamTMM_expression.pdf", width=7, height=4)  # CHANGE TITLE
pdf("output/ChIPseeker/THOR_qval10_WT_NPCvsESC_UniqueBamTMM_expression_FC05.pdf", width=7, height=4)  # CHANGE TITLE
pdf("output/ChIPseeker/THOR_qval10_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval10_FC15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval10_FC2_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval10_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC1.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval10_FC4_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval10_WT_NPCvsESC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval20_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval20_WT_NPCvsESC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 

pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_WT_NPCvsESC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 

pdf("output/ChIPseeker/THOR_qval30_WT_NPCvsESC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 

pdf("output/ChIPseeker/THOR_qval25_HET_NPCvsESC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_KO_NPCvsESC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_WT_2dNvsNPC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_HET_2dNvsNPC_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval40_WT_ESCvsNPC_uniqueBAMhousekeep_expression.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval60_WT_ESCvsNPC_uniqueBAMhousekeep_expression_promoterAnd5.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval50_HET_ESCvsNPC_uniqueBAMhousekeep_expression_promoterAnd5.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval50_KO_ESCvsNPC_uniqueBAMhousekeep_expression_promoterAnd5.pdf", width=7, height=4) # CHANGE TITLE 


pdf("output/ChIPseeker/THOR_qval20_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_FC2_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_FC4_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_FC8_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC1.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_qval0.01.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_qval0.001.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_qval0.0001.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_qval0.00001.pdf", width=7, height=4)  # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval15_WT_NPCvsESC_UniqueBamTMM_expression_promoterAnd5_FC1_qval0.001.pdf", width=7, height=4)  # CHANGE TITLE 



## GENOTYPE-EFFECT
pdf("output/ChIPseeker/THOR_qval25_ESC_WTvsHET_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_ESC_WTvsKO_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_NPC_WTvsHET_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_NPC_WTvsKO_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval25_2dN_WTvsHET_TMM_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval50_NPC_WTvsHET_housekeep_expression_promoterAnd5.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval50_NPC_WTvsKO_housekeep_expression_promoterAnd5.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval50_2dN_WTvsHET_housekeep_expression_promoterAnd5.pdf", width=7, height=4) # CHANGE TITLE 
pdf("output/ChIPseeker/THOR_qval50_2dN_WTvsKO_housekeep_expression_promoterAnd5.pdf", width=7, height=4) # CHANGE TITLE 



ESCvsNPC_annot_gain_lost_RNA %>%
    ggplot(aes(x = log2FoldChange, y = baseMean, color = significance)) +
        geom_point(alpha = 0.8, size = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        labs(title = "NPC vs ESC",
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



# Remove the genes that both gained and lost H3K27me3
## Identify genes that are present in both gain and lost categories
common_genes_HET <- intersect(WTvsHET_annot_gain$gene, WTvsHET_annot_lost$gene)
common_genes_KO <- intersect(WTvsKO_annot_gain$gene, WTvsKO_annot_lost$gene)

## Remove these genes from your gain and lost data frames
WTvsHET_annot_gain <- WTvsHET_annot_gain %>% filter(!(gene %in% common_genes_HET))
WTvsHET_annot_lost <- WTvsHET_annot_lost %>% filter(!(gene %in% common_genes_HET))

WTvsKO_annot_gain <- WTvsKO_annot_gain %>% filter(!(gene %in% common_genes_KO))
WTvsKO_annot_lost <- WTvsKO_annot_lost %>% filter(!(gene %in% common_genes_KO))

## Now bind the rows
WTvsHET_annot_gain_lost = WTvsHET_annot_gain %>% bind_rows(WTvsHET_annot_lost)
WTvsKO_annot_gain_lost = WTvsKO_annot_gain %>% bind_rows(WTvsKO_annot_lost)


# Merge files with RNA
WTvsHET_annot_gain_lost_RNA = WTvsHET_annot_gain_lost %>% 
    left_join(HET_vs_WT) %>%
    dplyr::select(gene, H3K27me3,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()


WTvsKO_annot_gain_lost_RNA = WTvsKO_annot_gain_lost %>% 
    left_join(KO_vs_WT) %>%
    dplyr::select(gene, H3K27me3,baseMean,log2FoldChange,padj) %>%
    filter(gene != "NA") %>%
    mutate(baseMean = replace_na(baseMean, 0),
           log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05) %>%  # add signif TRUE if 0.05
    unique()






# Volcano plot
count_data <- WTvsHET_annot_gain_lost_RNA %>%
    group_by(H3K27me3, significance) %>%
    summarise(up = sum(log2FoldChange > 0),
              down = sum(log2FoldChange < 0),
              total = n()) %>%
    ungroup() %>%
    group_by(H3K27me3) %>%
    mutate(total_panel = sum(total)) %>%
    ungroup()

pdf("output/ChIPseeker/THOR_qval5_HETvsWTunique_expression.pdf", width=7, height=4)  # CHANGE TITLE
pdf("output/ChIPseeker/THOR_qval15_HETvsWTunique_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE
pdf("output/ChIPseeker/THOR_qval10_HETvsWTunique_expression_promoterAnd5_FC05.pdf", width=7, height=4)  # CHANGE TITLE



WTvsHET_annot_gain_lost_RNA %>%
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

count_data <- WTvsKO_annot_gain_lost_RNA %>%
    group_by(H3K27me3, significance) %>%
    summarise(up = sum(log2FoldChange > 0),
              down = sum(log2FoldChange < 0),
              total = n()) %>%
    ungroup() %>%
    group_by(H3K27me3) %>%
    mutate(total_panel = sum(total)) %>%
    ungroup()

pdf("output/ChIPseeker/THOR_qval5_KOvsWTunique_expression.pdf", width=7, height=4) # CHANGE TITLE
pdf("output/ChIPseeker/THOR_qval5_KOvsWTunique_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE
pdf("output/ChIPseeker/THOR_qval10_KOvsWTunique_expression_promoterAnd5_FC05.pdf", width=7, height=4) # CHANGE TITLE


WTvsKO_annot_gain_lost_RNA %>%
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

**Elucidate what is the best paramter to use:**
--> Here for **RNaseq log2FC positive = more express in NPC**; log2FC negative = LESS express in NPC vs ESC. For ChIPseq **gain = H3K27me3 GAIN/incerase in NPC**

--> Uniquely aligned reads is MUCH better! 

--> Optimal parameter limited nb of false positive is: **uniquely aligned reads (uniqueBam), THOR qvalue25/30, expression log2FC 0.5, peak in promoter and 5' region**

**One-by-one comparison; Genotype and Time -effect**
--> Time-Effect, ESC to NPC; HET and KO; most sites lose H3K27me3!; very few gain! (in agreement with gene expresion...)
----> Thus, not clear how the genes are downregulated upon ESC to NPC in mutants... Seems NOT by H3K27me3, as compare to the WT where that is indeed the case!!!


--> Genotype-effect; probably MANY intergenic region affected; with LOST of H3K27me3 for HET and KO! Because, when filtering and keeping only the peak within the genes; we got comparable gain/lost genes affected


**Hypothesis**: HET mature faster than WT. To check this, check whether:
- Genes that gained H3K27me3 upon ESC to NPC are more closely related for WT / HET then WT / KO
- Genes that gained H3K27me3 in HET at ESC are the one that gained H3K27em3 in WT from ESC to NPC
- Export gene list
```bash
# Files
output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt
output/ChIPseeker/annotation_HET_ESCvsNPC_qval25_UniqueBamTMM.txt
output/ChIPseeker/annotation_KO_ESCvsNPC_qval25_UniqueBamTMM.txt

# Filter the gene that gain / lost H3K27me3 (signal in promoter or 5' only)
## Gain in WT
awk -F'\t' '($7 > 1) && ($12=="Promoter (<=1kb)" || $12=="Promoter (1-2kb)" || $12=="Promoter (2-3kb)" || $12=="5'\'' UTR")' output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt | cut -d $'\t' -f 21 | sort | uniq > output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Gain.txt
## Lose in WT
awk -F'\t' '($7 < 1) && ($12=="Promoter (<=1kb)" || $12=="Promoter (1-2kb)" || $12=="Promoter (2-3kb)" || $12=="5'\'' UTR")' output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM.txt | cut -d $'\t' -f 21 | sort | uniq > output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Lost.txt
## Gain in HET
awk -F'\t' '($7 > 1) && ($12=="Promoter (<=1kb)" || $12=="Promoter (1-2kb)" || $12=="Promoter (2-3kb)" || $12=="5'\'' UTR")' output/ChIPseeker/annotation_HET_ESCvsNPC_qval25_UniqueBamTMM.txt | cut -d $'\t' -f 21 | sort | uniq > output/ChIPseeker/annotation_HET_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Gain.txt
## Lose in HET
awk -F'\t' '($7 < 1) && ($12=="Promoter (<=1kb)" || $12=="Promoter (1-2kb)" || $12=="Promoter (2-3kb)" || $12=="5'\'' UTR")' output/ChIPseeker/annotation_HET_ESCvsNPC_qval25_UniqueBamTMM.txt | cut -d $'\t' -f 21 | sort | uniq > output/ChIPseeker/annotation_HET_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Lost.txt
## Gain in KO
awk -F'\t' '($7 > 1) && ($12=="Promoter (<=1kb)" || $12=="Promoter (1-2kb)" || $12=="Promoter (2-3kb)" || $12=="5'\'' UTR")' output/ChIPseeker/annotation_KO_ESCvsNPC_qval25_UniqueBamTMM.txt | cut -d $'\t' -f 21 | sort | uniq > output/ChIPseeker/annotation_KO_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Gain.txt
## Lose in KO
awk -F'\t' '($7 < 1) && ($12=="Promoter (<=1kb)" || $12=="Promoter (1-2kb)" || $12=="Promoter (2-3kb)" || $12=="5'\'' UTR")' output/ChIPseeker/annotation_KO_ESCvsNPC_qval25_UniqueBamTMM.txt | cut -d $'\t' -f 21 | sort | uniq > output/ChIPseeker/annotation_KO_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Lost.txt
## Gain in HET (at ESC; WT vs HET)
awk -F'\t' '($7 > 1) && ($12=="Promoter (<=1kb)" || $12=="Promoter (1-2kb)" || $12=="Promoter (2-3kb)" || $12=="5'\'' UTR")' output/ChIPseeker/annotation_ESC_WTvsHET_qval25_UniqueBamTMM.txt | cut -d $'\t' -f 21 | sort | uniq > output/ChIPseeker/annotation_ESC_WTvsHET_qval25_UniqueBamTMM_geneSymbol_Gain.txt
## Gain in KO (at ESC; WT vs KO)
awk -F'\t' '($7 > 1) && ($12=="Promoter (<=1kb)" || $12=="Promoter (1-2kb)" || $12=="Promoter (2-3kb)" || $12=="5'\'' UTR")' output/ChIPseeker/annotation_ESC_WTvsKO_qval25_UniqueBamTMM.txt | cut -d $'\t' -f 21 | sort | uniq > output/ChIPseeker/annotation_ESC_WTvsKO_qval25_UniqueBamTMM_geneSymbol_Gain.txt
```
- Venn diagram on [webtool](https://www.biovenn.nl/index.php)

--> Hard to conclude anything. More overlap with HET, but expected as more genes that LOST H3K27me3 in HET than KO...


Let's do **functional analysies on the genes that specifically (or not) LOST H3K27me3** in WT; and in HET; and in KO. Check their expression next:

```bash
conda activate deseq2
```

```R
# library
library("rtracklayer")
library("tidyverse")
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("DOSE")
library("pathview")
library("enrichplot")


# Import the gene list
## Speciifc LOST
gene_symbols <- readLines("output/ChIPseeker/ESCvsNPC_WTspecific_Lost.txt") # readLines transform line into vector
gene_symbols <- readLines("output/ChIPseeker/ESCvsNPC_HETspecific_Lost.txt") 
gene_symbols <- readLines("output/ChIPseeker/ESCvsNPC_KOspecific_Lost.txt") 
## WT changes TMM norm
gene_symbols <- readLines("output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Gain.txt") 
gene_symbols <- readLines("output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Lost.txt") 

genes <- mapIds(org.Hs.eg.db, 
                   keys = gene_symbols, 
                   column = "ENTREZID", 
                   keytype = "SYMBOL", 
                   multiVals = "first")



# Functional profiles
## KEGG
enrichKEGG <- enrichKEGG(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/ChIPseeker/functional_KEGG_ESCvsNPC_HETspecific_Lost.pdf", width=7, height=5)
pdf("output/ChIPseeker/functional_KEGG_ESCvsNPC_WT_Gain.pdf", width=7, height=6)
pdf("output/ChIPseeker/functional_KEGG_ESCvsNPC_WT_Lost.pdf", width=7, height=6)

dotplot(pairwise_termsim(enrichKEGG), showCategory = 15)
dev.off()
### NO ENRICHMENT ESCvsNPC_WTspecific_Lost
### NO ENRICHMENT ESCvsNPC_KOspecific_Lost 
### ENRICHMENT for HET; enuron activity


## GO
enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)

pdf("output/ChIPseeker/functional_GO_BP_ESCvsNPC_HETspecific_Lost_15.pdf", width=7, height=6)
pdf("output/ChIPseeker/functional_GO_BP_ESCvsNPC_KOspecific_Lost.pdf", width=7, height=4)
pdf("output/ChIPseeker/functional_GO_BP_ESCvsNPC_WT_Gain_15.pdf", width=7, height=6)
pdf("output/ChIPseeker/functional_GO_BP_ESCvsNPC_WT_Lost_15.pdf", width=7, height=6)
dotplot(enrichGO, showCategory = 15, title = "GO_Biological Process Enrichment Analysis")
dev.off()
### NO ENRICHMENT  ESCvsNPC_WTspecific_Lost
###  ENRICHMENT  ESCvsNPC_HETspecific_Lost >50
###  ENRICHMENT  ESCvsNPC_KOspecific_Lost 9
### ENRICHMENT ESCvsNPC_WT gain 322 hits
### ENRICHMENT ESCvsNPC_WT lost 1060 hits


enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "CC") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)

pdf("output/ChIPseeker/functional_GO_CC_ESCvsNPC_HETspecific_Lost_15.pdf", width=7, height=6)
pdf("output/ChIPseeker/functional_GO_CC_ESCvsNPC_KOspecific_Lost.pdf", width=7, height=2)
pdf("output/ChIPseeker/functional_GO_CC_ESCvsNPC_WT_Gain_15.pdf", width=7, height=6)
pdf("output/ChIPseeker/functional_GO_CC_ESCvsNPC_WT_Lost_15.pdf", width=7, height=6)
dotplot(enrichGO, showCategory = 15, title = "GO_Cellular Component Enrichment Analysis")
dev.off()
### ENRICHMENT  ESCvsNPC_HETspecific_Lost
### ENRICHMENT ESCvsNPC_WT 52 hits
### ENRICHMENT ESCvsNPC_WT 63 hits


enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "MF") # "BP" (Biological Process), "8" (Cellular Component), or "MF" (Molecular Function)

pdf("output/ChIPseeker/functional_GO_MF_ESCvsNPC_HETspecific_Lost_15.pdf", width=7, height=7)
pdf("output/ChIPseeker/functional_GO_MF_ESCvsNPC_KOspecific_Lost_15.pdf", width=7, height=3)
pdf("output/ChIPseeker/functional_GO_MF_ESCvsNPC_WT_Gain_15.pdf", width=7, height=6)
pdf("output/ChIPseeker/functional_GO_MF_ESCvsNPC_WT_Lost_15.pdf", width=7, height=6)
dotplot(enrichGO, showCategory = 15, title = "GO_Molecular Function Enrichment Analysis")
dev.off()
### ENRICHMENT  ESCvsNPC_HETspecific_Lost
### ENRICHMENT ESCvsNPC_WT lost 32 hits
### ENRICHMENT ESCvsNPC_WT gain 49 hits


## Disease
enrichDO <- enrichDO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
                         
pdf("output/ChIPseeker/functional_DO_ESCvsNPC_KOspecific_Lost_15.pdf", width=7, height=7)
dotplot(enrichDO, showCategory = 15, title = "Disease Ontology Enrichment Analysis")
dev.off()
### ESCvsNPC_HETspecific_Lost withdrawal disorder
### ESCvsNPC_KOspecific_Lost many! 30


# To retrieve entrezID to gene name:
ids <- c("1812", "3778", "4986", "84152", "5733")
mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID")



# CHECK EXPRESSION
## import RNAseq TIME-EFFECT
###  RNAseq ESC to NPC __ WT
WT_expr = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_WT_ESC_vs_WT_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene,log2FoldChange,padj) %>%
    add_column(genotype_expr = "WT")  %>%
    mutate(log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05 ) #  & (abs(log2FoldChange) > 0.5)
###  RNAseq ESC to NPC __ HET
HET_expr = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_HET_ESC_vs_HET_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene,log2FoldChange,padj)%>%
    add_column(genotype_expr = "HET")  %>%
    mutate(log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05  )
###  RNAseq ESC to NPC __ KO
KO_expr = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_KO_ESC_vs_KO_NPC.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene,log2FoldChange,padj)%>%
    add_column(genotype_expr = "KO") %>%
    mutate(log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05   )
## ADD the missing genes in each comparison:
#### Necessay so that each of the RNAseq comp contains ALL the genes
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)
#### Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name
#### Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()
all_genes <- gene_id_name %>% 
  separate(gene_id, into = c("gene", "version"), sep = "\\.") %>%
  dplyr::select(gene)
WT_expr <- all_genes %>%
  left_join(WT_expr, by = "gene") %>%
  mutate(log2FoldChange = replace_na(log2FoldChange, 0),
         padj = replace_na(padj, 1),
         genotype_expr = "WT",
         significance = padj <= 0.05)
HET_expr <- all_genes %>%
  left_join(HET_expr, by = "gene") %>%
  mutate(log2FoldChange = replace_na(log2FoldChange, 0),
         padj = replace_na(padj, 1),
         genotype_expr = "HET",
         significance = padj <= 0.05)
KO_expr <- all_genes %>%
  left_join(KO_expr, by = "gene") %>%
  mutate(log2FoldChange = replace_na(log2FoldChange, 0),
         padj = replace_na(padj, 1),
         genotype_expr = "KO",
         significance = padj <= 0.05)
### Combine
expr = WT_expr %>%
    bind_rows(HET_expr) %>%
    bind_rows(KO_expr)
### Covnert gene id to gene symbol
expr$gene_symbol <- mapIds(org.Hs.eg.db,
                         keys = expr$gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
expr_geneSymbol = expr %>% 
    filter(!is.na(gene_symbol))


## import RNAseq GENOTYPE-EFFECT
###  RNAseq WT vs HET __ ESC
ESC_expr = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_ESC_HET_vs_ESC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene,log2FoldChange,padj) %>%
    add_column(time_expr = "ESC")  %>%
    mutate(log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05)
###  RNAseq WT vs HET __ NPC
NPC_expr = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_NPC_HET_vs_NPC_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene,log2FoldChange,padj) %>%
    add_column(time_expr = "NPC")  %>%
    mutate(log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05)
###  RNAseq WT vs HET __ 2dN
X2dN_expr = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_2dN_HET_vs_2dN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene,log2FoldChange,padj) %>%
    add_column(time_expr = "2dN")  %>%
    mutate(log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05)
###  RNAseq WT vs HET __ 4wN
X4wN_expr = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_4wN_HET_R3R4_vs_4wN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene,log2FoldChange,padj) %>%
    add_column(time_expr = "4wN")  %>%
    mutate(log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05)
###  RNAseq WT vs HET __ 8wN
X8wN_expr = tibble(read.csv('../001__RNAseq/output/deseq2_hg38/raw_8wN_HET_vs_8wN_WT.txt')) %>%
    separate(gene, into = c("gene", "trash"), sep ="\\.") %>%
    dplyr::select(gene,log2FoldChange,padj) %>%
    add_column(time_expr = "8wN")  %>%
    mutate(log2FoldChange = replace_na(log2FoldChange, 0),
           padj = replace_na(padj, 1),  # replace baseMean of NA with 0 and padj of NA with 1 
           significance = padj <= 0.05)
## ADD the missing genes in each comparison:
#### Necessay so that each of the RNAseq comp contains ALL the genes
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)
#### Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name
#### Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()
all_genes <- gene_id_name %>% 
  separate(gene_id, into = c("gene", "version"), sep = "\\.") %>%
  dplyr::select(gene)
ESC_expr <- all_genes %>%
  left_join(ESC_expr, by = "gene") %>%
  mutate(log2FoldChange = replace_na(log2FoldChange, 0),
         padj = replace_na(padj, 1),
         time_expr = "ESC",
         significance = padj <= 0.05)
NPC_expr <- all_genes %>%
  left_join(NPC_expr, by = "gene") %>%
  mutate(log2FoldChange = replace_na(log2FoldChange, 0),
         padj = replace_na(padj, 1),
         time_expr = "NPC",
         significance = padj <= 0.05)
X2dN_expr <- all_genes %>%
  left_join(X2dN_expr, by = "gene") %>%
  mutate(log2FoldChange = replace_na(log2FoldChange, 0),
         padj = replace_na(padj, 1),
         time_expr = "2dN",
         significance = padj <= 0.05)
X4wN_expr <- all_genes %>%
  left_join(X4wN_expr, by = "gene") %>%
  mutate(log2FoldChange = replace_na(log2FoldChange, 0),
         padj = replace_na(padj, 1),
         time_expr = "4wN",
         significance = padj <= 0.05)
X8wN_expr <- all_genes %>%
  left_join(X8wN_expr, by = "gene") %>%
  mutate(log2FoldChange = replace_na(log2FoldChange, 0),
         padj = replace_na(padj, 1),
         time_expr = "8wN",
         significance = padj <= 0.05)
### Combine
expr = ESC_expr %>%
    bind_rows(NPC_expr) %>%
    bind_rows(X2dN_expr)%>%
    bind_rows(X4wN_expr)%>%
    bind_rows(X8wN_expr)
### Covnert gene id to gene symbol
expr$gene_symbol <- mapIds(org.Hs.eg.db,
                         keys = expr$gene,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
expr_geneSymbol = expr %>% 
    filter(!is.na(gene_symbol))






# import TPM
tpm_all_sample <- read_csv("../001__RNAseq/output/tpm_hg38/tpm_all_sample.txt") %>%
  dplyr::select(-"...1")
## Join with gene names for convenience
## Read GTF file
gtf_file <- "../../Master/meta/gencode.v19.annotation.gtf"
gtf_data <- import(gtf_file)
## Extract gene_id and gene_name
gene_data <- gtf_data[elementMetadata(gtf_data)$type == "gene"]
gene_id <- elementMetadata(gene_data)$gene_id
gene_name <- elementMetadata(gene_data)$gene_name
## Combine gene_id and gene_name into a data frame
gene_id_name <- data.frame(gene_id, gene_name) %>%
  unique() %>%
  as_tibble()
## Tidy df
tpm_all_sample_tidy = as_tibble(tpm_all_sample) %>%
  gather(key = "sample", value = "tpm", -Geneid) %>%
  mutate(tpm = log2(tpm+1)) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_")
tpm_all_sample_tidy$time <-
  factor(tpm_all_sample_tidy$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
tpm_all_sample_tidy$genotype <-
  factor(tpm_all_sample_tidy$genotype,
         c("WT", "KO", "HET", "iPSCWT","iPSCpatient"))
tpm_all_sample_tidy_gene_name = tpm_all_sample_tidy %>%
  dplyr::rename(gene_id=Geneid) %>%
  left_join(gene_id_name) %>%
  unique()
## Stat
tpm_all_sample_tidy_gene_name_stat <- tpm_all_sample_tidy_gene_name %>%
  dplyr::select(-replicate) %>%
  group_by(gene_id, gene_name, time, genotype) %>%
  summarise(mean=mean(tpm), median= median(tpm), SD=sd(tpm), n=n(), SE=SD/sqrt(n)) 	




## import gene list
### WT ESC to NPC TMM-norm GAIN
gene_symbols_WT <- as_tibble(read.table("output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Gain.txt", header = FALSE, stringsAsFactors = FALSE)) %>% 
    rename(gene_symbol = V1)
### WT ESC to NPC TMM-norm LOST
gene_symbols_WT <- as_tibble(read.table("output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Lost.txt", header = FALSE, stringsAsFactors = FALSE)) %>% 
    rename(gene_symbol = V1)


### control
gene_symbols_HET <- as_tibble(read.table("output/ChIPseeker/annotation_HET_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Lost.txt", header = FALSE, stringsAsFactors = FALSE)) %>% 
    rename(gene_symbol = V1)

### genotype-specific LOST
gene_symbols_HET <- as_tibble(read.table("output/ChIPseeker/ESCvsNPC_HETspecific_Lost.txt", header = FALSE, stringsAsFactors = FALSE)) %>% 
    rename(gene_symbol = V1)
gene_symbols_KO <- as_tibble(read.table("output/ChIPseeker/ESCvsNPC_KOspecific_Lost.txt", header = FALSE, stringsAsFactors = FALSE)) %>% 
    rename(gene_name = V1)
gene_symbols_HETandKO <- as_tibble(read.table("output/ChIPseeker/ESCvsNPC_HETandKOspecific_Lost.txt", header = FALSE, stringsAsFactors = FALSE)) %>% 
    rename(gene_symbol = V1)

gene_symbols_HET_HETandKO <- gene_symbols_HET %>%
  bind_rows(gene_symbols_HETandKO) %>%
  unique()

gene_symbols_GOneuronsBrain <- as_tibble(read.table("../004__IndirectEZH1TargetId/meta/neurons_brain_related_GO_geneList_geneSymbol.txt", header = FALSE, stringsAsFactors = FALSE)) %>% 
    rename(gene_name = V1)


# plot
## with logFC time comp
expr_geneSymbol$genotype_expr <-
  factor(expr_geneSymbol$genotype_expr,
         c("WT", "HET", "KO"))
         
count_data <- expr_geneSymbol %>%
  inner_join(gene_symbols_WT) %>%                # CHANGE HERE !!!!!!!!!
  group_by(genotype_expr, significance) %>%
  summarise(up = sum(log2FoldChange > 0),
            down = sum(log2FoldChange < 0),
            total = n()) %>%
  ungroup() %>%
  group_by(genotype_expr) %>%
  mutate(total_panel = sum(total)) %>%
  ungroup()

significant_data <- expr_geneSymbol %>%
  inner_join(gene_symbols_WT) %>%                  # CHANGE HERE !!!!!!!!!
  mutate(change_group = ifelse(log2FoldChange > 0, "Greater than 0", 
                               ifelse(log2FoldChange < 0, "Less than 0", "Equal to 0"))) %>%
  filter(significance)

all_data <- expr_geneSymbol %>%
  inner_join(gene_symbols_WT) %>%                # CHANGE HERE !!!!!!!!!
  mutate(change_group = ifelse(log2FoldChange > 0, "Greater than 0", 
                               ifelse(log2FoldChange < 0, "Less than 0", "Equal to 0")))

pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression.pdf", width=4, height=5)
pdf("output/ChIPseeker/ESCvsNPC_HETandKOspecific_Lost_expression.pdf", width=4, height=5)
pdf("output/ChIPseeker/ESCvsNPC_HET_HETandKOspecific_Lost_expression.pdf", width=6, height=4)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression.pdf", width=6, height=4)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression.pdf", width=6, height=4)

ggplot() +
  geom_boxplot(data = significant_data, aes(genotype_expr, log2FoldChange, fill = change_group)) +
  geom_jitter(data = all_data, aes(genotype_expr, log2FoldChange, color = significance), alpha = 0.8, size = 0.5, width = 0.1) +
  scale_fill_manual(values = c("darkred", "darkred")) + 
  scale_color_manual(values = c("grey", "red")) +
  geom_text(data = count_data %>% filter(significance),
            aes(x = genotype_expr, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
            hjust = 0.5, vjust = 1.1, size = 3, color = "black", check_overlap = TRUE) +
  geom_text(data = count_data %>% distinct(genotype_expr, .keep_all = TRUE),
            aes(x = genotype_expr, y = -Inf, label = paste("Total:", total_panel, "genes")),
            hjust = 0.5, vjust = -0.1, size = 3, color = "black", check_overlap = TRUE) +
  theme_bw()
dev.off()



## with logFC genotype comp
expr_geneSymbol$time_expr <-
  factor(expr_geneSymbol$time_expr,
         c("ESC", "NPC", "2dN","4wN","8wN"))
         
count_data <- expr_geneSymbol %>%
  inner_join(gene_symbols_HET_HETandKO) %>%   # CHANGE HERE !!!!!!!!!
  group_by(time_expr, significance) %>%
  summarise(up = sum(log2FoldChange > 0),
            down = sum(log2FoldChange < 0),
            total = n()) %>%
  ungroup() %>%
  group_by(time_expr) %>%
  mutate(total_panel = sum(total)) %>%
  ungroup()

significant_data <- expr_geneSymbol %>%
  inner_join(gene_symbols_HET_HETandKO) %>%
  mutate(change_group = ifelse(log2FoldChange > 0, "Greater than 0", 
                               ifelse(log2FoldChange < 0, "Less than 0", "Equal to 0"))) %>%
  filter(significance)

all_data <- expr_geneSymbol %>%
  inner_join(gene_symbols_HET_HETandKO) %>%
  mutate(change_group = ifelse(log2FoldChange > 0, "Greater than 0", 
                               ifelse(log2FoldChange < 0, "Less than 0", "Equal to 0")))

pdf("output/ChIPseeker/ESCvsNPC_HET_HETandKOspecific_Lost_expression_genotypeComp.pdf", width=7, height=4)
ggplot() +
  geom_boxplot(data = significant_data, aes(time_expr, log2FoldChange, fill = change_group)) +
  geom_jitter(data = all_data, aes(time_expr, log2FoldChange, color = significance), alpha = 0.8, size = 0.5, width = 0.1) +
  scale_fill_manual(values = c("darkred", "darkred")) + 
  scale_color_manual(values = c("grey", "red")) +
  geom_text(data = count_data %>% filter(significance),
            aes(x = time_expr, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
            hjust = 0.5, vjust = 1.1, size = 3, color = "black", check_overlap = TRUE) +
  geom_text(data = count_data %>% distinct(time_expr, .keep_all = TRUE),
            aes(x = time_expr, y = -Inf, label = paste("Total:", total_panel, "genes")),
            hjust = 0.5, vjust = -0.1, size = 3, color = "black", check_overlap = TRUE) +
  theme_bw()
dev.off()


## with TPM
pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_TPM.pdf", width=5, height=6)
tpm_all_sample_tidy_gene_name_stat %>% 
    inner_join(gene_symbols_HET) %>%
    filter(genotype %in% c("WT","HET","KO")) %>%
    ggplot(., aes(genotype,log2(mean))) +
    geom_boxplot(aes(fill=time))
dev.off()

pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_GOneuronsBrain_expression_TPM.pdf", width=5, height=6)
tpm_all_sample_tidy_gene_name_stat %>% 
    inner_join(gene_symbols_HET) %>%
    inner_join(gene_symbols_GOneuronsBrain) %>%
    filter(genotype %in% c("WT","HET","KO")) %>%
    ggplot(., aes(genotype,log2(mean))) +
    geom_boxplot(aes(fill=time))
dev.off()



```
***NOTE: keep in mind that with gene name conversion, I lose genes! For example I am not able to retrieve the 904 genes Lost in HET and HETandKO; only 800 of them***


--> With log2FC we do not see that the site that specifically LOST H3K27me3 in HET are more express... Thus let's try with TPM
----> same, does not show that the genes that lose H3K27me3 from ESC to NPC in HET are more express in HET vs WT or KO !!

--> Focus only on the one related to neurons activity; when joining with genes with GO neurons/brain (broadly, from the 004 project; `004__IndirectEZH1TargetId/meta/neurons_brain_related_GO_geneList_geneSymbol.txt`), it now works: in WT expression decrease; it increase in HET from ESC to NPC! 
----> It concern only 6 genes... Let's instead collect the genes with **GO related to neurons activity (identified from the analysis) and that showed decrease of H3K27me3 in HET**


Let's try to output the **gene expression level for each of the Gene Ontology categories**!

--> It seems that only the WT ChIPseq data from ESC to NPC make sense! In agreement with gene expression changes. 
----> Let's check genes that are diff. H3K27me3 regulated in WT and with GO neurons/brain; identify the one in the mutants that are express differentially
------> We would hypothesizie these genes are diff. H3K27me3 and validate this through ChIPqPCR. 

*NOTE: 1st part below suck, then check the part `# Scatter plot to check correlation`*

```R

# library
library("rtracklayer")
library("tidyverse")
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("DOSE")
library("pathview")
library("enrichplot")
library("ggrepel")

# Import the gene list
gene_symbols <- readLines("output/ChIPseeker/ESCvsNPC_HETspecific_Lost.txt") 
gene_symbols <- readLines("output/ChIPseeker/ESCvsNPC_KOspecific_Lost.txt") 
## WT changes TMM norm
gene_symbols <- readLines("output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Gain.txt") 
gene_symbols <- readLines("output/ChIPseeker/annotation_WT_ESCvsNPC_qval25_UniqueBamTMM_geneSymbol_Lost.txt") 



genes <- mapIds(org.Hs.eg.db, 
                   keys = gene_symbols, 
                   column = "ENTREZID", 
                   keytype = "SYMBOL", 
                   multiVals = "first")


## KEGG
enrichGO <- enrichKEGG(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
## GO
enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "BP") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "CC") 
enrichGO <- enrichGO(gene   = genes,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb         = "org.Hs.eg.db",
                         ont           = "MF") # "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
## Collect the gene list for each GO
go_genes <- enrichGO@result$geneID
go_genes <- go_genes[1:15] # focus on the 15 first element


## Convert to gene symbols
### extract the first 15 GO terms and corresponding genes
go_terms_genes <- head(enrichGO@result[, c("Description", "geneID")], 15) # can change to ID top have GO:00123 format
### split the gene IDs
go_terms_genes$geneID <- strsplit(go_terms_genes$geneID, "/")
### for each GO term, convert the gene IDs to symbols
go_terms_symbols <- lapply(go_terms_genes$geneID, function(genes) {
  mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
})
### add the symbols to the dataframe
go_terms_genes$gene_symbol <- go_terms_symbols
### Tidy the data
tidy_go_terms_genes <- as_tibble(go_terms_genes) %>%
  unnest(cols = c(geneID, gene_symbol))

# plot
## with logFC all cotegories together
expr_geneSymbol$genotype_expr <-
  factor(expr_geneSymbol$genotype_expr,
         c("WT", "HET", "KO"))
         
pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_GO_BP_15.pdf", width=3, height=4)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_GO_BP_15.pdf", width=3, height=4)

tidy_go_terms_genes %>% 
    dplyr::select(gene_symbol) %>%
    unique() %>%
    left_join(expr_geneSymbol) %>%
      ggplot(., aes(x = genotype_expr,y = log2FoldChange))  +
        geom_boxplot() +
        geom_jitter(aes(color = significance), alpha = 0.8, size = 0.5, width = 0.1) +
        scale_color_manual(values = c("grey", "red")) +
        theme_bw()
dev.off()


## with logFC keeping separated each GO category
expr_geneSymbol$genotype_expr <-
  factor(expr_geneSymbol$genotype_expr,
         c("WT", "HET", "KO"))
         
count_data <- expr_geneSymbol %>%
  inner_join(tidy_go_terms_genes) %>%                # CHANGE HERE !!!!!!!!!
  group_by(genotype_expr, significance, Description) %>%
  summarise(up = sum(log2FoldChange > 0),
            down = sum(log2FoldChange < 0),
            total = n()) %>%
  ungroup() %>%
  group_by(genotype_expr) %>%
  mutate(total_panel = sum(total)) %>%
  ungroup()

significant_data <- expr_geneSymbol %>%
  inner_join(tidy_go_terms_genes) %>%                  # CHANGE HERE !!!!!!!!!
  mutate(change_group = ifelse(log2FoldChange > 0, "Greater than 0", 
                               ifelse(log2FoldChange < 0, "Less than 0", "Equal to 0"))) %>%
  filter(significance)

all_data <- expr_geneSymbol %>%
  inner_join(tidy_go_terms_genes) %>%                # CHANGE HERE !!!!!!!!!
  mutate(change_group = ifelse(log2FoldChange > 0, "Greater than 0", 
                               ifelse(log2FoldChange < 0, "Less than 0", "Equal to 0")))


pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_GO_BP_15.pdf", width=10, height=15)
pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_GO_CC_15.pdf", width=10, height=15)
pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_GO_MF_15.pdf", width=10, height=15)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_GO_BP_15.pdf", width=20, height=15)


ggplot() +
  geom_boxplot(data = significant_data, aes(genotype_expr, log2FoldChange, fill = change_group)) +
  geom_jitter(data = all_data, aes(genotype_expr, log2FoldChange, color = significance), alpha = 0.8, size = 0.5, width = 0.1) +
  scale_fill_manual(values = c("darkred", "darkred")) + 
  scale_color_manual(values = c("grey", "red")) +
  geom_text(data = count_data %>% filter(significance),
            aes(x = genotype_expr, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
            hjust = 0.5, vjust = 1.1, size = 4, color = "black", check_overlap = TRUE) +
  geom_text(data = count_data %>% distinct(genotype_expr, .keep_all = TRUE),
            aes(x = genotype_expr, y = -Inf, label = paste("Total:", total_panel, "genes")),
            hjust = 0.5, vjust = -0.1, size = 4, color = "black", check_overlap = TRUE) +
  facet_wrap(~Description, nrow = 3) +
  theme_bw() +
  theme(strip.text = element_text(size = 14))
dev.off()




### Try to add signficiant test
pvalues_df <- lapply(split(significant_data, significant_data$Description), function(data) {
  # Filter the data for 'Less than 0' change group
  data_het <- data[data$change_group == "Greater than 0" & data$genotype_expr %in% c("WT", "HET"), ] ## CHANGE IF NEEDED
  data_ko <- data[data$change_group == "Greater than 0" & data$genotype_expr %in% c("WT", "KO"), ] ## CHANGE IF NEEDED
  
  # Perform t.test if there are two levels, else set pvalue as NA
  pvalue_het <- if (length(unique(data_het$genotype_expr)) == 2) 
                    t.test(data_het$log2FoldChange ~ data_het$genotype_expr)$p.value
                else NA
  pvalue_ko <- if (length(unique(data_ko$genotype_expr)) == 2) 
                    t.test(data_ko$log2FoldChange ~ data_ko$genotype_expr)$p.value
                else NA
  
  df <- data.frame(
    comparison = c("WT_HET", "WT_KO"),
    pvalue = c(pvalue_het, pvalue_ko),
    Description = unique(data$Description)
  )
  return(df)
})

pvalues_df <- do.call(rbind, pvalues_df)  # combine the list into a single data frame



pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_GO_BP_15.pdf", width=20, height=15)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_GO_CC_15.pdf", width=20, height=15)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_GO_MF_15.pdf", width=20, height=15)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_KEGG_15.pdf", width=20, height=15)

pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_GO_BP_15.pdf", width=20, height=15)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_GO_CC_15.pdf", width=20, height=15)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_GO_MF_15.pdf", width=20, height=15)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_KEGG_15.pdf", width=20, height=15)



ggplot() +
  geom_boxplot(data = significant_data, aes(genotype_expr, log2FoldChange, fill = change_group)) +
  geom_jitter(data = all_data, aes(genotype_expr, log2FoldChange, color = significance), alpha = 0.8, size = 0.5, width = 0.1) +
  scale_fill_manual(values = c("darkred", "darkred")) + 
  scale_color_manual(values = c("grey", "red")) +
  geom_text(data = count_data %>% filter(significance),
            aes(x = genotype_expr, y = Inf, label = paste(up, "genes up\n", down, "genes down")),
            hjust = 0.5, vjust = 1.1, size = 4, color = "black", check_overlap = TRUE) +
  geom_text(data = count_data %>% distinct(genotype_expr, .keep_all = TRUE),
            aes(x = genotype_expr, y = -Inf, label = paste("Total:", total_panel, "genes")),
            hjust = 0.5, vjust = -0.1, size = 4, color = "black", check_overlap = TRUE) +
  geom_text(data = pvalues_df[pvalues_df$comparison == "WT_HET",], aes(x = 1, y = -9, 
            label = paste("P value (", comparison, "): ", round(pvalue, 3), sep = "")), 
            size = 3.5, hjust = 0, color = "black") +
  geom_text(data = pvalues_df[pvalues_df$comparison == "WT_KO",], aes(x = 1, y = -10, 
            label = paste("P value (", comparison, "): ", round(pvalue, 3), sep = "")), 
            size = 3.5, hjust = 0, color = "black") +
  facet_wrap(~Description, nrow = 3) +
  theme_bw() +
  theme(strip.text = element_text(size = 14))
dev.off()




### Scatter plot to check correlation
###  check out which genes and try isolate the outlier; different between genotypes
all_data %>% filter(significance == TRUE, Description == "forebrain development") %>% dplyr::select(gene_symbol) %>% unique()

all_data_corr = all_data %>%
    filter(significance == TRUE, Description == "forebrain development") %>%  # !!! CHANGE HERE
    dplyr::select(gene_symbol, genotype_expr, log2FoldChange) %>% 
    spread(key = genotype_expr, value = log2FoldChange) %>%
    replace_na(list(WT = 0, HET = 0, KO = 0)) %>%  # replace NA in columns 'WT', 'HET', 'KO' with 0
    mutate(diff = abs(WT - KO)) # !!! CHANGE GNOTPYE


pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_forebrainDev_corr_WT_KO.pdf", width=20, height=15) # !!! CHANGE GNOTPYE
ggplot(all_data_corr, aes(x = WT, y = KO)) +    # !!! CHANGE GNOTPYE
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +  # horizontal line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +  # vertical line at 0
  geom_point(aes(shape = diff > 1, size = diff > 1, color = diff > 1)) + # change shape and size based on diff > 1
  scale_shape_manual(values = c(16, 17)) + # specify shapes
  scale_size_manual(values = c(3, 5)) +  # specify sizes
  geom_text(data = subset(all_data_corr, diff > 1), aes(label = gene_symbol), hjust = -0.3, vjust = -0.3, check_overlap = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") + # add linear regression line
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2 Fold Change (WT)", 
       y = "Log2 Fold Change (KO)",            # !!! CHANGE GNOTPYE
       title = "Correlation between Log2 Fold Changes for WT and HET Genotypes",
       shape = "Difference > 1",
       size = "Difference > 1",
       color = "Difference > 1") +
  theme_bw() +
  theme(legend.position = "bottom")  # place legend at bottom
dev.off()

pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_forebrainDev_corr_WT_HET_2.pdf", width=20, height=15)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_forebrainDev_corr_WT_KO_2.pdf", width=20, height=15) # !!! CHANGE GNOTPYE
ggplot(all_data_corr, aes(x = WT, y = KO)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +  # horizontal line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +  # vertical line at 0
  geom_point(aes(shape = diff > 2, size = diff > 2, color = diff > 2)) + # change shape and size based on diff > 1
  scale_shape_manual(values = c(16, 17)) + # specify shapes
  scale_size_manual(values = c(3, 5)) +  # specify sizes
  geom_text(data = subset(all_data_corr, diff > 2), aes(label = gene_symbol), hjust = -0.3, vjust = -0.3, check_overlap = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") + # add linear regression line
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2 Fold Change (WT)", 
       y = "Log2 Fold Change (KO)", 
       title = "Correlation between Log2 Fold Changes for WT and KO Genotypes",
       shape = "Difference > 2",
       size = "Difference > 2",
       color = "Difference > 2") +
  theme_bw() +
  theme(legend.position = "bottom")  # place legend at bottom
dev.off()
####


# Both genotypes scatter plot
## FOR LOST
all_data_corr <- all_data %>%
    filter(significance == TRUE, Description %in% c( "axon development", "axonogenesis", "regulation of nervous system development", "regulation of neurogenesis", "neuron projection guidance", "axon guidance", "central nervous system neuron differentiation") ) %>%   # !!! CHANGE HERE
    dplyr::select(gene_symbol, genotype_expr, log2FoldChange) %>%
    unique() %>%
    spread(key = genotype_expr, value = log2FoldChange) %>%
    replace_na(list(WT = 0, HET = 0, KO = 0)) %>%  
    mutate(diff_HET = abs(WT - HET),
           diff_KO = abs(WT - KO))
## FOR GAIN
all_data_corr <- all_data %>%
    filter(significance == TRUE, Description %in% c( "regulation of trans-synaptic signaling", "amine transport", "monoamine transport", "regulation of amine transport", "catecholamine transport", "dopamine transport") ) %>%   # !!! CHANGE HERE
    dplyr::select(gene_symbol, genotype_expr, log2FoldChange) %>%
    unique() %>%
    spread(key = genotype_expr, value = log2FoldChange) %>%
    replace_na(list(WT = 0, HET = 0, KO = 0)) %>%  
    mutate(diff_HET = abs(WT - HET),
           diff_KO = abs(WT - KO))
## paste list of GO: forebrain development, modulation of chemical synaptic transmission, regulation of trans-synaptic signaling, amine transport, monoamine transport, regulation of amine transport, catecholamine transport, dopamine transport, axon development, axonogenesis, regulation of nervous system development, regulation of neurogenesis, neuron projection guidance, axon guidance, central nervous system neuron differentiation


pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_forebrainDev_corr_WT_HET_KO.pdf", width=18, height=10)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_modulationOfChemicalSynapticTransmission_corr_WT_HET_KO.pdf", width=18, height=10)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_regulationOfTranssynapticSignaling_corr_WT_HET_KO.pdf", width=18, height=10)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_amine_corr_WT_HET_KO.pdf", width=18, height=10)
all_data_corr_long <- all_data_corr %>%
    gather(key = "Comparison", value = "log2FC", HET, KO)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_catecholamineDopamineTransport_corr_WT_HET_KO.pdf", width=18, height=10)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_axon_corr_WT_HET_KO.pdf", width=18, height=10)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_regulationNeuronDev_corr_WT_HET_KO.pdf", width=18, height=10)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_neuronGuidance_corr_WT_HET_KO.pdf", width=18, height=10)
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_CNSneuronDifferentiation_corr_WT_HET_KO.pdf", width=18, height=10)

all_data_corr_long <- all_data_corr %>%
    gather(key = "Comparison", value = "log2FC", HET, KO)
ggplot(all_data_corr_long, aes(x = WT, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(shape = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2, 
                   size = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2, 
                   color = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2)) +
    geom_text(data = all_data_corr_long %>% filter(ifelse(Comparison == "HET", diff_HET, diff_KO) > 2), 
              aes(label = gene_symbol), hjust = -0.3, vjust = -0.3, check_overlap = TRUE) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
    scale_shape_manual(values = c(16, 17)) +
    scale_size_manual(values = c(3, 5)) +
    scale_color_manual(values = c("black", "red")) +
    facet_wrap(~Comparison, scales = "free") +
    labs(x = "Log2 Fold Change (WT)", 
         y = "Log2 Fold Change (HET/KO)", 
         title = "Correlation between Log2 Fold Changes for WT and HET/KO Genotypes",
         shape = "Difference > 2",
         size = "Difference > 2",
         color = "Difference > 2") +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()


## All GO together with colored dot. --> each gene assigned only to one GO category (preferring the category that has more genes)
go_counts <- all_data %>%
    filter(Description %in% c("axon development", "axonogenesis", "regulation of nervous system development", "regulation of neurogenesis", "neuron projection guidance", "axon guidance", "central nervous system neuron differentiation")) %>%
    group_by(Description) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
go_counts <- all_data %>%
    filter(Description %in% c("regulation of trans-synaptic signaling", "amine transport", "monoamine transport", "regulation of amine transport", "catecholamine transport", "dopamine transport")) %>%
    group_by(Description) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
# Join back to original data and filter to keep only the most common GO category for each gene
all_data_single_go <- all_data %>%
    inner_join(go_counts, by = "Description") %>%
    group_by(gene_symbol) %>%
    filter(count == max(count)) %>%
    ungroup() %>%
    dplyr::select(-count)

all_data_corr <- all_data_single_go %>%
    filter(significance == TRUE) %>%
    dplyr::select(gene_symbol, genotype_expr, log2FoldChange, Description) %>%
    unique() %>%
    spread(key = genotype_expr, value = log2FoldChange) %>%
    replace_na(list(WT = 0, HET = 0, KO = 0)) %>%
    mutate(diff_HET = abs(WT - HET),
           diff_KO = abs(WT - KO))

all_data_corr_long <- all_data_corr %>%
    gather(key = "Comparison", value = "log2FC", HET, KO)

pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_all_corr_WT_HET_KO.pdf", width=18, height=10)
ggplot(all_data_corr_long, aes(x = WT, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(shape = ifelse(Comparison == "HET", diff_HET, diff_KO) > 1, 
                   size = ifelse(Comparison == "HET", diff_HET, diff_KO) > 1, 
                   color = Description)) + # Using Description for color
    geom_text(data = all_data_corr_long %>% filter(ifelse(Comparison == "HET", diff_HET, diff_KO) > 1), 
              aes(label = gene_symbol), hjust = -0.3, vjust = -0.3, check_overlap = FALSE) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
    scale_shape_manual(values = c(16, 17)) +
    scale_size_manual(values = c(3, 5)) +
    # You may need to define specific colors for your GO categories here
    facet_wrap(~Comparison, scales = "free") +
    labs(x = "Log2 Fold Change (WT)", 
         y = "Log2 Fold Change (HET/KO)", 
         title = "Correlation between Log2 Fold Changes for WT and HET/KO Genotypes",
         shape = "Difference > 1",
         size = "Difference > 1",
         color = "GO Category") +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()



pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_all_corr_WT_HET_KO.pdf", width=18, height=10)
ggplot(all_data_corr_long, aes(x = WT, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(shape = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2, 
                   size = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2, 
                   color = Description)) + # Using Description for color
    geom_text_repel(data = all_data_corr_long %>% filter(ifelse(Comparison == "HET", diff_HET, diff_KO) > 2), 
                    aes(label = gene_symbol), nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 10) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
    scale_shape_manual(values = c(16, 17)) +
    scale_size_manual(values = c(3, 5)) +
    # You may need to define specific colors for your GO categories here
    facet_wrap(~Comparison, scales = "free") +
    labs(x = "Log2 Fold Change (WT)", 
         y = "Log2 Fold Change (HET/KO)", 
         title = "Correlation between Log2 Fold Changes for WT and HET/KO Genotypes",
         shape = "Difference > 2",
         size = "Difference > 2",
         color = "GO Category") +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()



## All GO together with colored dot. --> each gene assigned only to one GO category (preferring the category that has more genes) + ONLY display genes that have a differential behavior exclusively in one of the genotypes (HET or KO), but not both


genotype_specific_genes <- all_data_corr %>%
  filter((diff_HET > 2 & diff_KO <= 2) | 
         (diff_KO > 2 & diff_HET <= 2) | 
         (diff_HET > 2 & diff_KO > 2 & sign(WT - HET) != sign(WT - KO))) %>%
  dplyr::select(gene_symbol)

# Filter data to keep only genotype-specific genes
all_data_corr_long_specific <- all_data_corr_long %>%
  filter(gene_symbol %in% genotype_specific_genes$gene_symbol)

# Create the plot
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_all_corr_WT_HET_KO_Specific.pdf", width=18, height=10)
ggplot(all_data_corr_long_specific, aes(x = WT, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(shape = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2, 
                   size = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2, 
                   color = Description)) +
    geom_text_repel(data = all_data_corr_long_specific %>% filter(ifelse(Comparison == "HET", diff_HET, diff_KO) > 2), 
                    aes(label = gene_symbol), nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 10) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
    scale_shape_manual(values = c(16, 17)) +
    scale_size_manual(values = c(3, 5)) +
    facet_wrap(~Comparison, scales = "free") +
    labs(x = "Log2 Fold Change (WT)", 
         y = "Log2 Fold Change (HET/KO)", 
         title = "Genotype-Specific Differential Genes: Correlation between Log2 Fold Changes for WT and HET/KO Genotypes",
         shape = "Difference > 2",
         size = "Difference > 2",
         color = "GO Category") +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()

## All GO together with colored dot. --> each gene assigned only to one GO category (preferring the category that has more genes) + ONLY display genes that have a differential behavior exclusively in one of the genotypes (HET or KO), but not both Colored in red
### Identify genotype-specific genes
genotype_specific_genes <- all_data_corr %>%
  filter((diff_HET > 2 & diff_KO <= 2) | 
         (diff_KO > 2 & diff_HET <= 2) | 
         (diff_HET > 2 & diff_KO > 2 & sign(WT - HET) != sign(WT - KO))) %>%
  dplyr::select(gene_symbol)

# Filter the genes to label
label_data <- all_data_corr_long %>%
    filter(ifelse(Comparison == "HET", diff_HET, diff_KO) > 2)


# Assign the colors
label_data$text_color <- ifelse(label_data$gene_symbol %in% genotype_specific_genes$gene_symbol, "red", "black")



### Add a column to specify text color for specific genes
all_data_corr_long <- all_data_corr_long %>%
  mutate(text_color = if_else(gene_symbol %in% genotype_specific_genes$gene_symbol, "red", "black"))




### Create the plot
pdf("output/ChIPseeker/ESCvsNPC_WT_Lost_expression_all_corr_WT_HET_KO_SpecificInRed.pdf", width=18, height=10)

ggplot(all_data_corr_long, aes(x = WT, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(shape = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2, 
                   size = ifelse(Comparison == "HET", diff_HET, diff_KO) > 2, 
                   color = Description)) +
    # Use label_data with color mapping
    geom_text_repel(data = label_data, 
                    aes(label = gene_symbol), 
                    color = label_data$text_color,
                    nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 10) +
    # Define specific colors for your GO categories if needed
    scale_color_manual(values = c("axon development" = "orange", "central nervous system neuron differentiation" = "green", "neuron projection guidance" = "blue", "regulation of nervous system development" = "purple")) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
    scale_shape_manual(values = c(16, 17)) +
    scale_size_manual(values = c(3, 5)) +
    facet_wrap(~Comparison, scales = "free") +
    labs(x = "Log2 Fold Change (WT)", 
         y = "Log2 Fold Change (HET/KO)", 
         title = "Correlation between Log2 Fold Changes for WT and HET/KO Genotypes",
         shape = "Difference > 2",
         size = "Difference > 2",
         color = "GO Category") +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()


pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_all_corr_WT_HET_KO_SpecificInRed.pdf", width=18, height=10)
pdf("output/ChIPseeker/ESCvsNPC_WT_Gain_expression_all_corr_WT_HET_KO_SpecificInRed_1.pdf", width=18, height=10)

ggplot(all_data_corr_long, aes(x = WT, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(shape = ifelse(Comparison == "HET", diff_HET, diff_KO) > 1, 
                   size = ifelse(Comparison == "HET", diff_HET, diff_KO) > 1, 
                   color = Description)) +
    # Use label_data with color mapping
    geom_text_repel(data = label_data, 
                    aes(label = gene_symbol), 
                    color = label_data$text_color,
                    nudge_x = 0.5, nudge_y = 0.5, max.overlaps = 10) +
    # Define specific colors for your GO categories if needed
    scale_color_manual(values = c("regulation of trans-synaptic signaling" = "yellow", "amine transport" = "brown", "monoamine transport" = "red")) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
    scale_shape_manual(values = c(16, 17)) +
    scale_size_manual(values = c(3, 5)) +
    facet_wrap(~Comparison, scales = "free") +
    labs(x = "Log2 Fold Change (WT)", 
         y = "Log2 Fold Change (HET/KO)", 
         title = "Correlation between Log2 Fold Changes for WT and HET/KO Genotypes",
         shape = "Difference > 1",
         size = "Difference > 1",
         color = "GO Category") +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()






# ChatGPT heatmap method with library("ComplexHeatmap") library("circlize")
## --> See Python Code
write.table(all_data_corr, file = "output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Gain_expression_modulationOfChemicalSynapticTransmission_corr_WT_HET_KO.txt", sep = "\t", row.names = FALSE)
write.table(all_data_corr, file = "output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Gain_expression_amine_corr_WT_HET_KO.txt", sep = "\t", row.names = FALSE)
write.table(all_data_corr, file = "output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Gain_expression_catecholamineDopamineTransport_corr_WT_HET_KO.txt", sep = "\t", row.names = FALSE)
write.table(all_data_corr, file = "output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_CNSneuronDifferentiation_corr_WT_HET_KO.txt", sep = "\t", row.names = FALSE)
write.table(all_data_corr, file = "output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_neuronGuidance_corr_WT_HET_KO.txt", sep = "\t", row.names = FALSE)
write.table(all_data_corr, file = "output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_regulationNeuronDev_corr_WT_HET_KO.txt", sep = "\t", row.names = FALSE)
write.table(all_data_corr, file = "output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_regulationNeuronDev_corr_WT_HET_KO.txt", sep = "\t", row.names = FALSE)
write.table(all_data_corr, file = "output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_axon_corr_WT_HET_KO.txt", sep = "\t", row.names = FALSE)


### Only keep the genes that are induced in WT 

pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_WTinduced_GO_BP_15.pdf", width=10, height=15)

expr_geneSymbol %>% 
    inner_join(tidy_go_terms_genes) %>% 
    group_by(gene_symbol) %>% 
    filter(any(genotype_expr == "WT" & log2FoldChange > 0)) %>%
    dplyr::select(gene_symbol,log2FoldChange,padj,genotype_expr,Description,significance) %>%
    unique() %>%
      ggplot(., aes(x = genotype_expr,y = log2FoldChange))  +
        geom_boxplot() +
        geom_jitter(aes(color = significance), alpha = 0.8, size = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        facet_wrap(~Description)
dev.off()

# plot with TPM


pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_TPM_GO_BP_15.pdf", width=25, height=26)
tpm_all_sample_tidy_gene_name_stat %>% 
    inner_join(tidy_go_terms_genes %>% dplyr::rename("gene_name"="gene_symbol")) %>%
    filter(genotype %in% c("WT","HET"), 
           time %in% c("ESC","NPC")) %>%
    ungroup() %>% unique() %>%
      ggplot(., aes(time,mean)) +
        geom_boxplot(aes(fill=genotype))+
        facet_wrap(~Description)
dev.off()


# over the entire time-course log_tpm

pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_TPM_timecourse_GO_BP_15.pdf", width=25, height=26)
tpm_all_sample_tidy_gene_name_stat %>% 
    inner_join(tidy_go_terms_genes %>% dplyr::rename("gene_name"="gene_symbol")) %>%
    filter(genotype %in% c("WT","HET","KO"), 
           time %in% c("ESC","NPC","2dN", "8wN")) %>%
ggplot(., aes(x = time, y = mean, color = genotype, group = genotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
  geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
  facet_wrap(~Description, scale = "free", nrow = 3) +
  scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()

# over the entire time-course with rlog counts
## log rlog counts
load("../001__RNAseq/output/deseq2_hg38/ddsTC_rld_filter.RData")
rlog_counts
rlog_counts_matrix <- assay(rlog_counts) 
## Make a clean table
rlog_counts_tidy <- as_tibble(rlog_counts_matrix, rownames = "gene_id") %>%
  gather(key = "sample", value = "rlog_counts", -gene_id) %>%
  separate(sample, into = c("time", "genotype", "replicate"), sep = "_") %>%
  left_join(gene_id_name)

## Compil with GO list
rlog_counts_tidy_GO <- rlog_counts_tidy %>%
  inner_join(tidy_go_terms_genes %>% dplyr::rename("gene_name"="gene_symbol") , by = "gene_name")


rlog_counts_tidy_GO$time <-
  factor(rlog_counts_tidy_GO$time,
         c("ESC", "NPC", "2dN", "4wN", "8wN"))
rlog_counts_tidy_GO$genotype <-
  factor(rlog_counts_tidy_GO$genotype,
         c("WT", "KO", "HET"))



pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_rlog_timecourse_GO_BP_15.pdf", width=25, height=26)
pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_rlog_timecourse_GO_CC_15.pdf", width=25, height=26)
pdf("output/ChIPseeker/ESCvsNPC_HETspecific_Lost_expression_rlog_timecourse_GO_MF_15.pdf", width=25, height=26)
ggplot(rlog_counts_tidy_GO, aes(x = time, y = rlog_counts, color = genotype, group = genotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, size = 1) +
  geom_point(stat = "summary", fun = mean, shape = 18, size = 3, stroke = 1.5) +
  facet_wrap(~Description, scale = "free", nrow = 3) +
  scale_color_manual(values=c("WT" = "black", "HET" = "blue", "KO" = "red")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),        # Increase facet_wrap title panel text size
    axis.title.x = element_text(size = 16)       # Increase x-axis legend text size
  )
dev.off()


```

--> Code is great to explore expression of the genes that are out from GO hits!

**Conclusion for ChIPseq with HET and KO:**

--> With log2FC we do NOT see higher log2FC in HET versus WT... Same using TPM...

--> Filtering for only the genes induced in the WT, also showed HET/KO less induced / WT!

--> Checking instead DEGs between genotypes at each time-point is even weirder; genes are mostly downregulated in HET; at each time-point...

--> Let's only focus on the WT ChIPseq ESC and NPC that make sens

The ChIPseq normalization using TMM with uniqueBAM does not work great (poorly in agreement with gene expression changes). Let's try using housekeeping genes for normalization with and without uniqueBAM. **Well, in the end; it worked great for WT only. And using housekeeping lead to similar results**



Dig in the code previously and output gene list and find these that may be worst investigating:
- Need to be in agreement with our HET/KO phenotype (could be more activtiy for HET; or more H3K27me3; or maturation faster. And the opposite for KO; more developmental related and slow/delay in maturation).
- for that purpose; let's create an excell table with all the genes with GO neuron/brain and check on IGV if changes of expression and H3K27me3 looks real; then dig in that short list. File is `output/ChIPseeker/candidate_ChIPqPCR.xlsx` in Google Drive. Let's use files where I generated Python_heatmaps (output of console) to investigate the gene candidates (`output/ChIPseeker/all_data_corr_ESCvsNPC_WT*`)

--> The IGV is pretty clean and expected when we check diff in H3K27me3 or expression; not sure that is the good strategy.. Let's instead directly find interesting gene!

--> ChatGPT has been used to dig into the gene list; to find HET candidates (61 unique genes in total `output/ChIPseeker/het_genotype_genes_with_categories.xlsx`):

- **Potential Gain of H3K27me3 in HET (H3K27me3 Lost)**:
WT log2FC is 0, and HET log2FC is negative (indicating that despite being lost in WT, the mark may be gained in HET as expression decreased). Genes that may have gained the H3K27me3 mark in HET despite being lost in WT.
- **Potential Less Loss or Maintenance of H3K27me3 in HET (H3K27me3 Lost)**:
HET log2FC is less than WT log2FC (indicating the reduction in expression is less important in HET, suggesting potential maintenance of the H3K27me3 mark). Genes that may have maintained or lost less of the H3K27me3 mark in HET.
- **More Gain of H3K27me3 in HET (H3K27me3 Gain)**:
HET log2FC is more strongly downregulated than WT log2FC (indicating more gain of the H3K27me3 mark in HET). Genes that have more gain of the H3K27me3 mark in HET.
- **Unique Downregulation in HET (H3K27me3 Gain)**:
WT log2FC is 0, and HET log2FC is negative (indicating genes that are uniquely more downregulated in HET, possibly due to gain of the H3K27me3 mark). Genes that are uniquely more downregulated in HET, possibly due to gain of the H3K27me3 mark.

--> Dig in the gene list and write description for each genes.

--> I did not found many HET genes related to activity. Or if so; they are LESS expressed in HET; so hard to conlude that HET neurons have higher neuronal activity...




### Heatmaps in Python
Let's not use the weird staticical comparison; instead; re-do GO analysis for each of Gain and Lost in WT. Then check directly what are the interesed GO to look at. And do the scatter plot representation and associated heatmaps of only the differential.

ChatGPT Python code for **heatmaps representation**, this is CLEAR and CONCISE. Let's isolate the genes and make a clean table as the one I saved for the other GO comparison:

```bash
conda activate deseq2
# conda install seaborn
python3
```

Enter python interpreter:

**Version WT-ordered only:**

```python
import pandas as pd
from matplotlib.patches import Rectangle
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv("output/ChIPseeker/all_data_corr_output.console", sep="\t", comment="#")

# Define a threshold
threshold = 2
filtered_data = data[(data['diff_HET'].abs() > threshold) | (data['diff_KO'].abs() > threshold)]

heatmap_data = filtered_data[['gene_symbol', 'WT', 'HET', 'KO']].set_index('gene_symbol')

heatmap_data.sort_values(by='WT', ascending=False, inplace=True)

opposite_genes_HET = heatmap_data[(heatmap_data['KO'] < heatmap_data['WT'])].index
opposite_genes_KO = heatmap_data[(heatmap_data['HET'] < heatmap_data['WT'])].index

fig, axs = plt.subplots(figsize=(10, 10))

sns.heatmap(heatmap_data[['WT', 'HET', 'KO']], cmap='coolwarm', center=0, annot=False, ax=axs)

for gene in heatmap_data.index:
    if gene in opposite_genes_HET:
        axs.add_patch(Rectangle((1, heatmap_data.index.get_loc(gene)), 1, 1, fill=False, edgecolor='green', lw=3))

    if gene in opposite_genes_KO:
        axs.add_patch(Rectangle((2, heatmap_data.index.get_loc(gene)), 1, 1, fill=False, edgecolor='green', lw=3))

axs.set_title(f'Heatmap')
axs.set_xlabel('Genotype')
axs.set_ylabel('Gene')

plt.tight_layout()

# Save the figure as PDF
plt.savefig("output/ChIPseeker/heatmap_Lost_ForebrainDev.pdf")
```


**Version splitting per categories:**
```python
import pandas as pd
from matplotlib.patches import Rectangle
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv("output/ChIPseeker/all_data_corr_output.console", sep="\t", comment="#")
data.head() 

# Define a threshold
threshold = 2
filtered_data = data[(data['diff_HET'].abs() > threshold) | (data['diff_KO'].abs() > threshold)]

heatmap_data = filtered_data[['gene_symbol', 'WT', 'HET', 'KO']].set_index('gene_symbol')

heatmap_data['category'] = 'More in HET'
heatmap_data.loc[heatmap_data['KO'] > heatmap_data['HET'], 'category'] = 'More in KO'
heatmap_data.loc[(heatmap_data['WT'] > heatmap_data['HET']) & (heatmap_data['WT'] > heatmap_data['KO']), 'category'] = 'Less in HET'
heatmap_data.loc[(heatmap_data['WT'] > heatmap_data['KO']) & (heatmap_data['WT'] > heatmap_data['HET']), 'category'] = 'Less in KO'

heatmap_data.sort_values(by=['category', 'WT'], ascending=[True, False], inplace=True)

opposite_genes = {
    'More in HET': heatmap_data[(heatmap_data['category'] == 'More in HET') & (heatmap_data['KO'] < heatmap_data['WT'])].index,
    'More in KO': heatmap_data[(heatmap_data['category'] == 'More in KO') & (heatmap_data['HET'] < heatmap_data['WT'])].index,
    'Less in KO': heatmap_data[(heatmap_data['category'] == 'Less in KO') & (heatmap_data['HET'] > heatmap_data['WT'])].index
}

categories = ['More in HET', 'More in KO', 'Less in HET', 'Less in KO']
fig, axs = plt.subplots(len(categories), figsize=(10, len(categories)*5))


for i, category in enumerate(categories):
    category_data = heatmap_data[heatmap_data['category'] == category]
    
    if not category_data.empty:  # Check if the DataFrame is not empty
        sns.heatmap(category_data[['WT', 'HET', 'KO']], cmap='coolwarm', center=0, annot=False, ax=axs[i])
        
        for gene in category_data.index:
            if gene in opposite_genes[category]:
                axs[i].add_patch(Rectangle((0, category_data.index.get_loc(gene)), 3, 1, fill=False, edgecolor='green', lw=3))
        
        axs[i].set_title(f'{category} Genes (n = {len(category_data)})')
        axs[i].set_xlabel('Genotype')
        axs[i].set_ylabel('Gene')

plt.tight_layout()

# Save the figure as PDF
plt.savefig("output/ChIPseeker/heatmap_Lost_ForebrainDev.pdf")
```



**Version splitting per categories; clean/beautifull:**

Data; Gain/Lost GO_related significant genes:
- output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Gain_expression_modulationOfChemicalSynapticTransmission_corr_WT_HET_KO.txt
- output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Gain_expression_amine_corr_WT_HET_KO.txt
- output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_axon_corr_WT_HET_KO.txt
- output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_regulationNeuronDev_corr_WT_HET_KO.txt
- output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_neuronGuidance_corr_WT_HET_KO.txt
- output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_CNSneuronDifferentiation_corr_WT_HET_KO.txt
```python
import pandas as pd
from matplotlib.patches import Rectangle
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv("output/ChIPseeker/all_data_corr_ESCvsNPC_WT_Lost_expression_CNSneuronDifferentiation_corr_WT_HET_KO.txt", sep="\t", comment="#")   # CHANGE HERE !!!!!!!!!!


# Define a threshold
threshold = 2
filtered_data = data[(data['diff_HET'].abs() > threshold) | (data['diff_KO'].abs() > threshold)]

heatmap_data = filtered_data[['gene_symbol', 'WT', 'HET', 'KO']].set_index('gene_symbol')

# Define genes with opposite behavior
opposite_genes = heatmap_data[
    ((heatmap_data['WT'] > 0) & (heatmap_data['KO'] < 0)) |
    ((heatmap_data['WT'] > 0) & (heatmap_data['HET'] < 0)) |
    ((heatmap_data['WT'] < 0) & (heatmap_data['KO'] > 0)) |
    ((heatmap_data['WT'] < 0) & (heatmap_data['HET'] > 0))
].index.tolist()

# Sort the data based on 'WT' column
heatmap_data.sort_values(by='WT', ascending=False, inplace=True)

# Plot the heatmap
fig, ax = plt.subplots(figsize=(4, 12))  # Adjust the width of the figure here

sns.heatmap(heatmap_data[['WT', 'HET', 'KO']], cmap='coolwarm', center=0, annot=False, ax=ax, cbar_kws={'orientation': 'horizontal', 'pad': 0.03})

# Highlight genes with opposite behavior
for gene in opposite_genes:
    ax.add_patch(Rectangle((0, heatmap_data.index.get_loc(gene)), 3, 1, fill=False, edgecolor='green', lw=3))

ax.set_title('Gene Expression Across Genotypes', loc='center', pad=20, fontsize=16)
ax.set_xlabel('Genotype', fontsize=14)
ax.set_ylabel('Gene', fontsize=14)
ax.yaxis.set_tick_params(rotation=0, labelsize=10)  # Rotate the gene names for better visibility

plt.tight_layout()

# Save the figure as PDF
plt.savefig("output/ChIPseeker/heatmap_Lost_CNSneuronDifferentiation.pdf")  # CHANGE HERE !!!!!!!!!!

```





### Clean code to do scatter plots and heatmaps from genes H3K27me3-dynamics ESC to NPC in WT with GO related to neurons/brain


- Perform GO analysis in WT gain/lost H3K27me3 from ESC to NPC
- Classify GO categories into higher order GO family (put together closely related categories) 
- Generate scatterplot graphs to show genes with a log2FC diff. 2 between genotypes (genes that behave differentially)
- Generate heatmaps for these genes and highlight the one with opposite behavior between genotypes
- Try to explain phenotype with key genes
- Propose a list of genes to validate by ChIPqPCR in mutants


XXX COPY PASTE FROM THE TOP XXX















#### THOR with housekeeping genes normalization

**Housekeeping gene list (bed) download** from [here](https://reg-gen.readthedocs.io/en/latest/thor/tool_usage.html): "We use regions 500 bps upstream of all housekeeping genes described by Eisenberg and Levanon (43) (C1orf43, CHMP2A, EMC7, GPI, PSMB2, PSMB4, RAB7A, REEP5, SNRPD3, VCP, VPS29) as control regions for the human genome" from the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5175345/).

--> In the file I download that is hg19 bed version; so I need re-generate bed file for hg38; moreover the genes are: EMC7, GPI, PSMB2, PSMB4, RAB7A, REEP5, SNRPD3, VCP, VPS29. 

- *NOTE: Check the 500bp upstream from TSS is taken (check on IGV with hg19): OK; not exactly, but I converted the bed from hg19 to hg38 with this [webtool](https://genome.ucsc.edu/cgi-bin/hgLiftOver).*
- *NOTE: Check on IGV on the CutRun data that these genes are similarly enriched (or not) in H3K27me3 in WT/HET/KO: OK*


--> This has been run within the previous part `## THOR with different SF (DiffBindTMM from uniquely/NON-uniquely aligned bam AND ChIPseqSpikeInFree ony)`; as `*housekeep*`

In the end, **only the time-effect WT ESC to NPC looks good**, replicates are clean and changes is in agreement with gene expression. Same as when using the TMM-method...






### Peak strenght and length (peak characteristic comparison) in WT ESC vs NPC

Let's compare the peak from ESC to NPC in WT (uniqueBAM_TMM normalzied with THOR works great):

1. Check the raw MACS2 peak:
- Use MACS2 pool peak to obtain bed file coordinate for either ESC or NPC peak
- Then use the THOR normalized bigwig to generate the deepTools heatmap


2. Check list of differential peak from THOR (UniqueBamTMM). For the heatmap plot cluster per gain/lost (just k = 2)
- Assign peak to genes and check
- And simply check peak coordinates




#### raw MACS2 peak for peak strenght/length comparison in WT ESC vs NPC

Need to decipher whether we used macs2 or macs2_unique (from uniquely aligned read BAM file); to answer this question; let's load our macs2 tracks to IGV and pick the one that fit well with the THOR bigwig files (`THOR_WT_ESCvsNPC_UniqueBamTMM`).

--> Let's pick unique qval 2.3 for MACS2; more in agreement with the bigwig (the other, non unique MACS2, seems to contain more false-positive peaks)

```bash
# 1. Raw MACS2 method
## Files location
002__ChIPseq/output/macs2_unique/broad_blacklist_qval2.30103/*pool*.broadPeak ## MACS2 unique
002__ChIPseq/output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/*bw ## THOR uniqueBAM
### concatenate ESC and NPC raw peaks
cat output/macs2_unique/broad_blacklist_qval2.30103/ESC_WT_H3K27me3_pool_peaks.broadPeak output/macs2_unique/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_pool_peaks.broadPeak > output/macs2_unique/broad_blacklist_qval2.30103/ESC_NPC_WT_H3K27me3_pool_peaks.broadPeak
bedtools sort -i output/macs2_unique/broad_blacklist_qval2.30103/ESC_NPC_WT_H3K27me3_pool_peaks.broadPeak > output/macs2_unique/broad_blacklist_qval2.30103/ESC_NPC_WT_H3K27me3_pool_peaks_sort.broadPeak

## Generate median bigwig files for THOR
conda activate BedToBigwig
sbatch scripts/bigwigmerge_THOR_WT_ESCvsNPC_UniqueBamTMM.sh # 3479177 ok

## deepTools plot
conda activate deeptools
sbatch scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_peak.sh # 3486731 ok
sbatch scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_NPC_peak.sh # 3486732 ok
sbatch scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak.sh # 3609625 ok

#### Filter out cluster 1 / 2 and 3 / 4 and plot heatmaps
awk '$13 == "cluster_1" || $13 == "cluster_2"' output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4.bed > output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_1_2.bed
bedtools sort -i output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_1_2.bed > output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_1_2_sort.bed


awk '$13 == "cluster_3" || $13 == "cluster_4"' output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4.bed > output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_3_4.bed
bedtools sort -i output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_3_4.bed > output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_3_4_sort.bed

sbatch scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_cluster_1_2_peak.sh # 3614310 ok
sbatch scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_cluster_3_4_peak.sh # 3614314
sbatch scripts/matrix_peak_5kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_cluster_3_4_peak.sh # 3614380 ok
sbatch scripts/matrix_peak_10kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_cluster_3_4_peak.sh # 3614369


# 2. List of differential peak from THOR method
## Files location
002__ChIPseq/output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/*bw
002__ChIPseq/output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25.bed

### Generate gain and lost bed files (filter the 16th FC column (eg. >1 more in NPC))
awk -F'\t' '$16 > 1' output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25.bed > output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25_positive.bed
awk -F'\t' '$16 < 1' output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25.bed > output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25_negative.bed


## deepTools plot
conda activate deeptools
sbatch scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25.sh # 3489952 ok
sbatch --dependency=afterany:3489952 scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_kmean2.sh # 3496104 ok
sbatch --dependency=afterany:3489952 scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_kmean5.sh # 3496105 ok

sbatch scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_positive.sh # 3497014 ok
sbatch scripts/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_negative.sh # 3497040 ok


```
*NOTE: `matrix_peak_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_peak.sh` this version contain genomic coordinates of ESC peak (MACS2 qval2.3 blacklist)*



--> 1. raw MACS2 peak method is showing that most peak that exist in ESC disapear/decrease, in NPC; similarly peak that exist in NPC where not present; or less, in ESC.
----> kmean is not informative as it cluster based on peak profile; but not if different between my bigiwg (eg. cluster all the peak that are enriched upstream the center for example)
----> interestingly, large peak (2kb) decrease; while smaller peak (500bp) increase from ESC to NPC! (kmeans analysis showed it); And smaller peaks are much more abundant than larger peaks.

--> 2. list of differential peak from THOR method is overall showing a decrease of H3K27me3; clearly visible when taking them all.
----> Separating the gain and lost; we indeed see more clearly the Lost or Gain; the Gain are much less clear.

Median size (calculated with `awk '{print $3 - $2}' your_file.bed | sort -n | awk '{a[NR] = $0} END {if (NR%2==1) print a[(NR+1)/2]; else print (a[NR/2] + a[NR/2+1])/2}'`):
- Lost 1,500bp
- Gain 950bp
- ESC 1,159bp
- NPC 573bp




Let's now check **peak location to feature using ChIPseeker + GO of associated genes**; check the following:
- ESC and NPC all peak (use MACS2 unique)
- The one that gain (use THOR qval25 positive)
- The one that lost (use THOR qval25 negative)



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

# Import peaks
peaks_ESC_macs2 =  read.table('output/macs2_unique/broad_blacklist_qval2.30103/ESC_WT_H3K27me3_pool_peaks.broadPeak') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9) 
peaks_NPC_macs2 =  read.table('output/macs2_unique/broad_blacklist_qval2.30103/NPC_WT_H3K27me3_pool_peaks.broadPeak') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, score=V5, strand=V6, signal_value=V7, pvalue=V8, qvalue=V9) 
peaks_lost =  read.table('output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25_negative.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, qvalue=V15, FC=V16) 
peaks_gain =  read.table('output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25_positive.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, name=V4, strand=V6, qvalue=V15, FC=V16) 

peaks_cl1_2 =  read.table('output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_1_2_sort.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, cluster=V13) 
peaks_cl3_4 =  read.table('output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_3_4_sort.bed') %>% dplyr::rename(Chr=V1, start=V2, end=V3, cluster=V13) 


# Tidy peaks
ESC_gr = makeGRangesFromDataFrame(peaks_ESC_macs2,keep.extra.columns=TRUE)
NPC_gr = makeGRangesFromDataFrame(peaks_NPC_macs2,keep.extra.columns=TRUE)
Lost_gr = makeGRangesFromDataFrame(peaks_lost,keep.extra.columns=TRUE)
peaks_cl1_2_gr = makeGRangesFromDataFrame(peaks_cl1_2,keep.extra.columns=TRUE)
peaks_cl3_4_gr = makeGRangesFromDataFrame(peaks_cl3_4,keep.extra.columns=TRUE)


gr_list <- list(ESC=ESC_gr, NPC=NPC_gr, Lost=Lost_gr, Gain=Gain_gr)
gr_list <- list(cl_1_2=peaks_cl1_2_gr, cl_3_4=peaks_cl3_4_gr)


## Genomic Annotation ALL TOGETHER
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

### Barplot
pdf("output/ChIPseeker/annotation_barplot_WT.pdf", width=14, height=5)
pdf("output/ChIPseeker/annotation_barplot_cl1_2_3_4.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()


## Get annotation data frame
cl_1_2_annot <- as.data.frame(peakAnnoList[["cl_1_2"]]@anno)
cl_3_4_annot <- as.data.frame(peakAnnoList[["cl_3_4"]]@anno)

## Convert entrez gene IDs to gene symbols
cl_1_2_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = cl_1_2_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
cl_3_4_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = cl_3_4_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")

cl_1_2_annot$gene <- mapIds(org.Hs.eg.db, keys = cl_1_2_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
cl_3_4_annot$gene <- mapIds(org.Hs.eg.db, keys = cl_3_4_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")


## Save output table
write.table(cl_1_2_annot, file="output/ChIPseeker/annotation_deepTools_cl_1_2.txt", sep="\t", quote=F, row.names=F)
write.table(cl_3_4_annot, file="output/ChIPseeker/annotation_deepTools_cl_3_4.txt", sep="\t", quote=F, row.names=F)

# GO-associated genes
cl_1_2_annot_gene = cl_1_2_annot %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))
cl_3_4_annot_gene = cl_3_4_annot %>%
  filter(annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR"))



ego <- enrichGO(gene = as.character(cl_3_4_annot_gene$geneSymbol), 
                keyType = "SYMBOL",     # Use ENSEMBL if want to use ENSG000XXXX format
                OrgDb = org.Hs.eg.db, 
                ont = "MF",          # “BP” (Biological Process), “MF” (Molecular Function), and “CC” (Cellular Component) 
                pAdjustMethod = "BH",   
                pvalueCutoff = 0.05, 
                readable = TRUE)

pdf("output/ChIPseeker/dotplot_MF_deepTools_cl_1_2.pdf", width=7, height=7)
pdf("output/ChIPseeker/dotplot_MF_deepTools_cl_3_4.pdf", width=7, height=7)
dotplot(ego, showCategory=20)
dev.off()

enrichKEGG <- enrichKEGG(gene   = as.character(cl_3_4_annot_gene$geneId),
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")

pdf("output/ChIPseeker/dotplot_KEGG_deepTools_cl_1_2.pdf", width=7, height=5)
pdf("output/ChIPseeker/dotplot_KEGG_deepTools_cl_3_4.pdf", width=7, height=7)
dotplot(enrichKEGG, showCategory=20)
dev.off()
```

--> At ESC; much more genes than intergenic regions are H3K27me3-enriched; as compare to NPC.
----> In agreement; most of the regions that lose H3K27me3 are gene regions (= become activated from ESC to NPC)


Let's count the number of genes assigned for cl1_2 and cl3_4 with `awk -F'\t' '($column_number ~ /Promoter \(<=1kb\)|Promoter \(1-2kb\)|Promoter \(2-3kb\)|5\' UTR/){print $gene_symbol_column_number}' input.txt | sort | uniq | wc -l`

```bash
# Files
output/ChIPseeker/annotation_deepTools_cl_1_2.txt
output/ChIPseeker/annotation_deepTools_cl_3_4.txt

awk -F'\t' '($16 ~ /Promoter \(<=1kb\)|Promoter \(1-2kb\)|Promoter \(2-3kb\)|5'\'' UTR/){print $25}' output/ChIPseeker/annotation_deepTools_cl_1_2.txt | sort | uniq | wc -l # 1,438
awk -F'\t' '($16 ~ /Promoter \(<=1kb\)|Promoter \(1-2kb\)|Promoter \(2-3kb\)|5'\'' UTR/){print $25}' output/ChIPseeker/annotation_deepTools_cl_3_4.txt | sort | uniq | wc -l # 3,649
```


#### Generate deepTools-heatmap of H3K27me3 in WT/HET/KO at ESC/NPC for (H3K27me3-WT + TC-regulated) clustered genes

The idea here is check whether one of our bigiwg normalization is in agreement with the clustered genes (H3K27me3-regulated in WT + TC-significant between genotypes):
- Convert list of clustered genes to gtf
- Generate deepTools heatmaps using various bigwig:
  - bigwig from THOR
  - bigwig_ChIPseqSpikeInFree (WT in agreement with THOR-diff peaks)
  - bigwig_ChIPseqSpikeInFree_BamToBedToBigwig (WT in agreement with THOR-diff peaks)
  - bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF (WT in agreement with THOR-diff peaks)
  - bigwig_DiffBind_TMM (WT in agreement with THOR-diff peaks)
  - bigwig_DiffBind_LIB ( *NO very bad* )
  - bigwig_UniqueBamUniqueSF_DiffBind_TMM (WT in agreement with THOR-diff peaks)

***NOTE: Would be great to troubleshoot the THOR-bigwig that uses ChIPseqSpikeInFree scaling factor...***
--> Related to the NOTE; all file THOR with ChIPseqSPikeInFree are very bad (usually all peaks are LOST from ESC to NPC). One file; the one with DiffBindTMM norm looks OK (`UniqueBamDiffBindTMM`); let's see how it looks

##### HET

Generate **GTF from Bed (convert bed to gtf)**:

```bash
conda activate BedToBigwig
# File with clustered genes
../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt
## Filter the different clusters (from 1 to 5)
awk -F, '$3 == 1 {gsub("\"", ""); print $1}' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster1.txt
awk -F, '$3 == 2 {gsub("\"", ""); print $1}' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster2.txt
awk -F, '$3 == 3 {gsub("\"", ""); print $1}' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster3.txt
awk -F, '$3 == 4 {gsub("\"", ""); print $1}' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster4.txt
awk -F, '$3 == 5 {gsub("\"", ""); print $1}' ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering.txt > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster5.txt


## Filter the gtf to keep only gene from cluster (order preserved; so no need sort; unique)
grep -F -f ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster1.txt meta/ENCFF159KBI.gtf > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster1.gtf
grep -F -f ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster2.txt meta/ENCFF159KBI.gtf > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster2.gtf
grep -F -f ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster3.txt meta/ENCFF159KBI.gtf > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster3.gtf
grep -F -f ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster4.txt meta/ENCFF159KBI.gtf > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster4.gtf
grep -F -f ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster5.txt meta/ENCFF159KBI.gtf > ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster5.gtf
```

Now generate deepTools plot:

```bash

## deepTools plot
conda activate deeptools
sbatch scripts/matrix_gene_2kb_bigwig_DiffBind_TMM_HETcluster1.sh # 3822300
sbatch scripts/matrix_gene_2kb_bigwig_DiffBind_TMM_HETcluster3.sh # 3822310
sbatch scripts/matrix_gene_2kb_bigwig_DiffBind_TMM_HETcluster4.sh # 3822320
sbatch scripts/matrix_gene_2kb_bigwig_DiffBind_TMM_HETcluster5.sh # 3827747

sbatch scripts/matrix_gene_2kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_HETcluster1.sh # 3822671
sbatch scripts/matrix_gene_2kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_HETcluster3.sh # 3822902
sbatch scripts/matrix_gene_2kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_HETcluster4.sh # 3822913
sbatch scripts/matrix_gene_2kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_HETcluster5.sh # 3827755

sbatch scripts/matrix_gene_2kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_HETcluster1.sh # 3822991
sbatch scripts/matrix_gene_2kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_HETcluster3.sh # 3822999
sbatch scripts/matrix_gene_2kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_HETcluster4.sh # 3823013

sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamTMM_HETcluster1.sh # 3823034
sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamTMM_HETcluster3.sh # 3823037
sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamTMM_HETcluster4.sh # 3823048
sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamTMM_HETcluster5.sh # 3827801

sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster1.sh # 3845499
sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster2.sh # 3845501
sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster3.sh # 3845503
sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster4.sh # 3845504
sbatch scripts/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster5.sh
```

--> Conclusion:
  - bigwig from THOR; works for cluster 3 only, not with cluster1 and 4. But only one that works for cluster5
  - bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF (WT in agreement with THOR-diff peaks); do not work (huge increase of H3K27me3 at NPC for all; mess up everythings)
  - bigwig_DiffBind_TMM (WT in agreement with THOR-diff peaks); work with cluster1/3, not with cluster4. Show much more H3K27me3 at ESC in HET. Do not work for cluster5
  - bigwig_UniqueBamUniqueSF_DiffBind_TMM (WT in agreement with THOR-diff peaks) work with cluster1/3, not with cluster4. Show much more H3K27me3 at ESC in HET. Do not work for cluster5

----> bigwig_DiffBind_TMM and bigwig_UniqueBamUniqueSF_DiffBind_TMM; perform the best, the more in agreement with gene expression. However; profile HET at ESC looks a bit weird; it s like super H3K27me3 as compare to WT... Like everywhere!
------> If we are able to show that there is overall MORE H3K27me3 in HET at ESC, by WB?; we re good!

----> After more consideration; THOR is the only option in agreement with cluster5; overall it seems to perform better...

--> For the WT `THOR_ESCvsNPC_UniqueBamDiffBindTMM` seems to work pretty well. Let's see how it goes with the HET! We could use THOR TMM-method to find diff. peaks; and then use THOR normalize with ChIPseqSpikeInFree-TMM normalize scaling factors to generate the bigwig for all files.
----> Seems to work well for HET as well; even though look overall VERY more H3K27me3 at ESC; not sure that is relevant; but lets say yes lol...


































