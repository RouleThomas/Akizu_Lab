# Import meta (genome) files
Files format to follow (ENCODE):\
**chr**:
`chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrM,chrX,chrY`\
**gene**:
`ENSG00000228253.1`

Files download from [ENCODE](https://www.encodeproject.org/data-standards/reference-sequences/) and cp to cluster: 
```bash
cp /home/roulet/tsclient/roule/Google\ Drive\ Streaming/Shared\ drives/akizulab/Personal\ Folders/Thomas/meta/* \
/scr1/users/roulet/Akizu_Lab/Master/meta`
```
ENCODE genome files in `/scr1/users/roulet/Akizu_Lab/Master/meta`:
- **GRCh38 fasta genome** = GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
- **GRCh38 gtf file** = ENCFF159KBI.gtf
- **hg19 fasta genome** = male.hg19.fasta
- **hg19 gtf file** = gencode.v19.annotation.gtf

# pipeline used from previous analyses
Infos from `Method_RNAseq, DEG, volcano, GSEA, heatmap.SZ.docx`:
- Raw fastq cleaned with *fastp* (adaptor removed, >10% poly-N sequences removed, low quality removed)
	- Q20, Q30, and GC content of the clean data were calculated.
    - Downstream analyses on the clean data with high quality.
		- Seems to be the default *fastp* parameter
 - Reads map on hg19 with *STAR* 
 - Reads count on gene feature with *featureCounts*
 - DEG with *DESEq2* (padj<0.05)
 - GSEA performed and plotted with *ClusterProfiler*


# Import files from Google drive to the cluster
##### 20230308, 20230309, 20230310
I cannot use a bash script to import files as if disconnection the transfer will fail. So cp the old-fashion way, let computer running o/n.\
**ESC, NPC, 2 days-neurons**
```bash
cp /home/roulet/tsclient/roule/Google\ Drive\ Streaming/Shared\ drives/akizulab/Primary\ Data/RNAseqs/EZH1\ RNAseq/1\ and\ 2\ month\ neuron\ RNAseq\ Aug2022/01.RawData/* \
/scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq/input
``` 
**1, 2 months-neurons**
```bash
cp -r /home/roulet/tsclient/roule/Google\ Drive\ Streaming/Shared\ drives/akizulab/Primary\ Data/RNAseqs/EZH1\ RNAseq/1\ and\ 2\ month\ neuron\ RNAseq\ Aug2022/01.RawData/ \
/scr1/users/roulet/Akizu_Lab/001_EZH1_Project/001__RNAseq/input
``` 
# File renaiming
##### 20230309, 20230310
Made a custom bash script to rename each files (files are already compressed so I modified script `organize_raw.sh` to keep only renaming function). 
```bash
# example for 1 file:
outdir="input"

x="NPC_WT_R1_1"
raw_f="P_WT_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"

if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi
# Run command time-per-time (ESC, then 2dN, then PNC):
sbatch rename_raw_ESC.sh
sbatch rename_raw_2dN.sh
sbatch rename_raw_NPC.sh
```

# Quality control with FASTQC on raw fastq
##### 20230310
FASTQC is not an available module. Let's download it:
```bash
# Download in Master/Software/
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
# Unzip
unzip fastqc_v0.12.1.zip
# Add execution right to the file
chmod +x FastQC/fastqc
# Add shortcut so that we only use *fastqc* to run it
## backup .bashrc file in case
cp ~/.bashrc ~/.bashrc.backup
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software/FastQC
# Restart terminal
```
Run fastqc
```bash
# example for 1 file:
fastqc -o output/fastqc input/ESC_HET_R1_1.fq.gz

# Run time-per-time (ESC, then 2dN, then PNC):
sbatch fastqc_raw_ESC.sh # 10979372 complete
sbatch fastqc_raw_2dN.sh # 10980363 complete
sbatch fastqc_raw_NPC.sh # 10982864 complete
sbatch fastqc_raw_4wN.sh # 10983894 complete
sbatch fastqc_raw_8wN.sh # 10984036 complete
# By mistake, I replace output/fastqc by the fastp-fastqc. 
# Let's re-generate fastqc for raw files in output/fastqc/raw
sbatch fastqc_raw_ESC.sh # 11062004
sbatch fastqc_raw_2dN.sh # 11062006
sbatch fastqc_raw_NPC.sh # 11062007
sbatch fastqc_raw_4wN.sh # 11062008
sbatch fastqc_raw_8wN.sh # 11062009
```

# Quality control with FASTP (trim)
##### 20210310
Install [fastp](https://github.com/OpenGene/fastp).
```bash
# Download in Master/Software/
wget https://opengene.org/fastp/fastp
chmod a+x ./fastp
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software
# Restart terminal
```
Run fastp
```bash
# example for 1 file:
fastp -i input/ESC_WT_R1_1.fq.gz -I input/ESC_WT_R1_2.fq.gz \
      -o output/fastp/ESC_WT_R1_1.fq.gz -O output/fastp/ESC_WT_R1_2.fq.gz \
	  -h output/fastp/ESC_WT_R1 -j output/fastp/ESC_WT_R1

# Run time-per-time (ESC, then 2dN, then PNC):
sbatch scripts/fastp_raw_ESC.sh # 10980247 complete (important infos)
sbatch scripts/fastp_raw_NPC2dN.sh # 10984692 complete (important infos)
sbatch scripts/fastp_raw_4wN8wN.sh # 10984925 complete (important infos)
sbatch scripts/fastp_raw_8wN_miss.sh # 11033843 complete (important infos)
```
Run fastqc on fastp-trimmed files
```bash
sbatch scripts/fastqc_fastp.sh # 11034677 (cancelled "due to time limit")
sbatch scripts/fastqc_fastp_1.sh # 11060976 complete
```

# Mapping with STAR
##### 20230310, 20230313, 20230314
Data are unstranded.
## Index the genome
*NOTE: theorically optimal size for `--sjdbOverhang` is [max read length]-1, thus need create specific index for specific read size. But the effect is marginal according to the [creator](https://github.com/alexdobin/STAR/issues/931). So let's keep it default.*

hg19 genome with 12CPU and 50G mem (time=<1.5day)
```bash
module load STAR/2.7.3a-GCC-9.3.0
# command
STAR --runThreadN 12 \
	--runMode genomeGenerate \
	--genomeDir /scr1/users/roulet/Akizu_Lab/Master/meta/STAR_hg19 \
	--genomeFastaFiles /scr1/users/roulet/Akizu_Lab/Master/meta/male.hg19.fasta \
	--sjdbGTFfile /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf 

# Run in slurm
sbatch STAR_index_hg19.sh # 10982789 complete
```


### Untrimmed fastq
Keep standard parameter as adapted for mammalian genome. Some examples [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html) and [here](https://biocorecrg.github.io/RNAseq_course_2019/alnpractical.html).

Mapping
```bash
module load STAR/2.7.3a-GCC-9.3.0
# example for 1 file:
STAR --genomeDir ../../Master/meta/STAR_hg19/ \
	--runThreadN 12 \
	--readFilesCommand zcat \
	--readFilesIn input/NPC_WT_R1_1.fq.gz input/NPC_WT_R1_2.fq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix output/STAR/NPC_WT_

# Run time-per-time:
sbatch scripts/STAR_raw_NPC.sh # 11034673, slight test, output ok
sbatch scripts/STAR_raw_NPC_1.sh # 11061777, 11065401
sbatch scripts/STAR_raw_ESC.sh # 11061825, 11065404
sbatch scripts/STAR_raw_2dN.sh # 11061828, 11065418
sbatch scripts/STAR_raw_4wN.sh # 11061917, 11065420
sbatch scripts/STAR_raw_8wN.sh # 11061920, 11065422
```
### Fastp-trimmed fastq
Mapping
```bash
module load STAR/2.7.3a-GCC-9.3.0
# example for 1 file:
STAR --genomeDir ../../Master/meta/STAR_hg19/ \
	--runThreadN 12 \
	--readFilesCommand zcat \
	--readFilesIn output/fastp/NPC_WT_R1_1.fq.gz output/fastp/NPC_WT_R1_2.fq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix output/STAR/fastp/NPC_WT_R1_

# Run time-per-time:
sbatch scripts/STAR_fastp_ESC.sh # 11062736, 11065467
sbatch scripts/STAR_fastp_NPC.sh # 11062742, 11065493
sbatch scripts/STAR_fastp_2dN.sh # 11062743, 11065515
sbatch scripts/STAR_fastp_4wN.sh # 11062744, 11065547
sbatch scripts/STAR_fastp_8wN.sh # 11062745, 11065561
```
Mapping indexation
```bash
# example for 1 file:
module load sam-bcf-tools
samtools index output/STAR/NPC_WT_Aligned.sortedByCoord.out.bam

# time per time (raw and fastp together):
sbatch STAR_index_NPC.sh # 11075019
sbatch STAR_index_ESC.sh # 11075041
sbatch STAR_index_2dN.sh # 11075067
sbatch STAR_index_4wN.sh # 11075079
sbatch STAR_index_8wN.sh # 11075108
```
*NOTE: next time do indexation after mapping; same script*

--> Let's compil the number of uniquely mapped reads for all files (add it in the Google Drive `RNAseq_infos.xlsx` file)
```bash
# Print nb of uniq map reads for raw mapping
for file in output/STAR/raw/*Log.final.out; do
    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $file | awk '{print $NF}')
    echo "$file: Number of uniquely mapped reads: $uniquely_mapped_reads"
done > output/STAR/raw/uniq_map_reads_counts_raw.txt

# Print nb of uniq map reads for fastp mapping
for file in output/STAR/fastp/*Log.final.out; do
    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $file | awk '{print $NF}')
    echo "$file: Number of uniquely mapped reads: $uniquely_mapped_reads"
done > output/STAR/fastp/uniq_map_reads_counts_fastp.txt
```
**More number of uniquely mapped reads for the fastp trimmed reads**, thus from now on, data analyses with the fastp-trimmed.
# Install IGV for vizualization
Go to [IGV](https://software.broadinstitute.org/software/igv/download) and copy link for Linux download.\
Go to `Master/software`
```bash
# download and install
wget https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_2.16.0_WithJava.zip
unzip IGV*
# Add shortcut so that we only use *igv.sh* to run it
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software/IGV_Linux_2.16.0
# Restart terminal
```

# Install Anaconda
```bash
module load Python/3.9.6*
```
Go to `Master/software` 
```bash
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh`
bash Anaconda3-2022.10-Linux-x86_64.sh
# Follow installation; accept all
```
Anaconda3 installed to `/home/roulet/anaconda3`; no need to module load it.
# Generate coverage (wiggle) files
Create a **deeptools; conda environment**
```bash
conda create -n deeptools -c bioconda deeptools
```
Generate few files RPKM-normalized for comparison with novogene analyses:
- 4wN_WT_R1 (H9)
- 4wN_KO_R1 (8del)
- 4wN_HET_R1 (het5)
```bash
conda activate deeptools
# example for 1 file:
bamCoverage --bam output/STAR/fastp/4wN_WT_R1_Aligned.sortedByCoord.out.bam \
	--outFileName output/temp/4wN_WT_R1_Aligned.sortedByCoord.out.bigwig \
	--outFileFormat bigwig \ 
	--normalizeUsing RPKM \ 
	--binSize 10
# Run for all 3 files:
sbatch RPKM_wig.sh # 11075632
```
Generate coverage files **TPM-normalized** (=BPM) for all samples (fastp-trimmmed reads)

*NOTE: RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb)). BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions)*

```bash
conda activate deeptools
# example for 1 file with correct BPM-TPM normalization:
bamCoverage --bam output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam \
	--outFileName output/bigwig/${x}_Aligned.sortedByCoord.out.bw \
	--outFileFormat bigwig \
	--normalizeUsing BPM \
	--binSize 10   
# run time-per-time:
sbatch scripts/TPM_bw_NPC.sh # 11075894
sbatch scripts/TPM_bw_ESC.sh # 11075901
sbatch scripts/TPM_bw_2dN.sh # 11075962
sbatch scripts/TPM_bw_4wN.sh # 11075965
sbatch scripts/TPM_bw_8wN.sh # 11075967
```

# Count the reads on gene feature
##### 20230314
Create a **featurecounts; conda environment**
```bash
conda create -c bioconda -n featurecounts subread
conda activate featurecounts
```

Let's first do a comparison raw vs fastp to confirm we end up with more counts using fastp
```bash
# ESC_WT_R1 raw:
featureCounts -p -C -O \
	-P -B -d 30 -D 1000 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/featurecounts/ESC_WT_R1.txt output/STAR/raw/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 53%
# ESC_WT_R1 fastp:
featureCounts -p -C -O \
	-P -B -d 30 -D 1000 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/temp/ESC_WT_R1.txt output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 52%
```
--> fastp is again better.

Many fragment are unassigned for Fragment_Lenght. Samples not enough sonicated maybe? Let's try some fine-tuning:

```bash
# Both ends mapped required but fragment short
featureCounts -p -C -O \
	-P -B -d 50 -D 600\
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/temp/ESC_WT_R1.txt output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 46.9%
# Both ends mapped required but fragment long
featureCounts -p -C -O \
	-P -B -d 30 -D 2500 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/temp/ESC_WT_R1.txt output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 66.6%
# Both ends mapped required but fragment even long
featureCounts -p -C -O \
	-P -B -d 30 -D 10000 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/temp/ESC_WT_R1.txt output/STAR/fastp/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 81.4%
# Both ends mapped not required 
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	o output/featurecounts/ESC_WT_R1.txt output/STAR/raw/ESC_WT_R1_Aligned.sortedByCoord.out.bam
## 87.2%
```
--> Seems that the longer fragment length allowed the better it is. (expected as no sonication in RNAseq; even though library prep focus for 300bp fragment)
- `-O` to count on meta feature (gene)
- `-p` (paired-end) with `-C` (not count paired-reads on 2 diff. chr.)
- NOT THIS: `-P -B -d 30 -D 1000` (set min and max paired reads to 30 to 1000bp)

Count on gene features with parameter
```bash
# example for 1 file:
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/featurecounts/${x}.txt output/STAR/raw/${x}_Aligned.sortedByCoord.out.bam

# time per time:
sbatch scripts/featurecounts_ESC.sh # 11076424
sbatch scripts/featurecounts_NPC.sh # 11076424
sbatch scripts/featurecounts_2dN.sh # 11076427
sbatch scripts/featurecounts_4wN.sh # 11076427
sbatch scripts/featurecounts_8wN.sh # 11076430
```
# Quality control metrics
Print number of succesfully assigned alignments for each sample (add to drive `RNAseq_infos.xlsx`)
```bash
for file in output/featurecounts/*.summary; do
    assigned_reads=$(grep "Assigned" $file | awk '{print $NF}')
    echo "$file: Assigned: $assigned_reads"
done > output/featurecounts/assigned_reads_counts.txt
```
Print the total number of reads
```bash
for file in output/STAR/fastp/*.final.out; do
    input_reads=$(grep "Number of input reads" $file | awk '{print $NF}')
    echo "$file: Number of input reads: $input_reads"
done > output/STAR/fastp/input_reads_counts.txt
```
Add these values to the `RNAseq_infos.xlsx`\
Then in R; see `/home/roulet/001_EZH1_project/001_EZH1_project.R`.

--> Overall >80% input reads as been assigned to (gene) features

# Install Bioconductor
Create a **deseq2; conda environment**
```bash
conda create -n deseq2 -c bioconda bioconductor-deseq2
```

see help here: "port all count sample (WT Rep 1 and 3"
https://github.com/RouleThomas/Wagner_Lab/blob/main/Tian_2022_RNAseq.md 









