# pipeline used from previous analyses
Seems Shuo used same method as for the ChIP: 
- bowtie2: `-q --local --no-mixed --no-unal --dovetail --phred33 -x`
- samtools: `-bS -q 20 -f 0x2`

[Here](https://hbctraining.github.io/Intro-to-ChIPseq-flipped/CUT&RUN/Intro_to_CUT&RUN.html) different parameters are used. But whatever works is ok!

Compared to ChIP-seq experiments, CUT&RUN data has low background and low sequence depth, making peak calling more susceptible to false positives. Thus, peak calling from CUT&RUN data requires high specificity for true positive peaks instead of high sensitivity. Recommended not to use MACS2 but [SEACR](https://github.com/FredHutch/SEACR)
# Import files from Google drive to the cluster

Using `cp`; check integrity (volume) of files!! (`ls -lh`)

# File renaiming

**Move all fastq files within the input folder** (now [file] are in input/folder1/[file]) using a loop:
```bash
for file in input/*/;
    do mv "$file"* input/;
done

rmdir input/*
```
- `input/*/` this go 1 folder downstream; whatever their names
- we move file to `input/` directory 

**Rename files**
Detail about file renaiming can be found in `CutRun_infos.xlsx` in my Google drive folder.


```bash
# example for 1 file:
outdir="input"

x="R1_Het5_Igg_1"
raw_f="R1_Het5_Igg_CKDL220021572-1A_HF3THDSX5_L1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"

if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# Run slurm for all:
sbatch rename_raw.sh # 11452878 ok
```


# Spike-in quality control

Script from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) with quality control check has been modified to work with zipped fastq files and adapted to my specific nomenclature; also now specified samples analyzed...
```bash
sbatch scripts/SNAP-CUTANA_K-MetStat_Panle_ShellScript.sh # 11454380
```
Then use the xlsx file from [EpiCypher](https://www.epicypher.com/products/nucleosomes/snap-cutana-k-metstat-panel) to generate quality control plot.


# Fastp trimming and fastqc
## Fastqc on raw reads
```bash
sbatch scripts/fastqc_WT.sh # 11454608
sbatch scripts/fastqc_HET.sh # 11454610
sbatch scripts/fastqc_KO.sh # 11454609
sbatch scripts/fastqc_patient.sh # 11454604
```

## Fastp to clean reads
```bash
sbatch scripts/fastp_WT.sh # 11454715
sbatch scripts/fastp_HET.sh # 11454712
sbatch scripts/fastp_KO.sh # 11454714
sbatch scripts/fastp_patient.sh # 11454699
```

## Fastqc on clean reads
```bash
sbatch scripts/fastqc_fastp_WT.sh # 
sbatch scripts/fastqc_fastp_HET.sh # 
sbatch scripts/fastqc_fastp_KO.sh # 
sbatch scripts/fastqc_fastp_patient.sh # 
```
