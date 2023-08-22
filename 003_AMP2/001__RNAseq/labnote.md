# Project

Function of AMP2 in mice.
- 3 organs: CB cerebellum, CT cortex, HP hippocampus

# Reference genome

Reference genome and annotation downloaded from [ENCODE](https://www.encodeproject.org/data-standards/reference-sequences/):
- mm10 XY reference genome (ENCODE3 used only one reference genome for analysis)
- mm10 GENCODE VM21 merged annotations gtf file 

--> Genome Files in `/Akizu_Lab/Master/meta_mice`

# FASTQC on raw

```bash
sbatch scripts/fastqc_raw_CB.sh # 4611595
sbatch scripts/fastqc_raw_CT.sh # 4611594
sbatch scripts/fastqc_raw_HP.sh # 4611593
```

--> XXX

# Quality control with FASTP (trim)

```bash
sbatch scripts/fastp_CB.sh # 4611647
sbatch scripts/fastp_CT.sh # 4611648
sbatch scripts/fastp_HP.sh # 4611649

```

Run fastqc on fastp-trimmed files


