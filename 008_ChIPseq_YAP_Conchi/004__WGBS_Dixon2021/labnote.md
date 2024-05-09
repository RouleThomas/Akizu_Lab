# Project and goals 

Re-analysis of WGBS (Whole Genome Bisulfit Sequencing) dataset from (Dixon et al)[10.1126/science.abd0875].

--> 2 rep WT H1 hESC


# Download data


- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject  PRJNA631028 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 18603395 ok
```

--> Seems like no input samples needed for WGBS (indeed it is just presence absence of mC)



## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename_008004.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_008004.txt
```

--> All good 



# Data analysis

## msPIPE 

Let's try [msPIPE](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04925-2) for the analysis. [Github](https://github.com/jkimlab/msPIPE). Seems it do everyhting from preprocessing to downstream analysis! *Paper from 2022 so quite recent*


### Run through Docker

#### Install Docker

Need to do sudo but I do not have privilege...


#### Install in a conda env

Below the requirements:
- R >= 4.1.3, python2 and python3
- Trim Galore
- Samtools
- cutadapt
- bowtie, bowtie2
- Bismark
- BS-Seeker2
- Other required R packages can be installed via msPIPE/bin/script/Package_install.R.


Create a conda env msPIPE

```bash
# attempt1
conda create -n msPIPE r-base=4.1.3 python

conda activate msPIPE
conda install -c conda-forge/label/cf202003 python=2.7 # fail here, impossible to install 2 python!!

conda env remove --name msPIPE

# attempt2
conda create -n msPIPE r-base=4.1.3 python=3.8
conda install trim-galore samtools cutadapt bowtie bowtie2 bismark bs-seeker2 # fail here

# attempt3 - copy bowtie2 env and install the remaining stuff

conda create --name msPIPE --clone bowtie2
conda activate msPIPE
conda install -c bioconda -c conda-forge -c anacond trim-galore # fail
conda install trim-galore 


FAIL AGAIN!!! XXXXXXXXXXX HERE

```
























XXXXXXXX HERE





# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 18390922 ok
sbatch scripts/fastp_QSER1.sh # 18544351 ok
```
