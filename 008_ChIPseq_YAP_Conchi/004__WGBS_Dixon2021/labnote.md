# Project and goals 

Re-analysis of WGBS (Whole Genome Bisulfit Sequencing) dataset from (Dixon et al)[10.1126/science.abd0875].

--> 2 rep WT H1 hESC


# Download data


- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject  PRJNA631028 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 18603395 xxx
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


Installation through Docker, let's try it! XXX


























XXXXXXXX HERE





# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 18390922 ok
sbatch scripts/fastp_QSER1.sh # 18544351 ok
```
