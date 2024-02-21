# Project and goals 

Re-analysis of RNAseq dataset from (Ciceri et al)[https://www.nature.com/articles/s41586-023-06984-8].

- Generate TPM expression per genes
--> Add this to RNAseqDataViewer shiny App
- Do DEGs, time point per time point



# Download data


- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject PRJNA803216 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 13293911 ok

```

## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename.txt
```

--> All good 




# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 13366413 xxx
```


## mapping fastp trim

```bash
sbatch --dependency=afterany:13366413 scripts/STAR_mapping_fastp.sh # 13367330 xxx
```

--> Around xx% aligned reads



# Count with featureCounts


--> **Whatever that is stranded or not, lets count with non stranded parameters to have same parameter as our data**


Count on gene features with parameter
```bash
conda activate featurecounts

# all samples:
sbatch --dependency=afterany:13367330 scripts/featurecounts.sh # 13368068 xxx
```

--> More than XXX% of succesfully assigned alignments



## Calculate TPM and RPKM

XXXXXXXXX Modify this below if all GOOD before





Use custom R script `RPKM_TPM_featurecounts.R` as follow:
```bash
conda activate deseq2
# Rscript scripts/RPKM_TPM_featurecounts.R INPUT OUTPUT_PREFIX
sbatch scripts/featurecounts_TPM.sh #
# mv all output to output/tpm or rpkm folder
mv output/featurecounts/*tpm* output/tpm/
mv output/featurecounts/*rpkm* output/rpkm/
```

All good. 













