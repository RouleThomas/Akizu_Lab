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
sbatch scripts/download_urls.sh # 13293911 XXX

```

















