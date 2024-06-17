# Background

BAP1 interact with some member of PRC2/1. Notably involved in histone de-ubiquitination (chromatin opening; as H2A Ub = compact).
BAP1 part of a big complex. KO is embryo lethal. 11 patients with HT BAP1 mutation = NDD / ASD syndrome --> 1.5 Fold more ubuiquitnation H2AK119Ub (as BAP1 not functional). 
Generate mouse model using CRE system to deplete BAP1. Strong phenotype in mice, hyperactive and burying assay (they do not burry anything, super strong behaviorial phenotype!) + dentate gyruse (part of the hippocampal formation in the temporal lobe of the brain) decrease along dev (CC3 expression (cell death marker) increase so could be a partial contributor of decreased size but not only that..)
During dentate gyrus development, Pax6 > Tbr2 > Prox1 gene expression switch (express + cell migrate) = marker of development (could be used as proxy for pseudotime analysis!). Found that Pax6 and Tbr2 decrease in the Mt (at e18). 

# Project

Casey Lim, MD/PhD student did the experiment.
Lab of Seon-Hee Kim 

WT vs *Bap1KO* mouse scRNA/ATACseq

BAP1 involved in NDD



# Docs

Integration of **scRNAseq/ATACseq** tutorial/infos [here](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/scAtacSeq/lab-sc_atac_seq.html) and [here](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette) and [here](https://satijalab.org/seurat/articles/multimodal_vignette.html). Also interesting discussion [here](https://github.com/satijalab/seurat/issues/5346):

- Start with RNA and perform QC
- Do cell clustering on RNA
- Integrate the ATAC data



What I could do 1st:
--> Use scRNASeq data to annotate cell types and then integrate ATACseq information
--> Identify motif from the ATACseq peak
--> DEGs WT vs KO
--> Diff. Acc. Regions WT vs KO; annotate diff. peak to genes (nearest TSS)
--> Then meeting from there what we do next

Next:
- Which promoters and enhancers become active in different cell types and conditions?
- Check if any TF binding sites are active in different cell types and conditions?
- Are some genes primed for expression, (e.g. the promoters show an open chromatin state, but the gene is not expressed yet)?



# file

- Transfer data from Hard drive to Google Drive then to HPC Cluster (transfer the aggregated cellranger counts output (`input_raw/A1B1_A2B2`) and the raw fastq (`input_raw/F001 and F002`))

A1 is WT snATACseq = ATAC_WT
A2 is Bap1 cKO snATACseq = ATAC_Bap1KO
B1 is WT snRNAseq = RNA_WT
B2 is Bap1 cKO snRNAseq = RNA_Bap1KO




## Counting with cellranger count

Within each folder in `/snRNAseq_Kcnc1_R320H/snRNAseq_Kcnc1_reorganized/*` I have two lanes L001 and L002 with I1/I2 and R1/R2 fastq. 

- *Option1*: Count separately the scRNAseq / scATACseq data; I need individual scRNAseq count file to use scrublet (doublet) and soupX (RNA contamination)
- *Option2*: Count together, using multiome kit special command (`cellranger-arc count`); let's see whether I can individualize the scRNAseq data for QC...; info [here](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/single-library-analysis)

--> Prefer option1


### Install cellranger-atac and -arc

For both options, need to install [cellranger-arc](https://kb.10xgenomics.com/hc/en-us/articles/360059656912-Can-I-analyze-only-the-Gene-Expression-data-from-my-single-cell-multiome-experiment) and [cellranger-atac](https://support.10xgenomics.com/single-cell-atac/software/downloads/latest) 



```bash
# Cellranger-arc
cd /scr1/users/roulet/Akizu_Lab/Master/software

curl -o cellranger-arc-2.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.0.2.tar.gz?Expires=1718694544&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hcmMvY2VsbHJhbmdlci1hcmMtMi4wLjIudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNzE4Njk0NTQ0fX19XX0_&Signature=Uy6F7La3bYC9yArvNY0PVBXAkD~dvozM53LBq3hwIlAuhtJBXqIjkm9IZ18YH8Mm95DrHrLJlepWamgY9YXR8jKBc6Ku16UEDLHEtakCr2oDdyQVTH3kjxBsiLt5vyu4CmtFqyzIyfSmOSK6bkg1V12J~6x-MI6O0z00f-io-oGFRzGTGYUZ0Fap-EerqDzbzBecazFz0WVxhXkk5MWVy3dz9xTNQCij7B5Ebv0ZsaNp4OLQ09WgjV5l938n32QfYzmth08kO3IPgVQd24dIMEYUnmBNC89d55S4hGb5cwmia5q7lhSW~GawhYpIqve5oAN98Cwu-mhOuu2zDTNs6A__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar -zxvf cellranger-arc-2.0.2.tar.gz

## add cellranger to our PATH
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software/cellranger-arc-2.0.2
## Restart terminal
which cellranger-arc



# cellranger-atac
cd /scr1/users/roulet/Akizu_Lab/Master/software

curl -o cellranger-atac-2.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-2.1.0.tar.gz?Expires=1718695272&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hdGFjL2NlbGxyYW5nZXItYXRhYy0yLjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE3MTg2OTUyNzJ9fX1dfQ__&Signature=D6x9msRKD308ae2Ylu9KheiWkE0CKmZGlNpKzA2WmH4~UrdR9v9TOjRsOqzl5j4YNhkY~H5HeS3tgfzuIFFMkgacXKEXiSDCFc~ex1pEvu5xplYS0ELQ6cOQaR0Aivy6vFGGMYrvelFeoK5w~22Tw~TtsHzvO3EoJRjKDRdzkKHqBdVy8BzOYnjXauHW448~lCyvcd0rlm-CDqgn~uRA6yijRZjTWTJgUFdwOPkAtxuAaL5lEYI1Zia0bl5mHh7tzenKzJlx~0tFEtukf8u07xbNnR9-QJsJQ-eVj4hPbAOnqUo5BnKjuCHJfneSwHnv8cQN7dyN-e5za9onAiNkxg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar -zxvf cellranger-atac-2.1.0.tar.gz

## add cellranger to our PATH
nano ~/.bashrc # add: export PATH=$PATH:/scr1/users/roulet/Akizu_Lab/Master/software/cellranger-atac-2.1.0
## Restart terminal
which cellranger-atac


# Add genome data
## cellranger-arc
cd 002_scRNAseq/meta
curl -O https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz # human for later just in case
curl -O https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz # mice for this project

tar -zxvf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
tar -zxvf refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz


## cellranger-atac
#--> TRY USING THE arc ones, looks like file are the same
# curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz 
# curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz

```





### Count

--> Decide later with Seurat, which option is better to use. 

```bash 
conda activate scRNAseq
which cellranger

# Run count using mice genome



## For option1
sbatch scripts/cellranger_count_RNA_WT.sh # 20409628 fail chemistry; 20409800 xxx
sbatch scripts/cellranger_count_RNA_Bap1KO.sh # 20409629 fail chemistry; 20409813 xxx
sbatch scripts/cellranger_count_ATAC_WT.sh # 20409600 fail chemistry; 20409818 fail need cellranger-atac; 20414461 xxx
sbatch scripts/cellranger_count_ATAC_Bap1KO.sh # 20409604 fail chemistry; 20409822 fail need cellranger-atac; 20414466 xxx

## For option2
# --> To do if option1 fail
#XXX sbatch scripts/cellranger_count_RNA_ATAC_WT.sh # 
#XXX sbatch scripts/cellranger_count_RNA_ATAC_Bap1KO.sh #  

```


--> As sc data generated with a multiome kit:
    - need to add `--chemistry ARC-v1` in counting for the RNA; solution found [here](https://bioinformatics.stackexchange.com/questions/18186/10x-low-rate-of-correct-barcodes-was-observed-for-the-candidate-chemistry-choice)
    - need to use `cellranger-atac count` to count ATAC exp solely; solution found [here](https://kb.10xgenomics.com/hc/en-us/articles/360061165691-Can-I-analyze-only-the-ATAC-data-from-my-single-cell-multiome-experiment)






