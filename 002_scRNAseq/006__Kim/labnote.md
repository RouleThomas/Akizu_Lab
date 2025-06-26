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




# Counting with cellranger count

Within each folder in `/snRNAseq_Kcnc1_R320H/snRNAseq_Kcnc1_reorganized/*` I have two lanes L001 and L002 with I1/I2 and R1/R2 fastq. 

- *Option1*: Count separately the scRNAseq / scATACseq data; I need individual scRNAseq count file to use scrublet (doublet) and soupX (RNA contamination)
- *Option2*: Count together, using multiome kit special command (`cellranger-arc count`); let's see whether I can individualize the scRNAseq data for QC...; info [here](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/single-library-analysis)

--> Prefer option1


## Install cellranger-atac and -arc

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





## Count

- **Option1**: Count RNA and ATAC separately
- **Option2**: Count RNA and ATAC together with [cell ranger arc](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/single-library-analysis)

```bash 
conda activate scRNAseq
which cellranger

# Run count using mice genome



## For option1
sbatch scripts/cellranger_count_RNA_WT.sh # 20409628 fail chemistry; 20409800 ok
sbatch scripts/cellranger_count_RNA_Bap1KO.sh # 20409629 fail chemistry; 20409813 ok
sbatch scripts/cellranger_count_ATAC_WT.sh # 20409600 fail chemistry; 20409818 fail need cellranger-atac; 20414461 too long; re run more ressource (750go) 20631094 ok
sbatch scripts/cellranger_count_ATAC_Bap1KO.sh # 20409604 fail chemistry; 20409822 fail need cellranger-atac; 20414466 too long; re run more ressource 20631095 ok

## For option2
# --> To do if option1 fail; option1 is not good as it make cell name different between RNA and ATAC
### Run cellranger arc count
sbatch scripts/cellranger_count_multiome_WT.sh #  23546723 ok
sbatch scripts/cellranger_count_multiome_Bap1KO.sh # 23546793 ok 

```


--> As sc data generated with a multiome kit:
    - need to add `--chemistry ARC-v1` in counting for the RNA; solution found [here](https://bioinformatics.stackexchange.com/questions/18186/10x-low-rate-of-correct-barcodes-was-observed-for-the-candidate-chemistry-choice)
    - need to use `cellranger-atac count` to count ATAC exp solely; solution found [here](https://kb.10xgenomics.com/hc/en-us/articles/360061165691-Can-I-analyze-only-the-ATAC-data-from-my-single-cell-multiome-experiment)



--> RNA samples count succesfully

--> ATAC samples very long to count (>3 days wit 400Go mem)

--> **If multiome10x, use option2: cellcount-arc (to maintain cellname identity)!!!**
- Need run cellranger [mkfastq](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/generating-fastqs-mkfastq). Need generate CSV samplesheet to indicate fastq metrics.
  NOTE: check the L00*, if same number we are running on a single lane, so `1` need to be used.
  NOTE: *SI-TT-A1* code for *RNA*; *SI-NA-A1* code for *ATAC*

```bash
# sample information
RNA_WT (input_raw/F001/B1/)
Lane,Sample,Index
5,B1_CKDL230036293-1A_22FCTTLT3,SI-TT-A1

RNA_Bap1KO (input_raw/F001/B2/)
Lane,Sample,Index
5,B2_CKDL230036294-1A_22FCTTLT3_S2,SI-TT-A1

ATAC_WT (input_raw/F002/A1/)
Lane,Sample,Index
8,A1SubLib_CKDL230036285-1A_22FF27LT3_S1,SI-NA-A1

ATAC_Bap1KO (input_raw/F002/A2/)
Lane,Sample,Index
8,A2SubLib_CKDL230036286-1A_22FF27LT3_S2,SI-NA-A1

# csv file:
##WT /006__Kim/libraries_WT.csv
fastqs,sample,library_type
/scr1/users/roulet/Akizu_Lab/002_scRNAseq/006__Kim/input_raw/F002/A1/,A1SubLib_CKDL230036285-1A_22FF27LT3,Chromatin Accessibility
/scr1/users/roulet/Akizu_Lab/002_scRNAseq/006__Kim/input_raw/F001/B1/,B1_CKDL230036293-1A_22FCTTLT3,Gene Expression
##Bap1KO /006__Kim/libraries_Bap1KO.csv
fastqs,sample,library_type
/scr1/users/roulet/Akizu_Lab/002_scRNAseq/006__Kim/input_raw/F002/A2/,A2SubLib_CKDL230036286-1A_22FF27LT3,Chromatin Accessibility
/scr1/users/roulet/Akizu_Lab/002_scRNAseq/006__Kim/input_raw/F001/B2/,B2_CKDL230036294-1A_22FCTTLT3,Gene Expression
```



## RNA contamination and doublet detection
- doublet detection using [scrublet](https://github.com/swolock/scrublet) **on the filtered matrix**
- ambient RNA correction using `soupX` in R before generating the Seurat object


```bash
conda deactivate # base environment needed

sbatch scripts/scrublet_RNA_WT.sh # 20606254 ok
sbatch scripts/scrublet_RNA_Bap1KO.sh # 20606261 ok

```


Doublet detection score (YAP1scRNAseq ranged between 0.1 to 42%); here ~27%
--> Table in `006__Kim/QC_metrics_all.xlsx`

--> Successfully assigned doublet


Now ambiant RNA contamination correction with soupX


**SoupX RNA cleaning**


```bash
conda activate scRNAseq
```



```R
# install.packages('SoupX')
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

# soupX decontamination
## Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('RNA_Bap1KO/outs')   # CHANGE FILE NAME HERE

## Assess % of conta
pdf("output/soupX/autoEstCont_RNA_Bap1KO.pdf", width=10, height=10)   # CHANGE FILE NAME HEREs
sc = autoEstCont(sc) 
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/RNA_Bap1KO.RData") # CHANGE FILE NAME HERE

```

--> 1-3% conta; success


# Seurat analysis

Let's analyse the RNA samples 1st:
- QC
- Define cluster
- Integrate ATACseq



```bash
conda activate scRNAseqV2
```

 
```R
# install.packages('SoupX')
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

## Load the matrix and Create SEURAT object
samples <- list(
  RNA_WT = "output/soupX/RNA_WT.RData",
  RNA_Bap1KO = "output/soupX/RNA_Bap1KO.RData"
)

seurat_objects <- list()

for (sample_name in names(samples)) {
  load(samples[[sample_name]])
  seurat_objects[[sample_name]] <- CreateSeuratObject(counts = out, project = sample_name)
}

## Function to assign Seurat objects to variables (unlist the list)
assign_seurat_objects <- function(seurat_objects_list) {
  for (sample_name in names(seurat_objects_list)) {
    assign(sample_name, seurat_objects_list[[sample_name]], envir = .GlobalEnv)
  }
}
assign_seurat_objects(seurat_objects) # This apply the function


# QUALITY CONTROL
# Function to add mitochondrial and ribosomal content
add_quality_control <- function(seurat_object) {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^Rp[sl]")
  return(seurat_object)
}
seurat_objects <- lapply(seurat_objects, add_quality_control) # This apply the function
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list


# Function to add doublet information
add_doublet_information <- function(sample_name, seurat_object) {
  doublet_file <- paste0("output/doublets/", sample_name, ".tsv")
  doublets <- read.table(doublet_file, header = FALSE, row.names = 1)
  colnames(doublets) <- c("Doublet_score", "Is_doublet")
  seurat_object <- AddMetaData(seurat_object, doublets)
  seurat_object$Doublet_score <- as.numeric(seurat_object$Doublet_score)
  return(seurat_object)
}
## Apply the function to each Seurat object in the list
for (sample_name in names(seurat_objects)) {
    seurat_objects[[sample_name]] <- add_doublet_information(sample_name, seurat_objects[[sample_name]])
  }

assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list



# Loop to generate QC plots for each sample
# Define the output directory
output_dir <- "output/seurat/"

# Loop through each Seurat object to generate the plots
for (sample_name in names(seurat_objects)) {
  seurat_object <- seurat_objects[[sample_name]]
  
  # Violin plot
  pdf(paste0(output_dir, "VlnPlot_QC_", sample_name, ".pdf"), width = 10, height = 6)
  print(VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) & 
    theme(plot.title = element_text(size = 10)))
  dev.off()
  
  # Scatter plots
  pdf(paste0(output_dir, "FeatureScatter_QC_RNAcount_mt_", sample_name, ".pdf"), width = 5, height = 5)
  print(FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  dev.off()
  
  pdf(paste0(output_dir, "FeatureScatter_QC_RNAcount_RNAfeature_", sample_name, ".pdf"), width = 5, height = 5)
  print(FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  dev.off()
  
  pdf(paste0(output_dir, "FeatureScatter_QC_rb_mt_", sample_name, ".pdf"), width = 5, height = 5)
  print(FeatureScatter(seurat_object, feature1 = "percent.rb", feature2 = "percent.mt"))
  dev.off()
  
  pdf(paste0(output_dir, "FeatureScatter_QC_mt_doublet_", sample_name, ".pdf"), width = 5, height = 5)
  print(FeatureScatter(seurat_object, feature1 = "percent.mt", feature2 = "Doublet_score"))
  dev.off()
  
  pdf(paste0(output_dir, "FeatureScatter_QC_RNAfeature_doublet_", sample_name, ".pdf"), width = 5, height = 5)
  print(FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "Doublet_score"))
  dev.off()
}

## some more manually
pdf("output/seurat/VlnPlot_QC_nCount_RNA_WT.pdf")
VlnPlot(seurat_object, features = c("nCount_RNA"), ncol = 4, pt.size = 0.1, y.max = 200000) & 
    theme(plot.title = element_text(size = 10) )
dev.off()


# QC filtering _ V1

### V1 not super stringeant; mit > 10 and RNAfeature 50
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 50 & seurat_object@meta.data$QC == 'Pass', 'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 50 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 10 & seurat_object@meta.data$QC == 'Pass', 'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 50 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list




# QC filtering _ V2

### V2 not super stringeant; mit > 5 and RNAfeature 50; rb >10
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 50 & seurat_object@meta.data$QC == 'Pass', 
                                  'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 50 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', 
                                  paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 5 & seurat_object@meta.data$QC == 'Pass', 
                                  'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 5 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', 
                                  paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 10 & seurat_object@meta.data$QC == 'Pass', 
                                  'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 10 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', 
                                  paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list




#### Write QC summary
qc_summary_list <- list()
# Collect QC summary for each sample
for (sample_name in names(seurat_objects)) {
  qc_summary <- table(seurat_objects[[sample_name]][['QC']])
  qc_summary_df <- as.data.frame(qc_summary)
  qc_summary_df$Sample <- sample_name
  qc_summary_list[[sample_name]] <- qc_summary_df
}
qc_summary_combined <- do.call(rbind, qc_summary_list)
# Write the data frame to a tab-separated text file
write.table(qc_summary_combined, file = "output/seurat/QC_summary_V2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



## subset seurat object to keep cells that pass the QC
subset_qc <- function(seurat_object) {
  seurat_object <- subset(seurat_object, subset = QC == 'Pass')
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- subset_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

# Normalize and scale data, then run cell cycle sorting
set.seed(42)
## Load gene marker of cell type
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

# Function to normalize, scale data, and perform cell cycle scoring
process_seurat_object <- function(seurat_object, mmus_s, mmus_g2m) {
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)  # zero-centres and scales it
  seurat_object <- CellCycleScoring(seurat_object, s.features = mmus_s, g2m.features = mmus_g2m)  # cell cycle sorting
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- process_seurat_object(seurat_objects[[sample_name]], mmus_s, mmus_g2m)
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

# write output summary phase
phase_summary_list <- list()
# Collect QC phase summary for each sample
for (sample_name in names(seurat_objects)) {
  phase_summary <- table(seurat_objects[[sample_name]][[]]$Phase)
  phase_summary_df <- as.data.frame(phase_summary)
  phase_summary_df$Sample <- sample_name
  phase_summary_list[[sample_name]] <- phase_summary_df
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

# Combine all summaries into one data frame
phase_summary_combined <- do.call(rbind, phase_summary_list)
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## plot cell cycle
# Calculate proportions
phase_summary_combined_tidy <- phase_summary_combined %>%
  group_by(Sample) %>%
  mutate(Prop = Freq / sum(Freq) * 100)

phase_summary_combined_tidy$Sample <- factor(phase_summary_combined_tidy$Sample, levels = c("RNA_WT", "RNA_Bap1KO")) # Reorder untreated 1st


# Plot
pdf("output/seurat/barPlot_CellCyclePhase_V2.pdf", width=5, height=6)
ggplot(phase_summary_combined_tidy, aes(x = Sample, y = Prop, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Genotype", y = "Proportion (%)", fill = "Cell Cycle Phase") +
  theme_bw() +
  scale_fill_manual(values = c("G1" = "#1f77b4", "G2M" = "#ff7f0e", "S" = "#2ca02c")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
##

###########################################################################
# saveRDS(RNA_WT, file = "output/seurat/RNA_WT_QCV2.rds") 
# saveRDS(RNA_Bap1KO, file = "output/seurat/RNA_Bap1KO_QCV2.rds") 
###########################################################################

# Cell type annotation

## Work on the WT sample 1st

################### RNA_WT
RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6650, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000)
RNA_WT <- RunPCA(RNA_WT, npcs = 50, verbose = FALSE)

## ELBOW ###########################################################################
pdf("output/seurat/elbow_RNA_WT.pdf", width=10, height=10)
ElbowPlot(RNA_WT) # 8 or 18
dev.off()
########################################################################### USELESS...


RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6650, verbose = TRUE, variable.features.n = 3000)
RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6650, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000, vst.flavor = 'v2')
RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6650, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000)
RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6650, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000)



RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6650, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = TRUE, variable.features.n = 3000)


## Best parameter V1 (with) QCV2_ QCV2_dim50kparam70res14algo3_noCellCycleRegression)
RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6637, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 3000)

RNA_WT <- RunPCA(RNA_WT, npcs = 50, verbose = FALSE)
RNA_WT <- RunUMAP(RNA_WT, reduction = "pca", dims = 1:50, verbose = FALSE)
RNA_WT <- FindNeighbors(RNA_WT, reduction = "pca", k.param = 70, dims = 1:50)
RNA_WT <- FindClusters(RNA_WT, resolution = 1.4, verbose = FALSE, algorithm = 3)

# pdf("output/seurat/UMAP_RNA_WT-dim10kparam15res04_noRegression.pdf", width=8, height=5)
# pdf("output/seurat/UMAP_RNA_WT-dim30kparam15res04_allRegression_vstv2.pdf", width=8, height=5)
# pdf("output/seurat/UMAP_RNA_WT-dim18kparam15res04_allRegression.pdf", width=8, height=5)
# pdf("output/seurat/UMAP_RNA_WT-dim18kparam15res04_noCellCycleRegression.pdf", width=8, height=5)
# pdf("output/seurat/UMAP_RNA_WT-dim30kparam15res04_noCellCycleNoRiboRegression.pdf", width=8, height=5)

pdf("output/seurat/UMAP_RNA_WT-QCV2_dim50kparam70res14algo3_noCellCycleRegression.pdf", width=8, height=5)
DimPlot(RNA_WT, reduction = "umap", label=TRUE)
dev.off()


# Check QC metrics

pdf("output/seurat/VlnPlot_QCmetrics_RNA_WT-dim50kparam30res07_noCellCycleRegression.pdf", width=20, height=5)
VlnPlot(RNA_WT,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_RNA_WT_Phase-dim50kparam30res07_noCellCycleRegression.pdf", width=8, height=5)
DimPlot(RNA_WT, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  

pdf("output/seurat/FeaturePlot_QCmetrics_RNA_WT_mt-allMarkersList4-dim50kparam30res07_noCellCycleRegression.pdf", width=15, height=10)
FeaturePlot(RNA_WT, features = c("percent.mt"), cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_RNA_WT_rb-allMarkersList4-dim50kparam30res07_noCellCycleRegression.pdf", width=15, height=10)
FeaturePlot(RNA_WT, features = c("percent.rb"), cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_RNA_WT_RNAfeat-allMarkersList4-dim50kparam30res07_noCellCycleRegression.pdf", width=15, height=10)
FeaturePlot(RNA_WT, features = c("nFeature_RNA"), cols = c("grey", "red"), min.cutoff = 4000)
dev.off()


pdf("output/seurat/FeaturePlot_QCmetrics_RNA_WT_RNAcount-allMarkersList4-dim50kparam30res07_noCellCycleRegression.pdf", width=15, height=10)
FeaturePlot(RNA_WT, features = c("nCount_RNA"), cols = c("grey", "red"), min.cutoff = 100)
dev.off()

# Check some genes

DefaultAssay(RNA_WT) <- "SCT" # For vizualization either use SCT or norm RNA




# pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList1-dim10kparam15res04_noRegression.pdf", width=15, height=20)
# pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList3-dim18kparam15res04_allRegression.pdf", width=15, height=30)
# pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList3-dim30kparam15res04_allRegression_vstv2.pdf", width=15, height=30)
# pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList3-dim18kparam15res04_noCellCycleRegression.pdf", width=15, height=30)


pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList4-QCV2_dim50kparam70res14algo3_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT, features = c("Pax6", "Eomes", "Prox1", "Neurod1", "Cck", "Crym", "Snca", "Tac2", "Pantr1", "Satb2", "Gad1", "Lhx1", "Nts"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()








pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList1-dim50kparam15res04_allRegression.pdf", width=15, height=20)

FeaturePlot(RNA_WT, features = c("Pax6","Eomes", "Prox1", "Pdzd2", "Tox3", "Dkk3", "Calb2", "Ociad2", "Fibcd1", "Pou3f1", "Wfs1", "Dcn", "Cacng5", "Ly6c1", "Coch", "Gad1", "Lhx1", "Hbb-bs", "Aif1", "Gfap", "Sox5", "Nts", "Satb2", "Cdh5"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList2-dim50kparam15res04_allRegression.pdf", width=15, height=25)

FeaturePlot(RNA_WT, features = c("Calb2", "Prox1", "Neurod1", "Mfap4", "Ppp1r14c", "Cck", "Lmo1", "Bcl11b", "Opcml", "Nrp1", "Crym", "S100a10", "Snca", "Grp", "Zbtb20", "Tac2", "Pcp4", "Mdk", "Lrpap1", "Mfap4", "Meis2", "Pantr1", "Pou3f2", "Cited2", "Pou3f3", "Mef2c", "Satb2", "Dab1", "Ntm", "Ptprz1", "Meg3", "Lmo3", "Dync1i1", "Ssbp2", "B3gat1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()





#################### LogNormalize ####################
RNA_WT <- NormalizeData(RNA_WT, normalization.method = "LogNormalize", scale.factor = 10000)
## Discover the 2000 first more variable genes
RNA_WT <- FindVariableFeatures(RNA_WT, selection.method = "vst", nfeatures = 2000)
## scale data to Z score (value centered around 0 and +/- 1)
all.genes <- rownames(RNA_WT)
RNA_WT <- ScaleData(RNA_WT, features = all.genes,vars.to.regress = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score"))

RNA_WT <- ScaleData(RNA_WT, features = all.genes)

## PCA
RNA_WT <- RunPCA(RNA_WT, features = VariableFeatures(object = RNA_WT))

# clustering
RNA_WT <- FindNeighbors(RNA_WT, dims = 1:18, k.param = 15)
RNA_WT <- FindClusters(RNA_WT, resolution = 0.4)
RNA_WT <- RunUMAP(RNA_WT, dims = 1:18, verbose = F)

pdf("output/seurat/UMAP_RNA_WT-dim18kparam15res04_allRegression_LogNormalize.pdf", width=10, height=10)
pdf("output/seurat/UMAP_RNA_WT-dim18kparam15res04_noRegression_LogNormalize.pdf", width=10, height=10)

DimPlot(RNA_WT,label.size = 4,repel = T,label = T)
dev.off()



pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList3-dim18kparam15res04_allRegression_LogNormalize.pdf", width=15, height=30)
pdf("output/seurat/FeaturePlot_SCT_RNA_WT-allMarkersList3-dim18kparam15res04_noRegression_LogNormalize.pdf", width=15, height=30)


FeaturePlot(RNA_WT, features = c("Pax6", "Eomes", "Prox1", "Neurod1", "Mfap4", "Ppp1r14c", "Pdzd2", "Tox3", "Dkk3", "Cck", "Lmo1", "Bcl11b", "Opcml", "Nrp1", "Crym", "Snca", "Zbtb20", "Tac2", "Mdk", "Lrpap1", "Meis2", "Pantr1", "Pou3f2", "Pou3f3", "Mef2c", "Satb2", "Dab1", "Ntm", "Ptprz1", "Gad1", "Lhx1", "Sox5", "Nts", "Meg3", "Lmo3", "Dync1i1", "Ssbp2", "B3gat1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()
################################################################################





# Data integration WT and Bap1KO


RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6637, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 2000) %>% 
    RunPCA(npcs = 20, verbose = FALSE)
RNA_Bap1KO <- SCTransform(RNA_Bap1KO, method = "glmGamPoi", ncells = 6938, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 2000) %>% 
    RunPCA(npcs = 20, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(RNA_WT = RNA_WT, RNA_Bap1KO = RNA_Bap1KO)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 2000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

embryo.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
RNA_WT_Bap1KO.sct <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT")

set.seed(42)

DefaultAssay(RNA_WT_Bap1KO.sct) <- "integrated"

RNA_WT_Bap1KO.sct <- RunPCA(RNA_WT_Bap1KO.sct, verbose = FALSE, npcs = 20)
RNA_WT_Bap1KO.sct <- RunUMAP(RNA_WT_Bap1KO.sct, reduction = "pca", dims = 1:20, verbose = FALSE)
RNA_WT_Bap1KO.sct <- FindNeighbors(RNA_WT_Bap1KO.sct, reduction = "pca", k.param = 70, dims = 1:20)
RNA_WT_Bap1KO.sct <- FindClusters(RNA_WT_Bap1KO.sct, resolution = 0.9, verbose = FALSE, algorithm = 4)

RNA_WT_Bap1KO.sct$orig.ident <- factor(RNA_WT_Bap1KO.sct$orig.ident, levels = c("RNA_WT", "RNA_Bap1KO")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Bap1KO_noSplit-QCV2_dim20kparam70res09algo4feat2500_noCellCycleRegression.pdf", width=8, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", label=TRUE)
dev.off()



DefaultAssay(RNA_WT_Bap1KO.sct) <- "SCT" # For vizualization either use SCT or norm RNA


pdf("output/seurat/FeaturePlot_QCmetrics_RNA_WT_Phase-QCV2_dim20kparam70res09algo4feat2500_noCellCycleRegression.pdf", width=6, height=5)
DimPlot(RNA_WT_Bap1KO.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  



pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersList4-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pax6", "Eomes", "Prox1", "Neurod1", "Cck", "Crym", "Snca", "Tac2", "Pantr1", "Satb2", "Gad1", "Lhx1", "Nts"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersListTopCA1-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Lmo1","Bcl11b","Opcml","Nrp1","Slc16a2","Lmo7","Nrp2","Xpr1","Crym","Zeb2","Insm1","Itm2b"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersListTopCA3-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("S100a10","Snca","Grp","Zbtb20","Stmn2","Mapt","Nrp2","Vim","Selm","Gap43","Caly","Slc16a2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersListTopDGGC-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Mfap4","Ppp1r14c","Kcnk1","Zbtb20","Btbd17","Gadd45g","Slc1a2","Cpe","Tmem178","Pcp4","Sema5a","Cd63"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersListTopDL-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pcp4","Mdk","Lrpap1","Mfap4","Nrn1","Hs3st1","Rnd2","Ulk4","Unc5d","Sstr2","Eomes","Cited2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersListTopML-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pou3f2","Cited2","Pou3f3","Gria2","Gm17750","Tsc22d1","Cttnbp2","Frmd4b","Fabp7","Igfbpl1","Nkain3","Sox11"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersListTopUL-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Dab1","Ntm","Ptprz1","Pou3f1","Tenm3","Itpr1","Flrt3","Tenm2","Sobp","Nefm","Limch1","Id2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersListTopSubC-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Sox5","Meg3","Lmo3","Dync1i1","Ssbp2","B3gat1","Nrgn","Tle4","Fezf2","Nr4a2","Bcl11b","Rprm"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-microgliaMarkers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Aif1", "Itgam","Cx3cr1", "Gpr34", "Gpr183", "P2ry12", "P2ry13", "Csf1r", "Tmem119"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-astrocyteMarkers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Gfap", "Slc1a2", "Acsl6", "Agt", "Aqp4", "Apoe", "S100b", "Sox9", "Gsta4"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-glutaminergicMarkers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Slc17a7", "Slc17a6", "Slc17a8", "Slc1a1", "Gls", "Meis2", "Slc1a6", "Grin1", "Grin2b"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-OPCMarkers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pdgfra", "Aldoc", "Olig1", "Epn2", "Neu4", "Nkx6-2", "Fyn", "Tnr", "Pcdh15"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-oligodendrocyteMarkers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Mog", "Mbp", "Mag", "Cldn14", "Klk6", "Eml1", "Nipal4", "Plp1", "Ermn"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-mossyCellsMarkers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Calb1", "Reln", "Gad1", "Slc17a6", "Npy", "Gria2", "Gria3"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-MeningealCellsMarkers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Cldn11", "Fn1", "Fbn1", "Col1a1", "Col3a1", "Pdgfrb", "Igf2", "Vtn", "Gjb6", "Aldh1a2" ), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-MeningealCellsPangloa-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Igfbp2", "Slc47a1", "Nov", "Nnat", "Ptgds", "Il33", "Pdgfra", "Lum", "Dc", "Foxc1", "Vtn", "Igf2", "Alx4", "Aldh1a2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-MeningealCellsPangloaTop-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=10, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Lum", "Foxc1", "Vtn", "Igf2", "Aldh1a2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-NSCQuiescActivChatGPT-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=10, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Id1", "Hes1", "Mki67", "Pcna", "Vim" ), max.cutoff = 2, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-INChatGPT-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=10, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Gad1", "Grin2d", "Sst", "Tac1", "Calb1", "Npy", "Gria3", "Lhx6", "Vip" ), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-interneuronMarkers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pvalb", "Cck", "Npy", "Nos1", "Sst", "Vip", "Fgf12", "Nxph11", "Kcnc2", "Calb2", "Calb1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-interneuron1Markers-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Vsx2", "En1", "Etv1", "Evx1", "Gad1", "Isl1", "Lhx1", "Lhx5", "Lhx3", "Lhx6", "Grm1", "Grin2d"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-top1Markers-QCV2_dim20kparam70res09algo4feat2500_noCellCycleRegression
.pdf", width=13, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pax6", "Eomes", "Sema5a", "Hs3st1", "Igfbpl1", "Satb2", "Nts", "Cck", "Snca", "Gad1", "Lhx1", "Pdgfra", "Csf1r"), cols = c("grey", "red"), max.cutoff = 1)
dev.off()


pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-Pax6Tbr2Prox1-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pax6", "Eomes", "Prox1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-Bap1Pax6Eomes-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=10, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Bap1", "Pax6", "Eomes"), max.cutoff = 1, cols = c("grey", "red"), split.by = "orig.ident")
dev.off()

pdf("output/seurat/RidgePlot_SCT_RNA_WT_Bap1KO-Bap1Pax6Eomes-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression.pdf", width=5, height=5)
RidgePlot(RNA_WT_Bap1KO.sct, features = c("Bap1", "Pax6", "Eomes"), cols = c("blue","red"), group.by = 'orig.ident', ncol =1, assay = "SCT")
dev.off()

pdf("output/seurat/VlnPlot_SCT_RNA_WT_Bap1KO-Bap1Pax6Eomes-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression.pdf", width=3, height=10)
VlnPlot(RNA_WT_Bap1KO.sct, features = c("Bap1", "Pax6", "Eomes"), cols = c("blue","red"), group.by = 'orig.ident', ncol =1, assay = "RNA")
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-interneuronGenes-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression.pdf", width=10, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Gad1", "Grin2d", "Calb1", "Npy", "Gria3", "Lhx6"),cols = c("grey", "red"), max.cutoff = 1 )
dev.off()


pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-Bap1complex-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=10, height=35)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Ogt", "Hcfc1", "Foxk1", "Foxk2", "Yy1", "Bap1", "Mbd5", "Mbd6", "Asxl1", "Asxl2", "Asxl3"),cols = c("grey", "red"), split.by = "orig.ident", max.cutoff = 2 )
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-Bap1-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=10, height=5)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Bap1"),cols = c("grey", "red"), split.by = "orig.ident", max.cutoff = 1 )
dev.off()

pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-DentateGyrusDev-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=10, height=15)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pax6", "Eomes", "Prox1"),cols = c("grey", "red"), split.by = "orig.ident")
dev.off()


pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-Phf7-QCV2_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=7, height=5)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Phf7"),cols = c("grey", "red"), split.by = "orig.ident")
dev.off()

pdf("output/seurat/UMAP_WT_Bap1KO_numeric_V1.pdf", width=12, height=6)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", split.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 6)
dev.off()


pdf("output/seurat/UMAP_WT_Bap1KO_noSplit_numeric_V1.pdf", width=7, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 6)
dev.off()

#overlapping orig.ident
pdf("output/seurat/UMAP_WT_Bap1KO_numeric_overlap_V1.pdf", width=6, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", group.by = "orig.ident", pt.size = 0.000001, cols = c("blue","red"))
dev.off()



#### Check raw/norm Bap1 count #####################
RNA_WT_Bap1KO.sct_WT <- subset(RNA_WT_Bap1KO.sct, subset = orig.ident == 'RNA_WT')
RNA_WT_Bap1KO.sct_Bap1KO <- subset(RNA_WT_Bap1KO.sct, subset = orig.ident == 'RNA_Bap1KO')


AverageExpression(RNA_WT_Bap1KO.sct, features = "Bap1", group.by = "orig.ident")


AverageExpression(RNA_WT_Bap1KO.sct, features = "Bap1", group.by = "seurat_clusters")


AverageExpression(RNA_WT_Bap1KO.sct_WT , features = "Bap1", group.by = "seurat_clusters")
AverageExpression(RNA_WT_Bap1KO.sct_Bap1KO , features = "Bap1", group.by = "seurat_clusters")

####################################


############ V1 naming

Cluster1 = PyNs_SubC_CA1 (subiculum PyNs)
Cluster2 = PyNs_SubC_CA23_1 (subiculum PyNs)
Cluster3 = PyNs_RSC_ML (Retrosplenial Cortical Pyramidal neurons, middle layer)
Cluster4 = IN_1 (interneuron)
Cluster5 = PyNs_RSC_UL (Retrosplenial Cortical Pyramidal neurons, upper layer)
Cluster6 = PyNs_SubC_CA23_2 (subiculum PyNs)
Cluster7 = DG_GC (Dentate Gyrus granule cells)
Cluster8 = NSC_1 (Neural Stem Cells)
Cluster9 = IN_2 (interneuron)
Cluster10 = SubC_1 (subiculum)
Cluster11 = NSC_2 (Neural Stem Cells)
Cluster12 = IP (Intermediate Progenitors)
Cluster13 = NSC_3 (Neural Stem Cells)
Cluster14 = Unknown_1
Cluster15 = SubC_2 (subiculum)
Cluster16 = CR (Cajal Retzius)
Cluster17 = PyNs_RSC_DL (Retrosplenial Cortical Pyramidal neurons, deep layer)
Cluster18 = Unknown_2
Cluster19 = Unknown_3


new.cluster.ids <- c(
  "PyNs_SubC_CA1" ,
  "PyNs_SubC_CA23_1" ,
  "PyNs_RSC_ML" ,
  "IN_1" ,
  "PyNs_RSC_UL" ,
  "PyNs_SubC_CA23_2" ,
  "DG_GC" ,
  "NSC_1" ,
  "IN_2" ,
  "SubC_1" ,
  "NSC_2" ,
  "IP" ,
  "NSC_3" ,
  "Unknown_1" ,
  "SubC_2" ,
  "CR" ,
  "PyNs_RSC_DL" ,
  "Unknown_2",
  "Unknown_3" 
)

names(new.cluster.ids) <- levels(RNA_WT_Bap1KO.sct)
RNA_WT_Bap1KO.sct <- RenameIdents(RNA_WT_Bap1KO.sct, new.cluster.ids)

RNA_WT_Bap1KO.sct$cluster.annot <- Idents(RNA_WT_Bap1KO.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Bap1KO_label_V1.pdf", width=12, height=6)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", split.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/test.pdf", width=12, height=6)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", split.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_WT_Bap1KO_noSplit_label_V1.pdf", width=7, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 4)
dev.off()

#overlapping orig.ident
pdf("output/seurat/UMAP_WT_Bap1KO_label_overlap_V1.pdf", width=6, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", group.by = "orig.ident", pt.size = 0.000001, cols = c("blue","red"))
dev.off()


# All in dotplot
DefaultAssay(RNA_WT_Bap1KO.sct) <- "SCT"

CA1 = Cck
CA3 = Crym, Snca
Pyramidal neurons middle layer (ML) = Pantr1
Interneurons (IN) = Gad1
Pyramidal neurons upper layer (UL) = Satb2
Dentate Gyrus Granule Cells (DG) = Prox1, Neurod1
Neural Stem Cells (NSC) = Pax6 (should have more cell types)
Subiculum (SubC) = Nts
Intermediate Progenitors (IP) = Eomes
Cajal Retzius (CR) = Lhx1
Pyramidal neurons deep layer (DL) = Tac2

all_markers <- c(
  "Pax6" ,
  "Eomes",
  "Prox1", "Neurod1",
  "Tac2",
  "Pantr1",
  "Satb2",
  "Nts",
  "Cck",
  "Crym", "Snca",
  "Gad1",
  "Lhx1"
)



levels(RNA_WT_Bap1KO.sct) <- c(
  "NSC_1" ,
  "NSC_2" ,
  "NSC_3" ,
  "IP" ,
  "DG_GC" ,
  "PyNs_RSC_DL" ,
  "PyNs_RSC_ML" ,
  "PyNs_RSC_UL" ,
  "SubC_1" ,
  "SubC_2" ,
  "PyNs_SubC_CA1" ,
  "PyNs_SubC_CA23_1" ,
  "PyNs_SubC_CA23_2" ,
  "IN_1" ,
  "IN_2" ,
  "CR" ,
  "Unknown_1" ,
  "Unknown_2",
  "Unknown_3" 
)



pdf("output/seurat/DotPlot_SCT_WT_Bap1KO_label_V1.pdf", width=7, height=4.5)
DotPlot(RNA_WT_Bap1KO.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()


########################################################



## Downsampling with bootstrap to compare the nb of cell per cell types

library("tidyverse")

### Identify the unique clusters
unique_clusters <- unique(Idents(RNA_WT_Bap1KO.sct))

### Create empty matrices to store cell counts
control_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
Bap1KO_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(control_clusters_counts) <- unique_clusters
colnames(Bap1KO_clusters_counts) <- unique_clusters

### Loop through 100 iterations
RNA_WT_Bap1KO.sct_WT <- which(RNA_WT_Bap1KO.sct$orig.ident == 'RNA_WT')
RNA_WT_Bap1KO.sct_Bap1KO <- which(RNA_WT_Bap1KO.sct$orig.ident == 'RNA_Bap1KO')

for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  RNA_WT_Bap1KO.sct_Bap1KO_downsample <- sample(RNA_WT_Bap1KO.sct_Bap1KO, 6637)
  RNA_WT_Bap1KO.sct_integrated_downsample <- RNA_WT_Bap1KO.sct[,c(RNA_WT_Bap1KO.sct_Bap1KO_downsample, RNA_WT_Bap1KO.sct_WT)]

  # Count nb of cells in each cluster
  control_clusters <- table(Idents(RNA_WT_Bap1KO.sct_integrated_downsample)[RNA_WT_Bap1KO.sct_integrated_downsample$orig.ident == "RNA_WT"])
  Bap1KO_clusters <- table(Idents(RNA_WT_Bap1KO.sct_integrated_downsample)[RNA_WT_Bap1KO.sct_integrated_downsample$orig.ident == "RNA_Bap1KO"])

  # Align the counts with the unique clusters
  control_clusters_counts[i, names(control_clusters)] <- as.numeric(control_clusters)
  Bap1KO_clusters_counts[i, names(Bap1KO_clusters)] <- as.numeric(Bap1KO_clusters)
}


### Calculate mean and standard error
mean_control_clusters <- colMeans(control_clusters_counts)
mean_Bap1KO_clusters <- colMeans(Bap1KO_clusters_counts)
std_error_WT_clusters <- apply(control_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("WT", "Bap1KO")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- control_clusters_counts[j,i]
    contingency_table[1,2] <- Bap1KO_clusters_counts[j,i]
    contingency_table[2,1] <- sum(control_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(Bap1KO_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_control_clusters),
  untreated = mean_control_clusters,
  dasatinib = mean_Bap1KO_clusters,
  std_error_WT = std_error_WT_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_WT, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "WT", "Bap1KO"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.05, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("WT", "Bap1KO")) # Reorder untreated 1st
plot_data$cluster <- factor(plot_data$cluster, levels = c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")) 



# Plotting using ggplot2
pdf("output/seurat/Cluster_cell_counts_BootstrapDownsampling100_WT_Bap1KO_numeric.pdf", width=9, height=4)
ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = filter(plot_data, condition == "Bap1KO"),
    aes(label = significance, y = value + std_error_WT_clusters),
    vjust = -0.8,
    position = position_dodge(0.9), size = 5
  ) +
  scale_fill_manual(values = c("WT" = "#4365AE", "Bap1KO" = "#981E33")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13)) +
  theme(axis.text.y = element_text(size = 13)) +
  ylim(0,900)
dev.off()




# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs

## DEGs keeping ALL genes
RNA_WT_Bap1KO.sct$celltype.stim <- paste(RNA_WT_Bap1KO.sct$seurat_clusters, RNA_WT_Bap1KO.sct$orig.ident,
    sep = "-")
Idents(RNA_WT_Bap1KO.sct) <- "celltype.stim"

cluster1 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "1-RNA_Bap1KO", ident.2 = "1-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
cluster2 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "2-RNA_Bap1KO", ident.2 = "2-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
cluster3 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "3-RNA_Bap1KO", ident.2 = "3-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
cluster4 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "4-RNA_Bap1KO", ident.2 = "4-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
cluster5 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "5-RNA_Bap1KO", ident.2 = "5-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
cluster6 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "6-RNA_Bap1KO", ident.2 = "6-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster7 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "7-RNA_Bap1KO", ident.2 = "7-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster8 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "8-RNA_Bap1KO", ident.2 = "8-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster9 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "9-RNA_Bap1KO", ident.2 = "9-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster10 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "10-RNA_Bap1KO", ident.2 = "10-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster11 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "11-RNA_Bap1KO", ident.2 = "11-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster12 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "12-RNA_Bap1KO", ident.2 = "12-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster13 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "13-RNA_Bap1KO", ident.2 = "13-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster14 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "14-RNA_Bap1KO", ident.2 = "14-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster15 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "15-RNA_Bap1KO", ident.2 = "15-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster16 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "16-RNA_Bap1KO", ident.2 = "16-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster17 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "17-RNA_Bap1KO", ident.2 = "17-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster18 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "18-RNA_Bap1KO", ident.2 = "18-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster19 <- FindMarkers(RNA_WT_Bap1KO.sct, ident.1 = "19-RNA_Bap1KO", ident.2 = "19-RNA_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  


### save output
## write.table(cluster1, file = "output/seurat/cluster1-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster1, file = "output/seurat/cluster1-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster2, file = "output/seurat/cluster2-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster3, file = "output/seurat/cluster3-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster4, file = "output/seurat/cluster4-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster5, file = "output/seurat/cluster5-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster6, file = "output/seurat/cluster6-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster7, file = "output/seurat/cluster7-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster8, file = "output/seurat/cluster8-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster9, file = "output/seurat/cluster9-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster10, file = "output/seurat/cluster10-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster11, file = "output/seurat/cluster11-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster12, file = "output/seurat/cluster12-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster13, file = "output/seurat/cluster13-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster14, file = "output/seurat/cluster14-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster15, file = "output/seurat/cluster15-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster16, file = "output/seurat/cluster16-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster17, file = "output/seurat/cluster17-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster18, file = "output/seurat/cluster18-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster19, file = "output/seurat/cluster19-Bap1KO_response_V1_allGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)


#### import all clsuter DEGs output :
cluster_types <- c("cluster1", "cluster2", "cluster3", 
                   "cluster4", "cluster5", "cluster6", 
                   "cluster7", "cluster8", 
                   "cluster9", "cluster10", "cluster11", 
                   "cluster12", "cluster13", "cluster14", "cluster15", 
                   "cluster16", "cluster17", "cluster18", "cluster19")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Bap1KO_response_V1_allGenes.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}







# DEGs nb dotplot
## Initialize an empty data frame to store the summary
DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())

## List of cell types
cell_types <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19")

## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/seurat/", cell_type, "-Bap1KO_response_V1_allGenes.txt", sep = "")
  deg_data <- read.table(file_name, header = TRUE, sep = "\t") ## Read the DEGs data
  num_degs <- sum(deg_data$p_val_adj < 0.05) ## Count the number of significant DEGs
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))  ## Append to the summary table
}


DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19")) 

DEG_count$Cluster_Number <- as.numeric(sub("cluster", "", DEG_count$Cell_Type))

DEG_count$Cluster_Number <- factor(DEG_count$Cluster_Number, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")) 


## Dotplot

cell_type_colors <- c(
  "cluster1" = "#FFB6C1", "cluster2" = "#FFA07A", "cluster3" = "#FFD700", "cluster4" = "#ADFF2F", 
  "cluster5" = "#7FFF00", "cluster6" = "#32CD32", "cluster7" = "#3CB371", "cluster8" = "#00FA9A", 
  "cluster9" = "#00CED1", "cluster10" = "#4682B4", "cluster11" = "#1E90FF", "cluster12" = "#6495ED", 
  "cluster13" = "#4169E1", "cluster14" = "#BA55D3", "cluster15" = "#DA70D6", "cluster16" = "#EE82EE", 
  "cluster17" = "#FF69B4", "cluster18" = "#FF1493", "cluster19" = "#DB7093"
)

# Generate the dot plot
pdf("output/seurat/Dotplot_Bap1KO_DEG_count_V1.pdf", width=9, height=4)
ggplot(DEG_count, aes(x = Cluster_Number, y = 1) )+
  geom_point(aes(size = ifelse(Num_DEGs == 0, 1, Num_DEGs), fill = Cell_Type) , shape = 21, color = "black") +
  scale_size_continuous(range = c(1, 15)) +
  scale_fill_manual(values = cell_type_colors, guide = "none") +  # Remove the legend for cell types
  theme_void() +
  labs(title = "Number of DEGs per Cell Type", x = "Cell Type", y = "", size = "Number of DEGs") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 50, size = 15),  # Adjust the position of x-axis labels
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(),
    legend.position = "right"
  ) +
  guides(size = guide_legend(title = "Number of DEGs", title.position = "top", title.hjust = 0.5))  # Keep only the size legend
dev.off()




### Find all markers 
all_markers <- FindAllMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(all_markers, file = "output/seurat/srat_WT_Bap1KO_all_markers_V1.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Display the top 10 CONSERVED marker genes of each cluster
Idents(RNA_WT_Bap1KO.sct) <- "seurat_clusters"

## DEGs cluster versus all other
cluster1.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "1", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster1")
cluster2.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "2", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster2")
cluster3.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "3", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster3")
cluster4.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "4", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster4")
cluster5.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "5", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster5")
cluster6.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "6", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster6")
cluster7.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "7", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster7")
cluster8.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "8", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster8")
cluster9.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "9", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster9")
cluster10.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "10", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster10")
cluster11.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "11", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster11")
cluster12.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "12", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster12")
cluster13.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "13", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster13")
cluster14.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "14", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster14")
cluster15.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "15", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster15")
cluster16.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "16", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster16")
cluster17.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "17", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster17")
cluster18.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "18", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster18")
cluster19.conserved <- FindConservedMarkers(RNA_WT_Bap1KO.sct, assay = "RNA", ident.1 = "19", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster19")


## Combine all conserved markers into one data frame
all_conserved <- bind_rows(cluster1.conserved,cluster2.conserved,cluster3.conserved,cluster4.conserved,cluster5.conserved,cluster6.conserved,cluster7.conserved,cluster8.conserved,cluster9.conserved,cluster10.conserved,cluster11.conserved,cluster12.conserved,cluster13.conserved,cluster14.conserved,cluster15.conserved,cluster16.conserved,cluster17.conserved,cluster18.conserved,cluster19.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/seurat/srat_all_conserved_markers_WT_Bap1KO_V1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## Find the top 5 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 5) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19")))


## Write the top 10 conserved markers for each cluster to a file

## Visualize the top 10/3 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(RNA_WT_Bap1KO.sct) <- c("1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7",
  "8",
  "9",
  "10",
  "11",
  "12",
  "13",
  "14",
  "15",
  "16",
  "17",
  "18",
  "19")

pdf("output/seurat/DotPlot_SCT_top5_conserved_WT_Bap1KO_V1.pdf", width=19, height=5)
DotPlot(RNA_WT_Bap1KO.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# save
## saveRDS(RNA_WT_Bap1KO.sct, file = "output/seurat/RNA_WT_Bap1KO.sct_V1_numeric.rds") 



RNA_WT_Bap1KO.sct <- readRDS(file = "output/seurat/RNA_WT_Bap1KO.sct_V1_numeric.rds")





############################ EasyCellType automatic annotation ##########################################
# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_WT_Bap1KO_all_markers_V1.txt", header = TRUE, row.names = 1)
### Filter either WT or cYAPKO
all_markers <- all_markers[grepl("RNA_WT$", all_markers$cluster), ]
all_markers <- all_markers[grepl("RNA_Bap1KO$", all_markers$cluster), ]

## Convert geneSymbol to EntrezID
all_markers$entrezid <- mapIds(org.Mm.eg.db,
                           keys=all_markers$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
all_markers <- na.omit(all_markers)

## Sort the datafram (data frame containing Entrez IDs, clusters and expression scores)

all_markers_sort <- data.frame(gene=all_markers$entrezid, cluster=all_markers$cluster, 
                      score=all_markers$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(all_markers_sort[, 1:3])

## Run the enrihcment analysis
annot.GSEA <- easyct(input.d, db="panglao", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=NULL, p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain", "Cerebellum", "Hippocampus"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?



annot.GSEA <- easyct(input.d, db="clustermole", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?

## plots

pdf("output/seurat/EasyCellType_dotplot_SCT_WT-cellmarker_brainCerebellumHippocampus.pdf", width=6, height=8)
pdf("output/seurat/EasyCellType_dotplot_SCT_WT-clustermole_brain.pdf", width=6, height=8)

plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()





##########################################################################################





############ V2 naming

Cluster1 = PyNs_SubC_CA1 (subiculum PyNs)
Cluster2 = PyNs_SubC_CA23_1 (subiculum PyNs)
Cluster3 = PyNs_RSC_ML (Retrosplenial Cortical Pyramidal neurons, middle layer)
Cluster4 = IN_1 (interneuron)
Cluster5 = PyNs_RSC_UL (Retrosplenial Cortical Pyramidal neurons, upper layer)
Cluster6 = PyNs_SubC_CA23_2 (subiculum PyNs)
Cluster7 = DG_GC (Dentate Gyrus granule cells)
Cluster8 = NSC_1 (Neural Stem Cells)
Cluster9 = IN_2 (interneuron)
Cluster10 = SubC_1 (subiculum)
Cluster11 = NSC_2 (Neural Stem Cells)
Cluster12 = IP (Intermediate Progenitors)
Cluster13 = NSC_3 (Neural Stem Cells)
Cluster14 = OPC (Oligodendrocyte progenitor cells)
Cluster15 = SubC_2 (subiculum)
Cluster16 = CR (Cajal Retzius)
Cluster17 = PyNs_RSC_DL (Retrosplenial Cortical Pyramidal neurons, deep layer)
Cluster18 = Microglia
Cluster19 = Unknown


new.cluster.ids <- c(
  "PyNs_SubC_CA1" ,
  "PyNs_SubC_CA23_1" ,
  "PyNs_RSC_ML" ,
  "IN_1" ,
  "PyNs_RSC_UL" ,
  "PyNs_SubC_CA23_2" ,
  "DG_GC" ,
  "NSC_1" ,
  "IN_2" ,
  "SubC_1" ,
  "NSC_2" ,
  "IP" ,
  "NSC_3" ,
  "OPC" ,
  "SubC_2" ,
  "CR" ,
  "PyNs_RSC_DL" ,
  "Microglia",
  "Unknown" 
)

names(new.cluster.ids) <- levels(RNA_WT_Bap1KO.sct)
RNA_WT_Bap1KO.sct <- RenameIdents(RNA_WT_Bap1KO.sct, new.cluster.ids)

RNA_WT_Bap1KO.sct$cluster.annot <- Idents(RNA_WT_Bap1KO.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Bap1KO_label_V2.pdf", width=12, height=6)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", split.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()


pdf("output/seurat/UMAP_WT_Bap1KO_noSplit_label_V2.pdf", width=7, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 4)
dev.off()

#overlapping orig.ident
pdf("output/seurat/UMAP_WT_Bap1KO_label_overlap_V1.pdf", width=6, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", group.by = "orig.ident", pt.size = 0.000001, cols = c("blue","red"))
dev.off()


# All in dotplot
DefaultAssay(RNA_WT_Bap1KO.sct) <- "SCT"

Neural Stem Cells (NSC) = Pax6 (should have more cell types)
Intermediate Progenitors (IP) = Eomes
Dentate Gyrus Granucle Cells (DG) = Prox1, Neurod1, Sema5a
CA1 = Cck, Insm1
CA3 = Crym, Snca, Nrp2
Pyramidal neurons deep layer (DL) = Tac2, Hs3st1, Nrn1
Pyramidal neurons middle layer (ML) = Pantr1, Igfbpl1, Frmd4b
Pyramidal neurons upper layer (UL) = Satb2, Itpr1
Interneurons (IN) = Gad1, Grin2d, Reln, Calb1, Npy, Gria3, Lhx6
Cajal Retzius (CR) = Lhx1
Subiculum (SubC) = Nts, Nr4a2, Lmo3, B3gat1
Microglia = Csf1r, Gpr34, Gpr183, Cx3cr1
OPC = Pdgfra, Olig1


all_markers <- c(
  "Pax6" ,
  "Eomes",
  "Prox1", "Neurod1", "Sema5a",
  "Tac2", "Hs3st1", "Nrn1",
  "Pantr1", "Igfbpl1", "Frmd4b",
  "Satb2", "Itpr1",
  "Nts", "Nr4a2", "Lmo3", "B3gat1",
  "Cck", "Insm1",
  "Crym", "Snca", "Nrp2",
  "Gad1", "Grin2d", "Calb1", "Npy", "Gria3", "Lhx6", # Reln removed
  "Lhx1",
  "Pdgfra", "Olig1",
  "Csf1r", "Gpr34", "Gpr183", "Cx3cr1"
)



levels(RNA_WT_Bap1KO.sct) <- c(
  "NSC_1" ,
  "NSC_2" ,
  "NSC_3" ,
  "IP" ,
  "DG_GC" ,
  "PyNs_RSC_DL" ,
  "PyNs_RSC_ML" ,
  "PyNs_RSC_UL" ,
  "SubC_1" ,
  "SubC_2" ,
  "PyNs_SubC_CA1" ,
  "PyNs_SubC_CA23_1" ,
  "PyNs_SubC_CA23_2" ,
  "IN_1" ,
  "IN_2" ,
  "CR" ,
  "OPC" ,
  "Microglia",
  "Unknown" 
)



pdf("output/seurat/DotPlot_SCT_WT_Bap1KO_label_V2.pdf", width=11, height=4.5)
DotPlot(RNA_WT_Bap1KO.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_WT_Bap1KO_label_V2vertical.pdf", width=11, height=4.5)
DotPlot(RNA_WT_Bap1KO.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()


########################################################


# Cell type proportion


pt <- table(Idents(RNA_WT_Bap1KO.sct), RNA_WT_Bap1KO.sct$orig.ident)
pt <- as.data.frame(pt)
pt <- pt %>%
  group_by(Var2) %>%
  mutate(Proportion = Freq / sum(Freq))

pt$Var1 <- as.character(pt$Var1)

pdf("output/seurat/cellTypeProp_SCT_WT_Bap1KO_V2.pdf", width=5, height=5)
ggplot(pt, aes(x = Var2, y = Proportion, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)), 
            position = position_fill(vjust = 0.5), size = 3) +
  theme_bw()
dev.off()


# Cell cycle proportion per cluster
## Using numeric cluster annotation

plot_cell_cycle_per_cluster <- function(RNA_WT_Bap1KO.sct, output_dir) {
  clusters <- unique(RNA_WT_Bap1KO.sct$seurat_clusters)
  for (cluster in clusters) {
    data <- RNA_WT_Bap1KO.sct@meta.data %>%
      filter(seurat_clusters == cluster) %>%
      group_by(orig.ident, Phase) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      group_by(orig.ident) %>%
      mutate(proportion = count / sum(count)) %>%
      ungroup()
    plot <- ggplot(data, aes(x = orig.ident, y = proportion, fill = Phase)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent) +
      labs(title = paste("Cluster", cluster), x = "Genotype", y = "Proportion (%)") +
      theme_bw() +
      scale_fill_manual(values = c("G1" = "#1f77b4", "G2M" = "#ff7f0e", "S" = "#2ca02c"))
    # Save plot to PDF
    pdf(paste0(output_dir, "cellCycle_Cluster_", cluster, ".pdf"), width = 5, height = 6)
    print(plot)
    dev.off()
  }
}
plot_cell_cycle_per_cluster(RNA_WT_Bap1KO.sct, output_dir = "output/seurat/")



```



# Shinny app

## RNA only

Created with [ShinyCell](https://github.com/SGDDNB/ShinyCell); and follow [shinyapps](https://www.shinyapps.io/) to put it online 

```bash
conda activate scRNAseqV2
# conda install -c anaconda hdf5
```

```R
# installation
## devtools::install_github("SGDDNB/ShinyCell")
## install.packages('rsconnect')


# Packages
library("Seurat")
library("ShinyCell")
library("rsconnect")

# Data import EMBRYO
RNA_WT_Bap1KO.sct <- readRDS(file = "output/seurat/embryo.combined.sct.rds")
DefaultAssay(RNA_WT_Bap1KO.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(RNA_WT_Bap1KO.sct)

makeShinyApp(RNA_WT_Bap1KO.sct, scConf, gene.mapping = TRUE,
             shiny.title = "RNA_QCV2_clusterV2",
             shiny.dir = "shinyApp_RNA_QCV2_clusterV2/") 

rsconnect::deployApp('shinyApp_RNA_QCV2_clusterV2')





```



--> https://roulethomas.shinyapps.io/shinyapp_rna_qcv2_clusterv2/ 





## RNA and ATAC 

Created with [ShinyCell](https://github.com/SGDDNB/ShinyCell); and follow [shinyapps](https://www.shinyapps.io/) to put it online 

```bash
conda activate SignacV5
module load hdf5 # to read `.h5` files
```

```R
# installation
## install.packages("devtools")
## devtools::install_github("SGDDNB/ShinyCell")
## install.packages('rsconnect')
## install.packages("DT")
## install.packages("ggdendro")
## install.packages("shinyhelper")

# Packages
library("Seurat")
library("ShinyCell")
library("rsconnect")

# Data import 1st version integration ATAC RNA
multiome_WT_Bap1KO_QCV3.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV3.sct_numeric_label.rds")
DefaultAssay(multiome_WT_Bap1KO_QCV3.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(multiome_WT_Bap1KO_QCV3.sct)

makeShinyApp(multiome_WT_Bap1KO_QCV3.sct, scConf, gene.mapping = TRUE,
             shiny.title = "multiome_WT_Bap1KO_QCV3",
             shiny.dir = "shinyApp_multiome_WT_Bap1KO_QCV3/") 

rsconnect::deployApp('shinyApp_multiome_WT_Bap1KO_QCV3')


# DAta import 2nd version ATAC RNA
## saveRDS(multiome_WT_Bap1KO_QCV2vC1.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric_label.rds") ##
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric_label.rds")
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA" # 
# Generate Shiny app V2
scConf = createConfig(multiome_WT_Bap1KO_QCV2vC1.sct)
makeShinyApp(multiome_WT_Bap1KO_QCV2vC1.sct, scConf, gene.mapping = TRUE,
             shiny.title = "multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000",
             shiny.dir = "shinyApp_multiome_WT_Bap1KO_QCV2vC1/") 
rsconnect::deployApp('shinyApp_multiome_WT_Bap1KO_QCV2vC1')



# Correct1, including RNA and ATAC assay
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds")
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA" # 
# Generate Shiny app V2
scConf = createConfig(multiome_WT_Bap1KO_QCV2vC1.sct)
makeShinyApp(multiome_WT_Bap1KO_QCV2vC1.sct, scConf, gene.mapping = TRUE,
             shiny.title = "multiome_WT_Bap1KO_QCV2vC1_correct1GeneActivityLinkPeaks",
             shiny.dir = "shinyApp_multiome_WT_Bap1KO_correct1/") 
rsconnect::deployApp('shinyApp_multiome_WT_Bap1KO_correct1')



```



# Generate bigwig

Generate coverage bigwig files to check for Bap1 exon7 deletion in Bap1KO.

--> Let's just generate the bigwig from the `possorted_genome_bam.bam` generated by 10X cellranger count 



```bash
conda activate deeptools

# raw
sbatch scripts/bamtobigwig_RNA_WT.sh # 21960401 ok
sbatch scripts/bamtobigwig_RNA_Bap1KO.sh # 21960405 ok


# BPM norm (=TPM)
sbatch scripts/bamtobigwig_BPMnorm_RNA_WT.sh # 21963563 ok
sbatch scripts/bamtobigwig_BPMnorm_RNA_Bap1KO.sh # 21963630 ok


```

--> Looking at the bam, Bap1 level is same WT and Bap1KO

--> Looking at the raw bigwig, Bap1 is NOT decrease...




# Seurat analysis on all cells

To decipher the Bap1 expression issue (eg. no clear decrease of Bap1 in the KO)

- Check Bap1 expression in all cells (including the one filtered out) = doublet.
- check expression in non RNA contamination corrected files



```bash
conda activate scRNAseqV2
```

 
```R
# install.packages('SoupX')
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

## Load the matrix and Create SEURAT object
samples <- list(
  RNA_WT = "output/soupX/RNA_WT.RData",
  RNA_Bap1KO = "output/soupX/RNA_Bap1KO.RData"
)

seurat_objects <- list()

for (sample_name in names(samples)) {
  load(samples[[sample_name]])
  seurat_objects[[sample_name]] <- CreateSeuratObject(counts = out, project = sample_name)
}

## Function to assign Seurat objects to variables (unlist the list)
assign_seurat_objects <- function(seurat_objects_list) {
  for (sample_name in names(seurat_objects_list)) {
    assign(sample_name, seurat_objects_list[[sample_name]], envir = .GlobalEnv)
  }
}
assign_seurat_objects(seurat_objects) # This apply the function


# QUALITY CONTROL
# Function to add mitochondrial and ribosomal content
add_quality_control <- function(seurat_object) {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^Rp[sl]")
  return(seurat_object)
}
seurat_objects <- lapply(seurat_objects, add_quality_control) # This apply the function
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list


# Function to add doublet information
add_doublet_information <- function(sample_name, seurat_object) {
  doublet_file <- paste0("output/doublets/", sample_name, ".tsv")
  doublets <- read.table(doublet_file, header = FALSE, row.names = 1)
  colnames(doublets) <- c("Doublet_score", "Is_doublet")
  seurat_object <- AddMetaData(seurat_object, doublets)
  seurat_object$Doublet_score <- as.numeric(seurat_object$Doublet_score)
  return(seurat_object)
}
## Apply the function to each Seurat object in the list
for (sample_name in names(seurat_objects)) {
    seurat_objects[[sample_name]] <- add_doublet_information(sample_name, seurat_objects[[sample_name]])
  }

assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list






# QC filtering _ V2
## NOT APPLY !!!



# Normalize and scale data, then run cell cycle sorting
set.seed(42)
## Load gene marker of cell type
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

# Function to normalize, scale data, and perform cell cycle scoring
process_seurat_object <- function(seurat_object, mmus_s, mmus_g2m) {
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)  # zero-centres and scales it
  seurat_object <- CellCycleScoring(seurat_object, s.features = mmus_s, g2m.features = mmus_g2m)  # cell cycle sorting
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- process_seurat_object(seurat_objects[[sample_name]], mmus_s, mmus_g2m)
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

# write output summary phase
phase_summary_list <- list()
# Collect QC phase summary for each sample
for (sample_name in names(seurat_objects)) {
  phase_summary <- table(seurat_objects[[sample_name]][[]]$Phase)
  phase_summary_df <- as.data.frame(phase_summary)
  phase_summary_df$Sample <- sample_name
  phase_summary_list[[sample_name]] <- phase_summary_df
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

# Combine all summaries into one data frame
phase_summary_combined <- do.call(rbind, phase_summary_list)
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V2_noQC.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## plot cell cycle
# Calculate proportions
phase_summary_combined_tidy <- phase_summary_combined %>%
  group_by(Sample) %>%
  mutate(Prop = Freq / sum(Freq) * 100)

phase_summary_combined_tidy$Sample <- factor(phase_summary_combined_tidy$Sample, levels = c("RNA_WT", "RNA_Bap1KO")) # Reorder untreated 1st


# Plot
pdf("output/seurat/barPlot_CellCyclePhase_V2_noQC.pdf", width=5, height=6)
ggplot(phase_summary_combined_tidy, aes(x = Sample, y = Prop, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Genotype", y = "Proportion (%)", fill = "Cell Cycle Phase") +
  theme_bw() +
  scale_fill_manual(values = c("G1" = "#1f77b4", "G2M" = "#ff7f0e", "S" = "#2ca02c")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
##




# Data integration WT and Bap1KO


RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 8698, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 2500) %>% 
    RunPCA(npcs = 20, verbose = FALSE)
RNA_Bap1KO <- SCTransform(RNA_Bap1KO, method = "glmGamPoi", ncells = 8904, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 2500) %>% 
    RunPCA(npcs = 20, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(RNA_WT = RNA_WT, RNA_Bap1KO = RNA_Bap1KO)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 2500)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

embryo.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
RNA_WT_Bap1KO.sct <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT")

set.seed(42)

DefaultAssay(RNA_WT_Bap1KO.sct) <- "integrated"

RNA_WT_Bap1KO.sct <- RunPCA(RNA_WT_Bap1KO.sct, verbose = FALSE, npcs = 20)
RNA_WT_Bap1KO.sct <- RunUMAP(RNA_WT_Bap1KO.sct, reduction = "pca", dims = 1:20, verbose = FALSE)
RNA_WT_Bap1KO.sct <- FindNeighbors(RNA_WT_Bap1KO.sct, reduction = "pca", k.param = 70, dims = 1:20)
RNA_WT_Bap1KO.sct <- FindClusters(RNA_WT_Bap1KO.sct, resolution = 0.9, verbose = FALSE, algorithm = 4)

RNA_WT_Bap1KO.sct$orig.ident <- factor(RNA_WT_Bap1KO.sct$orig.ident, levels = c("RNA_WT", "RNA_Bap1KO")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Bap1KO_noSplit-noQC_dim20kparam70res09algo4feat2500_noCellCycleRegression.pdf", width=8, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", label=TRUE)
dev.off()



DefaultAssay(RNA_WT_Bap1KO.sct) <- "SCT" # For vizualization either use SCT or norm RNA



pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-Bap1Pax6Eomes-noQC_dim16kparam30res07algo4feat2000_noCellCycleRegression
.pdf", width=10, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Bap1", "Pax6", "Eomes"), max.cutoff = 1, cols = c("grey", "red"), split.by = "orig.ident")
dev.off()

pdf("output/seurat/RidgePlot_SCT_RNA_WT_Bap1KO-Bap1Pax6Eomes-noQC_dim16kparam30res07algo4feat2000_noCellCycleRegression.pdf", width=5, height=5)
RidgePlot(RNA_WT_Bap1KO.sct, features = c("Bap1", "Pax6", "Eomes"), cols = c("blue","red"), group.by = 'orig.ident', ncol =1, assay = "SCT")
dev.off()

pdf("output/seurat/VlnPlot_SCT_RNA_WT_Bap1KO-Bap1Pax6Eomes-noQC_dim16kparam30res07algo4feat2000_noCellCycleRegression.pdf", width=3, height=10)
VlnPlot(RNA_WT_Bap1KO.sct, features = c("Bap1", "Pax6", "Eomes"), cols = c("blue","red"), group.by = 'orig.ident', ncol =1, assay = "RNA")
dev.off()



#### Check raw/norm Bap1 count #####################
RNA_WT_Bap1KO.sct_WT <- subset(RNA_WT_Bap1KO.sct, subset = orig.ident == 'RNA_WT')
RNA_WT_Bap1KO.sct_Bap1KO <- subset(RNA_WT_Bap1KO.sct, subset = orig.ident == 'RNA_Bap1KO')


AverageExpression(RNA_WT_Bap1KO.sct, features = "Bap1", group.by = "orig.ident")


AverageExpression(RNA_WT_Bap1KO.sct, features = "Bap1", group.by = "seurat_clusters")


AverageExpression(RNA_WT_Bap1KO.sct_WT , features = "Bap1", group.by = "seurat_clusters")
AverageExpression(RNA_WT_Bap1KO.sct_Bap1KO , features = "Bap1", group.by = "seurat_clusters")

####################################

```






## Condiments workflow to compare condition - WT vs Bap1KO 

- *Issue*: First test show hard to identify expected pseudotime trajectories when keeping all cell types: Some cell types are very different, spatially separated, leading to noise in the trajectory. Separating UL/ML/DP versus CA1/CA23 was challenging
--> *Solution*: Let's subset the UMAP to keep cell type within their same pseudotime trajectories (eg. when studying pseudotime traj. for CA1/CA23, we will remove cells from UL/DL)


Partitions of interest (start to end point):
- Part_DG_GC = `seurat_clusters %in% c(8,11,13,12,7)`
- Part_PyNs_RSC_UL = `seurat_clusters %in% c(8,11,13,12,17,3,5)`
- Part_SubC1 = `seurat_clusters %in% c(8,11,13,12,2,6,1,15,10)`



```bash
conda activate condiments_V6
```

```R
# package installation 
## install.packages("remotes")
## remotes::install_github("cran/spatstat.core")
## remotes::install_version("Seurat", "4.0.3")
## install.packages("magrittr")
## install.packages("magrittr")
## install.packages("dplyr")
## BiocManager::install("DelayedMatrixStats")
## BiocManager::install("tradeSeq")


# packages
library("condiments")
library("Seurat")
library("magrittr") # to use pipe
library("dplyr") # to use bind_cols and sample_frac
library("SingleCellExperiment") # for reducedDims
library("ggplot2")
library("slingshot")
library("DelayedMatrixStats")
library("tidyr")
library("tradeSeq")
library("cowplot")
library("scales")
library("pheatmap")

# Data import
RNA_WT_Bap1KO.sct <- readRDS(file = "output/seurat/RNA_WT_Bap1KO.sct_V1_numeric.rds")

DefaultAssay(RNA_WT_Bap1KO.sct) <- "RNA" # According to condiments workflow


# convert to SingleCellExperiment
RNA_WT_Bap1KO <- as.SingleCellExperiment(RNA_WT_Bap1KO.sct, assay = "RNA")

# Separate SCE object for each partitions:
Part_DG_GC <- RNA_WT_Bap1KO[, RNA_WT_Bap1KO$seurat_clusters %in% c(8,11,13,12,7)]
Part_PyNs_RSC_UL <- RNA_WT_Bap1KO[, RNA_WT_Bap1KO$seurat_clusters %in% c(8,11,13,12,17,3,5)]
Part_SubC1 <- RNA_WT_Bap1KO[, RNA_WT_Bap1KO$seurat_clusters %in% c(8,11,13,12,2,6,1,15,10)]
table(Part_PyNs_RSC_UL$seurat_clusters) # to double check



#############################################################################################
################### Time course effect, trajectory per trajectory ##############################
#############################################################################################

##$ Part_DG_GC ##$###########################################
#############################################################################################



scripts/fitGAM_6knots_Part_DG_GC_RNA_WT_Bap1KO_noCondition.sh 

set.seed(42)
traj_Part_DG_GC_noCondition <- readRDS("output/condiments/traj_Part_DG_GC_noCondition.rds")


## Genes that change with pseudotime

pseudotime_association <- associationTest(traj_Part_DG_GC_noCondition) # statistical test to check whether gene expression is constant across pseudotime within a lineage
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$fdr), ]
pseudotime_association$gene <- rownames(pseudotime_association)

pseudotime_association = as_tibble(pseudotime_association) # 9,587 fdr 0.05 DEG



# save output: write.table(pseudotime_association, file = c("output/condiments/pseudotime_association_traj_Part_DG_GC_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)
#--> Can do clustering on these genes if needed


## Genes that change between two pseudotime points (start vs end)
pseudotime_start_end_association <- startVsEndTest(traj_Part_DG_GC_noCondition, pseudotimeValues = c(0, 1))
pseudotime_start_end_association$gene <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association$fdr <- p.adjust(pseudotime_start_end_association$pvalue, method = "fdr")
pseudotime_start_end_association <- pseudotime_start_end_association[order(pseudotime_start_end_association$fdr), ]


##--> log2FC = end - start: negative log2fc means start point higher average expr than end point
# save output: write.table(pseudotime_start_end_association, file = c("output/condiments/pseudotime_start_end_association_traj_Part_DG_GC_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)

sce_cells <- colnames(traj_Part_DG_GC_noCondition) # collect cells of traj2
subset_traj_Part_DG_GC_noCondition.sct <- subset(RNA_WT_Bap1KO.sct, cells = sce_cells) # Create a seurat object with only cells from traj2

### plot gene per gene

top1_posFC = pseudotime_start_end_association %>% filter(fdr < 0.05, logFClineage1>0) %>% top_n(1, waldStat) %>% pull(gene)
top1_negFC = pseudotime_start_end_association %>% filter(fdr < 0.05, logFClineage1<0) %>% top_n(1, waldStat) %>% pull(gene)

pdf("output/condiments/plotSmoothers-traj_Part_DG_GC_noCondition_top1_posFC.pdf", width=5, height=4)
plotSmoothers(traj_Part_DG_GC_noCondition, subset_traj_Part_DG_GC_noCondition.sct[["RNA"]]@counts, gene = top1_posFC )
dev.off()

pdf("output/condiments/plotSmoothers-traj_Part_DG_GC_noCondition_top1_negFC.pdf",  width=5, height=4)
plotSmoothers(traj_Part_DG_GC_noCondition, subset_traj_Part_DG_GC_noCondition.sct[["RNA"]]@counts, gene = top1_negFC )
dev.off()


### plot top 25 genes
#### Select the top 25 genes with positive logFC and the top 25 with negative logFC
top25_posFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 > 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

top25_negFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 < 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

## Function to plot and save in PDF
plot_and_save <- function(genes_df, file_name) {
  pdf(file_name, width=5, height=4)
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$gene[i]
    logFC <- genes_df$logFClineage1[i]
    plot_title <- paste0(gene, " (logFC: ", round(logFC, 2), ")")
    p <- plotSmoothers(traj_Part_DG_GC_noCondition, subset_traj_Part_DG_GC_noCondition.sct[["RNA"]]@counts, gene = gene)
    p <- p + ggtitle(plot_title)
    print(p)
  }
  dev.off()
}



# Generate PDFs
plot_and_save(top25_posFC_genes, "output/condiments/plotSmoothers-top25_posFC_traj_Part_DG_GC_noCondition.pdf")
plot_and_save(top25_negFC_genes, "output/condiments/plotSmoothers-top25_negFC_traj_Part_DG_GC_noCondition.pdf")


## Identify Activation point = peak (maximum expression) of each gene along the pseudotime trajectory
### Identify peak of expression (max expr) of these Time-course DEG
traj_Part_DG_GC_noCondition
#### Extract pseudotime values
pseudotime <- colData(traj_Part_DG_GC_noCondition)$crv$pseudotime
#### Extract the expression matrix
expr_matrix <- assays(traj_Part_DG_GC_noCondition)$counts
#### Ensure the pseudotime values are named with the same cell names as the expression matrix columns
names(pseudotime) <- colnames(expr_matrix)
#### Function to find the peak pseudotime for each gene (raw and smoothed)
find_max_pseudotime <- function(gene_expr, pseudotime) {
  # Raw peak pseudotime
  raw_peak_pseudotime <- pseudotime[which.max(gene_expr)]
  # Smooth gene expression using loess
  smooth_model <- loess(gene_expr ~ pseudotime)
  smooth_expr <- predict(smooth_model)
  # Smooth peak pseudotime
  smooth_peak_pseudotime <- pseudotime[which.max(smooth_expr)]
  return(list(raw_peak_pseudotime = raw_peak_pseudotime, 
              smooth_peak_pseudotime = smooth_peak_pseudotime))
}
#### Apply the function to all genes
peak_values <- apply(expr_matrix, 1, function(x) find_max_pseudotime(as.numeric(x), pseudotime))
#### Convert the results to a data frame
peak_df <- data.frame(
  gene = rownames(expr_matrix),
  raw_peak_pseudotime = sapply(peak_values, `[[`, "raw_peak_pseudotime"),
  smooth_peak_pseudotime = sapply(peak_values, `[[`, "smooth_peak_pseudotime")
) %>% as_tibble()

# save output: write.table(peak_df, file = c("output/condiments/traj_Part_DG_GC_noCondition_ActivationPoint.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## heatmap activate/induced genes along pseudotime
### DEG Start End
pseudotime_start_end_association # filter log2fc >0 >1
pseudotime_start_end_association = read_tsv("output/condiments/pseudotime_start_end_association_traj_Part_DG_GC_noCondition.txt")
pseudotime_start_end_association_logFC0 = pseudotime_start_end_association %>% 
  filter(logFClineage1 > 0.5, fdr <0.05) %>%
  dplyr::select(gene) %>%
  unique()

#pdf("output/condiments/heatmap_pseudotime_start_end_association_logFClineageOver0.pdf", width=8, height=10)
pdf("output/condiments/heatmap_pseudotime_start_end_association_traj_Part_DG_GC_noCondition_logFClineageOver05fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_DG_GC_noCondition, gene = pseudotime_start_end_association_logFC0$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG time course
pseudotime_association = read_tsv("output/condiments/pseudotime_association_traj_Part_DG_GC_noCondition.txt")
pseudotime_association_deg = pseudotime_association %>%
  filter(fdr == 0)%>%
  dplyr::select(gene) %>%
  unique()


pdf("output/condiments/heatmap_pseudotime_association_traj_Part_DG_GC_noCondition_deg0pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_DG_GC_noCondition, gene = pseudotime_association_deg$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG Start End & time course
pseudotime_start_end_association_logFC0TCdeg05 = pseudotime_start_end_association %>% 
  filter(logFClineage1 > 0.5) %>%
  dplyr::select(gene) %>%
  unique() %>%
  left_join(pseudotime_association) %>%
  filter(fdr <0.05)


pdf("output/condiments/heatmap_pseudotime_start_end_association_traj_Part_DG_GC_noCondition_logFClineageOver05_DEG05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_DG_GC_noCondition, gene = pseudotime_start_end_association_logFC0TCdeg05$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()


## Identify which pseudotime value corespond to which cluster ################
pseudotime <- colData(traj_Part_DG_GC_noCondition)$crv$pseudotime
sce_cells <- colnames(traj_Part_DG_GC_noCondition)
subset_seurat <- subset(RNA_WT_Bap1KO.sct, cells = sce_cells) # Subset cell from traj2
clusters <- subset_seurat$seurat_clusters # Extract cluster information
### Combine pseudotime and cluster information into a data frame
pseudotime_cluster_df <- data.frame(
  cell = colnames(traj_Part_DG_GC_noCondition),
  pseudotime = pseudotime,
  cluster = clusters
)  %>%
  arrange(pseudotime)

switch_points <- which(diff(as.numeric(factor(pseudotime_cluster_df$cluster))) != 0) # Find the indices where the cluster changes
switch_pseudotimes <- pseudotime_cluster_df$pseudotime[switch_points] # Extract the pseudotime values at these switch points
switch_clusters_from <- pseudotime_cluster_df$cluster[switch_points]
switch_clusters_to <- pseudotime_cluster_df$cluster[switch_points + 1]
switch_df <- data.frame(
  switch_pseudotime = switch_pseudotimes,
  cluster_from = switch_clusters_from,
  cluster_to = switch_clusters_to 
) %>%
  group_by(cluster_from, cluster_to) %>%
  summarize(median_switch_pseudotime = median(switch_pseudotime), .groups = 'drop')
write.table(switch_df, file = c("output/condiments/switch_df_traj_Part_DG_GC_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)
##################################




##$ Part_PyNs_RSC_UL ##$###########################################
#############################################################################################

set.seed(42)
traj_Part_PyNs_RSC_UL_noCondition <- readRDS("output/condiments/traj_Part_PyNs_RSC_UL_noCondition.rds")


## Genes that change with pseudotime

pseudotime_association <- associationTest(traj_Part_PyNs_RSC_UL_noCondition) # statistical test to check whether gene expression is constant across pseudotime within a lineage
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$fdr), ]
pseudotime_association$gene <- rownames(pseudotime_association)

pseudotime_association = as_tibble(pseudotime_association) # 5,201 fdr 0.05 DEG



# save output: write.table(pseudotime_association, file = c("output/condiments/pseudotime_association_traj_Part_PyNs_RSC_UL_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)
#--> Can do clustering on these genes if needed


## Genes that change between two pseudotime points (start vs end)
pseudotime_start_end_association <- startVsEndTest(traj_Part_PyNs_RSC_UL_noCondition, pseudotimeValues = c(0, 1))
pseudotime_start_end_association$gene <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association$fdr <- p.adjust(pseudotime_start_end_association$pvalue, method = "fdr")
pseudotime_start_end_association <- pseudotime_start_end_association[order(pseudotime_start_end_association$fdr), ]


##--> log2FC = end - start: negative log2fc means start point higher average expr than end point
# save output: write.table(pseudotime_start_end_association, file = c("output/condiments/pseudotime_start_end_association_traj_Part_PyNs_RSC_UL_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)

sce_cells <- colnames(traj_Part_PyNs_RSC_UL_noCondition) # collect cells of traj
subset_traj_Part_PyNs_RSC_UL_noCondition.sct <- subset(RNA_WT_Bap1KO.sct, cells = sce_cells) # Create a seurat object with only cells from traj

### plot gene per gene

top1_posFC = pseudotime_start_end_association %>% filter(fdr < 0.05, logFClineage1>0) %>% top_n(1, waldStat) %>% pull(gene)
top1_negFC = pseudotime_start_end_association %>% filter(fdr < 0.05, logFClineage1<0) %>% top_n(1, waldStat) %>% pull(gene)

pdf("output/condiments/plotSmoothers-traj_Part_PyNs_RSC_UL_noCondition_top1_posFC.pdf", width=5, height=4)
plotSmoothers(traj_Part_PyNs_RSC_UL_noCondition, subset_traj_Part_PyNs_RSC_UL_noCondition.sct[["RNA"]]@counts, gene = top1_posFC )
dev.off()

pdf("output/condiments/plotSmoothers-traj_Part_PyNs_RSC_UL_noCondition_top1_negFC.pdf",  width=5, height=4)
plotSmoothers(traj_Part_PyNs_RSC_UL_noCondition, subset_traj_Part_PyNs_RSC_UL_noCondition.sct[["RNA"]]@counts, gene = top1_negFC )
dev.off()


### plot top 25 genes
#### Select the top 25 genes with positive logFC and the top 25 with negative logFC
top25_posFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 > 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

top25_negFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 < 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

## Function to plot and save in PDF
plot_and_save <- function(genes_df, file_name) {
  pdf(file_name, width=5, height=4)
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$gene[i]
    logFC <- genes_df$logFClineage1[i]
    plot_title <- paste0(gene, " (logFC: ", round(logFC, 2), ")")
    p <- plotSmoothers(traj_Part_PyNs_RSC_UL_noCondition, subset_traj_Part_PyNs_RSC_UL_noCondition.sct[["RNA"]]@counts, gene = gene)
    p <- p + ggtitle(plot_title)
    print(p)
  }
  dev.off()
}



# Generate PDFs
plot_and_save(top25_posFC_genes, "output/condiments/plotSmoothers-top25_posFC_traj_Part_PyNs_RSC_UL_noCondition.pdf")
plot_and_save(top25_negFC_genes, "output/condiments/plotSmoothers-top25_negFC_traj_Part_PyNs_RSC_UL_noCondition.pdf")


## Identify Activation point = peak (maximum expression) of each gene along the pseudotime trajectory
### Identify peak of expression (max expr) of these Time-course DEG
traj_Part_PyNs_RSC_UL_noCondition
#### Extract pseudotime values
pseudotime <- colData(traj_Part_PyNs_RSC_UL_noCondition)$crv$pseudotime
#### Extract the expression matrix
expr_matrix <- assays(traj_Part_PyNs_RSC_UL_noCondition)$counts
#### Ensure the pseudotime values are named with the same cell names as the expression matrix columns
names(pseudotime) <- colnames(expr_matrix)
#### Function to find the peak pseudotime for each gene (raw and smoothed)
find_max_pseudotime <- function(gene_expr, pseudotime) {
  # Raw peak pseudotime
  raw_peak_pseudotime <- pseudotime[which.max(gene_expr)]
  # Smooth gene expression using loess
  smooth_model <- loess(gene_expr ~ pseudotime)
  smooth_expr <- predict(smooth_model)
  # Smooth peak pseudotime
  smooth_peak_pseudotime <- pseudotime[which.max(smooth_expr)]
  return(list(raw_peak_pseudotime = raw_peak_pseudotime, 
              smooth_peak_pseudotime = smooth_peak_pseudotime))
}
#### Apply the function to all genes
peak_values <- apply(expr_matrix, 1, function(x) find_max_pseudotime(as.numeric(x), pseudotime))
#### Convert the results to a data frame
peak_df <- data.frame(
  gene = rownames(expr_matrix),
  raw_peak_pseudotime = sapply(peak_values, `[[`, "raw_peak_pseudotime"),
  smooth_peak_pseudotime = sapply(peak_values, `[[`, "smooth_peak_pseudotime")
) %>% as_tibble()

# save output: write.table(peak_df, file = c("output/condiments/traj_Part_PyNs_RSC_UL_noCondition_ActivationPoint.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## heatmap activate/induced genes along pseudotime
### DEG Start End
pseudotime_start_end_association # filter log2fc >0 >1
pseudotime_start_end_association = read_tsv("output/condiments/pseudotime_start_end_association_traj_Part_PyNs_RSC_UL_noCondition.txt")
pseudotime_start_end_association_logFC0 = pseudotime_start_end_association %>% 
  filter(logFClineage1 > 0, fdr <0.05) %>%
  dplyr::select(gene) %>%
  unique()

#pdf("output/condiments/heatmap_pseudotime_start_end_association_logFClineageOver0.pdf", width=8, height=10)
pdf("output/condiments/heatmap_pseudotime_start_end_association_traj_Part_PyNs_RSC_UL_noCondition_logFClineageOver0fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_PyNs_RSC_UL_noCondition, gene = pseudotime_start_end_association_logFC0$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG time course
pseudotime_association = read_tsv("output/condiments/pseudotime_association_traj_Part_PyNs_RSC_UL_noCondition.txt")
pseudotime_association_deg = pseudotime_association %>%
  filter(fdr == 0)%>%
  dplyr::select(gene) %>%
  unique()


pdf("output/condiments/heatmap_pseudotime_association_traj_Part_PyNs_RSC_UL_noCondition_deg0pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_PyNs_RSC_UL_noCondition, gene = pseudotime_association_deg$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG Start End & time course
pseudotime_start_end_association_logFC0TCdeg05 = pseudotime_start_end_association %>% 
  filter(logFClineage1 > 0.5) %>%
  dplyr::select(gene) %>%
  unique() %>%
  left_join(pseudotime_association) %>%
  filter(fdr <0.05)


pdf("output/condiments/heatmap_pseudotime_start_end_association_traj_Part_PyNs_RSC_UL_noCondition_logFClineageOver05_DEG05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_PyNs_RSC_UL_noCondition, gene = pseudotime_start_end_association_logFC0TCdeg05$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()


## Identify which pseudotime value corespond to which cluster ################
pseudotime <- colData(traj_Part_PyNs_RSC_UL_noCondition)$crv$pseudotime
sce_cells <- colnames(traj_Part_PyNs_RSC_UL_noCondition)
subset_seurat <- subset(RNA_WT_Bap1KO.sct, cells = sce_cells) # Subset cell from traj2
clusters <- subset_seurat$seurat_clusters # Extract cluster information
### Combine pseudotime and cluster information into a data frame
pseudotime_cluster_df <- data.frame(
  cell = colnames(traj_Part_PyNs_RSC_UL_noCondition),
  pseudotime = pseudotime,
  cluster = clusters
)  %>%
  arrange(pseudotime)

switch_points <- which(diff(as.numeric(factor(pseudotime_cluster_df$cluster))) != 0) # Find the indices where the cluster changes
switch_pseudotimes <- pseudotime_cluster_df$pseudotime[switch_points] # Extract the pseudotime values at these switch points
switch_clusters_from <- pseudotime_cluster_df$cluster[switch_points]
switch_clusters_to <- pseudotime_cluster_df$cluster[switch_points + 1]
switch_df <- data.frame(
  switch_pseudotime = switch_pseudotimes,
  cluster_from = switch_clusters_from,
  cluster_to = switch_clusters_to 
) %>%
  group_by(cluster_from, cluster_to) %>%
  summarize(median_switch_pseudotime = median(switch_pseudotime), .groups = 'drop')
write.table(switch_df, file = c("output/condiments/switch_df_traj_Part_PyNs_RSC_UL_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)
##################################







##$ Part_SubC1 ##$###########################################
#############################################################################################

set.seed(42)
traj_Part_SubC1_noCondition <- readRDS("output/condiments/traj_Part_SubC1_noCondition.rds")


## Genes that change with pseudotime

pseudotime_association <- associationTest(traj_Part_SubC1_noCondition) # statistical test to check whether gene expression is constant across pseudotime within a lineage
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$fdr), ]
pseudotime_association$gene <- rownames(pseudotime_association)

pseudotime_association = as_tibble(pseudotime_association) # xxx fdr 0.05 DEG



# save output: write.table(pseudotime_association, file = c("output/condiments/pseudotime_association_traj_Part_SubC1_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)
#--> Can do clustering on these genes if needed


## Genes that change between two pseudotime points (start vs end)
pseudotime_start_end_association <- startVsEndTest(traj_Part_SubC1_noCondition, pseudotimeValues = c(0, 1))
pseudotime_start_end_association$gene <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association$fdr <- p.adjust(pseudotime_start_end_association$pvalue, method = "fdr")
pseudotime_start_end_association <- pseudotime_start_end_association[order(pseudotime_start_end_association$fdr), ]


##--> log2FC = end - start: negative log2fc means start point higher average expr than end point
# save output: write.table(pseudotime_start_end_association, file = c("output/condiments/pseudotime_start_end_association_traj_Part_SubC1_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)

sce_cells <- colnames(traj_Part_SubC1_noCondition) # collect cells of traj
subset_traj_Part_SubC1_noCondition.sct <- subset(RNA_WT_Bap1KO.sct, cells = sce_cells) # Create a seurat object with only cells from traj

### plot gene per gene

top1_posFC = pseudotime_start_end_association %>% filter(fdr < 0.05, logFClineage1>0) %>% top_n(1, waldStat) %>% pull(gene)
top1_negFC = pseudotime_start_end_association %>% filter(fdr < 0.05, logFClineage1<0) %>% top_n(1, waldStat) %>% pull(gene)

pdf("output/condiments/plotSmoothers-traj_Part_SubC1_noCondition_top1_posFC.pdf", width=5, height=4)
plotSmoothers(traj_Part_SubC1_noCondition, subset_traj_Part_SubC1_noCondition.sct[["RNA"]]@counts, gene = top1_posFC )
dev.off()

pdf("output/condiments/plotSmoothers-traj_Part_SubC1_noCondition_top1_negFC.pdf",  width=5, height=4)
plotSmoothers(traj_Part_SubC1_noCondition, subset_traj_Part_SubC1_noCondition.sct[["RNA"]]@counts, gene = top1_negFC )
dev.off()


### plot top 25 genes
#### Select the top 25 genes with positive logFC and the top 25 with negative logFC
top25_posFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 > 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

top25_negFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 < 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

## Function to plot and save in PDF
plot_and_save <- function(genes_df, file_name) {
  pdf(file_name, width=5, height=4)
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$gene[i]
    logFC <- genes_df$logFClineage1[i]
    plot_title <- paste0(gene, " (logFC: ", round(logFC, 2), ")")
    p <- plotSmoothers(traj_Part_SubC1_noCondition, subset_traj_Part_SubC1_noCondition.sct[["RNA"]]@counts, gene = gene)
    p <- p + ggtitle(plot_title)
    print(p)
  }
  dev.off()
}



# Generate PDFs
plot_and_save(top25_posFC_genes, "output/condiments/plotSmoothers-top25_posFC_traj_Part_SubC1_noCondition.pdf")
plot_and_save(top25_negFC_genes, "output/condiments/plotSmoothers-top25_negFC_traj_Part_SubC1_noCondition.pdf")


## Identify Activation point = peak (maximum expression) of each gene along the pseudotime trajectory
### Identify peak of expression (max expr) of these Time-course DEG
traj_Part_SubC1_noCondition
#### Extract pseudotime values
pseudotime <- colData(traj_Part_SubC1_noCondition)$crv$pseudotime
#### Extract the expression matrix
expr_matrix <- assays(traj_Part_SubC1_noCondition)$counts
#### Ensure the pseudotime values are named with the same cell names as the expression matrix columns
names(pseudotime) <- colnames(expr_matrix)
#### Function to find the peak pseudotime for each gene (raw and smoothed)
find_max_pseudotime <- function(gene_expr, pseudotime) {
  # Raw peak pseudotime
  raw_peak_pseudotime <- pseudotime[which.max(gene_expr)]
  # Smooth gene expression using loess
  smooth_model <- loess(gene_expr ~ pseudotime)
  smooth_expr <- predict(smooth_model)
  # Smooth peak pseudotime
  smooth_peak_pseudotime <- pseudotime[which.max(smooth_expr)]
  return(list(raw_peak_pseudotime = raw_peak_pseudotime, 
              smooth_peak_pseudotime = smooth_peak_pseudotime))
}
#### Apply the function to all genes
peak_values <- apply(expr_matrix, 1, function(x) find_max_pseudotime(as.numeric(x), pseudotime))
#### Convert the results to a data frame
peak_df <- data.frame(
  gene = rownames(expr_matrix),
  raw_peak_pseudotime = sapply(peak_values, `[[`, "raw_peak_pseudotime"),
  smooth_peak_pseudotime = sapply(peak_values, `[[`, "smooth_peak_pseudotime")
) %>% as_tibble()

# save output: write.table(peak_df, file = c("output/condiments/traj_Part_SubC1_noCondition_ActivationPoint.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## heatmap activate/induced genes along pseudotime
### DEG Start End
pseudotime_start_end_association # filter log2fc >0 >1
pseudotime_start_end_association = read_tsv("output/condiments/pseudotime_start_end_association_traj_Part_SubC1_noCondition.txt")
pseudotime_start_end_association_logFC0 = pseudotime_start_end_association %>% 
  filter(logFClineage1 > 0, fdr <0.05) %>%
  dplyr::select(gene) %>%
  unique()

#pdf("output/condiments/heatmap_pseudotime_start_end_association_logFClineageOver0.pdf", width=8, height=10)
pdf("output/condiments/heatmap_pseudotime_start_end_association_traj_Part_SubC1_noCondition_logFClineageOver0fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_SubC1_noCondition, gene = pseudotime_start_end_association_logFC0$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG time course
pseudotime_association = read_tsv("output/condiments/pseudotime_association_traj_Part_SubC1_noCondition.txt")
pseudotime_association_deg = pseudotime_association %>%
  filter(fdr == 0)%>%
  dplyr::select(gene) %>%
  unique()


pdf("output/condiments/heatmap_pseudotime_association_traj_Part_SubC1_noCondition_deg0pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_SubC1_noCondition, gene = pseudotime_association_deg$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG Start End & time course
pseudotime_start_end_association_logFC0TCdeg05 = pseudotime_start_end_association %>% 
  filter(logFClineage1 > 0.5) %>%
  dplyr::select(gene) %>%
  unique() %>%
  left_join(pseudotime_association) %>%
  filter(fdr <0.05)


pdf("output/condiments/heatmap_pseudotime_start_end_association_traj_Part_SubC1_noCondition_logFClineageOver05_DEG05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj_Part_SubC1_noCondition, gene = pseudotime_start_end_association_logFC0TCdeg05$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()


## Identify which pseudotime value corespond to which cluster ################
pseudotime <- colData(traj_Part_SubC1_noCondition)$crv$pseudotime
sce_cells <- colnames(traj_Part_SubC1_noCondition)
subset_seurat <- subset(RNA_WT_Bap1KO.sct, cells = sce_cells) # Subset cell from traj2
clusters <- subset_seurat$seurat_clusters # Extract cluster information
### Combine pseudotime and cluster information into a data frame
pseudotime_cluster_df <- data.frame(
  cell = colnames(traj_Part_SubC1_noCondition),
  pseudotime = pseudotime,
  cluster = clusters
)  %>%
  arrange(pseudotime)

switch_points <- which(diff(as.numeric(factor(pseudotime_cluster_df$cluster))) != 0) # Find the indices where the cluster changes
switch_pseudotimes <- pseudotime_cluster_df$pseudotime[switch_points] # Extract the pseudotime values at these switch points
switch_clusters_from <- pseudotime_cluster_df$cluster[switch_points]
switch_clusters_to <- pseudotime_cluster_df$cluster[switch_points + 1]
switch_df <- data.frame(
  switch_pseudotime = switch_pseudotimes,
  cluster_from = switch_clusters_from,
  cluster_to = switch_clusters_to 
) %>%
  group_by(cluster_from, cluster_to) %>%
  summarize(median_switch_pseudotime = median(switch_pseudotime), .groups = 'drop')
write.table(switch_df, file = c("output/condiments/switch_df_traj_Part_SubC1_noCondition.txt"),sep="\t", quote=FALSE, row.names=FALSE)
##################################

















































#############################################################################################
################### Genotype effect, trajectory per trajectory ##############################
#############################################################################################

##$ Part_DG_GC ##$###########################################
#############################################################################################
# tidy
df <- bind_cols(
  as.data.frame(reducedDims(Part_DG_GC)$UMAP),
  as.data.frame(colData(Part_DG_GC)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
pdf("output/condiments/UMAP_RNA_WT_Bap1KO_Part_DG_GC.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = seurat_clusters)) +
  geom_point(size = .7) +
  labs(col = "seurat_clusters") +
  theme_classic()
dev.off()


## genotype overlap

pdf("output/condiments/UMAP_condition_RNA_WT_Bap1KO_Part_DG_GC.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = orig.ident)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("blue", "red")) + # Specify colors here
  labs(col = "Condition") +
  theme_classic()
dev.off()

## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$orig.ident,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

pdf("output/condiments/UMAP_imbalance_score_RNA_WT_Bap1KO_Part_DG_GC.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores") +
  theme_classic()
dev.off()


#  Trajectory Inference and Differential Topology
set.seed(42)

## PLOT with Separate trajectories
################ Paramater testing ########################## 
Part_DG_GC <- slingshot(Part_DG_GC, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_DG_GC)$seurat_clusters,
                 start.clus = '8', approx_points = 100)
Part_DG_GC <- slingshot(Part_DG_GC, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_DG_GC)$seurat_clusters,
                 start.clus = "8",
                 end.clus = c("7"),
                 approx_points = 200,
                 extend = 'n')
##########################################

Part_DG_GC <- slingshot(Part_DG_GC, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_DG_GC)$seurat_clusters,
                 start.clus = "8",
                 end.clus = c("7"),
                 approx_points = 100,
                 extend = 'n')


set.seed(42)
topologyTest(SlingshotDataSet(Part_DG_GC), Part_DG_GC$orig.ident) # KS_mean / 0.01 / 0.007581268 / 0.7581611


sdss <- slingshot_conditions(SlingshotDataSet(Part_DG_GC), Part_DG_GC$orig.ident)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "orig.ident")


# pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_Part_DG_GC_endApprox200extendn.pdf", width=5, height=5)

pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_Part_DG_GC_endApprox100extendnpdf", width=5, height=5)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = orig.ident)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(orig.ident, Lineage, Order),
            aes(group = interaction(Lineage, orig.ident)), size = 1.5) +
  theme_bw()
dev.off()






########## Code to label lineages if more than 1 #####################
# Add custom labels for each trajectory based on the Lineage
curves$label <- with(curves, ifelse(Lineage == 1, "Trajectory 1",
                               ifelse(Lineage == 2, "Trajectory 2", "Trajectory 3")))

pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_trajLabel.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  geom_text(data = curves %>% group_by(Lineage) %>% top_n(1, Order),
            aes(label = label, x = UMAP_1, y = UMAP_2, group = Lineage),
            size = 4, vjust = -1, hjust = 0.5) +
  theme_classic()
dev.off()

############### NEED TO MODIFY THE CODE BELOW TO ANNOTATE THE DIFFERENT TRAJECTORIES ###############
pdf("output/condiments/UMAP_trajectory_separated_trajAnnotated_RNA_WT_Bap1KO.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  annotate("text", x = -10, y = 6, label = "Lineage1", size = 5) +
  annotate("text", x = -7, y = -2.7, label = "Lineage2", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) +
  NULL
dev.off()
####### ################################################################################



## PLOT with common trajectories
df_2 <- bind_cols(
  as.data.frame(reducedDim(Part_DG_GC, "UMAP")),
  slingPseudotime(Part_DG_GC) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(Part_DG_GC) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
       #  Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
       #  pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst),
        # pst = max(pst) - pst)
)
curves <- slingCurves(Part_DG_GC, as.df = TRUE)

pdf("output/condiments/UMAP_trajectory_common_RNA_WT_Bap1KO_Part_DG_GC.pdf", width=5, height=5)

ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = Lineage1_pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1.5) +
  theme_classic()
dev.off()

########## Code to label lineages if more than 1 #####################
### With label
#### Calculate midpoint for each trajectory to place the label
curves_midpoints <- curves %>%
  group_by(Lineage) %>%
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))
pdf("output/condiments/UMAP_trajectory_common_label_RNA_WT_Bap1KO.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_midpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Add labels
  theme_classic()
dev.off()
curves_endpoints <- curves %>%
  group_by(Lineage) %>%
  arrange(Order) %>%
  top_n(1, Order) # Get the top/last ordered point for each group
pdf("output/condiments/UMAP_trajectory_common_label_RNA_WT_Bap1KO2.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_endpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
  theme_classic()
dev.off()
######################################################################


# Differential Progression
progressionTest(Part_DG_GC, conditions = Part_DG_GC$orig.ident, lineages = TRUE)

prog_res <- progressionTest(Part_DG_GC, conditions = Part_DG_GC$orig.ident, lineages = TRUE)


df_3 <-  slingPseudotime(Part_DG_GC) %>% as.data.frame() 

df_3$orig.ident <- Part_DG_GC$orig.ident
df_3 <- df_3 %>% 
  pivot_longer(-(orig.ident), names_to = "Lineage",
               values_to = "pst") %>%
  filter(!is.na(pst))

pdf("output/condiments/densityPlot_trajectory_lineages_RNA_WT_Bap1KO_Part_DG_GC.pdf", width=10, height=5)

ggplot(df_3, aes(x = pst)) +
  geom_density(alpha = .8, aes(fill = orig.ident), col = "transparent") +
  geom_density(aes(col = orig.ident), fill = "transparent", size = 1.5) +
  labs(x = "Pseudotime", fill = "orig.ident") +
  facet_wrap(~Lineage, scales = "free", nrow=2) +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("blue", "red"))
  )) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()
dev.off()


#### ->  save.image(file="output/condiments/condiments_RNA_WT_Bap1KO_Part_DG_GC.RData")
### load("output/condiments/condiments_RNA_WT_Bap1KO_Part_DG_GC.RData")
set.seed(42)

#  Differential expression

## Identify the needed number of knots 


# Run the DEGs trajectory per trajectory
### FROM THIS https://github.com/statOmics/tradeSeq/issues/64 :

#### Let's try to run the DEGs trajectory per trajectory
counts <- RNA_WT_Bap1KO.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(Part_DG_GC$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(Part_DG_GC, na = FALSE) [,1]
cellweights <- slingCurveWeights(Part_DG_GC) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

## Genotype DEG time course
Part_DG_GC <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, # part for DEG between genotype
     nknots = 6,
     sce = TRUE
   )

## DEG time course genotype together
Part_DG_GC_noConditions <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = NULL, # part for DEG between genotype
     nknots = 6,
     sce = TRUE
   )


### Worked, Estimated run 3hours which is OK !!! Let's run this in slurm job partition per partition




##$ Part_PyNs_RSC_UL #############################################
# tidy
df <- bind_cols(
  as.data.frame(reducedDims(Part_PyNs_RSC_UL)$UMAP),
  as.data.frame(colData(Part_PyNs_RSC_UL)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
pdf("output/condiments/UMAP_RNA_WT_Bap1KO_Part_PyNs_RSC_UL.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = seurat_clusters)) +
  geom_point(size = .7) +
  labs(col = "seurat_clusters") +
  theme_classic()
dev.off()


## genotype overlap

pdf("output/condiments/UMAP_condition_RNA_WT_Bap1KO_Part_PyNs_RSC_UL.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = orig.ident)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("blue", "red")) + # Specify colors here
  labs(col = "Condition") +
  theme_classic()
dev.off()

## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$orig.ident,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

pdf("output/condiments/UMAP_imbalance_score_RNA_WT_Bap1KO_Part_PyNs_RSC_UL.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores") +
  theme_classic()
dev.off()


#  Trajectory Inference and Differential Topology
set.seed(42)

## PLOT with Separate trajectories
################ Paramater testing ########################## 
Part_PyNs_RSC_UL <- slingshot(Part_PyNs_RSC_UL, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_PyNs_RSC_UL)$seurat_clusters,
                 start.clus = '8', approx_points = 100)
Part_PyNs_RSC_UL <- slingshot(Part_PyNs_RSC_UL, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_PyNs_RSC_UL)$seurat_clusters,
                 start.clus = "8",
                 end.clus = c("5"),
                 approx_points = 200,
                 extend = 'n')
##########################################

Part_PyNs_RSC_UL <- slingshot(Part_PyNs_RSC_UL, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_PyNs_RSC_UL)$seurat_clusters,
                 start.clus = "8",
                 end.clus = c("5"),
                 approx_points = 100,
                 extend = 'n')


set.seed(42)
topologyTest(SlingshotDataSet(Part_PyNs_RSC_UL), Part_PyNs_RSC_UL$orig.ident) # KS_mean / 0.01 / 0.007581268 / 0.7581611


sdss <- slingshot_conditions(SlingshotDataSet(Part_PyNs_RSC_UL), Part_PyNs_RSC_UL$orig.ident)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "orig.ident")


# pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_Part_PyNs_RSC_UL_endApprox200extendn.pdf", width=5, height=5)

pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_Part_PyNs_RSC_UL_endApprox100extendnpdf", width=5, height=5)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = orig.ident)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(orig.ident, Lineage, Order),
            aes(group = interaction(Lineage, orig.ident)), size = 1.5) +
  theme_bw()
dev.off()






########## Code to label lineages if more than 1 #####################
# Add custom labels for each trajectory based on the Lineage
curves$label <- with(curves, ifelse(Lineage == 1, "Trajectory 1",
                               ifelse(Lineage == 2, "Trajectory 2", "Trajectory 3")))

pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_trajLabel.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  geom_text(data = curves %>% group_by(Lineage) %>% top_n(1, Order),
            aes(label = label, x = UMAP_1, y = UMAP_2, group = Lineage),
            size = 4, vjust = -1, hjust = 0.5) +
  theme_classic()
dev.off()

############### NEED TO MODIFY THE CODE BELOW TO ANNOTATE THE DIFFERENT TRAJECTORIES ###############
pdf("output/condiments/UMAP_trajectory_separated_trajAnnotated_RNA_WT_Bap1KO.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  annotate("text", x = -10, y = 6, label = "Lineage1", size = 5) +
  annotate("text", x = -7, y = -2.7, label = "Lineage2", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) +
  NULL
dev.off()
####### ################################################################################



## PLOT with common trajectories
df_2 <- bind_cols(
  as.data.frame(reducedDim(Part_PyNs_RSC_UL, "UMAP")),
  slingPseudotime(Part_PyNs_RSC_UL) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(Part_PyNs_RSC_UL) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
       #  Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
       #  pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst),
        # pst = max(pst) - pst)
)
curves <- slingCurves(Part_PyNs_RSC_UL, as.df = TRUE)

pdf("output/condiments/UMAP_trajectory_common_RNA_WT_Bap1KO_Part_PyNs_RSC_UL.pdf", width=5, height=5)

ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = Lineage1_pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1.5) +
  theme_classic()
dev.off()

########## Code to label lineages if more than 1 #####################
### With label
#### Calculate midpoint for each trajectory to place the label
curves_midpoints <- curves %>%
  group_by(Lineage) %>%
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))
pdf("output/condiments/UMAP_trajectory_common_label_RNA_WT_Bap1KO.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_midpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Add labels
  theme_classic()
dev.off()
curves_endpoints <- curves %>%
  group_by(Lineage) %>%
  arrange(Order) %>%
  top_n(1, Order) # Get the top/last ordered point for each group
pdf("output/condiments/UMAP_trajectory_common_label_RNA_WT_Bap1KO2.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_endpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
  theme_classic()
dev.off()
######################################################################


# Differential Progression
progressionTest(Part_PyNs_RSC_UL, conditions = Part_PyNs_RSC_UL$orig.ident, lineages = TRUE)

prog_res <- progressionTest(Part_PyNs_RSC_UL, conditions = Part_PyNs_RSC_UL$orig.ident, lineages = TRUE)


df_3 <-  slingPseudotime(Part_PyNs_RSC_UL) %>% as.data.frame() 

df_3$orig.ident <- Part_PyNs_RSC_UL$orig.ident
df_3 <- df_3 %>% 
  pivot_longer(-(orig.ident), names_to = "Lineage",
               values_to = "pst") %>%
  filter(!is.na(pst))

pdf("output/condiments/densityPlot_trajectory_lineages_RNA_WT_Bap1KO_Part_PyNs_RSC_UL.pdf", width=10, height=5)

ggplot(df_3, aes(x = pst)) +
  geom_density(alpha = .8, aes(fill = orig.ident), col = "transparent") +
  geom_density(aes(col = orig.ident), fill = "transparent", size = 1.5) +
  labs(x = "Pseudotime", fill = "orig.ident") +
  facet_wrap(~Lineage, scales = "free", nrow=2) +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("blue", "red"))
  )) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()
dev.off()


#### ->  save.image(file="output/condiments/condiments_RNA_WT_Bap1KO_Part_PyNs_RSC_UL.RData")
### load("output/condiments/condiments_RNA_WT_Bap1KO.RData")
set.seed(42)

#  Differential expression

## Identify the needed number of knots 


# Run the DEGs trajectory per trajectory
### FROM THIS https://github.com/statOmics/tradeSeq/issues/64 :

#### Let's try to run the DEGs trajectory per trajectory
counts <- RNA_WT_Bap1KO.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(Part_PyNs_RSC_UL$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(Part_PyNs_RSC_UL, na = FALSE) [,1]
cellweights <- slingCurveWeights(Part_PyNs_RSC_UL) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj1 <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, 
     nknots = 6,
     sce = TRUE
   )



##$ Part_SubC1 #############################################
# tidy
df <- bind_cols(
  as.data.frame(reducedDims(Part_SubC1)$UMAP),
  as.data.frame(colData(Part_SubC1)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
pdf("output/condiments/UMAP_RNA_WT_Bap1KO_Part_SubC1.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = seurat_clusters)) +
  geom_point(size = .7) +
  labs(col = "seurat_clusters") +
  theme_classic()
dev.off()


## genotype overlap

pdf("output/condiments/UMAP_condition_RNA_WT_Bap1KO_Part_SubC1.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = orig.ident)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("blue", "red")) + # Specify colors here
  labs(col = "Condition") +
  theme_classic()
dev.off()

## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$orig.ident,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

pdf("output/condiments/UMAP_imbalance_score_RNA_WT_Bap1KO_Part_SubC1.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores") +
  theme_classic()
dev.off()


#  Trajectory Inference and Differential Topology
set.seed(42)

## PLOT with Separate trajectories
################ Paramater testing ########################## 
Part_SubC1 <- slingshot(Part_SubC1, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_SubC1)$seurat_clusters,
                 start.clus = '8', approx_points = 100)
Part_SubC1 <- slingshot(Part_SubC1, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_SubC1)$seurat_clusters,
                 start.clus = "8",
                 end.clus = c("10"),
                 approx_points = 200,
                 extend = 'n')
##########################################

Part_SubC1 <- slingshot(Part_SubC1, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_SubC1)$seurat_clusters,
                 start.clus = "8",
                 end.clus = c("10"),
                 approx_points = 100,
                 extend = 'n')


set.seed(42)
topologyTest(SlingshotDataSet(Part_SubC1), Part_SubC1$orig.ident) # KS_mean / 0.01 / 0.007581268 / 0.7581611


sdss <- slingshot_conditions(SlingshotDataSet(Part_SubC1), Part_SubC1$orig.ident)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "orig.ident")


# pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_Part_SubC1_endApprox200extendn.pdf", width=5, height=5)

pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_Part_SubC1_endApprox100extendnpdf", width=5, height=5)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = orig.ident)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(orig.ident, Lineage, Order),
            aes(group = interaction(Lineage, orig.ident)), size = 1.5) +
  theme_bw()
dev.off()






########## Code to label lineages if more than 1 #####################
# Add custom labels for each trajectory based on the Lineage
curves$label <- with(curves, ifelse(Lineage == 1, "Trajectory 1",
                               ifelse(Lineage == 2, "Trajectory 2", "Trajectory 3")))

pdf("output/condiments/UMAP_trajectory_separated_RNA_WT_Bap1KO_trajLabel.pdf", width=6, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  geom_text(data = curves %>% group_by(Lineage) %>% top_n(1, Order),
            aes(label = label, x = UMAP_1, y = UMAP_2, group = Lineage),
            size = 4, vjust = -1, hjust = 0.5) +
  theme_classic()
dev.off()

############### NEED TO MODIFY THE CODE BELOW TO ANNOTATE THE DIFFERENT TRAJECTORIES ###############
pdf("output/condiments/UMAP_trajectory_separated_trajAnnotated_RNA_WT_Bap1KO.pdf", width=5, height=5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(condition, Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  annotate("text", x = -10, y = 6, label = "Lineage1", size = 5) +
  annotate("text", x = -7, y = -2.7, label = "Lineage2", size = 5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank()) +
  NULL
dev.off()
####### ################################################################################



## PLOT with common trajectories
df_2 <- bind_cols(
  as.data.frame(reducedDim(Part_SubC1, "UMAP")),
  slingPseudotime(Part_SubC1) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(Part_SubC1) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
       #  Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
       #  pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst),
        # pst = max(pst) - pst)
)
curves <- slingCurves(Part_SubC1, as.df = TRUE)

pdf("output/condiments/UMAP_trajectory_common_RNA_WT_Bap1KO_Part_SubC1.pdf", width=5, height=5)

ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = Lineage1_pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1.5) +
  theme_classic()
dev.off()

########## Code to label lineages if more than 1 #####################
### With label
#### Calculate midpoint for each trajectory to place the label
curves_midpoints <- curves %>%
  group_by(Lineage) %>%
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))
pdf("output/condiments/UMAP_trajectory_common_label_RNA_WT_Bap1KO.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_midpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Add labels
  theme_classic()
dev.off()
curves_endpoints <- curves %>%
  group_by(Lineage) %>%
  arrange(Order) %>%
  top_n(1, Order) # Get the top/last ordered point for each group
pdf("output/condiments/UMAP_trajectory_common_label_RNA_WT_Bap1KO2.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1) +
  geom_text(data = curves_endpoints, aes(label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
  theme_classic()
dev.off()
######################################################################


# Differential Progression
progressionTest(Part_SubC1, conditions = Part_SubC1$orig.ident, lineages = TRUE)

prog_res <- progressionTest(Part_SubC1, conditions = Part_SubC1$orig.ident, lineages = TRUE)


df_3 <-  slingPseudotime(Part_SubC1) %>% as.data.frame() 

df_3$orig.ident <- Part_SubC1$orig.ident
df_3 <- df_3 %>% 
  pivot_longer(-(orig.ident), names_to = "Lineage",
               values_to = "pst") %>%
  filter(!is.na(pst))

pdf("output/condiments/densityPlot_trajectory_lineages_RNA_WT_Bap1KO_Part_SubC1.pdf", width=10, height=5)

ggplot(df_3, aes(x = pst)) +
  geom_density(alpha = .8, aes(fill = orig.ident), col = "transparent") +
  geom_density(aes(col = orig.ident), fill = "transparent", size = 1.5) +
  labs(x = "Pseudotime", fill = "orig.ident") +
  facet_wrap(~Lineage, scales = "free", nrow=2) +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("blue", "red"))
  )) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()
dev.off()


#### ->  save.image(file="output/condiments/condiments_RNA_WT_Bap1KO_Part_SubC1.RData")
### load("output/condiments/condiments_RNA_WT_Bap1KO_Part_SubC1.RData")
set.seed(42)

#  Differential expression

## Identify the needed number of knots 


# Run the DEGs trajectory per trajectory
### FROM THIS https://github.com/statOmics/tradeSeq/issues/64 :

#### Let's try to run the DEGs trajectory per trajectory
counts <- RNA_WT_Bap1KO.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(Part_SubC1$orig.ident) # identify conditions
#### Extract the pseudotimes and cell weights for the first lineage
pseudotimes <- slingPseudotime(Part_SubC1, na = FALSE) [,1]
cellweights <- slingCurveWeights(Part_SubC1) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]


traj1 <- fitGAM(
     counts = sub_counts, 
     pseudotime = sub_pseudotimes,
     cellWeights = sub_weights,
     conditions = sub_cond, 
     nknots = 6,
     sce = TRUE
   )




XXXXXXXXXXXX HERE BELOW NOT MODIFIED --> MODIFIED IF WANT CHECK GENOTYPE EFFECT ALONG PSEUDOTIME !

### Once slurm jobs finished for each trajectory; load traj1 <- readRDS("output/condiments/traj1.rds") OK



traj1 <- readRDS("output/condiments/traj1_RNA_WT_Bap1KO.rds.rds")
traj2 <- readRDS("output/condiments/traj2_RNA_WT_Bap1KO.rds.rds")
traj3 <- readRDS("output/condiments/traj3_RNA_WT_Bap1KO.rds.rds")


## DEGs between condition
condRes_traj5 <- conditionTest(traj5)
condRes_traj8_l2fc2 <- conditionTest(traj8, l2fc = log2(2)) # 

condRes_traj1_l2fc4 <- conditionTest(traj1, l2fc = log2(4)) # let s prefer to use this one
condRes_traj2_l2fc4 <- conditionTest(traj2, l2fc = log2(4)) # let s prefer to use this one




# Correct the pvalue with fdr

condRes_traj2_l2fc4$padj <- p.adjust(condRes_traj2_l2fc4$pvalue, "fdr")



### Save output tables
condRes_traj2_l2fc4$gene <- rownames(condRes_traj2_l2fc4) # create new column l;abel gene; as matrix before
condRes_traj2_l2fc4 <- condRes_traj2_l2fc4[, c(ncol(condRes_traj2_l2fc4), 1:(ncol(condRes_traj2_l2fc4)-1))] # just to put gene column 1st
write.table(condRes_traj2_l2fc4, file = c("output/condiments/condRes_traj2_l2fc4.txt"),sep="\t", quote=FALSE, row.names=FALSE)





# Heatmap clutering DEGs per traj _ REVISED METHOD
## import DEGs
condRes_traj1 <- read.table("output/condiments/condRes_traj1_l2fc4.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
condRes_traj1 <- read.table("output/condiments/condRes_traj1_l2fc2.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

## Isolate significant DEGs and transform into a vector
conditionGenes_traj1_vector <- condRes_traj1 %>% 
  filter(padj <= 0.05) %>%
  pull(gene)
  
# Predict smoothed values
yhatSmooth <- 
  predictSmooth(traj1, gene = conditionGenes_traj1_vector, nPoints = 50, tidy = FALSE) %>%
  log1p()

yhatSmoothScaled <- t(apply(yhatSmooth, 1, scales::rescale))
combinedData <- yhatSmoothScaled[, c(51:100, 1:50)]
# Generate heatmap with clustering
# Perform hierarchical clustering
hc <- hclust(dist(combinedData))
clusters <- cutree(hc, k=10)
# Create an annotation data frame for the rows based on cluster assignments
annotation_row <- data.frame(Cluster = factor(clusters))
# Define colors for each cluster
# 20
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen", "grey", "darkred", "darkblue", "gold", "darkgray", "lightblue", "lightgreen", "lightcoral", "lightpink", "lightcyan"))(20),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 10
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen" ))(10),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 8
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown" ))(8),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 7
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink" ))(7),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 6
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange"))(6),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 5
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple"))(5),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 4
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow" ))(4),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Generate the heatmap
pdf("output/condiments/clustered_heatmap_traj1.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_l2fc4_cl10.pdf", width=8, height=10)

pdf("output/condiments/clustered_heatmap_traj4_l2fc2_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj5_l2fc2_cl10.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj3_l2fc2_cl7.pdf", width=8, height=10)
pdf("output/condiments/clustered_heatmap_traj1_l2fc2_cl10.pdf", width=8, height=10)

pheatmap(combinedData,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Trajectory 1 - Hierarchical Clustering",
  legend = TRUE,
  cutree_rows = 10,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
)
dev.off()


# Line plots
library("reshape2")
library("stringr")
# Assuming yhatSmoothScaled contains your smoothed gene expression data
# Convert the yhatSmoothScaled data to a dataframe
df <- as.data.frame(yhatSmoothScaled)
df$Gene <- rownames(df)
# Transform the data into a long format
df_long <- melt(df, id.vars = "Gene", variable.name = "Pseudotime", value.name = "Expression")
# Attach the cluster information to the data frame
df$Cluster <- factor(clusters[df$Gene])
df_long$Cluster <- df$Cluster[match(df_long$Gene, df$Gene)]

# Extract condition column
df_long$Condition <- str_extract(df_long$Pseudotime, "condition[[:alnum:]]+")
df_long$Condition <- ifelse(str_detect(df_long$Condition, "WT"), "WT", "KO")

# Extract the point value and convert it to numeric
df_long$Updated_Pseudotime <- as.numeric(str_extract(df_long$Pseudotime, "(?<=point)\\d+"))

# Define colors for the conditions
color_map <- c("WT" = "blue", "KO" = "red")

# Plot using ggplot
pdf("output/condiments/clustered_linePlot_traj1_l2fc4_cl10.pdf", width=10, height=5)

pdf("output/condiments/clustered_linePlot_traj4_l2fc2_cl4.pdf", width=10, height=5)

ggplot(df_long, aes(x = as.numeric(Updated_Pseudotime), y = Expression, group = Gene)) + 
  geom_line(data = subset(df_long, Condition == "WT"), aes(color = Condition), alpha = 0.5) +
  geom_line(data = subset(df_long, Condition == "KO"), aes(color = Condition), alpha = 0.5) +
  scale_color_manual(values = color_map) + 
  facet_wrap(~Cluster, scales = "free_y", nrow=2) +
  theme_bw() +
  labs(title = "Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")

dev.off()

# Plot using ggplot
pdf("output/condiments/smoothed_linePlot_traj1_l2fc4_cl10_smooth.pdf", width=10, height=5)
pdf("output/condiments/smoothed_linePlot_traj2_l2fc2_cl10_smooth.pdf", width=10, height=5)

pdf("output/condiments/smoothed_linePlot_traj4_l2fc2_cl4_smooth.pdf", width=10, height=5)

ggplot(df_long, aes(x = Updated_Pseudotime, y = Expression, color = Condition)) + 
  geom_smooth(method = "loess", se = TRUE, span = 0.5) + 
  scale_color_manual(values = color_map) + 
  facet_wrap(~Cluster, scales = "free_y", nrow=2) +
  theme_bw() +
  labs(title = "Smoothed Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")
dev.off()




### Export gene list from each cluster
## Create a data frame with gene names and their respective cluster assignments
output_df <- data.frame(
  gene = rownames(combinedData),
  cluster = clusters
)

# Write the data frame to a .txt file
write.table(output_df, 
            file = "output/condiments/gene_clusters_traj1_l2fc2_cl10.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

# Check some genes individually
## FOR LINEAGE 1
counts <- RNA_WT_Bap1KO.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(RNA_WT_Bap1KO.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(RNA_WT_Bap1KO, na = FALSE) [,1]
cellweights <- slingCurveWeights(RNA_WT_Bap1KO) [,1]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj1_Aff3.pdf", width=8, height=4)
plotSmoothers(traj1, sub_counts, gene = "Aff3", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()


## FOR LINEAGE 2
counts <- RNA_WT_Bap1KO.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(RNA_WT_Bap1KO.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(RNA_WT_Bap1KO, na = FALSE) [,2]
cellweights <- slingCurveWeights(RNA_WT_Bap1KO) [,2]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj2_Gm28376.pdf", width=8, height=4)
plotSmoothers(traj2, sub_counts, gene = "Gm28376", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()

## FOR LINEAGE 3
counts <- RNA_WT_Bap1KO.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(RNA_WT_Bap1KO.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(RNA_WT_Bap1KO, na = FALSE) [,3]
cellweights <- slingCurveWeights(RNA_WT_Bap1KO) [,3]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj3_Fcrla.pdf", width=8, height=4)
plotSmoothers(traj3, sub_counts, gene = "Fcrla", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()

## FOR LINEAGE 4
counts <- RNA_WT_Bap1KO.combined.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(RNA_WT_Bap1KO.combined.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(RNA_WT_Bap1KO, na = FALSE) [,4]
cellweights <- slingCurveWeights(RNA_WT_Bap1KO) [,4]
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj4_Hist2h2ac.pdf", width=8, height=4)
plotSmoothers(traj4, sub_counts, gene = "Hist2h2ac", curvesCol = c("red","blue") ) +
scale_color_manual(values =c("red","blue"))
dev.off()



```

Let's run **fitGAM** to identify :
- Genes that show differential expression along the pseudotime trajectory (pseudotime-dependent DEGs).
- Genes that show significant differential expression between genotypes along the pseudotime trajectory (genotype-dependent DEGs).


```bash
conda activate condiments_V5

# trajectory per trajectory (all features, no parralelization) - genotype-dependent DEGs
sbatch scripts/fitGAM_6knots_Part_DG_GC_RNA_WT_Bap1KO.sh # 22719043 xxx
sbatch scripts/fitGAM_6knots_Part_PyNs_RSC_UL_RNA_WT_Bap1KO.sh # 22719116 xxx
sbatch scripts/fitGAM_6knots_Part_PyNs_SubC1_RNA_WT_Bap1KO.sh # 22719143 xxx

# trajectory per trajectory (all features, no parralelization) - pseudotime-dependent DEGs
sbatch scripts/fitGAM_6knots_Part_DG_GC_RNA_WT_Bap1KO_noCondition.sh # 22832785 ok
sbatch scripts/fitGAM_6knots_Part_PyNs_RSC_UL_RNA_WT_Bap1KO_noCondition.sh # 22833306 ok
sbatch scripts/fitGAM_6knots_Part_PyNs_SubC1_RNA_WT_Bap1KO_noCondition.sh # 22833684 fail; 22938128 ok
```

--> To remove genotype effect, specify `conditions = NULL` to `fitGam()`

- [Workshop](https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html) to identify DEG time-course: .
  - Many options but let's look 1st at *Genes that change between two pseudotime points*






## monocle3 with RNA assay 

```bash
conda activate monocle3_V1
```



```R
# Installation
library("igraph") # important to load it first so that v1.4.3 and not v1.5 is used (v1.5 is loaded via Seurat per default!)
library("monocle3")
library("Seurat")
library("SeuratWrappers")

# Data import EMBRYO

RNA_WT_Bap1KO.sct <- readRDS(file = "output/seurat/RNA_WT_Bap1KO.sct_V1_numeric.rds")
DefaultAssay(RNA_WT_Bap1KO.sct) <- "RNA" # 

# convert data to seurat object to cell_data_set
cds <- as.cell_data_set(RNA_WT_Bap1KO.sct)
cds <- cluster_cells(cds, resolution=1e-3) # Too many cluster; lead to too many trajectories
cds <- cluster_cells(cds, resolution=1e-2) # Too many cluster; lead to too many trajectories
cds <- cluster_cells(cds, resolution=1e-4) # Look good



pdf("output/monocle3/plot_cells_RNA_WT_Bap1KO_V1_numeric.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
dev.off()

pdf("output/monocle3/plot_cells_RNA_partition_WT_Bap1KO_V1_numeric.pdf", width=5, height=5)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
dev.off()


# subsetting partition
integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

# Trajectory analysis
## RAW
cds <- learn_graph(cds, use_partition = TRUE, verbose = TRUE)

pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_WT_Bap1KO_V1_numeric.pdf", width=5, height=5)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()

XXX NOT MODIFIED BELOW XXX

## COLORED BY PSEUDOTIME

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 3]))
pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_Root1_embryo_V2clust.pdf", width=5, height=5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")
dev.off()


## PLOT as seurat object
integrated.sub <- as.Seurat(cds, assay = NULL)
pdf("output/monocle3/FeaturePlot_RNA_trajectory_partition1_Root1_embryo_V2clust.pdf", width=5, height=5)
FeaturePlot(integrated.sub, "monocle3_pseudotime")
dev.off()


# Pseudotime differential genes (1 hour!!!)
cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)


## Save the output table
# write.table(cds_graph_test_results, file = "output/monocle3/cds_graph_test_results_RNA_V2clust.txt", sep="\t", row.names=TRUE, quote=FALSE)
# cds_graph_test_results <- read.table("output/monocle3/cds_graph_test_results_RNA.txt", header = TRUE, sep = "\t", row.names = 1, quote = "")

# Check expression of some genes
rowData(cds)$gene_short_name <- row.names(rowData(cds))
head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

pdf("output/monocle3/plot_cells_RNA_trajectory_partition1_Root1_TopDEGsGenes_embryo.pdf", width=5, height=5)
plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)
dev.off()

pdf("output/monocle3/FeaturePlot_RNA_trajectory_partition1_Root1_TopDEGsGenes_embryo.pdf", width=10, height=12)
FeaturePlot(integrated.sub, features = head(deg_ids), max.cutoff = 10, cols = c("grey", "red"))
dev.off()  


## Plot gene as a function of pseudotime
#### Import DEGs
cds_graph_test_results <- read.table("output/monocle3/cds_graph_test_results_RNA.txt", header = TRUE, sep = "\t", row.names = 1, quote = "")
rowData(cds)$gene_short_name <- row.names(rowData(cds))
head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))
### Identify top 10 morans_I genes
ordered_results <- cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE), ]
top_morans_I <- rownames(ordered_results)[1:10]

partition_1_cds <- as.cell_data_set(integrated.sub, assay = "RNA")
partition_1_cds <- estimate_size_factors(partition_1_cds)


partition_1_cds <- cds[rowData(cds)$gene_short_name %in% top_morans_I,
                       colData(cds)$monocle3_partitions %in% c("1")]
partition_1_cds <- order_cells(partition_1_cds)

partition_1_cds <- partition_1_cds[,Matrix::colSums(exprs(partition_1_cds)) != 0]
partition_1_cds <- estimate_size_factors(partition_1_cds)


pdf("output/monocle3/plot_genes_in_pseudotime_trajectory_partition1_top_morans_I_embryo.pdf", width=10, height=12)
plot_genes_in_pseudotime(partition_1_cds,
                         color_cells_by="cluster.annot",
                         min_expr=0.5)
dev.off()

### Identify top 10 qvalue genes
ordered_results <- cds_graph_test_results[order(cds_graph_test_results$q_value, decreasing = FALSE), ]
top_qvalue <- rownames(ordered_results)[1:10]

partition_1_cds <- as.cell_data_set(integrated.sub, assay = "RNA")
partition_1_cds <- estimate_size_factors(partition_1_cds)


partition_1_cds <- cds[rowData(cds)$gene_short_name %in% top_qvalue,
                       colData(cds)$monocle3_partitions %in% c("1")]
partition_1_cds <- order_cells(partition_1_cds)

partition_1_cds <- partition_1_cds[,Matrix::colSums(exprs(partition_1_cds)) != 0]
partition_1_cds <- estimate_size_factors(partition_1_cds)


pdf("output/monocle3/plot_genes_in_pseudotime_trajectory_partition1_top_qvalue_embryo.pdf", width=10, height=12)
plot_genes_in_pseudotime(partition_1_cds,
                         color_cells_by="cluster.annot",
                         min_expr=0.5)
dev.off()


# GEnerate heatmap



```

--> See notes from `002*/003*` for information about the *Monocle3* tool

--> Let's prefer using Slingshot within Condiments, to identify and fine tune pseudotime trajectories. It will allow us to keep the same trajectory for the analysis to identify milestone (=major cell stage during cell type diff. progression) AND differences between genotype.




# ATACseq integration


Let's use the [WNN](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac) integration method from seurat


## Install Signac for ATACseq

Create a new conda env as scRNAseqV3 to install *Signac*: `conda create --name scRNAseqV4 --clone scRNAseqV2`

```bash
conda activate scRNAseqV4
```


```R
install.packages("Signac")
# Fail


```

Try installing through Conda:

```bash
conda activate scRNAseqV4
conda install bioconda::r-signac 
```

--> Also fail... Env deleted `conda env remove --name scRNAseqV4`

Let's create a new conda environment specifically for using Signac. Follow recommendation from [Signac](https://stuartlab.org/signac/articles/install#current-release). Let's use R v4.0.3 as in the [paper](https://www.nature.com/articles/s41592-021-01282-5#Sec9)


```bash
conda create --name Signac r-base=4.0.3
conda activate Signac
```

```R
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
#--> many errors... of package dependendy, let's isntall them 1by1

install.packages("Matrix") #  Fail tried anaconda and work!


```

--> Matrix R package need R>4.4... So let's install it thorugh conda as it should use an older verison of Matrix: `conda install conda-forge::r-matrix`
  --> Give another try on R Signac; Fail again... Need SeuratObject
  --> Let's tyr installing Seurat through conda `conda install bioconda::r-seurat`
  --> Fail again: Env deleted `conda env remove --name Signac`


Let's give another try with starting from the most recent version of R...
 

 
```bash
conda create --name SignacV1 r-base=4.4.1
conda activate SignacV1
```

```R
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
#--> work!!! But need igraph and leiden to use isntall Seurat...

install.packages('Seurat')
#--> fail for igraph and leiden dependencies

install.packages('igraph')
#--> fail...

install.packages("remotes") # work
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
#--> Fail, try a previous Seurat version
remotes::install_version("Seurat", version = "4.3.0")
#--> Fail, same with igraph and leiden dependencies

```

--> Let's try install Seurat through conda `conda install bioconda::r-seurat`
  --> Fail... Let's try install igraph with conda  `conda install conda-forge::r-igraph`
  --> Fail... 



Let's give another try installing Seurat and Signac together with remotes, from [here](https://satijalab.org/seurat/articles/install.html). Maybe Seurat 1st then Signac


 
```bash
conda create --name SignacV2 r-base=4.4.1
conda activate SignacV2
```

```R
install.packages("remotes")
remotes::install_github("stuart-lab/signac", "seurat4", quiet = TRUE) 
#--> Fail, try install Seurat then Signac
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
#--> Fail, same igraph and leiden...
install.packages("Seurat")
```
--> Fail again....


Let's give another try installing Signac from anaconda, discuss [here](https://github.com/stuart-lab/signac/issues/637)
```bash
conda activate SignacV2
conda install bioconda::r-signac 
```
--> Fail again....

Let's try clone my deseq2 env and install Signac, then Seurat (eg. cloning deseq2 env, is what unlock me to install Seurat initially...)

```bash
conda create --name SignacV3 --clone deseq2

conda activate SignacV3
```

```R
install.packages("devtools")
install.packages("Seurat")
# --> Fail, SeuratObject not avail
install.packages("SeuratObject")
# --> Fail, Matrix version..., lets try install seurat object trought Anaconda... `conda install conda-forge::r-seuratobject`: FAIL


```


Let's try clone `scRNAseq` and install Signac: 


```bash
conda create --name SignacV4 --clone scRNAseq
conda activate SignacV4
```


```R
install.packages("Signac")
```

Fail.... Remvoe all previous env

--> Try install Signac and Seurat together at same time with Conda

```bash
conda create -n SignacV5 -c conda-forge -c bioconda r-base r-signac r-seurat
conda install conda-forge::r-reticulate # needed for RNA clustering

conda install -c conda-forge leidenalg python-igraph pandas umap-learn # to use algorithm4 Leiden/Louvain clustering

conda install conda-forge::r-metap # to run FindConservedMarkers()
```
--> WORK!!!! But need tidyverse and biocmanager for genome annotations...

```R
install.packages("hdf5r") # to read `.h5` files
install.packages("tidyverse") # fail, ragg dependency
install.packages("ragg") # fail, try thorugh conda  `conda install conda-forge::r-ragg`
install.packages("tidyverse") #  OK

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("biovizBase")
```
--> WORK!!!! 


## Run Signac


Step:
- Import RNA seurat object
- Import ATAC
- Clean ATACseq with [Signac](https://stuartlab.org/signac/articles/pbmc_vignette.html#pre-processing-workflow)
- Add in the ATAC-seq data as a second assay usnig [WWN](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac)


Great discussion on data integration [here](https://github.com/satijalab/seurat/issues/5346)


```bash
conda activate SignacV5
module load hdf5_18/1.8.21 # to read `.h5` files
```


```R
# library
library("Signac")
library("Seurat")
#library("hdf5r")
library("tidyverse")
library("EnsDb.Hsapiens.v86") # hg38
library("EnsDb.Mmusculus.v79") # mm10

# PRepare genome annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# Import RNA 
RNA_WT_Bap1KO.sct <- readRDS(file = "output/seurat/RNA_WT_Bap1KO.sct_V1_numeric.rds")

## WT QC cleaning ##
# Import ATAC
ATAC_WT_counts <- Read10X_h5("ATAC_WT/outs/filtered_peak_bc_matrix.h5")
ATAC_WT_metada = read.csv(
  file = "ATAC_WT/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

# Create seurat chromatin assay
chrom_assay_WT <- CreateChromatinAssay(
   counts = ATAC_WT_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = "ATAC_WT/outs/fragments.tsv.gz",
   min.cells = 10,
   annotation = annotations
 )

ATAC_WT <- CreateSeuratObject(
  counts = chrom_assay_WT,
  assay = "peaks",
  meta.data = ATAC_WT_metada
)


# QC metrics
## compute nucleosome signal score per cell
ATAC_WT <- NucleosomeSignal(object = ATAC_WT)
## compute TSS enrichment score per cell
ATAC_WT <- TSSEnrichment(object = ATAC_WT, fast = FALSE)
## add blacklist ratio and fraction of reads in peaks
ATAC_WT$pct_reads_in_peaks <- ATAC_WT$peak_region_fragments / ATAC_WT$passed_filters * 100
ATAC_WT$blacklist_ratio <- ATAC_WT$blacklist_region_fragments / ATAC_WT$peak_region_fragments


### QC plots
pdf("output/Signac/DensityScatter_ATAC_WT.pdf", width=5, height=5)
DensityScatter(ATAC_WT, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

pdf("output/Signac/VlnPlot_QC_ATAC_WT.pdf", width=12, height=6)
VlnPlot(
  object = ATAC_WT,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()


ATAC_WT$high.tss <- ifelse(ATAC_WT$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/Signac/TSSPlot_QC_ATAC_WT.pdf", width=12, height=6)
TSSPlot(ATAC_WT, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

ATAC_WT$nucleosome_group <- ifelse(ATAC_WT$nucleosome_signal > 1.2, 'NS > 1.2', 'NS < 1.2') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/Signac/FragmentHistogram_QC_ATAC_WT.pdf", width=12, height=6)
FragmentHistogram(object = ATAC_WT, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).

pdf("output/Signac/VlnPlot_QC_ATAC_WT_nCount_peaks.pdf", width=12, height=6)
VlnPlot(
  object = ATAC_WT,
  features = c('nCount_peaks'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 50000)
dev.off()

## subset cells that pass QC
ATAC_WT_QCV1 <- subset(
  x = ATAC_WT,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 50000 &
    pct_reads_in_peaks > 20 &
    nucleosome_signal < 1.2 &
    TSS.enrichment > 2
)

ATAC_WT_QCV1
ATAC_WT_QCV1$orig.ident <- "ATAC_WT"


## Bap1ko QC cleaning ##
# Import ATAC
ATAC_Bap1KO_counts <- Read10X_h5("ATAC_Bap1KO/outs/filtered_peak_bc_matrix.h5")
ATAC_Bap1KO_metada = read.csv(
  file = "ATAC_Bap1KO/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

# Create seurat chromatin assay
chrom_assay_Bap1KO <- CreateChromatinAssay(
   counts = ATAC_Bap1KO_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = "ATAC_Bap1KO/outs/fragments.tsv.gz",
   min.cells = 10,
   annotation = annotations
 )

ATAC_Bap1KO <- CreateSeuratObject(
  counts = chrom_assay_Bap1KO,
  assay = "peaks",
  meta.data = ATAC_Bap1KO_metada
)


# QC metrics
## compute nucleosome signal score per cell
ATAC_Bap1KO <- NucleosomeSignal(object = ATAC_Bap1KO)
## compute TSS enrichment score per cell
ATAC_Bap1KO <- TSSEnrichment(object = ATAC_Bap1KO, fast = FALSE)
## add blacklist ratio and fraction of reads in peaks
ATAC_Bap1KO$pct_reads_in_peaks <- ATAC_Bap1KO$peak_region_fragments / ATAC_Bap1KO$passed_filters * 100
ATAC_Bap1KO$blacklist_ratio <- ATAC_Bap1KO$blacklist_region_fragments / ATAC_Bap1KO$peak_region_fragments


### QC plots
pdf("output/Signac/DensityScatter_ATAC_Bap1KO.pdf", width=5, height=5)
DensityScatter(ATAC_Bap1KO, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

pdf("output/Signac/VlnPlot_QC_ATAC_Bap1KO.pdf", width=12, height=6)
VlnPlot(
  object = ATAC_Bap1KO,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()


ATAC_Bap1KO$high.tss <- ifelse(ATAC_Bap1KO$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/Signac/TSSPlot_QC_ATAC_Bap1KO.pdf", width=12, height=6)
TSSPlot(ATAC_Bap1KO, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

ATAC_Bap1KO$nucleosome_group <- ifelse(ATAC_Bap1KO$nucleosome_signal > 1.2, 'NS > 1.2', 'NS < 1.2') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/Signac/FragmentHistogram_QC_ATAC_Bap1KO.pdf", width=12, height=6)
FragmentHistogram(object = ATAC_Bap1KO, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).

pdf("output/Signac/VlnPlot_QC_ATAC_Bap1KO_nCount_peaks.pdf", width=12, height=6)
VlnPlot(
  object = ATAC_Bap1KO,
  features = c('nCount_peaks'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 50000)
dev.off()

## subset cells that pass QC
ATAC_Bap1KO_QCV1 <- subset(
  x = ATAC_Bap1KO,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 50000 &
    pct_reads_in_peaks > 20 &
    nucleosome_signal < 1.2 &
    TSS.enrichment > 2
)

ATAC_Bap1KO_QCV1
ATAC_Bap1KO_QCV1$orig.ident <- "ATAC_Bap1KO"

############## saveRDS #################################################################
# saveRDS(ATAC_WT_QCV1, file = "output/Signac/ATAC_WT_QCV1.rds") 
# saveRDS(ATAC_Bap1KO_QCV1, file = "output/Signac/ATAC_Bap1KO_QCV1.rds") 
##############################################################################
set.seed(42)

ATAC_WT_QCV1 = readRDS(file = "output/Signac/ATAC_WT_QCV1.rds")
ATAC_Bap1KO_QCV1 = readRDS(file = "output/Signac/ATAC_Bap1KO_QCV1.rds")


# Data integration RNA ATAC



# Method1: Combine WT and Bap1KO using the shared DNA accessibility assay
## compute LSI
ATAC_WT_QCV1 <- FindTopFeatures(ATAC_WT_QCV1, min.cutoff = 10)
ATAC_WT_QCV1 <- RunTFIDF(ATAC_WT_QCV1)
ATAC_WT_QCV1 <- RunSVD(ATAC_WT_QCV1)
ATAC_Bap1KO_QCV1 <- FindTopFeatures(ATAC_Bap1KO_QCV1, min.cutoff = 10)
ATAC_Bap1KO_QCV1 <- RunTFIDF(ATAC_Bap1KO_QCV1)
ATAC_Bap1KO_QCV1 <- RunSVD(ATAC_Bap1KO_QCV1)
## merge/combined dataset
ATAC_WT_Bap1KO_combined.method1 <- merge(ATAC_WT_QCV1, ATAC_Bap1KO_QCV1) # merge dataset
## process the combined dataset
ATAC_WT_Bap1KO_combined.method1 <- FindTopFeatures(ATAC_WT_Bap1KO_combined.method1, min.cutoff = 10)
ATAC_WT_Bap1KO_combined.method1 <- RunTFIDF(ATAC_WT_Bap1KO_combined.method1)
ATAC_WT_Bap1KO_combined.method1 <- RunSVD(ATAC_WT_Bap1KO_combined.method1)
ATAC_WT_Bap1KO_combined.method1 <- RunUMAP(ATAC_WT_Bap1KO_combined.method1, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(ATAC_WT_Bap1KO_combined.method1, group.by = "orig.ident")
## integration WT Bap1KO
### find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(ATAC_WT_QCV1, ATAC_Bap1KO_QCV1),
  anchor.features = rownames(ATAC_WT_QCV1),
  reduction = "rlsi",
  dims = 2:30
)
### integrate LSI embeddings
ATAC_WT_Bap1KO_integrated.method1 <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = ATAC_WT_Bap1KO_combined.method1[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30,
  k.weight = 32    # Default 100; Added after an error message that recommended me to reduce below 33
)

### create a new UMAP using the integrated embeddings
ATAC_WT_Bap1KO_integrated.method1 <- RunUMAP(ATAC_WT_Bap1KO_integrated.method1, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(ATAC_WT_Bap1KO_integrated.method1, group.by = "orig.ident")

pdf("output/Signac/ATAC_WT_Bap1KO_mergedVsintegrated.method1.pdf", width=12, height=6)
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
dev.off()

############## saveRDS #################################################################
# saveRDS(ATAC_WT_Bap1KO_integrated.method1, file = "output/Signac/ATAC_WT_Bap1KO_integrated.method1.rds")
# saveRDS(ATAC_WT_Bap1KO_combined.method1, file = "output/Signac/ATAC_WT_Bap1KO_combined.method1.rds")
##############################################################################
ATAC_WT_Bap1KO_integrated.method1 = readRDS(file = "output/Signac/ATAC_WT_Bap1KO_integrated.method1.rds")
ATAC_WT_Bap1KO_combined.method1 = readRDS(file = "output/Signac/ATAC_WT_Bap1KO_combined.method1.rds")


# RNA integration_method1 #####
## Make sure activate assay to RNA and Peaks respectively for RNA and ATAC
DefaultAssay(RNA_WT_Bap1KO.sct) <- "RNA"
DefaultAssay(ATAC_WT_Bap1KO_integrated.method1) <- "peaks"


RNA_ATAC_WT_Bap1KO <- FindMultiModalNeighbors(reduction.list = list(RNA_WT_Bap1KO.sct, ATAC_WT_Bap1KO_integrated.method1), dims.list = list(1:20, 2:30)) # Indicate nb dims to use (for RNA we used 20; 30 for ATAC)
# --> Too long to run in interactive. Let's run a slurm job FindMultiModalNeighbors_RNA_ATAC_WT_Bap1KO_V1.sh


# Load rds objet
XXX HERE slurm job XXX

RNA_ATAC_WT_Bap1KO = readRDS(file = "output/Signac/FindMultiModalNeighbors_RNA_ATAC_WT_Bap1KO_V1.rds")

#
RNA_ATAC_WT_Bap1KO <- RunUMAP(RNA_ATAC_WT_Bap1KO, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
RNA_ATAC_WT_Bap1KO <- FindClusters(RNA_ATAC_WT_Bap1KO, graph.name = "wsnn", algorithm = 3, verbose = TRUE)







## Testing Area ########

DefaultAssay(ATAC_Bap1KO_QCV1) <- "peaks"
ATAC_Bap1KO_QCV1 <- RunTFIDF(ATAC_Bap1KO_QCV1)
ATAC_Bap1KO_QCV1 <- FindTopFeatures(ATAC_Bap1KO_QCV1, min.cutoff = 'q0')
ATAC_Bap1KO_QCV1 <- RunSVD(ATAC_Bap1KO_QCV1)
ATAC_Bap1KO_QCV1 <- RunUMAP(ATAC_Bap1KO_QCV1, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # We exclude the first dimension as this is typically correlated with sequencing depth

pdf("output/Signac/DepthCor_QC_ATAC_Bap1KO.pdf", width=12, height=6)
DepthCor(ATAC_Bap1KO_QCV1)
dev.off()
####################




```

- *NOTE: RunTFIDF normalizes across cells to correct for differences in cellular sequencing depth, and across peaks to give higher values to more rare peaks.*


Integration of scATACseq conditions:
  - Method1: using the shared DNA accessibility assay. Discuss [here](https://stuartlab.org/signac/articles/integrate_atac.html)
  - Method2: simply merge using  `merge(x = ATAC_Bap1KO_QCV1, y = ATAC_WT_QCV1, add.cell.ids = c("Bap1KO", "WT"))` and process the combined dataset and put with RNA
https://stuartlab.org/signac/articles/integrate_atac.html

--> Test show that Method1 is better



Slurm job for `FindMultiModalNeighbors()`

```bash

conda activate SignacV5
module load hdf5_18/1.8.21 # to read `.h5` files


sbatch scripts/FindMultiModalNeighbors_RNA_ATAC_WT_Bap1KO_V1.sh # 23485732 xxx

```


--> I should put *RNA and ATAC assay in the same Seurat object*, using cell name. I cannot do it now, as I perform condition integration and data pre-processing separately on RNA and ATAC assay. Some of the step rename the cells, thus no more common cell name...
  --> I thus need to *re-perform all pre-processing by immediatly putting together RNA and ATAC* assay. I need to use the same clustering parameter as I used earlier.
  --> As reccomended in Seurat [WWN pipeline](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac)




# RNA and ATAC integration - Start from preprocess files

Treat RNA and ATAC assay into the same Seurat object; pre-process individually RNA and ATAC by indicating the DefaultAssay to use.

- Import RNA (import all parts from previous code at `# Seurat analysis` from data soupX import to Cell cycle normalization)
- Import ATAC
- Combine into same Seurat object
- Filter cell that pass ATAC QC (use same parameter as defined earlier)
- Integrate RNA data and perform clustering step
- Integrate ATAC data and perform clustering step
- Put together RNA and ATAC


Follow the [WNN Seurat vignette](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac)




```bash
conda activate SignacV5
```


```R
set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session...
library("tidyverse")
library("EnsDb.Hsapiens.v86") # hg38
library("EnsDb.Mmusculus.v79") # mm10




# PRepare genome annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


################################################
############# Import ATAC ########################
################################################
## WT
ATAC_WT_counts <- Read10X_h5("ATAC_WT/outs/filtered_peak_bc_matrix.h5")
ATAC_WT_metada = read.csv(
  file = "ATAC_WT/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
chrom_assay_WT <- CreateChromatinAssay(
   counts = ATAC_WT_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = "ATAC_WT/outs/fragments.tsv.gz",
   min.cells = 10,
   annotation = annotations
 )
ATAC_WT <- CreateSeuratObject(
  counts = chrom_assay_WT,
  assay = "peaks",
  meta.data = ATAC_WT_metada
)
## Bap1KO



########################################################################
############# Import RNA (from `# Seurat analysis`) ############
########################################################################
RNA_WT <- readRDS(file = "output/seurat/RNA_WT_QCV2.rds")
RNA_Bap1KO <- readRDS(file = "output/seurat/RNA_Bap1KO_QCV2.rds")

atac_wt_names = colnames(ATAC_WT)
rna_wt_names = colnames(RNA_WT)


intersect(atac_wt_names, rna_wt_names)


```


--> Fail... Cell names are different from the begining!! Need to use cellranger-arc to maintain cell name...




# RNA and ATAC integration _ Start from counting

Follow the great discussion on data integration [here](https://github.com/satijalab/seurat/issues/5346).

To maintain cell name identity. We need to re-perform the counting using `cellranger-arc count`
- Count individually for each sample --> `### Count`
- Run doublet and soupX on multiome output
- Import RNA and ATAC and create a single Seurat object
- Confirm that cell name are identical between RNA and ATAC
- For RNA, Perform QC step and genotype integration as in `# Seurat analysis` 
- For ATAC, Perform QC step and genotype integration as in `## Run Signac` 
- Follow the WWN integration method to put together ATAC and RNA and generate WNN map


## RNA contamination and doublet detection

**doublet detection**
- doublet detection using [scrublet](https://github.com/swolock/scrublet) **on the filtered matrix**
- ambient RNA correction using `soupX` in R before generating the Seurat object


```bash
conda deactivate # base environment needed

sbatch scripts/scrublet_multiome_WT.sh # 23918901 ok
sbatch scripts/scrublet_multiome_Bap1KO.sh # 23918922 ok
```


Doublet detection score (YAP1scRNAseq ranged between 0.1 to 42%); here ~27%
--> Table in `006__Kim/QC_metrics_all.xlsx`

--> Successfully assigned doublet, but weird very few doublet for WT this time... Not normal.. Lets change the treshold. For homogeneity let's change it manually for both WT and Bap1KO, cheking at the histogram plot; 0.15 look best, clear decline for both samples.


```bash
conda deactivate # base environment needed

sbatch scripts/scrublet_multiome_WT_tresh015.sh # 23932232 ok
sbatch scripts/scrublet_multiome_Bap1KO_tresh015.sh # 23932250 ok

awk '$3 == "True"' output/doublets/multiome_WT_tresh015.tsv | head
```





Now ambiant RNA contamination correction with soupX



**SoupX RNA cleaning**


```bash
conda activate scRNAseq
```



```R
# install.packages('SoupX')
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2") # for human mouse gene conversion for cell cycle genes

# soupX decontamination
## Decontaminate one channel of 10X data mapped with cellranger
sc = load10X('multiome_Bap1KO/outs')   # CHANGE FILE NAME HERE

## Assess % of conta
pdf("output/soupX/autoEstCont_multiome_Bap1KO.pdf", width=10, height=10)   # CHANGE FILE NAME HEREs
sc = autoEstCont(sc) 
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/multiome_Bap1KO.RData") # CHANGE FILE NAME HERE

```

--> 0.01% conta



## Seurat/Signac analysis

### Seurat part


**Option1**:
- Import RNA and perform QC, work in `scRNAseqV2` conda env
- Import ATAC and put with RNA in same seurat object, work in `SignacV5` conda env



```bash
conda activate scRNAseqV2
```



```R
set.seed(42)

# packages
library("SoupX")
library("Seurat")
library("tidyverse")
library("dplyr")
library("Seurat")
library("patchwork")
library("sctransform")
library("glmGamPoi")
library("celldex")
library("SingleR")
library("gprofiler2")


################################################
############# Import RNA ########################
################################################

## Load the matrix and Create SEURAT object
samples <- list(
  multiome_WT = "output/soupX/multiome_WT.RData",
  multiome_Bap1KO = "output/soupX/multiome_Bap1KO.RData"
)

seurat_objects <- list()

for (sample_name in names(samples)) {
  load(samples[[sample_name]])
  seurat_objects[[sample_name]] <- CreateSeuratObject(counts = out, project = sample_name)
}

## Function to assign Seurat objects to variables (unlist the list)
assign_seurat_objects <- function(seurat_objects_list) {
  for (sample_name in names(seurat_objects_list)) {
    assign(sample_name, seurat_objects_list[[sample_name]], envir = .GlobalEnv)
  }
}
assign_seurat_objects(seurat_objects) # This apply the function


# QUALITY CONTROL
# Function to add mitochondrial and ribosomal content
add_quality_control <- function(seurat_object) {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^Rp[sl]")
  return(seurat_object)
}
seurat_objects <- lapply(seurat_objects, add_quality_control) # This apply the function
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list


# Function to add doublet information
add_doublet_information <- function(sample_name, seurat_object) {
  doublet_file <- paste0("output/doublets/", sample_name,"_tresh015", ".tsv")
  doublets <- read.table(doublet_file, header = FALSE, row.names = 1)
  colnames(doublets) <- c("Doublet_score", "Is_doublet")
  seurat_object <- AddMetaData(seurat_object, doublets)
  seurat_object$Doublet_score <- as.numeric(seurat_object$Doublet_score)
  return(seurat_object)
}
## Apply the function to each Seurat object in the list
for (sample_name in names(seurat_objects)) {
    seurat_objects[[sample_name]] <- add_doublet_information(sample_name, seurat_objects[[sample_name]])
  }

assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list



# QC filtering _ V2

### V2 not super stringeant; mit > 5 and RNAfeature 50; rb >10
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 50 & seurat_object@meta.data$QC == 'Pass', 
                                  'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 50 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', 
                                  paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 5 & seurat_object@meta.data$QC == 'Pass', 
                                  'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 5 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', 
                                  paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 10 & seurat_object@meta.data$QC == 'Pass', 
                                  'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 10 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', 
                                  paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list





### V3 more  stringeant; mit > 3 and RNAfeature 400; rb >5; nCount 400
apply_qc <- function(seurat_object) {
  # Mark doublets
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  # Mark low nFeature_RNA
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC == 'Pass', 
                                  'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', 
                                  paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  # Mark high nCount_RNA
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nCount_RNA < 400 & seurat_object@meta.data$QC == 'Pass', 
                                  'Low_nCount_RNA', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nCount_RNA < 400 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nCount_RNA', 
                                  paste('Low_nCount_RNA', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  # Mark high mitochondrial percentage
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 3 & seurat_object@meta.data$QC == 'Pass', 
                                  'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 3 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', 
                                  paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  # Mark high ribosomal percentage
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 5 & seurat_object@meta.data$QC == 'Pass', 
                                  'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 5 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', 
                                  paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}

for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list








#### Write QC summary
qc_summary_list <- list()
# Collect QC summary for each sample
for (sample_name in names(seurat_objects)) {
  qc_summary <- table(seurat_objects[[sample_name]][['QC']])
  qc_summary_df <- as.data.frame(qc_summary)
  qc_summary_df$Sample <- sample_name
  qc_summary_list[[sample_name]] <- qc_summary_df
}
qc_summary_combined <- do.call(rbind, qc_summary_list)
# Write the data frame to a tab-separated text file
write.table(qc_summary_combined, file = "output/seurat/QC_summary_V2_multiomeOption1tresh015_QCV3.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



## subset seurat object to keep cells that pass the QC
subset_qc <- function(seurat_object) {
  seurat_object <- subset(seurat_object, subset = QC == 'Pass')
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- subset_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

# Normalize and scale data, then run cell cycle sorting
set.seed(42)
## Load gene marker of cell type

mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

# Function to normalize, scale data, and perform cell cycle scoring
process_seurat_object <- function(seurat_object, mmus_s, mmus_g2m) {
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)  # zero-centres and scales it
  seurat_object <- CellCycleScoring(seurat_object, s.features = mmus_s, g2m.features = mmus_g2m)  # cell cycle sorting
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- process_seurat_object(seurat_objects[[sample_name]], mmus_s, mmus_g2m)
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

# write output summary phase
phase_summary_list <- list()
# Collect QC phase summary for each sample
for (sample_name in names(seurat_objects)) {
  phase_summary <- table(seurat_objects[[sample_name]][[]]$Phase)
  phase_summary_df <- as.data.frame(phase_summary)
  phase_summary_df$Sample <- sample_name
  phase_summary_list[[sample_name]] <- phase_summary_df
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

# Combine all summaries into one data frame
phase_summary_combined <- do.call(rbind, phase_summary_list)
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V2_multiomeOption1_QCV3.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## plot cell cycle
# Calculate proportions
phase_summary_combined_tidy <- phase_summary_combined %>%
  group_by(Sample) %>%
  mutate(Prop = Freq / sum(Freq) * 100)

phase_summary_combined_tidy$Sample <- factor(phase_summary_combined_tidy$Sample, levels = c("multiome_WT", "multiome_Bap1KO")) # Reorder untreated 1st


# Plot
pdf("output/seurat/barPlot_CellCyclePhase_V2_multiomeOPtion1_QCV3.pdf", width=5, height=6)
ggplot(phase_summary_combined_tidy, aes(x = Sample, y = Prop, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Genotype", y = "Proportion (%)", fill = "Cell Cycle Phase") +
  theme_bw() +
  scale_fill_manual(values = c("G1" = "#1f77b4", "G2M" = "#ff7f0e", "S" = "#2ca02c")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
##


###########################################################################
# saveRDS(multiome_WT, file = "output/seurat/multiome_WT_QCV3_Option1.rds") 
# saveRDS(multiome_Bap1KO, file = "output/seurat/multiome_Bap1KO_QCV3_Option1.rds") 
###########################################################################





```

### Signac part


Now let's change conda env to import ATAC file onto our RNA seurat object


```bash
conda activate SignacV5
module load hdf5
```

```R

set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()

# remotes::install_github('immunogenomics/presto')


# Load RNA seurat object
multiome_WT <- readRDS("output/seurat/multiome_WT_QCV3_Option1.rds")
multiome_Bap1KO <- readRDS("output/seurat/multiome_Bap1KO_QCV3_Option1.rds")


########################
#### Add ATAC ######
########################



# WT #######################


multiome_WT_h5 <- Read10X_h5("multiome_WT/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 


# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- multiome_WT$RNA # Here we load the soupX scrublet corrected seurat
atac_counts <- multiome_WT_h5$Peaks # Here we load the h5 file

rna_names<-colnames(multiome_WT$RNA)
atac_names<-colnames(multiome_WT_h5$Peaks)

intersect <- intersect(atac_names, rna_names)

intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "multiome_WT/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
multiome_WT[["ATAC"]] <- chrom_assay




# Bap1KO #######################


multiome_Bap1KO_h5 <- Read10X_h5("multiome_Bap1KO/outs/filtered_feature_bc_matrix.h5") # the 10x hdf5 file contains both data types. 


# extract RNA and ATAC data and keep only the cells found in both assay
rna_counts <- multiome_Bap1KO$RNA # Here we load the soupX scrublet corrected seurat
atac_counts <- multiome_Bap1KO_h5$Peaks # Here we load the h5 file

rna_names<-colnames(multiome_Bap1KO$RNA)
atac_names<-colnames(multiome_Bap1KO_h5$Peaks)

intersect <- intersect(atac_names, rna_names)

intersect_atac_counts <- atac_counts[, intersect]

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(intersect_atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
intersect_atac_counts <- intersect_atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "multiome_Bap1KO/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = intersect_atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
multiome_Bap1KO[["ATAC"]] <- chrom_assay



##################
## QC ATAC #########
##################

DefaultAssay(multiome_WT) <- "ATAC" # 
DefaultAssay(multiome_Bap1KO) <- "ATAC" #

# WT ############

# QC metrics
## compute nucleosome signal score per cell
multiome_WT <- NucleosomeSignal(object = multiome_WT)
## compute TSS enrichment score per cell
multiome_WT <- TSSEnrichment(object = multiome_WT, fast = FALSE)


### QC plots
pdf("output/Signac/DensityScatter_multiome_WT_Option1.pdf", width=5, height=5)
DensityScatter(multiome_WT, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

pdf("output/Signac/VlnPlot_QC_multiome_WT_Option1.pdf", width=12, height=6)
VlnPlot(
  object = multiome_WT,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()


multiome_WT$high.tss <- ifelse(multiome_WT$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/Signac/TSSPlot_QC_multiome_WT_Option1.pdf", width=12, height=6)
TSSPlot(multiome_WT, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

multiome_WT$nucleosome_group <- ifelse(multiome_WT$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/Signac/FragmentHistogram_QC_multiome_WT_Option.pdf", width=12, height=6)
FragmentHistogram(object = multiome_WT, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).

pdf("output/Signac/VlnPlot_QC_multiome_WT_Option_nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = multiome_WT,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 50000)
dev.off()

## subset cells that pass QC
multiome_WT_QCV2 <- subset(
  x = multiome_WT,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 50000 &
    nucleosome_signal < 1.25 &
    TSS.enrichment > 2
)

multiome_WT_QCV2




# Bap1KO ############

# QC metrics
## compute nucleosome signal score per cell
multiome_Bap1KO <- NucleosomeSignal(object = multiome_Bap1KO)
## compute TSS enrichment score per cell
multiome_Bap1KO <- TSSEnrichment(object = multiome_Bap1KO, fast = FALSE)
## add blacklist ratio and fraction of reads in peaks


### QC plots
pdf("output/Signac/DensityScatter_multiome_Bap1KO_Option1.pdf", width=5, height=5)
DensityScatter(multiome_Bap1KO, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

pdf("output/Signac/VlnPlot_QC_multiome_Bap1KO_Option1.pdf", width=12, height=6)
VlnPlot(
  object = multiome_Bap1KO,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5 )
dev.off()


multiome_Bap1KO$high.tss <- ifelse(multiome_Bap1KO$TSS.enrichment > 2, 'High', 'Low') # 3 is coonly used but could be adjusted; based on violin and below plot
pdf("output/Signac/TSSPlot_QC_multiome_Bap1KO_Option1.pdf", width=12, height=6)
TSSPlot(multiome_Bap1KO, group.by = 'high.tss') + NoLegend()
dev.off()
#--> Here aim to have a clean sharp peak for High, flatter one for Low TSS.enrichment

multiome_Bap1KO$nucleosome_group <- ifelse(multiome_Bap1KO$nucleosome_signal > 1.25, 'NS > 1.25', 'NS < 1.25') # Adjust value based on nucleosome_signal VlnPlot
pdf("output/Signac/FragmentHistogram_QC_multiome_Bap1KO_Option.pdf", width=12, height=6)
FragmentHistogram(object = multiome_Bap1KO, group.by = 'nucleosome_group', region = "chr1-1-20000000") # region = "chr1-1-20000000" this need to be added to avoid bug discuss here: https://github.com/stuart-lab/signac/issues/199
dev.off()
#--> Here aim to have strong nucleosome-free region peak (around 50 bp), mono-nucleosome peak (around 200 bp), and sometimes di-nucleosome peak (around 400 bp).

pdf("output/Signac/VlnPlot_QC_multiome_Bap1KO_Option_nCount_ATAC.pdf", width=12, height=6)
VlnPlot(
  object = multiome_Bap1KO,
  features = c('nCount_ATAC'),
  pt.size = 0.1,
  ncol = 5 ) +
  ylim(0, 50000)
dev.off()

## subset cells that pass QC
multiome_Bap1KO_QCV2 <- subset(
  x = multiome_Bap1KO,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 50000 &
    nucleosome_signal < 1.25 &
    TSS.enrichment > 2
)

multiome_Bap1KO_QCV2





###########################################################################
# saveRDS(multiome_WT_QCV1, file = "output/seurat/multiome_WT_QCV2.rds") 
# saveRDS(multiome_Bap1KO_QCV1, file = "output/seurat/multiome_Bap1KO_QCV2.rds") 
###########################################################################
multiome_WT_QCV2 <- readRDS(file = "output/seurat/multiome_WT_QCV2.rds")
multiome_Bap1KO_QCV2 <- readRDS(file = "output/seurat/multiome_Bap1KO_QCV2.rds")


# Integrate WT and Bap1KO

######################################################################################################
## Pre-processing RNA  - QCV2 parameters ####################################################################
######################################################################################################

DefaultAssay(multiome_WT_QCV2) <- "RNA"
DefaultAssay(multiome_Bap1KO_QCV2) <- "RNA"



multiome_WT_QCV2 <- SCTransform(multiome_WT_QCV2, method = "glmGamPoi", ncells = 5949, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 2000)
multiome_Bap1KO_QCV2 <- SCTransform(multiome_Bap1KO_QCV2, method = "glmGamPoi", ncells = 6517, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 2000)


# Data integration (check active assay is 'SCT')
srat.list <- list(multiome_WT_QCV2 = multiome_WT_QCV2, multiome_Bap1KO_QCV2 = multiome_Bap1KO_QCV2)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 2000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

srat.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
multiome_WT_Bap1KO_QCV2.sct <- IntegrateData(anchorset = srat.anchors, normalization.method = "SCT")

set.seed(42)

###########################################################################
# saveRDS(multiome_WT_Bap1KO_QCV2.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2.sct.rds") 
###########################################################################
multiome_WT_Bap1KO_QCV2.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2.sct.rds")


DefaultAssay(multiome_WT_Bap1KO_QCV2.sct) <- "integrated"

multiome_WT_Bap1KO_QCV2.sct <- RunPCA(multiome_WT_Bap1KO_QCV2.sct, verbose = FALSE, npcs = 40)
multiome_WT_Bap1KO_QCV2.sct <- RunUMAP(multiome_WT_Bap1KO_QCV2.sct, reduction = "pca", dims = 1:40, verbose = FALSE)
multiome_WT_Bap1KO_QCV2.sct <- FindNeighbors(multiome_WT_Bap1KO_QCV2.sct, reduction = "pca", k.param = 35, dims = 1:40)
multiome_WT_Bap1KO_QCV2.sct <- FindClusters(multiome_WT_Bap1KO_QCV2.sct, resolution = 0.5, verbose = FALSE, algorithm = 4) # 

multiome_WT_Bap1KO_QCV2.sct$orig.ident <- factor(multiome_WT_Bap1KO_QCV2.sct$orig.ident, levels = c("multiome_WT", "multiome_Bap1KO")) # Reorder untreated 1st

pdf("output/Signac/UMAP_multiome_WT_Bap1KO-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression-numeric_V1.pdf", width=6, height=6)
DimPlot(multiome_WT_Bap1KO_QCV2.sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 6)
dev.off()

DefaultAssay(multiome_WT_Bap1KO_QCV2.sct) <- "SCT"


pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersList4-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=30)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c(  "Pax6" ,  "Eomes",  "Prox1", "Neurod1", "Sema5a",  "Tac2", "Hs3st1", "Nrn1",  "Pantr1", "Igfbpl1", "Frmd4b",  "Satb2", "Itpr1",  "Nts", "Nr4a2", "Lmo3", "B3gat1",  "Cck", "Insm1",  "Crym", "Snca", "Nrp2",  "Gad1", "Grin2d", "Calb1", "Npy", "Gria3", "Lhx6",  "Lhx1",  "Pdgfra", "Olig1",  "Csf1r", "Gpr34", "Gpr183", "Cx3cr1", "Aldh1a2", "Vtn", "Foxc1", "Id1", "Hes1", "Mki67", "Pcna", "Vim"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()








######################################################################################################
######################################################################################################
######################################################################################################
## Pre-processing RNA  - QCV2 parameters (after 20240830 Casey (C) meeting) = QCV2vC1 ##############################
######################################################################################################

DefaultAssay(multiome_WT_QCV2) <- "RNA"
DefaultAssay(multiome_Bap1KO_QCV2) <- "RNA"



multiome_WT_QCV2 <- SCTransform(multiome_WT_QCV2, method = "glmGamPoi", ncells = 5949, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 2000)
multiome_Bap1KO_QCV2 <- SCTransform(multiome_Bap1KO_QCV2, method = "glmGamPoi", ncells = 6517, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 2000)


# Data integration (check active assay is 'SCT')
srat.list <- list(multiome_WT_QCV2 = multiome_WT_QCV2, multiome_Bap1KO_QCV2 = multiome_Bap1KO_QCV2)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 2000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

srat.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
multiome_WT_Bap1KO_QCV2vC1.sct <- IntegrateData(anchorset = srat.anchors, normalization.method = "SCT")

set.seed(42)

###########################################################################
# saveRDS(multiome_WT_Bap1KO_QCV2.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1.sct.rds")
# multiome_WT_Bap1KO_QCV2.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2.sct.rds") 
###########################################################################


multiome_WT_Bap1KO_QCV2vC1.sct = multiome_WT_Bap1KO_QCV2.sct
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1.sct.rds")


DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "integrated"

multiome_WT_Bap1KO_QCV2vC1.sct <- RunPCA(multiome_WT_Bap1KO_QCV2vC1.sct, verbose = FALSE, npcs = 40)
multiome_WT_Bap1KO_QCV2vC1.sct <- RunUMAP(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "pca", dims = 1:40, verbose = FALSE)
multiome_WT_Bap1KO_QCV2vC1.sct <- FindNeighbors(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "pca", k.param = 42, dims = 1:40)
multiome_WT_Bap1KO_QCV2vC1.sct <- FindClusters(multiome_WT_Bap1KO_QCV2vC1.sct, resolution = 0.65, verbose = FALSE, algorithm = 4) # 

multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident, levels = c("multiome_WT", "multiome_Bap1KO")) # Reorder untreated 1st

pdf("output/Signac/UMAP_multiome_WT_Bap1KO-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-numeric_V1.pdf", width=6, height=6)

pdf("output/Signac/UMAP_test.pdf", width=6, height=6)
DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 6)
dev.off()


DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "SCT"


pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersList4-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression
.pdf", width=15, height=30)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = c(  "Pax6" ,  "Eomes",  "Prox1", "Neurod1", "Sema5a",  "Tac2", "Hs3st1", "Nrn1",  "Pantr1", "Igfbpl1", "Frmd4b",  "Satb2", "Itpr1",  "Nts", "Nr4a2", "Lmo3", "B3gat1",  "Cck", "Insm1",  "Crym", "Snca", "Nrp2",  "Gad1", "Grin2d", "Calb1", "Npy", "Gria3", "Lhx6",  "Lhx1",  "Pdgfra", "Olig1",  "Csf1r", "Gpr34", "Gpr183", "Cx3cr1", "Aldh1a2", "Vtn", "Foxc1", "Id1", "Hes1", "Mki67", "Pcna", "Vim"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-NSCtypesmarkers-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression
.pdf", width=10, height=10)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = c( "Id1", "Hes1" , "Mki67", "Pcna", "Vim" ), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

######################################################################################################
######################################################################################################
######################################################################################################



pdf("output/Signac/VlnPlot_QCmetrics_RNA_WT-QCV3_dim40kparam50res07algo4feat2000_noCellCycleRegression.pdf", width=20, height=5)
VlnPlot(multiome_WT_Bap1KO_QCV2.sct,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()


#### QC metrics investigation ####################################
pdf("output/Signac/VlnPlot_QCmetrics_RNA_WT_nFeature_RNA-QCV2_dim28kparam15res02algo4feat2000_noCellCycleRegression.pdf", width=5, height=5)
VlnPlot(multiome_WT_Bap1KO_QCV1.sct,features = c("nFeature_RNA")) +
  ylim(0,2000)
dev.off()
pdf("output/Signac/VlnPlot_QCmetrics_RNA_WT_nCount_RNA-QCV2_dim28kparam15res02algo4feat2000_noCellCycleRegression.pdf", width=5, height=5)
VlnPlot(multiome_WT_Bap1KO_QCV1.sct,features = c("nCount_RNA")) +
  ylim(0,5000)
dev.off()
############################################################


pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000_numeric_V1.pdf", width=12, height=6)
DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap", split.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 6)
dev.off()


pdf("output/Signac/UMAP_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000_noSplit_numeric_V1.pdf", width=5, height=5)
DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 6)
dev.off()

pdf("output/Signac/UMAP_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000_numeric_overlap_V1.pdf", width=6, height=5)
DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap", group.by = "orig.ident", pt.size = 0.000001, cols = c("blue","red"))
dev.off()





## Downsampling with bootstrap to compare the nb of cell per cell types

library("tidyverse")

### Identify the unique clusters
unique_clusters <- unique(Idents(multiome_WT_Bap1KO_QCV2vC1.sct))

### Create empty matrices to store cell counts
control_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
Bap1KO_clusters_counts <- matrix(0, nrow=100, ncol=length(unique_clusters))
colnames(control_clusters_counts) <- unique_clusters
colnames(Bap1KO_clusters_counts) <- unique_clusters

### Loop through 100 iterations
multiome_WT_Bap1KO_QCV2vC1.sct_WT <- which(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident == 'multiome_WT')
multiome_WT_Bap1KO_QCV2vC1.sct_Bap1KO <- which(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident == 'multiome_Bap1KO')

for (i in 1:100) { # Change this to 100 for the final run
  # Downsampling
  multiome_WT_Bap1KO_QCV2vC1.sct_Bap1KO_downsample <- sample(multiome_WT_Bap1KO_QCV2vC1.sct_Bap1KO, 5949)
  multiome_WT_Bap1KO_QCV2vC1.sct_integrated_downsample <- multiome_WT_Bap1KO_QCV2vC1.sct[,c(multiome_WT_Bap1KO_QCV2vC1.sct_Bap1KO_downsample, multiome_WT_Bap1KO_QCV2vC1.sct_WT)]

  # Count nb of cells in each cluster
  control_clusters <- table(Idents(multiome_WT_Bap1KO_QCV2vC1.sct_integrated_downsample)[multiome_WT_Bap1KO_QCV2vC1.sct_integrated_downsample$orig.ident == "multiome_WT"])
  Bap1KO_clusters <- table(Idents(multiome_WT_Bap1KO_QCV2vC1.sct_integrated_downsample)[multiome_WT_Bap1KO_QCV2vC1.sct_integrated_downsample$orig.ident == "multiome_Bap1KO"])

  # Align the counts with the unique clusters
  control_clusters_counts[i, names(control_clusters)] <- as.numeric(control_clusters)
  Bap1KO_clusters_counts[i, names(Bap1KO_clusters)] <- as.numeric(Bap1KO_clusters)
}


### Calculate mean and standard error
mean_control_clusters <- colMeans(control_clusters_counts)
mean_Bap1KO_clusters <- colMeans(Bap1KO_clusters_counts)
std_error_WT_clusters <- apply(control_clusters_counts, 2, sd) / sqrt(100)

# Chi-squared test
p_values <- numeric(length(unique_clusters))

for (i in 1:length(unique_clusters)) {
  # Create a matrix to store the counts for the chi-squared test
  contingency_table <- matrix(0, nrow=2, ncol=2)
  colnames(contingency_table) <- c("WT", "Bap1KO")
  rownames(contingency_table) <- c("Cluster", "NotCluster")
  
  for (j in 1:100) { # Number of bootstrap iterations
    contingency_table[1,1] <- control_clusters_counts[j,i]
    contingency_table[1,2] <- Bap1KO_clusters_counts[j,i]
    contingency_table[2,1] <- sum(control_clusters_counts[j,-i])
    contingency_table[2,2] <- sum(Bap1KO_clusters_counts[j,-i])
    
    # Perform the chi-squared test on the contingency table
    chi_test <- chisq.test(contingency_table)
    
    # Store the p-value
    p_values[i] <- p_values[i] + chi_test$p.value
  }
  
  # Average the p-values across all bootstrap iterations
  p_values[i] <- p_values[i] / 100
}

# Adjust the p-values
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Create a tidy data frame for plotting
plot_data <- data.frame(
  cluster = names(mean_control_clusters),
  untreated = mean_control_clusters,
  dasatinib = mean_Bap1KO_clusters,
  std_error_WT = std_error_WT_clusters,
  p_value = adjusted_p_values
) %>%
  gather(key = "condition", value = "value", -cluster, -std_error_WT, -p_value) %>%
  mutate(
    condition = if_else(condition == "untreated", "WT", "Bap1KO"),
    significance = ifelse(p_value < 0.0001, "***",
                       ifelse(p_value < 0.001, "**",
                              ifelse(p_value < 0.05, "*", "")))
  )

plot_data$condition <- factor(plot_data$condition, levels = c("WT", "Bap1KO")) # Reorder untreated 1st
plot_data$cluster <- factor(plot_data$cluster, levels = c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")) 



# Plotting using ggplot2
pdf("output/Signac/Cluster_cell_counts_BootstrapDownsampling100_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000_numeric.pdf", width=9, height=4)
ggplot(plot_data, aes(x = cluster, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    data = dplyr::filter(plot_data, condition == "Bap1KO"),
    aes(label = significance, y = value + std_error_WT_clusters),
    vjust = -0.8,
    position = position_dodge(0.9), size = 5
  ) +
  scale_fill_manual(values = c("WT" = "#4365AE", "Bap1KO" = "#981E33")) +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13)) +
  theme(axis.text.y = element_text(size = 13)) +
  ylim(0,900)
dev.off()






# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA"

# PART ADDED as *_correct ###########
multiome_WT_Bap1KO_QCV2vC1.sct <- NormalizeData(multiome_WT_Bap1KO_QCV2vC1.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(multiome_WT_Bap1KO_QCV2vC1.sct)
multiome_WT_Bap1KO_QCV2vC1.sct <- ScaleData(multiome_WT_Bap1KO_QCV2vC1.sct, features = all.genes) # zero-centres and scales it
#########


## DEGs keeping ALL genes
multiome_WT_Bap1KO_QCV2vC1.sct$celltype.stim <- paste(multiome_WT_Bap1KO_QCV2vC1.sct$seurat_clusters, multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident,
    sep = "-")
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "celltype.stim"

cluster1 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "1-multiome_Bap1KO", ident.2 = "1-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
cluster2 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "2-multiome_Bap1KO", ident.2 = "2-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
cluster3 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "3-multiome_Bap1KO", ident.2 = "3-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
cluster4 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "4-multiome_Bap1KO", ident.2 = "4-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")    
cluster5 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "5-multiome_Bap1KO", ident.2 = "5-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
cluster6 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "6-multiome_Bap1KO", ident.2 = "6-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster7 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "7-multiome_Bap1KO", ident.2 = "7-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster8 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "8-multiome_Bap1KO", ident.2 = "8-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster9 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "9-multiome_Bap1KO", ident.2 = "9-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster10 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "10-multiome_Bap1KO", ident.2 = "10-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster11 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "11-multiome_Bap1KO", ident.2 = "11-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster12 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "12-multiome_Bap1KO", ident.2 = "12-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA", min.cells.group = 2)   # min.cells.group = 2 added to avoid `error in validatecellgroups cell group 2 has fewer than 3 cells
cluster13 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "13-multiome_Bap1KO", ident.2 = "13-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster14 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "14-multiome_Bap1KO", ident.2 = "14-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster15 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "15-multiome_Bap1KO", ident.2 = "15-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster16 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "16-multiome_Bap1KO", ident.2 = "16-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster17 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "17-multiome_Bap1KO", ident.2 = "17-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster18 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "18-multiome_Bap1KO", ident.2 = "18-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
cluster19 <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "19-multiome_Bap1KO", ident.2 = "19-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA")  
    

### save output
write.table(cluster1, file = "output/Signac/cluster1-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster2, file = "output/Signac/cluster2-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster3, file = "output/Signac/cluster3-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster4, file = "output/Signac/cluster4-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster5, file = "output/Signac/cluster5-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster6, file = "output/Signac/cluster6-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster7, file = "output/Signac/cluster7-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster8, file = "output/Signac/cluster8-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster9, file = "output/Signac/cluster9-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster10, file = "output/Signac/cluster10-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster11, file = "output/Signac/cluster11-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster12, file = "output/Signac/cluster12-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster13, file = "output/Signac/cluster13-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster14, file = "output/Signac/cluster14-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster15, file = "output/Signac/cluster15-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster16, file = "output/Signac/cluster16-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster17, file = "output/Signac/cluster17-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster18, file = "output/Signac/cluster18-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cluster19, file = "output/Signac/cluster19-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#--> I added `*_correct` = the one I perfrom lognorm and scale prior DEG...



#### import all clsuter DEGs output :
cluster_types <- c("cluster1", "cluster2", "cluster3", 
                   "cluster4", "cluster5", "cluster6", 
                   "cluster7", "cluster8", 
                   "cluster9", "cluster10", "cluster11", 
                   "cluster12", "cluster13", "cluster14", "cluster15", 
                   "cluster16", "cluster17", "cluster18", "cluster19")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/Signac/", cluster, "-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}



# DEGs nb dotplot
## Initialize an empty data frame to store the summary
DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())

## List of cell types
cell_types <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19")

## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/Signac/", cell_type, "-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "")
  # Check if file exists
  if (!file.exists(file_name)) {
    print(paste("File not found:", file_name))
    next
  }
  # Read the DEGs data
  deg_data <- read.table(file_name, header = TRUE, sep = "\t")
  # Count the number of significant DEGs, defaulting to 0 if none are found
  num_degs <- sum(deg_data$p_val_adj < 0.05, na.rm = TRUE)
  # Ensure num_degs is not NA
  if (is.na(num_degs)) {
    num_degs <- 0
  }
  # Append to the summary table
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))
}


DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19")) 

DEG_count$Cluster_Number <- as.numeric(sub("cluster", "", DEG_count$Cell_Type))

DEG_count$Cluster_Number <- factor(DEG_count$Cluster_Number, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")) 


## Dotplot

cell_type_colors <- c(
  "cluster1" = "#FFB6C1", "cluster2" = "#FFA07A", "cluster3" = "#FFD700", "cluster4" = "#ADFF2F", 
  "cluster5" = "#7FFF00", "cluster6" = "#32CD32", "cluster7" = "#3CB371", "cluster8" = "#00FA9A", 
  "cluster9" = "#00CED1", "cluster10" = "#4682B4", "cluster11" = "#1E90FF", "cluster12" = "#6495ED", 
  "cluster13" = "#4169E1", "cluster14" = "#BA55D3", "cluster15" = "#DA70D6", "cluster16" = "#EE82EE", 
  "cluster17" = "#FF69B4", "cluster18" = "#FF1493", "cluster19" = "#DB7093"
)

# Generate the dot plot
pdf("output/Signac/Dotplot_Bap1KO_DEG_count_multiome_QCV2vC1_dim40kparam42res065algo4feat2000.pdf", width=9, height=4)
ggplot(DEG_count, aes(x = Cluster_Number, y = 1) )+
  geom_point(aes(size = ifelse(Num_DEGs == 0, 1, Num_DEGs), fill = Cell_Type) , shape = 21, color = "black") +
  scale_size_continuous(range = c(1, 15)) +
  scale_fill_manual(values = cell_type_colors, guide = "none") +  # Remove the legend for cell types
  theme_void() +
  labs(title = "Number of DEGs per Cell Type", x = "Cell Type", y = "", size = "Number of DEGs") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 50, size = 15),  # Adjust the position of x-axis labels
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(),
    legend.position = "right"
  ) +
  guides(size = guide_legend(title = "Number of DEGs", title.position = "top", title.hjust = 0.5))  # Keep only the size legend
dev.off()






## DEGs keeping ALL genes for NSC subcluster
#### in WT NSC vs NSC
multiome_WT_Bap1KO_QCV2vC1.sct$celltype.stim <- paste(multiome_WT_Bap1KO_QCV2vC1.sct$cluster.annot, multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident,
    sep = "-")
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "celltype.stim"

NSC_proliferative_2__NSC_proliferative_1_WTonly <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "NSC_proliferative_2-multiome_WT", ident.2 = "NSC_proliferative_1-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
NSC_proliferative_2__NSC_quiescent_WTonly <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "NSC_proliferative_2-multiome_WT", ident.2 = "NSC_quiescent-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
NSC_proliferative_1__NSC_quiescent_WTonly <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "NSC_proliferative_1-multiome_WT", ident.2 = "NSC_quiescent-multiome_WT",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
write.table(NSC_proliferative_2__NSC_proliferative_1_WTonly, file = "output/Signac/srat-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-NSC_proliferative_2__NSC_proliferative_1_WTonly.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(NSC_proliferative_2__NSC_quiescent_WTonly, file = "output/Signac/srat-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-NSC_proliferative_2__NSC_quiescent_WTonly.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(NSC_proliferative_1__NSC_quiescent_WTonly, file = "output/Signac/srat-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-NSC_proliferative_1__NSC_quiescent_WTonly.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#### in W and Bap1KO NSC vs NSC

Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "cluster.annot"

NSC_proliferative_2__NSC_proliferative_1_WTBap1KO <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "NSC_proliferative_2", ident.2 = "NSC_proliferative_1",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
NSC_proliferative_2__NSC_quiescent_WTBap1KO <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "NSC_proliferative_2", ident.2 = "NSC_quiescent",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
NSC_proliferative_1__NSC_quiescent_WTBap1KO <- FindMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, ident.1 = "NSC_proliferative_1", ident.2 = "NSC_quiescent",
    verbose = TRUE,
    test.use = "wilcox",
    logfc.threshold = -Inf,
    min.pct = -Inf,
    min.diff.pct = -Inf, # 
    assay = "RNA") 
write.table(NSC_proliferative_2__NSC_proliferative_1_WTBap1KO, file = "output/Signac/srat-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-NSC_proliferative_2__NSC_proliferative_1_WTBap1KO.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(NSC_proliferative_2__NSC_quiescent_WTBap1KO, file = "output/Signac/srat-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-NSC_proliferative_2__NSC_quiescent_WTBap1KO.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(NSC_proliferative_1__NSC_quiescent_WTBap1KO, file = "output/Signac/srat-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-NSC_proliferative_1__NSC_quiescent_WTBap1KO.txt", sep = "\t", quote = FALSE, row.names = TRUE)










### Find all markers 
##### With FindAllMarkers()
all_markers <- FindAllMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.table(all_markers, file = "output/Signac/srat_multiome_WT_Bap1KO_QCV3_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(all_markers, file = "output/Signac/srat_multiome_WT_Bap1KO-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(all_markers, file = "output/Signac/srat_multiome_WT_Bap1KO-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-all_markers_clusterName.txt", sep = "\t", quote = FALSE, row.names = TRUE)







# Display the CONSERVED marker genes of each cluster

## NEW version using multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric_label.rds")  (correct1)
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "cluster.annot"

## DEGs cluster versus all other
cluster1.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "1", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster1")
cluster2.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "2", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster2")
cluster3.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "3", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster3")
cluster4.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "4", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster4")
cluster5.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "5", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster5")
cluster6.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "6", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster6")
cluster7.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "7", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster7")
cluster8.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "8", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster8")
cluster9.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "9", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster9")
cluster10.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "10", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster10")
cluster11.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "11", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster11")
cluster12.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "12", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster12")
cluster13.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "13", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster13")
cluster14.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "14", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster14")
cluster15.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "15", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster15")
cluster16.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "16", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster16")
cluster17.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "17", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster17")
cluster18.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "18", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster18")
cluster19.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "19", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster19")
cluster20.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", ident.1 = "20", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster20")

## Combine all conserved markers into one data frame
all_conserved <- bind_rows(cluster1.conserved,cluster2.conserved,cluster3.conserved,cluster4.conserved,cluster5.conserved,cluster6.conserved,cluster7.conserved,cluster8.conserved,cluster9.conserved,cluster10.conserved,cluster11.conserved,cluster12.conserved,cluster13.conserved,cluster14.conserved,cluster15.conserved,cluster16.conserved,cluster17.conserved,cluster18.conserved,cluster19.conserved,cluster20.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/Signac/srat_all_conserved_markers_multiome_WT_Bap1KO_QCV3.txt", sep = "\t", quote = FALSE, row.names = TRUE)
##########################


## PAST version
Idents(multiome_WT_Bap1KO_QCV2.sct) <- "seurat_clusters"

## DEGs cluster versus all other
cluster1.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "1", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster1")
cluster2.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "2", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster2")
cluster3.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "3", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster3")
cluster4.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "4", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster4")
cluster5.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "5", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster5")
cluster6.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "6", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster6")
cluster7.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "7", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster7")
cluster8.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "8", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster8")
cluster9.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "9", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster9")
cluster10.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "10", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster10")
cluster11.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "11", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster11")
cluster12.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "12", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster12")
cluster13.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "13", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster13")
cluster14.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "14", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster14")
cluster15.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "15", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster15")
cluster16.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "16", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster16")
cluster17.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "17", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster17")
cluster18.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "18", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster18")
cluster19.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "19", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster19")
cluster20.conserved <- FindConservedMarkers(multiome_WT_Bap1KO_QCV2.sct, assay = "RNA", ident.1 = "20", grouping.var = "orig.ident", verbose = TRUE) %>% mutate(cluster = "cluster20")

## Combine all conserved markers into one data frame
all_conserved <- bind_rows(cluster1.conserved,cluster2.conserved,cluster3.conserved,cluster4.conserved,cluster5.conserved,cluster6.conserved,cluster7.conserved,cluster8.conserved,cluster9.conserved,cluster10.conserved,cluster11.conserved,cluster12.conserved,cluster13.conserved,cluster14.conserved,cluster15.conserved,cluster16.conserved,cluster17.conserved,cluster18.conserved,cluster19.conserved,cluster20.conserved)

all_conserved$gene <- rownames(all_conserved)
## Write all conserved markers to a file
write.table(all_conserved, file = "output/Signac/srat_all_conserved_markers_multiome_WT_Bap1KO_QCV3.txt", sep = "\t", quote = FALSE, row.names = TRUE)

all_conserved <- read_delim("output/Signac/srat_all_conserved_markers_multiome_WT_Bap1KO_QCV3.txt", delim = "\t", col_names = TRUE)
colnames(all_conserved) <- c("gene", colnames(all_conserved)[-ncol(all_conserved)]) # shift all column name to the right

## Find the top 5 conserved markers for each cluster
top10_conserved <- all_conserved %>%
  mutate(cluster = factor(cluster, levels = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19", "cluster20"))) %>% 
  separate(gene, into = c("gene", "suffix"), sep = "\\.\\.\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  group_by(cluster) %>% 
  arrange((max_pval)) %>% 
  slice_head(n = 5) %>% 
  ungroup() %>% 
  arrange(match(cluster, c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19", "cluster20")))


## Write the top 10 conserved markers for each cluster to a file

## Visualize the top 10/3 conserved markers for each cluster
marker_genes_conserved <- unique(top10_conserved$gene)
levels(multiome_WT_Bap1KO_QCV2.sct) <- c("1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7",
  "8",
  "9",
  "10",
  "11",
  "12",
  "13",
  "14",
  "15",
  "16",
  "17",
  "18",
  "19",
  "20")

DefaultAssay(multiome_WT_Bap1KO_QCV2.sct) <- "SCT"


pdf("output/Signac/DotPlot_SCT_top5_conserved_multiome_WT_Bap1KO_QCV3.pdf", width=19, height=5)
DotPlot(multiome_WT_Bap1KO_QCV2.sct, features = marker_genes_conserved, cols = c("grey", "red")) + RotatedAxis()
dev.off()


# SAVE #########################################################################################
## saveRDS(multiome_WT_Bap1KO_QCV2.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV3.sct_numeric.rds") 
## saveRDS(multiome_WT_Bap1KO_QCV2vC1.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric.rds") 
# I used multiome_WT_Bap1KO_QCV2 but was already multiome_WT_Bap1KO_QCV3... So here I update the file name
## saveRDS(multiome_WT_Bap1KO_QCV2vC1.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_correct1.sct_numeric.rds")  # THIS ONE IS GOOD FOR THE DEG! logNorm/ScaleData
################################################################################################
multiome_WT_Bap1KO_QCV3.sct = multiome_WT_Bap1KO_QCV2.sct

multiome_WT_Bap1KO_QCV3.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV3.sct_numeric.rds")





############################ EasyCellType automatic annotation ##########################################
# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/Signac/srat_multiome_WT_Bap1KO_QCV3_all_markers.txt", header = TRUE, row.names = 1)
### Filter either WT or cYAPKO
all_markers <- all_markers[grepl("multiome_WT$", all_markers$cluster), ]
all_markers <- all_markers[grepl("multiome_Bap1KO$", all_markers$cluster), ]

## Convert geneSymbol to EntrezID
all_markers$entrezid <- mapIds(org.Mm.eg.db,
                           keys=all_markers$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
all_markers <- na.omit(all_markers)

## Sort the datafram (data frame containing Entrez IDs, clusters and expression scores)

all_markers_sort <- data.frame(gene=all_markers$entrezid, cluster=all_markers$cluster, 
                      score=all_markers$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(all_markers_sort[, 1:3])

## Run the enrihcment analysis
annot.GSEA <- easyct(input.d, db="panglao", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=NULL, p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain", "Cerebellum", "Hippocampus"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?



annot.GSEA <- easyct(input.d, db="clustermole", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?

## plots


pdf("output/Signac/EasyCellType_dotplot_SCT_WT-cellmarker_brainCerebellumHippocampus.pdf", width=6, height=8)
pdf("output/Signac/EasyCellType_dotplot_SCT_WT-clustermole_brain.pdf", width=6, height=8)
plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


## check some genes


pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-EpendymalCell-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c(  "Rabl2", "Cfap54", "Ccdc153", "Foxj1", "Pifo", "Dynlrb2", "Rsph1", "Cfap44"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-RadialGliaCell-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c(  "Pax6", "Slc1a3", "Pdgfd", "Gli3", "Notch3", "Vcam1", "Hes5", "Olig2", "Gfap", "Emx2", "Cdh4", "Spry1", "Axin2", "Riiad1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-RadialGliaCell2-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c(  ), max.cutoff = 1, cols = c("grey", "red"))
dev.off()



pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-Astrocytes-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c( "Gfap", "Slc1a2", "Acsl6", "Agt", "Aqp4", "Apoe", "S100b", "Sox9", "Gsta4", "Srr", "Aldh1l1", "Slc39a12"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-ChoroidPlexusCells-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c( "Ttr", "Kl", "Clic6", "Prlr", "Chmp1a", "Slc26a11", "Slc23a2", "Wfikkn2", "Slc2a12", "Cldn1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-BasketCellChatGPT-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c( "Pvalb", "Gad1", "Gad2", "Cnr1", "Nkx2-1", "Syt2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-ChandelierCell-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c( "Ntf3" , "Sntb1", "Unc5b", "Rasgrp1", "Prkg1", "Vipr2" ), cols = c("grey", "red"))
dev.off()

pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-ChandelierCellChatGPT-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=10, height=10)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c( "Gad1", "Gad2","Pvalb", "Slc12a5" , "Ank3", "Syt2" ), cols = c("grey", "red"))
dev.off()



pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-TopMarkerGenes-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV3.sct, features = c( "Pax6", "Eomes", "Sema5a", "Hs3st1", "Igfbpl1", "Satb2", "Nts", "Cck", "Snca", "Gad1", "Lhx1", "Pdgfra", "Csf1r", "Foxj1", "Aqp4", "Pdgfd", "Slc12a5", "Aldh1a2"), cols = c("grey", "red"), max.cutoff = 1)
dev.off()


##########################################################################################


############ V2 naming

Cluster1 = PyNs_SubC_CA1 (subiculum PyNs)
Cluster2 = PyNs_SubC_CA23 (subiculum PyNs)
Cluster3 = IN_1 (interneuron)
Cluster4 = PyNs_RSC_UL (Retrosplenial Cortical Pyramidal neurons, upper layer)
Cluster5 = SubC_1 (subiculum)
Cluster6 = DG_GC (Dentate Gyrus granule cells)
Cluster7 = Astrocyte
Cluster8 = SubC_2 (subiculum)
Cluster9 = IN_2 (interneuron)
Cluster10 = NSC (Neural Stem Cells)
Cluster11 = IP (Intermediate Progenitors)
Cluster12 = PyNs_RSC_ML (Retrosplenial Cortical Pyramidal neurons, middle layer)
Cluster13 = CR (Cajal Retzius)
Cluster14 = PyNs_RSC_DL (Retrosplenial Cortical Pyramidal neurons, deep layer)
Cluster15 = OPC (Oligodendrocyte progenitor cells)
Cluster16 = Chandelier_Cells
Cluster17 = Meningeal_Cells
Cluster18 = Microglia
Cluster19 = Radial_Glia_Cells
Cluster20 = Ependymal_Cells



new.cluster.ids <- c(
  "PyNs_SubC_CA1" ,
  "PyNs_SubC_CA23" ,
  "IN_1",
  "PyNs_RSC_UL",
  "SubC_1",
  "DG_GC",
  "Astrocyte",
  "SubC_2",
  "IN_2",
  "NSC",
  "IP",
  "PyNs_RSC_ML",
  "CR",
  "PyNs_RSC_DL",
  "OPC",
  "Chandelier_Cells",
  "Meningeal_Cells",
  "Microglia",
  "Radial_Glia_Cells",
  "Ependymal_Cells"
)

names(new.cluster.ids) <- levels(multiome_WT_Bap1KO_QCV3.sct)
multiome_WT_Bap1KO_QCV3.sct <- RenameIdents(multiome_WT_Bap1KO_QCV3.sct, new.cluster.ids)

multiome_WT_Bap1KO_QCV3.sct$cluster.annot <- Idents(multiome_WT_Bap1KO_QCV3.sct) # create a new slot in my seurat object


pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV3_label.pdf", width=12, height=6)
DimPlot(multiome_WT_Bap1KO_QCV3.sct, reduction = "umap", split.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()


pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV3_noSplit_label.pdf", width=7, height=5)
DimPlot(multiome_WT_Bap1KO_QCV3.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 4)
dev.off()



# All in dotplot
DefaultAssay(multiome_WT_Bap1KO_QCV3.sct) <- "SCT"

Neural Stem Cells (NSC) = Pax6 (should have more cell types)
Intermediate Progenitors (IP) = Eomes
Dentate Gyrus Granucle Cells (DG) = Prox1, Neurod1, Sema5a
CA1 = Cck, Insm1
CA3 = Crym, Snca, Nrp2
Pyramidal neurons deep layer (DL) = Tac2, Hs3st1, Nrn1
Pyramidal neurons middle layer (ML) = Pantr1, Igfbpl1, Frmd4b
Pyramidal neurons upper layer (UL) = Satb2, Itpr1
Interneurons (IN) = Gad1, Grin2d, Reln, Calb1, Npy, Gria3, Lhx6
Cajal Retzius (CR) = Lhx1
Subiculum (SubC) = Nts, Nr4a2, Lmo3, B3gat1
Microglia = Csf1r, Gpr34, Gpr183, Cx3cr1
OPC = Pdgfra, Olig1
Ependymal_Cells = Foxj1, Cfap44, Dynlrb2
Astrocyte = Aqp4, Gfap
Radial Glia Cells = Gli3, Pdgfd
Chandellier cell = Ntf3, Prkg1, Slc12a5
Meningeal cells = Aldh1a2, Vtn, Lum, Foxc1, Igf2

all_markers <- c(
"Pax6",
"Eomes",
"Gli3", "Pdgfd",
"Pdgfra", "Olig1",
"Aqp4", "Gfap",
"Foxj1", "Cfap44", "Dynlrb2",
"Csf1r", "Gpr34", "Gpr183", "Cx3cr1",
"Aldh1a2", "Vtn", "Lum", "Foxc1", "Igf2",
"Lhx1",
"Prox1", "Neurod1", "Sema5a",
"Cck", "Insm1",
"Crym", "Snca", "Nrp2",
"Satb2", "Itpr1",
"Pantr1", "Igfbpl1", "Frmd4b",
"Tac2", "Hs3st1", "Nrn1",
"Nts", "Nr4a2", "Lmo3", "B3gat1",
"Gad1", "Grin2d", "Reln", "Calb1", "Npy", "Gria3", "Lhx6",
 "Ntf3", "Prkg1", "Slc12a5"
)



levels(multiome_WT_Bap1KO_QCV3.sct) <- c(
"NSC",
"IP",
"Radial_Glia_Cells",
"OPC",
"Astrocyte",
"Ependymal_Cells",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_ML",
"PyNs_RSC_DL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"Chandelier_Cells"
)



pdf("output/Signac/DotPlot_SCT_multiome_WT_Bap1KO_QCV3_label.pdf", width=11, height=4.5)
DotPlot(multiome_WT_Bap1KO_QCV3.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/Signac/DotPlot_SCT_multiome_WT_Bap1KO_QCV3_label_vertical.pdf", width=11, height=4.5)
DotPlot(multiome_WT_Bap1KO_QCV3.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()




############ V3 naming

Cluster1 = PyNs_SubC_CA23 (subiculum PyNs)
Cluster2 = IN_1 (interneuron)
Cluster3 = SubC_1 (subiculum)
Cluster4 = PyNs_SubC_CA1 (subiculum PyNs)
Cluster5 = PyNs_RSC_UL (Retrosplenial Cortical Pyramidal neurons, upper layer)
Cluster6 = DG_GC D(entate Gyrus granule cells)
Cluster7 = PyNs_RSC_MDL (Retrosplenial Cortical Pyramidal neurons, middle/deep layer)
Cluster8 = NSC_proliferative_1 (Neural Stem Cells)
Cluster9 = SubC_2 (subiculum)
Cluster10 = IN_2 (interneuron)
Cluster11 = NSC_quiescent (Neural Stem Cells)
Cluster12 = IN_SubC (subiculum interneurons) ???
Cluster13 = IP (Intermediate Progenitors)
Cluster14 = NSC_proliferative_2 (Neural Stem Cells)
Cluster15 = CR (Cajal Retzius)
Cluster16 = OPC (Oligodendrocyte progenitor cells)
Cluster17 = Meningeal_Cells
Cluster18 = Radial_Glia_Cells
Cluster19 = Microglia




new.cluster.ids <- c(
  "PyNs_SubC_CA23",
  "IN_1",
  "SubC_1",
  "PyNs_SubC_CA1",
  "PyNs_RSC_UL",
  "DG_GC",
  "PyNs_RSC_MDL",
  "NSC_proliferative_1",
  "SubC_2",
  "IN_2",
  "NSC_quiescent",
  "IN_SubC",
  "IP",
  "NSC_proliferative_2",
  "CR",
  "OPC",
  "Meningeal_Cells",
  "Radial_Glia_Cells",
  "Microglia"
)

DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA"

# put the level of my seurat object as the seurat cluster nb
new_levels <- unique(multiome_WT_Bap1KO_QCV2vC1.sct$seurat_clusters)
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- multiome_WT_Bap1KO_QCV2vC1.sct$seurat_clusters


# rename clsuter
names(new.cluster.ids) <- levels(multiome_WT_Bap1KO_QCV2vC1.sct)
multiome_WT_Bap1KO_QCV2vC1.sct <- RenameIdents(multiome_WT_Bap1KO_QCV2vC1.sct, new.cluster.ids)

multiome_WT_Bap1KO_QCV2vC1.sct$cluster.annot <- Idents(multiome_WT_Bap1KO_QCV2vC1.sct) # create a new slot in my seurat object


pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000_label.pdf", width=12, height=6)
DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap", split.by = "orig.ident", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()


pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000_noSplit_label.pdf", width=7, height=5)
DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 4)
dev.off()


multiome_WT_Bap1KO_QCV2vC1.sct$seurat_clusters <- Idents(multiome_WT_Bap1KO_QCV2vC1.sct) # create a new slot in my seurat object

pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000_noSplit_labelNumeric.pdf", width=7, height=5)
DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.3, label.size = 4)
dev.off()



# All in dotplot
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "SCT"

Neural Stem Cells (NSC) = Pax6 (should have more cell types)
Intermediate Progenitors (IP) = Eomes
Dentate Gyrus Granucle Cells (DG) = Prox1, Neurod1, Sema5a
CA1 = Cck, Insm1
CA3 = Crym, Snca, Nrp2
Pyramidal neurons deep layer (DL) = Tac2, Hs3st1, Nrn1
Pyramidal neurons middle layer (ML) = Pantr1, Igfbpl1, Frmd4b
Pyramidal neurons upper layer (UL) = Satb2, Itpr1
Interneurons (IN) = Gad1, Grin2d, Reln, Calb1, Npy, Gria3, Lhx6
Cajal Retzius (CR) = Lhx1
Subiculum (SubC) = Nts, Nr4a2, Lmo3, B3gat1
Microglia = Csf1r, Gpr34, Gpr183, Cx3cr1
OPC = Pdgfra, Olig1
Ependymal_Cells = Foxj1, Cfap44, Dynlrb2
Astrocyte = Aqp4, Gfap
Radial Glia Cells = Gli3, Pdgfd
Chandellier cell = Ntf3, Prkg1, Slc12a5
Meningeal cells = Aldh1a2, Vtn, Lum, Foxc1, Igf2

all_markers <- c(
"Pax6",
"Eomes",
"Gli3", "Pdgfd",
"Pdgfra", "Olig1",
"Foxj1", "Cfap44", "Dynlrb2",
"Csf1r", "Gpr34", "Gpr183", "Cx3cr1",
"Aldh1a2", "Vtn", "Lum", "Foxc1", "Igf2",
"Lhx1",
"Prox1", "Neurod1", "Sema5a",
"Cck", "Insm1",
"Crym", "Snca", "Nrp2",
"Satb2", "Itpr1",
"Pantr1", "Igfbpl1", "Frmd4b",
"Tac2", "Hs3st1", "Nrn1",
"Nts", "Nr4a2", "Lmo3", "B3gat1",
"Gad1", "Grin2d", "Reln", "Calb1", "Npy", "Gria3", "Lhx6",
 "Ntf3", "Prkg1", "Slc12a5"
)





levels(multiome_WT_Bap1KO_QCV2vC1.sct) <- c(
"NSC_quiescent",
"NSC_proliferative_1",
"NSC_proliferative_2",
"IP",
"Radial_Glia_Cells",
"OPC",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_MDL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"IN_SubC"
)



pdf("output/Signac/DotPlot_SCT_multiome_WT_Bap1KO_QCV2vC1_label.pdf", width=11, height=4.5)
DotPlot(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/Signac/DotPlot_SCT_multiome_WT_Bap1KO_QCV2vC1_label_vertical.pdf", width=11, height=4.5)
DotPlot(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()


########################################################


# Cell type proportion

Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- multiome_WT_Bap1KO_QCV2vC1.sct$cluster.annot

pt <- table(Idents(multiome_WT_Bap1KO_QCV2vC1.sct), multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident)
pt <- as.data.frame(pt)
pt <- pt %>%
  group_by(Var2) %>%
  mutate(Proportion = Freq / sum(Freq))

pt$Var1 <- as.character(pt$Var1)

pt$Var1 <- factor(pt$Var1, levels = c("NSC_quiescent",
"NSC_proliferative_1",
"NSC_proliferative_2",
"IP",
"Radial_Glia_Cells",
"OPC",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_MDL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"IN_SubC")) 


#pdf("output/Signac/cellTypeProp_SCT_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.pdf", width=5, height=5)
pdf("output/Signac/cellTypeProp_SCT_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.pdf", width=4, height=5)

ggplot(pt, aes(x = Var2, y = Proportion, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)), 
            position = position_fill(vjust = 0.5), size = 2.5) +
  theme_bw()
dev.off()


# Cell cycle proportion per cluster
## Using numeric cluster annotation

plot_cell_cycle_per_cluster <- function(multiome_WT_Bap1KO_QCV2vC1.sct, output_dir) {
  clusters <- unique(multiome_WT_Bap1KO_QCV2vC1.sct$seurat_clusters)
  for (cluster in clusters) {
    data <- multiome_WT_Bap1KO_QCV2vC1.sct@meta.data %>%
      dplyr::filter(seurat_clusters == cluster) %>%
      group_by(orig.ident, Phase) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      group_by(orig.ident) %>%
      mutate(proportion = count / sum(count)) %>%
      ungroup()

    plot <- ggplot(data, aes(x = orig.ident, y = proportion, fill = Phase)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent) +
      labs(title = paste("Cluster", cluster), x = "Genotype", y = "Proportion (%)") +
      theme_bw() +
      scale_fill_manual(values = c("G1" = "#1f77b4", "G2M" = "#ff7f0e", "S" = "#2ca02c")) +
      geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
                position = position_fill(vjust = 0.5), size = 5)

    # Save plot to PDF
    pdf(paste0(output_dir, "cellCycle_Cluster_QCV2vC1_dim40kparam42res065algo4feat2000_", cluster, ".pdf"), width = 5, height = 6)
    print(plot)
    dev.off()
  }
}
plot_cell_cycle_per_cluster(multiome_WT_Bap1KO_QCV2vC1.sct, output_dir = "output/Signac/")




# SAVE #########################################################################################
## saveRDS(multiome_WT_Bap1KO_QCV3.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV3.sct_numeric_label.rds") 
## saveRDS(multiome_WT_Bap1KO_QCV2vC1.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric_label.rds") 
##
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric_label.rds") 
################################################################################################


#overlapping orig.ident
pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV3_numeric_overlap.pdf", width=6, height=5)
DimPlot(multiome_WT_Bap1KO_QCV3.sct, reduction = "umap", group.by = "orig.ident", pt.size = 0.000001, cols = c("blue","red"))
dev.off()




######################################################################################################
## Pre-processing ATAC  ####################################################################
######################################################################################################

multiome_WT_Bap1KO_QCV3.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV3.sct_numeric_label.rds")

multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric_label.rds")

###  pre-processing and dimensional reductio
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "ATAC"
multiome_WT_Bap1KO_QCV2vC1.sct <- RunTFIDF(multiome_WT_Bap1KO_QCV2vC1.sct)
multiome_WT_Bap1KO_QCV2vC1.sct <- FindTopFeatures(multiome_WT_Bap1KO_QCV2vC1.sct, min.cutoff = 'q0')
multiome_WT_Bap1KO_QCV2vC1.sct <- RunSVD(multiome_WT_Bap1KO_QCV2vC1.sct)
multiome_WT_Bap1KO_QCV2vC1.sct <- RunUMAP(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = 'lsi', dims = 2:40, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # We exclude the first dimension as this is typically correlated with sequencing depth


### WNN graph, representing a weighted combination of RNA and ATAC-seq modalities

multiome_WT_Bap1KO_QCV2vC1.sct <- FindMultiModalNeighbors(multiome_WT_Bap1KO_QCV2vC1.sct, reduction.list = list("pca", "lsi"), dims.list = list(1:40, 2:40))
multiome_WT_Bap1KO_QCV2vC1.sct <- RunUMAP(multiome_WT_Bap1KO_QCV2vC1.sct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
multiome_WT_Bap1KO_QCV2vC1.sct <- FindClusters(multiome_WT_Bap1KO_QCV2vC1.sct, graph.name = "wsnn", algorithm = 3, verbose = TRUE)

# plot gene expression, ATAC-seq, or WNN analysis
p1 <- DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap", group.by = "cluster.annot", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "umap.atac", group.by = "cluster.annot", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "wnn.umap", group.by = "cluster.annot", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV2vC1_RNAATACWNN_ATACdim240.pdf", width=12, height=5)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident, levels = c("multiome_WT", "multiome_Bap1KO")) # Reorder untreated 1st

pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV2vC1_WNN_ATACdim240.pdf", width=5, height=5)
p = DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "wnn.umap", group.by = "cluster.annot", label = FALSE, label.size = 3, repel = TRUE) + ggtitle("Cell type") + NoLegend()
LabelClusters(p, id = "cluster.annot", fontface = "bold", color = "black", size = 3)
dev.off()
pdf("output/Signac/UMAP_multiome_WT_Bap1KO_QCV2vC1_WNN_ATACdim240_genotype.pdf", width=5, height=5)
DimPlot(multiome_WT_Bap1KO_QCV2vC1.sct, reduction = "wnn.umap", group.by = "orig.ident", label = FALSE, cols = c("blue", "red")) + ggtitle("Genotype")  + NoLegend()
dev.off()


#--> Different dim for ATAC tested; 2:40 best


# Calculate gene activity (count ATAC peak within gene and promoter)
GeneActivity = GeneActivity(
  multiome_WT_Bap1KO_QCV2vC1.sct,
  assay = "ATAC",
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = NULL,
  max.width = NULL,
  process_n = 2000,
  gene.id = FALSE,
  verbose = TRUE
)
## Add gene activity as new assay

multiome_WT_Bap1KO_QCV2vC1.sct[["GeneActivity"]] <- CreateAssayObject(counts = GeneActivity)

# SAVE ##########################################################################################
## saveRDS(multiome_WT_Bap1KO_QCV2vC1.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivity.sct_numeric_label.rds") 
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivity.sct_numeric_label.rds")
# In case I re run with correct1 (logNorm/ScaleData DEGs)
## saveRDS(multiome_WT_Bap1KO_QCV2vC1.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivity.sct_numeric_label.rds") 
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivity.sct_numeric_label.rds")
##########################################################################################



## Plot genes RNA ATAC

DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "SCT"
p1 <- FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "Eomes", reduction = "wnn.umap", max.cutoff = 1, cols = c("grey", "red")) + 
  ggtitle("RNA: Eomes") + theme(plot.title = element_text(hjust = 0.5))
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "GeneActivity"
p2 <- FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "Eomes", reduction = "wnn.umap", max.cutoff = 100, cols = c("grey", "red")) + 
  ggtitle("ATAC: Pax6") + theme(plot.title = element_text(hjust = 0.5))
pdf("output/Signac/FeaturePlot_WNNreduction_WT_Bap1KO_Eomes.pdf", width=10, height=5)
(p1 | p2) & NoLegend()  # The | operator ensures the plots are in the same rowdev.off()
dev.off()

DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "SCT"
p1 <- FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "Eomes", reduction = "umap", max.cutoff = 1, cols = c("grey", "red")) + 
  ggtitle("RNA: Eomes") + theme(plot.title = element_text(hjust = 0.5))
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "GeneActivity"
p2 <- FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "Eomes", reduction = "umap", max.cutoff = 100, cols = c("grey", "red")) + 
  ggtitle("ATAC: Eomes") + theme(plot.title = element_text(hjust = 0.5))
pdf("output/Signac/FeaturePlot_UMAPreduction_WT_Bap1KO_Eomes.pdf", width=10, height=5)
(p1 | p2) & NoLegend()  # The | operator ensures the plots are in the same rowdev.off()
dev.off()



# All in dotplot

## RNA
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- multiome_WT_Bap1KO_QCV2vC1.sct$cluster.annot
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "SCT"


all_markers <- c(
  "Slco1c1", "Gli3","Gm29260", # cl11 NSC_quiescent remove "Tnc", "8030451A03Rik"
  "Tnc", "8030451A03Rik","Egfr","Bcan", "Lrig1", # cl8 NSC_proliferative_1
  "Adamts19", "Gm29521","Adgrv1","Gm47520","Lef1", # cl14 NSC_proliferative_2
  "Eomes","Mfap4","Rmst","Chd7","Gm11266", # cl13 IP
  "Lmx1a","Ttc21a","Ccdc170","Gm19301","Wdr63", # cl18 Radial Glia Cells
  "Pdgfra","Olig2","Sox10","Bcas1","Smoc1","Gm10863", #cl16 OPC
  "Fli1","Mir142hg","Csf1r","Fyb","Inpp5d", # cl19 microglia
  "Arhgap29","Cped1","Slc6a20a","Col3a1","Col1a2", # cl17 Meningeal cells
  "Ndnf","Ebf3","Cacna2d2","Plekha7","Tmem163", # cl15 CR
  "Prox1","Sema5a","Adamts18","Cntnap5a","Ptpro",  # cl6 DG_GC
  "Fam189a1","Grin2a","Cck","Gm20754","Hcn1", # cl4 PyNs_SubC_CA1
  "Gm32647","6530403H02Rik","Rgs6","Zfp385b","Trpc7", # cl1 PyNs_SubC_CA23
  "Satb2","Tshz2","9130024F11Rik","Itpr1","Tshz3", # cl5 PyNs_RSC_UL
  "Sema3c","Meis2","Adgrb3","Auts2","Kcnq3", # cl7 PyNs_RSC_MDL
  "Ntng1","Cntn3","Ldb2","Schip1","A230006K03Rik", # cl3 SubC_1 
  "Hs3st4","Cobll1","Lrrtm3","Ctnna3","Grm8", # cl9 SubC_2
  "Nxph1","Maf","Gad2","Kcnc2","Erbb4", # cl2 IN_1
  "Adarb2","Dlx6os1","Igf1","Sorcs3", # cl10 IN_2 remove Erbb4
  "Camk2a","Cpne4","St18","Cpne7","Htr2c" # cl12 IN_SubC
)


# --> Filter the top 3 
all_markers <- c(
  "Slco1c1", "Gli3","Gm29260", # cl11 NSC_quiescent remove "Tnc", "8030451A03Rik"
  "Tnc", "8030451A03Rik","Egfr", # cl8 NSC_proliferative_1
  "Adamts19", "Gm29521","Adgrv1", # cl14 NSC_proliferative_2
  "Eomes","Mfap4","Rmst", # cl13 IP
  "Lmx1a","Ttc21a","Ccdc170", # cl18 Radial Glia Cells
  "Pdgfra","Olig2","Sox10", #cl16 OPC
  "Fli1","Mir142hg","Csf1r", # cl19 microglia
  "Arhgap29","Cped1","Slc6a20a", # cl17 Meningeal cells
  "Ndnf","Ebf3","Cacna2d2", # cl15 CR
  "Prox1","Sema5a","Adamts18",  # cl6 DG_GC
  "Fam189a1","Grin2a","Cck", # cl4 PyNs_SubC_CA1
  "Gm32647","6530403H02Rik","Rgs6", # cl1 PyNs_SubC_CA23
  "Satb2","Tshz2","9130024F11Rik", # cl5 PyNs_RSC_UL
  "Sema3c","Meis2","Kcnq3",# cl7 PyNs_RSC_MDL
  "Ntng1","Cntn3","Ldb2", # cl3 SubC_1 
  "Hs3st4","Cobll1","Lrrtm3", # cl9 SubC_2
  "Nxph1","Maf","Gad2", # cl2 IN_1
  "Adarb2","Dlx6os1","Igf1", # cl10 IN_2 remove Erbb4
  "Camk2a","Cpne4","St18" # cl12 IN_SubC
)


levels(multiome_WT_Bap1KO_QCV2vC1.sct) <- c(
"NSC_quiescent",
"NSC_proliferative_1",
"NSC_proliferative_2",
"IP",
"Radial_Glia_Cells",
"OPC",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_MDL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"IN_SubC"
)



pdf("output/Signac/DotPlot_SCT_multiome_WT_Bap1KO_QCV2vC1_label_vertical_top3.pdf", width=12, height=4.5)
DotPlot(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()


## GeneActivity

Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- multiome_WT_Bap1KO_QCV2vC1.sct$cluster.annot
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "SCT"


all_markers <- c(
  "Slco1c1", "Gli3","Gm29260", # cl11 NSC_quiescent remove "Tnc", "8030451A03Rik"
  "Tnc", "8030451A03Rik","Egfr","Bcan", "Lrig1", # cl8 NSC_proliferative_1
  "Adamts19", "Gm29521","Adgrv1","Gm47520","Lef1", # cl14 NSC_proliferative_2
  "Eomes","Mfap4","Rmst","Chd7","Gm11266", # cl13 IP
  "Lmx1a","Ttc21a","Ccdc170","Gm19301","Wdr63", # cl18 Radial Glia Cells
  "Pdgfra","Olig2","Sox10","Bcas1","Smoc1","Gm10863", #cl16 OPC
  "Fli1","Mir142hg","Csf1r","Fyb","Inpp5d", # cl19 microglia
  "Arhgap29","Cped1","Slc6a20a","Col3a1","Col1a2", # cl17 Meningeal cells
  "Ndnf","Ebf3","Cacna2d2","Plekha7","Tmem163", # cl15 CR
  "Prox1","Sema5a","Adamts18","Cntnap5a","Ptpro",  # cl6 DG_GC
  "Fam189a1","Grin2a","Cck","Gm20754","Hcn1", # cl4 PyNs_SubC_CA1
  "Gm32647","6530403H02Rik","Rgs6","Zfp385b","Trpc7", # cl1 PyNs_SubC_CA23
  "Satb2","Tshz2","9130024F11Rik","Itpr1","Tshz3", # cl5 PyNs_RSC_UL
  "Sema3c","Meis2","Adgrb3","Auts2","Kcnq3", # cl7 PyNs_RSC_MDL
  "Ntng1","Cntn3","Ldb2","Schip1","A230006K03Rik", # cl3 SubC_1 
  "Hs3st4","Cobll1","Lrrtm3","Ctnna3","Grm8", # cl9 SubC_2
  "Nxph1","Maf","Gad2","Kcnc2","Erbb4", # cl2 IN_1
  "Adarb2","Dlx6os1","Igf1","Sorcs3", # cl10 IN_2 remove Erbb4
  "Camk2a","Cpne4","St18","Cpne7","Htr2c" # cl12 IN_SubC
)


# --> Filter the top 3 
all_markers <- c(
  "Slco1c1", "Gli3","Gm29260", # cl11 NSC_quiescent remove "Tnc", "8030451A03Rik"
  "Tnc", "8030451A03Rik","Egfr", # cl8 NSC_proliferative_1
  "Adamts19", "Gm29521","Adgrv1", # cl14 NSC_proliferative_2
  "Eomes","Mfap4","Rmst", # cl13 IP
  "Lmx1a","Ttc21a","Ccdc170", # cl18 Radial Glia Cells
  "Pdgfra","Olig2","Sox10", #cl16 OPC
  "Fli1","Mir142hg","Csf1r", # cl19 microglia
  "Arhgap29","Cped1","Slc6a20a", # cl17 Meningeal cells
  "Ndnf","Ebf3","Cacna2d2", # cl15 CR
  "Prox1","Sema5a","Adamts18",  # cl6 DG_GC
  "Fam189a1","Grin2a","Cck", # cl4 PyNs_SubC_CA1
  "Gm32647","6530403H02Rik","Rgs6", # cl1 PyNs_SubC_CA23
  "Satb2","Tshz2","9130024F11Rik", # cl5 PyNs_RSC_UL
  "Sema3c","Meis2","Kcnq3",# cl7 PyNs_RSC_MDL
  "Ntng1","Cntn3","Ldb2", # cl3 SubC_1 
  "Hs3st4","Cobll1","Lrrtm3", # cl9 SubC_2
  "Nxph1","Maf","Gad2", # cl2 IN_1
  "Adarb2","Dlx6os1","Igf1", # cl10 IN_2 remove Erbb4
  "Camk2a","Cpne4","St18" # cl12 IN_SubC
)


levels(multiome_WT_Bap1KO_QCV2vC1.sct) <- c(
"NSC_quiescent",
"NSC_proliferative_1",
"NSC_proliferative_2",
"IP",
"Radial_Glia_Cells",
"OPC",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_MDL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"IN_SubC"
)



pdf("output/Signac/DotPlot_GeneActivity_multiome_WT_Bap1KO_QCV2vC1_label_vertical_top3.pdf", width=12, height=4.5)
DotPlot(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()




# some plots



pdf("output/Signac/CoveragePlot-Pdgfra.pdf", width=5, height=5)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Pdgfra', features = 'Pdgfra', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, group.by = "cluster.annot" )
dev.off()


pdf("output/Signac/CoveragePlot-Gad1.pdf", width=5, height=5)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Gad1', features = 'Gad1', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, group.by = "cluster.annot" )
dev.off()

pdf("output/Signac/CoveragePlot-Lhx1.pdf", width=5, height=5)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Lhx1', features = 'Lhx1', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, group.by = "cluster.annot" )
dev.off()

pdf("output/Signac/CoveragePlot-Grin2a.pdf", width=5, height=5)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Grin2a', features = 'Grin2a', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, group.by = "cluster.annot" )
dev.off()


pdf("output/Signac/CoveragePlot-Trhde.pdf", width=5, height=5)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Trhde', features = 'Trhde', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, group.by = "cluster.annot" )
dev.off()

pdf("output/Signac/CoveragePlot-Bap1.pdf", width=5, height=5)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Bap1', features = 'Bap1', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, group.by = "cluster.annot" )
dev.off()


pdf("output/Signac/CoveragePlot-Bap1.pdf", width=5, height=5)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Bap1', features = 'Bap1', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, group.by = "cluster.annot" )
dev.off()


# WT vs Bap1KO
## Reorder metrics
multiome_WT_Bap1KO_QCV2vC1.sct$combined_ident <- paste(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident, multiome_WT_Bap1KO_QCV2vC1.sct$cluster.annot, sep = "_")

Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "combined_ident"

multiome_WT_Bap1KO_QCV2vC1.sct$combined_ident <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$combined_ident, 
                                                     levels = c("multiome_WT_NSC_quiescent", "multiome_Bap1KO_NSC_quiescent",
                                                                "multiome_WT_NSC_proliferative_1" ,"multiome_Bap1KO_NSC_proliferative_1",
                                                                "multiome_WT_NSC_proliferative_2" ,"multiome_Bap1KO_NSC_proliferative_2",
                                                                "multiome_WT_IP" ,"multiome_Bap1KO_IP",
                                                                "multiome_WT_Radial_Glia_Cells" ,"multiome_Bap1KO_Radial_Glia_Cells",
                                                                "multiome_WT_OPC" ,"multiome_Bap1KO_OPC",
                                                                "multiome_WT_Microglia" ,"multiome_Bap1KO_Microglia",
                                                                "multiome_WT_Meningeal_Cells" ,"multiome_Bap1KO_Meningeal_Cells",
                                                                "multiome_WT_CR" ,"multiome_Bap1KO_CR",
                                                                "multiome_WT_DG_GC" ,"multiome_Bap1KO_DG_GC",
                                                                "multiome_WT_PyNs_SubC_CA1" ,"multiome_Bap1KO_PyNs_SubC_CA1",
                                                                "multiome_WT_PyNs_SubC_CA23" ,"multiome_Bap1KO_PyNs_SubC_CA23",
                                                                "multiome_WT_PyNs_RSC_UL" ,"multiome_Bap1KO_PyNs_RSC_UL",
                                                                "multiome_WT_PyNs_RSC_MDL" ,"multiome_Bap1KO_PyNs_RSC_MDL",
                                                                "multiome_WT_SubC_1" ,"multiome_Bap1KO_SubC_1",
                                                                "multiome_WT_SubC_2" ,"multiome_Bap1KO_SubC_2",
                                                                "multiome_WT_IN_1" ,"multiome_Bap1KO_IN_1",
                                                                "multiome_WT_IN_2" ,"multiome_Bap1KO_IN_2",
                                                                "multiome_WT_IN_SubC" ,"multiome_Bap1KO_IN_SubC"))


pdf("output/Signac/CoveragePlot_condition-Epha3.pdf", width=6, height=9)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Epha3', features = 'Epha3', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, group.by = "combined_ident")
dev.off()

# Mdga2
pdf("output/Signac/CoveragePlot_condition-Mdga2-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Mdga2', features = 'Mdga2', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Mdga2_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr12-67217092-67224465"), features = 'Mdga2', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()

# Lrrc7
pdf("output/Signac/CoveragePlot_condition-Lrrc7-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Lrrc7', features = 'Lrrc7', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Lrrc7_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr3-158559019-158563930"), features = 'Lrrc7', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()

# Tenm4
pdf("output/Signac/CoveragePlot_condition-Tenm4-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Tenm4', features = 'Tenm4', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Tenm4_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr7-96208102-96217029"), features = 'Tenm4', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Tenm4_TSS2-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr7-96170754-96172985"), features = 'Tenm4', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()

# Rmst
pdf("output/Signac/CoveragePlot_condition-Rmst-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Rmst', features = 'Rmst', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Lrrc7_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr3-158559019-158563930"), features = 'Lrrc7', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()


# Dscam
pdf("output/Signac/CoveragePlot_condition-Dscam-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Dscam', features = 'Dscam', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Dscam_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr16-97169399-97172371"), features = 'Dscam', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()


# Tenm2
pdf("output/Signac/CoveragePlot_condition-Tenm2-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Tenm2', features = 'Tenm2', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Dscam_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr16-97169399-97172371"), features = 'Dscam', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()

# Slc24a2
pdf("output/Signac/CoveragePlot_condition-Slc24a2-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Slc24a2', features = 'Slc24a2', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Slc24a2_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr4-87226581-87231384"), features = 'Slc24a2', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()

# Itpr1
pdf("output/Signac/CoveragePlot_condition-Itpr1-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Itpr1', features = 'Itpr1', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Itpr1_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr6-108212957-108214105"), features = 'Itpr1', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()


# Kirrel3
pdf("output/Signac/CoveragePlot_condition-Kirrel3-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Kirrel3', features = 'Kirrel3', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Itpr1_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr6-108212957-108214105"), features = 'Kirrel3', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()

# Sema6a
pdf("output/Signac/CoveragePlot_condition-Sema6a-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region = 'Sema6a', features = 'Sema6a', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"),   extend.upstream = 1000)
dev.off()
pdf("output/Signac/CoveragePlot_condition-Sema6a_TSS-DG_GC.pdf", width=5, height=2)
CoveragePlot(multiome_WT_Bap1KO_QCV2vC1.sct, region =  c("chr18-47364261-47369828"), features = 'Sema6a', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, links = FALSE, group.by = "combined_ident", idents = c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC"))
dev.off()




# Identify diff access regions (DAR)
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- 'ATAC'

# wilcox is the default option for test.use
NSC_quiescent_WTvsBap1KO <- FindMarkers(
  object = multiome_WT_Bap1KO_QCV2vC1.sct,
  ident.1 = "multiome_WT_NSC_quiescent",
  ident.2 = "multiome_Bap1KO_NSC_quiescent",
  test.use = 'wilcox',
  min.pct = 0.1
) %>%
  add_column(cluster = "NSC_quiescent")
NSC_proliferative_1_WTvsBap1KO <- FindMarkers(
  object = multiome_WT_Bap1KO_QCV2vC1.sct,
  ident.1 = "multiome_WT_NSC_proliferative_1",
  ident.2 = "multiome_Bap1KO_NSC_proliferative_1",
  test.use = 'wilcox',
  min.pct = 0.1
)

### Here is to automatize the process clsuter per cluster:
#### Define the cluster pairs for comparison
cluster_pairs <- list(
  c("multiome_WT_NSC_quiescent", "multiome_Bap1KO_NSC_quiescent", "NSC_quiescent"),
  c("multiome_WT_NSC_proliferative_1", "multiome_Bap1KO_NSC_proliferative_1", "NSC_proliferative_1"),
  c("multiome_WT_NSC_proliferative_2", "multiome_Bap1KO_NSC_proliferative_2", "NSC_proliferative_2"),
  c("multiome_WT_IP", "multiome_Bap1KO_IP", "IP"),
  c("multiome_WT_Radial_Glia_Cells", "multiome_Bap1KO_Radial_Glia_Cells", "Radial_Glia_Cells"),
  c("multiome_WT_OPC", "multiome_Bap1KO_OPC", "OPC"),
  c("multiome_WT_Microglia", "multiome_Bap1KO_Microglia", "Microglia"),
  c("multiome_WT_Meningeal_Cells", "multiome_Bap1KO_Meningeal_Cells", "Meningeal_Cells"),
  c("multiome_WT_CR", "multiome_Bap1KO_CR", "CR"),
  c("multiome_WT_DG_GC", "multiome_Bap1KO_DG_GC", "DG_GC"),
  c("multiome_WT_PyNs_SubC_CA1", "multiome_Bap1KO_PyNs_SubC_CA1", "PyNs_SubC_CA1"),
  c("multiome_WT_PyNs_SubC_CA23", "multiome_Bap1KO_PyNs_SubC_CA23", "PyNs_SubC_CA23"),
  c("multiome_WT_PyNs_RSC_UL", "multiome_Bap1KO_PyNs_RSC_UL", "PyNs_RSC_UL"),
  c("multiome_WT_PyNs_RSC_MDL", "multiome_Bap1KO_PyNs_RSC_MDL", "PyNs_RSC_MDL"),
  c("multiome_WT_SubC_1", "multiome_Bap1KO_SubC_1", "SubC_1"),
  c("multiome_WT_SubC_2", "multiome_Bap1KO_SubC_2", "SubC_2"),
  c("multiome_WT_IN_1", "multiome_Bap1KO_IN_1", "IN_1"),
  c("multiome_WT_IN_2", "multiome_Bap1KO_IN_2", "IN_2"),
  c("multiome_WT_IN_SubC", "multiome_Bap1KO_IN_SubC", "IN_SubC") )

## Function to run FindMarkers and return a tibble with cluster and gene information
run_find_markers <- function(pair) {
  ident_1 <- pair[1]
  ident_2 <- pair[2]
  cluster_name <- pair[3]
  
  # Run FindMarkers for each pair of clusters
  markers <- FindMarkers(
    object = multiome_WT_Bap1KO_QCV2vC1.sct,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = 'wilcox',
    min.pct = 0.1
  ) %>%
  rownames_to_column(var = "query_region") %>%  # Add gene names as a column
  add_column(cluster = cluster_name)    # Add cluster name as a new column
  
  return(markers)
}

## Combine all FindMarkers results into one tibble
all_markers_tibble <- map_dfr(cluster_pairs, run_find_markers) %>%
  as_tibble()
## save output
write.table(all_markers_tibble, file = "output/Signac/DAR_peaks_QCV2vC1_dim40kparam42res065algo4feat2000.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


## Find closest gene to DAR
DAR_closestGene <- ClosestFeature(multiome_WT_Bap1KO_QCV2vC1.sct, regions = all_markers_tibble$query_region)
DAR_genes = all_markers_tibble %>% 
  left_join(DAR_closestGene) %>%
  as_tibble() %>%
  unique()

## save output
write.table(DAR_genes, file = "output/Signac/DAR_genes_QCV2vC1_dim40kparam42res065algo4feat2000.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########### COUNT UP DOWN PEAK PER CLUSTER ##################################################################
all_markers_tibble = read_tsv("output/Signac/DAR_peaks_QCV2vC1_dim40kparam42res065algo4feat2000.txt")

all_markers_tibble %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarize(
    nb_upregulated = sum(avg_log2FC > 0),
    nb_downregulated = sum(avg_log2FC < 0)
  )

DAR_genes = read_tsv("output/Signac/DAR_genes_QCV2vC1_dim40kparam42res065algo4feat2000.txt")


DAR_genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select(gene_name, avg_log2FC, p_val_adj, cluster) %>%
  unique() %>%
  group_by(cluster) %>%
  summarize(
    nb_upregulated = sum(avg_log2FC > 0),
    nb_downregulated = sum(avg_log2FC < 0)
  )


##############################################################################################################

# Link peaks to genes - Find peaks that are correlated with the expression of nearby genes
## Install the BS genome (GC content, region lengths, and dinucleotide base frequencies for regions in the assay and add to the feature metadata)
library("BSgenome.Mmusculus.UCSC.mm10") # BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

## Add genomic information to Seurat object

### Convert rownames of ATAC counts to GRanges
grange.counts <- StringToGRanges(rownames(multiome_WT_Bap1KO_QCV2vC1.sct[["ATAC"]]), sep = c(":", "-"))
### Filter for standard chromosomes
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- multiome_WT_Bap1KO_QCV2vC1.sct[as.vector(grange.use), ]
### Get annotations for the mouse genome
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
### Adjust chromosome naming style
seqlevelsStyle(annotations) <- 'UCSC'
### Set genome to mm10 (mouse)
genome = 'mm10'

## Run  RegionStats
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "ATAC"
multiome_WT_Bap1KO_QCV2vC1.sct <- RegionStats(
  object = multiome_WT_Bap1KO_QCV2vC1.sct,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  assay = "ATAC",
  verbose = TRUE
)

## Run  LinkPeaks
multiome_WT_Bap1KO_QCV2vC1.sct = LinkPeaks(
  multiome_WT_Bap1KO_QCV2vC1.sct,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  genes.use = NULL,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  gene.id = FALSE,
  verbose = TRUE
)



# SAVE ##########################################################################################
## saveRDS(multiome_WT_Bap1KO_QCV2vC1.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_XXXTOLOADANDRUNLINKPEAKS.rds") 
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_XXXTOLOADANDRUNLINKPEAKS.rds")
# --> JUST A VERIFICATION THAT correct1 did not messedup linkpeaks; run this and check whether the Values are similar.
## saveRDS(LinkPeaks, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds") 
# --> Value is not exactly, the same. Maybe because of random calculation, maybe not, as we don't know lets prefer use correct1 version.
# PAST version:
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")
# NEW correct1:
#saveRDS(multiome_WT_Bap1KO_QCV2vC1.sct, file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds") 
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds")
##########################################################################################



Links = as_tibble(Links(multiome_WT_Bap1KO_QCV2vC1.sct))



### Adjust the pvalue and select positive corr as in the https://www.nature.com/articles/s41467-024-45199-x#Fig2 paper
Links$adjusted_pvalue <- p.adjust(Links$pvalue, method = "BH")
#write.table(Links, file = c("output/Signac/LinkPeaks_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1.txt"),sep="\t", quote=FALSE, row.names=FALSE)

### Isolate signif genes
Links_signif = Links %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%  # 28,201 Link Signif
  dplyr::select(gene) %>%
  unique()
### Reorder the gene based on their cluster max expression

Links %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%  
  dplyr::select(peak) %>% # 22,747 Link peak Signif
  unique()

## log norm and Scale GeneActivity assay
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "GeneActivity"

multiome_WT_Bap1KO_QCV2vC1.sct <- NormalizeData(multiome_WT_Bap1KO_QCV2vC1.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(multiome_WT_Bap1KO_QCV2vC1.sct)
multiome_WT_Bap1KO_QCV2vC1.sct <- ScaleData(multiome_WT_Bap1KO_QCV2vC1.sct, features = all.genes) # zero-centres and scales it



###### Find all markers 
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "cluster.annot"
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA"

Links_markers <- FindAllMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, features = Links_signif$gene, assay = "RNA", only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
###### Identify in which cluster the Links_markers gene is highly express
Links_markers_pval= as_tibble(Links_markers) %>%
  group_by(gene) %>%
  dplyr::filter(p_val == min(p_val)) %>%
  dplyr::select(gene, cluster)
#write.table(Links_markers, file = "output/Signac/srat_multiome_WT_Bap1KO-QCV2vC1_dim40kparam42res065algo4feat2000correct1_noCellCycleRegression-Links_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)



# plot heatmap
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_SCTscaledata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "SCT", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_SCTdata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "SCT", slot= "data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()

pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_RNAdata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", slot= "data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_RNArawdata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", slot= "raw.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()


pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_GeneActivityScaldata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_GeneActivitydata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_GeneActivityrawdata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "raw.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()

# testing aesthetics
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_GeneActivityScaldata-dispminmax1.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , disp.min = -1, disp.max = 1, angle = 0, hjust = 0.5)
dev.off()

pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_GeneActivityScaldata-dispmin25max05.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , disp.min = -2.5, disp.max = 0.5, angle = 0, hjust = 0.5)
dev.off()

pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_GeneActivityScaldata-bluewhitered.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, 
          assay = "GeneActivity", 
          slot= "scale.data", 
          features = Links_markers_pval$gene, 
          group.by = "cluster.annot", 
          angle = 0, 
          hjust = 0.5) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))  # Adjust colors for higher contrast
dev.off()

pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0_GeneActivityScaldata-bluewhiteredBreaks.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, 
          assay = "GeneActivity", 
          slot= "scale.data", 
          features = Links_markers_pval$gene, 
          group.by = "cluster.annot", 
          angle = 0, 
          hjust = 0.5) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"), 
                       breaks = seq(-2, 2, by = 0.5))  # Modify breaks to increase color contrast
dev.off()


# final aesthetics - pretty

pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000correct1_LinksPadj05Score0_SCTscaledata_pretty.pdf", width=6, height=3)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "SCT", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE)
dev.off()

pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000correct1_LinksPadj05Score0_GeneActivtiyscaledata_pretty.pdf", width=6, height=3)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE)
dev.off()
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000correct1_LinksPadj05Score0_GeneActivtiyscaledata_pretty1.pdf", width=6, height=3)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE, disp.max = 2)
dev.off()
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000correct1_LinksPadj05Score0_GeneActivtiyscaledata_pretty2.pdf", width=6, height=3)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5, draw.lines = FALSE, label = FALSE, disp.max = 2, disp.min = -2)
dev.off()



# more stringeant filtering

### Adjust the pvalue and select positive corr as in the https://www.nature.com/articles/s41467-024-45199-x#Fig2 paper
Links$adjusted_pvalue <- p.adjust(Links$pvalue, method = "BH")
### Isolate signif genes
Links_signif = Links %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%
  dplyr::select(gene) %>%
  unique()
### Reorder the gene based on their cluster max expression


###### Find all markers 
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "cluster.annot"
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA"
Links_markers <- FindAllMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, features = Links_signif$gene, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) # V1 min.pct = 0.01, logfc.threshold = 0.1; V2 min.pct = 0.25, logfc.threshold = 0.25
###### Identify in which cluster the Links_markers gene is highly express
Links_markers_pval= as_tibble(Links_markers) %>%
  group_by(gene) %>%
  dplyr::filter(p_val == min(p_val)) %>%
  dplyr::select(gene, cluster)
#V1: write.table(Links_markers, file = "output/Signac/srat_multiome_WT_Bap1KO-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-Links_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#V2: write.table(Links_markers, file = "output/Signac/srat_multiome_WT_Bap1KO-QCV2vC1_dim40kparam42res065algo4feat2000_noCellCycleRegression-Links_markersV2.txt", sep = "\t", quote = FALSE, row.names = TRUE)

pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0FindAllMarkersminpct025logfc025_SCTscaledata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "SCT", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()
pdf("output/Signac/DoHeatmap_QCV2vC1_dim40kparam42res065algo4feat2000_LinksPadj05Score0FindAllMarkersminpct025logfc025_GeneActivityscaledata.pdf", width=8, height=4)
DoHeatmap(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", slot= "scale.data", features = Links_markers_pval$gene, group.by = "cluster.annot" , angle = 0, hjust = 0.5)
dev.off()







# Select top 2-5 Link genes for each cluster in dotplot


top_genes <- Links_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 3) %>%
  ungroup()

# View the result
top_genes




all_markers <- c(
  "Hgf", "Pdgfrb","Clu", # cl11 NSC_quiescent remove "Tnc", "8030451A03Rik"
  "Aqp4", "Aldh1l1","Acsbg1", # cl8 NSC_proliferative_1
  "Adamts19", "Gm29521","Slc2a9", # cl14 NSC_proliferative_2
  "Eomes","Mfap4","Rmst", # cl13 IP
  "Ccdc170","Stpg1","Dynlrb2", # cl18 Radial Glia Cells
  "Sox10","Bcas1","Pdgfra", #cl16 OPC
  "Vav1","Ikzf1","Mir142hg", # cl19 microglia
  "Tbx18","Slc6a13","Foxd1", # cl17 Meningeal cells
  "Lhx1","Lhx5","Lhx1os", # cl15 CR
  "Sema5a","Prox1","Rarb",  # cl6 DG_GC
  "Cck","Fam189a1","Grin2a", # cl4 PyNs_SubC_CA1
  "Trpc7","6530403H02Rik","Tgfbr3", # cl1 PyNs_SubC_CA23
  "9130024F11Rik","Satb2","Ntf3", # cl5 PyNs_RSC_UL
  "Sema3c","Meis2","Kcnq3",# cl7 PyNs_RSC_MDL
  "Cntn3","Nrgn","Ntng1", # cl3 SubC_1 
  "Cobll1","Nr4a2","Grm8", # cl9 SubC_2
  "Lhx6","Nxph2","Cntnap3", # cl2 IN_1
  "Adarb2","Igf1","Npas1", # cl10 IN_2 remove Erbb4
  "Cpne7","St18","Htr2c" # cl12 IN_SubC
)


levels(multiome_WT_Bap1KO_QCV2vC1.sct) <- c(
"NSC_quiescent",
"NSC_proliferative_1",
"NSC_proliferative_2",
"IP",
"Radial_Glia_Cells",
"OPC",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_MDL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"IN_SubC"
)



pdf("output/Signac/DotPlot_GeneActivity_multiome_WT_Bap1KO_QCV2vC1_label_vertical_top3Links.pdf", width=12, height=4.5)
DotPlot(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "GeneActivity", features = all_markers)   +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") + # option="magma"
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()




pdf("output/Signac/DotPlot_SCT_multiome_WT_Bap1KO_QCV2vC1_label_vertical_top3Links.pdf", width=12, height=4.5)
DotPlot(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "SCT", features = all_markers)    +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis() + # option="magma"
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()





# DEGs number colored in a UMAP
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "cluster.annot"

DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())
## List of cell types
cell_types <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19")
## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/Signac/", cell_type, "-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt", sep = "")
  deg_data <- read.table(file_name, header = TRUE, sep = "\t") ## Read the DEGs data
  num_degs <- sum(deg_data$p_val_adj < 0.05 & abs(deg_data$avg_log2FC) > 0.5) ## Count the number of significant DEGs
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))  ## Append to the summary table
}
DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9", "cluster10", "cluster11", "cluster12", "cluster13", "cluster14", "cluster15", "cluster16", "cluster17", "cluster18", "cluster19")) 
DEG_count$Cluster_Number <- as.numeric(sub("cluster", "", DEG_count$Cell_Type))
DEG_count$Cluster_Number <- factor(DEG_count$Cluster_Number, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")) 

cluster_to_cell_name <- c(
  "cluster1" = "PyNs_SubC_CA23",
  "cluster2" = "IN_1",
  "cluster3" = "SubC_1",
  "cluster4" = "PyNs_SubC_CA1",
  "cluster5" = "PyNs_RSC_UL",
  "cluster6" = "DG_GC",
  "cluster7" = "PyNs_RSC_MDL",
  "cluster8" = "NSC_proliferative_1",
  "cluster9" = "SubC_2",
  "cluster10" = "IN_2",
  "cluster11" = "NSC_quiescent",
  "cluster12" = "IN_SubC",
  "cluster13" = "IP",
  "cluster14" = "NSC_proliferative_2",
  "cluster15" = "CR",
  "cluster16" = "OPC",
  "cluster17" = "Meningeal_Cells",
  "cluster18" = "Radial_Glia_Cells",
  "cluster19" = "Microglia"
)

# Add the new column by mapping the Cell_Type values
DEG_count$Cell_Name <- cluster_to_cell_name[DEG_count$Cell_Type]






# Add DEG information to my seurat object - DEG_count
cell_clusters <- multiome_WT_Bap1KO_QCV2vC1.sct@meta.data$cluster.annot
names(cell_clusters) <- rownames(multiome_WT_Bap1KO_QCV2vC1.sct@meta.data)
DEG_named_vector <- DEG_count$Num_DEGs[match(cell_clusters, DEG_count$Cell_Name)]
names(DEG_named_vector) <- names(cell_clusters)
# Integrate DEG values into the Seurat object
multiome_WT_Bap1KO_QCV2vC1.sct <- AddMetaData(multiome_WT_Bap1KO_QCV2vC1.sct, metadata = DEG_named_vector, col.name = "DEG")
# Create a UMAP plot colored by qval values
pdf("output/Signac/FeaturePlot_WNN_multiome_WT_Bap1KO_QCV2vC1_DEG.pdf", width=6, height=6)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "DEG", pt.size = 0.5, reduction = "wnn.umap") +
  scale_colour_viridis(option="magma") # 
dev.off()

# WNN REDUCTION ############
# Add values on the heatmap
## Extract UMAP coordinates
umap_coordinates <- as.data.frame(multiome_WT_Bap1KO_QCV2vC1.sct@reductions$wnn.umap@cell.embeddings)
umap_coordinates$cluster <- multiome_WT_Bap1KO_QCV2vC1.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(wnnUMAP_1, wnnUMAP_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DEG_count %>% dplyr::rename( "cluster"="Cell_Name"))
## Create a UMAP plot colored by DEG values, with cluster DEG counts as text annotations
pdf("output/Signac/FeaturePlot_WNN_multiome_WT_Bap1KO_QCV2vC1_DEG_numeric.pdf", width=6, height=6)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "DEG", pt.size = 0.5, reduction = "wnn.umap") +
  scale_colour_viridis(option="mako") + # 
  geom_text(data = cluster_centers, aes(x = wnnUMAP_1, y = wnnUMAP_2, label = Num_DEGs), 
            size = 5, color = "red", fontface = "bold") 
dev.off()

# RNA UMAP REDUCTION ############ 
umap_coordinates <- as.data.frame(multiome_WT_Bap1KO_QCV2vC1.sct@reductions$umap@cell.embeddings)
umap_coordinates$cluster <- multiome_WT_Bap1KO_QCV2vC1.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DEG_count %>% dplyr::rename( "cluster"="Cell_Name"))
## Create a UMAP plot colored by DEG values, with cluster DEG counts as text annotations
pdf("output/Signac/FeaturePlot_RNAUMAP_multiome_WT_Bap1KO_QCV2vC1_DEGpadj05fc05_numeric_correct1.pdf", width=6, height=6)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "DEG", pt.size = 0.5, reduction = "umap") +
  scale_colour_viridis(option="mako") + # 
  geom_text(data = cluster_centers, aes(x = umap_1, y = umap_2, label = Num_DEGs), 
            size = 5, color = "red", fontface = "bold") 
dev.off()







# DARs number colored in a UMAP - raw DARs
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "cluster.annot"

DAR_peaks_count <- read_tsv("output/Signac/DAR_peaks_QCV2vC1_dim40kparam42res065algo4feat2000.txt") %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DARs=n())
# Add DAR information to my seurat object - DAR_peaks_count
cell_clusters <- multiome_WT_Bap1KO_QCV2vC1.sct@meta.data$cluster.annot
names(cell_clusters) <- rownames(multiome_WT_Bap1KO_QCV2vC1.sct@meta.data)
DAR_named_vector <- DAR_peaks_count$Num_DARs[match(cell_clusters, DAR_peaks_count$cluster)]
names(DAR_named_vector) <- names(cell_clusters)
# Integrate DAR values into the Seurat object
multiome_WT_Bap1KO_QCV2vC1.sct <- AddMetaData(multiome_WT_Bap1KO_QCV2vC1.sct, metadata = DAR_named_vector, col.name = "DAR")
# Create a UMAP plot colored by qval values
pdf("output/Signac/FeaturePlot_WNN_multiome_WT_Bap1KO_QCV2vC1_DAR.pdf", width=6, height=6)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "DAR", pt.size = 0.5, reduction = "wnn.umap") +
  scale_colour_viridis(option="magma") # 
dev.off()

# WNN REDUCTION ############
# Add values on the heatmap
## Extract UMAP coordinates
umap_coordinates <- as.data.frame(multiome_WT_Bap1KO_QCV2vC1.sct@reductions$wnn.umap@cell.embeddings)
umap_coordinates$cluster <- multiome_WT_Bap1KO_QCV2vC1.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(wnnUMAP_1, wnnUMAP_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DAR_peaks_count)
## Create a UMAP plot colored by DAR values, with cluster DEG counts as text annotations
pdf("output/Signac/FeaturePlot_WNN_multiome_WT_Bap1KO_QCV2vC1_DAR_numeric.pdf", width=6, height=6)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "DAR", pt.size = 0.5, reduction = "wnn.umap") +
  scale_colour_viridis(option="mako") + # 
  geom_text(data = cluster_centers, aes(x = wnnUMAP_1, y = wnnUMAP_2, label = Num_DARs), 
            size = 5, color = "red", fontface = "bold") 
dev.off()

# RNA UMAP REDUCTION ############

# Add values on the heatmap
## Extract UMAP coordinates
umap_coordinates <- as.data.frame(multiome_WT_Bap1KO_QCV2vC1.sct@reductions$umap@cell.embeddings)
umap_coordinates$cluster <- multiome_WT_Bap1KO_QCV2vC1.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DAR_peaks_count)
## Create a UMAP plot colored by DAR values, with cluster DEG counts as text annotations
pdf("output/Signac/FeaturePlot_RNAUMAP_multiome_WT_Bap1KO_QCV2vC1_DAR_numeric.pdf", width=6, height=6)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "DAR", pt.size = 0.5, reduction = "umap") +
  scale_colour_viridis(option="mako") + # 
  geom_text(data = cluster_centers, aes(x = umap_1, y = umap_2, label = Num_DARs), 
            size = 5, color = "red", fontface = "bold") 
dev.off()





# DARs number colored in a UMAP - gene with a DARs
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "cluster.annot"

DAR_gene_count <- read_tsv("output/Signac/DAR_genes_QCV2vC1_dim40kparam42res065algo4feat2000.txt") %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select(gene_name, cluster) %>%
  unique() %>%
  group_by(cluster) %>%
  summarise(Num_DARs=n())

# WNN REDUCTION ############
# Add DAR information to my seurat object - DAR_gene_count
cell_clusters <- multiome_WT_Bap1KO_QCV2vC1.sct@meta.data$cluster.annot
names(cell_clusters) <- rownames(multiome_WT_Bap1KO_QCV2vC1.sct@meta.data)
DAR_named_vector <- DAR_gene_count$Num_DARs[match(cell_clusters, DAR_gene_count$cluster)]
names(DAR_named_vector) <- names(cell_clusters)
# Integrate DAR values into the Seurat object
multiome_WT_Bap1KO_QCV2vC1.sct <- AddMetaData(multiome_WT_Bap1KO_QCV2vC1.sct, metadata = DAR_named_vector, col.name = "DARgene")
# Create a UMAP plot colored by qval values
pdf("output/Signac/FeaturePlot_WNN_multiome_WT_Bap1KO_QCV2vC1_DARgene.pdf", width=6, height=6)
FeaturePlot(multiome_WT_Bap1KO_QCV2vC1.sct, features = "DARgene", pt.size = 0.5, reduction = "wnn.umap") +
  scale_colour_viridis(option="magma") # 
dev.off()



# Count nb of up down reg genes

cluster_types <- c("cluster1",
"cluster2",
"cluster3",
"cluster4",
"cluster5",
"cluster6",
"cluster7",
"cluster8",
"cluster9",
"cluster10",
"cluster11",
"cluster12",
"cluster13",
"cluster14",
"cluster15",
"cluster16",
"cluster17",
"cluster18",
"cluster19")
## Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/Signac/", cluster, "-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}
## Define cluster names and corresponding data frames
clusters <- list(
  "PyNs_SubC_CA23" = cluster1,
  "IN_1" = cluster2,
  "SubC_1" = cluster3,
  "PyNs_SubC_CA1" = cluster4,
  "PyNs_RSC_UL" = cluster5,
  "DG_GC" = cluster6,
  "PyNs_RSC_MDL" = cluster7,
  "NSC_proliferative_1" = cluster8,
  "SubC_2" = cluster9,
  "IN_2" = cluster10,
  "NSC_quiescent" = cluster11,
  "IN_SubC" = cluster12,
  "IP" = cluster13,
  "NSC_proliferative_2" = cluster14,
  "CR" = cluster15,
  "OPC" = cluster16,
  "Meningeal_Cells" = cluster17,
  "Radial_Glia_Cells" = cluster18,
  "Microglia" = cluster19
)

## Function to count up- and down-regulated genes in each cluster
count_up_down <- function(cluster_data) {
  significant_genes <- cluster_data %>% dplyr::filter(p_val_adj < 0.05)
  nb_upregulated <- sum(significant_genes$avg_log2FC > 0)
  nb_downregulated <- sum(significant_genes$avg_log2FC < 0)
  
  return(data.frame(nb_upregulated = nb_upregulated, nb_downregulated = nb_downregulated))
}
results <- lapply(clusters, count_up_down)
final_results <- do.call(rbind, results) %>%
  as_tibble(rownames = "cluster") %>%
  arrange(cluster)
row.names(final_results) <- names(clusters)
print(final_results)





# proportion of DEG genes also DAR - Overlap DAR DEG same cell type
## import DAR genes
DAR_genes = read_tsv("output/Signac/DAR_genes_QCV2vC1_dim40kparam42res065algo4feat2000.txt") %>%
  dplyr::select(gene_name, p_val_adj, avg_log2FC, cluster)
DAR_genes_signif = DAR_genes %>%
  dplyr::filter(p_val_adj <0.05) %>%
  dplyr::rename("gene" = "gene_name")
## Function to count up- and down-DAR genes in each cluster
count_up_down <- function(data) {
  significant_genes <- data %>% dplyr::filter(p_val_adj < 0.05)
  results <- significant_genes %>%
    group_by(cluster) %>%
    summarize(
      nb_up = sum(avg_log2FC > 0),
      nb_down = sum(avg_log2FC < 0)
    )
  return(results)
}
final_results <- count_up_down(DAR_genes_signif)
print(final_results)


## import DEG genes
cluster_types <- c("cluster1",
"cluster2",
"cluster3",
"cluster4",
"cluster5",
"cluster6",
"cluster7",
"cluster8",
"cluster9",
"cluster10",
"cluster11",
"cluster12",
"cluster13",
"cluster14",
"cluster15",
"cluster16",
"cluster17",
"cluster18",
"cluster19")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/Signac/", cluster, "-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}

PyNs_SubC_CA23 = cluster1  %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "PyNs_SubC_CA23")
IN_1 = cluster2 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "IN_1")
SubC_1 = cluster3 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "SubC_1")
PyNs_SubC_CA1 = cluster4 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "PyNs_SubC_CA1")
PyNs_RSC_UL = cluster5 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble()  %>%
  add_column(cluster= "PyNs_RSC_UL")
DG_GC = cluster6 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "DG_GC")
PyNs_RSC_MDL = cluster7 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "PyNs_RSC_MDL")
NSC_proliferative_1 = cluster8 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "NSC_proliferative_1")
SubC_2 = cluster9 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "SubC_2")
IN_2 = cluster10 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "IN_2")
NSC_quiescent = cluster11 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "NSC_quiescent")
IN_SubC = cluster12 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "IN_SubC")
IP = cluster13 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "IP")
NSC_proliferative_2 = cluster14 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "NSC_proliferative_2")
CR = cluster15 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "CR")
OPC = cluster16 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "OPC")
Meningeal_Cells = cluster17 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "Meningeal_Cells")
Radial_Glia_Cells = cluster18 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "Radial_Glia_Cells")
Microglia = cluster19 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "Microglia")


DEG_genes = PyNs_SubC_CA23 %>%
  bind_rows(IN_1) %>%
  bind_rows(SubC_1) %>%
  bind_rows(PyNs_SubC_CA1) %>%
  bind_rows(PyNs_RSC_UL) %>%
  bind_rows(DG_GC) %>%
  bind_rows(PyNs_RSC_MDL) %>%
  bind_rows(NSC_proliferative_1) %>%
  bind_rows(SubC_2) %>%
  bind_rows(IN_2) %>%
  bind_rows(NSC_quiescent) %>%
  bind_rows(IN_SubC) %>%
  bind_rows(IP) %>%
  bind_rows(NSC_proliferative_2) %>%
  bind_rows(CR) %>%
  bind_rows(OPC) %>%
  bind_rows(Meningeal_Cells) %>%
  bind_rows(Radial_Glia_Cells) %>%
  bind_rows(Microglia) 


DEG_genes_signif = DEG_genes %>%
  dplyr::filter(p_val_adj < 0.05 ) # & abs(avg_log2FC) > 0.5
   
   
# plot
DEG_genes_signif
DAR_genes_signif
   
### output save table DEG DAR ###
DEG_genes_signif %>%
  dplyr::rename("avg_log2FC_DEG" = "avg_log2FC", "p_val_adj_DEG" = "p_val_adj") %>%
  dplyr::select(gene, avg_log2FC_DEG,p_val_adj_DEG) %>%
  left_join(DAR_genes_signif) %>%
  dplyr::rename("avg_log2FC_DAR" = "avg_log2FC", "p_val_adj_DAR" = "p_val_adj") %>%
  dplyr::select(gene, avg_log2FC_DEG,p_val_adj_DEG, avg_log2FC_DAR, p_val_adj_DAR) %>%
  unique() %>%
  mutate(overlap_with_DAR = ifelse(gene %in% DAR_genes_signif$gene, "Overlap", "No Overlap")) %>%
  write_tsv(., "output/Signac/DEG_genes_signif_DAR.tsv")
###

DEG_DAR_genes_signif <- DEG_genes_signif %>%
  mutate(overlap_with_DAR = ifelse(gene %in% DAR_genes_signif$gene, "Overlap", "No Overlap"))
   
   
DEG_DAR_genes_signif$cluster <- factor(DEG_DAR_genes_signif$cluster, levels = c("NSC_quiescent",
"NSC_proliferative_1",
"NSC_proliferative_2",
"IP",
"Radial_Glia_Cells",
"OPC",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_MDL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"IN_SubC")) 

pdf("output/Signac/barplot_DEGpadj05fc05_DAR_prop.pdf", width=6, height=4)
DEG_DAR_genes_signif %>%
  group_by(cluster, overlap_with_DAR) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(proportion = count / sum(count) * 100) %>% # Calculate proportions
ggplot(., aes(x = cluster, y = proportion, fill = overlap_with_DAR)) +
  geom_bar(stat = "identity", position = "fill")  +
  geom_text(aes(label = count), 
            position = position_fill(vjust = 0.5), # Place text within the bar
            size = 3, color = "white", fontface = "bold") +
  labs(title = "Proportion of DEGs Overlapping with DARs by Cell Type",
       x = "Cell Type",
       y = "Proportion (%)",
       fill = "Overlap Status") +
  theme_bw() +
  scale_fill_manual(values = c("Overlap" = "blue", "No Overlap" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off() 
  
   
   
### only up regulated genes with DAR opening


DEG_genes_signif_pos = DEG_genes_signif %>% 
  dplyr::filter(avg_log2FC >0)

DAR_genes_signif_pos = DAR_genes_signif %>% 
  dplyr::filter(avg_log2FC >0)
   

DEG_DAR_genes_signif_pos <- DEG_genes_signif_pos %>%
  rowwise() %>%
  mutate(overlap_with_DAR = ifelse(
    any(gene == DAR_genes_signif_pos$gene & cluster == DAR_genes_signif_pos$cluster), 
    "Overlap", 
    "No Overlap"
  )) %>%
  ungroup()

### output save table DEG DAR ###
DEG_DAR_genes_signif_pos %>%
  unique() %>%
  write_tsv(., "output/Signac/DEG_DAR_genes_signif_pos.tsv")
###

   
DEG_DAR_genes_signif_pos$cluster <- factor(DEG_DAR_genes_signif_pos$cluster, levels = c("NSC_quiescent",
"NSC_proliferative_1",
"NSC_proliferative_2",
"IP",
"Radial_Glia_Cells",
"OPC",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_MDL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"IN_SubC")) 

pdf("output/Signac/barplot_DEGpadj05fc05_DAR_pos_prop.pdf", width=6, height=4)
DEG_DAR_genes_signif_pos %>%
  group_by(cluster, overlap_with_DAR) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(proportion = count / sum(count) * 100) %>% # Calculate proportions
ggplot(., aes(x = cluster, y = proportion, fill = overlap_with_DAR)) +
  geom_bar(stat = "identity", position = "fill")  +
  geom_text(aes(label = count), 
            position = position_fill(vjust = 0.5), # Place text within the bar
            size = 3, color = "white", fontface = "bold") +
  labs(title = "Proportion of upregulated DEGs Overlapping with open DARs by Cell Type",
       x = "Cell Type",
       y = "Proportion (%)",
       fill = "Overlap Status") +
  theme_bw() +
  scale_fill_manual(values = c("Overlap" = "blue", "No Overlap" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off() 
  

### only down regulated genes with DAR closing


DEG_genes_signif_neg = DEG_genes_signif %>% 
  dplyr::filter(avg_log2FC <0)

DAR_genes_signif_neg = DAR_genes_signif %>% 
  dplyr::filter(avg_log2FC <0)
   

DEG_DAR_genes_signif_neg <- DEG_genes_signif_neg %>%
  rowwise() %>%
  mutate(overlap_with_DAR = ifelse(
    any(gene == DAR_genes_signif_neg$gene & cluster == DAR_genes_signif_neg$cluster), 
    "Overlap", 
    "No Overlap"
  )) %>%
  ungroup()

### output save table DEG DAR ###
DEG_DAR_genes_signif_neg %>%
  unique() %>%
  write_tsv(., "output/Signac/DEG_DAR_genes_signif_neg.tsv")
###
   
DEG_DAR_genes_signif_neg$cluster <- factor(DEG_DAR_genes_signif_neg$cluster, levels = c("NSC_quiescent",
"NSC_proliferative_1",
"NSC_proliferative_2",
"IP",
"Radial_Glia_Cells",
"OPC",
"Microglia",
"Meningeal_Cells",
"CR",
"DG_GC",
"PyNs_SubC_CA1",
"PyNs_SubC_CA23",
"PyNs_RSC_UL",
"PyNs_RSC_MDL",
"SubC_1",
"SubC_2",
"IN_1",
"IN_2",
"IN_SubC")) 

pdf("output/Signac/barplot_DEGpadj05fc05_DAR_neg_prop.pdf", width=6, height=4)
DEG_DAR_genes_signif_neg %>%
  group_by(cluster, overlap_with_DAR) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(proportion = count / sum(count) * 100) %>% # Calculate proportions
ggplot(., aes(x = cluster, y = proportion, fill = overlap_with_DAR)) +
  geom_bar(stat = "identity", position = "fill")  +
  geom_text(aes(label = count), 
            position = position_fill(vjust = 0.5), # Place text within the bar
            size = 3, color = "white", fontface = "bold") +
  labs(title = "Proportion of downregulated DEGs Overlapping with closing DARs by Cell Type",
       x = "Cell Type",
       y = "Proportion (%)",
       fill = "Overlap Status") +
  theme_bw() +
  scale_fill_manual(values = c("Overlap" = "blue", "No Overlap" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off() 
  



```


**Usefull tutorial:**
- *Diff. access. identification*: https://stuartlab.org/signac/articles/pbmc_vignette#find-differentially-accessible-peaks-between-cell-types
- *Motif enrichment*: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac
- *Find peaks that are correlated with the expression of nearby genes*: https://stuartlab.org/signac/reference/linkpeaks and need to compute genomic information: https://stuartlab.org/signac/reference/regionstats


### Pathway analysis

Let's do some pathway analysis for *WNT signaling*. Pathway of interest collected from [msigdb](https://www.gsea-msigdb.org/gsea/msigdb/mouse/genesets.jsp?collection=M2):
**List1:**
  - BIOCARTA_WNT_PATHWAY
  - REACTOME_SIGNALING_BY_WNT 
  - REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING 
  - WP_WNT_SIGNALING 
  - WP_WNT_SIGNALING_PATHWAY
  - WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY
**List2:**
  - GOBP_WNT_SIGNALING_PATHWAY
  - GOBP_CANONICAL_WNT_SIGNALING_PATHWAY
  - GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY
  - GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY
  - BIOCARTA_SHH_PATHWAY
  - allCombined_WNT
**List3:**
  - WP_HEDGEHOG_SIGNALING_PATHWAY
  - REACTOME_HEDGEHOG_ON_STATE
  - REACTOME_HEDGEHOG_OFF_STATE
  - REACTOME_HEDGEHOG_LIGAND_BIOGENESIS
  - REACTOME_SIGNALING_BY_HEDGEHOG
  - allCombined_SHH (154 genes)

Gene names downloaded and available at `output/Pathway/geneList_[PATHWAY].txt`.

```bash
conda activate SignacV5
module load hdf5
```

```R

set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()
# remotes::install_github('immunogenomics/presto')

# import seurat object
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")


# GSEA plot
library("fgsea")
#### import all cluster DEGs output:
cluster_types <- c("cluster1",
"cluster2",
"cluster3",
"cluster4",
"cluster5",
"cluster6",
"cluster7",
"cluster8",
"cluster9",
"cluster10",
"cluster11",
"cluster12",
"cluster13",
"cluster14",
"cluster15",
"cluster16",
"cluster17",
"cluster18",
"cluster19")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/Signac/", cluster, "-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}
##############################################################################################
################################# NB OF UP DOWN GENES ###############################################
# Define cluster names and corresponding data frames
clusters <- list(
  "PyNs_SubC_CA23" = cluster1,
  "IN_1" = cluster2,
  "SubC_1" = cluster3,
  "PyNs_SubC_CA1" = cluster4,
  "PyNs_RSC_UL" = cluster5,
  "DG_GC" = cluster6,
  "PyNs_RSC_MDL" = cluster7,
  "NSC_proliferative_1" = cluster8,
  "SubC_2" = cluster9,
  "IN_2" = cluster10,
  "NSC_quiescent" = cluster11,
  "IN_SubC" = cluster12,
  "IP" = cluster13,
  "NSC_proliferative_2" = cluster14,
  "CR" = cluster15,
  "OPC" = cluster16,
  "Meningeal_Cells" = cluster17,
  "Radial_Glia_Cells" = cluster18,
  "Microglia" = cluster19
)
# Function to count up- and down-regulated genes in each cluster
count_up_down <- function(cluster_data) {
  significant_genes <- cluster_data %>% dplyr::filter(p_val_adj < 0.05)
  nb_upregulated <- sum(significant_genes$avg_log2FC > 0)
  nb_downregulated <- sum(significant_genes$avg_log2FC < 0)
  return(data.frame(nb_upregulated = nb_upregulated, nb_downregulated = nb_downregulated))
}
# Iterate through each cluster and get the count of up- and down-regulated genes
results <- lapply(clusters, count_up_down)
# Combine results into a data frame
final_results <- do.call(rbind, results)
# Add cluster names as row names
row.names(final_results) <- names(clusters)
# Output the final results
print(final_results)
##############################################################################################
##############################################################################################


# Rename cluster with cluster name
PyNs_SubC_CA23 = cluster1 
IN_1 = cluster2
SubC_1 = cluster3
PyNs_SubC_CA1 = cluster4
PyNs_RSC_UL = cluster5
DG_GC = cluster6
PyNs_RSC_MDL = cluster7
NSC_proliferative_1 = cluster8
SubC_2 = cluster9
IN_2 = cluster10
NSC_quiescent = cluster11
IN_SubC = cluster12
IP = cluster13
NSC_proliferative_2 = cluster14
CR = cluster15
OPC = cluster16
Meningeal_Cells = cluster17
Radial_Glia_Cells = cluster18
Microglia = cluster19


# GSEA plot one by one
## load list of genes to test
### List1
fgsea_sets <- list(
  BIOCARTA_WNT_PATHWAY = read_table(file = c("output/Pathway/geneList_BIOCARTA_WNT_PATHWAY.txt"))$Genes,
  REACTOME_SIGNALING_BY_WNT = read_table(file = c("output/Pathway/geneList_REACTOME_SIGNALING_BY_WNT.txt"))$Genes,
  REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING = read_table(file = c("output/Pathway/geneList_REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING.txt"))$Genes,
  WP_WNT_SIGNALING = read_table(file = c("output/Pathway/geneList_WP_WNT_SIGNALING.txt"))$Genes,
  WP_WNT_SIGNALING_PATHWAY = read_table(file = c("output/Pathway/geneList_WP_WNT_SIGNALING_PATHWAY.txt"))$Genes,
  WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY = read_table(file = c("output/Pathway/geneList_WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY.txt"))$Genes
)

### List2
fgsea_sets <- list(
  GOBP_WNT_SIGNALING_PATHWAY = read_table(file = c("output/Pathway/geneList_GOBP_WNT_SIGNALING_PATHWAY.txt"))$Genes,
  GOBP_CANONICAL_WNT_SIGNALING_PATHWAY = read_table(file = c("output/Pathway/geneList_GOBP_CANONICAL_WNT_SIGNALING_PATHWAY.txt"))$Genes,
  GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY = read_table(file = c("output/Pathway/geneList_GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY.txt"))$Genes,
  GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY = read_table(file = c("output/Pathway/geneList_GOBP_REGULATION_OF_WNT_SIGNALING_PATHWAY.txt"))$Genes,
  BIOCARTA_SHH_PATHWAY = read_table(file = c("output/Pathway/geneList_BIOCARTA_SHH_PATHWAY.txt"))$Genes,
  allCombined_WNT = read_table(file = c("output/Pathway/geneList_allCombined_WNT.txt"))$Genes
)

### List3
fgsea_sets <- list(
  WP_HEDGEHOG_SIGNALING_PATHWAY = read_table(file = c("output/Pathway/geneList_WP_HEDGEHOG_SIGNALING_PATHWAY.txt"))$Genes,
  REACTOME_HEDGEHOG_ON_STATE = read_table(file = c("output/Pathway/geneList_REACTOME_HEDGEHOG_ON_STATE.txt"))$Genes,
  REACTOME_HEDGEHOG_OFF_STATE = read_table(file = c("output/Pathway/geneList_REACTOME_HEDGEHOG_OFF_STATE.txt"))$Genes,
  REACTOME_HEDGEHOG_LIGAND_BIOGENESIS = read_table(file = c("output/Pathway/geneList_REACTOME_HEDGEHOG_LIGAND_BIOGENESIS.txt"))$Genes,
  REACTOME_SIGNALING_BY_HEDGEHOG = read_table(file = c("output/Pathway/geneList_REACTOME_SIGNALING_BY_HEDGEHOG.txt"))$Genes,
  allCombined_SHH = read_table(file = c("output/Pathway/geneList_allCombined_SHH.txt"))$Genes
)





# FC only ##############################################
# GSEA all together
## List of clusters to analyze
clusters <- list("PyNs_SubC_CA23", "IN_1", "SubC_1", "PyNs_SubC_CA1", "PyNs_RSC_UL", "DG_GC", "PyNs_RSC_MDL",
                 "NSC_proliferative_1", "SubC_2", "IN_2", "NSC_quiescent", "IN_SubC", "IP", 
                 "NSC_proliferative_2", "CR", "OPC", "Meningeal_Cells", "Radial_Glia_Cells", "Microglia")

## Automate the GSEA and plot generation
for (cluster_name in clusters) {
  # Rank genes based on avg_log2FC for the current cluster
  genes <- get(cluster_name) %>%   # Ensure you have a list/dataframe for each cluster, e.g., PyNs_SubC_CA23, etc.
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC)
  ranks <- deframe(genes)


  # Run GSEA for each pathway
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 100000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(ES))
  # Save results
  write_tsv(fgseaResTidy, paste0("output/Pathway/GSEA_Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1_", cluster_name, ".tsv"))
  # Generate plots for each pathway
  for (pathway in names(fgsea_sets)) {
    # Create the plot
    gsea_plot <- plotEnrichment(fgsea_sets[[pathway]], ranks) +
      labs(title = paste0(pathway, "-", cluster_name)) +
      theme_bw()
    
    # Save the plot using ggsave
    ggsave(
      filename = paste0("output/Pathway/GSEA_Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1_", pathway, "_", cluster_name, ".pdf"),
      plot = gsea_plot,
      width = 5, height = 3, device = "pdf"
    )
  }
}

 

#--> 100000 time is OK;  1000000 = still warning; tested 10000000; too long!



# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c( "PyNs_SubC_CA23" ,
"IN_1" ,
"SubC_1" ,
"PyNs_SubC_CA1" ,
"PyNs_RSC_UL" ,
"DG_GC" ,
"PyNs_RSC_MDL" ,
"NSC_proliferative_1" ,
"SubC_2" ,
"IN_2" ,
"NSC_quiescent" ,
"IN_SubC" ,
"IP" ,
"NSC_proliferative_2" ,
"CR" ,
"OPC" ,
"Meningeal_Cells" ,
"Radial_Glia_Cells" ,
"Microglia" )

## Initialize an empty list to store the results for each cluster type
all_results <- list()
## Loop over each cluster type
for (cluster in cluster_types) {
  
  # Extract genes for the current cluster
  genes <- get(cluster) %>% 
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC)
  
  ranks <- deframe(genes)
  
  # Run GSEA for the current cluster
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 100000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(ES))
  
  # Extract summary table and add cluster column
  fgseaResTidy_summary = fgseaResTidy %>% 
    dplyr::select(pathway, pval, padj, ES, size, NES) %>%
    mutate(cluster = cluster) %>%
    arrange(padj) %>% 
    head()
  
  # Store results in the list
  all_results[[cluster]] <- fgseaResTidy_summary
}
## Combine results from all cluster types into one table
final_results <- bind_rows(all_results, .id = "cluster")



write.table(final_results, file = c("output/Pathway/gsea_output_Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1-List3.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap all GSEA
pdf("output/Pathway/heatmap_gsea_padj-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1-List3.pdf", width=8, height=4)
ggplot(final_results, aes(x=cluster, y=pathway, fill=NES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Norm. Enrichment\nScore") +
  geom_text(aes(label=sprintf("%.2f", NES)), 
            color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()







# rank_FCpval: sign(log2(FC))*-log(p value) ##############################################
# GSEA all together
## List of clusters to analyze
clusters <- list("PyNs_SubC_CA23", "IN_1", "SubC_1", "PyNs_SubC_CA1", "PyNs_RSC_UL", "DG_GC", "PyNs_RSC_MDL",
                 "NSC_proliferative_1", "SubC_2", "IN_2", "NSC_quiescent", "IN_SubC", "IP", 
                 "NSC_proliferative_2", "CR", "OPC", "Meningeal_Cells", "Radial_Glia_Cells", "Microglia")

## Automate the GSEA and plot generation
for (cluster_name in clusters) {
  # Rank genes based on avg_log2FC for the current cluster
  genes <- get(cluster_name) %>%   # Ensure you have a list/dataframe for each cluster, e.g., PyNs_SubC_CA23, etc.
    rownames_to_column(var = "gene") %>%
    mutate(rank_FCpval = sign(avg_log2FC) * -log(p_val)) %>%  # Create rank column
    arrange(desc(rank_FCpval)) %>%
    dplyr::select(gene, rank_FCpval)
  ranks <- deframe(genes)


  # Run GSEA for each pathway
  fgseaRes <- fgsea(fgsea_sets, stats = ranks)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(ES))
  # Save results
  write_tsv(fgseaResTidy, paste0("output/Pathway/GSEA_Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_rank_FCpval_correct1_", cluster_name, ".tsv"))
  # Generate plots for each pathway
  for (pathway in names(fgsea_sets)) {
    # Create the plot
    gsea_plot <- plotEnrichment(fgsea_sets[[pathway]], ranks) +
      labs(title = paste0(pathway, "-", cluster_name)) +
      theme_bw()
    
    # Save the plot using ggsave
    ggsave(
      filename = paste0("output/Pathway/GSEA_Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_rank_FCpval_correct1_", pathway, "_", cluster_name, ".pdf"),
      plot = gsea_plot,
      width = 5, height = 3, device = "pdf"
    )
  }
}




# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c( "PyNs_SubC_CA23" ,
"IN_1" ,
"SubC_1" ,
"PyNs_SubC_CA1" ,
"PyNs_RSC_UL" ,
"DG_GC" ,
"PyNs_RSC_MDL" ,
"NSC_proliferative_1" ,
"SubC_2" ,
"IN_2" ,
"NSC_quiescent" ,
"IN_SubC" ,
"IP" ,
"NSC_proliferative_2" ,
"CR" ,
"OPC" ,
"Meningeal_Cells" ,
"Radial_Glia_Cells" ,
"Microglia" )

## Initialize an empty list to store the results for each cluster type
all_results <- list()
## Loop over each cluster type
for (cluster in cluster_types) {
  
  # Extract genes for the current cluster
  genes <- get(cluster) %>% 
    rownames_to_column(var = "gene")  %>%
    mutate(rank_FCpval = sign(avg_log2FC) * -log(p_val)) %>%  # Create rank column
    arrange(desc(rank_FCpval)) %>%
    dplyr::select(gene, rank_FCpval)
  
  ranks <- deframe(genes)
  
  # Run GSEA for the current cluster
  fgseaRes <- fgsea(fgsea_sets, stats = ranks)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(ES))
  
  # Extract summary table and add cluster column
  fgseaResTidy_summary = fgseaResTidy %>% 
    dplyr::select(pathway, pval, padj, ES, size, NES) %>%
    mutate(cluster = cluster) %>%
    arrange(padj) %>% 
    head()
  
  # Store results in the list
  all_results[[cluster]] <- fgseaResTidy_summary
}
## Combine results from all cluster types into one table
final_results <- bind_rows(all_results, .id = "cluster")

write.table(final_results, file = c("output/Pathway/gsea_output_Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_rank_FCpval_correct1-List3.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap all GSEA
pdf("output/Pathway/heatmap_gsea_padj-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_rank_FCpval_correct1-List3.pdf", width=8, height=4)
ggplot(final_results, aes(x=cluster, y=pathway, fill=NES)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0, name="Norm. Enrichment\nScore") +
  geom_text(aes(label=sprintf("%.2f", NES)), 
            color = ifelse(final_results$padj <= 0.05, "black", "grey50"),  # change btween pvalue, qvalue,p.adjust
            size=2) +
  coord_fixed()  # Force aspect ratio of the plot to be 1:1
dev.off()



```

- *NOTE: if the fgsea list contain too many samples, it fail and missed some caclulation!*



### RNA and ATAC

I check whether RNA DEG are associated with DAR; not so much; only 20-25% are, in the right direction (eg. induced RNA with open chromatin). Let's see whether looking at raw RNA and ATAC level improve this. Let's try to generate heatmap of RNA ATAC on the DEG genes




```bash
conda activate SignacV5
module load hdf5
```

```R
set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
library("ggbeeswarm")
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()
# remotes::install_github('immunogenomics/presto')

# import seurat object
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")

### Test with scaling data ###
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA" # According to condiments workflow
multiome_WT_Bap1KO_QCV2vC1.sct <- NormalizeData(multiome_WT_Bap1KO_QCV2vC1.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(multiome_WT_Bap1KO_QCV2vC1.sct)
multiome_WT_Bap1KO_QCV2vC1.sct <- ScaleData(multiome_WT_Bap1KO_QCV2vC1.sct, features = all.genes) # zero-centres and scales it
#########

# import DEGs
cluster_types <- c("cluster1",
"cluster2",
"cluster3",
"cluster4",
"cluster5",
"cluster6",
"cluster7",
"cluster8",
"cluster9",
"cluster10",
"cluster11",
"cluster12",
"cluster13",
"cluster14",
"cluster15",
"cluster16",
"cluster17",
"cluster18",
"cluster19")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/Signac/", cluster, "-Bap1KO_response_multiome_QCV2vC1_dim40kparam42res065algo4feat2000_allGenes_correct1.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}

PyNs_SubC_CA23 = cluster1  %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "PyNs_SubC_CA23")
IN_1 = cluster2 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "IN_1")
SubC_1 = cluster3 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "SubC_1")
PyNs_SubC_CA1 = cluster4 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "PyNs_SubC_CA1")
PyNs_RSC_UL = cluster5 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble()  %>%
  add_column(cluster= "PyNs_RSC_UL")
DG_GC = cluster6 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "DG_GC")
PyNs_RSC_MDL = cluster7 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "PyNs_RSC_MDL")
NSC_proliferative_1 = cluster8 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "NSC_proliferative_1")
SubC_2 = cluster9 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "SubC_2")
IN_2 = cluster10 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "IN_2")
NSC_quiescent = cluster11 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "NSC_quiescent")
IN_SubC = cluster12 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "IN_SubC")
IP = cluster13 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "IP")
NSC_proliferative_2 = cluster14 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "NSC_proliferative_2")
CR = cluster15 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "CR")
OPC = cluster16 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "OPC")
Meningeal_Cells = cluster17 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "Meningeal_Cells")
Radial_Glia_Cells = cluster18 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "Radial_Glia_Cells")
Microglia = cluster19 %>% 
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  add_column(cluster= "Microglia")

DEG_genes = PyNs_SubC_CA23 %>%
  bind_rows(IN_1) %>%
  bind_rows(SubC_1) %>%
  bind_rows(PyNs_SubC_CA1) %>%
  bind_rows(PyNs_RSC_UL) %>%
  bind_rows(DG_GC) %>%
  bind_rows(PyNs_RSC_MDL) %>%
  bind_rows(NSC_proliferative_1) %>%
  bind_rows(SubC_2) %>%
  bind_rows(IN_2) %>%
  bind_rows(NSC_quiescent) %>%
  bind_rows(IN_SubC) %>%
  bind_rows(IP) %>%
  bind_rows(NSC_proliferative_2) %>%
  bind_rows(CR) %>%
  bind_rows(OPC) %>%
  bind_rows(Meningeal_Cells) %>%
  bind_rows(Radial_Glia_Cells) %>%
  bind_rows(Microglia) 

DEG_genes_signif = DEG_genes %>%
  dplyr::filter(p_val_adj < 0.05 ) # & abs(avg_log2FC) > 0.5

DEG_genes_signif_pos = DEG_genes_signif %>% 
  dplyr::filter(avg_log2FC >0)

DEG_genes_signif_pos__PyNs_SubC_CA23 = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "PyNs_SubC_CA23")
DEG_genes_signif_pos__IN_1 = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "IN_1")
DEG_genes_signif_pos__SubC_1 = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "SubC_1")
DEG_genes_signif_pos__PyNs_SubC_CA1 = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "PyNs_SubC_CA1")
DEG_genes_signif_pos__PyNs_RSC_UL = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "PyNs_RSC_UL")
DEG_genes_signif_pos__DG_GC = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "DG_GC")
DEG_genes_signif_pos__PyNs_RSC_MDL = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "PyNs_RSC_MDL")
DEG_genes_signif_pos__NSC_proliferative_1 = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "NSC_proliferative_1")
DEG_genes_signif_pos__SubC_2 = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "SubC_2")
DEG_genes_signif_pos__IN_2 = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "IN_2")
DEG_genes_signif_pos__NSC_quiescent = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "NSC_quiescent")
DEG_genes_signif_pos__IN_SubC = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "IN_SubC")
DEG_genes_signif_pos__IP = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "IP")
DEG_genes_signif_pos__NSC_proliferative_2 = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "NSC_proliferative_2")
DEG_genes_signif_pos__CR = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "CR")
DEG_genes_signif_pos__OPC = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "OPC")
DEG_genes_signif_pos__Meningeal_Cells = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "Meningeal_Cells")
DEG_genes_signif_pos__Radial_Glia_Cells = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "Radial_Glia_Cells")
DEG_genes_signif_pos__Microglia = DEG_genes_signif_pos %>%
  dplyr::filter(cluster == "Microglia")




DEG_genes_signif_neg = DEG_genes_signif %>% 
  dplyr::filter(avg_log2FC <0)

DEG_genes_signif_neg__PyNs_SubC_CA23 = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "PyNs_SubC_CA23")
DEG_genes_signif_neg__IN_1 = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "IN_1")
DEG_genes_signif_neg__SubC_1 = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "SubC_1")
DEG_genes_signif_neg__PyNs_SubC_CA1 = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "PyNs_SubC_CA1")
DEG_genes_signif_neg__PyNs_RSC_UL = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "PyNs_RSC_UL")
DEG_genes_signif_neg__DG_GC = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "DG_GC")
DEG_genes_signif_neg__PyNs_RSC_MDL = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "PyNs_RSC_MDL")
DEG_genes_signif_neg__NSC_proliferative_1 = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "NSC_proliferative_1")
DEG_genes_signif_neg__SubC_2 = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "SubC_2")
DEG_genes_signif_neg__IN_2 = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "IN_2")
DEG_genes_signif_neg__NSC_quiescent = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "NSC_quiescent")
DEG_genes_signif_neg__IN_SubC = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "IN_SubC")
DEG_genes_signif_neg__IP = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "IP")
DEG_genes_signif_neg__NSC_proliferative_2 = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "NSC_proliferative_2")
DEG_genes_signif_neg__CR = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "CR")
DEG_genes_signif_neg__OPC = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "OPC")
DEG_genes_signif_neg__Meningeal_Cells = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "Meningeal_Cells")
DEG_genes_signif_neg__Radial_Glia_Cells = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "Radial_Glia_Cells")
DEG_genes_signif_neg__Microglia = DEG_genes_signif_neg %>%
  dplyr::filter(cluster == "Microglia")







#### PLOT RNA

genes_of_interest <- DEG_genes_signif_pos$gene




genes_of_interest <- DEG_genes_signif_neg__Microglia$gene

genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(multiome_WT_Bap1KO_QCV2vC1.sct@assays$RNA@data)]
### extract WT and cYAPKO gene expression values from RNA assay
WT_expression <- multiome_WT_Bap1KO_QCV2vC1.sct@assays$RNA@scale.data[genes_of_interest, colnames(multiome_WT_Bap1KO_QCV2vC1.sct)[multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident == "multiome_WT"]]
Bap1KO_expression <- multiome_WT_Bap1KO_QCV2vC1.sct@assays$RNA@scale.data[genes_of_interest, colnames(multiome_WT_Bap1KO_QCV2vC1.sct)[multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident == "multiome_Bap1KO"]]
### mean expression values for each gene
WT_mean <- rowMeans(WT_expression)
Bap1KO_mean <- rowMeans(Bap1KO_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  WT = WT_mean,
  Bap1KO = Bap1KO_mean
) %>% 
pivot_longer(cols = c(WT, Bap1KO), names_to = "Condition", values_to = "Expression")
## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(WT_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)
data_for_plot$Condition <- factor(data_for_plot$Condition, levels = c("WT", "Bap1KO")) 

# pdf("output/Signac/boxplot_RNAexpression_DEG_genes_signif_pos_multiome_WT_Bap1KO_QCV2vC1.pdf", width=2, height=4)
# pdf("output/Signac/boxplot_RNAexpression_DEG_genes_signif_pos_multiome_WT_Bap1KO_QCV2vC1_scaleData.pdf", width=2, height=4)
# pdf("output/Signac/boxplot_RNAexpression_multiome_WT_Bap1KO_QCV2vC1-DEG_genes_signif_neg__DG_GC.pdf", width=2, height=4)

pdf("output/Signac/boxplot_RNAexpression_multiome_WT_Bap1KO_QCV2vC1_scaleData-DEG_genes_signif_neg__Microglia.pdf", width=2, height=4)
ggplot(data = data_for_plot, aes(x = Condition, y = Expression, fill = Condition)) +
  scale_fill_manual(values = c("WT" = "blue", "Bap1KO" = "red")) +
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.75), size = 1, color = NA) +
  geom_boxplot(notch = TRUE, outlier.size = -1, color = "black", lwd = 1, alpha = 0.7, show.legend = FALSE) +
  ggbeeswarm::geom_quasirandom(shape = 21, size = 2, dodge.width = 0.75, color = "black", alpha = 0.3, show.legend = FALSE) +
  theme_bw() +
  ylab("Expression") +
  xlab("Genotype") +
  theme(
    legend.position = "none", # Remove legend
    text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15, color = "black")
  )
dev.off()



pdf("output/Signac/heatmap_RNAexpression_DEG_genes_signif_pos_multiome_WT_Bap1KO_QCV2vC1.pdf", width=3, height=2)
ggplot(data_for_plot, aes(x=Gene, y=Condition, fill=Expression)) + 
  geom_tile(color = "black") +  # Add black contour to each tile
  theme_bw() +  # Use black-white theme for cleaner look
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  ) +
  scale_fill_viridis(direction = 1, option = "viridis", name="Expression") 
dev.off()



#### PLOT ATAC

genes_of_interest <- DEG_genes_signif_pos$gene   # CHANGE HERE neg pos




genes_of_interest <- DEG_genes_signif_neg__Microglia$gene

genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(multiome_WT_Bap1KO_QCV2vC1.sct@assays$GeneActivity@data)]
### extract WT and cYAPKO gene expression values from GeneActivity assay
WT_expression <- multiome_WT_Bap1KO_QCV2vC1.sct@assays$GeneActivity@scale.data[genes_of_interest, colnames(multiome_WT_Bap1KO_QCV2vC1.sct)[multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident == "multiome_WT"]]
Bap1KO_expression <- multiome_WT_Bap1KO_QCV2vC1.sct@assays$GeneActivity@scale.data[genes_of_interest, colnames(multiome_WT_Bap1KO_QCV2vC1.sct)[multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident == "multiome_Bap1KO"]]
### mean expression values for each gene
WT_mean <- rowMeans(WT_expression)
Bap1KO_mean <- rowMeans(Bap1KO_expression)
data_for_plot <- data.frame(
  Gene = genes_of_interest,
  WT = WT_mean,
  Bap1KO = Bap1KO_mean
) %>% 
pivot_longer(cols = c(WT, Bap1KO), names_to = "Condition", values_to = "Expression")
## ORder from low to high express
### Reordering the genes based on their mean expression in WT in ascending order
ordered_genes <- names(sort(WT_mean))
### Extracting unique gene names from the ordered list
unique_ordered_genes <- unique(ordered_genes)
### Filter out rows from data_for_plot that don't have their genes in unique_ordered_genes
data_for_plot <- data_for_plot[data_for_plot$Gene %in% unique_ordered_genes, ]
### Set the factor levels for 'Gene' according to the unique ordered list
data_for_plot$Gene <- factor(data_for_plot$Gene, levels = unique_ordered_genes)
data_for_plot$Condition <- factor(data_for_plot$Condition, levels = c("WT", "Bap1KO")) 

#pdf("output/Signac/boxplot_GeneActivity_DEG_genes_signif_neg_multiome_WT_Bap1KO_QCV2vC1.pdf", width=2, height=4) # CHANGE HERE neg pos
#pdf("output/Signac/boxplot_GeneActivity_multiome_WT_Bap1KO_QCV2vC1-DEG_genes_signif_neg__Microglia.pdf", width=2, height=4) # CHANGE HERE neg pos
#pdf("output/Signac/boxplot_GeneActivity_DEG_genes_signif_pos_multiome_WT_Bap1KO_QCV2vC1_scaleData.pdf", width=2, height=4) 
pdf("output/Signac/boxplot_GeneActivity_multiome_WT_Bap1KO_QCV2vC1_scaleData-DEG_genes_signif_neg__Microglia.pdf", width=2, height=4) 
ggplot(data = data_for_plot, aes(x = Condition, y = Expression, fill = Condition)) +
  scale_fill_manual(values = c("WT" = "blue", "Bap1KO" = "red")) +
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.75), size = 1, color = NA) +
  geom_boxplot(notch = TRUE, outlier.size = -1, color = "black", lwd = 1, alpha = 0.7, show.legend = FALSE) +
  ggbeeswarm::geom_quasirandom(shape = 21, size = 2, dodge.width = 0.75, color = "black", alpha = 0.3, show.legend = FALSE) +
  theme_bw() +
  ylab("Expression") +
  xlab("Genotype") +
  theme(
    legend.position = "none", # Remove legend
    text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15, color = "black")
  )
dev.off()





# Gene List:

DEG_genes_signif_pos__PyNs_SubC_CA23 
DEG_genes_signif_pos__IN_1
DEG_genes_signif_pos__SubC_1 
DEG_genes_signif_pos__PyNs_SubC_CA1 
DEG_genes_signif_pos__PyNs_RSC_UL 
DEG_genes_signif_pos__DG_GC 
DEG_genes_signif_pos__PyNs_RSC_MDL
DEG_genes_signif_pos__NSC_proliferative_1 
DEG_genes_signif_pos__SubC_2 
DEG_genes_signif_pos__IN_2 
DEG_genes_signif_pos__NSC_quiescent 
DEG_genes_signif_pos__IN_SubC
DEG_genes_signif_pos__IP 
DEG_genes_signif_pos__NSC_proliferative_2 
DEG_genes_signif_pos__CR 
DEG_genes_signif_pos__OPC 
DEG_genes_signif_pos__Meningeal_Cells 
DEG_genes_signif_pos__Radial_Glia_Cells
DEG_genes_signif_pos__Microglia 



DEG_genes_signif_neg__PyNs_SubC_CA23
DEG_genes_signif_neg__IN_1 
DEG_genes_signif_neg__SubC_1 
DEG_genes_signif_neg__PyNs_SubC_CA1 
DEG_genes_signif_neg__PyNs_RSC_UL
DEG_genes_signif_neg__DG_GC 
DEG_genes_signif_neg__PyNs_RSC_MDL 
DEG_genes_signif_neg__NSC_proliferative_1 
DEG_genes_signif_neg__SubC_2 
DEG_genes_signif_neg__IN_2 
DEG_genes_signif_neg__NSC_quiescent 
DEG_genes_signif_neg__IN_SubC 
DEG_genes_signif_neg__IP 
DEG_genes_signif_neg__NSC_proliferative_2 
DEG_genes_signif_neg__CR 
DEG_genes_signif_neg__OPC 
DEG_genes_signif_neg__Meningeal_Cells
DEG_genes_signif_neg__Radial_Glia_Cells
DEG_genes_signif_neg__Microglia 




# Light testing - Random
## Verification that marker genes are good
multiome_WT_Bap1KO_QCV2vC1.sct
Idents(multiome_WT_Bap1KO_QCV2vC1.sct) <- "cluster.annot"

all_markers <- FindAllMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA" # According to condiments workflow

multiome_WT_Bap1KO_QCV2vC1.sct <- NormalizeData(multiome_WT_Bap1KO_QCV2vC1.sct, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(multiome_WT_Bap1KO_QCV2vC1.sct)
multiome_WT_Bap1KO_QCV2vC1.sct <- ScaleData(multiome_WT_Bap1KO_QCV2vC1.sct, features = all.genes) # zero-centres and scales it
#########


all_markers <- FindAllMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#--> ALL GOODD!!!! But RNA needed to be logNorm and scale!

```


--> RNA is good; it seems I am using counts, non scale.data, but result are good.

--> ATAC is weird; always show more signal in Bap1; even for genes downregulated... I tried using scale.data; and it is even worst. Even more enriched toward Bap1KO.





### Pseudotime condiments


#### Install Signac with condiments


```bash
# try clone condiments_V6 and isntall Signac in R with `install.packages("Signac")`
conda create --name condiments_Signac --clone condiments_V6
#--> Fail
conda env remove --name condiments_Signac

# try clone Signac and isntall condiments in R with `BiocManager::install("condiments")`
conda create --name condiments_Signac --clone SignacV5
#--> package Cairo, ggrastr, scater, distinct, fail; try install them individually
module load cairo/1.17.4-GCCcore-11.3.0
#--> Re try `BiocManager::install("condiments")`, fail same; try Anaconda install Cairo:  conda install conda-forge::r-cairo
#--> Then re try `BiocManager::install("condiments")`
```

--> WORK!!!






#### Run condiments RNA assay - V1 common condition

```bash
conda activate condiments_Signac
```



```R
# package installation 
## install.packages("remotes")
## remotes::install_github("cran/spatstat.core")
## remotes::install_version("Seurat", "4.0.3")
## install.packages("magrittr")
## install.packages("magrittr")
## install.packages("dplyr")
## BiocManager::install("DelayedMatrixStats")
## BiocManager::install("tradeSeq")


# packages
library("condiments")
library("Seurat")
library("magrittr") # to use pipe
library("dplyr") # to use bind_cols and sample_frac
library("SingleCellExperiment") # for reducedDims
library("ggplot2")
library("slingshot")
library("DelayedMatrixStats")
library("tidyr")
library("tradeSeq") 
library("cowplot")
library("scales")
library("pheatmap")

# Data import

multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric_label.rds")


DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA" # According to condiments workflow


# convert to SingleCellExperiment
RNA_WT_Bap1KO <- as.SingleCellExperiment(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA")


# tidy
df <- bind_cols(
  as.data.frame(reducedDims(RNA_WT_Bap1KO)$UMAP),
  as.data.frame(colData(RNA_WT_Bap1KO)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
## genotype overlap
pdf("output/condiments/UMAP_treatment_multiome_WT_Bap1KO_QCV2vC1.pdf", width=6, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col = orig.ident)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("blue", "red")) + # Specify colors here
  labs(col = "Treatment") +
  theme_classic()
dev.off()

## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(umap_1, umap_2) %>% as.matrix(), 
  conditions = df$orig.ident,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

pdf("output/condiments/UMAP_imbalance_score_multiome_WT_Bap1KO_QCV2vC1.pdf", width=5, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores") +
  theme_classic()
dev.off()


#  Trajectory Inference and Differential Topology
set.seed(42)

## PLOT with Separate trajectories
### Testing Area ############
humangastruloid2472hrs <- slingshot(humangastruloid2472hrs, reducedDim = 'UMAP',
                 clusterLabels = colData(humangastruloid2472hrs)$seurat_clusters,
                 start.clus = '3', approx_points = 100)

humangastruloid2472hrs <- slingshot(humangastruloid2472hrs, reducedDim = 'UMAP',
                 clusterLabels = colData(humangastruloid2472hrs)$seurat_clusters,
                 start.clus = '3', end.clus = c("8","1","4") ,approx_points = 100)

#                 extend = 'n', stretch = 0)


##########################################
RNA_WT_Bap1KO <- slingshot(RNA_WT_Bap1KO, reducedDim = 'UMAP',
                 clusterLabels = colData(RNA_WT_Bap1KO)$cluster.annot,
                 start.clus = 'NSC_quiescent', end.clus = c("PyNs_RSC_UL","DG_GC") ,approx_points = 100, extend = 'n')




#test reduceDim PCA or subset endoderm

set.seed(42)
topologyTest(SlingshotDataSet(RNA_WT_Bap1KO), RNA_WT_Bap1KO$orig.ident) #  


sdss <- slingshot_conditions(SlingshotDataSet(RNA_WT_Bap1KO), RNA_WT_Bap1KO$orig.ident)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "orig.ident")

#  
pdf("output/condiments/UMAP_trajectory_separated_multiome_WT_Bap1KO_QCV2vC1_STARTNSCquiescentENDPyNsRSCULDGGCextendN.pdf", width=6, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col = orig.ident)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(orig.ident, Lineage, Order),
            aes(group = interaction(Lineage, orig.ident)), size = 1.5) +
  theme_classic()
dev.off()



# Add custom labels for each trajectory based on the Lineage
curves$label <- with(curves, ifelse(Lineage == 1, "Trajectory 1",
                               ifelse(Lineage == 2, "Trajectory 2",
                               ifelse(Lineage == 3, "Trajectory 3",
                               ifelse(Lineage == 4, "Trajectory 4",
                               ifelse(Lineage == 5, "Trajectory 5",
                               ifelse(Lineage == 6, "Trajectory 6",
                               ifelse(Lineage == 7, "Trajectory 7",
                               ifelse(Lineage == 8, "Trajectory 8", "Trajectory 9" )))))))))
                               

# 
pdf("output/condiments/UMAP_trajectory_separated_multiome_WT_Bap1KO_QCV2vC1_STARTNSCquiescentENDPyNsRSCULDGGCextendN_trajLabel.pdf", width=6, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col = orig.ident)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(orig.ident, Lineage, Order),
            aes(group = interaction(Lineage, orig.ident)), size = 1.5) +
  geom_text(data = curves %>% group_by(Lineage) %>% top_n(1, Order),
            aes(label = label, x = umap_1, y = umap_2, group = Lineage),
            size = 4, vjust = -1, hjust = 0.5) +
  theme_classic()
dev.off()


## PLOT with common trajectories - Individually
df_2 <- bind_cols(
  as.data.frame(reducedDim(RNA_WT_Bap1KO, "UMAP")),
  slingPseudotime(RNA_WT_Bap1KO) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(RNA_WT_Bap1KO) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         Lineage3_pst = if_else(is.na(Lineage3_pst), 0, Lineage3_pst),
         Lineage4_pst = if_else(is.na(Lineage4_pst), 0, Lineage4_pst),
         Lineage5_pst = if_else(is.na(Lineage5_pst), 0, Lineage5_pst),
         Lineage6_pst = if_else(is.na(Lineage6_pst), 0, Lineage6_pst),
         Lineage7_pst = if_else(is.na(Lineage7_pst), 0, Lineage7_pst),
         Lineage8_pst = if_else(is.na(Lineage8_pst), 0, Lineage8_pst),
         Lineage9_pst = if_else(is.na(Lineage9_pst), 0, Lineage9_pst))
curves <- slingCurves(RNA_WT_Bap1KO, as.df = TRUE)
### Function to create the plot for each lineage
create_plot <- function(lineage_number) {
  df_2 <- df_2 %>%
    mutate(pst = case_when(
      !!sym(paste0("Lineage", lineage_number, "_pst")) > 0 ~ !!sym(paste0("Lineage", lineage_number, "_pst")),
      TRUE ~ 0
    ),
    group = if_else(pst > 0, paste0("lineage", lineage_number), "other"))
  curves_filtered <- curves %>% filter(Lineage == lineage_number)
  curves_endpoints <- curves_filtered %>%
    group_by(Lineage) %>%
    arrange(Order) %>%
    top_n(1, Order) # Get the top/last ordered point for each group
  df_2_lineage <- df_2 %>% filter(group == paste0("lineage", lineage_number))
  df_2_other <- df_2 %>% filter(group == "other")
  p <- ggplot() +
    geom_point(data = df_2_other, aes(x = umap_1, y = umap_2), size = .7, color = "grey85") +
    geom_point(data = df_2_lineage, aes(x = umap_1, y = umap_2, col = pst), size = .7) +
    scale_color_viridis_c() +
    labs(col = "Pseudotime", title = paste("Lineage", lineage_number)) +
    geom_path(data = curves_filtered %>% arrange(Order),
              aes(x = umap_1, y = umap_2, group = Lineage), col = "black", size = 1) +
    geom_text(data = curves_endpoints, aes(x = umap_1, y = umap_2, label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
    theme_classic()
  return(p)
}
### Generate the plots for each lineage
plots <- list()
for (i in 1:9) {
  plots[[i]] <- create_plot(i)
}
pdf("output/condiments/UMAP_trajectory_common_label_multiome_WT_Bap1KO_QCV2vC1_STARTNSCquiescentENDPyNsRSCULDGGCextendN_Lineage123456789.pdf", width=25, height=10)
gridExtra::grid.arrange(grobs = plots, ncol = 5)
dev.off()



## PLOT with separate trajectories - Individually
### WT
RNA_WT <- RNA_WT_Bap1KO[, RNA_WT_Bap1KO$orig.ident == "multiome_WT"]

df_2 <- bind_cols(
  as.data.frame(reducedDim(RNA_WT, "UMAP")),
  slingPseudotime(RNA_WT) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(RNA_WT) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         Lineage3_pst = if_else(is.na(Lineage3_pst), 0, Lineage3_pst),
         Lineage4_pst = if_else(is.na(Lineage4_pst), 0, Lineage4_pst),
         Lineage5_pst = if_else(is.na(Lineage5_pst), 0, Lineage5_pst),
         Lineage6_pst = if_else(is.na(Lineage6_pst), 0, Lineage6_pst),
         Lineage7_pst = if_else(is.na(Lineage7_pst), 0, Lineage7_pst),
         Lineage8_pst = if_else(is.na(Lineage8_pst), 0, Lineage8_pst),
         Lineage9_pst = if_else(is.na(Lineage9_pst), 0, Lineage9_pst))
curves <- slingCurves(RNA_WT, as.df = TRUE)
### Function to create the plot for each lineage
create_plot <- function(lineage_number) {
  df_2 <- df_2 %>%
    mutate(pst = case_when(
      !!sym(paste0("Lineage", lineage_number, "_pst")) > 0 ~ !!sym(paste0("Lineage", lineage_number, "_pst")),
      TRUE ~ 0
    ),
    group = if_else(pst > 0, paste0("lineage", lineage_number), "other"))
  curves_filtered <- curves %>% filter(Lineage == lineage_number)
  curves_endpoints <- curves_filtered %>%
    group_by(Lineage) %>%
    arrange(Order) %>%
    top_n(1, Order) # Get the top/last ordered point for each group
  df_2_lineage <- df_2 %>% filter(group == paste0("lineage", lineage_number))
  df_2_other <- df_2 %>% filter(group == "other")
  p <- ggplot() +
    geom_point(data = df_2_other, aes(x = umap_1, y = umap_2), size = .7, color = "grey85") +
    geom_point(data = df_2_lineage, aes(x = umap_1, y = umap_2, col = pst), size = .7) +
    scale_color_viridis_c() +
    labs(col = "Pseudotime", title = paste("Lineage", lineage_number)) +
    geom_path(data = curves_filtered %>% arrange(Order),
              aes(x = umap_1, y = umap_2, group = Lineage), col = "black", size = 1) +
    geom_text(data = curves_endpoints, aes(x = umap_1, y = umap_2, label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
    theme_classic()
  return(p)
}
### Generate the plots for each lineage
plots <- list()
for (i in 1:9) {
  plots[[i]] <- create_plot(i)
}
pdf("output/condiments/UMAP_trajectory_common_label_multiome_WT_Bap1KO_QCV2vC1_STARTNSCquiescentENDPyNsRSCULDGGCextendN_Lineage123456789_WTonly.pdf", width=25, height=10)
gridExtra::grid.arrange(grobs = plots, ncol = 5)
dev.off()


### Bap1KO
RNA_Bap1KO <- RNA_WT_Bap1KO[, RNA_WT_Bap1KO$orig.ident == "multiome_Bap1KO"]

df_2 <- bind_cols(
  as.data.frame(reducedDim(RNA_Bap1KO, "UMAP")),
  slingPseudotime(RNA_Bap1KO) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(RNA_Bap1KO) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         Lineage3_pst = if_else(is.na(Lineage3_pst), 0, Lineage3_pst),
         Lineage4_pst = if_else(is.na(Lineage4_pst), 0, Lineage4_pst),
         Lineage5_pst = if_else(is.na(Lineage5_pst), 0, Lineage5_pst),
         Lineage6_pst = if_else(is.na(Lineage6_pst), 0, Lineage6_pst),
         Lineage7_pst = if_else(is.na(Lineage7_pst), 0, Lineage7_pst),
         Lineage8_pst = if_else(is.na(Lineage8_pst), 0, Lineage8_pst),
         Lineage9_pst = if_else(is.na(Lineage9_pst), 0, Lineage9_pst))
curves <- slingCurves(RNA_Bap1KO, as.df = TRUE)
### Function to create the plot for each lineage
create_plot <- function(lineage_number) {
  df_2 <- df_2 %>%
    mutate(pst = case_when(
      !!sym(paste0("Lineage", lineage_number, "_pst")) > 0 ~ !!sym(paste0("Lineage", lineage_number, "_pst")),
      TRUE ~ 0
    ),
    group = if_else(pst > 0, paste0("lineage", lineage_number), "other"))
  curves_filtered <- curves %>% filter(Lineage == lineage_number)
  curves_endpoints <- curves_filtered %>%
    group_by(Lineage) %>%
    arrange(Order) %>%
    top_n(1, Order) # Get the top/last ordered point for each group
  df_2_lineage <- df_2 %>% filter(group == paste0("lineage", lineage_number))
  df_2_other <- df_2 %>% filter(group == "other")
  p <- ggplot() +
    geom_point(data = df_2_other, aes(x = umap_1, y = umap_2), size = .7, color = "grey85") +
    geom_point(data = df_2_lineage, aes(x = umap_1, y = umap_2, col = pst), size = .7) +
    scale_color_viridis_c() +
    labs(col = "Pseudotime", title = paste("Lineage", lineage_number)) +
    geom_path(data = curves_filtered %>% arrange(Order),
              aes(x = umap_1, y = umap_2, group = Lineage), col = "black", size = 1) +
    geom_text(data = curves_endpoints, aes(x = umap_1, y = umap_2, label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
    theme_classic()
  return(p)
}
### Generate the plots for each lineage
plots <- list()
for (i in 1:9) {
  plots[[i]] <- create_plot(i)
}
pdf("output/condiments/UMAP_trajectory_common_label_multiome_WT_Bap1KO_QCV2vC1_STARTNSCquiescentENDPyNsRSCULDGGCextendN_Lineage123456789_Bap1KOonly.pdf", width=25, height=10)
gridExtra::grid.arrange(grobs = plots, ncol = 5)
dev.off()









# Differential Progression
progressionTest(RNA_WT_Bap1KO, conditions = RNA_WT_Bap1KO$orig.ident, lineages = TRUE)

prog_res <- progressionTest(RNA_WT_Bap1KO, conditions = RNA_WT_Bap1KO$orig.ident, lineages = TRUE)

df_3 <-  slingPseudotime(RNA_WT_Bap1KO) %>% as.data.frame() 

df_3$condition <- RNA_WT_Bap1KO$orig.ident
df_3 <- df_3 %>% 
  pivot_longer(-condition, names_to = "Lineage",
               values_to = "pst") %>%
  filter(!is.na(pst))

pdf("output/condiments/densityPlot_trajectory_lineages_multiome_WT_Bap1KO_QCV2vC1_STARTNSCquiescentENDPyNsRSCULDGGCextendN.pdf", width=15, height=5)
ggplot(df_3, aes(x = pst)) +
  geom_density(alpha = .8, aes(fill = condition), col = "transparent") +
  geom_density(aes(col = condition), fill = "transparent", size = 1.5) +
  labs(x = "Pseudotime", fill = "condition") +
  facet_wrap(~Lineage, scales = "free", nrow=2) +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("blue", "red"))
  )) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()
dev.off()


#### ->  save.image(file="output/condiments/condiments_multiome_WT_Bap1KO_QCV2vC1.RData")
### load("output/condiments/condiments_multiome_WT_Bap1KO_QCV2vC1.RData")
set.seed(42)

#  Differential expression
# --> Run fitGam() through Slurm



################### Time Course effect COMMON CONDITIONS ######################################################
## TRAJECTORY9 ##################
set.seed(42)
traj9_RNA_common <- readRDS("output/condiments/traj9_RNA_common.rds")


## DEGs between condition
traj9_RNA_common_l2fc0 <- conditionTest(traj9_RNA_common, l2fc = 0) # 


# Correct the pvalue with fdr
traj9_RNA_common_l2fc0$padj <- p.adjust(traj9_RNA_common_l2fc0$pvalue, "fdr")


### Save output tables
traj9_RNA_common_l2fc0$gene <- rownames(traj9_RNA_common_l2fc0) # create new column label gene; as matrix before
condRes_traj9_RNA_common_l2fc0 <- traj9_RNA_common_l2fc0[, c(ncol(traj9_RNA_common_l2fc0), 1:(ncol(traj9_RNA_common_l2fc0)-1))] # just to put gene column 1st
write.table(condRes_traj9_RNA_common_l2fc0, file = c("output/condiments/condRes_traj9_RNA_common_l2fc0.txt"),sep="\t", quote=FALSE, row.names=FALSE)

# Heatmap clutering DEGs per traj _ REVISED METHOD
## import DEGs
condRes_traj9_RNA_common_l2fc0 <- read.table("output/condiments/condRes_traj9_RNA_common_l2fc0.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) 

## Isolate significant DEGs and transform into a vector
conditionGenes_traj9_l2fc0_vector <- condRes_traj9_RNA_common_l2fc0 %>% 
  filter(padj <= 0.05) %>%
  pull(gene)

# Predict smoothed values
yhatSmooth <- 
  predictSmooth(traj9_RNA_common, gene = conditionGenes_traj9_l2fc0_vector, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth, 1, scales::rescale))
combinedData <- yhatSmoothScaled[, c(51:100, 1:50)]
# Generate heatmap with clustering
# Perform hierarchical clustering
hc <- hclust(dist(combinedData))
clusters <- cutree(hc, k=10) # !!!!!!!!!!!!!!!!!! CHANGE CLUSTER NB HERE !!!!!!!!!!!!!!!!!!
# Create an annotation data frame for the rows based on cluster assignments
annotation_row <- data.frame(Cluster = factor(clusters))



# Line plots
library("reshape2")
library("stringr")
# Assuming yhatSmoothScaled contains your smoothed gene expression data
# Convert the yhatSmoothScaled data to a dataframe
df <- as.data.frame(yhatSmoothScaled)
df$Gene <- rownames(df)
# Transform the data into a long format
df_long <- melt(df, id.vars = "Gene", variable.name = "Pseudotime", value.name = "Expression")
# Attach the cluster information to the data frame
df$Cluster <- factor(clusters[df$Gene])
df_long$Cluster <- df$Cluster[match(df_long$Gene, df$Gene)]

# Extract condition column
df_long$Condition <- ifelse(str_detect(df_long$Pseudotime, "WT"), "WT", "Bap1KO")

# Extract the point value and convert it to numeric
df_long$Updated_Pseudotime <- as.numeric(str_extract(df_long$Pseudotime, "(?<=point)\\d+"))

# Define colors for the conditions
color_map <- c("WT" = "blue", "Bap1KO" = "red")

gene_counts <- df_long %>%
  group_by(Cluster) %>%
  summarise(GeneCount = n_distinct(Gene))
df_long <- df_long %>%
  left_join(gene_counts, by = "Cluster") %>%
  mutate(ClusterLabel = paste0("Cluster ", Cluster, " (", GeneCount, " genes)"))

# Plot using ggplot
pdf("output/condiments/clustered_linePlot_traj9_RNA_common_l2fc0_cl10.pdf", width=10, height=5)
ggplot(df_long, aes(x = as.numeric(Updated_Pseudotime), y = Expression, group = Gene)) + 
  geom_line(data = subset(df_long, Condition == "WT"), aes(color = Condition), alpha = 0.5) +
  geom_line(data = subset(df_long, Condition == "Bap1KO"), aes(color = Condition), alpha = 0.5) +
  scale_color_manual(values = color_map) + 
  facet_wrap(~ClusterLabel, scales = "free_y", nrow = 2) +  # Use the updated ClusterLabel column
  theme_bw() +
  labs(title = "Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")

dev.off()

# Plot using ggplot
pdf("output/condiments/smoothed_linePlot_traj9_RNA_common_l2fc0_cl10.pdf", width=10, height=5)
ggplot(df_long, aes(x = Updated_Pseudotime, y = Expression, color = Condition)) + 
  geom_smooth(method = "loess", se = TRUE, span = 0.5) + 
  scale_color_manual(values = color_map) + 
  facet_wrap(~ClusterLabel, scales = "free_y", nrow = 2) +  # Use the updated ClusterLabel column
  theme_bw() +
  labs(title = "Smoothed Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")
dev.off()



### Export gene list from each cluster
## Create a data frame with gene names and their respective cluster assignments
output_df <- data.frame(
  gene = rownames(combinedData),
  cluster = clusters
)

# Write the data frame to a .txt file
write.table(output_df, 
            file = "output/condiments/gene_clusters_traj9_RNA_common_l2fc0_cl10.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

# Check some genes individually - RNA
## FOR LINEAGE 9
counts <- multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(RNA_WT_Bap1KO, na = FALSE) [,9] # HERE INDICATE TRAJ
cellweights <- slingCurveWeights(RNA_WT_Bap1KO) [,9] # HERE INDICATE TRAJ
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj9_RNA_common-Eomes.pdf", width=4, height=2)
plotSmoothers(traj9_RNA_common, sub_counts, gene = "Eomes", curvesCol = c("blue","red") ) +
scale_color_manual(values =c("blue","red"))
dev.off()


# Check some genes individually - ATAC
xxxxxx NOT WORKING xxxxxx
multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds")
## FOR LINEAGE 9
counts <- multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct[["GeneActivity"]]@counts # Collect the counts from seurat
cond <- factor(multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(RNA_WT_Bap1KO, na = FALSE) [,9] # HERE INDICATE TRAJ
cellweights <- slingCurveWeights(RNA_WT_Bap1KO) [,9] # HERE INDICATE TRAJ
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj9_ATAC_common-Prox1.pdf", width=4, height=2)
plotSmoothers(traj9_RNA_common, sub_counts, gene = "Prox1", curvesCol = c("blue","red") ) +
scale_color_manual(values =c("blue","red"))
dev.off()
xxxxxxxxxxxxxxxxxxxxxxxx





# Heatmap representation
# Define colors for each cluster
# 20
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen", "grey", "darkred", "darkblue", "gold", "darkgray", "lightblue", "lightgreen", "lightcoral", "lightpink", "lightcyan"))(20),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 10
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen" ))(10),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 8
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown" ))(8),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 7
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink" ))(7),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 6
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange"))(6),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 5
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple"))(5),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 4
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow" ))(4),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Generate the heatmap

col_order <- order(grepl("WT", colnames(combinedData)), decreasing = TRUE)
combinedData <- combinedData[, col_order]
pdf("output/condiments/heatmap_traj9_RNA_common_l2fc0_cl10.pdf", width=5, height=5)
pheatmap(combinedData,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  legend = TRUE,
  cutree_rows = 10,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
)
dev.off()






## TRAJECTORY2 ##################
set.seed(42)
traj2_RNA_common <- readRDS("output/condiments/traj2_RNA_common.rds")


## DEGs between condition
traj2_RNA_common_l2fc0 <- conditionTest(traj2_RNA_common, l2fc = 0) # 


# Correct the pvalue with fdr
traj2_RNA_common_l2fc0$padj <- p.adjust(traj2_RNA_common_l2fc0$pvalue, "fdr")


### Save output tables
traj2_RNA_common_l2fc0$gene <- rownames(traj2_RNA_common_l2fc0) # create new column label gene; as matrix before
condRes_traj2_RNA_common_l2fc0 <- traj2_RNA_common_l2fc0[, c(ncol(traj2_RNA_common_l2fc0), 1:(ncol(traj2_RNA_common_l2fc0)-1))] # just to put gene column 1st
write.table(condRes_traj2_RNA_common_l2fc0, file = c("output/condiments/condRes_traj2_RNA_common_l2fc0.txt"),sep="\t", quote=FALSE, row.names=FALSE)

# Heatmap clutering DEGs per traj _ REVISED METHOD
## import DEGs
condRes_traj2_RNA_common_l2fc0 <- read.table("output/condiments/condRes_traj2_RNA_common_l2fc0.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) 

## Isolate significant DEGs and transform into a vector
conditionGenes_traj2_l2fc0_vector <- condRes_traj2_RNA_common_l2fc0 %>% 
  filter(padj <= 0.05) %>%
  pull(gene)

# Predict smoothed values
yhatSmooth <- 
  predictSmooth(traj2_RNA_common, gene = conditionGenes_traj2_l2fc0_vector, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth, 1, scales::rescale))
combinedData <- yhatSmoothScaled[, c(51:100, 1:50)]
# Generate heatmap with clustering
# Perform hierarchical clustering
hc <- hclust(dist(combinedData))
clusters <- cutree(hc, k=10) # !!!!!!!!!!!!!!!!!! CHANGE CLUSTER NB HERE !!!!!!!!!!!!!!!!!!
# Create an annotation data frame for the rows based on cluster assignments
annotation_row <- data.frame(Cluster = factor(clusters))



# Line plots
library("reshape2")
library("stringr")
# Assuming yhatSmoothScaled contains your smoothed gene expression data
# Convert the yhatSmoothScaled data to a dataframe
df <- as.data.frame(yhatSmoothScaled)
df$Gene <- rownames(df)
# Transform the data into a long format
df_long <- melt(df, id.vars = "Gene", variable.name = "Pseudotime", value.name = "Expression")
# Attach the cluster information to the data frame
df$Cluster <- factor(clusters[df$Gene])
df_long$Cluster <- df$Cluster[match(df_long$Gene, df$Gene)]

# Extract condition column
df_long$Condition <- ifelse(str_detect(df_long$Pseudotime, "WT"), "WT", "Bap1KO")

# Extract the point value and convert it to numeric
df_long$Updated_Pseudotime <- as.numeric(str_extract(df_long$Pseudotime, "(?<=point)\\d+"))

# Define colors for the conditions
color_map <- c("WT" = "blue", "Bap1KO" = "red")

gene_counts <- df_long %>%
  group_by(Cluster) %>%
  summarise(GeneCount = n_distinct(Gene))
df_long <- df_long %>%
  left_join(gene_counts, by = "Cluster") %>%
  mutate(ClusterLabel = paste0("Cluster ", Cluster, " (", GeneCount, " genes)"))

# Plot using ggplot
pdf("output/condiments/clustered_linePlot_traj2_RNA_common_l2fc0_cl10.pdf", width=10, height=5)
ggplot(df_long, aes(x = as.numeric(Updated_Pseudotime), y = Expression, group = Gene)) + 
  geom_line(data = subset(df_long, Condition == "WT"), aes(color = Condition), alpha = 0.5) +
  geom_line(data = subset(df_long, Condition == "Bap1KO"), aes(color = Condition), alpha = 0.5) +
  scale_color_manual(values = color_map) + 
  facet_wrap(~ClusterLabel, scales = "free_y", nrow = 2) +  # Use the updated ClusterLabel column
  theme_bw() +
  labs(title = "Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")

dev.off()

# Plot using ggplot
pdf("output/condiments/smoothed_linePlot_traj2_RNA_common_l2fc0_cl10.pdf", width=10, height=5)
ggplot(df_long, aes(x = Updated_Pseudotime, y = Expression, color = Condition)) + 
  geom_smooth(method = "loess", se = TRUE, span = 0.5) + 
  scale_color_manual(values = color_map) + 
  facet_wrap(~ClusterLabel, scales = "free_y", nrow = 2) +  # Use the updated ClusterLabel column
  theme_bw() +
  labs(title = "Smoothed Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")
dev.off()



### Export gene list from each cluster
## Create a data frame with gene names and their respective cluster assignments
output_df <- data.frame(
  gene = rownames(combinedData),
  cluster = clusters
)

# Write the data frame to a .txt file
write.table(output_df, 
            file = "output/condiments/gene_clusters_traj2_RNA_common_l2fc0_cl10.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)

# Check some genes individually
## FOR LINEAGE 2
counts <- multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(RNA_WT_Bap1KO, na = FALSE) [,2] # HERE INDICATE TRAJ
cellweights <- slingCurveWeights(RNA_WT_Bap1KO) [,2] # HERE INDICATE TRAJ
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]

pdf("output/condiments/plotSmoothers_traj2_RNA_common-Eomes.pdf", width=4, height=2)
plotSmoothers(traj2_RNA_common, sub_counts, gene = "Eomes", curvesCol = c("blue","red") ) +
scale_color_manual(values =c("blue","red"))
dev.off()


# Heatmap representation
# Define colors for each cluster
# 20
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen", "grey", "darkred", "darkblue", "gold", "darkgray", "lightblue", "lightgreen", "lightcoral", "lightpink", "lightcyan"))(20),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 10
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen" ))(10),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 8
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown" ))(8),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 7
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink" ))(7),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 6
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange"))(6),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 5
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple"))(5),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 4
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow" ))(4),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Generate the heatmap

col_order <- order(grepl("WT", colnames(combinedData)), decreasing = TRUE)
combinedData <- combinedData[, col_order]
pdf("output/condiments/heatmap_traj2_RNA_common_l2fc0_cl10.pdf", width=5, height=5)
pheatmap(combinedData,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  legend = TRUE,
  cutree_rows = 10,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
)
dev.off()






```


Let's do GO/pathway analysis on:
- clustered pseudotime -DEG (traj9 and 2)
- DAR NSC overlapping with pseudotime -DEG (traj9 and 2)


```bash
conda activate deseq2
```

```R

# library
library("tidyverse")
library("enrichR")
library("ggrepel")



# Code for quick overview; all cluster together; ie. multiple gene list at one ####################


# Define the GO database to use
dbs <- c("GO_Biological_Process_2023")

# Read the gene cluster file
clusters <- read_tsv("output/condiments/gene_clusters_traj9_RNA_common_l2fc0_cl10.txt")
clusters <- read_tsv("output/condiments/gene_clusters_traj2_RNA_common_l2fc0_cl10.txt")
clusters <- read_tsv("output/condiments/vennDiagram_DARNSC_DEGtraj9.txt") %>% # added a column cluster to make code work
  add_column(cluster = "1")
clusters <- read_tsv("output/condiments/vennDiagram_DARNSC_DEGtraj2.txt") %>% # added a column cluster to make code work
  add_column(cluster = "1")

# Iterate over each cluster, run enrichment analysis, and store results
go_results_list <- list()

# Iterate over each cluster
for (cluster_id in unique(clusters$cluster)) {
  # Extract genes for the cluster
  cluster_genes <- clusters %>%
    filter(cluster == cluster_id) %>%
    pull(gene) %>%
    unique()
  
  cat("\nProcessing Cluster:", cluster_id, "with", length(cluster_genes), "genes\n")
  print(head(cluster_genes, 5))  # Print first 5 genes for this cluster

  # Run enrichment
  enriched <- enrichr(cluster_genes, dbs)$GO_Biological_Process_2023
  
  # Introduce a 2-second pause between API calls
  Sys.sleep(2)

  # Check the enrichment result
  if (!is.null(enriched) && nrow(enriched) > 0) {
    enriched <- enriched %>%
      filter(Adjusted.P.value <= 0.05) %>%
      mutate(cluster = cluster_id, logAdjP = -log10(Adjusted.P.value))
    go_results_list[[as.character(cluster_id)]] <- enriched
  } else {
    cat("No enrichment results for cluster:", cluster_id, "\n")
  }
}

go_results <- bind_rows(go_results_list)

# Simplify GO terms and prepare for visualization
go_results <- go_results %>%
  separate(Overlap, into = c("count", "total"), sep = "/", convert = TRUE) %>%
  mutate(
    Term = gsub("\\(GO:[0-9]+\\)", "", Term),  # Remove GO IDs from terms
    Term = str_trunc(Term, 50),  # Truncate long terms for readability
    logAdjP = -log10(Adjusted.P.value),  # Convert adjusted p-value to -log10 scale
    geneRatio = (count / total) * 100
  ) 


# Select the top terms per cluster for visualization
top_go_results <- go_results %>%
  group_by(cluster) %>%
  slice_max(logAdjP, n = 10)  # Top 10 terms per cluster

# Reorder terms for better plotting
top_go_results <- top_go_results %>%
  mutate(Term = factor(Term, levels = unique(Term[order(logAdjP)]))) %>%
  mutate(cluster = factor(cluster))

top_go_results <- top_go_results %>%
  mutate(cluster = factor(cluster, levels = 1:10))  # Ensure correct ordering of clusters

# Define custom colors for clusters
cluster_colors <- setNames(
  colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen"))(10),
  as.character(1:10)  # Map colors to cluster levels
)

# Generate the plot
pdf("output/condiments/enrichR_GO_Biological_Process_2023_traj2_RNA_common_l2fc0_cl10.pdf", width = 15, height = 5)
ggplot(top_go_results, aes(x = logAdjP, y = Term, color = cluster, size = Combined.Score)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ cluster, scales = "free", nrow = 3) +
  scale_color_manual(values = cluster_colors) +  # Use custom color mapping
  labs(
    title = "GO Biological Process Enrichment by Cluster",
    x = "-log10(Adjusted P-value)",
    y = "GO Biological Process",
    size = "Combined Score",
    color = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  )
dev.off()


# enrichR_GO_Biological_Process_2023_traj2_RNA_common_l2fc0_cl10_cluster1 
pdf("output/condiments/enrichR_GO_Biological_Process_2023_vennDiagram_DARNSC_DEGtraj2.pdf", width = 5, height = 3)
top_go_results %>%
  dplyr::filter(cluster == "1") %>%
ggplot(., aes(x = logAdjP, y = Term, color = cluster, size = geneRatio)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ cluster, scales = "free", nrow = 3) +
  scale_color_manual(values = cluster_colors) +  # Use custom color mapping
  labs(
    title = "GO Biological Process Enrichment by Cluster",
    x = "-log10(Adjusted P-value)",
    y = "GO Biological Process",
    size = "gene ratio",
    color = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  )
dev.off()






# Define the GO database to use
dbs <- c("KEGG_2019_Mouse")

# Read the gene cluster file
clusters <- read_tsv("output/condiments/gene_clusters_traj9_RNA_common_l2fc0_cl10.txt")
clusters <- read_tsv("output/condiments/gene_clusters_traj2_RNA_common_l2fc0_cl10.txt")

# Iterate over each cluster, run enrichment analysis, and store results
go_results_list <- list()

# Iterate over each cluster
for (cluster_id in unique(clusters$cluster)) {
  # Extract genes for the cluster
  cluster_genes <- clusters %>%
    filter(cluster == cluster_id) %>%
    pull(gene) %>%
    unique()
  
  cat("\nProcessing Cluster:", cluster_id, "with", length(cluster_genes), "genes\n")
  print(head(cluster_genes, 5))  # Print first 5 genes for this cluster

  # Run enrichment
  enriched <- enrichr(cluster_genes, dbs)$KEGG_2019_Mouse
  
  # Introduce a 2-second pause between API calls
  Sys.sleep(2)

  # Check the enrichment result
  if (!is.null(enriched) && nrow(enriched) > 0) {
    enriched <- enriched %>%
      filter(Adjusted.P.value <= 0.05) %>%
      mutate(cluster = cluster_id, logAdjP = -log10(Adjusted.P.value))
    go_results_list[[as.character(cluster_id)]] <- enriched
  } else {
    cat("No enrichment results for cluster:", cluster_id, "\n")
  }
}

go_results <- bind_rows(go_results_list)

# Simplify GO terms and prepare for visualization
go_results <- go_results %>%
  mutate(
    Term = gsub("\\(GO:[0-9]+\\)", "", Term),  # Remove GO IDs from terms
    Term = str_trunc(Term, 50),  # Truncate long terms for readability
    logAdjP = -log10(Adjusted.P.value)  # Convert adjusted p-value to -log10 scale
  )

# Select the top terms per cluster for visualization
top_go_results <- go_results %>%
  group_by(cluster) %>%
  slice_max(logAdjP, n = 10)  # Top 10 terms per cluster

# Reorder terms for better plotting
top_go_results <- top_go_results %>%
  mutate(Term = factor(Term, levels = unique(Term[order(logAdjP)]))) %>%
  mutate(cluster = factor(cluster))

top_go_results <- top_go_results %>%
  mutate(cluster = factor(cluster, levels = 1:10))  # Ensure correct ordering of clusters

# Define custom colors for clusters
cluster_colors <- setNames(
  colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen"))(10),
  as.character(1:10)  # Map colors to cluster levels
)

# Generate the plot
pdf("output/condiments/enrichR_KEGG_2019_Mouse_traj2_RNA_common_l2fc0_cl10.pdf", width = 15, height = 5)
ggplot(top_go_results, aes(x = logAdjP, y = Term, color = cluster, size = Combined.Score)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ cluster, scales = "free", nrow = 3) +
  scale_color_manual(values = cluster_colors) +  # Use custom color mapping
  labs(
    title = "GO Biological Process Enrichment by Cluster",
    x = "-log10(Adjusted P-value)",
    y = "GO Biological Process",
    size = "Combined Score",
    color = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  )
dev.off()






```









#### Run condiments RNA assay - V2 separated condition

Identify pseudotime traj for WT and Bap1KO independently; Focus on previous common trajectory 2 and 9 (DG_GC, CA123)



```bash
conda activate condiments_Signac
```



```R
# package installation 
## install.packages("remotes")
## remotes::install_github("cran/spatstat.core")
## remotes::install_version("Seurat", "4.0.3")
## install.packages("magrittr")
## install.packages("magrittr")
## install.packages("dplyr")
## BiocManager::install("DelayedMatrixStats")
## BiocManager::install("tradeSeq")


# packages
library("condiments")
library("Seurat")
library("magrittr") # to use pipe
library("dplyr") # to use bind_cols and sample_frac
library("SingleCellExperiment") # for reducedDims
library("ggplot2")
library("slingshot")
library("DelayedMatrixStats")
library("tidyr")
library("tradeSeq") 
library("cowplot")
library("scales")
library("pheatmap")

# Data import and separate conditons

multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000.sct_numeric_label.rds")


DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA" # According to condiments workflow


multiome_WT_QCV2vC1.sct <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, orig.ident == "multiome_WT")
multiome_Bap1KO_QCV2vC1.sct <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, orig.ident == "multiome_Bap1KO")


################################################################################################
### WT ################################################
################################################################################################

DefaultAssay(multiome_WT_QCV2vC1.sct) <- "RNA" # According to condiments workflow

# convert to SingleCellExperiment
RNA_WT <- as.SingleCellExperiment(multiome_WT_QCV2vC1.sct, assay = "RNA")

#  Trajectory Inference and Differential Topology
set.seed(42)


### Testing area ############


##########################################

RNA_WT <- slingshot(RNA_WT, reducedDim = 'UMAP',
                 clusterLabels = colData(RNA_WT)$cluster.annot,
                 start.clus = 'NSC_quiescent', end.clus = c("PyNs_RSC_UL","DG_GC") ,approx_points = 100, extend = 'n')




df_2 <- bind_cols(
  as.data.frame(reducedDim(RNA_WT, "UMAP")),
  slingPseudotime(RNA_WT) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(RNA_WT) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst),
        # pst = max(pst) - pst)
)
curves <- slingCurves(RNA_WT, as.df = TRUE)




pdf("output/condiments/UMAP_trajectory_common_RNA_WT_StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn.pdf", width=5, height=5)
ggplot(df_2, aes(x = umap_1, y = umap_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1.5) +
  theme_classic()
dev.off()




## PLOT with trajectories - Individually
df_2 <- bind_cols(
  as.data.frame(reducedDim(RNA_WT, "UMAP")),
  slingPseudotime(RNA_WT) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(RNA_WT) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         Lineage3_pst = if_else(is.na(Lineage3_pst), 0, Lineage3_pst),
         Lineage4_pst = if_else(is.na(Lineage4_pst), 0, Lineage4_pst),
         Lineage5_pst = if_else(is.na(Lineage5_pst), 0, Lineage5_pst),
         Lineage6_pst = if_else(is.na(Lineage6_pst), 0, Lineage6_pst),
         Lineage7_pst = if_else(is.na(Lineage7_pst), 0, Lineage7_pst),
         Lineage8_pst = if_else(is.na(Lineage8_pst), 0, Lineage8_pst),
         Lineage9_pst = if_else(is.na(Lineage9_pst), 0, Lineage9_pst))
curves <- slingCurves(RNA_WT, as.df = TRUE)
### Function to create the plot for each lineage
create_plot <- function(lineage_number) {
  df_2 <- df_2 %>%
    mutate(pst = case_when(
      !!sym(paste0("Lineage", lineage_number, "_pst")) > 0 ~ !!sym(paste0("Lineage", lineage_number, "_pst")),
      TRUE ~ 0
    ),
    group = if_else(pst > 0, paste0("lineage", lineage_number), "other"))
  curves_filtered <- curves %>% filter(Lineage == lineage_number)
  curves_endpoints <- curves_filtered %>%
    group_by(Lineage) %>%
    arrange(Order) %>%
    top_n(1, Order) # Get the top/last ordered point for each group
  df_2_lineage <- df_2 %>% filter(group == paste0("lineage", lineage_number))
  df_2_other <- df_2 %>% filter(group == "other")
  p <- ggplot() +
    geom_point(data = df_2_other, aes(x = umap_1, y = umap_2), size = .7, color = "grey85") +
    geom_point(data = df_2_lineage, aes(x = umap_1, y = umap_2, col = pst), size = .7) +
    scale_color_viridis_c() +
    labs(col = "Pseudotime", title = paste("Lineage", lineage_number)) +
    geom_path(data = curves_filtered %>% arrange(Order),
              aes(x = umap_1, y = umap_2, group = Lineage), col = "black", size = 1) +
    geom_text(data = curves_endpoints, aes(x = umap_1, y = umap_2, label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
    theme_classic()
  return(p)
}
### Generate the plots for each lineage
plots <- list()
for (i in 1:9) {
  plots[[i]] <- create_plot(i)
}
pdf("output/condiments/UMAP_trajectory_common_label_RNA_WT_StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn_Lineage123456789.pdf", width=25, height=7)
gridExtra::grid.arrange(grobs = plots, ncol = 5)
dev.off()




#### ->  save.image(file="output/condiments/condiments_RNA_WT_StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn.RData")
### load("output/condiments/condiments_RNA_WT_StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn.RData")
set.seed(42)
##
# RUN fitGAM
##


################### Time Course effect (UNTREATED CONDITION) - TRAJECTORY 9 - ##############################################
# NSC --> DG_GC
set.seed(42)
traj9_RNA_WT <- readRDS("output/condiments/traj9_RNA_WT.rds")
## Genes that change with pseudotime

pseudotime_association <- associationTest(traj9_RNA_WT) # statistical test to check whether gene expression is constant across pseudotime within a lineage
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$fdr), ]
pseudotime_association$gene <- rownames(pseudotime_association)

pseudotime_association = as_tibble(pseudotime_association) 
# save output: write.table(pseudotime_association, file = c("output/condiments/pseudotime_association_traj9_RNA_WT.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## Genes that change between two pseudotime points (start vs end)
pseudotime_start_end_association <- startVsEndTest(traj9_RNA_WT, pseudotimeValues = NULL)
pseudotime_start_end_association$gene <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association$fdr <- p.adjust(pseudotime_start_end_association$pvalue, method = "fdr")
pseudotime_start_end_association <- pseudotime_start_end_association[order(pseudotime_start_end_association$fdr), ]
##--> log2FC = end - start: negative log2fc means start point higher average expr than end point
# save output: write.table(pseudotime_start_end_association, file = c("output/condiments/pseudotime_start_end_association_traj9_RNA_WT.txt"),sep="\t", quote=FALSE, row.names=FALSE)


sce_cells <- colnames(traj9_RNA_WT) # collect cells of traj1
subset_traj9_RNA_WT_multiome_WT_Bap1KO_QCV2vC1.sct <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, cells = sce_cells) # Create a seurat object with only cells from traj1
subset_traj9_RNA_WT_multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct, cells = sce_cells) # Create a seurat object with only cells from traj1



### plot top 25 genes
#### Select the top 25 genes with positive logFC and the top 25 with negative logFC
top25_posFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 > 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

top25_negFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 < 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

## Function to plot and save in PDF
plot_and_save <- function(genes_df, file_name) {
  pdf(file_name, width=5, height=4)
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$gene[i]
    logFC <- genes_df$logFClineage1[i]
    plot_title <- paste0(gene, " (logFC: ", round(logFC, 2), ")")
    p <- plotSmoothers(traj9_RNA_WT, subset_traj9_RNA_WT_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts, gene = gene)
    p <- p + ggtitle(plot_title)
    print(p)
  }
  dev.off()
}

# Generate PDFs
plot_and_save(top25_posFC_genes, "output/condiments/plotSmoothers-top25_posFC_traj9_RNA_WT.pdf")
plot_and_save(top25_negFC_genes, "output/condiments/plotSmoothers-top25_negFC_traj9_RNA_WT.pdf")


### plot unique genes (Pax6, Tbr2/Eomes, Prox1)
pdf("output/condiments/plotSmoothers-Prox1-traj9_RNA_WT.pdf",  width=5, height=4)
plotSmoothers(traj9_RNA_WT, subset_traj9_RNA_WT_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts, gene = "Prox1" )
dev.off()


### plot average of multiple genes - multiple gene lists
#### Define the gene list
gene_list_1 <- c("Pax6")
gene_list_2 <- c("Eomes")
gene_list_3 <- c("Prox1")
# Extract the counts for the genes of interest
counts_matrix <- subset_traj9_RNA_WT_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts
counts_matrix <- subset_traj9_RNA_WT_multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct[["GeneActivity"]]@counts

# Subset the data for each gene list and calculate the average expression
subset_data_1 <- counts_matrix[gene_list_1, , drop = FALSE]
subset_data_2 <- counts_matrix[gene_list_2, , drop = FALSE]
subset_data_3 <- counts_matrix[gene_list_3, , drop = FALSE]
# Extract pseudotime
pseudotime_vector <- colData(traj9_RNA_WT)$crv
# Create data frames for plotting
plot_data_1 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_1)), GeneList = "GeneList 1")
plot_data_2 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_2)), GeneList = "GeneList 2")
plot_data_3 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_3)), GeneList = "GeneList 3")
# Combine the data frames
plot_data <- rbind(plot_data_1, plot_data_2, plot_data_3)
plot_data <- plot_data %>%
  group_by(GeneList) %>%
  mutate(SmoothedExpression = predict(loess(Expression ~ Pseudotime.pseudotime, span = 0.5)))
# Combine the smoothed data frames for plotting
smoothed_plot_data <- rbind(plot_data_1, plot_data_2, plot_data_3)
# Plot the smoothed trajectories
pdf("output/condiments/plotSmoothers-Pax6EomesProx1-traj9_RNA_WT.pdf", width = 5, height = 4)
pdf("output/condiments/plotSmoothers-Pax6EomesProx1-traj9_ATAC_WT.pdf", width = 5, height = 4)

ggplot(plot_data, aes(x = Pseudotime.pseudotime, y = Expression, color = GeneList)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = SmoothedExpression), size = 1) +
  scale_color_manual(values = c("GeneList 1" = "green", "GeneList 2" = "blue", "GeneList 3" = "purple")) + # Custom colors
  labs(
    x = "Pseudotime",
    y = "Log(Expression + 1)",
    color = "Gene List"
  ) +
  theme_bw()
dev.off()






## Identify Activation point = peak (maximum expression) of each gene along the pseudotime trajectory
### Identify peak of expression (max expr) of these Time-course DEG
traj9_RNA_WT
#### Extract pseudotime values
pseudotime <- colData(traj9_RNA_WT)$crv$pseudotime
#### Extract the expression matrix
expr_matrix <- assays(traj9_RNA_WT)$counts
#### Ensure the pseudotime values are named with the same cell names as the expression matrix columns
names(pseudotime) <- colnames(expr_matrix)
#### Function to find the peak pseudotime for each gene (raw and smoothed)
find_max_pseudotime <- function(gene_expr, pseudotime) {
  # Raw peak pseudotime
  raw_peak_pseudotime <- pseudotime[which.max(gene_expr)]
  # Smooth gene expression using loess
  smooth_model <- loess(gene_expr ~ pseudotime)
  smooth_expr <- predict(smooth_model)
  # Smooth peak pseudotime
  smooth_peak_pseudotime <- pseudotime[which.max(smooth_expr)]
  return(list(raw_peak_pseudotime = raw_peak_pseudotime, 
              smooth_peak_pseudotime = smooth_peak_pseudotime))
}
#### Apply the function to all genes
peak_values <- apply(expr_matrix, 1, function(x) find_max_pseudotime(as.numeric(x), pseudotime))
#### Convert the results to a data frame
peak_df <- data.frame(
  gene = rownames(expr_matrix),
  raw_peak_pseudotime = sapply(peak_values, `[[`, "raw_peak_pseudotime"),
  smooth_peak_pseudotime = sapply(peak_values, `[[`, "smooth_peak_pseudotime")
) %>% as_tibble()


# save output: write.table(peak_df, file = c("output/condiments/traj9_RNA_WT_ActivationPoint.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## heatmap activate/induced genes along pseudotime
### DEG Start End
pseudotime_start_end_association # filter log2fc >0 >1
pseudotime_start_end_association = read_tsv("output/condiments/pseudotime_start_end_association_traj9_RNA_WT.txt")
pseudotime_start_end_association_logFC0 = pseudotime_start_end_association %>% 
  filter(fdr <0.05) %>%
  dplyr::select(gene) %>%
  unique()

pdf("output/condiments/heatmap_pseudotime_start_end_association_traj9_RNA_WT_fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj9_RNA_WT, gene = pseudotime_start_end_association_logFC0$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()




### DEG time course
pseudotime_association = read_tsv("output/condiments/pseudotime_association_traj9_RNA_WT.txt")
pseudotime_association_deg = pseudotime_association %>%
  filter(fdr <0.05)%>%
  dplyr::select(gene) %>%
  unique()


pdf("output/condiments/heatmap_pseudotime_association_traj9_RNA_WT_fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj9_RNA_WT, gene = pseudotime_association_deg$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

#### Gene Activity (ATAC) pseudotime
##### Load Seurat object with GeneActivity information
multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds")
#### import LinkGenes
Gene_LinkPeak_signif = read_tsv("output/Signac/LinkPeaks_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1.txt") %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%
  dplyr::select(gene) %>%
  unique()
#### Extract pseudotime and cell information
pseudotime_data <- colData(traj9_RNA_WT)$crv
pseudotime_values <- pseudotime_data$pseudotime
names(pseudotime_values) <- rownames(pseudotime_data)
cell_ids <- colnames(traj9_RNA_WT)
#### Ensure matching cells between pseudotime and Seurat object
common_cells <- intersect(cell_ids, colnames(multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct))
#### Subset pseudotime and GeneActivity for the common cells
pseudotime_values <- pseudotime_values[common_cells]
#### Subset to only available genes
valid_genes <- intersect(pseudotime_association_deg$gene, available_genes)

# valid_genes <- intersect(valid_genes, Gene_LinkPeak_signif$gene) # HERE FILTER TO KEEP LINK PEAK ONLY

gene_activity_data <- GetAssayData(multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct, 
                                   assay = "GeneActivity", 
                                   slot = "scale.data")[valid_genes, common_cells]
cell_order <- order(pseudotime_values)
gene_activity_ordered <- gene_activity_data[, cell_order]
#### Order genes by their RNA expression peak
#### Predict RNA expression smoothness for the genes
yhatSmooth <- predictSmooth(traj9_RNA_WT, gene = valid_genes, nPoints = 25, tidy = FALSE)
#### Order genes by maximum expression in pseudotime
gene_order <- order(apply(yhatSmooth, 1, which.max))  # Order by RNA max along pseudotime
gene_activity_final <- gene_activity_ordered[gene_order, ]  # Reorder ATAC signals to match RNA ordering
#### Scale and plot the gene activity heatmap
gene_activity_scaled <- t(scale(t(gene_activity_final)))

n_bins <- 25
# Generate 25 evenly spaced bins across pseudotime
bin_edges <- seq(1, ncol(gene_activity_scaled), length.out = n_bins + 1)
bin_centers <- round((bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2)
# Aggregate gene activity data into bins by averaging
gene_activity_binned <- sapply(seq_along(bin_centers), function(i) {
  cols_to_average <- seq(floor(bin_edges[i]), ceiling(bin_edges[i + 1]) - 1)
  if (length(cols_to_average) > 1) {
    rowMeans(gene_activity_scaled[, cols_to_average, drop = FALSE])
  } else {
    gene_activity_scaled[, cols_to_average, drop = FALSE]
  }
})

gene_activity_binned_scaled <- t(scale(t(gene_activity_binned)))

#_LinkPeaks
pdf("output/condiments/heatmap_pseudotime_association_deg_traj9_ATAC_WT_fdr05.pdf", width = 8, height = 10)
pheatmap(gene_activity_binned_scaled,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = FALSE, 
         show_colnames = FALSE,
         breaks = seq(-1, 2, length.out = 101),
         labels_col = paste0("Bin ", seq_len(n_bins))) # Add bin labels for clarity
dev.off()


pdf("output/condiments/heatmap_pseudotime_association_traj9_RNA_WT_fdr05_ATACgenes.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj9_RNA_WT, gene = valid_genes, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()





## Identify which pseudotime value corespond to which cluster ################
pseudotime <- colData(traj9_RNA_WT)$crv$pseudotime
sce_cells <- colnames(traj9_RNA_WT)
subset_seurat <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, cells = sce_cells) # Subset cell from traj2
clusters <- subset_seurat$cluster.annot # Extract cluster information
### Combine pseudotime and cluster information into a data frame
pseudotime_cluster_df <- data.frame(
  cell = colnames(traj9_RNA_WT),
  pseudotime = pseudotime,
  cluster = clusters
)  %>%
  arrange(pseudotime)

switch_points <- which(diff(as.numeric(factor(pseudotime_cluster_df$cluster))) != 0) # Find the indices where the cluster changes
switch_pseudotimes <- pseudotime_cluster_df$pseudotime[switch_points] # Extract the pseudotime values at these switch points
switch_clusters_from <- pseudotime_cluster_df$cluster[switch_points]
switch_clusters_to <- pseudotime_cluster_df$cluster[switch_points + 1]
switch_df <- data.frame(
  switch_pseudotime = switch_pseudotimes,
  cluster_from = switch_clusters_from,
  cluster_to = switch_clusters_to 
) %>%
  group_by(cluster_from, cluster_to) %>%
  summarize(median_switch_pseudotime = median(switch_pseudotime), .groups = 'drop')
write.table(switch_df, file = c("output/condiments/switch_df_traj9_RNA_WT.txt"),sep="\t", quote=FALSE, row.names=FALSE)
##################################








################### Time Course effect (UNTREATED CONDITION) - TRAJECTORY 2 - ##############################################
# NSC --> CA123
set.seed(42)
traj2_RNA_WT <- readRDS("output/condiments/traj2_RNA_WT.rds")
## Genes that change with pseudotime

pseudotime_association <- associationTest(traj2_RNA_WT) # statistical test to check whether gene expression is constant across pseudotime within a lineage
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$fdr), ]
pseudotime_association$gene <- rownames(pseudotime_association)
pseudotime_association = as_tibble(pseudotime_association) 
# save output: write.table(pseudotime_association, file = c("output/condiments/pseudotime_association_traj2_RNA_WT.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## Genes that change between two pseudotime points (start vs end)
pseudotime_start_end_association <- startVsEndTest(traj2_RNA_WT, pseudotimeValues = NULL)
pseudotime_start_end_association$gene <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association$fdr <- p.adjust(pseudotime_start_end_association$pvalue, method = "fdr")
pseudotime_start_end_association <- pseudotime_start_end_association[order(pseudotime_start_end_association$fdr), ]
##--> log2FC = end - start: negative log2fc means start point higher average expr than end point
# save output: write.table(pseudotime_start_end_association, file = c("output/condiments/pseudotime_start_end_association_traj2_RNA_WT.txt"),sep="\t", quote=FALSE, row.names=FALSE)


sce_cells <- colnames(traj2_RNA_WT) # collect cells of traj1
subset_traj2_RNA_WT_multiome_WT_Bap1KO_QCV2vC1.sct <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, cells = sce_cells) # Create a seurat object with only cells from traj1
subset_traj2_RNA_WT_multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct, cells = sce_cells) # Create a seurat object with only cells from traj1



### plot top 25 genes
#### Select the top 25 genes with positive logFC and the top 25 with negative logFC
top25_posFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 > 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

top25_negFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 < 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

## Function to plot and save in PDF
plot_and_save <- function(genes_df, file_name) {
  pdf(file_name, width=5, height=4)
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$gene[i]
    logFC <- genes_df$logFClineage1[i]
    plot_title <- paste0(gene, " (logFC: ", round(logFC, 2), ")")
    p <- plotSmoothers(traj2_RNA_WT, subset_traj2_RNA_WT_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts, gene = gene)
    p <- p + ggtitle(plot_title)
    print(p)
  }
  dev.off()
}

# Generate PDFs
plot_and_save(top25_posFC_genes, "output/condiments/plotSmoothers-top25_posFC_traj2_RNA_WT.pdf")
plot_and_save(top25_negFC_genes, "output/condiments/plotSmoothers-top25_negFC_traj2_RNA_WT.pdf")


### plot unique genes (Pax6, Tbr2/Eomes, Prox1)
pdf("output/condiments/plotSmoothers-Prox1-traj2_RNA_WT.pdf",  width=5, height=4)
plotSmoothers(traj2_RNA_WT, subset_traj2_RNA_WT_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts, gene = "Prox1" )
dev.off()


### plot average of multiple genes - multiple gene lists
#### Define the gene list
gene_list_1 <- c("Pax6")
gene_list_2 <- c("Eomes")
gene_list_3 <- c("Prox1")
# Extract the counts for the genes of interest
counts_matrix <- subset_traj2_RNA_WT_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts
counts_matrix <- subset_traj2_RNA_WT_multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct[["GeneActivity"]]@counts
# Subset the data for each gene list and calculate the average expression
subset_data_1 <- counts_matrix[gene_list_1, , drop = FALSE]
subset_data_2 <- counts_matrix[gene_list_2, , drop = FALSE]
subset_data_3 <- counts_matrix[gene_list_3, , drop = FALSE]
# Extract pseudotime
pseudotime_vector <- colData(traj2_RNA_WT)$crv
# Create data frames for plotting
plot_data_1 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_1)), GeneList = "GeneList 1")
plot_data_2 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_2)), GeneList = "GeneList 2")
plot_data_3 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_3)), GeneList = "GeneList 3")

# Combine the data frames
plot_data <- rbind(plot_data_1, plot_data_2, plot_data_3)

plot_data <- plot_data %>%
  group_by(GeneList) %>%
  mutate(SmoothedExpression = predict(loess(Expression ~ Pseudotime.pseudotime, span = 0.5)))


# Combine the smoothed data frames for plotting
smoothed_plot_data <- rbind(plot_data_1, plot_data_2, plot_data_3)
# Plot the smoothed trajectories
pdf("output/condiments/plotSmoothers-Pax6EomesProx1-traj2_RNA_WT.pdf", width = 5, height = 4)
ggplot(plot_data, aes(x = Pseudotime.pseudotime, y = Expression, color = GeneList)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = SmoothedExpression), size = 1) +
  scale_color_manual(values = c("GeneList 1" = "green", "GeneList 2" = "blue", "GeneList 3" = "purple")) + # Custom colors
  labs(
    x = "Pseudotime",
    y = "Log(Expression + 1)",
    color = "Gene List"
  ) +
  theme_bw()
dev.off()


## Identify Activation point = peak (maximum expression) of each gene along the pseudotime trajectory
### Identify peak of expression (max expr) of these Time-course DEG
traj2_RNA_WT
#### Extract pseudotime values
pseudotime <- colData(traj2_RNA_WT)$crv$pseudotime
#### Extract the expression matrix
expr_matrix <- assays(traj2_RNA_WT)$counts
#### Ensure the pseudotime values are named with the same cell names as the expression matrix columns
names(pseudotime) <- colnames(expr_matrix)
#### Function to find the peak pseudotime for each gene (raw and smoothed)
find_max_pseudotime <- function(gene_expr, pseudotime) {
  # Raw peak pseudotime
  raw_peak_pseudotime <- pseudotime[which.max(gene_expr)]
  # Smooth gene expression using loess
  smooth_model <- loess(gene_expr ~ pseudotime)
  smooth_expr <- predict(smooth_model)
  # Smooth peak pseudotime
  smooth_peak_pseudotime <- pseudotime[which.max(smooth_expr)]
  return(list(raw_peak_pseudotime = raw_peak_pseudotime, 
              smooth_peak_pseudotime = smooth_peak_pseudotime))
}
#### Apply the function to all genes
peak_values <- apply(expr_matrix, 1, function(x) find_max_pseudotime(as.numeric(x), pseudotime))
#### Convert the results to a data frame
peak_df <- data.frame(
  gene = rownames(expr_matrix),
  raw_peak_pseudotime = sapply(peak_values, `[[`, "raw_peak_pseudotime"),
  smooth_peak_pseudotime = sapply(peak_values, `[[`, "smooth_peak_pseudotime")
) %>% as_tibble()
# save output: write.table(peak_df, file = c("output/condiments/traj2_RNA_WT_ActivationPoint.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## heatmap activate/induced genes along pseudotime
### DEG Start End
pseudotime_start_end_association # filter log2fc >0 >1
pseudotime_start_end_association = read_tsv("output/condiments/pseudotime_start_end_association_traj2_RNA_WT.txt")
pseudotime_start_end_association_logFC0 = pseudotime_start_end_association %>% 
  filter(fdr <0.05) %>%
  dplyr::select(gene) %>%
  unique()

pdf("output/condiments/heatmap_pseudotime_start_end_association_traj2_RNA_WT_fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj2_RNA_WT, gene = pseudotime_start_end_association_logFC0$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG time course
pseudotime_association = read_tsv("output/condiments/pseudotime_association_traj2_RNA_WT.txt")
pseudotime_association_deg = pseudotime_association %>%
  filter(fdr <0.05)%>%
  dplyr::select(gene) %>%
  unique()


pdf("output/condiments/heatmap_pseudotime_association_traj2_RNA_WT_fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj2_RNA_WT, gene = pseudotime_association_deg$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()


#### Gene Activity (ATAC) pseudotime
##### Load Seurat object with GeneActivity information
multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds")
#### import LinkGenes
Gene_LinkPeak_signif = read_tsv("output/Signac/LinkPeaks_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1.txt") %>%
  dplyr::filter(adjusted_pvalue < 0.05, score >0) %>%
  dplyr::select(gene) %>%
  unique()
#### Extract pseudotime and cell information
pseudotime_data <- colData(traj2_RNA_WT)$crv
pseudotime_values <- pseudotime_data$pseudotime
names(pseudotime_values) <- rownames(pseudotime_data)
cell_ids <- colnames(traj2_RNA_WT)
#### Ensure matching cells between pseudotime and Seurat object
common_cells <- intersect(cell_ids, colnames(multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct))
#### Subset pseudotime and GeneActivity for the common cells
pseudotime_values <- pseudotime_values[common_cells]
#### Subset to only available genes
valid_genes <- intersect(pseudotime_association_deg$gene, available_genes)

# valid_genes <- intersect(valid_genes, Gene_LinkPeak_signif$gene) # HERE FILTER TO KEEP LINK PEAK ONLY

gene_activity_data <- GetAssayData(multiome_WT_Bap1KO_QCV2vC1_GeneActivity.sct, 
                                   assay = "GeneActivity", 
                                   slot = "scale.data")[valid_genes, common_cells]
cell_order <- order(pseudotime_values)
gene_activity_ordered <- gene_activity_data[, cell_order]
#### Order genes by their RNA expression peak
#### Predict RNA expression smoothness for the genes
yhatSmooth <- predictSmooth(traj2_RNA_WT, gene = valid_genes, nPoints = 25, tidy = FALSE)
#### Order genes by maximum expression in pseudotime
gene_order <- order(apply(yhatSmooth, 1, which.max))  # Order by RNA max along pseudotime
gene_activity_final <- gene_activity_ordered[gene_order, ]  # Reorder ATAC signals to match RNA ordering
#### Scale and plot the gene activity heatmap
gene_activity_scaled <- t(scale(t(gene_activity_final)))

n_bins <- 25
# Generate 25 evenly spaced bins across pseudotime
bin_edges <- seq(1, ncol(gene_activity_scaled), length.out = n_bins + 1)
bin_centers <- round((bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2)
# Aggregate gene activity data into bins by averaging
gene_activity_binned <- sapply(seq_along(bin_centers), function(i) {
  cols_to_average <- seq(floor(bin_edges[i]), ceiling(bin_edges[i + 1]) - 1)
  if (length(cols_to_average) > 1) {
    rowMeans(gene_activity_scaled[, cols_to_average, drop = FALSE])
  } else {
    gene_activity_scaled[, cols_to_average, drop = FALSE]
  }
})

gene_activity_binned_scaled <- t(scale(t(gene_activity_binned)))

#_LinkPeaks
pdf("output/condiments/heatmap_pseudotime_association_deg_traj2_ATAC_WT_fdr05.pdf", width = 8, height = 10)
pheatmap(gene_activity_binned_scaled,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = FALSE, 
         show_colnames = FALSE,
         breaks = seq(-1, 2, length.out = 101),
         labels_col = paste0("Bin ", seq_len(n_bins))) # Add bin labels for clarity
dev.off()


pdf("output/condiments/heatmap_pseudotime_association_traj2_RNA_WT_fdr05_ATACgenes.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj2_RNA_WT, gene = valid_genes, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()











## Identify which pseudotime value corespond to which cluster ################
pseudotime <- colData(traj2_RNA_WT)$crv$pseudotime
sce_cells <- colnames(traj2_RNA_WT)
subset_seurat <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, cells = sce_cells) # Subset cell from traj2
clusters <- subset_seurat$cluster.annot # Extract cluster information
### Combine pseudotime and cluster information into a data frame
pseudotime_cluster_df <- data.frame(
  cell = colnames(traj2_RNA_WT),
  pseudotime = pseudotime,
  cluster = clusters
)  %>%
  arrange(pseudotime)

switch_points <- which(diff(as.numeric(factor(pseudotime_cluster_df$cluster))) != 0) # Find the indices where the cluster changes
switch_pseudotimes <- pseudotime_cluster_df$pseudotime[switch_points] # Extract the pseudotime values at these switch points
switch_clusters_from <- pseudotime_cluster_df$cluster[switch_points]
switch_clusters_to <- pseudotime_cluster_df$cluster[switch_points + 1]
switch_df <- data.frame(
  switch_pseudotime = switch_pseudotimes,
  cluster_from = switch_clusters_from,
  cluster_to = switch_clusters_to 
) %>%
  group_by(cluster_from, cluster_to) %>%
  summarize(median_switch_pseudotime = median(switch_pseudotime), .groups = 'drop')
write.table(switch_df, file = c("output/condiments/switch_df_traj2_RNA_WT.txt"),sep="\t", quote=FALSE, row.names=FALSE)
##################################













################################################################################################
### Bap1KO ################################################
################################################################################################

DefaultAssay(multiome_Bap1KO_QCV2vC1.sct) <- "RNA" # According to condiments workflow

# convert to SingleCellExperiment
RNA_Bap1KO <- as.SingleCellExperiment(multiome_Bap1KO_QCV2vC1.sct, assay = "RNA")

#  Trajectory Inference and Differential Topology
set.seed(42)


### Param used for WT ############
RNA_Bap1KO <- slingshot(RNA_Bap1KO, reducedDim = 'UMAP',
                 clusterLabels = colData(RNA_Bap1KO)$cluster.annot,
                 start.clus = 'NSC_quiescent', end.clus = c("PyNs_RSC_UL","DG_GC") ,approx_points = 100, extend = 'n')

##########################################

RNA_Bap1KO <- slingshot(RNA_Bap1KO, reducedDim = 'UMAP',
                 clusterLabels = colData(RNA_Bap1KO)$cluster.annot,
                 start.clus = 'NSC_quiescent', end.clus = c("PyNs_RSC_UL","DG_GC") ,approx_points = 100, extend = 'n')




df_2 <- bind_cols(
  as.data.frame(reducedDim(RNA_Bap1KO, "UMAP")),
  slingPseudotime(RNA_Bap1KO) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(RNA_Bap1KO) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst),
        # pst = max(pst) - pst)
)
curves <- slingCurves(RNA_Bap1KO, as.df = TRUE)




pdf("output/condiments/UMAP_trajectory_common_RNA_Bap1KO_StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn.pdf", width=5, height=5)
ggplot(df_2, aes(x = umap_1, y = umap_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1.5) +
  theme_classic()
dev.off()




## PLOT with trajectories - Individually
df_2 <- bind_cols(
  as.data.frame(reducedDim(RNA_Bap1KO, "UMAP")),
  slingPseudotime(RNA_Bap1KO) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(RNA_Bap1KO) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         Lineage3_pst = if_else(is.na(Lineage3_pst), 0, Lineage3_pst),
         Lineage4_pst = if_else(is.na(Lineage4_pst), 0, Lineage4_pst),
         Lineage5_pst = if_else(is.na(Lineage5_pst), 0, Lineage5_pst),
         Lineage6_pst = if_else(is.na(Lineage6_pst), 0, Lineage6_pst),
         Lineage7_pst = if_else(is.na(Lineage7_pst), 0, Lineage7_pst),
         Lineage8_pst = if_else(is.na(Lineage8_pst), 0, Lineage8_pst),
         Lineage9_pst = if_else(is.na(Lineage9_pst), 0, Lineage9_pst),
         Lineage10_pst = if_else(is.na(Lineage10_pst), 0, Lineage10_pst))
curves <- slingCurves(RNA_Bap1KO, as.df = TRUE)
### Function to create the plot for each lineage
create_plot <- function(lineage_number) {
  df_2 <- df_2 %>%
    mutate(pst = case_when(
      !!sym(paste0("Lineage", lineage_number, "_pst")) > 0 ~ !!sym(paste0("Lineage", lineage_number, "_pst")),
      TRUE ~ 0
    ),
    group = if_else(pst > 0, paste0("lineage", lineage_number), "other"))
  curves_filtered <- curves %>% filter(Lineage == lineage_number)
  curves_endpoints <- curves_filtered %>%
    group_by(Lineage) %>%
    arrange(Order) %>%
    top_n(1, Order) # Get the top/last ordered point for each group
  df_2_lineage <- df_2 %>% filter(group == paste0("lineage", lineage_number))
  df_2_other <- df_2 %>% filter(group == "other")
  p <- ggplot() +
    geom_point(data = df_2_other, aes(x = umap_1, y = umap_2), size = .7, color = "grey85") +
    geom_point(data = df_2_lineage, aes(x = umap_1, y = umap_2, col = pst), size = .7) +
    scale_color_viridis_c() +
    labs(col = "Pseudotime", title = paste("Lineage", lineage_number)) +
    geom_path(data = curves_filtered %>% arrange(Order),
              aes(x = umap_1, y = umap_2, group = Lineage), col = "black", size = 1) +
    geom_text(data = curves_endpoints, aes(x = umap_1, y = umap_2, label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
    theme_classic()
  return(p)
}
### Generate the plots for each lineage
plots <- list()
for (i in 1:10) {
  plots[[i]] <- create_plot(i)
}
pdf("output/condiments/UMAP_trajectory_common_label_RNA_Bap1KO_StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn_Lineage12345678910.pdf", width=25, height=7)
gridExtra::grid.arrange(grobs = plots, ncol = 5)
dev.off()




#### ->  save.image(file="output/condiments/condiments_RNA_Bap1KO_StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn.RData")
### load("output/condiments/condiments_RNA_Bap1KO_StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn.RData")
set.seed(42)
##
# RUN fitGAM
##



################### Time Course effect (BAP1KO CONDITION) - TRAJECTORY 8 - ##############################################
# NSC --> DG_GC
set.seed(42)
traj8_RNA_Bap1KO <- readRDS("output/condiments/traj8_RNA_Bap1KO.rds")
## Genes that change with pseudotime

pseudotime_association <- associationTest(traj8_RNA_Bap1KO) # statistical test to check whether gene expression is constant across pseudotime within a lineage
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$fdr), ]
pseudotime_association$gene <- rownames(pseudotime_association)
pseudotime_association = as_tibble(pseudotime_association) 
# save output: write.table(pseudotime_association, file = c("output/condiments/pseudotime_association_traj8_RNA_Bap1KO.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## Genes that change between two pseudotime points (start vs end)
pseudotime_start_end_association <- startVsEndTest(traj8_RNA_Bap1KO, pseudotimeValues = NULL)
pseudotime_start_end_association$gene <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association$fdr <- p.adjust(pseudotime_start_end_association$pvalue, method = "fdr")
pseudotime_start_end_association <- pseudotime_start_end_association[order(pseudotime_start_end_association$fdr), ]
##--> log2FC = end - start: negative log2fc means start point higher average expr than end point
# save output: write.table(pseudotime_start_end_association, file = c("output/condiments/pseudotime_start_end_association_traj8_RNA_Bap1KO.txt"),sep="\t", quote=FALSE, row.names=FALSE)


sce_cells <- colnames(traj8_RNA_Bap1KO) # collect cells of traj1
subset_traj8_RNA_Bap1KO_multiome_WT_Bap1KO_QCV2vC1.sct <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, cells = sce_cells) # Create a seurat object with only cells from traj1



### plot top 25 genes
#### Select the top 25 genes with positive logFC and the top 25 with negative logFC
top25_posFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 > 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

top25_negFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 < 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

## Function to plot and save in PDF
plot_and_save <- function(genes_df, file_name) {
  pdf(file_name, width=5, height=4)
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$gene[i]
    logFC <- genes_df$logFClineage1[i]
    plot_title <- paste0(gene, " (logFC: ", round(logFC, 2), ")")
    p <- plotSmoothers(traj8_RNA_Bap1KO, subset_traj8_RNA_Bap1KO_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts, gene = gene)
    p <- p + ggtitle(plot_title)
    print(p)
  }
  dev.off()
}

# Generate PDFs
plot_and_save(top25_posFC_genes, "output/condiments/plotSmoothers-top25_posFC_traj8_RNA_Bap1KO.pdf")
plot_and_save(top25_negFC_genes, "output/condiments/plotSmoothers-top25_negFC_traj8_RNA_Bap1KO.pdf")


### plot unique genes (Pax6, Tbr2/Eomes, Prox1)
pdf("output/condiments/plotSmoothers-Eomes-traj8_RNA_Bap1KO.pdf",  width=5, height=4)
plotSmoothers(traj8_RNA_Bap1KO, subset_traj8_RNA_Bap1KO_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts, gene = "Eomes" )
dev.off()


### plot average of multiple genes - multiple gene lists
#### Define the gene list
gene_list_1 <- c("Pax6")
gene_list_2 <- c("Eomes")
gene_list_3 <- c("Prox1")
# Extract the counts for the genes of interest
counts_matrix <- subset_traj8_RNA_Bap1KO_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts
# Subset the data for each gene list and calculate the average expression
subset_data_1 <- counts_matrix[gene_list_1, , drop = FALSE]
subset_data_2 <- counts_matrix[gene_list_2, , drop = FALSE]
subset_data_3 <- counts_matrix[gene_list_3, , drop = FALSE]
# Extract pseudotime
pseudotime_vector <- colData(traj8_RNA_Bap1KO)$crv
# Create data frames for plotting
plot_data_1 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_1)), GeneList = "GeneList 1")
plot_data_2 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_2)), GeneList = "GeneList 2")
plot_data_3 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_3)), GeneList = "GeneList 3")

# Combine the data frames
plot_data <- rbind(plot_data_1, plot_data_2, plot_data_3)

plot_data <- plot_data %>%
  group_by(GeneList) %>%
  mutate(SmoothedExpression = predict(loess(Expression ~ Pseudotime.pseudotime, span = 0.5)))


# Combine the smoothed data frames for plotting
smoothed_plot_data <- rbind(plot_data_1, plot_data_2, plot_data_3)
# Plot the smoothed trajectories
pdf("output/condiments/plotSmoothers-Pax6EomesProx1-traj8_RNA_Bap1KO.pdf", width = 5, height = 4)
ggplot(plot_data, aes(x = Pseudotime.pseudotime, y = Expression, color = GeneList)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = SmoothedExpression), size = 1) +
  scale_color_manual(values = c("GeneList 1" = "green", "GeneList 2" = "blue", "GeneList 3" = "purple")) + # Custom colors
  labs(
    x = "Pseudotime",
    y = "Log(Expression + 1)",
    color = "Gene List"
  ) +
  theme_bw()
dev.off()






## Identify Activation point = peak (maximum expression) of each gene along the pseudotime trajectory
### Identify peak of expression (max expr) of these Time-course DEG
traj8_RNA_Bap1KO
#### Extract pseudotime values
pseudotime <- colData(traj8_RNA_Bap1KO)$crv$pseudotime
#### Extract the expression matrix
expr_matrix <- assays(traj8_RNA_Bap1KO)$counts
#### Ensure the pseudotime values are named with the same cell names as the expression matrix columns
names(pseudotime) <- colnames(expr_matrix)
#### Function to find the peak pseudotime for each gene (raw and smoothed)
find_max_pseudotime <- function(gene_expr, pseudotime) {
  # Raw peak pseudotime
  raw_peak_pseudotime <- pseudotime[which.max(gene_expr)]
  # Smooth gene expression using loess
  smooth_model <- loess(gene_expr ~ pseudotime)
  smooth_expr <- predict(smooth_model)
  # Smooth peak pseudotime
  smooth_peak_pseudotime <- pseudotime[which.max(smooth_expr)]
  return(list(raw_peak_pseudotime = raw_peak_pseudotime, 
              smooth_peak_pseudotime = smooth_peak_pseudotime))
}
#### Apply the function to all genes
peak_values <- apply(expr_matrix, 1, function(x) find_max_pseudotime(as.numeric(x), pseudotime))
#### Convert the results to a data frame
peak_df <- data.frame(
  gene = rownames(expr_matrix),
  raw_peak_pseudotime = sapply(peak_values, `[[`, "raw_peak_pseudotime"),
  smooth_peak_pseudotime = sapply(peak_values, `[[`, "smooth_peak_pseudotime")
) %>% as_tibble()


# save output: write.table(peak_df, file = c("output/condiments/traj8_RNA_Bap1KO_ActivationPoint.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## heatmap activate/induced genes along pseudotime
### DEG Start End
pseudotime_start_end_association # filter log2fc >0 >1
pseudotime_start_end_association = read_tsv("output/condiments/pseudotime_start_end_association_traj8_RNA_Bap1KO.txt")
pseudotime_start_end_association_logFC0 = pseudotime_start_end_association %>% 
  filter(fdr <0.05) %>%
  dplyr::select(gene) %>%
  unique()

pdf("output/condiments/heatmap_pseudotime_start_end_association_traj8_RNA_Bap1KO_fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj8_RNA_Bap1KO, gene = pseudotime_start_end_association_logFC0$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG time course
pseudotime_association = read_tsv("output/condiments/pseudotime_association_traj8_RNA_Bap1KO.txt")
pseudotime_association_deg = pseudotime_association %>%
  filter(fdr <0.05)%>%
  dplyr::select(gene) %>%
  unique()


pdf("output/condiments/heatmap_pseudotime_association_traj8_RNA_Bap1KO_fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj8_RNA_Bap1KO, gene = pseudotime_association_deg$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()



## Identify which pseudotime value corespond to which cluster ################
pseudotime <- colData(traj8_RNA_Bap1KO)$crv$pseudotime
sce_cells <- colnames(traj8_RNA_Bap1KO)
subset_seurat <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, cells = sce_cells) # Subset cell from traj2
clusters <- subset_seurat$cluster.annot # Extract cluster information
### Combine pseudotime and cluster information into a data frame
pseudotime_cluster_df <- data.frame(
  cell = colnames(traj8_RNA_Bap1KO),
  pseudotime = pseudotime,
  cluster = clusters
)  %>%
  arrange(pseudotime)

switch_points <- which(diff(as.numeric(factor(pseudotime_cluster_df$cluster))) != 0) # Find the indices where the cluster changes
switch_pseudotimes <- pseudotime_cluster_df$pseudotime[switch_points] # Extract the pseudotime values at these switch points
switch_clusters_from <- pseudotime_cluster_df$cluster[switch_points]
switch_clusters_to <- pseudotime_cluster_df$cluster[switch_points + 1]
switch_df <- data.frame(
  switch_pseudotime = switch_pseudotimes,
  cluster_from = switch_clusters_from,
  cluster_to = switch_clusters_to 
) %>%
  group_by(cluster_from, cluster_to) %>%
  summarize(median_switch_pseudotime = median(switch_pseudotime), .groups = 'drop')
write.table(switch_df, file = c("output/condiments/switch_df_traj8_RNA_Bap1KO.txt"),sep="\t", quote=FALSE, row.names=FALSE)
##################################

















################### Time Course effect (BAP1KO CONDITION) - TRAJECTORY 2 - ##############################################
# NSC --> CA123
set.seed(42)
traj2_RNA_Bap1KO <- readRDS("output/condiments/traj2_RNA_Bap1KO.rds")
## Genes that change with pseudotime

pseudotime_association <- associationTest(traj2_RNA_Bap1KO) # statistical test to check whether gene expression is constant across pseudotime within a lineage
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$fdr), ]
pseudotime_association$gene <- rownames(pseudotime_association)
pseudotime_association = as_tibble(pseudotime_association) 
# save output: write.table(pseudotime_association, file = c("output/condiments/pseudotime_association_traj2_RNA_Bap1KO.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## Genes that change between two pseudotime points (start vs end)
pseudotime_start_end_association <- startVsEndTest(traj2_RNA_Bap1KO, pseudotimeValues = NULL)
pseudotime_start_end_association$gene <- rownames(pseudotime_start_end_association)
pseudotime_start_end_association$fdr <- p.adjust(pseudotime_start_end_association$pvalue, method = "fdr")
pseudotime_start_end_association <- pseudotime_start_end_association[order(pseudotime_start_end_association$fdr), ]
##--> log2FC = end - start: negative log2fc means start point higher average expr than end point
# save output: write.table(pseudotime_start_end_association, file = c("output/condiments/pseudotime_start_end_association_traj2_RNA_Bap1KO.txt"),sep="\t", quote=FALSE, row.names=FALSE)


sce_cells <- colnames(traj2_RNA_Bap1KO) # collect cells of traj1
subset_traj2_RNA_Bap1KO_multiome_WT_Bap1KO_QCV2vC1.sct <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, cells = sce_cells) # Create a seurat object with only cells from traj1



### plot top 25 genes
#### Select the top 25 genes with positive logFC and the top 25 with negative logFC
top25_posFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 > 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

top25_negFC_genes <- pseudotime_start_end_association %>%
  filter(fdr < 0.05, logFClineage1 < 0) %>%
  top_n(25, waldStat) %>%
  arrange(desc(waldStat)) %>%
  select(gene, logFClineage1)

## Function to plot and save in PDF
plot_and_save <- function(genes_df, file_name) {
  pdf(file_name, width=5, height=4)
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$gene[i]
    logFC <- genes_df$logFClineage1[i]
    plot_title <- paste0(gene, " (logFC: ", round(logFC, 2), ")")
    p <- plotSmoothers(traj2_RNA_Bap1KO, subset_traj2_RNA_Bap1KO_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts, gene = gene)
    p <- p + ggtitle(plot_title)
    print(p)
  }
  dev.off()
}

# Generate PDFs
plot_and_save(top25_posFC_genes, "output/condiments/plotSmoothers-top25_posFC_traj2_RNA_Bap1KO.pdf")
plot_and_save(top25_negFC_genes, "output/condiments/plotSmoothers-top25_negFC_traj2_RNA_Bap1KO.pdf")


### plot unique genes (Pax6, Tbr2/Eomes, Prox1)
pdf("output/condiments/plotSmoothers-Eomes-traj2_RNA_Bap1KO.pdf",  width=5, height=4)
plotSmoothers(traj2_RNA_Bap1KO, subset_traj2_RNA_Bap1KO_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts, gene = "Eomes" )
dev.off()


### plot average of multiple genes - multiple gene lists
#### Define the gene list
gene_list_1 <- c("Pax6")
gene_list_2 <- c("Eomes")
gene_list_3 <- c("Prox1")
# Extract the counts for the genes of interest
counts_matrix <- subset_traj2_RNA_Bap1KO_multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts
# Subset the data for each gene list and calculate the average expression
subset_data_1 <- counts_matrix[gene_list_1, , drop = FALSE]
subset_data_2 <- counts_matrix[gene_list_2, , drop = FALSE]
subset_data_3 <- counts_matrix[gene_list_3, , drop = FALSE]
# Extract pseudotime
pseudotime_vector <- colData(traj2_RNA_Bap1KO)$crv
# Create data frames for plotting
plot_data_1 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_1)), GeneList = "GeneList 1")
plot_data_2 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_2)), GeneList = "GeneList 2")
plot_data_3 <- data.frame(Pseudotime = pseudotime_vector, Expression = log1p(as.numeric(subset_data_3)), GeneList = "GeneList 3")

# Combine the data frames
plot_data <- rbind(plot_data_1, plot_data_2, plot_data_3)

plot_data <- plot_data %>%
  group_by(GeneList) %>%
  mutate(SmoothedExpression = predict(loess(Expression ~ Pseudotime.pseudotime, span = 0.5)))


# Combine the smoothed data frames for plotting
smoothed_plot_data <- rbind(plot_data_1, plot_data_2, plot_data_3)
# Plot the smoothed trajectories
pdf("output/condiments/plotSmoothers-Pax6EomesProx1-traj2_RNA_Bap1KO.pdf", width = 5, height = 4)
ggplot(plot_data, aes(x = Pseudotime.pseudotime, y = Expression, color = GeneList)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = SmoothedExpression), size = 1) +
  scale_color_manual(values = c("GeneList 1" = "green", "GeneList 2" = "blue", "GeneList 3" = "purple")) + # Custom colors
  labs(
    x = "Pseudotime",
    y = "Log(Expression + 1)",
    color = "Gene List"
  ) +
  theme_bw()
dev.off()






## Identify Activation point = peak (maximum expression) of each gene along the pseudotime trajectory
### Identify peak of expression (max expr) of these Time-course DEG
traj2_RNA_Bap1KO
#### Extract pseudotime values
pseudotime <- colData(traj2_RNA_Bap1KO)$crv$pseudotime
#### Extract the expression matrix
expr_matrix <- assays(traj2_RNA_Bap1KO)$counts
#### Ensure the pseudotime values are named with the same cell names as the expression matrix columns
names(pseudotime) <- colnames(expr_matrix)
#### Function to find the peak pseudotime for each gene (raw and smoothed)
find_max_pseudotime <- function(gene_expr, pseudotime) {
  # Raw peak pseudotime
  raw_peak_pseudotime <- pseudotime[which.max(gene_expr)]
  # Smooth gene expression using loess
  smooth_model <- loess(gene_expr ~ pseudotime)
  smooth_expr <- predict(smooth_model)
  # Smooth peak pseudotime
  smooth_peak_pseudotime <- pseudotime[which.max(smooth_expr)]
  return(list(raw_peak_pseudotime = raw_peak_pseudotime, 
              smooth_peak_pseudotime = smooth_peak_pseudotime))
}
#### Apply the function to all genes
peak_values <- apply(expr_matrix, 1, function(x) find_max_pseudotime(as.numeric(x), pseudotime))
#### Convert the results to a data frame
peak_df <- data.frame(
  gene = rownames(expr_matrix),
  raw_peak_pseudotime = sapply(peak_values, `[[`, "raw_peak_pseudotime"),
  smooth_peak_pseudotime = sapply(peak_values, `[[`, "smooth_peak_pseudotime")
) %>% as_tibble()


# save output: write.table(peak_df, file = c("output/condiments/traj2_RNA_Bap1KO_ActivationPoint.txt"),sep="\t", quote=FALSE, row.names=FALSE)


## heatmap activate/induced genes along pseudotime
### DEG Start End
pseudotime_start_end_association # filter log2fc >0 >1
pseudotime_start_end_association = read_tsv("output/condiments/pseudotime_start_end_association_traj2_RNA_Bap1KO.txt")
pseudotime_start_end_association_logFC0 = pseudotime_start_end_association %>% 
  filter(fdr <0.05) %>%
  dplyr::select(gene) %>%
  unique()

pdf("output/condiments/heatmap_pseudotime_start_end_association_traj2_RNA_Bap1KO_fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj2_RNA_Bap1KO, gene = pseudotime_start_end_association_logFC0$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()

### DEG time course
pseudotime_association = read_tsv("output/condiments/pseudotime_association_traj2_RNA_Bap1KO.txt")
pseudotime_association_deg = pseudotime_association %>%
  filter(fdr <0.05)%>%
  dplyr::select(gene) %>%
  unique()


pdf("output/condiments/heatmap_pseudotime_association_traj2_RNA_Bap1KO_fdr05.pdf", width=8, height=10)
yhatSmooth <- predictSmooth(traj2_RNA_Bap1KO, gene = pseudotime_association_deg$gene, nPoints = 25, tidy = FALSE)
yhatSmooth <- yhatSmooth[order(apply(yhatSmooth,1,which.max)), ]
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:25]))),
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
dev.off()



## Identify which pseudotime value corespond to which cluster ################
pseudotime <- colData(traj2_RNA_Bap1KO)$crv$pseudotime
sce_cells <- colnames(traj2_RNA_Bap1KO)
subset_seurat <- subset(multiome_WT_Bap1KO_QCV2vC1.sct, cells = sce_cells) # Subset cell from traj2
clusters <- subset_seurat$cluster.annot # Extract cluster information
### Combine pseudotime and cluster information into a data frame
pseudotime_cluster_df <- data.frame(
  cell = colnames(traj2_RNA_Bap1KO),
  pseudotime = pseudotime,
  cluster = clusters
)  %>%
  arrange(pseudotime)

switch_points <- which(diff(as.numeric(factor(pseudotime_cluster_df$cluster))) != 0) # Find the indices where the cluster changes
switch_pseudotimes <- pseudotime_cluster_df$pseudotime[switch_points] # Extract the pseudotime values at these switch points
switch_clusters_from <- pseudotime_cluster_df$cluster[switch_points]
switch_clusters_to <- pseudotime_cluster_df$cluster[switch_points + 1]
switch_df <- data.frame(
  switch_pseudotime = switch_pseudotimes,
  cluster_from = switch_clusters_from,
  cluster_to = switch_clusters_to 
) %>%
  group_by(cluster_from, cluster_to) %>%
  summarize(median_switch_pseudotime = median(switch_pseudotime), .groups = 'drop')
write.table(switch_df, file = c("output/condiments/switch_df_traj2_RNA_Bap1KO.txt"),sep="\t", quote=FALSE, row.names=FALSE)
##################################


```




#### Run condiments RNA assay - V3 NSC prol2 --> DG_GC (Task post meeting 20250327)

Starting cluster = NSC_proliferative_2; express Hopx genes / end cluster = DG_GC


Lets **isolate the cell types** first; follow same method as Kcnc1 project (`002*/005*`)




```bash
conda activate condiments_Signac
```



```R
# package installation 
## install.packages("remotes")
## remotes::install_github("cran/spatstat.core")
## remotes::install_version("Seurat", "4.0.3")
## install.packages("magrittr")
## install.packages("magrittr")
## install.packages("dplyr")
## BiocManager::install("DelayedMatrixStats")
## BiocManager::install("tradeSeq")


# packages
library("condiments")
library("Seurat")
library("magrittr") # to use pipe
library("dplyr") # to use bind_cols and sample_frac
library("SingleCellExperiment") # for reducedDims
library("ggplot2")
library("slingshot")
library("DelayedMatrixStats")
library("tidyr")
library("tradeSeq") 
library("cowplot")
library("scales")
library("pheatmap")

# Data import


multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds")


DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA" # According to condiments workflow


# convert to SingleCellExperiment
WT_Bap1KO <- as.SingleCellExperiment(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA")



## SEPARATE CELLS for each trajectory ############################

########################################################
############################ DG_GC ############################
########################################################

# First filter based on cell type
Part_DG_GC <- WT_Bap1KO[, WT_Bap1KO$cluster.annot %in% c("NSC_proliferative_2", "IP", "DG_GC")]
table(Part_DG_GC$cluster.annot) # to double check


# tidy
df <- bind_cols(
  as.data.frame(reducedDims(Part_DG_GC)$UMAP),
  as.data.frame(colData(Part_DG_GC)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
pdf("output/condiments/UMAP_WT_Bap1KO-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-Part_DG_GC.pdf", width=6, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col = cluster.annot)) +
  geom_point(size = .7) +
  labs(col = "orig.ident") +
  theme_classic()
dev.off()


## Second filter based on UMAP coordinate
umap_coords <- reducedDims(Part_DG_GC)$UMAP

# Filter conditions based on your description:
# Keep cells with UMAP_1 > -3 and UMAP_2 < 2.5
selected_cells <- umap_coords[,1] < 4 &  umap_coords[,2] < -2.5 &  umap_coords[,2] > -9.5

# Subset your SCE object
Part_DG_GC_subset <- Part_DG_GC[, selected_cells]

# Check resulting subset
dim(Part_DG_GC_subset)

df <- bind_cols(
  as.data.frame(reducedDims(Part_DG_GC_subset)$UMAP),
  as.data.frame(colData(Part_DG_GC_subset)[, -3])
  ) %>%
  sample_frac(1)

pdf("output/condiments/UMAP_WT_Bap1KO-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-Part_DG_GC_subset.pdf", width=7, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col =  cluster.annot)) +
  geom_point(size = .7) +
  theme_classic()
dev.off()



## genotype overlap
pdf("output/condiments/UMAP_genotype_WT_Bap1KO-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-Part_DG_GC_subset.pdf", width=7, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col = orig.ident)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("blue", "red")) + # Specify colors here
  labs(col = "Genotype") +
  theme_classic()
dev.off()

## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(umap_1, umap_2) %>% as.matrix(), 
  conditions = df$orig.ident,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

pdf("output/condiments/UMAP_imbalance_score_WT_Bap1KO-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-Part_DG_GC_subset.pdf", width=7, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores") +
  theme_classic()
dev.off()





## PLOT with Separate trajectories
### Testing Area ############


##########################################



Part_DG_GC_subset <- slingshot(Part_DG_GC_subset, reducedDim = 'UMAP',
                 clusterLabels = colData(Part_DG_GC_subset)$cluster.annot,
                 start.clus = c("NSC_proliferative_2"), end.clus = c("DG_GC") ,approx_points = 100, extend = 'pc1', stretch = 1)




#test reduceDim PCA or subset endoderm
topologyTest(SlingshotDataSet(Part_DG_GC_subset), Part_DG_GC_subset$orig.ident) #  


sdss <- slingshot_conditions(SlingshotDataSet(Part_DG_GC_subset), Part_DG_GC_subset$orig.ident)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "orig.ident")



#  

pdf("output/condiments/UMAP_trajectory_separated_WT_Bap1KO-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-Part_DG_GC_subset-STARTNSCprol2_ENDDGGC_points100extendpc1stretch1.pdf", width=7, height=5)
ggplot(df, aes(x = umap_1, y = umap_2, col = orig.ident)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(orig.ident, Lineage, Order),
            aes(group = interaction(Lineage, orig.ident)), size = 1.5) +
  theme_classic()
dev.off()










## PLOT with separate trajectories - Individually
### WT
Part_DG_GC_subset_WT <- Part_DG_GC_subset[, Part_DG_GC_subset$orig.ident == "multiome_WT"]

df_2 <- bind_cols(
  as.data.frame(reducedDim(Part_DG_GC_subset_WT, "UMAP")),
  slingPseudotime(Part_DG_GC_subset_WT) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(Part_DG_GC_subset_WT) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst))
curves <- slingCurves(Part_DG_GC_subset_WT, as.df = TRUE)
### Function to create the plot for each lineage
create_plot <- function(lineage_number) {
  df_2 <- df_2 %>%
    mutate(pst = case_when(
      !!sym(paste0("Lineage", lineage_number, "_pst")) > 0 ~ !!sym(paste0("Lineage", lineage_number, "_pst")),
      TRUE ~ 0
    ),
    group = if_else(pst > 0, paste0("lineage", lineage_number), "other"))
  curves_filtered <- curves %>% filter(Lineage == lineage_number)
  curves_endpoints <- curves_filtered %>%
    group_by(Lineage) %>%
    arrange(Order) %>%
    top_n(1, Order) # Get the top/last ordered point for each group
  df_2_lineage <- df_2 %>% filter(group == paste0("lineage", lineage_number))
  df_2_other <- df_2 %>% filter(group == "other")
  p <- ggplot() +
    geom_point(data = df_2_other, aes(x = umap_1, y = umap_2), size = .7, color = "grey85") +
    geom_point(data = df_2_lineage, aes(x = umap_1, y = umap_2, col = pst), size = .7) +
    scale_color_viridis_c() +
    labs(col = "Pseudotime", title = paste("Lineage", lineage_number)) +
    geom_path(data = curves_filtered %>% arrange(Order),
              aes(x = umap_1, y = umap_2, group = Lineage), col = "black", size = 1) +
    geom_text(data = curves_endpoints, aes(x = umap_1, y = umap_2, label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
    theme_classic()
  return(p)
}
### Generate the plots for each lineage
plots <- list()
for (i in 1:1) {
  plots[[i]] <- create_plot(i)
}
pdf("output/condiments/UMAP_trajectory_common_label_Part_DG_GC_subset_WT-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-STARTNSCprol2_ENDDGGC_points100extendpc1stretch1_WTonly.pdf", width=7, height=5)
gridExtra::grid.arrange(grobs = plots, ncol = 1)
dev.off()


### Bap1KO
Part_DG_GC_subset_Bap1KO <- Part_DG_GC_subset[, Part_DG_GC_subset$orig.ident == "multiome_Bap1KO"]

df_2 <- bind_cols(
  as.data.frame(reducedDim(Part_DG_GC_subset_Bap1KO, "UMAP")),
  slingPseudotime(Part_DG_GC_subset_Bap1KO) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(Part_DG_GC_subset_Bap1KO) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst))
curves <- slingCurves(Part_DG_GC_subset_Bap1KO, as.df = TRUE)
### Function to create the plot for each lineage
create_plot <- function(lineage_number) {
  df_2 <- df_2 %>%
    mutate(pst = case_when(
      !!sym(paste0("Lineage", lineage_number, "_pst")) > 0 ~ !!sym(paste0("Lineage", lineage_number, "_pst")),
      TRUE ~ 0
    ),
    group = if_else(pst > 0, paste0("lineage", lineage_number), "other"))
  curves_filtered <- curves %>% filter(Lineage == lineage_number)
  curves_endpoints <- curves_filtered %>%
    group_by(Lineage) %>%
    arrange(Order) %>%
    top_n(1, Order) # Get the top/last ordered point for each group
  df_2_lineage <- df_2 %>% filter(group == paste0("lineage", lineage_number))
  df_2_other <- df_2 %>% filter(group == "other")
  p <- ggplot() +
    geom_point(data = df_2_other, aes(x = umap_1, y = umap_2), size = .7, color = "grey85") +
    geom_point(data = df_2_lineage, aes(x = umap_1, y = umap_2, col = pst), size = .7) +
    scale_color_viridis_c() +
    labs(col = "Pseudotime", title = paste("Lineage", lineage_number)) +
    geom_path(data = curves_filtered %>% arrange(Order),
              aes(x = umap_1, y = umap_2, group = Lineage), col = "black", size = 1) +
    geom_text(data = curves_endpoints, aes(x = umap_1, y = umap_2, label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
    theme_classic()
  return(p)
}
### Generate the plots for each lineage
plots <- list()
for (i in 1:1) {
  plots[[i]] <- create_plot(i)
}
pdf("output/condiments/UMAP_trajectory_common_label_Part_DG_GC_subset_WT-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-STARTNSCprol2_ENDDGGC_points100extendpc1stretch1_Bap1KOonly.pdf", width=7, height=5)
gridExtra::grid.arrange(grobs = plots, ncol = 1)
dev.off()



# Differential Progression
prog_res <- progressionTest(Part_DG_GC_subset, conditions = Part_DG_GC_subset$orig.ident, lineages = TRUE)

df_3 <-  slingPseudotime(Part_DG_GC_subset) %>% as.data.frame() 

df_3$orig.ident <- Part_DG_GC_subset$orig.ident
df_3 <- df_3 %>% 
  pivot_longer(-orig.ident, names_to = "Lineage",
               values_to = "pst") %>%
  filter(!is.na(pst))

pdf("output/condiments/densityPlot_trajectory_lineage_Part_DG_GC_subset-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-STARTNSCprol2_ENDDGGC_points100extendpc1stretch1.pdf", width=6, height=3)
ggplot(df_3, aes(x = pst)) +
  geom_density(alpha = .8, aes(fill = orig.ident), col = "transparent") +
  geom_density(aes(col = orig.ident), fill = "transparent", size = 1.5) +
  labs(x = "Pseudotime", fill = "orig.ident") +
  facet_wrap(~Lineage, scales = "free", nrow=2) +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("blue", "red"))
  )) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()
dev.off()


#### ->  save.image(file="output/condiments/condiments-Part_DG_GC_subset_STARTNSCprol2_ENDDGGC_points100extendpc1stretch1-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.RData")
### load("output/condiments/condiments-Part_DG_GC_subset_STARTNSCprol2_ENDDGGC_points100extendpc1stretch1-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.RData")
set.seed(42)

#  Differential expression
# --> Run fitGam() through Slurm


################### Time Course effect COMMON CONDITIONS ######################################################
## TRAJECTORY1 - DG_GC ##################
set.seed(42)
traj1 <- readRDS("output/condiments/traj1_Part_DG_GC_subset-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.rds")



## DEGs between condition
traj1_l2fc0 <- conditionTest(traj1, l2fc = 0) # 


# Correct the pvalue with fdr
traj1_l2fc0$padj <- p.adjust(traj1_l2fc0$pvalue, "fdr")


### Save output tables
traj1_l2fc0$gene <- rownames(traj1_l2fc0) # create new column label gene; as matrix before
condRes_traj1_l2fc0 <- traj1_l2fc0[, c(ncol(traj1_l2fc0), 1:(ncol(traj1_l2fc0)-1))] # just to put gene column 1st
write.table(condRes_traj1_l2fc0, file = c("output/condiments/condRes-traj1_DG_GC-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0.txt"),sep="\t", quote=FALSE, row.names=FALSE)

# Heatmap clutering DEGs per traj _ REVISED METHOD
## import DEGs
condRes_traj1_l2fc0 <- read.table("output/condiments/condRes-traj1_DG_GC-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) 

## Isolate significant DEGs and transform into a vector
conditionGenes_traj1_l2fc0 <- condRes_traj1_l2fc0 %>% 
  filter(padj <= 0.05) %>%
  pull(gene)

# Predict smoothed values
yhatSmooth <- 
  predictSmooth(traj1, gene = conditionGenes_traj1_l2fc0, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth, 1, scales::rescale))
combinedData <- yhatSmoothScaled[, c(51:100, 1:50)]
# Generate heatmap with clustering
# Perform hierarchical clustering
hc <- hclust(dist(combinedData))
clusters <- cutree(hc, k=6) # !!!!!!!!!!!!!!!!!! CHANGE CLUSTER NB HERE !!!!!!!!!!!!!!!!!!
# Create an annotation data frame for the rows based on cluster assignments
annotation_row <- data.frame(Cluster = factor(clusters))


# Line plots
library("reshape2")
library("stringr")
# Assuming yhatSmoothScaled contains your smoothed gene expression data
# Convert the yhatSmoothScaled data to a dataframe
df <- as.data.frame(yhatSmoothScaled)
df$Gene <- rownames(df)
# Transform the data into a long format
df_long <- melt(df, id.vars = "Gene", variable.name = "Pseudotime", value.name = "Expression")
# Attach the cluster information to the data frame
df$Cluster <- factor(clusters[df$Gene])
df_long$Cluster <- df$Cluster[match(df_long$Gene, df$Gene)]

# Extract condition column
df_long$Condition <- ifelse(str_detect(df_long$Pseudotime, "WT"), "WT", "Bap1KO")

# Extract the point value and convert it to numeric
df_long$Updated_Pseudotime <- as.numeric(str_extract(df_long$Pseudotime, "(?<=point)\\d+"))

# Define colors for the conditions
color_map <- c("WT" = "blue", "Bap1KO" = "red")

gene_counts <- df_long %>%
  group_by(Cluster) %>%
  summarise(GeneCount = n_distinct(Gene))
df_long <- df_long %>%
  left_join(gene_counts, by = "Cluster") %>%
  mutate(ClusterLabel = paste0("Cluster ", Cluster, " (", GeneCount, " genes)"))

# Plot using ggplot
#pdf("output/condiments/clustered_linePlot_traj1_DG_GC-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0-cl10.pdf", width=10, height=5)
pdf("output/condiments/clustered_linePlot_traj1_DG_GC-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0-cl6.pdf", width=12, height=5)
ggplot(df_long, aes(x = as.numeric(Updated_Pseudotime), y = Expression, group = Gene)) + 
  geom_line(data = subset(df_long, Condition == "WT"), aes(color = Condition), alpha = 0.5) +
  geom_line(data = subset(df_long, Condition == "Bap1KO"), aes(color = Condition), alpha = 0.5) +
  scale_color_manual(values = color_map) + 
  facet_wrap(~ClusterLabel, scales = "free_y", nrow = 2) +  # Use the updated ClusterLabel column
  theme_bw() +
  labs(title = "Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")
dev.off()

# Plot using ggplot
pdf("output/condiments/smoothed_linePlot_traj1_DG_GC-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0-cl6.pdf", width=10, height=5)
ggplot(df_long, aes(x = Updated_Pseudotime, y = Expression, color = Condition)) + 
  geom_smooth(method = "loess", se = TRUE, span = 0.5) + 
  scale_color_manual(values = color_map) + 
  facet_wrap(~ClusterLabel, scales = "free_y", nrow = 2) +  # Use the updated ClusterLabel column
  theme_bw() +
  labs(title = "Smoothed Gene Expression Dynamics Across Pseudotime by Cluster",
       x = "Pseudotime",
       y = "Expression Level")
dev.off()



### Export gene list from each cluster
## Create a data frame with gene names and their respective cluster assignments
output_df <- data.frame(
  gene = rownames(combinedData),
  cluster = clusters
)

# Write the data frame to a .txt file
write.table(output_df, 
            file = "output/condiments/gene_clusters-traj1_DG_GC-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0-cl6.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)


# Check some genes individually - RNA
## FOR LINEAGE 1
counts <- multiome_WT_Bap1KO_QCV2vC1.sct[["RNA"]]@counts # Collect the counts from seurat
cond <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident) # identify conditions
pseudotimes <- slingPseudotime(Part_DG_GC_subset, na = FALSE) [,1] # HERE INDICATE TRAJ
cellweights <- slingCurveWeights(Part_DG_GC_subset) [,1] # HERE INDICATE TRAJ
#### Subset the counts, pseudotimes, and cell weights for non-zero weights:
sub_weights <- cellweights[cellweights != 0]
sub_pseudotimes <- pseudotimes[names(pseudotimes) %in% names(sub_weights)]
sub_counts <- counts[, colnames(counts) %in% names(sub_weights)]
sub_cond <- cond[colnames(counts) %in% names(sub_weights)]
pdf("output/condiments/plotSmoothers-traj1_DG_GC_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-RNA_common-Prox1.pdf", width=5, height=2)
plotSmoothers(traj1, sub_counts, gene = "Prox1", curvesCol = c("blue","red"), lwd = 1, size = 0.2) +
scale_color_manual(values =c("blue","red"))
dev.off()




# Heatmap representation
# Define colors for each cluster
# 20
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen", "grey", "darkred", "darkblue", "gold", "darkgray", "lightblue", "lightgreen", "lightcoral", "lightpink", "lightcyan"))(20),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 10
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown", "cyan", "darkgreen" ))(10),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 8
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink", "brown" ))(8),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 7
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange", "pink" ))(7),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 6
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple", "orange"))(6),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 5
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow", "purple"))(5),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# 4
cluster_colors <- setNames(colorRampPalette(c("red", "blue", "green", "yellow" ))(4),
                           unique(annotation_row$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Generate the heatmap

col_order <- order(grepl("WT", colnames(combinedData)), decreasing = TRUE)
combinedData <- combinedData[, col_order]
pdf("output/condiments/heatmap-traj1_DG_GC_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0_cl6.pdf", width=5, height=5)
pheatmap(combinedData,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  legend = TRUE,
  cutree_rows = 6,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors
)
dev.off()



```

--> 8 clusters seems optimal for pseudotime DEG DG_GC



##### GO pseudotime V3 NSC prol2 --> DG_GC (8 clusters)

Let's do GO analysis for all the 8 clusters of pseudotime DEGs


```bash
conda activate deseq2
```


```R
# Required packages
library("clusterProfiler")
library("org.Mm.eg.db")  
library("enrichplot")
library("tidyverse")
library("patchwork")


gene_clusters_traj1_DG_GC <- read.table("output/condiments/gene_clusters-traj1_DG_GC-dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0-cl8.txt", 
                            header = TRUE, 
                            sep = "\t", 
                            stringsAsFactors = FALSE)



# GO BP
pdf("output/Pathway/dotplot_BP-traj1_DG_GC_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0_cl8.pdf", width = 8, height = 6)
# Loop through clusters 1 to 10
for (cluster_id in sort(unique(gene_clusters_traj1_DG_GC$cluster))) {
  message("Processing cluster: ", cluster_id)
  gene_list <- gene_clusters_traj1_DG_GC %>%
    filter(cluster == cluster_id) %>%
    pull(gene)
  ego <- enrichGO(gene = gene_list,
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  readable = TRUE)
  if (!is.null(ego) && nrow(ego) > 0) {
    print(dotplot(ego, showCategory = 20) + ggtitle(paste("Cluster", cluster_id)))
  } else {
    print(ggplot() + ggtitle(paste("Cluster", cluster_id, "- No Enrichment")) + theme_void())
  }
}
dev.off()



# KEGG
pdf("output/Pathway/dotplot_KEGG-traj1_DG_GC_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0_cl8.pdf", width = 8, height = 6)
# Loop through clusters 1 to 10
for (cluster_id in sort(unique(gene_clusters_traj1_DG_GC$cluster))) {
  message("Processing KEGG cluster: ", cluster_id)
  gene_symbols <- gene_clusters_traj1_DG_GC %>%
    filter(cluster == cluster_id) %>%
    pull(gene)
  # Convert to ENTREZ IDs
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = gene_symbols,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first") %>%
    na.omit() %>% as.character()
  # Perform KEGG enrichment if gene list not empty
  ekegg <- if (length(entrez_ids) > 0) {
    enrichKEGG(gene = entrez_ids,
               organism = "mmu",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")
  } else { NULL }
  # Plot
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    print(dotplot(ekegg, showCategory = 20) + ggtitle(paste("Cluster", cluster_id)))
  } else {
    print(ggplot() + ggtitle(paste("Cluster", cluster_id, "- No KEGG Enrichment")) + theme_void())
  }
}
dev.off()




######## Specific case ################

genes_cluster <- gene_clusters_traj1_DG_GC %>%
  filter(cluster == 10) %>%
  pull(gene)

# Convert SYMBOLs to ENTREZ IDs
entrez_cluster <- mapIds(org.Mm.eg.db,
                          keys = genes_cluster,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
                          multiVals = "first") %>%
  na.omit() %>% as.character()

# KEGG enrichment
ekegg_cluster <- enrichKEGG(gene = entrez_cluster,
                             organism = "mmu",
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH")

# Plot
pdf("output/Pathway/dotplot_KEGG-traj1_DG_GC_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-l2fc0_cl8_cluster10.pdf", width = 6, height = 6)
if (!is.null(ekegg_cluster) && nrow(ekegg_cluster) > 0) {
  print(dotplot(ekegg_cluster, showCategory = 10) )
} else {
  print(ggplot() + ggtitle("No KEGG Enrichment") + theme_void())
}
dev.off()




```













#### FitGam slurm jobs



```bash
conda activate condiments_Signac

# trajectory per condition together (for DEG condition, condiments) - pseudotime-condition DEG
### traj of interest NSC --> DG_GC traj9
sbatch scripts/fitGAM_6knots_traj9_RNA_common.sh # 31299442 ok
### traj of interest NSC --> CA123 traj2
sbatch scripts/fitGAM_6knots_traj2_RNA_common.sh # 31299487 ok
### traj of interest NSC_prol_2 --> DG_GC isolated cells
sbatch scripts/fitGAM_6knots_traj1_RNA_common_DG_GC.sh # 45797474 ok





# trajectory per CONDITION SEP (all features, no parralelization) - pseudotime-dependent DEGs
## Following parameters: *StartNSCquiescentEndPyNsRSCULDGGCapprox100extendn*

### traj of interest NSC --> DG_GC: COMMON=9 WT=traj9; Bap1KO=traj8
sbatch scripts/fitGAM_6knots_traj9_RNA_WT.sh # 30439293 ok
sbatch scripts/fitGAM_6knots_traj8_RNA_Bap1KO.sh # 30458721 ok

### traj of interest NSC --> CA123: COMMON=2 WT=traj2; Bap1KO=traj2
sbatch scripts/fitGAM_6knots_traj2_RNA_WT.sh # 30439292 ok
sbatch scripts/fitGAM_6knots_traj2_RNA_Bap1KO.sh # 30458706 ok
```









#### Run condiments WNN reduction

```bash
conda activate condiments_Signac
```



```R
# package installation 
## install.packages("remotes")
## remotes::install_github("cran/spatstat.core")
## remotes::install_version("Seurat", "4.0.3")
## install.packages("magrittr")
## install.packages("magrittr")
## install.packages("dplyr")
## BiocManager::install("DelayedMatrixStats")
## BiocManager::install("tradeSeq")


# packages
library("condiments")
library("Seurat")
library("magrittr") # to use pipe
library("dplyr") # to use bind_cols and sample_frac
library("SingleCellExperiment") # for reducedDims
library("ggplot2")
library("slingshot")
library("DelayedMatrixStats")
library("tidyr")
# library("tradeSeq") # NOT INSTALL MAY NEED TO IF BUG!! I ddi not try yet
library("cowplot")
library("scales")
library("pheatmap")

# Data import

multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")


DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA" # According to condiments workflow


# convert to SingleCellExperiment
RNA_WT_Bap1KO <- as.SingleCellExperiment(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA")


# tidy
df <- bind_cols(
  as.data.frame(reducedDims(RNA_WT_Bap1KO)$WNN.UMAP),
  as.data.frame(colData(RNA_WT_Bap1KO)[, -3])
  ) %>%
  sample_frac(1)

# PLOT
## genotype overlap
pdf("output/condiments/UMAP_treatment_multiome_WT_Bap1KO_QCV2vC1_WNN.pdf", width=6, height=5)
ggplot(df, aes(x = wnnUMAP_1, y = wnnUMAP_2, col = orig.ident)) +
  geom_point(size = .7) +
  scale_color_manual(values = c("blue", "red")) + # Specify colors here
  labs(col = "Treatment") +
  theme_classic()
dev.off()

## imbalance score
scores <- condiments::imbalance_score(
  Object = df %>% select(wnnUMAP_1, wnnUMAP_2) %>% as.matrix(), 
  conditions = df$orig.ident,
  k = 20, smooth = 40)
df$scores <- scores$scaled_scores

pdf("output/condiments/UMAP_imbalance_score_multiome_WT_Bap1KO_QCV2vC1_WNN.pdf", width=5, height=5)
ggplot(df, aes(x = wnnUMAP_1, y = wnnUMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores") +
  theme_classic()
dev.off()


#  Trajectory Inference and Differential Topology
set.seed(42)

## PLOT with Separate trajectories
### Testing Area ############
humangastruloid2472hrs <- slingshot(humangastruloid2472hrs, reducedDim = 'UMAP',
                 clusterLabels = colData(humangastruloid2472hrs)$seurat_clusters,
                 start.clus = '3', approx_points = 100)

humangastruloid2472hrs <- slingshot(humangastruloid2472hrs, reducedDim = 'UMAP',
                 clusterLabels = colData(humangastruloid2472hrs)$seurat_clusters,
                 start.clus = '3', end.clus = c("8","1","4") ,approx_points = 100)

#                 extend = 'n', stretch = 0)


##########################################
RNA_WT_Bap1KO <- slingshot(RNA_WT_Bap1KO, reducedDim = 'UMAP',
                 clusterLabels = colData(RNA_WT_Bap1KO)$cluster.annot,
                 start.clus = 'NSC_quiescent', end.clus = c("PyNs_RSC_UL","DG_GC") ,approx_points = 100, extend = 'n')




#test reduceDim PCA or subset endoderm

set.seed(42)
topologyTest(SlingshotDataSet(RNA_WT_Bap1KO), RNA_WT_Bap1KO$orig.ident) #  


sdss <- slingshot_conditions(SlingshotDataSet(RNA_WT_Bap1KO), RNA_WT_Bap1KO$orig.ident)
curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "orig.ident")
# Rename columns in curves
colnames(curves)[which(colnames(curves) == "umap_1")] <- "wnnUMAP_1"
colnames(curves)[which(colnames(curves) == "umap_2")] <- "wnnUMAP_2"

#  
pdf("output/condiments/UMAP_trajectory_separated_multiome_WT_Bap1KO_QCV2vC1_STARTNSCquiescentENDPyNsRSCULDGGCextendN_WNN.pdf", width=6, height=5)
ggplot(df, aes(x = wnnUMAP_1, y = wnnUMAP_2, col = orig.ident)) +
  geom_point(size = .7, alpha = .2) +
  scale_color_brewer(palette = "Accent") +
  geom_path(data = curves %>% arrange(orig.ident, Lineage, Order),
            aes(group = interaction(Lineage, orig.ident)), size = 1.5) +
  theme_classic()
dev.off()


```



--> Not good to use WNN reduction for pseudotime. Pseudotime is temporal dynamics of gene expression. So focus on gene expression; so need top use the RNA UMAP version.

# Task post meeting 20250327

- Import correct1 .Rds file
- Identify unbiased marker genes in each cluster
- Generate Vln plot for Hopx (do statistics NSC_proliferative_2 vs all other clusters)
- Perform pseudotime analysis for DG_GC (starting cluster = NSC_proliferative_2; express Hopx genes / end cluster = DG_GC)



```bash
conda activate SignacV5
module load hdf5
```

```R
set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
library("ggbeeswarm")
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()
# remotes::install_github('immunogenomics/presto')

# import seurat object
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds")

multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident <- factor(multiome_WT_Bap1KO_QCV2vC1.sct$orig.ident, levels = c("multiome_WT", "multiome_Bap1KO"))



###############################################################
# VLN PLOTS with STATISTICS #####################
###############################################################


# Check some genes
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "RNA"

## perform DEG cluster vs all cluster
all_markers <- FindAllMarkers(multiome_WT_Bap1KO_QCV2vC1.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.01)
### output file
write.table(all_markers, file = "output/Signac/srat_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# map marker gene to its specific cluster
gene_cluster_map <- c(
  Hopx = "NSC_proliferative_2"
)
gene_cluster_df <- tibble::tibble(
  gene = names(gene_cluster_map),
  cluster = unname(gene_cluster_map)
)

# collect padj value
gene_cluster_df %>%
  left_join(all_markers, by = c("gene", "cluster")) 

# Prepare the data to annotate each gene with its p_val_adj
sig_data <- gene_cluster_df %>%
  left_join(all_markers, by = c("gene", "cluster")) %>%
  mutate(
    Identity = cluster,
    y_pos = avg_log2FC + 0.25,  # adjust based on expression range
    significance = ifelse(!is.na(p_val_adj), paste0("p.adj = ", signif(p_val_adj, 2)), "not found")
  ) %>%
  dplyr::select(gene, Identity, y_pos, significance, p_val_adj)

pdf("output/Signac/VlnPlot_RNA_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-Hopx_RNA.pdf", width=7, height=3.5)
for (gene in names(gene_cluster_map)) {
  print(paste("Generating plot for:", gene))
  gene_sig_data <- sig_data %>% dplyr::filter(gene == !!gene)
  # Define title with gene and p.adj
  if (nrow(gene_sig_data) == 1 && !is.na(gene_sig_data$p_val_adj)) {
    title_text <- paste0(gene, " (p.adj = ", signif(gene_sig_data$p_val_adj, 2), ")")
  } else {
    title_text <- paste0(gene, " (p.adj = not found)")
  }
  # Generate violin plot
  p <- VlnPlot(multiome_WT_Bap1KO_QCV2vC1.sct, 
               features = gene, 
               pt.size = 0.1) +
    ggtitle(title_text) +
    theme(plot.title = element_text(size=10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(p)
}
dev.off()


pdf("output/Signac/VlnPlot_RNA_multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks-Hopx_SCT.pdf", width=7, height=3.5)
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1.sct) <- "SCT"
for (gene in names(gene_cluster_map)) {
  print(paste("Generating plot for:", gene))
  gene_sig_data <- sig_data %>% dplyr::filter(gene == !!gene)
  # Define title with gene and p.adj
  if (nrow(gene_sig_data) == 1 && !is.na(gene_sig_data$p_val_adj)) {
    title_text <- paste0(gene, " (p.adj = ", signif(gene_sig_data$p_val_adj, 2), ")")
  } else {
    title_text <- paste0(gene, " (p.adj = not found)")
  }
  # Generate violin plot
  p <- VlnPlot(multiome_WT_Bap1KO_QCV2vC1.sct, 
               features = gene, 
               pt.size = 0.1) +
    ggtitle(title_text) +
    theme(plot.title = element_text(size=10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(p)
}
dev.off()


```



# GRN with Pando

Gene regulatory network using chromatin accessibility and gene expression data. Follow Pando guideline [here](https://quadbio.github.io/Pando/index.html).

## install Pando


```bash
# clone SignacV5 and install Pando in R
conda create --name Signac_Pando --clone SignacV5 # devtools::install_github('quadbiolab/Pando')
```



--> Work!!!


## Run Pando

### All cells

```bash
conda activate Signac_Pando
module load hdf5
```



```R
set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()
library("Pando")
library("chromVARmotifs") # devtools::install_github("GreenleafLab/chromVARmotifs")



# import and rename Seurat obj
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")


# Get gene names and filter out 'Rik' genes
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1_GRN.sct) <- "RNA" 
genes <- rownames(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
genes_filtered <- genes[!grepl("Rik", genes)]

# Find variable features
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1_GRN.sct) <- "RNA" 
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- FindVariableFeatures(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, 
                                                            selection.method = "vst", 
                                                            nfeatures = 3000)

# Get variable features and filter out 'Rik' genes
variable_genes <- multiome_WT_Bap1KO_QCV2vC1_GRN.sct[["RNA"]]@var.features
variable_genes_filtered <- variable_genes[!grepl("Rik", variable_genes)]
#--> 2859 genes instead of 3000

# import the phastConsElements60way.mm10 bed file and convert it to GR range
phastConsElements60way.mm10 = read_tsv("../meta/phastConsElements60way.mm10", col_names  = FALSE) %>%
  dplyr::rename("chr" = "X1", "start" = "X2", "end" = "X3")
## Convert to GRange object
phastConsElements60way_GRanges =  GRanges(seqnames = phastConsElements60way.mm10$chr,
                          ranges = IRanges(start = phastConsElements60way.mm10$start, 
                                           end = phastConsElements60way.mm10$end),
                          score = phastConsElements60way.mm10$X5)

# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA",
  regions = phastConsElements60way_GRanges  # Optional but recommended, see notes
  )



# Scan for TF motifs
#################### FAIL ########################################
data(motifs)

multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN.sct,
    pfm = motifs,
    genome = BSgenome.Mmusculus.UCSC.mm10
)
# --> Long, 2hrs; shorter 30mnin if  using  phastConsElements60way_GRanges
# --> This fail as data(motifs) is for mouse!
#####################################################################


## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)


######################################################################################################
# Inferring the GRN - Variable features only ####################################################################
## Variable Feature
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = variable_genes_filtered 
)

# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

# some QC plots
pdf("output/Pando/plot_gof_VariableFeatures.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_VariableFeatures.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)

pdf("output/Pando/plot_network_graph_VariableFeatures.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_VariableFeatures.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, layout='fr')
dev.off()

## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN.Ezh2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN.TF.sct, tf='Ezh2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Ezh2-VariableFeatures.pdf", width=10, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN.Ezh2.sct, tf='Ezh2')
dev.off()



######################################################################################################
# Inferring the GRN - All genes ####################################################################
## 
multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered 
)
### Save the Seurat object with the inferred GRN
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct.rds")



# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

# some QC plots
pdf("output/Pando/plot_gof_allGenes.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_allGenes.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)

pdf("output/Pando/plot_network_graph_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, layout='fr')
dev.off()

## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN.Ezh2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN.TF.sct, tf='Ezh2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Ezh2-allGenes.pdf", width=10, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN.Ezh2.sct, tf='Ezh2')
dev.off()




```

*Notes:*
- Make sure to run variableFeature on RNA assay prior running GRN code!! [Issue](https://github.com/quadbio/Pando/issues/29)
- `initiate_grn()` optional but better to add `regions`.  --> constraining the set of peaks to more confident regions cuts down on runtime and makes the resulting GRN more robust; discuss [here](https://quadbio.github.io/Pando/articles/getting_started.html). Here is where I found the mice equivalent: https://support.bioconductor.org/p/96226/. Issue discuss [here](https://github.com/quadbio/Pando/issues/62):
  - *phastConsElements60way* mice element data downloaded [here](https://genome.ucsc.edu/cgi-bin/hgTables). `60 Vert. El` selected; adn output format as bed: `phastConsElements60way.mm10` --> File transfer to `002*/meta/` folder
- For mice fin_motifs() step, issue discuss [here](https://github.com/quadbio/Pando/issues/26). Need to use a different motifs database, with nomenclature adapted for mice.
  - motif database install with [chromVARmotifs](https://github.com/GreenleafLab/chromVARmotifs)
- For mice need to remove all `*Rik` genes otherwise run into a [bug](https://github.com/quadbio/Pando/issues/1) at `infer_grn()`
  - I try subset seurat object ; remove Rik genes, but it fail, ATAC assay disapear... So instead I will mention to use non Rik genes at `infer_grn()`
  - I remove the `*Rik` genes only from the variable features, as per default the `infer_grn()` run on the variable feature RNA
- Not clear whether `infer_grn()` to be run on all genes or not; I think yes and running on VariableFeateure is for computational time reason. Mention [here](https://quadbio.github.io/Pando/articles/getting_started.html)
- `GREAT` seems better than `Signac` for `peak_to_gene_method = ` at `infer_grn()`, do more [stuff](https://quadbio.github.io/Pando/articles/association.html)
  - Not for us, it gave very few Tf-gene candidate as compare to Signac. *Let's use Signac method*
- Another option to test is only to work with peaks correlated with other genes, identified after running LinkPeaks(). Recommended [here](https://quadbio.github.io/Pando/articles/regions.html). **But maybe our target of interest will not be present...!**


--> None Bap1 complex TF are found in the tf dataset after `infer_grn()`...


Let's re-run Pando, but subset UMAP WT; Bap1KO. Do GRN for Bap1 and compare GRN map between WT and Bap1KO.








### WT and Bap1KO separated with regions to phastConsElements60way

```bash
srun --mem=500g --cpus-per-task=8 --pty bash -l

conda activate Signac_Pando
module load hdf5
```



```R
set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()
library("Pando")
library("chromVARmotifs") # devtools::install_github("GreenleafLab/chromVARmotifs")

# to run job in parrallel
library(doParallel) # install.packages("doParallel")
registerDoParallel(8) # update nb of core here!!


# import and rename Seurat obj
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")


# Get gene names and filter out 'Rik' genes
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1_GRN.sct) <- "RNA" 
genes <- rownames(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
genes_filtered <- genes[!grepl("Rik", genes)]

# Separate Seurat into WT and Bap1KO

multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, subset = orig.ident == "multiome_WT")
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, subset = orig.ident == "multiome_Bap1KO")



##################################################################################
# WT ##################################################################################
##################################################################################

# import the phastConsElements60way.mm10 bed file and convert it to GR range
phastConsElements60way.mm10 = read_tsv("../meta/phastConsElements60way.mm10", col_names  = FALSE) %>%
  dplyr::rename("chr" = "X1", "start" = "X2", "end" = "X3")
## Convert to GRange object
phastConsElements60way_GRanges =  GRanges(seqnames = phastConsElements60way.mm10$chr,
                          ranges = IRanges(start = phastConsElements60way.mm10$start, 
                                           end = phastConsElements60way.mm10$end),
                          score = phastConsElements60way.mm10$X5)

# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA",
  regions = phastConsElements60way_GRanges  # Optional but recommended, see notes
  )



# Scan for TF motifs

## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)

######################################################################################################
# Inferring the GRN - All genes ####################################################################

# TESTING THE MODEL TO HAVE OUR TF of interest in the GRN ##################
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn()
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct.rds")
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); only 14 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 415 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_cvglmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'cv.glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 396 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_bagging_ridge.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'bagging_ridge', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 0 Tfs
# --> brms (library(brms) and xgb (library(xgboost) not install
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_bagging_ridge.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'bagging_ridge', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 0 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_bagging_cv.glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'cv.glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 17 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 18 Tfs

############################################################################################################

## PARAMETER TO  USE:
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T,
    tf_cor = 0.05 # default 0.1
)
#--> YES Foxk1, Foxk2, Yy1 



# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)



# some QC plots
pdf("output/Pando/plot_gof_allGenes_WT.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_allGenes_WT.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)

pdf("output/Pando/plot_network_graph_allGenes_WT.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_allGenes_WT.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, layout='fr')
dev.off()

## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Yy1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Yy1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Yy1-allGenes_WT.pdf", width=10, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Yy1.sct, tf='Yy1')
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Foxk1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk1-allGenes_WT.pdf", width=10, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk1.sct, tf='Foxk1')
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Foxk2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk2-allGenes_WT.pdf", width=10, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk2.sct, tf='Foxk2')
dev.off()



XXX BELOW NOT MOD

##################################################################################
# Bap1KO ##################################################################################
##################################################################################



# import the phastConsElements60way.mm10 bed file and convert it to GR range
phastConsElements60way.mm10 = read_tsv("../meta/phastConsElements60way.mm10", col_names  = FALSE) %>%
  dplyr::rename("chr" = "X1", "start" = "X2", "end" = "X3")
## Convert to GRange object
phastConsElements60way_GRanges =  GRanges(seqnames = phastConsElements60way.mm10$chr,
                          ranges = IRanges(start = phastConsElements60way.mm10$start, 
                                           end = phastConsElements60way.mm10$end),
                          score = phastConsElements60way.mm10$X5)

# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA",
  regions = phastConsElements60way_GRanges  # Optional but recommended, see notes
  )



# Scan for TF motifs



## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)





######################################################################################################
# Inferring the GRN - All genes ####################################################################
## 
multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered 
)



# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

# some QC plots
pdf("output/Pando/plot_gof_allGenes.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_allGenes.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)

pdf("output/Pando/plot_network_graph_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, layout='fr')
dev.off()

## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN.Ezh2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN.TF.sct, tf='Ezh2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Ezh2-allGenes.pdf", width=10, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN.Ezh2.sct, tf='Ezh2')
dev.off()

```



--> Foxk1, Foxk2, Yy1; hard to find signfic to construct GRN. Had to reduce `tf_cor` treshold to 0.05 (insted of 0.1) in `infer_grn()`
  --> Could play with `find_modules()` parameter, but  trehold already relaxed...
  --> Let's try not to restrict analysis to conserved peaks see if that improve






### WT and Bap1KO separated without regions

```bash
srun --mem=500g --cpus-per-task=8 --pty bash -l

conda activate Signac_Pando
module load hdf5
```


```R
set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()
library("Pando")
library("chromVARmotifs") # devtools::install_github("GreenleafLab/chromVARmotifs")

# to run job in parrallel
library(doParallel) # install.packages("doParallel")
registerDoParallel(8) # update nb of core here!!


# import and rename Seurat obj
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")


# Get gene names and filter out 'Rik' genes
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1_GRN.sct) <- "RNA" 
genes <- rownames(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
genes_filtered <- genes[!grepl("Rik", genes)]

# Separate Seurat into WT and Bap1KO

multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, subset = orig.ident == "multiome_WT")
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, subset = orig.ident == "multiome_Bap1KO")



##################################################################################
# WT ##################################################################################
##################################################################################

# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA"
  # regions = phastConsElements60way_GRanges  # Optional but recommended, see notes
  )



# Scan for TF motifs

## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)

######################################################################################################
# Inferring the GRN - All genes ####################################################################

# TESTING THE MODEL TO HAVE OUR TF of interest in the GRN ##################
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn()
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct.rds")
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); only 14 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 415 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_cvglmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'cv.glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 396 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_bagging_ridge.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'bagging_ridge', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 0 Tfs
# --> brms (library(brms) and xgb (library(xgboost) not install
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_bagging_ridge.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'bagging_ridge', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 0 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_bagging_cv.glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'cv.glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 17 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 18 Tfs

############################################################################################################

## PARAMETER TO  USE:
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T,
    tf_cor = 0.05 # default 0.1
)
#--> YES Foxk1, Foxk2, Yy1 ; 616 TFs
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_noRegions.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_noRegions.sct.rds")


# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)



# some QC plots
pdf("output/Pando/plot_gof_allGenes_WT_noRegions.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_allGenes_WT_noRegions.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)

pdf("output/Pando/plot_network_graph_allGenes_WT_noRegions.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_allGenes_WT_noRegions.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, layout='fr')
dev.off()

## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Yy1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Yy1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Yy1-allGenes_WT_noRegions.pdf", width=5, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Yy1.sct, tf='Yy1', circular=F, label_nodes = "all")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Foxk1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk1-allGenes_WT_noRegions.pdf", width=15, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk1.sct, tf='Foxk1', label_nodes = "tfs")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Foxk2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk2-allGenes_WT_noRegions.pdf", width=7, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk2.sct, tf='Foxk2', circular=F, label_nodes = "all")
dev.off()




##################################################################################
# Bap1KO ##################################################################################
##################################################################################




# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA"
  )



# Scan for TF motifs



## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)






######################################################################################################
# Inferring the GRN - All genes ####################################################################
## 
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T,
    tf_cor = 0.05 # default 0.1
)
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes_noRegions.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes_noRegions.sct.rds")



# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

# some QC plots
pdf("output/Pando/plot_gof_Bap1KO_allGenes.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_Bap1KO_allGenes.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)

pdf("output/Pando/plot_network_graph_Bap1KO_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_Bap1KO_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, layout='fr')
dev.off()


## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Yy1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Yy1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Yy1-allGenes_Bap1_noRegions.pdf", width=5, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Yy1.sct, tf='Yy1', circular=F, label_nodes = "all")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Foxk1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk1-allGenes_Bap1_noRegions.pdf", width=15, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk1.sct, tf='Foxk1', label_nodes = "tfs")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Foxk2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk2-allGenes_Bap1_noRegions.pdf", width=7, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk2.sct, tf='Foxk2', circular=F, label_nodes = "all")
dev.off()



```






### WT and Bap1KO separated with regions LinkPeaks() - gene peak correlated

As recommended [here](https://quadbio.github.io/Pando/articles/regions.html); let's use only peaks correlated with the expression of their nearby genes. --> Strict, but robust GRN

```bash
srun --mem=500g --cpus-per-task=8 --pty bash -l

conda activate Signac_Pando
module load hdf5
```


```R
set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()
library("Pando")
library("chromVARmotifs") # devtools::install_github("GreenleafLab/chromVARmotifs")

# to run job in parrallel
library(doParallel) # install.packages("doParallel")
registerDoParallel(8) # update nb of core here!!


# import and rename Seurat obj
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")

# Select LinkPeaks()
Links = as_tibble(Links(multiome_WT_Bap1KO_QCV2vC1_GRN.sct))
## Convert to GRange object
Links_GRanges =  GRanges(Links)

# Get gene names and filter out 'Rik' genes
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1_GRN.sct) <- "RNA" 
genes <- rownames(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
genes_filtered <- genes[!grepl("Rik", genes)]

# Separate Seurat into WT and Bap1KO

multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, subset = orig.ident == "multiome_WT")
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, subset = orig.ident == "multiome_Bap1KO")



##################################################################################
# WT ##################################################################################
##################################################################################

# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA",
  regions = Links_GRanges  # Optional but recommended, see notes
  )



# Scan for TF motifs

## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)

######################################################################################################
# Inferring the GRN - All genes ####################################################################

# TESTING THE MODEL TO HAVE OUR TF of interest in the GRN ##################
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> 508 TFs; No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn()
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct.rds")
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); only 14 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 415 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_cvglmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'cv.glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 396 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_bagging_ridge.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'bagging_ridge', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 0 Tfs
# --> brms (library(brms) and xgb (library(xgboost) not install
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_bagging_ridge.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'bagging_ridge', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 0 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_bagging_cv.glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'cv.glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 17 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 18 Tfs

############################################################################################################

## PARAMETER TO  USE:
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T,
    tf_cor = 0.05 # default 0.1
)
#--> YES Foxk1, Foxk2, Yy1 ; 617 TFs
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_RegionsLinkPeaks.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- readRDS(file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_RegionsLinkPeaks.sct.rds")



# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)



# some QC plots
pdf("output/Pando/plot_gof_allGenes_WT_RegionsLinkPeaks.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_allGenes_WT_RegionsLinkPeaks.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)

pdf("output/Pando/plot_network_graph_allGenes_WT_RegionsLinkPeaks.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_allGenes_WT_RegionsLinkPeaks.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, layout='fr')
dev.off()

## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Yy1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Yy1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Yy1-allGenes_WT_RegionsLinkPeaks.pdf", width=5, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Yy1.sct, tf='Yy1', circular=F, label_nodes = "all")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Foxk1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk1-allGenes_WT_RegionsLinkPeaks.pdf", width=15, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk1.sct, tf='Foxk1', label_nodes = "tfs")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Foxk2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk2-allGenes_WT_RegionsLinkPeaks.pdf", width=7, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk2.sct, tf='Foxk2', circular=F, label_nodes = "all")
dev.off()




##################################################################################
# Bap1KO ##################################################################################
##################################################################################

XXX BELOW NOT MODIFIED


# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA"
  )



# Scan for TF motifs



## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)






######################################################################################################
# Inferring the GRN - All genes ####################################################################
## 
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T,
    tf_cor = 0.05 # default 0.1
)
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes_noRegions.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes_noRegions.sct.rds")



# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

# some QC plots
pdf("output/Pando/plot_gof_Bap1KO_allGenes.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_Bap1KO_allGenes.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)

pdf("output/Pando/plot_network_graph_Bap1KO_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_Bap1KO_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, layout='fr')
dev.off()


## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Yy1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Yy1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Yy1-allGenes_Bap1_noRegions.pdf", width=5, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Yy1.sct, tf='Yy1', circular=F, label_nodes = "all")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Foxk1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk1-allGenes_Bap1_noRegions.pdf", width=15, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk1.sct, tf='Foxk1', label_nodes = "tfs")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Foxk2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk2-allGenes_Bap1_noRegions.pdf", width=7, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk2.sct, tf='Foxk2', circular=F, label_nodes = "all")
dev.off()



```



--> Not sure what to do with these GRN...






### WT and Bap1KO separated with regions LinkPeaks() - gene peak correlated - Check Pax6: ONLY DG_GC trajectory


XXXY HER E!!

As recommended [here](https://quadbio.github.io/Pando/articles/regions.html); let's use only peaks correlated with the expression of their nearby genes. --> Strict, but robust GRN.

--> Lets filter here to only keep cells from the  DG_GC trajectory; as objective is to check Pax6 gene

```bash
srun --mem=500g --cpus-per-task=8 --pty bash -l

conda activate Signac_Pando
module load hdf5
```

```R
set.seed(42)

# library
library("Signac")
library("Seurat")
#library("hdf5r") # need to reinstall it at each session... with install.packages("hdf5r")
library("tidyverse")
library("EnsDb.Mmusculus.v79") # mm10
library("reticulate") # needed to use FindClusters()
library("metap") # needed to use FindConservedMarkers()
use_python("~/anaconda3/envs/SignacV5/bin/python") # to specify which python to use... Needed for FindClusters()
library("Pando")
library("chromVARmotifs") # devtools::install_github("GreenleafLab/chromVARmotifs")

# to run job in parrallel
library(doParallel) # install.packages("doParallel")
registerDoParallel(8) # update nb of core here!!


# import and rename Seurat obj
#multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000GeneActivityLinkPeaks.sct_numeric_label.rds")
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- readRDS(file = "output/seurat/multiome_WT_Bap1KO_QCV2vC1_dim40kparam42res065algo4feat2000correct1GeneActivityLinkPeaks.sct_numeric_label.rds")




# Get gene names and filter out 'Rik' genes
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1_GRN.sct) <- "RNA" 
genes <- rownames(multiome_WT_Bap1KO_QCV2vC1_GRN.sct)
genes_filtered <- genes[!grepl("Rik", genes)]

# Find variable features
DefaultAssay(multiome_WT_Bap1KO_QCV2vC1_GRN.sct) <- "RNA" 
multiome_WT_Bap1KO_QCV2vC1_GRN.sct <- FindVariableFeatures(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, 
                                                            selection.method = "vst", 
                                                            nfeatures = 3000)

# Get variable features and filter out 'Rik' genes
variable_genes <- multiome_WT_Bap1KO_QCV2vC1_GRN.sct[["RNA"]]@var.features
variable_genes_filtered <- variable_genes[!grepl("Rik", variable_genes)]
#--> 2859 genes instead of 3000




# Select LinkPeaks()
Links = as_tibble(Links(multiome_WT_Bap1KO_QCV2vC1_GRN.sct))
## Convert to GRange object
Links_GRanges =  GRanges(Links)



# Only keep DG_GC trajectory cells
# Step 1: Filter based on cluster identity
Part_DG_GC <- subset(
  multiome_WT_Bap1KO_QCV2vC1_GRN.sct, 
  subset = cluster.annot %in% c("NSC_proliferative_2", "IP", "DG_GC")
)

# Step 2: Filter based on UMAP coordinates
umap_coords <- Embeddings(Part_DG_GC, "umap")

# Apply your custom UMAP boundaries
selected_cells <- which(
  umap_coords[,1] < 4 &
  umap_coords[,2] < -2.5 &
  umap_coords[,2] > -9.5
)

# Subset the Seurat object
Part_DG_GC_subset <- subset(Part_DG_GC, cells = colnames(Part_DG_GC)[selected_cells])

# Optional: Check dimensions
dim(Part_DG_GC_subset)



# Separate Seurat into WT and Bap1KO

multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, subset = orig.ident == "multiome_WT")
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- subset(multiome_WT_Bap1KO_QCV2vC1_GRN.sct, subset = orig.ident == "multiome_Bap1KO")



##################################################################################
# WT ##################################################################################
##################################################################################

# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA",
  regions = Links_GRanges  # Optional but recommended, see notes
  )



# Scan for TF motifs

## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)

######################################################################################################
# Inferring the GRN - All genes ####################################################################

# TESTING THE MODEL TO HAVE OUR TF of interest in the GRN ##################
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> 508 TFs; No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn()
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct.rds")
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); only 14 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 415 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_cvglmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'cv.glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 396 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_Signac_bagging_ridge.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'bagging_ridge', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 0 Tfs
# --> brms (library(brms) and xgb (library(xgboost) not install
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_bagging_ridge.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'bagging_ridge', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 0 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_bagging_cv.glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'cv.glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> xxx No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 17 Tfs
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_GREAT_glmnet.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'GREAT', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glmnet', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T 
)
#--> No Foxk1, Foxk2, Yy1 -----> TRY OTHER PARAMETER IN infer)grn(); 18 Tfs

############################################################################################################

## PARAMETER TO  USE:
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T,
    tf_cor = 0.05 # default 0.1
)
#--> YES Foxk1, Foxk2, Yy1 ; 617 TFs
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_RegionsLinkPeaks.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- readRDS(file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes_RegionsLinkPeaks.sct.rds")



# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)



# some QC plots
pdf("output/Pando/plot_gof_allGenes_WT_RegionsLinkPeaks.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_allGenes_WT_RegionsLinkPeaks.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)

pdf("output/Pando/plot_network_graph_allGenes_WT_RegionsLinkPeaks.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_allGenes_WT_RegionsLinkPeaks.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, layout='fr')
dev.off()

## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Yy1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Yy1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Yy1-allGenes_WT_RegionsLinkPeaks.pdf", width=5, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Yy1.sct, tf='Yy1', circular=F, label_nodes = "all")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Foxk1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk1-allGenes_WT_RegionsLinkPeaks.pdf", width=15, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk1.sct, tf='Foxk1', label_nodes = "tfs")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.TF.sct, tf='Foxk2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk2-allGenes_WT_RegionsLinkPeaks.pdf", width=7, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_WT_allGenes.Foxk2.sct, tf='Foxk2', circular=F, label_nodes = "all")
dev.off()




##################################################################################
# Bap1KO ##################################################################################
##################################################################################

XXX BELOW NOT MODIFIED


# create grn object
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- initiate_grn(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
  peak_assay = "ATAC",
  rna_assay = "RNA"
  )



# Scan for TF motifs



## Prep mouse TF motifs database
library(BSgenome.Mmusculus.UCSC.mm10)
data("mouse_pwms_v2")
x <- character()
for(i in 1:length(mouse_pwms_v2@listData)){
x[i] <- mouse_pwms_v2@listData[[i]]@name
}
motif2tf <- data.frame(motif = names(mouse_pwms_v2@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(mouse_pwms_v2@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>%
subset(gene_id != "XP" & gene_id != "NP")
mouse_pwms_v3 <- subset(mouse_pwms_v2, names(mouse_pwms_v2@listData) %in% motif2tf$motif)

## Find motifs
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct <- find_motifs(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
    pfm = mouse_pwms_v3,
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
)






######################################################################################################
# Inferring the GRN - All genes ####################################################################
## 
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- infer_grn(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO.sct,
    peak_to_gene_method = 'Signac', # or use 'GREAT' consider overlapping regulatory regions, lets keep Signac
    method = 'glm', # other model can be tested: ('glmnet', 'cv.glmnet', 'xgb')
    genes = genes_filtered,
    parallel = T,
    tf_cor = 0.05 # default 0.1
)
#saveRDS(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, file = "output/Pando/multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes_noRegions.sct.rds")
#load: multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- readRDS(file = "multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes_noRegions.sct.rds")



# Find modules = Genes regulated by each TF
GetNetwork(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
coef(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- find_modules(
    multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

# some QC plots
pdf("output/Pando/plot_gof_Bap1KO_allGenes.pdf", width=7, height=6)
plot_gof(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, point_size=3)
dev.off()
pdf("output/Pando/plot_module_metrics_Bap1KO_allGenes.pdf", width=7, height=4)
plot_module_metrics(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
dev.off()

# GRN plots
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)

pdf("output/Pando/plot_network_graph_Bap1KO_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct)
dev.off()
pdf("output/Pando/plot_network_graph_fr_Bap1KO_allGenes.pdf", width=10, height=10)
plot_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, layout='fr')
dev.off()


## 1 TF
multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct <- get_network_graph(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.sct, 
    graph_name = 'full_graph', 
    umap_method = 'none')



multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Yy1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Yy1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Yy1-allGenes_Bap1_noRegions.pdf", width=5, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Yy1.sct, tf='Yy1', circular=F, label_nodes = "all")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk1.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Foxk1', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk1-allGenes_Bap1_noRegions.pdf", width=15, height=10)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk1.sct, tf='Foxk1', label_nodes = "tfs")
dev.off()

multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk2.sct <- get_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.TF.sct, tf='Foxk2', graph='full_graph')
pdf("output/Pando/plot_tf_network-Foxk2-allGenes_Bap1_noRegions.pdf", width=7, height=2)
plot_tf_network(multiome_WT_Bap1KO_QCV2vC1_GRN_Bap1KO_allGenes.Foxk2.sct, tf='Foxk2', circular=F, label_nodes = "all")
dev.off()



```





