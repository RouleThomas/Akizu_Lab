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
sbatch scripts/cellranger_count_RNA_WT.sh # 20409628 fail chemistry; 20409800 ok
sbatch scripts/cellranger_count_RNA_Bap1KO.sh # 20409629 fail chemistry; 20409813 ok
sbatch scripts/cellranger_count_ATAC_WT.sh # 20409600 fail chemistry; 20409818 fail need cellranger-atac; 20414461 too long; re run more ressource (750go) 20631094 ok
sbatch scripts/cellranger_count_ATAC_Bap1KO.sh # 20409604 fail chemistry; 20409822 fail need cellranger-atac; 20414466 too long; re run more ressource 20631095 ok

## For option2
# --> To do if option1 fail
#XXX sbatch scripts/cellranger_count_RNA_ATAC_WT.sh # 
#XXX sbatch scripts/cellranger_count_RNA_ATAC_Bap1KO.sh #  

```


--> As sc data generated with a multiome kit:
    - need to add `--chemistry ARC-v1` in counting for the RNA; solution found [here](https://bioinformatics.stackexchange.com/questions/18186/10x-low-rate-of-correct-barcodes-was-observed-for-the-candidate-chemistry-choice)
    - need to use `cellranger-atac count` to count ATAC exp solely; solution found [here](https://kb.10xgenomics.com/hc/en-us/articles/360061165691-Can-I-analyze-only-the-ATAC-data-from-my-single-cell-multiome-experiment)



--> RNA samples count succesfully

--> ATAC samples very long to count (>3 days wit 400Go mem)






## RNA contamination and doublet detection
- doublet detection using [scrublet](https://github.com/swolock/scrublet) **on the filtered matrix**
- ambient RNA correction using `soupX` in R before generating the Seurat object

Too many samples so let's run scrublet in 6 separate bash jobs

```bash
conda deactivate # base environment needed

# Run doublet detection/scrublet genotype/time tissue (CB/CX) together

sbatch scripts/scrublet_RNA_WT.sh # 20606254 xxx
sbatch scripts/scrublet_RNA_Bap1KO.sh # 20606261 xxx



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




# Generate bigwig

Generate coverage bigwig files to check for Bap1 exon7 deletion in Bap1KO.

--> Let's just generate the bigwig from the `possorted_genome_bam.bam` generated by 10X cellranger count 



```bash
conda activate deeptools

# raw
sbatch scripts/bamtobigwig_RNA_WT.sh # 21960401 ok
sbatch scripts/bamtobigwig_RNA_Bap1KO.sh # 21960405 ok


# BPM norm (=TPM)
sbatch scripts/bamtobigwig_BPMnorm_RNA_WT.sh # 21963563 xxx
sbatch scripts/bamtobigwig_BPMnorm_RNA_Bap1KO.sh # 21963630 xxx


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


```bash
conda activate SignacV5
module load hdf5_18/1.8.21 # to read `.h5` files
```


```R
# library
library("Signac")
library("Seurat")
library("hdf5r")
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
##############################################################################



# RNA integration
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



