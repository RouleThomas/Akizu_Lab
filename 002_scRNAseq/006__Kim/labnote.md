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


## TEST RNA regression  --> THE WINNER PRO WINNER !!! RNA regression 
RNA_WT <- SCTransform(RNA_WT, method = "glmGamPoi", ncells = 6637, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 50, verbose = FALSE)
RNA_Bap1KO <- SCTransform(RNA_Bap1KO, method = "glmGamPoi", ncells = 6938, vars.to.regress = c("percent.mt","nCount_RNA","percent.rb"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 50, verbose = FALSE)
# Data integration (check active assay is 'SCT')
srat.list <- list(RNA_WT = RNA_WT, RNA_Bap1KO = RNA_Bap1KO)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

embryo.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
RNA_WT_Bap1KO.sct <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT")

set.seed(42)

DefaultAssay(RNA_WT_Bap1KO.sct) <- "integrated"

RNA_WT_Bap1KO.sct <- RunPCA(RNA_WT_Bap1KO.sct, verbose = FALSE, npcs = 50)
RNA_WT_Bap1KO.sct <- RunUMAP(RNA_WT_Bap1KO.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
RNA_WT_Bap1KO.sct <- FindNeighbors(RNA_WT_Bap1KO.sct, reduction = "pca", k.param = 30, dims = 1:50)
RNA_WT_Bap1KO.sct <- FindClusters(RNA_WT_Bap1KO.sct, resolution = 0.7, verbose = FALSE, algorithm = 4)

RNA_WT_Bap1KO.sct$orig.ident <- factor(RNA_WT_Bap1KO.sct$orig.ident, levels = c("RNA_WT", "RNA_Bap1KO")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Bap1KO_noSplit-QCV2_dim50kparam30res07algo4_noCellCycleRegression.pdf", width=8, height=5)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", label=TRUE)
dev.off()



DefaultAssay(RNA_WT_Bap1KO.sct) <- "SCT" # For vizualization either use SCT or norm RNA


pdf("output/seurat/FeaturePlot_SCT_RNA_WT_Bap1KO-allMarkersList4-QCV2_dim50kparam30res07algo4_noCellCycleRegression
.pdf", width=15, height=10)
FeaturePlot(RNA_WT_Bap1KO.sct, features = c("Pax6", "Eomes", "Prox1", "Neurod1", "Cck", "Crym", "Snca", "Tac2", "Pantr1", "Satb2", "Gad1", "Lhx1", "Nts"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()






pdf("output/seurat/UMAP_WT_Bap1KO_split_V1.pdf", width=8, height=6)
DimPlot(RNA_WT_Bap1KO.sct, reduction = "umap", split.by = "orig.ident", label=TRUE)
dev.off()









```