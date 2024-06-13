# Project

Study KCNC1 mutation (Kcnc1-Arg320His). In human HT KCNC1 missense mutation (dominant negative=mutated non functional prot compete with the functional prot.) associated with epilepsy and cerebellar ataxia (Cerebellum/CB not working well, but histology do not show anything --> scRNAseq may). Generated mice model with same mutation; also got phentoype (HM mice die very young, HT are used here). So KCNC1 is a potassium channel and if it is not working properly, the excitatory neurons will not be slow-down, too much active and lead to epilepsy. KCNC work in tetramer with other KCNC proteins.

A drug has been design to promote potassium channel activation and it rescue the phenotype in mice (currently, tested in human): paper [here](https://pubmed.ncbi.nlm.nih.gov/38266642/)

KCNC1 express day12 in mice neuron (important for fast-spike interneuron (inhibitory); express in Cortex/Ctx and CB. In CB expressed in unipolar brush cell (UBC), purkinje layer interneuron (PLI), Golgi, Molecular Layer Interneuron 2 (MLI2), Purkinje). Histology no marker for MLI2. --> [webtool](https://portal.nemoarchive.org/) for cell type in CB. This [paper](https://www.sciencedirect.com/science/article/pii/S0896627324002484) show that MLI2 where KCNC1 is express, inhibit MLI1 who inhibit Purkinje cell activity. Loss of Purkinje activtiy or Purkinje cells leads to ataxia. 

2 genotype (WT, KCNC1-mutant); 2 brain regions (Ctx, CB); 3 ages (p14,35,180 = 14 is when KCNC1 start being express in intenreunon; p35= 1st phenotype very mild; p180= degeneration)



# Objectives

- Check cell population changes in WT vs mutant; of Ctx and CB. Then DEGs analysis
- check where KCNC1 is express (compare with CB [paper](https://pubmed.ncbi.nlm.nih.gov/34616064/)). Check trajectory of expression too.
- consequence of interneuron not working properly; may be increased excitatory neurons activity (leading to epilepsy); check if more excitatory neurons, or more activtiy (check marker genes of activity maybe? Channel?)
- Check KCNC1 expression over time and within the different cell types. 
    - The weird thing is that the KCNC1 mutation disease is progressive. Not clear why. Maybe where KCNC1 expression pattern change, or the cell types appearition themself change (Like MLI2 appear later?)

--> Focus the analysis on the cell types expressing KCNC1



# Paper for cell type annotation

- CB: [paper](https://pubmed.ncbi.nlm.nih.gov/34616064/); cell types here should be found at p35-p180
  supp file in `docs/41586_2021_3220_MOESM4_ESM.xlsx` with marker genes.
- Ctx: XXX



# Data acess

To access data folow email
```
Service Request : Fileshare_Access_2801 
Requester : goldberge 
Fileshare : goldberg_lab_scb 
AccessType : ReadWrite 
Fileshare Path : \\ressmb05.research.chop.edu\goldberg_lab_scb 
```

Open files app/Connect to server; add Fileshare path and add my credentials

Then, right click on the server and open with terminal

--> Copy all files to `002*/005*`

```bash
cp -r * /scr1/users/roulet/Akizu_Lab/002_scRNAseq/005__Goldberg/input_raw/

```

--> ALL GOOD

File used in `/snRNAseq_Kcnc1_R320H/snRNAseq_Kcnc1_reorganized/*`

For the analysis I will follow the YAP1 scRNAseq in `002003`


## Counting with cellranger count

Within each folder in `/snRNAseq_Kcnc1_R320H/snRNAseq_Kcnc1_reorganized/*` I have two lanes L001 and L002 with I1/I2 and R1/R2 fastq. 








```bash 
conda activate scRNAseq
which cellranger

# Run count using mice genome

## p14 _ WT _ CB and CX
sbatch scripts/cellranger_count_WT_p14_CB_Rep1.sh # 20178099 ok
sbatch scripts/cellranger_count_WT_p14_CB_Rep2.sh # 20178116 ok
sbatch scripts/cellranger_count_WT_p14_CB_Rep3.sh # 20178127 ok
sbatch scripts/cellranger_count_WT_p14_CX_Rep1.sh # 20178154 ok 
sbatch scripts/cellranger_count_WT_p14_CX_Rep2.sh # 20178162 ok
sbatch scripts/cellranger_count_WT_p14_CX_Rep3.sh # 20178170 ok


## p14 _ Kcnc1 _ CB and CX
sbatch scripts/cellranger_count_Kcnc1_p14_CB_Rep1.sh # 20177866 FAIL corrupted; 20178043 FAIL should have delete previous; 20203124 ok
sbatch scripts/cellranger_count_Kcnc1_p14_CB_Rep2.sh # 20177913 FAIL corrupted; 20178032 FAIL should have delete previous; 20203125 ok
sbatch scripts/cellranger_count_Kcnc1_p14_CB_Rep3.sh # 20177921 ok
sbatch scripts/cellranger_count_Kcnc1_p14_CX_Rep1.sh # 20177943 FAIL corrupted; 20178594 FAIL should have delete previous; 20203127 ok
sbatch scripts/cellranger_count_Kcnc1_p14_CX_Rep2.sh # 20177952 ok
sbatch scripts/cellranger_count_Kcnc1_p14_CX_Rep3.sh # 20177959 ok



## p35 _ WT _ CB and CX
sbatch scripts/cellranger_count_WT_p35_CB_Rep1.sh # 20178351 fail fastq path; 20203366 ok
sbatch scripts/cellranger_count_WT_p35_CB_Rep2.sh # 20178354 fail fastq path; 20203367 fail sample name 20250639 ok 
sbatch scripts/cellranger_count_WT_p35_CB_Rep3.sh # 20178357 fail fastq path; 20203368 ok
sbatch scripts/cellranger_count_WT_p35_CX_Rep1.sh # 20178372 ok 
sbatch scripts/cellranger_count_WT_p35_CX_Rep2.sh # 20178379 fail fastq path; 20203383 ok
sbatch scripts/cellranger_count_WT_p35_CX_Rep3.sh # 20178386 fail fastq path; 20203384 ok


## p35 _ Kcnc1 _ CB and CX
sbatch scripts/cellranger_count_Kcnc1_p35_CB_Rep1.sh # 20178221 fail fastq path; 20203342 ok
sbatch scripts/cellranger_count_Kcnc1_p35_CB_Rep2.sh # 20178240 fail fastq path; 20203343 ok
sbatch scripts/cellranger_count_Kcnc1_p35_CB_Rep3.sh # 20178250 fail fastq path; 20203344 ok
sbatch scripts/cellranger_count_Kcnc1_p35_CX_Rep1.sh # 20178258 fail fastq path; 20203351 ok 
sbatch scripts/cellranger_count_Kcnc1_p35_CX_Rep2.sh # 20178280 fail fastq path; 20203352 ok
sbatch scripts/cellranger_count_Kcnc1_p35_CX_Rep3.sh # 20178344 fail fastq path; 20203353 ok



## p180 _ WT _ CB and CX
sbatch scripts/cellranger_count_WT_p180_CB_Rep1.sh # 20178494 fail fastq path; 20203461 ok
sbatch scripts/cellranger_count_WT_p180_CB_Rep2.sh # 20178497 fail fastq path; 20203464 fail sample name 20251147 ok 
sbatch scripts/cellranger_count_WT_p180_CB_Rep3.sh # 20178501 fail fastq path; 20203465 ok
sbatch scripts/cellranger_count_WT_p180_CX_Rep1.sh # 20178518 fail fastq path; 20203505 ok
sbatch scripts/cellranger_count_WT_p180_CX_Rep2.sh # 20178525 fail fastq path; 20203508 ok
sbatch scripts/cellranger_count_WT_p180_CX_Rep3.sh # 20178587 fail fastq path; 20203509 ok


## p180 _ Kcnc1 _ CB and CX
sbatch scripts/cellranger_count_Kcnc1_p180_CB_Rep1.sh # 20178463 fail fastq path; 20203418 ok
sbatch scripts/cellranger_count_Kcnc1_p180_CB_Rep2.sh # 20178466 fail fastq path; 20203419 fail sample name 20250949 ok 
sbatch scripts/cellranger_count_Kcnc1_p180_CB_Rep3.sh # 20178471 fail fastq path; 20203420 ok
sbatch scripts/cellranger_count_Kcnc1_p180_CX_Rep1.sh # 20178478 fail fastq path; 20203432 ok
sbatch scripts/cellranger_count_Kcnc1_p180_CX_Rep2.sh # 20178481 fail fastq path; 20203437 ok
sbatch scripts/cellranger_count_Kcnc1_p180_CX_Rep3.sh # 20178485 fail fastq path; 20203438 ok


```


--> Bug for:
- `cellranger_count_Kcnc1_p14_CB_Rep1`: `Sequence and quality length mismatch: file: "/scr1/users/roulet/Akizu_Lab/002_scRNAseq/003__YAP1/input/24hgastruloidhumanUN/24hgastruloidhumanUN_S1_L001_R2_001.fastq.gz", line: 292599520`
    - Another file available in `8_31_23 core reanalyzed corrupted files`

-----> To avoid error, let's directly use the reanalyzed one for the concerned samples: 
- Kcnc1_p14_CB_Rep1
- Kcnc1_p14_CB_Rep2
- Kcnc1_p14_CX_Rep1
- WT_p35_CX_Rep1


--> All samples successfully counted








## RNA contamination and doublet detection
- doublet detection using [scrublet](https://github.com/swolock/scrublet) **on the filtered matrix**
- ambient RNA correction using `soupX` in R before generating the Seurat object

Too many samples so let's run scrublet in 6 separate bash jobs

```bash
conda deactivate # base environment needed

python3 scrublet.py [input_path] [output_path]

# Run doublet detection/scrublet genotype/time tissue (CB/CX) together

sbatch scripts/scrublet_WT_p14.sh # 20255539 fail miss cp script; 20259650 ok
sbatch scripts/scrublet_Kcnc1_p14.sh # 20255755 fail miss cp script; 20259651 ok
sbatch --dependency=afterany:20250639 scripts/scrublet_WT_p35.sh # 20255864 ok
sbatch scripts/scrublet_Kcnc1_p35.sh # 20255900 fail miss cp script; 20259652 ok
sbatch --dependency=afterany:20251147 scripts/scrublet_WT_p180.sh # 20255934 ok
sbatch --dependency=afterany:20250949 scripts/scrublet_Kcnc1_p180.sh # 20255974 Kcnc1_p180_CB_Rep3 fail in predicting but it output file; RAS ok

# Run manual treshold doublet identification for Kcnc1_p180_CB_Rep3 sample
## Run scripts that include histogram to define treshold value
python3 scripts/scrublet_doublets_histogram.py Kcnc1_p180_CB_Rep3/outs/filtered_feature_bc_matrix output/doublets/Kcnc1_p180_CB_Rep3.tsv
## Run scripts with manual treshold, test 2 just to see
python3 scripts/scrublet_doublets_manualTresh.py Kcnc1_p180_CB_Rep3/outs/filtered_feature_bc_matrix output/doublets/Kcnc1_p180_CB_Rep3_tresh02.tsv 0.2
python3 scripts/scrublet_doublets_manualTresh.py Kcnc1_p180_CB_Rep3/outs/filtered_feature_bc_matrix output/doublets/Kcnc1_p180_CB_Rep3_tresh025.tsv 0.25
python3 scripts/scrublet_doublets_manualTresh.py Kcnc1_p180_CB_Rep3/outs/filtered_feature_bc_matrix output/doublets/Kcnc1_p180_CB_Rep3_tresh04.tsv 0.4


```
Doublet detection score (YAP1scRNAseq ranged between 0.1 to 42%): now range from 0 to 30% (median ~10%)
--> Table in `QC/QC_metrics_all.xlsx`

- IMPORTANT NOTE: sample `Kcnc1_p180_CB_Rep3` "failed to automatically identify doublet score threshold". From [this](https://www.biostars.org/p/9572647/), let's assign trehsold to XXX + a new script as been created to handle manual treshold setting `scripts/scrublet_doublets_manualTresh.py` + generate a script that generate histogram `scripts/scrublet_doublets_histogram.py`
-----> Treshold of 0.2, 0.25 and 0.4 have been tested


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
sc = load10X('WT_p14_CB_Rep1/outs')   # CHANGE FILE NAME HERE

## Assess % of conta
pdf("output/soupX/autoEstCont_WT_p14_CB_Rep1.pdf", width=10, height=10)   # CHANGE FILE NAME HEREs
sc = autoEstCont(sc) 
dev.off()
## Generate the corrected matrix
out = adjustCounts(sc)
## Save the matrix
save(out, file = "output/soupX/WT_p14_CB_Rep1.RData") # CHANGE FILE NAME HERE

```

--> sample `Kcnc1_p35_CB_Rep2` is highly contaminated, add to add `forceAccept = TRUE` in `sc = autoEstCont(sc, forceAccept = TRUE)`; all other sampels are <30% (median 10%)


# Data integration 

**Step1_General Trend**: 
Investigate *genotype effect for each brain regions*
- Integrate all time points and genotype together within the same UMAP:
  - **1 UMAP per brain regions** --> highlight genotype/time

**Step2_Time specific**:
Investigate *genotype effect at each time point for each brain regions*
- Integrate genotype together for each time points
  - **1 UMAP per time and per brain regions** --> highlight genotype


How to deal with the Bio Rep for data integration is discuss [here](https://github.com/satijalab/seurat/issues/4753). Author mentioned to "run integration for all samples together, but set the control replicates as reference". Not clear... I read [here](https://github.com/satijalab/seurat/issues/2059) that integrating replicate together first, then integrating the condition works best:
- Step1: Integrate the 3 Bio rep for WT and for Mt separately for each time point > integrate the 3 time points together.
- Step2: Inegrate the 3 Bio Rep for WT and for Mt separately and integrate at each time point.

--> Follow [this](https://github.com/satijalab/seurat/issues/6003) or this [Workflow2](https://github.com/satijalab/seurat/discussions/5702):
- normalize all samples/replicate individually (eg. WT_p14_CB_R1 + WT_p14_CB_R1 + WT_p14_CB_R3 = WT_p14_CB and same for Kcnc1 at each time points)
- integrate then the other genotype (WT_p14_CB + Kcnc1_p14_CB = WT-Kcnc1_p14_CB) = **Step1_General Trend=2 step data integration**
- integrate then the time points (WT-Kcnc1_p14_CB + WT-Kcnc1_p35_CB + WT-Kcnc1_p180_CB = WT-Kcnc1_CB) = **Step2_Time specific=3 step data integration**


## General trend

Then `conda activate scRNAseqV2` and go into R and filtered out **RNA contamination and start with SEURAT**.

Tuto seurat [here](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html)


```bash
conda activate scRNAseqV2
```

For the sample `Kcnc1_p180_CB_Rep3`; I pick the scrublet Treshold showing the better corr. score for RNAfeature-doublet in the `FeatureScatter_*` plot. They all have the same score... So I pick 0.25 which look good regarding the histogram.

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
  WT_p14_CB_Rep1 = "output/soupX/WT_p14_CB_Rep1.RData",
  WT_p14_CB_Rep2 = "output/soupX/WT_p14_CB_Rep2.RData",
  WT_p14_CB_Rep3 = "output/soupX/WT_p14_CB_Rep3.RData",
  Kcnc1_p14_CB_Rep1 = "output/soupX/Kcnc1_p14_CB_Rep1.RData",
  Kcnc1_p14_CB_Rep2 = "output/soupX/Kcnc1_p14_CB_Rep2.RData",
  Kcnc1_p14_CB_Rep3 = "output/soupX/Kcnc1_p14_CB_Rep3.RData",

  WT_p35_CB_Rep1 = "output/soupX/WT_p35_CB_Rep1.RData",
  WT_p35_CB_Rep2 = "output/soupX/WT_p35_CB_Rep2.RData",
  WT_p35_CB_Rep3 = "output/soupX/WT_p35_CB_Rep3.RData",
  Kcnc1_p35_CB_Rep1 = "output/soupX/Kcnc1_p35_CB_Rep1.RData",
  Kcnc1_p35_CB_Rep2 = "output/soupX/Kcnc1_p35_CB_Rep2.RData",
  Kcnc1_p35_CB_Rep3 = "output/soupX/Kcnc1_p35_CB_Rep3.RData",

  WT_p180_CB_Rep1 = "output/soupX/WT_p180_CB_Rep1.RData",
  WT_p180_CB_Rep2 = "output/soupX/WT_p180_CB_Rep2.RData",
  WT_p180_CB_Rep3 = "output/soupX/WT_p180_CB_Rep3.RData",
  Kcnc1_p180_CB_Rep1 = "output/soupX/Kcnc1_p180_CB_Rep1.RData",
  Kcnc1_p180_CB_Rep2 = "output/soupX/Kcnc1_p180_CB_Rep2.RData",
  Kcnc1_p180_CB_Rep3 = "output/soupX/Kcnc1_p180_CB_Rep3.RData"
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
  if (sample_name != "Kcnc1_p180_CB_Rep3") {
    seurat_objects[[sample_name]] <- add_doublet_information(sample_name, seurat_objects[[sample_name]])
  }
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list

### Kcnc1_p180_CB_Rep3 Done separately as doublet treshold set manually
doublets <- read.table("output/doublets/Kcnc1_p180_CB_Rep3_tresh025.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
Kcnc1_p180_CB_Rep3 <- AddMetaData(Kcnc1_p180_CB_Rep3,doublets)
Kcnc1_p180_CB_Rep3$Doublet_score <- as.numeric(Kcnc1_p180_CB_Rep3$Doublet_score) # make score as numeric
## Re include  Kcnc1_p180_CB_Rep3 in our list

seurat_objects[["Kcnc1_p180_CB_Rep3"]] <- Kcnc1_p180_CB_Rep3


##################### Doublet testing for sample Kcnc1_p180_CB_Rep3 with manual treshold for doublet detection #########
Kcnc1_p180_CB_Rep3[["percent.mt"]] <- PercentageFeatureSet(Kcnc1_p180_CB_Rep3, pattern = "^mt-")
Kcnc1_p180_CB_Rep3[["percent.rb"]] <- PercentageFeatureSet(Kcnc1_p180_CB_Rep3, pattern = "^Rp[sl]")

doublets <- read.table("output/doublets/Kcnc1_p180_CB_Rep3_tresh04.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
Kcnc1_p180_CB_Rep3 <- AddMetaData(Kcnc1_p180_CB_Rep3,doublets)
Kcnc1_p180_CB_Rep3$Doublet_score <- as.numeric(Kcnc1_p180_CB_Rep3$Doublet_score) # make score as numeric
head(Kcnc1_p180_CB_Rep3[[]])

pdf("output/seurat/FeatureScatter_QC_mt_doublet_Kcnc1_p180_CB_Rep3_tresh04.pdf", width=5, height=5)
FeatureScatter(Kcnc1_p180_CB_Rep3, feature1 = "percent.mt", feature2 = "Doublet_score")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_Kcnc1_p180_CB_Rep3_tresh04.pdf", width=5, height=5)
FeatureScatter(Kcnc1_p180_CB_Rep3, feature1 = "nFeature_RNA", feature2 = "Doublet_score")
dev.off()
########################################################################################################################


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


# QC filtering _ V1

### V1 not super stringeant; mit > 20 and RNAfeaure 100
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 100 & seurat_object@meta.data$QC == 'Pass', 'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 100 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 20 & seurat_object@meta.data$QC == 'Pass', 'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 100 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
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
write.table(qc_summary_combined, file = "output/seurat/QC_summary_V1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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
# Combine all summaries into one data frame
phase_summary_combined <- do.call(rbind, phase_summary_list)
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Cell type annotation

############ SAVE workspace #################################################################
# save.image(file = "CB_QC_V1.RData")
load("CB_QC_V1.RData")
#############################################################################################



## Work on the cleanest WT sample for each time point; WT_p14_CB_Rep2, WT_p35_CB_Rep3, WT_p180_CB_Rep3
# elbow
WT_p14_CB_Rep2 <- SCTransform(WT_p14_CB_Rep2, method = "glmGamPoi", ncells = 13576, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000)
WT_p14_CB_Rep2 <- RunPCA(WT_p14_CB_Rep2, npcs = 50, verbose = FALSE)

## ELBOW ###########################################################################
pdf("output/seurat/Elbow_WT_p14_CB_Rep2.pdf", width=10, height=10)
ElbowPlot(WT_p14_CB_Rep2) # 6 or 10
dev.off()
########################################################################### USELESS...


WT_p14_CB_Rep2 <- RunPCA(WT_p14_CB_Rep2, npcs = 10, verbose = FALSE)
WT_p14_CB_Rep2 <- RunPCA(WT_p14_CB_Rep2, npcs = 10, verbose = FALSE)
WT_p14_CB_Rep2 <- RunUMAP(WT_p14_CB_Rep2, reduction = "pca", dims = 1:10, verbose = FALSE)
WT_p14_CB_Rep2 <- FindNeighbors(WT_p14_CB_Rep2, reduction = "pca", k.param = 15, dims = 1:10)
WT_p14_CB_Rep2 <- FindClusters(WT_p14_CB_Rep2, resolution = 0.4, verbose = FALSE, algorithm = 4)


pdf("output/seurat/UMAP_WT_p14_CB_Rep2-dim10kparam15res04.pdf", width=5, height=5)
DimPlot(WT_p14_CB_Rep2, reduction = "umap", label=TRUE)
dev.off()


# Check QC metrics

pdf("output/seurat/VlnPlot_QCmetrics_WT_p14_CB_Rep2.pdf", width=20, height=5)
VlnPlot(WT_p14_CB_Rep2,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()



# Check some genes

DefaultAssay(WT_p14_CB_Rep2) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kcnc1.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Kcnc1"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_Purkinje.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Ppp1r17"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_Granular.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Gabra6"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_Golgi.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Slc6a5", "Grm2", "Sst"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_MLI1.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Prkcd", "Sorcs3", "Ptprk"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_MLI2.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Prkcd", "Nxph1", "Cdh22"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_PLI_3.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Htr2a", "Edil3"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_PLI_1PLI_2.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Aldh1a3", "Slc6a5"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_UnipolarBrush.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Eomes"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-Kozareva2021_Bergmann.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB_Rep2, features = c("Gdf10"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()



























# Integrate the Bio Rep

# clustering
## Optimal parameter 50 dim
WT_p14_CB_Rep1 <- SCTransform(WT_p14_CB_Rep1, method = "glmGamPoi", ncells = 12877, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% 
    RunPCA(npcs = 50, verbose = FALSE)
WT_p14_CB_Rep2 <- SCTransform(WT_p14_CB_Rep2, method = "glmGamPoi", ncells = 13576, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 50, verbose = FALSE)
WT_p14_CB_Rep3 <- SCTransform(WT_p14_CB_Rep3, method = "glmGamPoi", ncells = 13420, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>%
    RunPCA(npcs = 50, verbose = FALSE)

# Data integration (check active assay is 'SCT')
srat.list <- list(WT_p14_CB_Rep1 = WT_p14_CB_Rep1, WT_p14_CB_Rep2 = WT_p14_CB_Rep2, WT_p14_CB_Rep3 = WT_p14_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_p14_CB.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_p14_CB.sct <- IntegrateData(anchorset = WT_p14_CB.anchors, normalization.method = "SCT")
set.seed(42)
DefaultAssay(WT_p14_CB.sct) <- "integrated"
WT_p14_CB.sct <- RunPCA(WT_p14_CB.sct, verbose = FALSE, npcs = 50)
WT_p14_CB.sct <- RunUMAP(WT_p14_CB.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_p14_CB.sct <- FindNeighbors(WT_p14_CB.sct, reduction = "pca", k.param = 15, dims = 1:50)
WT_p14_CB.sct <- FindClusters(WT_p14_CB.sct, resolution = 0.4, verbose = FALSE, algorithm = 4)


pdf("output/seurat/UMAP_WT_p14_CB_splitOrigIdent.pdf", width=12, height=4)
DimPlot(WT_p14_CB.sct, reduction = "umap", split.by = "orig.ident", label=TRUE)
dev.off()
pdf("output/seurat/UMAP_WT_p14_CB.pdf", width=12, height=4)
DimPlot(WT_p14_CB.sct, reduction = "umap", label=TRUE)
dev.off()
pdf("output/seurat/UMAP_WT_p14_CB_groupOrigIdent.pdf", width=12, height=4)
DimPlot(WT_p14_CB.sct, reduction = "umap", group.by = "orig.ident" , label=TRUE)
dev.off()

# Check QC metrics

pdf("output/seurat/VlnPlot_QCmetrics_WT_p14_CB.pdf", width=25, height=5)
VlnPlot(srat_WT,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()









# Check some genes

DefaultAssay(WT_p14_CB.sct) <- "SCT" # For vizualization either use SCT or norm RNA
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kcnc1.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Kcnc1"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_Purkinje.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Ppp1r17"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_Granular.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Gabra6"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_Golgi.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Slc6a5", "Grm2", "Sst"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_MLI1.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Prkcd", "Sorcs3", "Ptprk"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_MLI2.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Prkcd", "Nxph1", "Cdh22"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_PLI_3.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Htr2a", "Edil3"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_PLI_1PLI_2.pdf", width=7, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Aldh1a3", "Slc6a5"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_UnipolarBrush.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Eomes"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-Kozareva2021_Bergmann.pdf", width=5, height=5)
FeaturePlot(WT_p14_CB.sct, features = c("Gdf10"), max.cutoff = 5, cols = c("grey", "red"))
dev.off()














































```

