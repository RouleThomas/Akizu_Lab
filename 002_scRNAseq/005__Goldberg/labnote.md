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
  supp file in `docs/41586_2021_3220_MOESM4_ESM.xlsx` with marker genes + Fig1C


- Ctx: XXX

--> Cool [db](https://www.panglaodb.se/markers.html) for cell type marker genes!  

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


## General trend Cerebellum (CB)

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


##### QC filtering _ V1 ############################################

### V1 not super stringeant; mit > 20 and RNAfeaure 100 #############################################
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
#####################################################################################################



##### QC filtering _ V2 (used for p35) ############################################

### mit > 15, rb > 22 and RNAfeaure 100 #############################################
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 100 & seurat_object@meta.data$QC == 'Pass', 'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 100 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 15 & seurat_object@meta.data$QC == 'Pass', 'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 15 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 22 & seurat_object@meta.data$QC == 'Pass', 'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 22 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list
#####################################################################################################



##### QC filtering _ V3 (used for p14) ############################################

###  mit > 15, rb > 10 and RNAfeaure 400 #############################################
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC == 'Pass', 'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 15 & seurat_object@meta.data$QC == 'Pass', 'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 15 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 10 & seurat_object@meta.data$QC == 'Pass', 'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 10 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list
#####################################################################################################



##### QC filtering _ V4 (used for p180) ############################################

###  mit > 7, rb > 5 and RNAfeaure 400 #############################################
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC == 'Pass', 'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 7 & seurat_object@meta.data$QC == 'Pass', 'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 7 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 5 & seurat_object@meta.data$QC == 'Pass', 'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 5 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list
#####################################################################################################


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
write.table(qc_summary_combined, file = "output/seurat/QC_summary_V2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # for p35?
write.table(qc_summary_combined, file = "output/seurat/QC_summary_V3.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # for p14 and p35?
write.table(qc_summary_combined, file = "output/seurat/QC_summary_V4.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # for p180

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
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # for p35?
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V3.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # for p14 and p35?
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V4.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # for p180




############ SAVE samples #################################################################
save_seurat_objects <- function(seurat_objects_list, output_dir) {
  for (sample_name in names(seurat_objects_list)) {
    file_name <- paste0(output_dir, sample_name, "_V2_ReProcess_numeric.rds") # _V2_numeric \\\--> V3 for p14 ; v4 for p180; V2 for p35
    saveRDS(seurat_objects_list[[sample_name]], file = file_name) # V2 got a ReProcess for sample p35
  }
}
output_dir <- "output/seurat/"
## Call the function to save the Seurat objects
save_seurat_objects(seurat_objects, output_dir)
###############################################################################################
############# READ samples  (QC V1 V2 or V3 = V1_numeric...)################################################
# Function to load Seurat objects
load_seurat_objects <- function(file_paths) {
  seurat_objects <- list()
  for (file_path in file_paths) {
    sample_name <- gsub("_V2_ReProcess_numeric.rds", "", basename(file_path))
    seurat_objects[[sample_name]] <- readRDS(file_path)
  }
  return(seurat_objects)
}
output_dir <- "output/seurat/"
file_paths <- list.files(output_dir, pattern = "_V2_ReProcess_numeric.rds$", full.names = TRUE)
# Call the function to load the Seurat objects
seurat_objects <- load_seurat_objects(file_paths)
# Loop through the list and assign each Seurat object to a variable with the same name
for (sample_name in names(seurat_objects)) {
  assign(sample_name, seurat_objects[[sample_name]])
}
# 1 sample: WT_p14_CB_Rep2 <- readRDS(file = "output/seurat/WT_p14_CB_Rep2_V1_numeric.rds")
## Kcnc1_p14_CB_Rep1 <- readRDS(file = "output/seurat/Kcnc1_p14_CB_Rep1_V1_numeric.rds")
################################################################################################







## Work on the cleanest WT sample for each time point; WT_p14_CB_Rep2, WT_p35_CB_Rep3, WT_p180_CB_Rep3

################### WT_p14_CB_Rep2
set.seed(42)
WT_p14_CB_Rep2 <- SCTransform(WT_p14_CB_Rep2, method = "glmGamPoi", ncells = 13576, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000)
WT_p14_CB_Rep2 <- RunPCA(WT_p14_CB_Rep2, npcs = 30, verbose = FALSE)

## ELBOW ###########################################################################
pdf("output/seurat/Elbow_WT_p14_CB_Rep2.pdf", width=10, height=10)
ElbowPlot(WT_p14_CB_Rep2) # 6 or 10
dev.off()
########################################################################### USELESS...



WT_p14_CB_Rep2 <- RunPCA(WT_p14_CB_Rep2, npcs = 30, verbose = FALSE)
WT_p14_CB_Rep2 <- RunUMAP(WT_p14_CB_Rep2, reduction = "pca", dims = 1:30, verbose = FALSE)
WT_p14_CB_Rep2 <- FindNeighbors(WT_p14_CB_Rep2, reduction = "pca", k.param = 40, dims = 1:30)
WT_p14_CB_Rep2 <- FindClusters(WT_p14_CB_Rep2, resolution = 0.7, verbose = FALSE, algorithm = 4)


pdf("output/seurat/UMAP_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression.pdf", width=5, height=5)
DimPlot(WT_p14_CB_Rep2, reduction = "umap", label=TRUE)
dev.off()





# Check some genes

DefaultAssay(WT_p14_CB_Rep2) <- "SCT" # For vizualization either use SCT or norm RNA


pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-List1.pdf", width=25, height=25)
FeaturePlot(WT_p14_CB_Rep2, features = c("Kcnc1","Ppp1r17", "Gabra6", "Slc6a5", "Grm2", "Sst", "Prkcd", "Sorcs3", "Ptprk", "Nxph1", "Cdh22", "Htr2a", "Edil3", "Aldh1a3", "Slc6a5","Eomes", "Gdf10"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-purkinje1.pdf", width=25, height=25)
FeaturePlot(WT_p14_CB_Rep2, features = c("Calb1", "Gad2", "Grid2", "Gad1", "Slc1a6", "Pcp4", "Car8", "Camk2a", "Hcn1", "Gria3","Slc32a1", "Slc12a5"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()




pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-bergmann1.pdf", width=25, height=25)
FeaturePlot(WT_p14_CB_Rep2, features = c("Fabp7", "Ntsr2", "Egfr", "Ttc21b", "Hepacam", "Zeb2", "Ufl1", "Nrxn3", "Atp7a", "Bicd2", "Gja1", "Metrn"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-KimList4.pdf", width=25, height=25)
FeaturePlot(WT_p14_CB_Rep2, features = c("Pax6", "Eomes", "Prox1", "Neurod1", "Cck", "Crym", "Snca", "Tac2", "Pantr1", "Satb2", "Gad1", "Lhx1", "Nts"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-KimList4extended.pdf", width=25, height=25)
FeaturePlot(WT_p14_CB_Rep2, features = c("Sema5a","Grin2d", "Reln", "Calb1", "Npy", "Gria3", "Lhx6","Insm1","Crym","Nrp2","Hs3st1", "Nrn1","Igfbpl1", "Frmd4b","Itpr1","Lhx1","Nr4a2", "Lmo3", "B3gat1","Csf1r", "Gpr34", "Gpr183", "Cx3cr1","Pdgfra", "Olig1"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()


# PanglaoDB all brains top 12 genes
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoAdrenergicNeurons.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Pnmt", "Npff", "Ddc", "Dbh", "Slc18a2", "Slc12a7", "Syt1", "Th"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoAstrocyte.pdf", width=25, height=25)
FeaturePlot(WT_p14_CB_Rep2, features = c("Gfap", "Slc1a2", "Acsl6", "Agt", "Aqp4", "Apoe", "S100B", "Sox9", "Gsta4", "Srr", "Aldh1l1","Slc39a12"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoAstrocyte.pdf", width=25, height=25)
FeaturePlot(WT_p14_CB_Rep2, features = c("Gfap", "Slc1a2", "Acsl6", "Agt", "Aqp4", "Apoe", "S100B", "Sox9", "Gsta4", "Srr", "Aldh1l1","Slc39a12"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoCR.pdf", width=25, height=25) # CR only in Ctx, not CB
FeaturePlot(WT_p14_CB_Rep2, features = c("Car10", "Nrg1", "Edil3", "Flrt3", "Cdh4", "Srgap3", "Fat3", "Lingo2", "Lhx1", "Pnoc", "Cntnap2", "Zcchc12"), max.cutoff = 2, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoCholinNeurons.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Chat", "Slc5a7", "Slc18a3", "Ache", "Tac1", "Acly", "Brca1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoChoroidPlexCells.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Ttr", "Kl", "Clic6", "Prlr", "Chmp1a", "Slc26a11", "Slc23a2", "Wfikkn2", "Slc2a12", "Cldn1", "Slc29a4", "Slc13a4"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoDopaNeur.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Slc6a3", "Th", "Pitx3", "Smad3", "Neurod6", "Slc18a2", "Ddc", "Zim3", "Tenm1", "Scn2a", "Prkcg", "Mapk8ip2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoEpendymal.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Rabl2", "Cfap54", "Ccdc153", "Foxj1", "Pifo", "Dynlrb2", "Rsph1", "Cfap44", "Pcp4l1", "Ak8", "Tmem212", "Tm4sf1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoGABAneur.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Gad1", "Gad2", "Slc6a1", "Gabbr2", "Gadd45b", "Pax2", "Slc32a1", "Vip", "Pvalb", "Dlx1", "Tnfaip8l3", "Sema3c"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoGlutaNeur.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Slc17a7", "Slc17a6", "Slc17a8", "Slc1a1", "Gls", "Meis2", "Slc1a6", "Grin1", "Grin2b", "Slc1a2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoGlycNeur.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Slc6a9", "Slc32a1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoImmatureNeur.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Neurog2", "Tbr1", "Epha3", "Yy2", "Icam5", "Foxd2", "Slc6a5", "Bcl2", "Nf1", "Syn1", "Neurod1", "Creb1", "Notch1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoInterneuron.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Pvalb", "Cck", "Npy", "Nos1", "Sst", "Vip", "Fgf12", "Nxph1", "Eomes", "Kcnc2", "Calb2", "Calb1", "Vsx2", "En1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoMeningCell.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Igfbp2", "Slc47a1", "Nov", "Nnat", "Ptgds", "Il33", "Pdgfra", "Lum", "Dcn", "Foxc1", "Vtn", "Igf2", "Alx4"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoMicroglia.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Itgam", "Tmem119", "Cxcr1", "P2ry12", "Aif1", "Csf1r", "Aldh1a2", "Adrb2", "Sall1", "Sphk1", "Colec12", "Ccrl2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoMotorNeur.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Chat", "Mnx1", "Isl2", "Vsx2", "En1", "Evx1", "Evx2", "Fgf1", "Lhx3", "Nkx6-1", "Ngfr", "Reg2", "Sim1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoNSC.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Pbk", "Sox1", "Sp9", "Arx", "Hes5", "Eomes", "Ascl2", "Dynlt1c", "Barhl2", "Lhx9", "Shox2", "Dbx1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoNeuroblast.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Ntng1", "Dcx", "Dlx2", "Neurog1", "Ncam1", "Neurod1", "Dab1", "Cnr1", "Cckbr", "Grm5", "Ezh2", "Pros1", "Mark2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoNeuroendocrineCells.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Nkx2-2", "Npy", "Gpr88", "Ccdc85c", "Tac1", "Penk", "Dpf1", "Prss12", "Serpini1", "Cck", "Arhgig", "Dusp26", "Negr1", "Spock1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoNeurons.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Prph", "Disp2", "Tubb3", "Rbfox3", "Chat", "Isl1", "Nefm", "Nrgn", "Sgip1", "Eno2", "Th", "Calb1", "Epo", "Csf3"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoNoradr.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Th", "Slc18a2", "Ddc", "Dbh", "Slc6a2", "Slc39a11", "Slc9b2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoOPC.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Pdgfra", "Aldoc", "Olig1", "Epn2", "Neu4", "Nkx6-2", "Fyn", "Tnr", "Pcdh15", "Cspg5", "Nnat", "Etv5", "Slc25a47"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoOligodendrocyte.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Mog", "Mbp", "Mag", "Cldn14", "Klk6", "Eml1", "Nipal4", "Plp1", "Ermn", "Plekhh1", "Il33", "Cntf", "Pde8a", "Myo1d", "Itgb4"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoPinealocyte.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Sag", "Gngt2", "Rom1", "Adra1b", "Pmepa1", "Pde10a", "Crem", "Aanat", "Gng13", "Cacna1f", "Pde6b", "Neurod1", "Chrnb4", "Chrna3"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoPurkinjeNeurons.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Calb1", "Gad2", "Grid2", "Gad1", "Slc1a6", "Pcp4", "Car8", "Camk2a", "Hcn1", "Gria3","Slc32a1", "Slc12a5","Cnpy1", "Dlg4", "Slc24a3", "Clstn3", "Gabra1", "Slc24a2", "Kcnma1", "Penk", "Hivep2", "Nrip3", "Tshz1", "Pcbp3", "Gpsm1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoPyramidalCells.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Nrgn", "Dab1", "Dlg4", "Ccm2", "Map2", "Calb1", "Adrb2", "Rtn4", "Pde2a", "Kcnb1","Kcnb2", "Kcnq2", "Kcnq3", "Satb2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoRadialGliaCells.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Pax6", "Slc1a3", "Pdgfd", "Gli3", "Notch3", "Vcam1", "Hes5", "Olig2", "Gfap", "Emx2", "Cdh4", "Spry1", "Axin2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoRetinalGanglionCells.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Pou4f2", "Rbpms", "Pou4f1", "Isl1", "Grip1", "Cpne4", "Kctd8", "Mab21l2", "Dtx1", "Barhl2", "Narf", "Tox2", "Zfp580"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoSatelliteGliaCells.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Ptgds", "Il17a", "Sox2", "Glul", "Slc1a3", "P2rx7", "P2ry12", "Tlr4", "Cnga3", "Aqp4", "Kcnj10", "Cxcl8", "Cdh1", "Ednrb", "Gfap", "Slc6a1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoSchwannCells.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Mpz", "Sox10", "S100b", "Nov", "Cryab", "Matn2", "Aldoc", "Nf2", "Plp1", "Pmp22", "Cnp", "Nrn1", "Cadm4", "Kcna1", "Rxrg", "Apod"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoSerotonergicNeurons.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Fev", "Slc6a4", "Tph2", "Ddc", "Slc18a2", "Esm1", "Slc22a3", "Tph1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoTanycytes.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Dio2", "Scn7a", "Rax", "Adm", "Crym", "Gpr50", "Vcan", "Tgfb2", "Ctgf", "Rgs7bp", "Rgcc", "Igfbp5", "Fndc3"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-panglaoTrigeminalneurons.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Fgf13", "Crebrf", "Kcnab2", "Atl1", "Hspb8", "Fkbp1b", "Tln2", "Synm", "Gdap1", "Dgkz", "Prune2", "Slc17a6", "Phf24", "Dgkh"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-ChatGPTcluster3.pdf", width=15, height=15)
FeaturePlot(WT_p14_CB_Rep2, features = c("Sox2", "Nestin", "Dcx", "Neurod1", "Tubb3" , "Ascl1","P2ry12", "Tmem119", "Cx3cr1","Glast" , "Aldh1l1","Pdgfrb", "Anpep" ,"Cd45", "Emr1", "Cd11b"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-dim30kparam40res07-countMtRbRegression-List2.pdf", width=30, height=70)
FeaturePlot(WT_p14_CB_Rep2, features = c("Calb1", "Slc1a6", "Car8","Gabra6","Sorcs3", "Ptprk","Nxph1", "Cdh22","Edil3","Fabp7", "Zeb2", "Hepacam","Pax6","Gad1", "Gria3","Itpr1","Slc1a2", "Apoe", "Aqp4", "Slc39a12","Glul", "Slc1a3", "Ednrb","Ttr", "Clic6", "Slc13a4", "Kl","Cfap54", "Ccdc153", "Cfap44", "Tmem212","Gad2", "Slc6a1","Ptgds", "Dcn", "Colec12","Ntng1", "Grm5","Slc18a2", "Ddc","Aldoc", "Nnat" ,"Mbp", "Mag", "Plp1", "Pdgfd", "Gli3","Rbpms","S100b", "Cnp","Slc6a4", "Tph2" ), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


# Selected marker genes:
Purkinje= Calb1, Slc1a6, Car8
Granular= Gabra6
Golgi=
MLI1= Sorcs3, "Ptprk"
MLI2=  "Nxph1", "Cdh22",
PLI_3= Edil3
PLI1_2=
Unipolar_brush=
Bergmann= Fabp7, Zeb2, Hepacam

From 006__Kim:
NSC= Pax6
IN= Gad1, Gria3
PyNs_RSC_UL (Retrosplenial Cortical Pyramidal neurons, upper layer)= Itpr1

From Pangeo:
Astrocyte= Slc1a2, Apoe, Aqp4, Slc39a12
Astrocyte(but identfy wit hsatellite glia cell marker)= Glul, Slc1a3, Ednrb
Choroid plexus cells= Ttr, Clic6, Slc13a4, Kl # very small
Ependymal cells= Cfap54, Ccdc153, Cfap44, Tmem212 # very small
GABAergic neurons (to put with IN from Kim)= Gad2, Slc6a1
Interneurons= Nxph1
Meningeal cells= Ptgds, Dcn, 
Microglia= Colec12
Neuroblast= Ntng1, Grm5
Noradrenergic neurons= Slc18a2, Ddc # very small
OPC= Aldoc, Nnat #
Oligodendrocyte= Mbp, Mag, Plp1
Radial Glia Cells= Slc1a3, Pdgfd, Gli3
RetinalGanglionCells (not thiscell type but very specific)= Rbpms
Schwann Cells= S100b, Aldoc, Cnp
Serotonergic neurons= Slc6a4, Tph2, Ddc, Slc18a2 # very small



# Check QC metrics

pdf("output/seurat/VlnPlot_QCmetrics_WT_p14_CB_Rep2.pdf", width=20, height=5)
VlnPlot(WT_p14_CB_Rep2,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()




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

# --> OPTIMAL PARAMETERS NOT FOUND YET for WT_p14_CB_Rep2 !!!

## Automatic cell type annotation

### Find all markers 
all_markers <- FindAllMarkers(WT_p14_CB_Rep2, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(all_markers, file = "output/seurat/srat_WT_p14_CB_Rep2_dim30kparam40res07_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)


############################ EasyCellType automatic annotation ##########################################
# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_WT_p14_CB_Rep2_dim30kparam40res07_all_markers.txt", header = TRUE, row.names = 1)

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


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain", "Cerebellum", "Hippocampus"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?

annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Cerebellum"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?



annot.GSEA <- easyct(input.d, db="clustermole", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?

## plots



pdf("output/seurat/EasyCellType_dotplot_WT_p14_CB_Rep2_dim30kparam40res07-cellmarker_brainCerebellumHippocampus.pdf", width=6, height=8)
pdf("output/seurat/EasyCellType_dotplot_WT_p14_CB_Rep2_dim30kparam40res07-cellmarker_cerebellum.pdf", width=6, height=8)

pdf("output/seurat/EasyCellType_dotplot_WT_p14_CB_Rep2_dim30kparam40res07-clustermole_brain.pdf", width=6, height=8)

plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()










##########################################
## integration WT Kcnc1 p14 all replicates (1st replicate, then genotype) ######
 ##########################################


WT_p14_CB_Rep1$replicate <- "Rep1"
WT_p14_CB_Rep2$replicate <- "Rep2"
WT_p14_CB_Rep3$replicate <- "Rep3"

WT_p14_CB_Rep1$condition <- "WT"
WT_p14_CB_Rep2$condition <- "WT"
WT_p14_CB_Rep3$condition <- "WT"

Kcnc1_p14_CB_Rep1$replicate <- "Rep1"
Kcnc1_p14_CB_Rep2$replicate <- "Rep2"
Kcnc1_p14_CB_Rep3$replicate <- "Rep3"

Kcnc1_p14_CB_Rep1$condition <- "Kcnc1"
Kcnc1_p14_CB_Rep2$condition <- "Kcnc1"
Kcnc1_p14_CB_Rep3$condition <- "Kcnc1"

set.seed(42)

# Test Replicate integration then genotype integration (2 step integration - with regression repeated for replicates and genotypes integration = 2stepIntegrationRegressRepeated)
## WT Rep
WT_p14_CB_Rep1 <- SCTransform(WT_p14_CB_Rep1, method = "glmGamPoi", ncells = 12877, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 
WT_p14_CB_Rep2 <- SCTransform(WT_p14_CB_Rep2, method = "glmGamPoi", ncells = 13576, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 
WT_p14_CB_Rep3 <- SCTransform(WT_p14_CB_Rep3, method = "glmGamPoi", ncells = 13420, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 


srat.list <- list(WT_p14_CB_Rep1 = WT_p14_CB_Rep1, WT_p14_CB_Rep2 = WT_p14_CB_Rep2, WT_p14_CB_Rep3 = WT_p14_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_p14_CB.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_p14_CB.sct <- IntegrateData(anchorset = WT_p14_CB.anchors, normalization.method = "SCT")

## Kcnc1 Rep
Kcnc1_p14_CB_Rep1 <- SCTransform(Kcnc1_p14_CB_Rep1, method = "glmGamPoi", ncells = 10495, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 
Kcnc1_p14_CB_Rep2 <- SCTransform(Kcnc1_p14_CB_Rep2, method = "glmGamPoi", ncells = 12431, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 
Kcnc1_p14_CB_Rep3 <- SCTransform(Kcnc1_p14_CB_Rep3, method = "glmGamPoi", ncells = 16683, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 

srat.list <- list(Kcnc1_p14_CB_Rep1 = Kcnc1_p14_CB_Rep1, Kcnc1_p14_CB_Rep2 = Kcnc1_p14_CB_Rep2, Kcnc1_p14_CB_Rep3 = Kcnc1_p14_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
Kcnc1_p14_CB.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
Kcnc1_p14_CB.sct <- IntegrateData(anchorset = Kcnc1_p14_CB.anchors, normalization.method = "SCT")


## WT and Kcnc1 

WT_p14_CB.sct <- SCTransform(WT_p14_CB.sct, method = "glmGamPoi", ncells = 39873, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 
Kcnc1_p14_CB.sct <- SCTransform(Kcnc1_p14_CB.sct, method = "glmGamPoi", ncells = 39609, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 


srat.list <- list(WT_p14_CB.sct = WT_p14_CB.sct, Kcnc1_p14_CB.sct = Kcnc1_p14_CB.sct)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

WT_Kcnc1_p14_CB.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_p14_CB.sct <- IntegrateData(anchorset = WT_Kcnc1_p14_CB.anchors, normalization.method = "SCT")

#### UMAP
DefaultAssay(WT_Kcnc1_p14_CB.sct) <- "integrated"

WT_Kcnc1_p14_CB.sct <- RunPCA(WT_Kcnc1_p14_CB.sct, verbose = FALSE, npcs = 30)
WT_Kcnc1_p14_CB.sct <- RunUMAP(WT_Kcnc1_p14_CB.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
WT_Kcnc1_p14_CB.sct <- FindNeighbors(WT_Kcnc1_p14_CB.sct, reduction = "pca", k.param = 40, dims = 1:30)
WT_Kcnc1_p14_CB.sct <- FindClusters(WT_Kcnc1_p14_CB.sct, resolution = 0.7, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p14_CB.sct$condition <- factor(WT_Kcnc1_p14_CB.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1-2stepIntegrationRegressRepeated-dim30kparam40res07.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=TRUE)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_splitCondition-2stepIntegrationRegressRepeated-dim30kparam40res07.pdf", width=13, height=6)
DimPlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=TRUE, split.by = "condition")
dev.off()
pdf("output/seurat/UMAP_WT_Kcnc1_splitReplicate-2stepIntegrationRegressRepeated-dim30kparam40res07.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=TRUE, split.by = "replicate")
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_Phase-2stepIntegrationRegressRepeated-dim30kparam40res07.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p14_CB.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_nFeature_RNA-2stepIntegrationRegressRepeated-dim30kparam40res07.pdf", width=10, height=6)
FeaturePlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_percentmt-2stepIntegrationRegressRepeated-dim30kparam40res07.pdf", width=10, height=6)
FeaturePlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_percentrb-2stepIntegrationRegressRepeated-dim30kparam40res07.pdf", width=10, height=6)
FeaturePlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=FALSE, features = "percent.rb")
dev.off()  




# Test Replicate integration then genotype integration (2 step integration - with regression only for genotype integration = 2stepIntegrationRegressNotRepeated)
## WT Rep
WT_p14_CB_Rep1 <- SCTransform(WT_p14_CB_Rep1, method = "glmGamPoi", ncells = 12877, verbose = TRUE, variable.features.n = 3000) 
WT_p14_CB_Rep2 <- SCTransform(WT_p14_CB_Rep2, method = "glmGamPoi", ncells = 13576, verbose = TRUE, variable.features.n = 3000) 
WT_p14_CB_Rep3 <- SCTransform(WT_p14_CB_Rep3, method = "glmGamPoi", ncells = 13420, verbose = TRUE, variable.features.n = 3000) 


srat.list <- list(WT_p14_CB_Rep1 = WT_p14_CB_Rep1, WT_p14_CB_Rep2 = WT_p14_CB_Rep2, WT_p14_CB_Rep3 = WT_p14_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_p14_CB.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_p14_CB.sct <- IntegrateData(anchorset = WT_p14_CB.anchors, normalization.method = "SCT")

## Kcnc1 Rep
Kcnc1_p14_CB_Rep1 <- SCTransform(Kcnc1_p14_CB_Rep1, method = "glmGamPoi", ncells = 10495, verbose = TRUE, variable.features.n = 3000) 
Kcnc1_p14_CB_Rep2 <- SCTransform(Kcnc1_p14_CB_Rep2, method = "glmGamPoi", ncells = 12431, verbose = TRUE, variable.features.n = 3000) 
Kcnc1_p14_CB_Rep3 <- SCTransform(Kcnc1_p14_CB_Rep3, method = "glmGamPoi", ncells = 16683, verbose = TRUE, variable.features.n = 3000) 

srat.list <- list(Kcnc1_p14_CB_Rep1 = Kcnc1_p14_CB_Rep1, Kcnc1_p14_CB_Rep2 = Kcnc1_p14_CB_Rep2, Kcnc1_p14_CB_Rep3 = Kcnc1_p14_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
Kcnc1_p14_CB.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
Kcnc1_p14_CB.sct <- IntegrateData(anchorset = Kcnc1_p14_CB.anchors, normalization.method = "SCT")


## WT and Kcnc1 

WT_p14_CB.sct <- SCTransform(WT_p14_CB.sct, method = "glmGamPoi", ncells = 39873, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 
Kcnc1_p14_CB.sct <- SCTransform(Kcnc1_p14_CB.sct, method = "glmGamPoi", ncells = 39609, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb"), verbose = TRUE, variable.features.n = 3000) 


srat.list <- list(WT_p14_CB.sct = WT_p14_CB.sct, Kcnc1_p14_CB.sct = Kcnc1_p14_CB.sct)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)

WT_Kcnc1_p14_CB.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_p14_CB.sct <- IntegrateData(anchorset = WT_Kcnc1_p14_CB.anchors, normalization.method = "SCT")

#### UMAP
DefaultAssay(WT_Kcnc1_p14_CB.sct) <- "integrated"

WT_Kcnc1_p14_CB.sct <- RunPCA(WT_Kcnc1_p14_CB.sct, verbose = FALSE, npcs = 30)
WT_Kcnc1_p14_CB.sct <- RunUMAP(WT_Kcnc1_p14_CB.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
WT_Kcnc1_p14_CB.sct <- FindNeighbors(WT_Kcnc1_p14_CB.sct, reduction = "pca", k.param = 40, dims = 1:30)
WT_Kcnc1_p14_CB.sct <- FindClusters(WT_Kcnc1_p14_CB.sct, resolution = 0.7, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p14_CB.sct$condition <- factor(WT_Kcnc1_p14_CB.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1-2stepIntegrationRegressNotRepeated-dim30kparam40res07.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=TRUE)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_splitCondition-2stepIntegrationRegressNotRepeated-dim30kparam40res07.pdf", width=13, height=6)
DimPlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=TRUE, split.by = "condition")
dev.off()
pdf("output/seurat/UMAP_WT_Kcnc1_splitReplicate-2stepIntegrationRegressNotRepeated-dim30kparam40res07.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=TRUE, split.by = "replicate")
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_Phase-2stepIntegrationRegressNotRepeated-dim30kparam40res07.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p14_CB.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_nFeature_RNA-2stepIntegrationRegressNotRepeated-dim30kparam40res07.pdf", width=10, height=6)
FeaturePlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_percentmt-2stepIntegrationRegressNotRepeated-dim30kparam40res07.pdf", width=10, height=6)
FeaturePlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_percentrb-2stepIntegrationRegressNotRepeated-dim30kparam40res07.pdf", width=10, height=6)
FeaturePlot(WT_Kcnc1_p14_CB.sct, reduction = "umap", label=FALSE, features = "percent.rb")
dev.off()  




# Test Replicate and Genotype integration (1 step integration)
## WT Rep

### Regv2 
WT_p14_CB_Rep1 <- SCTransform(WT_p14_CB_Rep1, method = "glmGamPoi", ncells = 12520, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb", "nFeature_RNA")) 
WT_p14_CB_Rep2 <- SCTransform(WT_p14_CB_Rep2, method = "glmGamPoi", ncells = 13474, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb", "nFeature_RNA")) 
WT_p14_CB_Rep3 <- SCTransform(WT_p14_CB_Rep3, method = "glmGamPoi", ncells = 13231, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb", "nFeature_RNA")) 
Kcnc1_p14_CB_Rep1 <- SCTransform(Kcnc1_p14_CB_Rep1, method = "glmGamPoi", ncells = 10410, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb", "nFeature_RNA")) 
Kcnc1_p14_CB_Rep2 <- SCTransform(Kcnc1_p14_CB_Rep2, method = "glmGamPoi", ncells = 11053, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb", "nFeature_RNA")) 
Kcnc1_p14_CB_Rep3 <- SCTransform(Kcnc1_p14_CB_Rep3, method = "glmGamPoi", ncells = 15827, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb", "nFeature_RNA")) 



### Reg v1 _ better than Regv2
WT_p14_CB_Rep1 <- SCTransform(WT_p14_CB_Rep1, method = "glmGamPoi", ncells = 12369, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p14_CB_Rep2 <- SCTransform(WT_p14_CB_Rep2, method = "glmGamPoi", ncells = 13414, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p14_CB_Rep3 <- SCTransform(WT_p14_CB_Rep3, method = "glmGamPoi", ncells = 13181, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p14_CB_Rep1 <- SCTransform(Kcnc1_p14_CB_Rep1, method = "glmGamPoi", ncells = 10382, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p14_CB_Rep2 <- SCTransform(Kcnc1_p14_CB_Rep2, method = "glmGamPoi", ncells = 10934, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p14_CB_Rep3 <- SCTransform(Kcnc1_p14_CB_Rep3, method = "glmGamPoi", ncells = 15577, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 




srat.list <- list(WT_p14_CB_Rep1 = WT_p14_CB_Rep1, WT_p14_CB_Rep2 = WT_p14_CB_Rep2, WT_p14_CB_Rep3 = WT_p14_CB_Rep3, Kcnc1_p14_CB_Rep1 = Kcnc1_p14_CB_Rep1, Kcnc1_p14_CB_Rep2 = Kcnc1_p14_CB_Rep2, Kcnc1_p14_CB_Rep3 = Kcnc1_p14_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_Kcnc1_p14_CB_1step.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_p14_CB_1step.sct <- IntegrateData(anchorset = WT_Kcnc1_p14_CB_1step.anchors, normalization.method = "SCT")


#### UMAP
DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "integrated"

WT_Kcnc1_p14_CB_1step.sct <- RunPCA(WT_Kcnc1_p14_CB_1step.sct, verbose = FALSE, npcs = 30)
WT_Kcnc1_p14_CB_1step.sct <- RunUMAP(WT_Kcnc1_p14_CB_1step.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
WT_Kcnc1_p14_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p14_CB_1step.sct, reduction = "pca", k.param = 50, dims = 1:30)
WT_Kcnc1_p14_CB_1step.sct <- FindClusters(WT_Kcnc1_p14_CB_1step.sct, resolution = 0.35, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p14_CB_1step.sct$condition <- factor(WT_Kcnc1_p14_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim30kparam50res035.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()


# genes

DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "SCT"


pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim40kparam30res05-countMtRbRegression-List4.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = c("Calb1", "Slc1a6", "Car8", "Gabra6", "Pax6", "Sorcs3", "Ptprk", "Nxph1", "Cdh22", "Zeb2", "Hepacam", "Aqp4", "Slc39a12", "Kl", "Clic6", "Slc13a4", "Ttr", "Cfap54", "Ccdc153", "Cfap44", "Tmem212", "Ptgds", "Dcn", "Ntng1", "Grm5", "Aldoc", "Cnp", "Cspg5", "Mbp", "Mag", "Plp1", "Slc18a2", "Ddc", "Slc6a4", "Tph2", "Lef1", "Notum", "Apcdd1", "Nxph1", "Dynlt1c", "Otx1", "Rnd3", "Pvalb", "Cck", "Sst", "Myh11"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()








#### QC metrics investigation ####################################
pdf("output/seurat/VlnPlot_QCmetrics_SCT_WT_Kcnc1_p14_CB_1step-1stepIntegrationRegressNotRepeated-QCV3dim30kparam20res04-countMtRbRegression.pdf", width=20, height=5)
VlnPlot(WT_Kcnc1_p14_CB_1step.sct,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nFeature_RNA_SCT_WT_Kcnc1_p14_CB_1step-1stepIntegrationRegressNotRepeated-QCV3dim30kparam20res04-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p14_CB_1step.sct,features = c("nFeature_RNA")) +
  ylim(0,2000)
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nCount_RNA_SCT_WT_Kcnc1_p14_CB_1step-1stepIntegrationRegressNotRepeated-QCV3dim30kparam20res04-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p14_CB_1step.sct,features = c("nCount_RNA")) +
  ylim(0,10000)
dev.off()


pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p14_CB_nFeature_RNA-1stepIntegrationRegressNotRepeated-QCV3dim30kparam20res04.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p14_CB_percentmt-1stepIntegrationRegressNotRepeated-QCV3dim30kparam20res04.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p14_CB_percentrb-1stepIntegrationRegressNotRepeated-QCV3dim30kparam20res04.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.rb")
dev.off()  
############################################################









pdf("output/seurat/UMAP_WT_Kcnc1_splitCondition-1stepIntegrationRegressNotRepeated-QCV2dim13kparam30res05.pdf", width=13, height=6)
DimPlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "condition")
dev.off()
pdf("output/seurat/UMAP_WT_Kcnc1_splitReplicate-1stepIntegrationRegressNotRepeated-QCV2dim13kparam30res05.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "replicate")
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_Phase-1stepIntegrationRegressNotRepeated-QCV2dim13kparam30res05.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p14_CB_1step.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  


# genes

DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "SCT"

pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-1stepIntegrationRegressNotRepeated-QCV2dim30kparam20res05-countMtRbRegression-List3.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = c("Calb1", "Slc1a6", "Car8", "Gabra6", "Pax6", "Sorcs3", "Ptprk", "Nxph1", "Cdh22", "Zeb2", "Hepacam", "Aqp4", "Slc39a12", "Kl", "Clic6", "Slc13a4", "Ttr", "Cfap54", "Ccdc153", "Cfap44", "Tmem212", "Ptgds", "Dcn", "Ntng1", "Grm5","Aldoc", "Cnp", "Mbp", "Mag", "Plp1", "Slc18a2", "Ddc", "Slc6a4", "Tph2"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()



pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB_Rep2-1stepIntegrationRegressNotRepeated-QCV2dim30kparam40res07-countMtRbRegression-List2.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = c("Calb1", "Slc1a6", "Car8","Gabra6","Sorcs3", "Ptprk","Nxph1", "Cdh22","Edil3","Fabp7", "Zeb2", "Hepacam","Pax6","Gad1", "Gria3","Itpr1","Slc1a2", "Apoe", "Aqp4", "Slc39a12","Glul", "Slc1a3", "Ednrb","Ttr", "Clic6", "Slc13a4", "Kl","Cfap54", "Ccdc153", "Cfap44", "Tmem212","Gad2", "Slc6a1","Ptgds", "Dcn", "Colec12","Ntng1", "Grm5","Slc18a2", "Ddc","Aldoc", "Nnat" ,"Mbp", "Mag", "Plp1", "Pdgfd", "Gli3","Rbpms","S100b", "Cnp","Slc6a4", "Tph2" ), max.cutoff = 1, cols = c("grey", "red"))
dev.off()




# WT vs Kcnc1 gene expr ############
pdf("output/seurat/FeaturePlot_SCT_WT_p14_CB-1stepIntegrationRegressNotRepeated-QCV3dim30kparam50res035-PurkinjeDEGs.pdf", width=10, height=25)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = c("Grid2", "Frmpd4", "Kcnab1", "Gria3", "Garnl3" ,"Dgkb"),  cols = c("grey", "red"), split.by = "condition") #  max.cutoff = 10, min.cutoff = 1
dev.off()



###########################

## investigate to find optimal nb of dimension
### Vizualize the 1st 100 PC
pdf("output/seurat/DimHeatmap_1stepIntegrationRegressNotRepeated.pdf", width=10, height=100)
DimHeatmap(WT_Kcnc1_p14_CB_1step.sct, dims = 1:70, cells = 500, balanced = TRUE)
dev.off()
### Elbow
pdf("output/seurat/Elbow_1stepIntegrationRegressNotRepeated.pdf", width=10, height=10)
ElbowPlot(WT_Kcnc1_p14_CB_1step.sct, ndims = 70) # 8 
dev.off()
### Elbow quantification
#### Determine percent of variation associated with each PC
pct <- WT_Kcnc1_p14_CB_1step.sct[["pca"]]@stdev / sum(WT_Kcnc1_p14_CB_1step.sct[["pca"]]@stdev) * 100
#### Calculate cumulative percents for each PC
cumu <- cumsum(pct)
#### Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1] # 57
co1 # 57
##### Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
#### last point where change of % of variation is more than 0.1%.
co2 # 13
##### Minimum of the two calculation
pcs <- min(co1, co2)
pcs
###### Create a dataframe with values
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
# Elbow plot to visualize 
pdf("output/seurat/ElbowQuantif_1stepIntegrationRegressNotRepeated.pdf", width=5, height=4)
 ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()





# save ##################
## saveRDS(WT_p14_CB.sct, file = "output/seurat/WT_p14_CB.sct_V1_numeric.rds") 
## saveRDS(Kcnc1_p14_CB.sct, file = "output/seurat/Kcnc1_p14_CB.sct_V1_numeric.rds") 
## saveRDS(WT_Kcnc1_p14_CB.sct, file = "output/seurat/WT_Kcnc1_p14_CB.sct_V1_numeric.rds") 
## saveRDS(WT_Kcnc1_p14_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V1_numeric.rds") 
## saveRDS(WT_Kcnc1_p14_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V2_numeric.rds") # regMtRbFeaCount
## WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V1_numeric.rds") # regMtRbCount
## saveRDS(WT_Kcnc1_p14_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V3_numeric.rds") # regMtRbCount with QC_V3
## saveRDS(WT_Kcnc1_p14_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V4_numeric.rds") # regMtRbCount with QC_V3; after Naiara review Goldberg_V2.pptx; QCV3dim30kparam50res035
## saveRDS(WT_Kcnc1_p14_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V5_numeric.rds") # regMtRbCount with QC_V3; after Naiara review Goldberg_V2.pptx; QCV3dim30kparam50res035 with name V2 (corr PLI12 PLI23 instead of interneuron)
WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V4_numeric.rds") # regMtRbCount with QC_V3; after Naiara review Goldberg_V2.pptx; QCV3dim30kparam50res035 without name

WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V5_numeric.rds") # regMtRbCount with QC_V3; after Naiara review Goldberg_V2.pptx; QCV3dim30kparam50res035 with name V2 (corr PLI12 PLI23 instead of interneuron)
set.seed(42)
##########


## Let's work with 1step integration: 
# --> 1stepIntegrationRegressNotRepeatedregMtRbCou-QCV2dim30kparam30res05 = WT_Kcnc1_p14_CB_1step_V1

WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V1_numeric.rds")


############################ EasyCellType automatic annotation ##########################################

### Find all markers 
all_markers <- FindAllMarkers(WT_Kcnc1_p14_CB_1step.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(all_markers, file = "output/seurat/srat_WT_Kcnc1_p14_CB_1step_QCV3_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
## V1 = QCV1
## QCV3


# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_WT_Kcnc1_p14_CB_1step_QCV3_all_markers.txt", header = TRUE, row.names = 1)



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


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain", "Cerebellum", "Hippocampus"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?


annot.GSEA <- easyct(input.d, db="clustermole", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?


## plots

#pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p14_CB_1step_V1-cellmarker_CerebellumBrainHippocampus.pdf", width=6, height=8)
#pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p14_CB_1step_V1-clustermole_Brain.pdf", width=6, height=8)

pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p14_CB_1step_QCV3-clustermole_Brain.pdf", width=6, height=8)
plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# panglao --> NOT working... 'subscript out of bounds' error... I tried gene as ENSEMBL, entrezID and geneSymbo, human/mic, everything...


############ V1 naming

Cluster1: Granular_1
Cluster2: MLI1
Cluster3: Neuroblast_1
Cluster4: Granular_2
Cluster5: Granular_3
Cluster6: Granular_4
Cluster7: Granular_5
Cluster8: MLI2
Cluster9: Granular_6
Cluster10: Endothelial_Cells
Cluster11: Granular_7
Cluster12: Astrocyte
Cluster13: OPC
Cluster14: Interneuron_1 (mature?)
Cluster15: Bergmann_Glia
Cluster16: Interneuron_2
Cluster17: Oligodendrocyte
Cluster18: EpendymalMeningeal_Cells
Cluster19: Muscle_Cells
Cluster20: Purkinje_Cells
Cluster21: Neuroblast_2 (late?)
Cluster22: NPC
Cluster23: Choroid_Plexus_Cells



new.cluster.ids <- c(
"Granular_1",
"MLI1",
"Neuroblast_1",
"Granular_2",
"Granular_3",
"Granular_4",
"Granular_5",
"MLI2",
"Granular_6",
"Endothelial_Cells",
"Granular_7",
"Astrocyte",
"OPC",
"Interneuron_1",
"Bergmann_Glia",
"Interneuron_2",
"Oligodendrocyte",
"EpendymalMeningeal_Cells",
"Muscle_Cells",
"Purkinje_Cells",
"Neuroblast_2" ,
"NPC",
"Choroid_Plexus_Cells"
)

names(new.cluster.ids) <- levels(WT_Kcnc1_p14_CB_1step.sct)
WT_Kcnc1_p14_CB_1step.sct <- RenameIdents(WT_Kcnc1_p14_CB_1step.sct, new.cluster.ids)

WT_Kcnc1_p14_CB_1step.sct$cluster.annot <- Idents(WT_Kcnc1_p14_CB_1step.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam20res04_label.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()


pdf("output/seurat/UMAP_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam20res04_noSplit_label.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5)
dev.off()



# All in dotplot
DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "SCT"

List5:
1,4,5,6,7,9,11- Granular = Gabra6, Pax6
3,21- Neuroblast= Ntng1, Grm5
2- MLI1 = Sorcs3, Ptprk
8- MLI2 = Nxph1, Cdh22
10- Endothelial_Cells= Lef1, Notum, Apcdd1
12- Astrocyte = Zeb2, Hepacam
12- Astrocyte_2 = Aqp4, Slc39a12
13- OPC = Aldoc, Cnp, Cspg5
14- Interneuron= Rgs6, Kcnc2
17- Oligodendrocyte= Mbp, Mag, Plp1
18- Ependymal cells= Cfap54, Ccdc153, Cfap44, Tmem212
18- Meningeal cell = Ptgds, Dcn
19- Muscle= Dlc1, Pdgfrb
20- Purkinje cells = Calb1, Slc1a6, Car8
22- NPC= Lmx1a, Adcy2, Mdga1, Eomes
23- Choroid plexus cells = Kl, Clic6, Slc13a4, Ttr


all_markers <- c(
"Gabra6", "Pax6", "Ntng1", "Grm5", "Sorcs3", "Ptprk", "Nxph1", "Cdh22", "Lef1", "Notum", "Apcdd1", "Zeb2", "Hepacam", "Aqp4", "Slc39a12", "Aldoc", "Cnp", "Cspg5", "Rgs6", "Kcnc2", "Mbp", "Mag", "Plp1", "Cfap54", "Ccdc153", "Cfap44", "Tmem212", "Ptgds", "Dcn", "Dlc1", "Pdgfrb", "Calb1", "Slc1a6", "Car8", "Lmx1a", "Adcy2", "Mdga1", "Eomes", "Kl", "Clic6", "Slc13a4", "Ttr"
)



levels(WT_Kcnc1_p14_CB_1step.sct) <- c(
"Granular_1",
"Granular_2",
"Granular_3",
"Granular_4",
"Granular_5",
"Granular_6",
"Granular_7",
"Neuroblast_1",
"Neuroblast_2",
"MLI1",
"MLI2",
"Endothelial_Cells",
"Astrocyte",
"Bergmann_Glia",
"OPC",
"Interneuron_1",
"Interneuron_2",
"Oligodendrocyte",
"EpendymalMeningeal_Cells",
"Muscle_Cells",
"Purkinje_Cells",
"NPC",
"Choroid_Plexus_Cells"
)



pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam20res04_label.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p14_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam20res04_label_vertical.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p14_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()





pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p14_CB_1step-QCV3dim30kparam20res04-List5.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = c("Gabra6", "Pax6", "Ntng1", "Grm5", "Sorcs3", "Ptprk", "Nxph1", "Cdh22", "Lef1", "Notum", "Apcdd1", "Zeb2", "Hepacam", "Aqp4", "Slc39a12", "Aldoc", "Cnp", "Cspg5", "Rgs6", "Kcnc2", "Mbp", "Mag", "Plp1", "Cfap54", "Ccdc153", "Cfap44", "Tmem212", "Ptgds", "Dcn", "Dlc1", "Pdgfrb", "Calb1", "Slc1a6", "Car8", "Lmx1a", "Adcy2", "Mdga1", "Eomes", "Kl", "Clic6", "Slc13a4", "Ttr"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()








############ V2 naming (WT_Kcnc1_p14_CB_1step.sct_V4_numeric.rds = QCV3dim30kparam50res035)

Cluster1: Interneuron is instead PLI23
Cluster2: Granular
Cluster3: MLI1
Cluster4: Granular
Cluster5: Granular
Cluster6: Granular
Cluster7: MLI2
Cluster8: Endothelial
Cluster9: Granular
Cluster10: Astrocyte
Cluster11: OPC
Cluster12: Bergmann_Glia
Cluster13: PLI is instead PLI12
Cluster14: Oligodendrocyte
Cluster15: Mix_Microglia_Meningeal
Cluster16: Endothelial_Mural 
Cluster17: Purkinje
Cluster18: Golgi
Cluster19: Unipolar_Brush
Cluster20: Choroid_Plexus

new.cluster.ids <- c(
  "PLI23",
  "Granular_1",
  "MLI1",
  "Granular_2",
  "Granular_3",
  "Granular_4",
  "MLI2",
  "Endothelial",
  "Granular_5",
  "Astrocyte",
  "OPC",
  "Bergmann_Glia",
  "PLI12",
  "Oligodendrocyte",
  "Mix_Microglia_Meningeal",
  "Endothelial_Mural" ,
  "Purkinje",
  "Golgi",
  "Unipolar_Brush",
  "Choroid_Plexus"
)

names(new.cluster.ids) <- levels(WT_Kcnc1_p14_CB_1step.sct)
WT_Kcnc1_p14_CB_1step.sct <- RenameIdents(WT_Kcnc1_p14_CB_1step.sct, new.cluster.ids)
WT_Kcnc1_p14_CB_1step.sct$cluster.annot <- Idents(WT_Kcnc1_p14_CB_1step.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam50res035_label.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam50res035_noSplit_label.pdf", width=9, height=6)
DimPlot(WT_Kcnc1_p14_CB_1step.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5)
dev.off()



# All in dotplot
DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "SCT"


List6:
Granular = Gabra6, Pax6
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
PLI12 = Klhl1, Gfra2, Aldh1a3
PLI23 = Galntl6, Kcnc2
Golgi = Pax2
Unipolar_Brush = Eomes
Purkinje = Calb1, Slc1a6, Car8
Astrocyte = Zeb2
Bergman_Glia = Aqp4, Slc39a12
Oligodendrocyte= Mbp, Mag, Plp1
OPC = Aldoc, Cnp
Microglia = Itgam, Cx3cr1
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Endothelial_Mural =  Dlc1, Pdgfrb
Choroid plexus cells = Kl,  Ttr


all_markers <- c(
  "Gabra6","Pax6",
  "Sorcs3", "Ptprk",
  "Nxph1", "Cdh22",
  "Klhl1", "Gfra2", "Aldh1a3",
  "Galntl6", "Kcnc2",
  "Pax2",
  "Eomes",
  "Calb1", "Slc1a6", "Car8",
  "Zeb2",
  "Aqp4", "Slc39a12",
  "Mbp", "Mag", "Plp1",
  "Aldoc", "Cnp",
  "Itgam", "Cx3cr1",
  "Ptgds", "Dcn",
  "Lef1", "Notum", "Apcdd1",
  "Dlc1", "Pdgfrb",
  "Kl",  "Ttr"
)



levels(WT_Kcnc1_p14_CB_1step.sct) <- c(
  "Granular_1",
  "Granular_2",
  "Granular_3",
  "Granular_4",
  "Granular_5",
  "MLI1" ,
  "MLI2" ,
  "PLI12" ,
  "PLI23",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergmann_Glia",
  "Oligodendrocyte",
  "OPC" ,
  "Mix_Microglia_Meningeal",
  "Endothelial",
  "Endothelial_Mural",
  "Choroid_Plexus"
)



pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam50res035_label.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p14_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam50res035_label_vertical.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p14_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()




pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p14_CB_1step-QCV3dim30kparam50res035-List6.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = all_markers, max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p14_CB_1step-QCV3dim30kparam50res035-Kcnc1.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = "Kcnc1", cols = c("grey", "red"), max.cutoff = 1)
dev.off()

## p14 cell type proportion ###############################
### count nb of cells in each cluster
WT_p14_CB_Rep1 = table(Idents(WT_Kcnc1_p14_CB_1step.sct)[WT_Kcnc1_p14_CB_1step.sct$orig.ident == "WT_p14_CB_Rep1"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p14",
             replicate= "Rep1")
WT_p14_CB_Rep2 = table(Idents(WT_Kcnc1_p14_CB_1step.sct)[WT_Kcnc1_p14_CB_1step.sct$orig.ident == "WT_p14_CB_Rep2"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p14",
             replicate= "Rep2")
WT_p14_CB_Rep3 = table(Idents(WT_Kcnc1_p14_CB_1step.sct)[WT_Kcnc1_p14_CB_1step.sct$orig.ident == "WT_p14_CB_Rep3"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p14",
             replicate= "Rep3")

Kcnc1_p14_CB_Rep1 = table(Idents(WT_Kcnc1_p14_CB_1step.sct)[WT_Kcnc1_p14_CB_1step.sct$orig.ident == "Kcnc1_p14_CB_Rep1"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p14",
             replicate= "Rep1")
Kcnc1_p14_CB_Rep2 = table(Idents(WT_Kcnc1_p14_CB_1step.sct)[WT_Kcnc1_p14_CB_1step.sct$orig.ident == "Kcnc1_p14_CB_Rep2"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p14",
             replicate= "Rep2")
Kcnc1_p14_CB_Rep3 = table(Idents(WT_Kcnc1_p14_CB_1step.sct)[WT_Kcnc1_p14_CB_1step.sct$orig.ident == "Kcnc1_p14_CB_Rep3"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p14",
             replicate= "Rep3")


p14_CB = WT_p14_CB_Rep1 %>%
  bind_rows(WT_p14_CB_Rep2) %>%
  bind_rows(WT_p14_CB_Rep3) %>%
  bind_rows(Kcnc1_p14_CB_Rep1) %>%
  bind_rows(Kcnc1_p14_CB_Rep2) %>%
  bind_rows(Kcnc1_p14_CB_Rep3) %>%
  as_tibble()
  



### Keeping all replicates
p14_CB_prop = p14_CB %>%
  group_by(replicate, genotype) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = (count / total_count) * 100)

p14_CB_prop$genotype <-
  factor(p14_CB_prop$genotype,
         c("WT", "Kcnc1"))

pdf("output/seurat/histogramProp_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam50res035.pdf", width=7, height=4)
ggbarplot(p14_CB_prop, x = "cluster", y = "proportion", fill = "genotype",
                  color = "genotype", palette = c("black", "blue"),
                  position = position_dodge(0.8), # Separate bars by genotype
                  add = "mean_se", # Add error bars
                  lab.pos = "out", lab.size = 3) +
  stat_compare_means(aes(group = genotype), method = "t.test", label = "p.signif") +
  theme_bw() +
  labs(x = "Cell Type (Cluster)", y = "Cell Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




### Outlier removed (WT_Rep2)

p14_CB_filt = WT_p14_CB_Rep1 %>%
  bind_rows(WT_p14_CB_Rep3) %>%
  bind_rows(Kcnc1_p14_CB_Rep1) %>%
  bind_rows(Kcnc1_p14_CB_Rep2) %>%
  bind_rows(Kcnc1_p14_CB_Rep3) %>%
  as_tibble()
  


p14_CB_filt_prop = p14_CB_filt %>%
  group_by(replicate, genotype) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = (count / total_count) * 100)

p14_CB_filt_prop$genotype <-
  factor(p14_CB_filt_prop$genotype,
         c("WT", "Kcnc1"))

pdf("output/seurat/histogramProp_WT_Kcnc1_p14_CB_1step_filtOutlier_QCV3dim30kparam50res035.pdf", width=7, height=4)
ggbarplot(p14_CB_filt_prop, x = "cluster", y = "proportion", fill = "genotype",
                  color = "genotype", palette = c("black", "blue"),
                  position = position_dodge(0.8), # Separate bars by genotype
                  add = "mean_se", # Add error bars
                  lab.pos = "out", lab.size = 3) +
  stat_compare_means(aes(group = genotype), method = "t.test", label = "p.signif") + # or p.format
  theme_bw() +
  labs(x = "Cell Type (Cluster)", y = "Cell Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
# --> This is very long and has been run in slurm job




# DEGs number colored in a UMAP
Idents(WT_Kcnc1_p14_CB_1step.sct) <- "cluster.annot"

DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())
## List of cell types
cell_types <- c(    "PLI23",
  "Granular_1",
  "MLI1",
  "Granular_2",
  "Granular_3",
  "Granular_4",
  "MLI2",
  "Endothelial",
  "Granular_5",
  "Astrocyte",
  "OPC",
  "Bergmann_Glia",
  "PLI12",
  "Oligodendrocyte",
  "Mix_Microglia_Meningeal",
  "Endothelial_Mural" ,
  "Purkinje",
  "Golgi",
  "Unipolar_Brush",
  "Choroid_Plexus")
## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/seurat/", cell_type, "-Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes_correct.txt", sep = "")
  deg_data <- read.table(file_name, header = TRUE, sep = "\t") ## Read the DEGs data
  num_degs <- sum(deg_data$p_val_adj < 0.05 & (deg_data$avg_log2FC > 0.25 | deg_data$avg_log2FC < -0.25)) ## Count the number of significant DEGs
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))  ## Append to the summary table
}

DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = c(  "Granular_1",
  "Granular_2",
  "Granular_3",
  "Granular_4",
  "Granular_5",
  "MLI1" ,
  "MLI2" ,
  "PLI12" ,
  "PLI23",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergmann_Glia",
  "Oligodendrocyte",
  "OPC" ,
  "Mix_Microglia_Meningeal",
  "Endothelial",
  "Endothelial_Mural",
  "Choroid_Plexus")) 
  
  
# Add DEG information to my seurat object - DEG_count
cell_clusters <- WT_Kcnc1_p14_CB_1step.sct@meta.data$cluster.annot
names(cell_clusters) <- rownames(WT_Kcnc1_p14_CB_1step.sct@meta.data)
DEG_named_vector <- DEG_count$Num_DEGs[match(cell_clusters, DEG_count$Cell_Type)]
names(DEG_named_vector) <- names(cell_clusters)
# Integrate DEG values into the Seurat object
WT_Kcnc1_p14_CB_1step.sct <- AddMetaData(WT_Kcnc1_p14_CB_1step.sct, metadata = DEG_named_vector, col.name = "DEG")
# Create a UMAP plot colored by qval values
pdf("output/seurat/FeaturePlot_WT_Kcnc1_p14_CB_1step_DEG.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = "DEG", pt.size = 0.5, reduction = "umap", label = TRUE) +
  scale_colour_viridis() #  option="magma"
dev.off()
# Add values on the heatmap
## Extract UMAP coordinates
umap_coordinates <- as.data.frame(WT_Kcnc1_p14_CB_1step.sct@reductions$umap@cell.embeddings)
umap_coordinates$cluster <- WT_Kcnc1_p14_CB_1step.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DEG_count %>% dplyr::rename( "cluster"="Cell_Type"))
## Create a UMAP plot colored by DEG values, with cluster DEG counts as text annotations
pdf("output/seurat/FeaturePlot_WT_Kcnc1_p14_CB_1step_DEGpadj05fc025_numeric.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p14_CB_1step.sct, features = "DEG", pt.size = 0.5, reduction = "umap") +
  scale_colour_viridis(option="mako") + # option="magma"
  geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = Num_DEGs), 
            size = 4, color = "red", fontface = "bold") 
dev.off()



# SCPA

# Compare WT and cYAPKO using SCPA ##########################################

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
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")

DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "RNA" # Recommended 



# Test different Pathway collections and generate enrichment plot for each cell types (C2 = Pathway, C5 = ontology )
## import Pathways
pathways <- msigdbr("Mus musculus", "C5") %>%          # !!!!!! CHANGE HERE PATHWAYS !!!!!!
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
"Granular_1", "Granular_2", "Granular_3", "Granular_4", "Granular_5",
  "Interneuron", "MLI1", "MLI2", "PLI", "Golgi", "Unipolar_Brush",
  "Purkinje", "Astrocyte", "Bergmann_Glia", "Oligodendrocyte", "OPC",
  "Mix_Microglia_Meningeal", "Endothelial", "Endothelial_Mural", "Choroid_Plexus"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(WT_Kcnc1_p14_CB_1step.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  Kcnc1 <- seurat_extract(WT_Kcnc1_p14_CB_1step.sct,
                           meta1 = "condition", value_meta1 = "Kcnc1",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, Kcnc1),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_SB_p14_C5_", cluster, ".txt")       # !!!!!! CHANGE HERE PATHWAYS !!!!!!
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}
#--> Long ~2hrs



## load all the comparison for each cell type (FC qval information)
clusters = c(
"Granular_1", "Granular_2", "Granular_3", "Granular_4", "Granular_5",
  "Interneuron", "MLI1", "MLI2", "PLI", "Golgi", "Unipolar_Brush",
  "Purkinje", "Astrocyte", "Bergmann_Glia", "Oligodendrocyte", "OPC",
  "Mix_Microglia_Meningeal", "Endothelial", "Endothelial_Mural", "Choroid_Plexus"
)
## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_CB_p14_C2_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster)) %>% as_tibble()

## Filter pathway of interest
pathways <- c(
  "WP_NEURODEGENERATION_WITH_BRAIN_IRON_ACCUMULATION_NBIA_SUBTYPES_PATHWAY",
  "WP_NEUROINFLAMMATION_AND_GLUTAMATERGIC_SIGNALING",
  "KEGG_ALZHEIMERS_DISEASE",
  "WP_ALZHEIMERS_DISEASE",
  "KEGG_PARKINSONS_DISEASE",
  "WP_PARKINSONS_DISEASE_PATHWAY",
  "KEGG_HUNTINGTONS_DISEASE",
  "WP_ERK_PATHWAY_IN_HUNTINGTONS_DISEASE",
  "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
  "WP_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
  "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
  "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS",
  "REACTOME_POTASSIUM_CHANNELS",
  "REACTOME_ION_CHANNEL_TRANSPORT",
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  "REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS",
  "REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_DOPAMINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_SEROTONIN_NEUROTRANSMITTER_RELEASE_CYCLE",
  "ALCALA_APOPTOSIS",
  "KEGG_APOPTOSIS",
  "REACTOME_APOPTOSIS",
  "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
  "WP_APOPTOSIS",
  "WP_APOPTOSIS_MODULATION_AND_SIGNALING",
  "REACTOME_DEATH_RECEPTOR_SIGNALLING",
  "REACTOME_DISEASES_OF_PROGRAMMED_CELL_DEATH",
  "REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES",
  "REACTOME_PROGRAMMED_CELL_DEATH",
  "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
  "BIOCARTA_DEATH_PATHWAY",
  "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES",
  "WP_PKCGAMMA_CALCIUM_SIGNALING_PATHWAY_IN_ATAXIA"
)
classes <- c(
  rep("Neurodegeneration", 10),
  rep("Neuronal_Activity", 13),
  rep("Apoptosis", 12),
  "ROS",
  "Ataxia"
)
colors <- c(
  rep("Purple", 10),
  rep("Green", 13),
  rep("Dark Orange", 12),
  "Red",
  "Blue"
)
pathway_tibble <- tibble(
  Pathway = pathways,
  Class = classes,
  Color = colors
)

pathway_all_data = pathway_tibble %>%
  left_join(all_data)

pathway_all_data$Pathway <- factor(pathway_all_data$Pathway, levels = pathways)


### dotplot
pdf("output/Pathway/dotplot_SCPA_CB_p14_C2_selectedPathwayV1.pdf", width=12, height=8)
ggplot(pathway_all_data, aes(x = cluster, y = Pathway)) +
  geom_point(aes(size = ifelse(qval > 1.4, qval, NA), color = Color), na.rm = TRUE) +
  scale_size_continuous(range = c(3, 10), breaks = c(1.5, 2, 3, 4), name = "qval") +
  scale_color_identity() +
  theme_bw() +
  labs(x = "Cluster",
       y = "Pathway",
       color = "Class Color") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "right") +
  guides(color = guide_legend(title = "Class Color", override.aes = list(size = 5)),
         size = guide_legend(title = "qval"))
dev.off()




# GSEA plot
library("fgsea")


#### import all clsuter DEGs output :
cluster_types <- c("Granular_1",
  "Granular_2",
  "Granular_3",
  "Granular_4",
  "Granular_5",
  "MLI1" ,
  "MLI2" ,
  "PLI12" ,
  "PLI23",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergmann_Glia",
  "Oligodendrocyte",
  "OPC" ,
  "Mix_Microglia_Meningeal",
  "Endothelial",
  "Endothelial_Mural",
  "Choroid_Plexus")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes_correct.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}

## load list of genes to test
### List1
fgsea_sets <- list(
  PathwaysOfNeurodegeneration = read_table(file = "output/Pathway/geneList_mmu05022.txt")$Genes,
  GOBP_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS = read_table(file = "output/Pathway/geneList_GOBP_NEURON_INTRINSIC_APOPTOTIC.txt")$Genes,
  GOBP_NEUROINFLAMMATORY_RESPONSE = read_table(file = "output/Pathway/geneList_GOBP_NEUROINFLAMMATORY_RESPONSE.txt")$Genes
)
### List2
fgsea_sets <- list(
  GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION = read_table(file = "output/Pathway/geneList_GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION.txt")$Genes,
  GOBP_NEURONAL_ACTION_POTENTIAL = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_ACTION_POTENTIAL.txt")$Genes,
  GOBP_NEURONAL_SIGNAL_TRANSDUCTION = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_SIGNAL_TRANSDUCTION.txt")$Genes,
  REACTOME_POTASSIUM_CHANNELS = read_table(file = "output/Pathway/geneList_REACTOME_POTASSIUM_CHANNELS.txt")$Genes,
  REACTOME_ION_CHANNEL_TRANSPORT = read_table(file = "output/Pathway/geneList_REACTOME_ION_CHANNEL_TRANSPORT.txt")$Genes
)
### List3
fgsea_sets <- list(
  REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING = read_table(file = "output/Pathway/geneList_REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING.txt")$Genes,
  REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS = read_table(file = "output/Pathway/geneList_REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS.txt")$Genes,
  REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION = read_table(file = "output/Pathway/geneList_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION.txt")$Genes,
  REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS = read_table(file = "output/Pathway/geneList_REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS.txt")$Genes
)
### List4
fgsea_sets <- list(
  PathwaysOfNeurodegeneration = read_table(file = "output/Pathway/geneList_mmu05022.txt")$Genes,
  GOBP_NEURONAL_ACTION_POTENTIAL = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_ACTION_POTENTIAL.txt")$Genes,
  REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION = read_table(file = "output/Pathway/geneList_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION.txt")$Genes
)


## Rank genes based on FC
genes <- Purkinje %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(ES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -NES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()


## plot GSEA
pdf("output/Pathway/GSEA_Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes-mmu05022-Purkinje.pdf", width=5, height=3)
plotEnrichment(fgsea_sets[["PathwaysOfNeurodegeneration"]],
               ranks) + labs(title="PathwaysOfNeurodegeneration-Purkinje") +
               theme_bw()
dev.off()


# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c("Granular_1",
  "Granular_2",
  "Granular_3",
  "Granular_4",
  "Granular_5",
  "MLI1" ,
  "MLI2" ,
  "PLI12" ,
  "PLI23",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergmann_Glia",
  "Oligodendrocyte",
  "OPC" ,
  "Mix_Microglia_Meningeal",
  "Endothelial",
  "Endothelial_Mural",
  "Choroid_Plexus")

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
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 10000)
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

write.table(final_results, file = c("output/Pathway/gsea_output_Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes_correct-List4.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap all GSEA
pdf("output/Pathway/heatmap_gsea_output_Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes_correct-List4.pdf", width=10, height=3)
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

### if need reorder terms
final_results$pathway <- factor(final_results$pathway, levels = c("PathwaysOfNeurodegeneration", "GOBP_NEURONAL_ACTION_POTENTIAL", "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION")) 

pdf("output/Pathway/heatmap_gsea_output_Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes_correct-List4.pdf", width=10, height=3)
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






##########################################
## integration WT Kcnc1 p35 all replicates (1st replicate, then genotype) ######
 ##########################################


WT_p35_CB_Rep1$replicate <- "Rep1"
WT_p35_CB_Rep2$replicate <- "Rep2"
WT_p35_CB_Rep3$replicate <- "Rep3"

WT_p35_CB_Rep1$condition <- "WT"
WT_p35_CB_Rep2$condition <- "WT"
WT_p35_CB_Rep3$condition <- "WT"

Kcnc1_p35_CB_Rep1$replicate <- "Rep1"
Kcnc1_p35_CB_Rep2$replicate <- "Rep2"   #   Kcnc1_p35_CB_Rep2_nCountRNA500
Kcnc1_p35_CB_Rep3$replicate <- "Rep3"

Kcnc1_p35_CB_Rep1$condition <- "Kcnc1"
Kcnc1_p35_CB_Rep2$condition <- "Kcnc1" # Kcnc1_p35_CB_Rep2_nCountRNA500 Kcnc1_p35_CB_Rep2
Kcnc1_p35_CB_Rep3$condition <- "Kcnc1"

set.seed(42)


# Replicate and Genotype integration (1 step integration)
## WT Rep

### Reg v1 _ better than Regv2 (OG that worked but I cannot repeat...)
WT_p35_CB_Rep1 <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 8496, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2 <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 11577, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3 <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 14582, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1 <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 11642, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2 <- SCTransform(Kcnc1_p35_CB_Rep2, method = "glmGamPoi", ncells = 33698, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3 <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 18205, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 


### Reg v1 _ better than Regv2 - Correct for bug RNA quantity leading to downregulation of all genes _ QCV3
WT_p35_CB_Rep1 <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 7286, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2 <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 10660, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3 <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 13642, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1 <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 10211, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2 <- SCTransform(Kcnc1_p35_CB_Rep2, method = "glmGamPoi", ncells = 26802, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3 <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 16324, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 


### Reg v1 _ better than Regv2 - Correct for bug RNA quantity leading to downregulation of all genes _ QCV2 # GOOD ONE!!
WT_p35_CB_Rep1 <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 7299, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2 <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 10683, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3 <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 13664, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1 <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 10264, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2 <- SCTransform(Kcnc1_p35_CB_Rep2, method = "glmGamPoi", ncells = 33623, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3 <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 16447, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
#--> Kcnc1 Rep2 is problematic, many cells, but many with very low level of nCount_RNA, it disturb the normalization and integration... Let's try to subset this sample, to have approximately the same number of cells (increase QC filtering) with same nCount_RNA as the other for this specific sample.. (~7k cells is ideal)

Kcnc1_p35_CB_Rep2_nCountRNA1000 <- subset(Kcnc1_p35_CB_Rep2, subset = nCount_RNA > 1000) # 1533 cells
Kcnc1_p35_CB_Rep2_nCountRNA700 <- subset(Kcnc1_p35_CB_Rep2, subset = nCount_RNA > 700) # 3279 cells
Kcnc1_p35_CB_Rep2_nCountRNA400 <- subset(Kcnc1_p35_CB_Rep2, subset = nCount_RNA > 400) # 11535 cells

Kcnc1_p35_CB_Rep2_nCountRNA500 <- subset(Kcnc1_p35_CB_Rep2, subset = nCount_RNA > 500) # 6526 cells; lets use this one!!
#--> Repeat SCTrnsform for that sample and reperfrorm integration


### Reg v1 _ better than Regv2 - Correct for bug RNA quantity leading to downregulation of all genes _ QCV2 with Kcnc1 V2 filter
DefaultAssay(Kcnc1_p35_CB_Rep2_nCountRNA500) <- "RNA"


WT_p35_CB_Rep1 <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 7299, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2 <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 10683, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3 <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 13664, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1 <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 10264, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2 <- SCTransform(Kcnc1_p35_CB_Rep2_nCountRNA500, method = "glmGamPoi", ncells = 6526, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3 <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 16447, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 






srat.list <- list(WT_p35_CB_Rep1 = WT_p35_CB_Rep1, WT_p35_CB_Rep2 = WT_p35_CB_Rep2, WT_p35_CB_Rep3 = WT_p35_CB_Rep3, Kcnc1_p35_CB_Rep1 = Kcnc1_p35_CB_Rep1, Kcnc1_p35_CB_Rep2 = Kcnc1_p35_CB_Rep2, Kcnc1_p35_CB_Rep3 = Kcnc1_p35_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_Kcnc1_p35_CB_1step.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_p35_CB_1step.sct <- IntegrateData(anchorset = WT_Kcnc1_p35_CB_1step.anchors, normalization.method = "SCT")


#### UMAP
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "integrated"

WT_Kcnc1_p35_CB_1step.sct <- RunPCA(WT_Kcnc1_p35_CB_1step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_p35_CB_1step.sct <- RunUMAP(WT_Kcnc1_p35_CB_1step.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_Kcnc1_p35_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p35_CB_1step.sct, reduction = "pca", k.param = 20, dims = 1:50)
WT_Kcnc1_p35_CB_1step.sct <- FindClusters(WT_Kcnc1_p35_CB_1step.sct, resolution = 0.3, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p35_CB_1step.sct$condition <- factor(WT_Kcnc1_p35_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st


# pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim50kparam50res03.pdf", width=7, height=6)
# pdf("output/seurat/UMAP_WT_Kcnc1_p35_CBtest.pdf", width=7, height=6)

pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou_Kcnc1Rep2nCountRNA500-QCV2dim50kparam20res03.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()

################################################################################
## TRY HARMONY ########################################################################################
# install.packages("harmony")
library("harmony")
### Reg v1 _ better than Regv2 - Correct for bug RNA quantity leading to downregulation of all genes _ QCV2
WT_p35_CB_Rep1 <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 7299, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2 <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 10683, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3 <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 13664, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1 <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 10264, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2 <- SCTransform(Kcnc1_p35_CB_Rep2, method = "glmGamPoi", ncells = 33623, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3 <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 16447, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 

srat.list <- list(WT_p35_CB_Rep1 = WT_p35_CB_Rep1, WT_p35_CB_Rep2 = WT_p35_CB_Rep2, WT_p35_CB_Rep3 = WT_p35_CB_Rep3, Kcnc1_p35_CB_Rep1 = Kcnc1_p35_CB_Rep1, Kcnc1_p35_CB_Rep2 = Kcnc1_p35_CB_Rep2, Kcnc1_p35_CB_Rep3 = Kcnc1_p35_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)


srat.list <- list(WT_p35_CB_Rep1 = WT_p35_CB_Rep1, WT_p35_CB_Rep2 = WT_p35_CB_Rep2, WT_p35_CB_Rep3 = WT_p35_CB_Rep3, Kcnc1_p35_CB_Rep1 = Kcnc1_p35_CB_Rep1, Kcnc1_p35_CB_Rep3 = Kcnc1_p35_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
#-->Same shit downregulation all genes...
srat.list <- list(WT_p35_CB_Rep1 = WT_p35_CB_Rep1, Kcnc1_p35_CB_Rep1 = Kcnc1_p35_CB_Rep1)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
#-->Same shit downregulation all genes...

srat.list <- list(WT_p35_CB_Rep2 = WT_p35_CB_Rep2, WT_p35_CB_Rep3 = WT_p35_CB_Rep3, Kcnc1_p35_CB_Rep1 = Kcnc1_p35_CB_Rep1, Kcnc1_p35_CB_Rep3 = Kcnc1_p35_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)

## Merge dataset, and not integrate!
WT_Kcnc1_p35_CB_1step.sct <- merge(x = srat.list[[1]], y = srat.list[2:length(srat.list)], merge.data=TRUE)
VariableFeatures(WT_Kcnc1_p35_CB_1step.sct) <- features

WT_Kcnc1_p35_CB_1step.sct <- RunPCA(WT_Kcnc1_p35_CB_1step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_p35_CB_1step.sct <- RunHarmony(WT_Kcnc1_p35_CB_1step.sct, assay.use="SCT", group.by.vars = "condition") # group.by.vars is the experimental covariates to be corrected/harmonized by the algorithm.
WT_Kcnc1_p35_CB_1step.sct <- RunUMAP(WT_Kcnc1_p35_CB_1step.sct, reduction = "harmony", dims = 1:50)
WT_Kcnc1_p35_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p35_CB_1step.sct, reduction = "harmony", k.param = 20, dims = 1:50)
WT_Kcnc1_p35_CB_1step.sct <- FindClusters(WT_Kcnc1_p35_CB_1step.sct, resolution = 0.3, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p35_CB_1step.sct$condition <- factor(WT_Kcnc1_p35_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV2dim50kparam20res03Harmony.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()
# --> Same puttin orig.ident in RunHarmony(*, group.by.vars)... A lot of downreagulated genes for cluster1 (Granular)...
########################################################################################################################
########################################################################################################################




################################################################################
## TRY LOGNORM WITH SCT ########################################################################################
### Reg v1 _ better than Regv2 - Correct for bug RNA quantity leading to downregulation of all genes _ QCV2
## NormalizeData
WT_p35_CB_Rep1 <- NormalizeData(WT_p35_CB_Rep1, normalization.method = "LogNormalize", scale.factor = 10000) 
WT_p35_CB_Rep2 <- NormalizeData(WT_p35_CB_Rep2, normalization.method = "LogNormalize", scale.factor = 10000) 
WT_p35_CB_Rep3 <- NormalizeData(WT_p35_CB_Rep3, normalization.method = "LogNormalize", scale.factor = 10000) 
Kcnc1_p35_CB_Rep1 <- NormalizeData(Kcnc1_p35_CB_Rep1, normalization.method = "LogNormalize", scale.factor = 10000) 
Kcnc1_p35_CB_Rep2 <- NormalizeData(Kcnc1_p35_CB_Rep2, normalization.method = "LogNormalize", scale.factor = 10000) 
Kcnc1_p35_CB_Rep3 <- NormalizeData(Kcnc1_p35_CB_Rep3, normalization.method = "LogNormalize", scale.factor = 10000) 

## SCT
WT_p35_CB_Rep1 <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 7299, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2 <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 10683, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3 <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 13664, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1 <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 10264, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2 <- SCTransform(Kcnc1_p35_CB_Rep2, method = "glmGamPoi", ncells = 33623, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3 <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 16447, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 

srat.list <- list(WT_p35_CB_Rep1 = WT_p35_CB_Rep1, WT_p35_CB_Rep2 = WT_p35_CB_Rep2, WT_p35_CB_Rep3 = WT_p35_CB_Rep3, Kcnc1_p35_CB_Rep1 = Kcnc1_p35_CB_Rep1, Kcnc1_p35_CB_Rep2 = Kcnc1_p35_CB_Rep2, Kcnc1_p35_CB_Rep3 = Kcnc1_p35_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_Kcnc1_p35_CB_1step.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_p35_CB_1step.sct <- IntegrateData(anchorset = WT_Kcnc1_p35_CB_1step.anchors, normalization.method = "SCT")


DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "integrated"

WT_Kcnc1_p35_CB_1step.sct <- RunPCA(WT_Kcnc1_p35_CB_1step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_p35_CB_1step.sct <- RunUMAP(WT_Kcnc1_p35_CB_1step.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_Kcnc1_p35_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p35_CB_1step.sct, reduction = "pca", k.param = 20, dims = 1:50)
WT_Kcnc1_p35_CB_1step.sct <- FindClusters(WT_Kcnc1_p35_CB_1step.sct, resolution = 0.3, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p35_CB_1step.sct$condition <- factor(WT_Kcnc1_p35_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st


pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou_LogNormSCT-QCV2dim50kparam20res03.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()

########################################################################################################################
########################################################################################################################



# CORRECT VERSION (issue found for the global downregulation)
WT_p35_CB_Rep1 <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 7299, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2 <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 10683, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3 <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 13664, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1 <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 10264, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2 <- SCTransform(Kcnc1_p35_CB_Rep2, method = "glmGamPoi", ncells = 33623, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3 <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 16447, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 



srat.list <- list(WT_p35_CB_Rep1 = WT_p35_CB_Rep1, WT_p35_CB_Rep2 = WT_p35_CB_Rep2, WT_p35_CB_Rep3 = WT_p35_CB_Rep3, Kcnc1_p35_CB_Rep1 = Kcnc1_p35_CB_Rep1, Kcnc1_p35_CB_Rep2 = Kcnc1_p35_CB_Rep2, Kcnc1_p35_CB_Rep3 = Kcnc1_p35_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_Kcnc1_p35_CB_1step.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_p35_CB_1step.sct <- IntegrateData(anchorset = WT_Kcnc1_p35_CB_1step.anchors, normalization.method = "SCT")


#### UMAP
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "integrated"

WT_Kcnc1_p35_CB_1step.sct <- RunPCA(WT_Kcnc1_p35_CB_1step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_p35_CB_1step.sct <- RunUMAP(WT_Kcnc1_p35_CB_1step.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_Kcnc1_p35_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p35_CB_1step.sct, reduction = "pca", k.param = 20, dims = 1:50)
WT_Kcnc1_p35_CB_1step.sct <- FindClusters(WT_Kcnc1_p35_CB_1step.sct, resolution = 0.3, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p35_CB_1step.sct$condition <- factor(WT_Kcnc1_p35_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st



# pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV2dim50kparam20res03.pdf", width=7, height=6)
pdf("output/seurat/UMAP_test.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()


## THIS VERSION HAS BEEN SAVED AS:  output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V2_label_ReProcess.rds










# genes

DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "SCT"


# pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim50-List6.pdf", width=30, height=60)
pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV2dim50-List7.pdf", width=30, height=60)

FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = c(  "Gabra6", "Pax6",
  "Kcnc2",
  "Sorcs3", "Ptprk",
  "Nxph1", "Cdh22",
  "Aldh1a3",
  "Pax2",
  "Eomes", "Rgs6", "Tafa2",
  "Calb1", "Slc1a6", "Car8",
  "Zeb2",
  "Aqp4", "Slc39a12",
  "Vcan", "Sox6",
  "Ptgds", "Dcn",
  "Lef1", "Notum", "Apcdd1",
  "Kl",  "Ttr",
  "Actb", "Tmsb4x"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()







#### QC metrics investigation ####################################
pdf("output/seurat/VlnPlot_QCmetrics_SCT_WT_Kcnc1_p35_CB_1step-1stepIntegrationRegressNotRepeated-QCV2dim50kparam20res03-countMtRbRegression.pdf", width=20, height=5)
VlnPlot(WT_Kcnc1_p35_CB_1step.sct,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nFeature_RNA_SCT_WT_Kcnc1_p35_CB_1step-1stepIntegrationRegressNotRepeated-QCV2dim50kparam20res03-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p35_CB_1step.sct,features = c("nFeature_RNA")) +
  ylim(0,2000)
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nCount_RNA_SCT_WT_Kcnc1_p35_CB_1step-1stepIntegrationRegressNotRepeated-QCV2dim50kparam20res03-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p35_CB_1step.sct,features = c("nCount_RNA")) +
  ylim(0,10000)
dev.off()


pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p35_CB_nFeature_RNA-1stepIntegrationRegressNotRepeated-QCV2dim50.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p35_CB_percentmt-1stepIntegrationRegressNotRepeated-QCV2dim50.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p35_CB_percentrb-1stepIntegrationRegressNotRepeated-QCV2dim50.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.rb")
dev.off()  
############################################################



pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB_splitCondition-1stepIntegrationRegressNotRepeated-QCV3dim50kparam50res03.pdf", width=13, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "condition")
dev.off()
pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB_splitReplicate-1stepIntegrationRegressNotRepeated-QCV3dim50kparam50res03.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "replicate")
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_WT_p35_CB_Kcnc1_Phase-1stepIntegrationRegressNotRepeated-QCV3dim50kparam50res03.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  


# save ##################
## saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step.sct_V1_numeric.rds") # regMtRbCount with QC_V3
# WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step.sct_V1_numeric.rds")
## saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_numeric.rds") # regMtRbCount with QC_V3
## saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_label.rds") # regMtRbCount with QC_V3
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_numeric.rds")

WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_label.rds")
## Let's work with 1step integration: 
# --> 1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim50kparam50res03 = WT_Kcnc1_p35_CB_1step_V1
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_label.rds")

# Save the Kcnc1_p35_CB_Rep2 with additional QC filtering: too many cells with low nCount_RNA so I filtered to nCount_RNA > 500; and end up with 6526 cells instead of ~30k cells saveRDS(Kcnc1_p35_CB_Rep2_nCountRNA500, file = "output/seurat/Kcnc1_p35_CB_Rep2_nCountRNA500.rds")

Kcnc1_p35_CB_Rep2_nCountRNA500 <- readRDS(file = "output/seurat/Kcnc1_p35_CB_Rep2_nCountRNA500.rds")

# saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step_Kcnc1Rep2nCountRNA500-QCV2dim50kparam20res03.sct_V1_numeric_ReProcess.rds") #

WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step_Kcnc1Rep2nCountRNA500-QCV2dim50kparam20res03.sct_V1_numeric_ReProcess.rds")

set.seed(42)


## BUG I cannot retrieve previous version... Lets compare QC_V2 vs V3 vs V4 with same parameters
#saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_numeric_ReProcess.rds") # this is with QC_V3 all repeated from the start!!! --> NOT the one...
#saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V1_numeric_ReProcess.rds") # this is with QC_V2 all repeated from the start!!! --> NOT the one...
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V1_numeric_ReProcess.rds")
#saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V2_label_ReProcess.rds") 

WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V2_label_ReProcess.rds") # THIS IS THE LAST ONE!


##########




############################ EasyCellType automatic annotation ##########################################

### Find all markers 
all_markers <- FindAllMarkers(WT_Kcnc1_p35_CB_1step.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(all_markers, file = "output/seurat/srat_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03_all_markers.txt", header = TRUE, row.names = 1)



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


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain", "Cerebellum", "Hippocampus"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?


annot.GSEA <- easyct(input.d, db="clustermole", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?





## plots



#pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p35_CB_1step_V1-cellmarker_CerebellumBrainHippocampus.pdf", width=6, height=8)
#pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p35_CB_1step_V1-clustermole_Brain.pdf", width=6, height=8)

pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03-clustermole_Brain.pdf", width=6, height=8)
plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




# panglao --> NOT working... 'subscript out of bounds' error... I tried gene as ENSEMBL, entrezID and geneSymbo, human/mic, everything...


#




## check some genes



DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "SCT"


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV2dim50-SuppTab1Marker.pdf", width=30, height=40)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = c( "Ppp1r17", "Gabra6", "Grm2", "Sst","Prkcd", "Sorcs3", "Ptprk", "Nxph1", "Cdh22","Htr2a", "Edil3","Aldh1a3", "Slc6a5","Eomes","Gdf10"
), max.cutoff = 1, cols = c("grey", "red"))
dev.off()

pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV2dim50-oligodendrocyt.pdf", width=30, height=20)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = c( "Mbp", "Mag", "Plp1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV2dim50-PLImarkers.pdf", width=30, height=30)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = c("Klhl1", "Nrp1", "Gfra2", "Slc16a12", "Aldh1a3", "Slc6a5", "Galntl6", "Scg2", "Fstl5", "Zeb2", "Kcnc2", "Htr2a"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()




############ V1 naming (output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_numeric.rds = QCV3dim50kparam50res03)

Cluster1: Granular
Cluster2: MLI1
Cluster3: Interneuron
Cluster4: Endothelial_Stalk
Cluster5: MLI2
Cluster6: Astrocyte
Cluster7: PLI
Cluster8: Golgi
Cluster9: Bergman_Glia
Cluster10: UBC
Cluster11: Endothelial
Cluster12: Choroid_Plexus
Cluster13: Purkinje
Cluster14: Meningeal
Cluster15: OPC


new.cluster.ids <- c(
  "Granular",
  "MLI1",
  "Interneuron",
  "Endothelial_Stalk",
  "MLI2",
  "Astrocyte",
  "PLI",
  "Golgi",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Endothelial",
  "Choroid_Plexus",
  "Purkinje",
  "Meningeal",
  "OPC"
)

names(new.cluster.ids) <- levels(WT_Kcnc1_p35_CB_1step.sct)
WT_Kcnc1_p35_CB_1step.sct <- RenameIdents(WT_Kcnc1_p35_CB_1step.sct, new.cluster.ids)
WT_Kcnc1_p35_CB_1step.sct$cluster.annot <- Idents(WT_Kcnc1_p35_CB_1step.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03_label.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03_noSplit_label.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5)
dev.off()



# All in dotplot
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "SCT"


List7:
Granular = Gabra6, Pax6
Interneuron = Kcnc2
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
PLI = Aldh1a3
Golgi = Pax2
Unipolar_Brush = Eomes, Rgs6, Tafa2
Purkinje = Calb1, Slc1a6, Car8
Astrocyte = Zeb2
Bergman_Glia = Aqp4, Slc39a12
OPC = Vcan, Sox6
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Choroid plexus cells = Kl,  Ttr
Endothelial_Stalk = Actb, Tmsb4x



all_markers <- c(
  "Gabra6", "Pax6",
  "Kcnc2",
  "Sorcs3", "Ptprk",
  "Nxph1", "Cdh22",
  "Aldh1a3",
  "Pax2",
  "Eomes", "Rgs6", "Tafa2",
  "Calb1", "Slc1a6", "Car8",
  "Zeb2",
  "Aqp4", "Slc39a12",
  "Vcan", "Sox6",
  "Ptgds", "Dcn",
  "Lef1", "Notum", "Apcdd1",
  "Kl",  "Ttr",
  "Actb", "Tmsb4x"
)



levels(WT_Kcnc1_p35_CB_1step.sct) <- c(
  "Granular",
  "Interneuron",
  "MLI1",
  "MLI2",
  "PLI",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergman_Glia",
  "OPC",
  "Meningeal",
  "Endothelial",
  "Choroid_Plexus",
  "Endothelial_Stalk"
)



pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03_label.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p35_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03_label_vertical.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p35_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()




pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03-List7.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = all_markers, max.cutoff = 1, cols = c("grey", "red"))
dev.off()



pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03-Kcnc1.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = "Kcnc1", max.cutoff = 1, cols = c("grey", "red"))
dev.off()






############ V2 naming ("output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V1_numeric_ReProcess.rds = QCV2dim50kparam20res03))
Cluster1= Granule (Gabra6 , Pax6)
Cluster2= MLI1 (Sorcs3, Ptprk)
Cluster3= PLI23 (Galntl6, Kcnc2)
Cluster4= MLI2 (Nxph1, Cdh22)
Cluster5= Endothelial_Stalk (Actb, Tmsb4x)
Cluster6= Astrocyte (Zeb2)
Cluster7= PLI12 (Klhl1, Gfra2, Aldh1a3)
Cluster8= Golgi (Pax2)
Cluster9= Unipolar_Brush (Eomes, Rgs6, Tafa2)
Cluster10= Bergman_Glia (Aqp4, Slc39a12)
Cluster11= Endothelial (Lef1, Notum, Apcdd1)
Cluster12= Choroid_Plexus (Kl, Ttr)
Cluster13= Purkinje (Calb1, Slc1a6, Car8)
Cluster14= Meningeal (Ptgds, Dcn)
Cluster15= Unknown_Neuron_Subpop
Cluster16= OPC (Vcan, Sox6)


new.cluster.ids <- c(
  "Granule",
  "MLI1",
  "PLI23" ,
  "MLI2" ,
  "Endothelial_Stalk",
  "Astrocyte" ,
  "PLI12" ,
  "Golgi" ,
  "Unipolar_Brush" ,
  "Bergman_Glia",
  "Endothelial",
  "Choroid_Plexus",
  "Purkinje",
  "Meningeal",
  "Unknown_Neuron_Subpop",
  "OPC" 
)

names(new.cluster.ids) <- levels(WT_Kcnc1_p35_CB_1step.sct)
WT_Kcnc1_p35_CB_1step.sct <- RenameIdents(WT_Kcnc1_p35_CB_1step.sct, new.cluster.ids)
WT_Kcnc1_p35_CB_1step.sct$cluster.annot <- Idents(WT_Kcnc1_p35_CB_1step.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03_label.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03_noSplit_label.pdf", width=8, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5)
dev.off()



# All in dotplot
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "SCT"


List10:
Granule =Gabra6 , Pax6
MLI1 =Sorcs3, Ptprk
MLI2 =Nxph1, Cdh22
PLI12 =Klhl1, Gfra2, Aldh1a3
PLI23 =Galntl6, Kcnc2
Golgi =Pax2
Unipolar_Brush =Eomes, Rgs6, Tafa2
Purkinje =Calb1, Slc1a6, Car8
Astrocyte =Zeb2
Bergman_Glia =Aqp4, Slc39a12
OPC= Vcan, Sox6
Meningeal =Ptgds, Dcn
Choroid_Plexus =Kl, Ttr
Endothelial =Lef1, Notum, Apcdd1
Endothelial_Stalk =Actb, Tmsb4x
Unknown_Neuron_Subpop


all_markers <- c(
  "Gabra6" , "Pax6",
  "Sorcs3", "Ptprk",
  "Nxph1", "Cdh22",
  "Klhl1", "Gfra2", "Aldh1a3",
  "Galntl6", "Kcnc2",
  "Pax2",
  "Eomes", "Rgs6", "Tafa2",
  "Calb1", "Slc1a6", "Car8",
  "Zeb2",
  "Aqp4", "Slc39a12",
  "Vcan", "Sox6",
  "Ptgds", "Dcn",
  "Kl", "Ttr",
  "Lef1", "Notum", "Apcdd1",
  "Actb", "Tmsb4x"
)



levels(WT_Kcnc1_p35_CB_1step.sct) <- c(
  "Granule",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergman_Glia",
  "OPC",
  "Meningeal",
  "Choroid_Plexus",
  "Endothelial",
  "Endothelial_Stalk",
  "Unknown_Neuron_Subpop"
)



pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03_label.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p35_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03_label_vertical.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p35_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()




pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03-List10.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = all_markers, max.cutoff = 1, cols = c("grey", "red"))
dev.off()



pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03-Kcnc1.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = "Kcnc1", max.cutoff = 1, cols = c("grey", "red"))
dev.off()






## p35 cell type proportion ############################### 
### count nb of cells in each cluster
WT_p35_CB_Rep1 = table(Idents(WT_Kcnc1_p35_CB_1step.sct)[WT_Kcnc1_p35_CB_1step.sct$orig.ident == "WT_p35_CB_Rep1"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p35",
             replicate= "Rep1")
WT_p35_CB_Rep2 = table(Idents(WT_Kcnc1_p35_CB_1step.sct)[WT_Kcnc1_p35_CB_1step.sct$orig.ident == "WT_p35_CB_Rep2"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p35",
             replicate= "Rep2")
WT_p35_CB_Rep3 = table(Idents(WT_Kcnc1_p35_CB_1step.sct)[WT_Kcnc1_p35_CB_1step.sct$orig.ident == "WT_p35_CB_Rep3"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p35",
             replicate= "Rep3")

Kcnc1_p35_CB_Rep1 = table(Idents(WT_Kcnc1_p35_CB_1step.sct)[WT_Kcnc1_p35_CB_1step.sct$orig.ident == "Kcnc1_p35_CB_Rep1"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p35",
             replicate= "Rep1")
Kcnc1_p35_CB_Rep2 = table(Idents(WT_Kcnc1_p35_CB_1step.sct)[WT_Kcnc1_p35_CB_1step.sct$orig.ident == "Kcnc1_p35_CB_Rep2"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p35",
             replicate= "Rep2")
Kcnc1_p35_CB_Rep3 = table(Idents(WT_Kcnc1_p35_CB_1step.sct)[WT_Kcnc1_p35_CB_1step.sct$orig.ident == "Kcnc1_p35_CB_Rep3"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p35",
             replicate= "Rep3")


p35_CB = WT_p35_CB_Rep1 %>%
  bind_rows(WT_p35_CB_Rep2) %>%
  bind_rows(WT_p35_CB_Rep3) %>%
  bind_rows(Kcnc1_p35_CB_Rep1) %>%
  bind_rows(Kcnc1_p35_CB_Rep2) %>%
  bind_rows(Kcnc1_p35_CB_Rep3) %>%
  as_tibble()
  



### Keeping all replicates
p35_CB_prop = p35_CB %>%
  group_by(replicate, genotype) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = (count / total_count) * 100)

p35_CB_prop$genotype <-
  factor(p35_CB_prop$genotype,
         c("WT", "Kcnc1"))

pdf("output/seurat/histogramProp_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03.pdf", width=7, height=4)
ggbarplot(p35_CB_prop, x = "cluster", y = "proportion", fill = "genotype",
                  color = "genotype", palette = c("black", "blue"),
                  position = position_dodge(0.8), # Separate bars by genotype
                  add = "mean_se", # Add error bars
                  lab.pos = "out", lab.size = 3) +
  stat_compare_means(aes(group = genotype), method = "t.test", label = "p.signif") +
  theme_bw() +
  labs(x = "Cell Type (Cluster)", y = "Cell Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



### Outlier removed (WT and Kcnc1_Rep2)
XXX NOT MODIFFED

p35_CB_filt = WT_p35_CB_Rep1 %>%
  bind_rows(WT_p35_CB_Rep3) %>%
  bind_rows(Kcnc1_p35_CB_Rep1) %>%
  bind_rows(Kcnc1_p35_CB_Rep3) %>%
  as_tibble()
  


p35_CB_filt_prop = p35_CB_filt %>%
  group_by(replicate, genotype) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = (count / total_count) * 100)

p35_CB_filt_prop$genotype <-
  factor(p35_CB_filt_prop$genotype,
         c("WT", "Kcnc1"))

pdf("output/seurat/histogramProp_WT_Kcnc1_p35_CB_1step_filtOutlier_QCV3dim50kparam50res03.pdf", width=7, height=4)
ggbarplot(p35_CB_filt_prop, x = "cluster", y = "proportion", fill = "genotype",
                  color = "genotype", palette = c("black", "blue"),
                  position = position_dodge(0.8), # Separate bars by genotype
                  add = "mean_se", # Add error bars
                  lab.pos = "out", lab.size = 3) +
  stat_compare_means(aes(group = genotype), method = "t.test", label = "p.signif") + # or p.format
  theme_bw() +
  labs(x = "Cell Type (Cluster)", y = "Cell Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()







# differential expressed genes across conditions ##########################
############## Default method ##########
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "RNA"

WT_Kcnc1_p35_CB_1step.sct$celltype.stim <- paste(WT_Kcnc1_p35_CB_1step.sct$cluster.annot, WT_Kcnc1_p35_CB_1step.sct$condition,
    sep = "-")
Idents(WT_Kcnc1_p35_CB_1step.sct) <- "celltype.stim"

# Define the list of clusters for comparison
clusters <- c(
  "Granule",
  "MLI1",
  "PLI23" ,
  "MLI2" ,
  "Endothelial_Stalk",
  "Astrocyte" ,
  "PLI12" ,
  "Golgi" ,
  "Unipolar_Brush" ,
  "Bergman_Glia",
  "Endothelial",
  "Choroid_Plexus",
  "Purkinje",
  "Meningeal",
  "Unknown_Neuron_Subpop",
  "OPC" 
)



# Initialize an empty list to store the results for each cluster
cluster_markers <- list()
# Loop through each cluster and perform FindMarkers
for (cluster in clusters) {
  cat("Processing cluster:", cluster, "\n")
  # Run FindMarkers for each cluster comparing Kcnc1 vs WT condition
  markers <- FindMarkers(WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = paste(cluster, "Kcnc1", sep = "-"), 
                         ident.2 = paste(cluster, "WT", sep = "-"), 
                         verbose = TRUE, 
                         test.use = "wilcox",
                         logfc.threshold = -Inf,
                         min.pct = -Inf,
                         min.diff.pct = -Inf,
                         assay = "RNA")
  # Store the result in the list
  cluster_markers[[cluster]] <- markers
  # Save the result as a text file
  output_filename <- paste0("output/seurat/", cluster, "-Kcnc1_response_p35_CB_QCV3dim50kparam20res03Reprocess_allGenes.txt")
  write.table(markers, file = output_filename, sep = "\t", quote = FALSE, row.names = TRUE)
}
# --> Too long run in slurm job


## Slight test other DEG method (bimod, roc, t, negbinom, poisson, LR, MAST)
Idents(WT_Kcnc1_p35_CB_1step.sct) <- "celltype.stim"

Purkinje_bimod <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT",
                         test.use = "bimod")
OPC_bimod <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "OPC-Kcnc1", 
                         ident.2 = "OPC-WT",
                         test.use = "bimod")
Meningeal_bimod <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Meningeal-Kcnc1", 
                         ident.2 = "Meningeal-WT",
                         test.use = "bimod")
cat("Upregulated:", sum(Meningeal_bimod$p_val_adj < 0.05 & Meningeal_bimod$avg_log2FC > 0), 
    "\nDownregulated:", sum(Meningeal_bimod$p_val_adj < 0.05 & Meningeal_bimod$avg_log2FC < 0), "\n")
Purkinje_roc <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT",
                         test.use = "roc")
cat("Upregulated:", sum(Purkinje_roc$p_val_adj < 0.05 & Purkinje_roc$avg_log2FC > 0), 
    "\nDownregulated:", sum(Purkinje_roc$p_val_adj < 0.05 & Purkinje_roc$avg_log2FC < 0), "\n")
Purkinje_t <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT",
                         test.use = "t")
cat("Upregulated:", sum(Purkinje_t$p_val_adj < 0.05 & Purkinje_t$avg_log2FC > 0), 
    "\nDownregulated:", sum(Purkinje_t$p_val_adj < 0.05 & Purkinje_t$avg_log2FC < 0), "\n")
Purkinje_negbinom <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT",
                         test.use = "negbinom")
cat("Upregulated:", sum(Purkinje_negbinom$p_val_adj < 0.05 & Purkinje_negbinom$avg_log2FC > 0), 
    "\nDownregulated:", sum(Purkinje_negbinom$p_val_adj < 0.05 & Purkinje_negbinom$avg_log2FC < 0), "\n")
Purkinje_poisson <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT",
                         test.use = "poisson")
cat("Upregulated:", sum(Purkinje_poisson$p_val_adj < 0.05 & Purkinje_poisson$avg_log2FC > 0), 
    "\nDownregulated:", sum(Purkinje_poisson$p_val_adj < 0.05 & Purkinje_poisson$avg_log2FC < 0), "\n")
Purkinje_LR <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT",
                         test.use = "LR")
cat("Upregulated:", sum(Purkinje_LR$p_val_adj < 0.05 & Purkinje_LR$avg_log2FC > 0), 
    "\nDownregulated:", sum(Purkinje_LR$p_val_adj < 0.05 & Purkinje_LR$avg_log2FC < 0), "\n")
Purkinje_MAST <- FindMarkers(object = WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Purkinje-Kcnc1", 
                         ident.2 = "Purkinje-WT",
                         test.use = "MAST")
cat("Upregulated:", sum(Purkinje_MAST$p_val_adj < 0.05 & Purkinje_MAST$avg_log2FC > 0), 
    "\nDownregulated:", sum(Purkinje_MAST$p_val_adj < 0.05 & Purkinje_MAST$avg_log2FC < 0), "\n")
#--> Same.. all method give many downregulated genes..


############## Perform DE analysis after pseudobulking ##########

### Transform seurat to pseudobulk seurat
pseudo_WT_Kcnc1_p35_CB_1step.sct <- AggregateExpression(WT_Kcnc1_p35_CB_1step.sct, assays = "RNA", return.seurat = T, group.by = c("condition", "orig.ident", "cluster.annot"))
### each 'cell' is a genotype-condition-celltype pseudobulk profile
tail(Cells(pseudo_WT_Kcnc1_p35_CB_1step.sct))
head(Cells(pseudo_WT_Kcnc1_p35_CB_1step.sct))
### Extract cell names and update the pseudo seurat with $informations
#### Extract cell names
cell_names <- Cells(pseudo_WT_Kcnc1_p35_CB_1step.sct)
#### Split the cell names into parts
#### The split pattern below splits into 3 parts: condition, orig.ident, and cluster.annot
cell_parts <- strsplit(cell_names, "_", fixed = TRUE)
#### Extract condition, orig.ident, and cluster.annot from each cell name
condition <- sapply(cell_parts, function(x) x[1]) # The first part is the condition
orig_ident <- sapply(cell_parts, function(x) paste(x[2:length(x) - 1], collapse = "_")) # Middle parts are orig.ident, which may contain underscores
cluster_annot <- sapply(cell_parts, function(x) x[length(x)]) # The last part is the cluster annotation
#### Add extracted data as metadata to the Seurat object
pseudo_WT_Kcnc1_p35_CB_1step.sct$condition <- condition
pseudo_WT_Kcnc1_p35_CB_1step.sct$orig.ident <- orig_ident
pseudo_WT_Kcnc1_p35_CB_1step.sct$cluster.annot <- cluster_annot

### Prepare and perform pseudobulk DEG
pseudo_WT_Kcnc1_p35_CB_1step.sct$celltype.stim <- paste(pseudo_WT_Kcnc1_p35_CB_1step.sct$cluster.annot, pseudo_WT_Kcnc1_p35_CB_1step.sct$condition, sep = "_")

Idents(pseudo_WT_Kcnc1_p35_CB_1step.sct) <- "celltype.stim"
DefaultAssay(pseudo_WT_Kcnc1_p35_CB_1step.sct) <- "RNA"

### round the counts to be used by deseq2 (require integer)
#### Access the raw counts matrix in the RNA assay
counts_matrix <- GetAssayData(pseudo_WT_Kcnc1_p35_CB_1step.sct, slot = "counts", assay = "RNA")
#### Round the counts to the nearest integer
rounded_counts_matrix <- round(counts_matrix)
#### Replace the counts in the Seurat object with the rounded counts
pseudo_WT_Kcnc1_p35_CB_1step.sct <- SetAssayData(pseudo_WT_Kcnc1_p35_CB_1step.sct, slot = "counts", new.data = rounded_counts_matrix, assay = "RNA")



Purkinje_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Purkinje_Kcnc1", 
                         ident.2 = "Purkinje_WT",
                         test.use = "DESeq2")
Granule_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Granule_Kcnc1", 
                         ident.2 = "Granule_WT",
                         test.use = "DESeq2")
Astrocyte_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Astrocyte_Kcnc1", 
                         ident.2 = "Astrocyte_WT",
                         test.use = "DESeq2")
Meningeal_pseudobulk <- FindMarkers(object = pseudo_WT_Kcnc1_p35_CB_1step.sct, 
                         ident.1 = "Meningeal_Kcnc1", 
                         ident.2 = "Meningeal_WT",
                         test.use = "DESeq2")
cat("Upregulated:", sum(Meningeal_pseudobulk$p_val_adj < 0.05 & Meningeal_pseudobulk$avg_log2FC > 0), 
    "\nDownregulated:", sum(Meningeal_pseudobulk$p_val_adj < 0.05 & Meningeal_pseudobulk$avg_log2FC < 0), "\n")
#-=-> DESEQ2 to be run on the raw counts of RNA assay!                         
#--> Few to no DEGs









# DEGs number colored in a UMAP
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V2_label_ReProcess.rds")

Idents(WT_Kcnc1_p35_CB_1step.sct) <- "cluster.annot"

DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())
## List of cell types
cell_types <- c(     "Granule",
  "MLI1",
  "PLI23" ,
  "MLI2" ,
  "Endothelial_Stalk",
  "Astrocyte" ,
  "PLI12" ,
  "Golgi" ,
  "Unipolar_Brush" ,
  "Bergman_Glia",
  "Endothelial",
  "Choroid_Plexus",
  "Purkinje",
  "Meningeal",
  "Unknown_Neuron_Subpop",
  "OPC" )
## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/seurat/", cell_type, "-Kcnc1_response_p35_CB_QCV2dim50kparam20res03_allGenes_correct.txt", sep = "")
  deg_data <- read.table(file_name, header = TRUE, sep = "\t") ## Read the DEGs data
  num_degs <- sum(deg_data$p_val_adj < 0.05 & (deg_data$avg_log2FC > 0.25 | deg_data$avg_log2FC < -0.25) ) ## Count the number of significant DEGs
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))  ## Append to the summary table
}




DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = c(   "Granule",
  "MLI1",
  "PLI23" ,
  "MLI2" ,
  "Endothelial_Stalk",
  "Astrocyte" ,
  "PLI12" ,
  "Golgi" ,
  "Unipolar_Brush" ,
  "Bergman_Glia",
  "Endothelial",
  "Choroid_Plexus",
  "Purkinje",
  "Meningeal",
  "Unknown_Neuron_Subpop",
  "OPC")) 
# Add DEG information to my seurat object - DEG_count
cell_clusters <- WT_Kcnc1_p35_CB_1step.sct@meta.data$cluster.annot
names(cell_clusters) <- rownames(WT_Kcnc1_p35_CB_1step.sct@meta.data)
DEG_named_vector <- DEG_count$Num_DEGs[match(cell_clusters, DEG_count$Cell_Type)]
names(DEG_named_vector) <- names(cell_clusters)
# Integrate DEG values into the Seurat object
WT_Kcnc1_p35_CB_1step.sct <- AddMetaData(WT_Kcnc1_p35_CB_1step.sct, metadata = DEG_named_vector, col.name = "DEG")
# Create a UMAP plot colored by qval values
pdf("output/seurat/FeaturePlot_WT_Kcnc1_p35_CB_1step_DEG_correct.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = "DEG", pt.size = 0.5, reduction = "umap", label = TRUE) +
  scale_colour_viridis() #  option="magma"
dev.off()
# Add values on the heatmap
## Extract UMAP coordinates
umap_coordinates <- as.data.frame(WT_Kcnc1_p35_CB_1step.sct@reductions$umap@cell.embeddings)
umap_coordinates$cluster <- WT_Kcnc1_p35_CB_1step.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DEG_count %>% dplyr::rename( "cluster"="Cell_Type"))
## Create a UMAP plot colored by DEG values, with cluster DEG counts as text annotations
pdf("output/seurat/FeaturePlot_WT_Kcnc1_p35_CB_1step_DEGpadj05fc025_numeric.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = "DEG", pt.size = 0.5, reduction = "umap") +
  scale_colour_viridis(option="mako") + # option="magma"
  geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = Num_DEGs), 
            size = 4, color = "red", fontface = "bold") 
dev.off()




# SCPA

# Compare WT and cYAPKO using SCPA ##########################################

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
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")

DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "RNA" # Recommended 



# Test different Pathway collections and generate enrichment plot for each cell types (C2 = Pathway, C5 = ontology )
## import Pathways
pathways <- msigdbr("Mus musculus", "C5") %>%          # !!!!!! CHANGE HERE PATHWAYS !!!!!!
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
 "Granular",
  "Interneuron",
  "MLI1",
  "MLI2_1",
  "MLI2_2",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergman_Glia",
  "OPC",
  "Meningeal",
  "Endothelial",
  "Choroid_Plexus",
  "Endothelial_Stalk"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(WT_Kcnc1_p35_CB_1step.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  Kcnc1 <- seurat_extract(WT_Kcnc1_p35_CB_1step.sct,
                           meta1 = "condition", value_meta1 = "Kcnc1",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, Kcnc1),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_CB_p35_C5_", cluster, ".txt")       # !!!!!! CHANGE HERE PATHWAYS !!!!!!
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}
#--> Long ~2hrs



## load all the comparison for each cell type (FC qval information)
clusters = c(
 "Granular",
  "Interneuron",
  "MLI1",
  "MLI2_1",
  "MLI2_2",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergman_Glia",
  "OPC",
  "Meningeal",
  "Endothelial",
  "Choroid_Plexus",
  "Endothelial_Stalk"
)
## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_CB_p35_C2_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster)) %>% as_tibble()

## Filter pathway of interest
pathways <- c(
  "WP_NEURODEGENERATION_WITH_BRAIN_IRON_ACCUMULATION_NBIA_SUBTYPES_PATHWAY",
  "WP_NEUROINFLAMMATION_AND_GLUTAMATERGIC_SIGNALING",
  "KEGG_ALZHEIMERS_DISEASE",
  "WP_ALZHEIMERS_DISEASE",
  "KEGG_PARKINSONS_DISEASE",
  "WP_PARKINSONS_DISEASE_PATHWAY",
  "KEGG_HUNTINGTONS_DISEASE",
  "WP_ERK_PATHWAY_IN_HUNTINGTONS_DISEASE",
  "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
  "WP_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
  "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
  "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS",
  "REACTOME_POTASSIUM_CHANNELS",
  "REACTOME_ION_CHANNEL_TRANSPORT",
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  "REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS",
  "REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_DOPAMINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_SEROTONIN_NEUROTRANSMITTER_RELEASE_CYCLE",
  "ALCALA_APOPTOSIS",
  "KEGG_APOPTOSIS",
  "REACTOME_APOPTOSIS",
  "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
  "WP_APOPTOSIS",
  "WP_APOPTOSIS_MODULATION_AND_SIGNALING",
  "REACTOME_DEATH_RECEPTOR_SIGNALLING",
  "REACTOME_DISEASES_OF_PROGRAMMED_CELL_DEATH",
  "REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES",
  "REACTOME_PROGRAMMED_CELL_DEATH",
  "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
  "BIOCARTA_DEATH_PATHWAY",
  "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES",
  "WP_PKCGAMMA_CALCIUM_SIGNALING_PATHWAY_IN_ATAXIA"
)
classes <- c(
  rep("Neurodegeneration", 10),
  rep("Neuronal_Activity", 13),
  rep("Apoptosis", 12),
  "ROS",
  "Ataxia"
)
colors <- c(
  rep("Purple", 10),
  rep("Green", 13),
  rep("Dark Orange", 12),
  "Red",
  "Blue"
)
pathway_tibble <- tibble(
  Pathway = pathways,
  Class = classes,
  Color = colors
)

pathway_all_data = pathway_tibble %>%
  left_join(all_data)

pathway_all_data$Pathway <- factor(pathway_all_data$Pathway, levels = pathways)



### dotplot
pdf("output/Pathway/dotplot_SCPA_CB_p35_C2_selectedPathwayV1.pdf", width=12, height=8)
ggplot(pathway_all_data %>%
         mutate(cluster = ifelse(cluster == "MLI2_2", "PLI", cluster)) %>%
         mutate(cluster = ifelse(cluster == "MLI2_1", "MLI2", cluster)),
       aes(x = cluster, y = Pathway)) +
  geom_point(aes(size = ifelse(qval > 1.4, qval, NA), color = Color), na.rm = TRUE) +
  scale_size_continuous(range = c(3, 10), breaks = c(1.5, 2, 3, 4), name = "qval") +
  scale_color_identity() +
  theme_bw() +
  labs(x = "Cluster",
       y = "Pathway",
       color = "Class Color") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "right") +
  guides(color = guide_legend(title = "Class Color", override.aes = list(size = 5)),
         size = guide_legend(title = "qval")) 
dev.off()






# GSEA plot
library("fgsea")


#### import all clsuter DEGs output :
cluster_types <- c( "Granule",
  "MLI1",
  "MLI2" ,
  "PLI12" ,
  "PLI23" ,
  "Golgi" ,
  "Unipolar_Brush" ,
  "Purkinje",
  "Astrocyte" ,
  "Bergman_Glia",
  "OPC",
  "Meningeal",
  "Endothelial",
  "Endothelial_Stalk",
  "Choroid_Plexus",
  "Unknown_Neuron_Subpop")


# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p35_CB_QCV2dim50kparam20res03_allGenes_correct.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}

## load list of genes to test
### List1
fgsea_sets <- list(
  PathwaysOfNeurodegeneration = read_table(file = "output/Pathway/geneList_mmu05022.txt")$Genes,
  GOBP_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS = read_table(file = "output/Pathway/geneList_GOBP_NEURON_INTRINSIC_APOPTOTIC.txt")$Genes,
  GOBP_NEUROINFLAMMATORY_RESPONSE = read_table(file = "output/Pathway/geneList_GOBP_NEUROINFLAMMATORY_RESPONSE.txt")$Genes
)
### List2
fgsea_sets <- list(
  GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION = read_table(file = "output/Pathway/geneList_GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION.txt")$Genes,
  GOBP_NEURONAL_ACTION_POTENTIAL = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_ACTION_POTENTIAL.txt")$Genes,
  GOBP_NEURONAL_SIGNAL_TRANSDUCTION = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_SIGNAL_TRANSDUCTION.txt")$Genes,
  REACTOME_POTASSIUM_CHANNELS = read_table(file = "output/Pathway/geneList_REACTOME_POTASSIUM_CHANNELS.txt")$Genes,
  REACTOME_ION_CHANNEL_TRANSPORT = read_table(file = "output/Pathway/geneList_REACTOME_ION_CHANNEL_TRANSPORT.txt")$Genes
)
### List3
fgsea_sets <- list(
  REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING = read_table(file = "output/Pathway/geneList_REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING.txt")$Genes,
  REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS = read_table(file = "output/Pathway/geneList_REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS.txt")$Genes,
  REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION = read_table(file = "output/Pathway/geneList_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION.txt")$Genes,
  REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS = read_table(file = "output/Pathway/geneList_REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS.txt")$Genes
)
### List4
fgsea_sets <- list(
  PathwaysOfNeurodegeneration = read_table(file = "output/Pathway/geneList_mmu05022.txt")$Genes,
  GOBP_NEURONAL_ACTION_POTENTIAL = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_ACTION_POTENTIAL.txt")$Genes,
  REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION = read_table(file = "output/Pathway/geneList_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION.txt")$Genes
)


## Rank genes based on FC
genes <- Purkinje %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 10000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(ES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -NES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()


## plot GSEA
pdf("output/Pathway/GSEA_Kcnc1_response_p35_CB_QCV3dim50kparam50res03_allGenes-mmu05022-Purkinje.pdf", width=5, height=3)

plotEnrichment(fgsea_sets[["PathwaysOfNeurodegeneration"]],
               ranks) + labs(title="PathwaysOfNeurodegeneration-Purkinje") +
               theme_bw()
dev.off()


# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c( "Granule",
  "MLI1",
  "MLI2" ,
  "PLI12" ,
  "PLI23" ,
  "Golgi" ,
  "Unipolar_Brush" ,
  "Purkinje",
  "Astrocyte" ,
  "Bergman_Glia",
  "OPC",
  "Meningeal",
  "Endothelial",
  "Endothelial_Stalk",
  "Choroid_Plexus",
  "Unknown_Neuron_Subpop")

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
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 10000)
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

write.table(final_results, file = c("output/Pathway/gsea_output_Kcnc1_response_p35_CB_QCV2dim50kparam20res03_allGenes_correct-List4.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap all GSEA
pdf("output/Pathway/heatmap_gsea_padj-Kcnc1_response_p35_CB_QCV2dim50kparam20res03_allGenes_correct-List4.pdf", width=10, height=3)
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

### if need reorder terms
final_results$pathway <- factor(final_results$pathway, levels = c("PathwaysOfNeurodegeneration", "GOBP_NEURONAL_ACTION_POTENTIAL", "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION")) 

pdf("output/Pathway/heatmap_gsea_padj-Kcnc1_response_p35_CB_QCV2dim50kparam20res03_allGenes_correct-List4.pdf", width=10, height=3)
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






##########################################
## integration WT Kcnc1 p180 all replicates (1st replicate, then genotype) ######
 ##########################################


WT_p180_CB_Rep1$replicate <- "Rep1"
WT_p180_CB_Rep2$replicate <- "Rep2"
WT_p180_CB_Rep3$replicate <- "Rep3"

WT_p180_CB_Rep1$condition <- "WT"
WT_p180_CB_Rep2$condition <- "WT"
WT_p180_CB_Rep3$condition <- "WT"

Kcnc1_p180_CB_Rep1$replicate <- "Rep1"
Kcnc1_p180_CB_Rep2$replicate <- "Rep2"
Kcnc1_p180_CB_Rep3$replicate <- "Rep3"

Kcnc1_p180_CB_Rep1$condition <- "Kcnc1"
Kcnc1_p180_CB_Rep2$condition <- "Kcnc1"
Kcnc1_p180_CB_Rep3$condition <- "Kcnc1"

set.seed(42)


# Replicate and Genotype integration (1 step integration)
## WT Rep

### Reg v1 _ better than Regv2
WT_p180_CB_Rep1 <- SCTransform(WT_p180_CB_Rep1, method = "glmGamPoi", ncells = 10375, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p180_CB_Rep2 <- SCTransform(WT_p180_CB_Rep2, method = "glmGamPoi", ncells = 10468, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p180_CB_Rep3 <- SCTransform(WT_p180_CB_Rep3, method = "glmGamPoi", ncells = 11641, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep1 <- SCTransform(Kcnc1_p180_CB_Rep1, method = "glmGamPoi", ncells = 10864, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep2 <- SCTransform(Kcnc1_p180_CB_Rep2, method = "glmGamPoi", ncells = 11931, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep3 <- SCTransform(Kcnc1_p180_CB_Rep3, method = "glmGamPoi", ncells = 12882, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 




srat.list <- list(WT_p180_CB_Rep1 = WT_p180_CB_Rep1, WT_p180_CB_Rep2 = WT_p180_CB_Rep2, WT_p180_CB_Rep3 = WT_p180_CB_Rep3, Kcnc1_p180_CB_Rep1 = Kcnc1_p180_CB_Rep1, Kcnc1_p180_CB_Rep2 = Kcnc1_p180_CB_Rep2, Kcnc1_p180_CB_Rep3 = Kcnc1_p180_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_Kcnc1_p180_CB_1step.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_p180_CB_1step.sct <- IntegrateData(anchorset = WT_Kcnc1_p180_CB_1step.anchors, normalization.method = "SCT")


#### UMAP
DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "integrated"

WT_Kcnc1_p180_CB_1step.sct <- RunPCA(WT_Kcnc1_p180_CB_1step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_p180_CB_1step.sct <- RunUMAP(WT_Kcnc1_p180_CB_1step.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_Kcnc1_p180_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p180_CB_1step.sct, reduction = "pca", k.param = 20, dims = 1:50)
WT_Kcnc1_p180_CB_1step.sct <- FindClusters(WT_Kcnc1_p180_CB_1step.sct, resolution = 0.2, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p180_CB_1step.sct$condition <- factor(WT_Kcnc1_p180_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50kparam20res02.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()

# genes

DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "SCT"

pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50-List6_7.pdf", width=30, height=60)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = c("Gabra6", "Pax6", "Kcnc2", "Sorcs3", "Ptprk", "Nxph1", "Cdh22", "Aldh1a3", "Pax2", "Eomes", "Rgs6", "Tafa2", "Calb1", "Slc1a6", "Car8", "Zeb2", "Aqp4", "Slc39a12", "Vcan", "Sox6", "Mbp", "Mag", "Plp1", "Aldoc", "Cnp", "Itgam", "Cx3cr1", "Ptgds", "Dcn", "Lef1", "Notum", "Apcdd1", "Dlc1", "Pdgfrb", "Kl",  "Ttr", "Actb", "Tmsb4x"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50-Kcnc1.pdf", width=10, height=10)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = c("Kcnc1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()




List6:
Granular = Gabra6, Pax6
Interneuron = Kcnc2
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
PLI = Aldh1a3
Golgi = Pax2
Unipolar_Brush = Eomes
Purkinje = Calb1, Slc1a6, Car8
Astrocyte = Zeb2
Bergman_Glia = Aqp4, Slc39a12
Oligodendrocyte= Mbp, Mag, Plp1
OPC = Aldoc, Cnp
Microglia = Itgam, Cx3cr1
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Endothelial_Mural =  Dlc1, Pdgfrb
Choroid plexus cells = Kl,  Ttr

List7:
Granular = Gabra6, Pax6
Interneuron = Kcnc2
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
Golgi = Pax2
Unipolar_Brush = Eomes, Rgs6, Tafa2
Purkinje = Calb1, Slc1a6, Car8
Astrocyte = Zeb2
Bergman_Glia = Aqp4, Slc39a12
OPC = Vcan, Sox6
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Choroid plexus cells = Kl,  Ttr
Endothelial_Stalk = Actb, Tmsb4x




#### QC metrics investigation ####################################
pdf("output/seurat/VlnPlot_QCmetrics_SCT_WT_Kcnc1_p180_CB_1step-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02-countMtRbRegression.pdf", width=20, height=5)
VlnPlot(WT_Kcnc1_p180_CB_1step.sct,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nFeature_RNA_SCT_WT_Kcnc1_p180_CB_1step-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p180_CB_1step.sct,features = c("nFeature_RNA")) +
  ylim(0,2000)
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nCount_RNA_SCT_WT_Kcnc1_p180_CB_1step-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p180_CB_1step.sct,features = c("nCount_RNA")) +
  ylim(0,10000)
dev.off()


pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p180_CB_nFeature_RNA-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p180_CB_percentmt-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p180_CB_percentrb-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.rb")
dev.off()  
############################################################



pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB_splitCondition-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=13, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "condition")
dev.off()
pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB_splitReplicate-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "replicate")
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_WT_p180_CB_Kcnc1_Phase-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  


# save ######################################################
## saveRDS(WT_Kcnc1_p180_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p180_CB_1step.sct_V1_numeric.rds") # regMtRbCount with QC_V3
## saveRDS(WT_Kcnc1_p180_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p180_CB_1step.sct_V2_numeric.rds") # regMtRbCount with QC_V4
## saveRDS(WT_Kcnc1_p180_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V2_numeric.rds") # regMtRbCount with QC_V4
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V2_numeric.rds") # QC_V4

## saveRDS(WT_Kcnc1_p180_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V1_label.rds") # regMtRbCount with QC_V4
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V1_label.rds") # QC_V4 with PLI12 PLI23
########################################################################



set.seed(42)
##########


## Let's work with 1step integration: 
# --> 1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50kparam20res02 = WT_Kcnc1_p180_CB_1step_V1
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V2_numeric.rds")


############################ EasyCellType automatic annotation ##########################################

### Find all markers 
all_markers <- FindAllMarkers(WT_Kcnc1_p180_CB_1step.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(all_markers, file = "output/seurat/srat_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_all_markers.txt", header = TRUE, row.names = 1)



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


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain", "Cerebellum", "Hippocampus"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?


annot.GSEA <- easyct(input.d, db="clustermole", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?





## plots



#pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p180_CB_1step_V1-cellmarker_CerebellumBrainHippocampus.pdf", width=6, height=8)
#pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p180_CB_1step_V1-clustermole_Brain.pdf", width=6, height=8)

pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02-clustermole_Brain.pdf", width=6, height=8)
plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




# panglao --> NOT working... 'subscript out of bounds' error... I tried gene as ENSEMBL, entrezID and geneSymbo, human/mic, everything...


#




## check some genes



DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "SCT"


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50-SuppTab1Marker.pdf", width=30, height=40)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = c( "Ppp1r17", "Gabra6", "Grm2", "Sst","Prkcd", "Sorcs3", "Ptprk", "Nxph1", "Cdh22","Htr2a", "Edil3","Aldh1a3", "Slc6a5","Eomes","Gdf10"
), max.cutoff = 1, cols = c("grey", "red"))
dev.off()






############ V1 naming (output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V2_numeric.rds = QCV4dim50kparam20res02)

Cluster1: Granular_1 (Gabra6, Pax6)
Cluster2: Granular_2 (Gabra6, Pax6)
Cluster3: MLI1 (Sorcs3, Ptprk)
Cluster4: Granular_3 (Gabra6, Pax6)
Cluster5: MLI2_1 (Nxph1, Cdh22)
Cluster6: PLI23 (Galntl6, Kcnc2)
Cluster7: Astrocyte (Zeb2)
Cluster8: PLI12 (Klhl1, Gfra2, Aldh1a3)
Cluster9: Bergman_Glia (Aqp4, Slc39a12)
Cluster10: Unipolar_Brush (Eomes, Rgs6, Tafa2)
Cluster11: Mix_Endothelial_EndothelialMural (Lef1, Notum, Apcdd1, Dlc1, Pdgfrb)
Cluster12: Meningeal (Ptgds, Dcn)
Cluster13: Choroid_Plexus (Kl, Ttr)
Cluster14: Golgi (Pax2)
Cluster15: Purkinje (Calb1, Slc1a6, Car8)
Cluster16: Unknown_Neuron_Subpop
Cluster17: Oligodendrocyte (Mbp, Mag, Plp1)



new.cluster.ids <- c(
  "Granular_1",
  "Granular_2",
  "MLI1",
  "Granular_3",
  "MLI2",
  "PLI23",
  "Astrocyte",
  "PLI12",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)

names(new.cluster.ids) <- levels(WT_Kcnc1_p180_CB_1step.sct)
WT_Kcnc1_p180_CB_1step.sct <- RenameIdents(WT_Kcnc1_p180_CB_1step.sct, new.cluster.ids)
WT_Kcnc1_p180_CB_1step.sct$cluster.annot <- Idents(WT_Kcnc1_p180_CB_1step.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_label.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_noSplit_label.pdf", width=9, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5)
dev.off()



# All in dotplot
DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "SCT"


List8:
Cluster1: Granular_1 (Gabra6, Pax6)
Cluster3: MLI1 (Sorcs3, Ptprk)
Cluster5: MLI2 (Nxph1, Cdh22)
Cluster6: PLI23 (Galntl6, Kcnc2)
Cluster7: Astrocyte (Zeb2)
Cluster9: Bergman_Glia (Aqp4, Slc39a12)
Cluster10: Unipolar_Brush (Eomes, Rgs6, Tafa2)
Cluster11: Mix_Endothelial_EndothelialMural (Lef1, Notum, Apcdd1, Dlc1, Pdgfrb)
Cluster12: Meningeal (Ptgds, Dcn)
Cluster13: Choroid_Plexus (Kl, Ttr)
Cluster14: Golgi (Pax2)
Cluster15: Purkinje (Calb1, Slc1a6, Car8)
Cluster16: Unknown_Neuron_Subpop (Gm42397,Hcrtr2)
Cluster17: Oligodendrocyte (Mbp, Mag, Plp1)



all_markers <- c(
  "Gabra6", "Pax6",
  "Sorcs3", "Ptprk",
  "Nxph1", "Cdh22",
  "Klhl1", "Gfra2", "Aldh1a3",
  "Galntl6", "Kcnc2",
  "Zeb2",
  "Aqp4", "Slc39a12",
  "Eomes", "Rgs6", "Tafa2",
  "Lef1", "Notum", "Apcdd1", "Dlc1", "Pdgfrb",
  "Ptgds", "Dcn",
  "Kl", "Ttr",
  "Pax2",
  "Calb1", "Slc1a6", "Car8",
  "Gm42397","Hcrtr2",
  "Mbp", "Mag", "Plp1"
)



levels(WT_Kcnc1_p180_CB_1step.sct) <- c(
  "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)



pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_label.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p180_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_label_vertical.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p180_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()




pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02-List8.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = all_markers, max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02-Kcnc1.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = "Kcnc1", max.cutoff = 1, cols = c("grey", "red"))
dev.off()






## p180 cell type proportion ###############################
### count nb of cells in each cluster
WT_p180_CB_Rep1 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "WT_p180_CB_Rep1"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p180",
             replicate= "Rep1")
WT_p180_CB_Rep2 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "WT_p180_CB_Rep2"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p180",
             replicate= "Rep2")
WT_p180_CB_Rep3 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "WT_p180_CB_Rep3"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p180",
             replicate= "Rep3")

Kcnc1_p180_CB_Rep1 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "Kcnc1_p180_CB_Rep1"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p180",
             replicate= "Rep1")
Kcnc1_p180_CB_Rep2 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "Kcnc1_p180_CB_Rep2"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p180",
             replicate= "Rep2")
Kcnc1_p180_CB_Rep3 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "Kcnc1_p180_CB_Rep3"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p180",
             replicate= "Rep3")


p180_CB = WT_p180_CB_Rep1 %>%
  bind_rows(WT_p180_CB_Rep2) %>%
  bind_rows(WT_p180_CB_Rep3) %>%
  bind_rows(Kcnc1_p180_CB_Rep1) %>%
  bind_rows(Kcnc1_p180_CB_Rep2) %>%
  bind_rows(Kcnc1_p180_CB_Rep3) %>%
  as_tibble()
  



### Keeping all replicates
p180_CB_prop = p180_CB %>%
  group_by(replicate, genotype) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = (count / total_count) * 100)

p180_CB_prop$genotype <-
  factor(p180_CB_prop$genotype,
         c("WT", "Kcnc1"))

pdf("output/seurat/histogramProp_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02.pdf", width=7, height=4)
ggbarplot(p180_CB_prop, x = "cluster", y = "proportion", fill = "genotype",
                  color = "genotype", palette = c("black", "blue"),
                  position = position_dodge(0.8), # Separate bars by genotype
                  add = "mean_se", # Add error bars
                  lab.pos = "out", lab.size = 3) +
  stat_compare_means(aes(group = genotype), method = "t.test", label = "p.signif") +
  theme_bw() +
  labs(x = "Cell Type (Cluster)", y = "Cell Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
# --> Too long run in slurm job



# DEGs number colored in a UMAP
Idents(WT_Kcnc1_p180_CB_1step.sct) <- "cluster.annot"

DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())
## List of cell types
cell_types <- c(      "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte")
## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/seurat/", cell_type, "-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct.txt", sep = "")
  deg_data <- read.table(file_name, header = TRUE, sep = "\t") ## Read the DEGs data
  num_degs <- sum(deg_data$p_val_adj < 0.05 & (deg_data$avg_log2FC > 0.25 | deg_data$avg_log2FC < -0.25)) ## Count the number of significant DEGs
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))  ## Append to the summary table
}

DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = c(  "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte")) 
  
  
# Add DEG information to my seurat object - DEG_count
cell_clusters <- WT_Kcnc1_p180_CB_1step.sct@meta.data$cluster.annot
names(cell_clusters) <- rownames(WT_Kcnc1_p180_CB_1step.sct@meta.data)
DEG_named_vector <- DEG_count$Num_DEGs[match(cell_clusters, DEG_count$Cell_Type)]
names(DEG_named_vector) <- names(cell_clusters)
# Integrate DEG values into the Seurat object
WT_Kcnc1_p180_CB_1step.sct <- AddMetaData(WT_Kcnc1_p180_CB_1step.sct, metadata = DEG_named_vector, col.name = "DEG")
# Create a UMAP plot colored by qval values
pdf("output/seurat/FeaturePlot_WT_Kcnc1_p180_CB_1step_DEG.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = "DEG", pt.size = 0.5, reduction = "umap", label = TRUE) +
  scale_colour_viridis() #  option="magma"
dev.off()
# Add values on the heatmap
## Extract UMAP coordinates
umap_coordinates <- as.data.frame(WT_Kcnc1_p180_CB_1step.sct@reductions$umap@cell.embeddings)
umap_coordinates$cluster <- WT_Kcnc1_p180_CB_1step.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DEG_count %>% dplyr::rename( "cluster"="Cell_Type"))
## Create a UMAP plot colored by DEG values, with cluster DEG counts as text annotations
pdf("output/seurat/FeaturePlot_WT_Kcnc1_p180_CB_1step_DEGpadj05fc025_numeric.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = "DEG", pt.size = 0.5, reduction = "umap") +
  scale_colour_viridis(option="mako") + # option="magma"
  geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = Num_DEGs), 
            size = 4, color = "red", fontface = "bold") 
dev.off()




# SCPA

# Compare WT and cYAPKO using SCPA ##########################################

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
library("SCPA")
library("circlize")
library("magrittr")
library("msigdb")
library("msigdbr")
library("ComplexHeatmap")
library("ggrepel")
library("ggpubr")

DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "RNA" # Recommended 



# Test different Pathway collections and generate enrichment plot for each cell types (C2 = Pathway, C5 = ontology )
## import Pathways
pathways <- msigdbr("Mus musculus", "C5") %>%          # !!!!!! CHANGE HERE PATHWAYS !!!!!!
format_pathways()
names(pathways) <- sapply(pathways, function(x) x$Pathway[1]) # just to name the list, so easier to visualise

# Code to save output for each cell type comparison
clusters = c(
  "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2_1",
  "MLI2_2",
  "Interneuron",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)
### Loop through each value
for (cluster in clusters) {
  #### Extract data for WT and cYAPKO based on current value
  WT <- seurat_extract(WT_Kcnc1_p180_CB_1step.sct,
                       meta1 = "condition", value_meta1 = "WT",
                       meta2 = "cluster.annot", value_meta2 = cluster)

  Kcnc1 <- seurat_extract(WT_Kcnc1_p180_CB_1step.sct,
                           meta1 = "condition", value_meta1 = "Kcnc1",
                           meta2 = "cluster.annot", value_meta2 = cluster)

  ##### Compare pathways
  WT_cYAPKO <- compare_pathways(samples = list(WT, Kcnc1),
                                pathways = pathways,
                                parallel = TRUE, cores = 8)

  ##### Write to file using the current value in the filename
  output_filename <- paste0("output/Pathway/SCPA_CB_p180_C5_", cluster, ".txt")       # !!!!!! CHANGE HERE PATHWAYS !!!!!!
  write.table(WT_cYAPKO, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)
}
#--> Long ~2hrs



## load all the comparison for each cell type (FC qval information)
clusters = c(
  "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2_1",
  "MLI2_2",
  "Interneuron",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)
## import with a function
### A function to read and add the cluster column
read_and_add_cluster <- function(cluster) {
  path <- paste0("output/Pathway/SCPA_CB_p180_C2_", cluster, ".txt")
  df <- read.delim(path, header = TRUE) %>%
    add_column(cluster = cluster)
  return(df)
}
### Use lapply to apply the function on each cluster and bind all data frames together
all_data <- bind_rows(lapply(clusters, read_and_add_cluster)) %>% as_tibble()

## Filter pathway of interest
pathways <- c(
  "WP_NEURODEGENERATION_WITH_BRAIN_IRON_ACCUMULATION_NBIA_SUBTYPES_PATHWAY",
  "WP_NEUROINFLAMMATION_AND_GLUTAMATERGIC_SIGNALING",
  "KEGG_ALZHEIMERS_DISEASE",
  "WP_ALZHEIMERS_DISEASE",
  "KEGG_PARKINSONS_DISEASE",
  "WP_PARKINSONS_DISEASE_PATHWAY",
  "KEGG_HUNTINGTONS_DISEASE",
  "WP_ERK_PATHWAY_IN_HUNTINGTONS_DISEASE",
  "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
  "WP_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
  "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
  "REACTOME_NEUROTRANSMITTER_RELEASE_CYCLE",
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS",
  "REACTOME_POTASSIUM_CHANNELS",
  "REACTOME_ION_CHANNEL_TRANSPORT",
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  "REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS",
  "REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_DOPAMINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_GLUTAMATE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_NOREPINEPHRINE_NEUROTRANSMITTER_RELEASE_CYCLE",
  "REACTOME_SEROTONIN_NEUROTRANSMITTER_RELEASE_CYCLE",
  "ALCALA_APOPTOSIS",
  "KEGG_APOPTOSIS",
  "REACTOME_APOPTOSIS",
  "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
  "WP_APOPTOSIS",
  "WP_APOPTOSIS_MODULATION_AND_SIGNALING",
  "REACTOME_DEATH_RECEPTOR_SIGNALLING",
  "REACTOME_DISEASES_OF_PROGRAMMED_CELL_DEATH",
  "REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES",
  "REACTOME_PROGRAMMED_CELL_DEATH",
  "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
  "BIOCARTA_DEATH_PATHWAY",
  "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES",
  "WP_PKCGAMMA_CALCIUM_SIGNALING_PATHWAY_IN_ATAXIA"
)
classes <- c(
  rep("Neurodegeneration", 10),
  rep("Neuronal_Activity", 13),
  rep("Apoptosis", 12),
  "ROS",
  "Ataxia"
)
colors <- c(
  rep("Purple", 10),
  rep("Green", 13),
  rep("Dark Orange", 12),
  "Red",
  "Blue"
)
pathway_tibble <- tibble(
  Pathway = pathways,
  Class = classes,
  Color = colors
)

pathway_all_data = pathway_tibble %>%
  left_join(all_data)

pathway_all_data$Pathway <- factor(pathway_all_data$Pathway, levels = pathways)


### dotplot
pdf("output/Pathway/dotplot_SCPA_CB_p180_C2_selectedPathwayV1.pdf", width=12, height=8)
ggplot(pathway_all_data, aes(x = cluster, y = Pathway)) +
  geom_point(aes(size = ifelse(qval > 1.4, qval, NA), color = Color), na.rm = TRUE) +
  scale_size_continuous(range = c(3, 10), breaks = c(1.5, 2, 3, 4), name = "qval") +
  scale_color_identity() +
  theme_bw() +
  labs(x = "Cluster",
       y = "Pathway",
       color = "Class Color") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "right") +
  guides(color = guide_legend(title = "Class Color", override.aes = list(size = 5)),
         size = guide_legend(title = "qval"))
dev.off()

# --> NO SIGNIFICANT TERMS!





# GSEA plot
library("fgsea")


#### import all clsuter DEGs output :
cluster_types <- c(     "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}

## load list of genes to test
### List1
fgsea_sets <- list(
  PathwaysOfNeurodegeneration = read_table(file = "output/Pathway/geneList_mmu05022.txt")$Genes,
  GOBP_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS = read_table(file = "output/Pathway/geneList_GOBP_NEURON_INTRINSIC_APOPTOTIC.txt")$Genes,
  GOBP_NEUROINFLAMMATORY_RESPONSE = read_table(file = "output/Pathway/geneList_GOBP_NEUROINFLAMMATORY_RESPONSE.txt")$Genes
)
### List2
fgsea_sets <- list(
  GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION = read_table(file = "output/Pathway/geneList_GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION.txt")$Genes,
  GOBP_NEURONAL_ACTION_POTENTIAL = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_ACTION_POTENTIAL.txt")$Genes,
  GOBP_NEURONAL_SIGNAL_TRANSDUCTION = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_SIGNAL_TRANSDUCTION.txt")$Genes,
  REACTOME_POTASSIUM_CHANNELS = read_table(file = "output/Pathway/geneList_REACTOME_POTASSIUM_CHANNELS.txt")$Genes,
  REACTOME_ION_CHANNEL_TRANSPORT = read_table(file = "output/Pathway/geneList_REACTOME_ION_CHANNEL_TRANSPORT.txt")$Genes
)
### List3
fgsea_sets <- list(
  REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING = read_table(file = "output/Pathway/geneList_REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING.txt")$Genes,
  REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS = read_table(file = "output/Pathway/geneList_REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS.txt")$Genes,
  REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION = read_table(file = "output/Pathway/geneList_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION.txt")$Genes,
  REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS = read_table(file = "output/Pathway/geneList_REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS.txt")$Genes
)
### List4
fgsea_sets <- list(
  PathwaysOfNeurodegeneration = read_table(file = "output/Pathway/geneList_mmu05022.txt")$Genes,
  GOBP_NEURONAL_ACTION_POTENTIAL = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_ACTION_POTENTIAL.txt")$Genes,
  REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION = read_table(file = "output/Pathway/geneList_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION.txt")$Genes
)


## Rank genes based on FC
genes <- Purkinje %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(ES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -NES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()


## plot GSEA
pdf("output/Pathway/GSEA_Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes-mmu05022-Purkinje.pdf", width=5, height=3)

plotEnrichment(fgsea_sets[["PathwaysOfNeurodegeneration"]],
               ranks) + labs(title="PathwaysOfNeurodegeneration-Purkinje") +
               theme_bw()
dev.off()


# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c(    "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte")

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
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 10000)
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

write.table(final_results, file = c("output/Pathway/gsea_output_Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct-List4.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap all GSEA
pdf("output/Pathway/heatmap_gsea_padj-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct-List4.pdf", width=15, height=3)
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


### if need reorder terms
final_results$pathway <- factor(final_results$pathway, levels = c("PathwaysOfNeurodegeneration", "GOBP_NEURONAL_ACTION_POTENTIAL", "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION")) 

pdf("output/Pathway/heatmap_gsea_padj-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct-List4.pdf", width=10, height=3)
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

--> How many dims to use? Not clear, using Elbow and [quantification](https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html) say 13dims (but look very few!!). We do not observe significant changes of clustering by changes the numb of dims. Only small cluster are affected (Serotonergic neurons)

--> p35 CB sample is problematic. Most genes in most cluster are downregulated, hinting at a technical issue... From all samples, Kcnc1_Rep2 seems oputlier with a lot of cells and very low nCount_RNA. I tried removing this sample but same issue occur. I also tried using [Harmony](https://github.com/immunogenomics/harmony) instead SCT-seurat pipleline for data sample integration without success, same bias is occuring. I also testing working on Rep1 of WT and Kcnc1 and XXX
  - [Discussion on Harmony pipeline to integrate with SCTransform](https://github.com/satijalab/seurat/issues/4896): Pipeline I used
  - [Discussion about Harmony and seurat](https://github.com/immunogenomics/harmony/issues/41)


    --> Could try [MNN integration](https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.md) method? 

--> Let's try alternative DEG test; discuss [here](https://satijalab.org/seurat/articles/de_vignette)

--> **FALSE ALERT!!! The issue was that I forgetted to normalize and scale data prior doing DEG....**



## General trend Cortex (CX)

I copied beginning and p180 part from `## General trend Cerebellum (CB)`


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
  WT_p14_CX_Rep1 = "output/soupX/WT_p14_CX_Rep1.RData",
  WT_p14_CX_Rep2 = "output/soupX/WT_p14_CX_Rep2.RData",
  WT_p14_CX_Rep3 = "output/soupX/WT_p14_CX_Rep3.RData",
  Kcnc1_p14_CX_Rep1 = "output/soupX/Kcnc1_p14_CX_Rep1.RData",
  Kcnc1_p14_CX_Rep2 = "output/soupX/Kcnc1_p14_CX_Rep2.RData",
  Kcnc1_p14_CX_Rep3 = "output/soupX/Kcnc1_p14_CX_Rep3.RData",

  WT_p35_CX_Rep1 = "output/soupX/WT_p35_CX_Rep1.RData",
  WT_p35_CX_Rep2 = "output/soupX/WT_p35_CX_Rep2.RData",
  WT_p35_CX_Rep3 = "output/soupX/WT_p35_CX_Rep3.RData",
  Kcnc1_p35_CX_Rep1 = "output/soupX/Kcnc1_p35_CX_Rep1.RData",
  Kcnc1_p35_CX_Rep2 = "output/soupX/Kcnc1_p35_CX_Rep2.RData",
  Kcnc1_p35_CX_Rep3 = "output/soupX/Kcnc1_p35_CX_Rep3.RData",

  WT_p180_CX_Rep1 = "output/soupX/WT_p180_CX_Rep1.RData",
  WT_p180_CX_Rep2 = "output/soupX/WT_p180_CX_Rep2.RData",
  WT_p180_CX_Rep3 = "output/soupX/WT_p180_CX_Rep3.RData",
  Kcnc1_p180_CX_Rep1 = "output/soupX/Kcnc1_p180_CX_Rep1.RData",
  Kcnc1_p180_CX_Rep2 = "output/soupX/Kcnc1_p180_CX_Rep2.RData",
  Kcnc1_p180_CX_Rep3 = "output/soupX/Kcnc1_p180_CX_Rep3.RData"
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


##### QC filtering _ V1 ############################################
#--> Not enough stringeant
##### QC filtering _ V2 (used for CB p35) ############################################
### mit > 15, rb > 22 and RNAfeaure 100 #############################################
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 100 & seurat_object@meta.data$QC == 'Pass', 'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 100 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 15 & seurat_object@meta.data$QC == 'Pass', 'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 15 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 22 & seurat_object@meta.data$QC == 'Pass', 'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 22 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list
#####################################################################################################



##### QC filtering _ V3 (used for CB p14) ############################################

###  mit > 15, rb > 10 and RNAfeaure 400 #############################################
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC == 'Pass', 'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 15 & seurat_object@meta.data$QC == 'Pass', 'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 15 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 10 & seurat_object@meta.data$QC == 'Pass', 'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 10 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list
#####################################################################################################



##### QC filtering _ V4 (used for CB  p180) ############################################

###  mit > 7, rb > 5 and RNAfeaure 400 #############################################
apply_qc <- function(seurat_object) {
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$Is_doublet == 'True', 'Doublet', 'Pass')
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC == 'Pass', 'Low_nFeature', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$nFeature_RNA < 400 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'Low_nFeature', paste('Low_nFeature', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 7 & seurat_object@meta.data$QC == 'Pass', 'High_MT', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.mt > 7 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_MT', paste('High_MT', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 5 & seurat_object@meta.data$QC == 'Pass', 'High_RB', seurat_object@meta.data$QC)
  seurat_object[['QC']] <- ifelse(seurat_object@meta.data$percent.rb > 5 & seurat_object@meta.data$QC != 'Pass' & seurat_object@meta.data$QC != 'High_RB', paste('High_RB', seurat_object@meta.data$QC, sep = ','), seurat_object@meta.data$QC)
  return(seurat_object)
}
for (sample_name in names(seurat_objects)) {
  seurat_objects[[sample_name]] <- apply_qc(seurat_objects[[sample_name]])
}
assign_seurat_objects(seurat_objects) # This NEED to be reapply to apply the previous function to all individual in our list
#####################################################################################################


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
write.table(qc_summary_combined, file = "output/seurat/QC_summary_V2_CX.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # 


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
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V2_CX.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) # 



############ SAVE samples #################################################################
save_seurat_objects <- function(seurat_objects_list, output_dir) {
  for (sample_name in names(seurat_objects_list)) {
    file_name <- paste0(output_dir, sample_name, "_V2_numeric.rds") # CB: _V2_numeric \\\--> V3 for p14 ; v4 for p180; V2 for p35
    saveRDS(seurat_objects_list[[sample_name]], file = file_name) # V2 got a ReProcess for sample p35
  }
}
output_dir <- "output/seurat/"
## Call the function to save the Seurat objects
save_seurat_objects(seurat_objects, output_dir)
###############################################################################################
############# READ samples  (QC V1 V2 or V3 = V1_numeric...)################################################
# Function to load Seurat objects
load_seurat_objects <- function(file_paths) {
  seurat_objects <- list()
  for (file_path in file_paths) {
    sample_name <- gsub("_V2_numeric.rds", "", basename(file_path))
    seurat_objects[[sample_name]] <- readRDS(file_path)
  }
  return(seurat_objects)
}
output_dir <- "output/seurat/"
file_paths <- list.files(output_dir, pattern = "_V2_numeric.rds$", full.names = TRUE)
# Call the function to load the Seurat objects
seurat_objects <- load_seurat_objects(file_paths)
# Loop through the list and assign each Seurat object to a variable with the same name
for (sample_name in names(seurat_objects)) {
  assign(sample_name, seurat_objects[[sample_name]])
}
# 1 sample: WT_p14_CB_Rep2 <- readRDS(file = "output/seurat/WT_p14_CB_Rep2_V1_numeric.rds")
## Kcnc1_p14_CB_Rep1 <- readRDS(file = "output/seurat/Kcnc1_p14_CB_Rep1_V1_numeric.rds")
################################################################################################

XXXY HERE !!! LOAD V2





##########################################
## integration WT Kcnc1 p180 all replicates (1st replicate, then genotype) ######
 ##########################################


WT_p180_CB_Rep1$replicate <- "Rep1"
WT_p180_CB_Rep2$replicate <- "Rep2"
WT_p180_CB_Rep3$replicate <- "Rep3"

WT_p180_CB_Rep1$condition <- "WT"
WT_p180_CB_Rep2$condition <- "WT"
WT_p180_CB_Rep3$condition <- "WT"

Kcnc1_p180_CB_Rep1$replicate <- "Rep1"
Kcnc1_p180_CB_Rep2$replicate <- "Rep2"
Kcnc1_p180_CB_Rep3$replicate <- "Rep3"

Kcnc1_p180_CB_Rep1$condition <- "Kcnc1"
Kcnc1_p180_CB_Rep2$condition <- "Kcnc1"
Kcnc1_p180_CB_Rep3$condition <- "Kcnc1"

set.seed(42)


# Replicate and Genotype integration (1 step integration)
## WT Rep

### Reg v1 _ better than Regv2
WT_p180_CB_Rep1 <- SCTransform(WT_p180_CB_Rep1, method = "glmGamPoi", ncells = 10375, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p180_CB_Rep2 <- SCTransform(WT_p180_CB_Rep2, method = "glmGamPoi", ncells = 10468, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p180_CB_Rep3 <- SCTransform(WT_p180_CB_Rep3, method = "glmGamPoi", ncells = 11641, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep1 <- SCTransform(Kcnc1_p180_CB_Rep1, method = "glmGamPoi", ncells = 10864, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep2 <- SCTransform(Kcnc1_p180_CB_Rep2, method = "glmGamPoi", ncells = 11931, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p180_CB_Rep3 <- SCTransform(Kcnc1_p180_CB_Rep3, method = "glmGamPoi", ncells = 12882, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 




srat.list <- list(WT_p180_CB_Rep1 = WT_p180_CB_Rep1, WT_p180_CB_Rep2 = WT_p180_CB_Rep2, WT_p180_CB_Rep3 = WT_p180_CB_Rep3, Kcnc1_p180_CB_Rep1 = Kcnc1_p180_CB_Rep1, Kcnc1_p180_CB_Rep2 = Kcnc1_p180_CB_Rep2, Kcnc1_p180_CB_Rep3 = Kcnc1_p180_CB_Rep3)
features <- SelectIntegrationFeatures(object.list = srat.list, nfeatures = 3000)
srat.list <- PrepSCTIntegration(object.list = srat.list, anchor.features = features)
WT_Kcnc1_p180_CB_1step.anchors <- FindIntegrationAnchors(object.list = srat.list, normalization.method = "SCT",
    anchor.features = features)
WT_Kcnc1_p180_CB_1step.sct <- IntegrateData(anchorset = WT_Kcnc1_p180_CB_1step.anchors, normalization.method = "SCT")


#### UMAP
DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "integrated"

WT_Kcnc1_p180_CB_1step.sct <- RunPCA(WT_Kcnc1_p180_CB_1step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_p180_CB_1step.sct <- RunUMAP(WT_Kcnc1_p180_CB_1step.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_Kcnc1_p180_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p180_CB_1step.sct, reduction = "pca", k.param = 20, dims = 1:50)
WT_Kcnc1_p180_CB_1step.sct <- FindClusters(WT_Kcnc1_p180_CB_1step.sct, resolution = 0.2, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p180_CB_1step.sct$condition <- factor(WT_Kcnc1_p180_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50kparam20res02.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()

# genes

DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "SCT"

pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50-List6_7.pdf", width=30, height=60)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = c("Gabra6", "Pax6", "Kcnc2", "Sorcs3", "Ptprk", "Nxph1", "Cdh22", "Aldh1a3", "Pax2", "Eomes", "Rgs6", "Tafa2", "Calb1", "Slc1a6", "Car8", "Zeb2", "Aqp4", "Slc39a12", "Vcan", "Sox6", "Mbp", "Mag", "Plp1", "Aldoc", "Cnp", "Itgam", "Cx3cr1", "Ptgds", "Dcn", "Lef1", "Notum", "Apcdd1", "Dlc1", "Pdgfrb", "Kl",  "Ttr", "Actb", "Tmsb4x"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50-Kcnc1.pdf", width=10, height=10)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = c("Kcnc1"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()




List6:
Granular = Gabra6, Pax6
Interneuron = Kcnc2
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
PLI = Aldh1a3
Golgi = Pax2
Unipolar_Brush = Eomes
Purkinje = Calb1, Slc1a6, Car8
Astrocyte = Zeb2
Bergman_Glia = Aqp4, Slc39a12
Oligodendrocyte= Mbp, Mag, Plp1
OPC = Aldoc, Cnp
Microglia = Itgam, Cx3cr1
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Endothelial_Mural =  Dlc1, Pdgfrb
Choroid plexus cells = Kl,  Ttr

List7:
Granular = Gabra6, Pax6
Interneuron = Kcnc2
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
Golgi = Pax2
Unipolar_Brush = Eomes, Rgs6, Tafa2
Purkinje = Calb1, Slc1a6, Car8
Astrocyte = Zeb2
Bergman_Glia = Aqp4, Slc39a12
OPC = Vcan, Sox6
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Choroid plexus cells = Kl,  Ttr
Endothelial_Stalk = Actb, Tmsb4x




#### QC metrics investigation ####################################
pdf("output/seurat/VlnPlot_QCmetrics_SCT_WT_Kcnc1_p180_CB_1step-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02-countMtRbRegression.pdf", width=20, height=5)
VlnPlot(WT_Kcnc1_p180_CB_1step.sct,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nFeature_RNA_SCT_WT_Kcnc1_p180_CB_1step-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p180_CB_1step.sct,features = c("nFeature_RNA")) +
  ylim(0,2000)
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nCount_RNA_SCT_WT_Kcnc1_p180_CB_1step-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p180_CB_1step.sct,features = c("nCount_RNA")) +
  ylim(0,10000)
dev.off()


pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p180_CB_nFeature_RNA-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p180_CB_percentmt-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p180_CB_percentrb-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.rb")
dev.off()  
############################################################



pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB_splitCondition-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=13, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "condition")
dev.off()
pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB_splitReplicate-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "replicate")
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_WT_p180_CB_Kcnc1_Phase-1stepIntegrationRegressNotRepeated-QCV4dim50kparam20res02.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  


# save ######################################################
## saveRDS(WT_Kcnc1_p180_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p180_CB_1step.sct_V1_numeric.rds") # regMtRbCount with QC_V3
## saveRDS(WT_Kcnc1_p180_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p180_CB_1step.sct_V2_numeric.rds") # regMtRbCount with QC_V4
## saveRDS(WT_Kcnc1_p180_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V2_numeric.rds") # regMtRbCount with QC_V4
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V2_numeric.rds") # QC_V4

## saveRDS(WT_Kcnc1_p180_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V1_label.rds") # regMtRbCount with QC_V4
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V1_label.rds") # QC_V4 with PLI12 PLI23
########################################################################



set.seed(42)
##########


## Let's work with 1step integration: 
# --> 1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50kparam20res02 = WT_Kcnc1_p180_CB_1step_V1
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V2_numeric.rds")


############################ EasyCellType automatic annotation ##########################################

### Find all markers 
all_markers <- FindAllMarkers(WT_Kcnc1_p180_CB_1step.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(all_markers, file = "output/seurat/srat_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_all_markers.txt", header = TRUE, row.names = 1)



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


annot.GSEA <- easyct(input.d, db="cellmarker", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain", "Cerebellum", "Hippocampus"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?


annot.GSEA <- easyct(input.d, db="clustermole", # cellmarker or panglao or clustermole
                    species="Mouse", #  Human or Mouse
                    tissue=c("Brain"), p_cut=0.5,   # to see: data(cellmarker_tissue), data(clustermole_tissue), data(panglao_tissue)
                    test="GSEA")    # GSEA or fisher?





## plots



#pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p180_CB_1step_V1-cellmarker_CerebellumBrainHippocampus.pdf", width=6, height=8)
#pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p180_CB_1step_V1-clustermole_Brain.pdf", width=6, height=8)

pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02-clustermole_Brain.pdf", width=6, height=8)
plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




# panglao --> NOT working... 'subscript out of bounds' error... I tried gene as ENSEMBL, entrezID and geneSymbo, human/mic, everything...


#




## check some genes



DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "SCT"


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV4dim50-SuppTab1Marker.pdf", width=30, height=40)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = c( "Ppp1r17", "Gabra6", "Grm2", "Sst","Prkcd", "Sorcs3", "Ptprk", "Nxph1", "Cdh22","Htr2a", "Edil3","Aldh1a3", "Slc6a5","Eomes","Gdf10"
), max.cutoff = 1, cols = c("grey", "red"))
dev.off()






############ V1 naming (output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V2_numeric.rds = QCV4dim50kparam20res02)

Cluster1: Granular_1 (Gabra6, Pax6)
Cluster2: Granular_2 (Gabra6, Pax6)
Cluster3: MLI1 (Sorcs3, Ptprk)
Cluster4: Granular_3 (Gabra6, Pax6)
Cluster5: MLI2_1 (Nxph1, Cdh22)
Cluster6: PLI23 (Galntl6, Kcnc2)
Cluster7: Astrocyte (Zeb2)
Cluster8: PLI12 (Klhl1, Gfra2, Aldh1a3)
Cluster9: Bergman_Glia (Aqp4, Slc39a12)
Cluster10: Unipolar_Brush (Eomes, Rgs6, Tafa2)
Cluster11: Mix_Endothelial_EndothelialMural (Lef1, Notum, Apcdd1, Dlc1, Pdgfrb)
Cluster12: Meningeal (Ptgds, Dcn)
Cluster13: Choroid_Plexus (Kl, Ttr)
Cluster14: Golgi (Pax2)
Cluster15: Purkinje (Calb1, Slc1a6, Car8)
Cluster16: Unknown_Neuron_Subpop
Cluster17: Oligodendrocyte (Mbp, Mag, Plp1)



new.cluster.ids <- c(
  "Granular_1",
  "Granular_2",
  "MLI1",
  "Granular_3",
  "MLI2",
  "PLI23",
  "Astrocyte",
  "PLI12",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)

names(new.cluster.ids) <- levels(WT_Kcnc1_p180_CB_1step.sct)
WT_Kcnc1_p180_CB_1step.sct <- RenameIdents(WT_Kcnc1_p180_CB_1step.sct, new.cluster.ids)
WT_Kcnc1_p180_CB_1step.sct$cluster.annot <- Idents(WT_Kcnc1_p180_CB_1step.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_label.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_noSplit_label.pdf", width=9, height=6)
DimPlot(WT_Kcnc1_p180_CB_1step.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5)
dev.off()



# All in dotplot
DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "SCT"


List8:
Cluster1: Granular_1 (Gabra6, Pax6)
Cluster3: MLI1 (Sorcs3, Ptprk)
Cluster5: MLI2 (Nxph1, Cdh22)
Cluster6: PLI23 (Galntl6, Kcnc2)
Cluster7: Astrocyte (Zeb2)
Cluster9: Bergman_Glia (Aqp4, Slc39a12)
Cluster10: Unipolar_Brush (Eomes, Rgs6, Tafa2)
Cluster11: Mix_Endothelial_EndothelialMural (Lef1, Notum, Apcdd1, Dlc1, Pdgfrb)
Cluster12: Meningeal (Ptgds, Dcn)
Cluster13: Choroid_Plexus (Kl, Ttr)
Cluster14: Golgi (Pax2)
Cluster15: Purkinje (Calb1, Slc1a6, Car8)
Cluster16: Unknown_Neuron_Subpop (Gm42397,Hcrtr2)
Cluster17: Oligodendrocyte (Mbp, Mag, Plp1)



all_markers <- c(
  "Gabra6", "Pax6",
  "Sorcs3", "Ptprk",
  "Nxph1", "Cdh22",
  "Klhl1", "Gfra2", "Aldh1a3",
  "Galntl6", "Kcnc2",
  "Zeb2",
  "Aqp4", "Slc39a12",
  "Eomes", "Rgs6", "Tafa2",
  "Lef1", "Notum", "Apcdd1", "Dlc1", "Pdgfrb",
  "Ptgds", "Dcn",
  "Kl", "Ttr",
  "Pax2",
  "Calb1", "Slc1a6", "Car8",
  "Gm42397","Hcrtr2",
  "Mbp", "Mag", "Plp1"
)



levels(WT_Kcnc1_p180_CB_1step.sct) <- c(
  "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)



pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_label.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p180_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red")) + RotatedAxis()
dev.off()

pdf("output/seurat/DotPlot_SCT_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02_label_vertical.pdf", width=11, height=4.5)
DotPlot(WT_Kcnc1_p180_CB_1step.sct, assay = "SCT", features = all_markers, cols = c("grey", "red"))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
dev.off()




pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02-List8.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = all_markers, max.cutoff = 1, cols = c("grey", "red"))
dev.off()


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02-Kcnc1.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = "Kcnc1", max.cutoff = 1, cols = c("grey", "red"))
dev.off()






## p180 cell type proportion ###############################
### count nb of cells in each cluster
WT_p180_CB_Rep1 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "WT_p180_CB_Rep1"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p180",
             replicate= "Rep1")
WT_p180_CB_Rep2 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "WT_p180_CB_Rep2"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p180",
             replicate= "Rep2")
WT_p180_CB_Rep3 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "WT_p180_CB_Rep3"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "WT",
             time= "p180",
             replicate= "Rep3")

Kcnc1_p180_CB_Rep1 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "Kcnc1_p180_CB_Rep1"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p180",
             replicate= "Rep1")
Kcnc1_p180_CB_Rep2 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "Kcnc1_p180_CB_Rep2"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p180",
             replicate= "Rep2")
Kcnc1_p180_CB_Rep3 = table(Idents(WT_Kcnc1_p180_CB_1step.sct)[WT_Kcnc1_p180_CB_1step.sct$orig.ident == "Kcnc1_p180_CB_Rep3"]) %>%
  as.data.frame() %>%
  dplyr::rename("cluster"= "Var1" , "count" = "Freq") %>%
  add_column(genotype= "Kcnc1",
             time= "p180",
             replicate= "Rep3")


p180_CB = WT_p180_CB_Rep1 %>%
  bind_rows(WT_p180_CB_Rep2) %>%
  bind_rows(WT_p180_CB_Rep3) %>%
  bind_rows(Kcnc1_p180_CB_Rep1) %>%
  bind_rows(Kcnc1_p180_CB_Rep2) %>%
  bind_rows(Kcnc1_p180_CB_Rep3) %>%
  as_tibble()
  



### Keeping all replicates
p180_CB_prop = p180_CB %>%
  group_by(replicate, genotype) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = (count / total_count) * 100)

p180_CB_prop$genotype <-
  factor(p180_CB_prop$genotype,
         c("WT", "Kcnc1"))

pdf("output/seurat/histogramProp_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02.pdf", width=7, height=4)
ggbarplot(p180_CB_prop, x = "cluster", y = "proportion", fill = "genotype",
                  color = "genotype", palette = c("black", "blue"),
                  position = position_dodge(0.8), # Separate bars by genotype
                  add = "mean_se", # Add error bars
                  lab.pos = "out", lab.size = 3) +
  stat_compare_means(aes(group = genotype), method = "t.test", label = "p.signif") +
  theme_bw() +
  labs(x = "Cell Type (Cluster)", y = "Cell Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# differential expressed genes across conditions
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
# --> Too long run in slurm job



# DEGs number colored in a UMAP
Idents(WT_Kcnc1_p180_CB_1step.sct) <- "cluster.annot"

DEG_count <- data.frame(Cell_Type = character(), Num_DEGs = integer())
## List of cell types
cell_types <- c(      "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte")
## Loop through each cell type to count the number of significant DEGs
for (cell_type in cell_types) {
  file_name <- paste("output/seurat/", cell_type, "-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct.txt", sep = "")
  deg_data <- read.table(file_name, header = TRUE, sep = "\t") ## Read the DEGs data
  num_degs <- sum(deg_data$p_val_adj < 0.05 & (deg_data$avg_log2FC > 0.25 | deg_data$avg_log2FC < -0.25)) ## Count the number of significant DEGs
  DEG_count <- rbind(DEG_count, data.frame(Cell_Type = cell_type, Num_DEGs = num_degs))  ## Append to the summary table
}

DEG_count$Cell_Type <- factor(DEG_count$Cell_Type, levels = c(  "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte")) 
  
  
# Add DEG information to my seurat object - DEG_count
cell_clusters <- WT_Kcnc1_p180_CB_1step.sct@meta.data$cluster.annot
names(cell_clusters) <- rownames(WT_Kcnc1_p180_CB_1step.sct@meta.data)
DEG_named_vector <- DEG_count$Num_DEGs[match(cell_clusters, DEG_count$Cell_Type)]
names(DEG_named_vector) <- names(cell_clusters)
# Integrate DEG values into the Seurat object
WT_Kcnc1_p180_CB_1step.sct <- AddMetaData(WT_Kcnc1_p180_CB_1step.sct, metadata = DEG_named_vector, col.name = "DEG")
# Create a UMAP plot colored by qval values
pdf("output/seurat/FeaturePlot_WT_Kcnc1_p180_CB_1step_DEG.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = "DEG", pt.size = 0.5, reduction = "umap", label = TRUE) +
  scale_colour_viridis() #  option="magma"
dev.off()
# Add values on the heatmap
## Extract UMAP coordinates
umap_coordinates <- as.data.frame(WT_Kcnc1_p180_CB_1step.sct@reductions$umap@cell.embeddings)
umap_coordinates$cluster <- WT_Kcnc1_p180_CB_1step.sct@meta.data$cluster.annot
## Calculate cluster centers
cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = umap_coordinates, FUN = mean) %>%
  left_join(DEG_count %>% dplyr::rename( "cluster"="Cell_Type"))
## Create a UMAP plot colored by DEG values, with cluster DEG counts as text annotations
pdf("output/seurat/FeaturePlot_WT_Kcnc1_p180_CB_1step_DEGpadj05fc025_numeric.pdf", width=6, height=6)
FeaturePlot(WT_Kcnc1_p180_CB_1step.sct, features = "DEG", pt.size = 0.5, reduction = "umap") +
  scale_colour_viridis(option="mako") + # option="magma"
  geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = Num_DEGs), 
            size = 4, color = "red", fontface = "bold") 
dev.off()






# GSEA plot
library("fgsea")


#### import all clsuter DEGs output :
cluster_types <- c(     "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte")
# Loop over each cluster type to read data and assign to a variable
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  assign(cluster, data)
}

## load list of genes to test
### List1
fgsea_sets <- list(
  PathwaysOfNeurodegeneration = read_table(file = "output/Pathway/geneList_mmu05022.txt")$Genes,
  GOBP_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS = read_table(file = "output/Pathway/geneList_GOBP_NEURON_INTRINSIC_APOPTOTIC.txt")$Genes,
  GOBP_NEUROINFLAMMATORY_RESPONSE = read_table(file = "output/Pathway/geneList_GOBP_NEUROINFLAMMATORY_RESPONSE.txt")$Genes
)
### List2
fgsea_sets <- list(
  GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION = read_table(file = "output/Pathway/geneList_GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION.txt")$Genes,
  GOBP_NEURONAL_ACTION_POTENTIAL = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_ACTION_POTENTIAL.txt")$Genes,
  GOBP_NEURONAL_SIGNAL_TRANSDUCTION = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_SIGNAL_TRANSDUCTION.txt")$Genes,
  REACTOME_POTASSIUM_CHANNELS = read_table(file = "output/Pathway/geneList_REACTOME_POTASSIUM_CHANNELS.txt")$Genes,
  REACTOME_ION_CHANNEL_TRANSPORT = read_table(file = "output/Pathway/geneList_REACTOME_ION_CHANNEL_TRANSPORT.txt")$Genes
)
### List3
fgsea_sets <- list(
  REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING = read_table(file = "output/Pathway/geneList_REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING.txt")$Genes,
  REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS = read_table(file = "output/Pathway/geneList_REACTOME_NA_CL_DEPENDENT_NEUROTRANSMITTER_TRANSPORTERS.txt")$Genes,
  REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION = read_table(file = "output/Pathway/geneList_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION.txt")$Genes,
  REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS = read_table(file = "output/Pathway/geneList_REACTOME_VOLTAGE_GATED_POTASSIUM_CHANNELS.txt")$Genes
)
### List4
fgsea_sets <- list(
  PathwaysOfNeurodegeneration = read_table(file = "output/Pathway/geneList_mmu05022.txt")$Genes,
  GOBP_NEURONAL_ACTION_POTENTIAL = read_table(file = "output/Pathway/geneList_GOBP_NEURONAL_ACTION_POTENTIAL.txt")$Genes,
  REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION = read_table(file = "output/Pathway/geneList_REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION.txt")$Genes
)


## Rank genes based on FC
genes <- Purkinje %>%  ## CHANGE HERE GENE LIST !!!!!!!!!!!!!!!! ##
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene, avg_log2FC)

ranks <- deframe(genes)
head(ranks)
## Run GSEA

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(ES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -NES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()


## plot GSEA
pdf("output/Pathway/GSEA_Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes-mmu05022-Purkinje.pdf", width=5, height=3)

plotEnrichment(fgsea_sets[["PathwaysOfNeurodegeneration"]],
               ranks) + labs(title="PathwaysOfNeurodegeneration-Purkinje") +
               theme_bw()
dev.off()


# Save output table for all pathway and cluster
## Define the list of cluster types
cluster_types <- c(    "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2",
  "PLI12",
  "PLI23",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte")

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
  fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 10000)
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

write.table(final_results, file = c("output/Pathway/gsea_output_Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct-List4.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap all GSEA
pdf("output/Pathway/heatmap_gsea_padj-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct-List4.pdf", width=15, height=3)
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


### if need reorder terms
final_results$pathway <- factor(final_results$pathway, levels = c("PathwaysOfNeurodegeneration", "GOBP_NEURONAL_ACTION_POTENTIAL", "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION")) 

pdf("output/Pathway/heatmap_gsea_padj-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct-List4.pdf", width=10, height=3)
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














## Pseudotime condiments

Let's integrate the three time points together (SCT method, same as replicates/genotypes integration); and perform pseudotime analysis. Granules to focus on.

### Data integration all samples

Let's try:
- *two step integration*: Load our WT_Kcnc1 integrated sample for each time points. And integrate them. Do not re-regress variables.
- *one step integration*: Load all individual samples and integrate them in one step altogether.
- *two step integrateMerge: merge the 1step integrated dataset*: Load our WT_Kcnc1 integrated sample for each time points. Merge them; and then perform analysis. Recommended by Satijalab [here](https://github.com/satijalab/seurat/issues/2468). For the *variable features* I follow [this](https://github.com/satijalab/seurat/issues/5135) and [this](https://github.com/satijalab/seurat/issues/6185); so "set the VariableFeatures to the features in the scale.data slot of the merged SCT assay, which by default includes only the residuals from variable features."


Let's run in an R script the 1/2 step integration

```bash
conda activate scRNAseqV2

sbatch scripts/dataIntegration_CB_1step.sh # 29969813 fail
sbatch scripts/dataIntegration_CB_2step.sh # 29969881 ok

conda activate monocle3
# --> Remvoed some packages loaded; only include the necessary
sbatch scripts/dataIntegration_CB_1step_monocle3Conda.sh # 30052867 fail 

conda activate scRNAseqV2
sbatch scripts/dataIntegration_CB_1step_reductionRPCA.sh # 30202253 error 'pca not found in this seurat object'; added an extra line from seurat tutorial
sbatch scripts/dataIntegration_CB_1step_reductionRPCA.sh # 30223575 ok



sbatch scripts/dataIntegration_CB_2step_integrateMerge.sh # 

```

--> Error `unable to coerce from TsparseMatrix to [CR]sparseMatrixwhen length of 'i' slot exceeds 2^31-1` discussed [here](https://github.com/satijalab/seurat/discussions/7311). Solution discussed:
  - merge replicates, then do integration
  - test using `scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,anchor.features = features, reduction="rpca", dims = 1:30)`. Also show in tutorial [here](https://satijalab.org/seurat/articles/integration_rpca.html)
  - either test Seurat v5; use `conda activate monocle3` with seuratv5 --> FAIL


#### 2step integration



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
library("RColorBrewer")


# Load 2step integration
WT_Kcnc1_CB_2step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_2step-dim50kparam20res02.rds")
## Add the time information (p14, p35, p180) from orig.ident
WT_Kcnc1_CB_2step.sct$time <- stringr::str_extract(WT_Kcnc1_CB_2step.sct$orig.ident, "p\\d+")


DefaultAssay(WT_Kcnc1_CB_2step.sct) <- "integrated"

WT_Kcnc1_CB_2step.sct <- RunPCA(WT_Kcnc1_CB_2step.sct, verbose = FALSE, npcs = 30)
WT_Kcnc1_CB_2step.sct <- RunUMAP(WT_Kcnc1_CB_2step.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
WT_Kcnc1_CB_2step.sct <- FindNeighbors(WT_Kcnc1_CB_2step.sct, reduction = "pca", k.param = 20, dims = 1:30)
WT_Kcnc1_CB_2step.sct <- FindClusters(WT_Kcnc1_CB_2step.sct, resolution = 0.2, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_CB_2step.sct$condition <- factor(WT_Kcnc1_CB_2step.sct$condition, levels = c("WT", "Kcnc1")) 
WT_Kcnc1_CB_2step.sct$time <- factor(WT_Kcnc1_CB_2step.sct$time, levels = c("p14", "p35", "p180")) 

pdf("output/seurat/UMAP_WT_Kcnc1_CB-2stepIntegrationRegressNotRepeatedregMtRbCou-dim30kparam20res02.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_CB_2step.sct, reduction = "umap", label=TRUE, raster = FALSE) # add `raster = FALSE` when more than 100k cells in the plot
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-2stepIntegrationRegressNotRepeatedregMtRbCou-dim30kparam20res02_splitOrig.pdf", width=60, height=6)
DimPlot(WT_Kcnc1_CB_2step.sct, reduction = "umap", label=TRUE, raster = FALSE, split.by = "orig.ident") # add `raster = FALSE` when more than 100k cells in the plot
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-2stepIntegrationRegressNotRepeatedregMtRbCou-dim30kparam20res02_splitTime.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_CB_2step.sct, reduction = "umap", label=TRUE, raster = FALSE, split.by = "time") # add `raster = FALSE` when more than 100k cells in the plot
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-2stepIntegrationRegressNotRepeatedregMtRbCou-dim30kparam20res02_groupTime.pdf", width=5, height=6)
my_cols = brewer.pal(3,"Dark2")
DimPlot(WT_Kcnc1_CB_2step.sct, reduction = "umap", label=TRUE, raster = FALSE, group.by = "time", cols=alpha(my_cols,0.3)) # add `raster = FALSE` when more than 100k cells in the plot
dev.off()



# List678: (*combine list6 (p14), list7(p35), list8(p180)*)
Granule = Gabra6, Pax6
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
PLI12 = Klhl1, Gfra2, Aldh1a3
PLI23 = Galntl6, Kcnc2
Golgi = Pax2
Unipolar_Brush = Eomes, Rgs6, Tafa2
Purkinje = Calb1, Slc1a6, Car8
Bergman_Glia = Aqp4, Slc39a12
Astrocyte = Zeb2
OPC = Aldoc, Cnp,     Vcan, Sox6
Oligodendrocyte= Mbp, Mag, Plp1
Microglia = Itgam, Cx3cr1
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Endothelial_Mural =  Dlc1, Pdgfrb
Endothelial_Stalk = Actb, Tmsb4x
Choroid plexus cells = Kl,  Ttr





DefaultAssay(WT_Kcnc1_CB_2step.sct) <- "SCT" # For vizualization either use SCT or norm RNA


pdf("output/seurat/FeaturePlot_SCT_WT_CB-dim30-List678.pdf", width=25, height=40)
FeaturePlot(WT_Kcnc1_CB_2step.sct, features = c("Gabra6", "Pax6", "Sorcs3", "Ptprk", "Nxph1", "Cdh22",  "Klhl1", "Gfra2", "Aldh1a3", "Galntl6", "Kcnc2", "Pax2", "Eomes", "Rgs6", "Tafa2", "Calb1", "Slc1a6", "Car8", "Aqp4", "Slc39a12", "Zeb2", "Aldoc", "Cnp", "Vcan", "Sox6", "Mbp", "Mag", "Plp1", "Itgam", "Cx3cr1", "Ptgds", "Dcn", "Lef1", "Notum", "Apcdd1", "Dlc1", "Pdgfrb", "Actb", "Tmsb4x", "Kl",  "Ttr"), max.cutoff = 2, cols = c("grey", "red"), raster = FALSE)
dev.off()



```


--> *2step integration*: seems all cell types are put together whatever their time; hard to observe pseudotime trajectories maintenance, notably in Granules...

Let's test *1step integration*


#### 1step integration




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
library("RColorBrewer")


# Load 1step integration
WT_Kcnc1_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_1step-dim50kparam20res02.rds")
## Add the time information (p14, p35, p180) from orig.ident
WT_Kcnc1_CB_1step.sct$time <- stringr::str_extract(WT_Kcnc1_CB_1step.sct$orig.ident, "p\\d+")


DefaultAssay(WT_Kcnc1_CB_1step.sct) <- "integrated"

WT_Kcnc1_CB_1step.sct <- RunPCA(WT_Kcnc1_CB_1step.sct, verbose = FALSE, npcs = 50)
WT_Kcnc1_CB_1step.sct <- RunUMAP(WT_Kcnc1_CB_1step.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
WT_Kcnc1_CB_1step.sct <- FindNeighbors(WT_Kcnc1_CB_1step.sct, reduction = "pca", k.param = 20, dims = 1:50)
WT_Kcnc1_CB_1step.sct <- FindClusters(WT_Kcnc1_CB_1step.sct, resolution = 0.2, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_CB_1step.sct$condition <- factor(WT_Kcnc1_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) 
WT_Kcnc1_CB_1step.sct$time <- factor(WT_Kcnc1_CB_1step.sct$time, levels = c("p14", "p35", "p180")) 

pdf("output/seurat/UMAP_WT_Kcnc1_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-dim50kparam20res02.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_CB_1step.sct, reduction = "umap", label=TRUE, raster = FALSE) # add `raster = FALSE` when more than 100k cells in the plot
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-dim50kparam20res02_splitOrig.pdf", width=60, height=6)
DimPlot(WT_Kcnc1_CB_1step.sct, reduction = "umap", label=TRUE, raster = FALSE, split.by = "orig.ident") # add `raster = FALSE` when more than 100k cells in the plot
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-dim50kparam20res02_splitTime.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_CB_1step.sct, reduction = "umap", label=TRUE, raster = FALSE, split.by = "time") # add `raster = FALSE` when more than 100k cells in the plot
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-dim50kparam20res02_groupTime.pdf", width=5, height=6)
my_cols = brewer.pal(3,"Dark2")
DimPlot(WT_Kcnc1_CB_1step.sct, reduction = "umap", label=TRUE, raster = FALSE, group.by = "time", cols=alpha(my_cols,0.3)) # add `raster = FALSE` when more than 100k cells in the plot
dev.off()



# List678: (*combine list6 (p14), list7(p35), list8(p180)*)
Granule = Gabra6, Pax6
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
PLI12 = Klhl1, Gfra2, Aldh1a3
PLI23 = Galntl6, Kcnc2
Golgi = Pax2
Unipolar_Brush = Eomes, Rgs6, Tafa2
Purkinje = Calb1, Slc1a6, Car8
Bergman_Glia = Aqp4, Slc39a12
Astrocyte = Zeb2
OPC = Aldoc, Cnp,     Vcan, Sox6
Oligodendrocyte= Mbp, Mag, Plp1
Microglia = Itgam, Cx3cr1
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Endothelial_Mural =  Dlc1, Pdgfrb
Endothelial_Stalk = Actb, Tmsb4x
Choroid plexus cells = Kl,  Ttr


DefaultAssay(WT_Kcnc1_CB_1step.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_WT_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-dim50-List678.pdf", width=25, height=40)
FeaturePlot(WT_Kcnc1_CB_1step.sct, features = c("Gabra6", "Pax6", "Sorcs3", "Ptprk", "Nxph1", "Cdh22",  "Klhl1", "Gfra2", "Aldh1a3", "Galntl6", "Kcnc2", "Pax2", "Eomes", "Rgs6", "Tafa2", "Calb1", "Slc1a6", "Car8", "Aqp4", "Slc39a12", "Zeb2", "Aldoc", "Cnp", "Vcan", "Sox6", "Mbp", "Mag", "Plp1", "Itgam", "Cx3cr1", "Ptgds", "Dcn", "Lef1", "Notum", "Apcdd1", "Dlc1", "Pdgfrb", "Actb", "Tmsb4x", "Kl",  "Ttr"), max.cutoff = 2, cols = c("grey", "red"), raster = FALSE)
dev.off()
```










#### integrateMerge


Let's test **integrateMerge option**

--> Not sure about working on SCT assay for the UMAP SCT... Before I used the integrated assay to generate UMAP. I do not think that is a big deal. 



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
library("RColorBrewer")
set.seed(42)

# load seurat object
WT_Kcnc1_CB_integrateMerge.sct <- readRDS(file = "output/seurat/WT_Kcnc1_CB_integrateMerge.sct.rds")
## Add the time information (p14, p35, p180) from orig.ident
WT_Kcnc1_CB_integrateMerge.sct$time <- stringr::str_extract(WT_Kcnc1_CB_integrateMerge.sct$orig.ident, "p\\d+")
# Load the variable features from the previous SCT individual objects
VariableFeatures(WT_Kcnc1_CB_integrateMerge.sct[["SCT"]]) <- rownames(WT_Kcnc1_CB_integrateMerge.sct[["SCT"]]@scale.data)





#### UMAP
DefaultAssay(WT_Kcnc1_CB_integrateMerge.sct) <- "SCT"

WT_Kcnc1_CB_integrateMerge.sct <- RunPCA(WT_Kcnc1_CB_integrateMerge.sct, verbose = FALSE, npcs = 40)
WT_Kcnc1_CB_integrateMerge.sct <- RunUMAP(WT_Kcnc1_CB_integrateMerge.sct, reduction = "pca", dims = 1:40, verbose = FALSE)
WT_Kcnc1_CB_integrateMerge.sct <- FindNeighbors(WT_Kcnc1_CB_integrateMerge.sct, reduction = "pca", k.param = 15, dims = 1:40)
WT_Kcnc1_CB_integrateMerge.sct <- FindClusters(WT_Kcnc1_CB_integrateMerge.sct, resolution = 0.3, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_CB_integrateMerge.sct$condition <- factor(WT_Kcnc1_CB_integrateMerge.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st
WT_Kcnc1_CB_integrateMerge.sct$time <- factor(WT_Kcnc1_CB_integrateMerge.sct$time, levels = c("p14", "p35", "p180")) 

pdf("output/seurat/UMAP_WT_Kcnc1_CB-integrateMerge-dim40kparam15res03.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_CB_integrateMerge.sct, reduction = "umap", label=TRUE, raster = FALSE)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-integrateMerge-dim40kparam15res03_splitOrig.pdf", width=60, height=6)
DimPlot(WT_Kcnc1_CB_integrateMerge.sct, reduction = "umap", label=TRUE, raster = FALSE, split.by = "orig.ident") # add `raster = FALSE` when more than 100k cells in the plot
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-integrateMerge-dim40kparam15res03_splitTime.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_CB_integrateMerge.sct, reduction = "umap", label=TRUE, raster = FALSE, split.by = "time") # add `raster = FALSE` when more than 100k cells in the plot
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-integrateMerge-dim40kparam15res03_groupTime.pdf", width=7, height=6)
my_cols = brewer.pal(3,"Dark2")
DimPlot(WT_Kcnc1_CB_integrateMerge.sct, reduction = "umap", label=TRUE, raster = FALSE, group.by = "time", cols=alpha(my_cols,0.3)) # add `raster = FALSE` when more than 100k cells in the plot
dev.off()



# List678: (*combine list6 (p14), list7(p35), list8(p180)*)
Granule = Gabra6, Pax6
MLI1 = Sorcs3, Ptprk
MLI2 = Nxph1, Cdh22
PLI12 = Klhl1, Gfra2, Aldh1a3
PLI23 = Galntl6, Kcnc2
Golgi = Pax2
Unipolar_Brush = Eomes, Rgs6, Tafa2
Purkinje = Calb1, Slc1a6, Car8
Bergman_Glia = Aqp4, Slc39a12
Astrocyte = Zeb2
OPC = Aldoc, Cnp,     Vcan, Sox6
Oligodendrocyte= Mbp, Mag, Plp1
Microglia = Itgam, Cx3cr1
Meningeal = Ptgds, Dcn
Endothelial = Lef1, Notum, Apcdd1
Endothelial_Mural =  Dlc1, Pdgfrb
Endothelial_Stalk = Actb, Tmsb4x
Choroid plexus cells = Kl,  Ttr


DefaultAssay(WT_Kcnc1_CB_integrateMerge.sct) <- "SCT" # For vizualization either use SCT or norm RNA

pdf("output/seurat/FeaturePlot_SCT_WT_CB-integrateMerge-dim40-List678.pdf", width=25, height=40)
FeaturePlot(WT_Kcnc1_CB_integrateMerge.sct, features = c("Gabra6", "Pax6", "Sorcs3", "Ptprk", "Nxph1", "Cdh22",  "Klhl1", "Gfra2", "Aldh1a3", "Galntl6", "Kcnc2", "Pax2", "Eomes", "Rgs6", "Tafa2", "Calb1", "Slc1a6", "Car8", "Aqp4", "Slc39a12", "Zeb2", "Aldoc", "Cnp", "Vcan", "Sox6", "Mbp", "Mag", "Plp1", "Itgam", "Cx3cr1", "Ptgds", "Dcn", "Lef1", "Notum", "Apcdd1", "Dlc1", "Pdgfrb", "Actb", "Tmsb4x", "Kl",  "Ttr"), max.cutoff = 2, cols = c("grey", "red"), raster = FALSE)
dev.off()



# V1 naming
cluster1: Granule_1
cluster2: PLI
cluster3: Granule_2
cluster4: Mix
cluster5: MLI1_1
cluster6: MLI1_2
cluster7: MLI2_1
cluster8: MLI1_3
cluster9: Granule_3
cluster10: Astrocyte
cluster11: Endothelial
cluster12: MLI2_2
cluster13: BergmanGlia
cluster14: MLI2_3
cluster15: UnipolarBrush
cluster16: Meningeal
cluster17: Golgi
cluster18: Purkinje
cluster19: ChoroidPlexusCells
cluster20: OPC
cluster21: Oligodendrocyte
cluster22: Unknown


new.cluster.ids <- c(
  "Granule_1",
  "PLI",
  "Granule_2",
  "Mix",
  "MLI1_1",
  "MLI1_2",
  "MLI2_1",
  "MLI1_3",
  "Granule_3",
  "Astrocyte",
  "Endothelial",
  "MLI2_2",
  "BergmanGlia",
  "MLI2_3",
  "UnipolarBrush",
  "Meningeal",
  "Golgi",
  "Purkinje",
  "ChoroidPlexusCells",
  "OPC",
  "Oligodendrocyte",
  "Unknown"
)

names(new.cluster.ids) <- levels(WT_Kcnc1_CB_integrateMerge.sct)
WT_Kcnc1_CB_integrateMerge.sct <- RenameIdents(WT_Kcnc1_CB_integrateMerge.sct, new.cluster.ids)
WT_Kcnc1_CB_integrateMerge.sct$cluster.annot <- Idents(WT_Kcnc1_CB_integrateMerge.sct) # create a new slot in my seurat object


pdf("output/seurat/UMAP_WT_Kcnc1_CB-integrateMerge-dim40kparam15res03_label.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_CB_integrateMerge.sct, reduction = "umap", split.by = "time", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3, raster = FALSE)
dev.off()

pdf("output/seurat/UMAP_WT_Kcnc1_CB-integrateMerge-dim40kparam15res03_noSplit_label.pdf", width=9, height=6)
DimPlot(WT_Kcnc1_CB_integrateMerge.sct, reduction = "umap",  label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5, raster = FALSE)
dev.off()

WT_Kcnc1_CB_integrateMerge.sct$condition <- factor(WT_Kcnc1_CB_integrateMerge.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st


pdf("output/seurat/UMAP_WT_Kcnc1_CB-integrateMerge-dim40kparam15res03_splitCondition_label.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_CB_integrateMerge.sct, reduction = "umap", split.by = "condition", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 3, raster = FALSE)
dev.off()



# save ######################################################
## saveRDS(WT_Kcnc1_CB_integrateMerge.sct, file = "output/seurat/WT_Kcnc1_CB_integrateMerge-dim40kparam15res03-labelv1.rds") # 
############################################################################################################

```

--> Very good! Clearly show pattern of cell types from>to p14>p35>p180; notably for Granule, MLI1, MLI2







### Pseudotime


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
library("readr")

# Separate WT Kcnc1 genotypes
humangastruloid_UNTREATED.combined.sct <- subset(humangastruloid.combined.sct, treatment == "UNTREATED")
humangastruloid_DASATINIB.combined.sct <- subset(humangastruloid.combined.sct, treatment == "DASATINIB")


################################################################################################
### CONTROL/UNTREATED ################################################
################################################################################################

DefaultAssay(humangastruloid_UNTREATED.combined.sct) <- "RNA" # According to condiments workflow

# convert to SingleCellExperiment
humangastruloidUNTREATED <- as.SingleCellExperiment(humangastruloid_UNTREATED.combined.sct, assay = "RNA")

#  Trajectory Inference and Differential Topology
set.seed(42)


### Testing area ############


##########################################
"Mixed_Epiblast_Ectoderm_PrimitiveStreak", "Endoderm", "Cardiomyocyte", "ProliferatingCardiacMesoderm", 

humangastruloidUNTREATED <- slingshot(humangastruloidUNTREATED, reducedDim = 'UMAP',
                 clusterLabels = colData(humangastruloidUNTREATED)$cluster.annot,
                 start.clus = 'Epiblast', end.clus = c("CPC2", "CPC1", "Endoderm", "Cardiomyocyte"), approx_points = 200, extend = "n", stretch = 1.6)





df_2 <- bind_cols(
  as.data.frame(reducedDim(humangastruloidUNTREATED, "UMAP")),
  slingPseudotime(humangastruloidUNTREATED) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(humangastruloidUNTREATED) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         pst = if_else(Lineage1 > Lineage2, Lineage1_pst, Lineage2_pst),
        # pst = max(pst) - pst)
)
curves <- slingCurves(humangastruloidUNTREATED, as.df = TRUE)




pdf("output/condiments/UMAP_trajectory_common_humangastruloidUNTREATED_StartEpiEndCPC12EndoCardioextendnstretch16approx200.pdf", width=5, height=5)
ggplot(df_2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = .7, aes(col = pst)) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage), col = "black", size = 1.5) +
  theme_classic()
dev.off()




## PLOT with trajectories - Individually
df_2 <- bind_cols(
  as.data.frame(reducedDim(humangastruloidUNTREATED, "UMAP")),
  slingPseudotime(humangastruloidUNTREATED) %>% as.data.frame() %>%
    dplyr::rename_with(paste0, "_pst", .cols = everything()),
  slingCurveWeights(humangastruloidUNTREATED) %>% as.data.frame(),
  ) %>%
  mutate(Lineage1_pst = if_else(is.na(Lineage1_pst), 0, Lineage1_pst),
         Lineage2_pst = if_else(is.na(Lineage2_pst), 0, Lineage2_pst),
         Lineage3_pst = if_else(is.na(Lineage3_pst), 0, Lineage3_pst),
         Lineage4_pst = if_else(is.na(Lineage4_pst), 0, Lineage4_pst),
         Lineage5_pst = if_else(is.na(Lineage5_pst), 0, Lineage5_pst),
         Lineage6_pst = if_else(is.na(Lineage6_pst), 0, Lineage6_pst))
curves <- slingCurves(humangastruloidUNTREATED, as.df = TRUE)
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
    geom_point(data = df_2_other, aes(x = UMAP_1, y = UMAP_2), size = .7, color = "grey85") +
    geom_point(data = df_2_lineage, aes(x = UMAP_1, y = UMAP_2, col = pst), size = .7) +
    scale_color_viridis_c() +
    labs(col = "Pseudotime", title = paste("Lineage", lineage_number)) +
    geom_path(data = curves_filtered %>% arrange(Order),
              aes(x = UMAP_1, y = UMAP_2, group = Lineage), col = "black", size = 1) +
    geom_text(data = curves_endpoints, aes(x = UMAP_1, y = UMAP_2, label = Lineage), size = 4, vjust = -1, hjust = -1, col = "red") +  # Use endpoints for labels
    theme_classic()
  return(p)
}
### Generate the plots for each lineage
plots <- list()
for (i in 1:6) {
  plots[[i]] <- create_plot(i)
}
pdf("output/condiments/UMAP_trajectory_common_label_humangastruloidUNTREATED_StartEpiEndCPC12EndoCardioextendnstretch16approx200_Lineage123456.pdf", width=25, height=5)
gridExtra::grid.arrange(grobs = plots, ncol = 6)
dev.off()






#### ->  save.image(file="output/condiments/condiments_humangastruloidUNTREATED2472hrs_StartEpiEndCPC12EndoCardioextendnstretch16approx200.RData")
### load("output/condiments/condiments_humangastruloidUNTREATED2472hrs_StartEpiEndCPC12EndoCardioextendnstretch16approx200.RData")
set.seed(42)
##
# RUN fitGAM
##

```



## Slurm jobs

### DEGs WT vs Kcnc1 cluster per cluster


```bash
conda activate scRNAseqV2

# p14 CB
sbatch scripts/DEG_allGenes_WT_Kcnc1_p14_CB.sh # 27801074 ok
sbatch scripts/DEG_allGenes_WT_Kcnc1_p14_CB_correct.sh # 29330177 ok

# p35 CB
sbatch scripts/DEG_allGenes_WT_Kcnc1_p35_CB.sh # 27801335 ok
sbatch scripts/DEG_allGenes_WT_Kcnc1_p35_CB_correct.sh # 29328955 ok

# p180 CB
sbatch scripts/DEG_allGenes_WT_Kcnc1_p180_CB.sh # 27801489 ok
sbatch scripts/DEG_allGenes_WT_Kcnc1_p180_CB_correct.sh # 29330555 ok
```



- *NOTE: I forget to LogNorm and Scale prior doing DEG; the `*_correct` are corrected + updated cluster name notbably PLI12 and PLI23 instead of MLI_1 MLI_2; all good.*



## Comparison effect logNorm/ScaleData after integration

I realized I did not apply logNorm/ScaleData normalization after data integration and prior doing the DGE. For p14 and p180, analysis showed a comparable number of DGEs up/down. However, all genes where downregulated in p35 sample, that's where I realized logNorm/ScaleData was not applied. Let's compare with and without logNorm/ScaleData how it affect the nb of DEGs and direction.



```R
# packages
library("tidyverse")


# import list of DEGs

## p14 ########
### without logNorm/ScaleData ################################
cluster_types <- c(
  "Granular_1", "Granular_2", "Granular_3", "Granular_4", "Granular_5",
  "Interneuron", "MLI1", "MLI2", "PLI", "Golgi", "Unipolar_Brush",
  "Purkinje", "Astrocyte", "Bergmann_Glia", "Oligodendrocyte", "OPC",
  "Mix_Microglia_Meningeal", "Endothelial", "Endothelial_Mural", "Choroid_Plexus"
)
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1) %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  data$cluster <- cluster
  assign(cluster, data)
}

cluster_data_list <- mget(cluster_types)
p14 <- bind_rows(cluster_data_list) # Bind all cluster data frames into a single data frame
#### Function to count up- and down-regulated genes in each cluster
count_up_down <- function(cluster_data) {
  significant_genes <- cluster_data %>% filter(p_val_adj < 0.05)
  nb_upregulated <- sum(significant_genes$avg_log2FC > 0)
  nb_downregulated <- sum(significant_genes$avg_log2FC < 0)
  return(data.frame(nb_upregulated = nb_upregulated, nb_downregulated = nb_downregulated))
}

clusters_split <- split(p14, p14$cluster)
results <- lapply(clusters_split, count_up_down)
final_results <- do.call(rbind, results)
row.names(final_results) <- names(clusters_split)
print(final_results)


### with logNorm/ScaleData (*_correct)
cluster_types <- c(
  "PLI23",
  "Granular_1",
  "MLI1",
  "Granular_2",
  "Granular_3",
  "Granular_4",
  "MLI2",
  "Endothelial",
  "Granular_5",
  "Astrocyte",
  "OPC",
  "Bergmann_Glia",
  "PLI12",
  "Oligodendrocyte",
  "Mix_Microglia_Meningeal",
  "Endothelial_Mural" ,
  "Purkinje",
  "Golgi",
  "Unipolar_Brush",
  "Choroid_Plexus"
)
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p14_CB_QCV3dim30kparam50res035_allGenes_correct.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1) %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  data$cluster <- cluster
  assign(cluster, data)
}

cluster_data_list <- mget(cluster_types)
p14_correct <- bind_rows(cluster_data_list) # Bind all cluster data frames into a single data frame
#### Function to count up- and down-regulated genes in each cluster
count_up_down <- function(cluster_data) {
  significant_genes <- cluster_data %>% filter(p_val_adj < 0.05)
  nb_upregulated <- sum(significant_genes$avg_log2FC > 0)
  nb_downregulated <- sum(significant_genes$avg_log2FC < 0)
  return(data.frame(nb_upregulated = nb_upregulated, nb_downregulated = nb_downregulated))
}

clusters_split <- split(p14_correct, p14_correct$cluster)
results <- lapply(clusters_split, count_up_down)
final_results <- do.call(rbind, results)
row.names(final_results) <- names(clusters_split)
print(final_results)






## p35 ########
### without logNorm/ScaleData ################################
cluster_types <- c(
  "Granular",
  "Interneuron",
  "MLI1",
  "MLI2_1",
  "MLI2_2",
  "Golgi",
  "Unipolar_Brush",
  "Purkinje",
  "Astrocyte",
  "Bergman_Glia",
  "OPC",
  "Meningeal",
  "Endothelial",
  "Choroid_Plexus",
  "Endothelial_Stalk"
)
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p35_CB_QCV3dim50kparam50res03_allGenes.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1) %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  data$cluster <- cluster
  assign(cluster, data)
}

cluster_data_list <- mget(cluster_types)
p35 <- bind_rows(cluster_data_list) # Bind all cluster data frames into a single data frame
#### Function to count up- and down-regulated genes in each cluster
count_up_down <- function(cluster_data) {
  significant_genes <- cluster_data %>% filter(p_val_adj < 0.05)
  nb_upregulated <- sum(significant_genes$avg_log2FC > 0)
  nb_downregulated <- sum(significant_genes$avg_log2FC < 0)
  return(data.frame(nb_upregulated = nb_upregulated, nb_downregulated = nb_downregulated))
}

clusters_split <- split(p35, p35$cluster)
results <- lapply(clusters_split, count_up_down)
final_results <- do.call(rbind, results)
row.names(final_results) <- names(clusters_split)
print(final_results)


### with logNorm/ScaleData (*_correct)
cluster_types <- c(
  "Granule",
  "MLI1",
  "PLI23" ,
  "MLI2" ,
  "Endothelial_Stalk",
  "Astrocyte" ,
  "PLI12" ,
  "Golgi" ,
  "Unipolar_Brush" ,
  "Bergman_Glia",
  "Endothelial",
  "Choroid_Plexus",
  "Purkinje",
  "Meningeal",
  "Unknown_Neuron_Subpop",
  "OPC" 
)
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p35_CB_QCV2dim50kparam20res03_allGenes_correct.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1) %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  data$cluster <- cluster
  assign(cluster, data)
}

cluster_data_list <- mget(cluster_types)
p35_correct <- bind_rows(cluster_data_list) # Bind all cluster data frames into a single data frame
#### Function to count up- and down-regulated genes in each cluster
count_up_down <- function(cluster_data) {
  significant_genes <- cluster_data %>% filter(p_val_adj < 0.05)
  nb_upregulated <- sum(significant_genes$avg_log2FC > 0)
  nb_downregulated <- sum(significant_genes$avg_log2FC < 0)
  return(data.frame(nb_upregulated = nb_upregulated, nb_downregulated = nb_downregulated))
}

clusters_split <- split(p35_correct, p35_correct$cluster)
results <- lapply(clusters_split, count_up_down)
final_results <- do.call(rbind, results)
row.names(final_results) <- names(clusters_split)
print(final_results)










## p180 ########
### without logNorm/ScaleData ################################
cluster_types <- c(
  "Granular_1",
  "Granular_2",
  "Granular_3",
  "MLI1",
  "MLI2_1",
  "MLI2_2",
  "Interneuron",
  "Astrocyte",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1) %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  data$cluster <- cluster
  assign(cluster, data)
}

cluster_data_list <- mget(cluster_types)
p180 <- bind_rows(cluster_data_list) # Bind all cluster data frames into a single data frame
#### Function to count up- and down-regulated genes in each cluster
count_up_down <- function(cluster_data) {
  significant_genes <- cluster_data %>% filter(p_val_adj < 0.05)
  nb_upregulated <- sum(significant_genes$avg_log2FC > 0)
  nb_downregulated <- sum(significant_genes$avg_log2FC < 0)
  return(data.frame(nb_upregulated = nb_upregulated, nb_downregulated = nb_downregulated))
}

clusters_split <- split(p180, p180$cluster)
results <- lapply(clusters_split, count_up_down)
final_results <- do.call(rbind, results)
row.names(final_results) <- names(clusters_split)
print(final_results)


### with logNorm/ScaleData (*_correct)
cluster_types <- c(
  "Granular_1",
  "Granular_2",
  "MLI1",
  "Granular_3",
  "MLI2",
  "PLI23",
  "Astrocyte",
  "PLI12",
  "Bergman_Glia",
  "Unipolar_Brush",
  "Mix_Endothelial_EndothelialMural",
  "Meningeal",
  "Choroid_Plexus",
  "Golgi",
  "Purkinje",
  "Unknown_Neuron_Subpop",
  "Oligodendrocyte"
)
for (cluster in cluster_types) {
  file_path <- paste0("output/seurat/", cluster, "-Kcnc1_response_p180_CB_QCV4dim50kparam20res02_allGenes_correct.txt")
  data <- read.delim(file_path, header = TRUE, row.names = 1) %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  data$cluster <- cluster
  assign(cluster, data)
}

cluster_data_list <- mget(cluster_types)
p180_correct <- bind_rows(cluster_data_list) # Bind all cluster data frames into a single data frame
#### Function to count up- and down-regulated genes in each cluster
count_up_down <- function(cluster_data) {
  significant_genes <- cluster_data %>% filter(p_val_adj < 0.05)
  nb_upregulated <- sum(significant_genes$avg_log2FC > 0)
  nb_downregulated <- sum(significant_genes$avg_log2FC < 0)
  return(data.frame(nb_upregulated = nb_upregulated, nb_downregulated = nb_downregulated))
}

clusters_split <- split(p180_correct, p180_correct$cluster)
results <- lapply(clusters_split, count_up_down)
final_results <- do.call(rbind, results)
row.names(final_results) <- names(clusters_split)
print(final_results)




```






# ShinyApp

Let's create a shiny app to make gene search easier!
- WT_Kcnc1_p14_CB_1step_V1 = 1stepIntegrationRegressNotRepeatedregMtRbCou-QCV2dim30kparam30res05



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

WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V1_numeric.rds")
DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "RNA" # 

# Generate Shiny app V1
scConf = createConfig(WT_Kcnc1_p14_CB_1step.sct)

makeShinyApp(WT_Kcnc1_p14_CB_1step.sct, scConf, gene.mapping = TRUE,
             shiny.title = "WT_Kcnc1_p14_CB_1step_V1",
             shiny.dir = "shinyApp_WT_Kcnc1_p14_CB_1step_V1/") 

rsconnect::deployApp('shinyApp_WT_Kcnc1_p14_CB_1step_V1')



# Generate Shiny app QCV3
DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "RNA" # 

scConf = createConfig(WT_Kcnc1_p14_CB_1step.sct)

makeShinyApp(WT_Kcnc1_p14_CB_1step.sct, scConf, gene.mapping = TRUE,
             shiny.title = "WT_Kcnc1_p14_CB_1step_QCV3",
             shiny.dir = "shinyApp_WT_Kcnc1_p14_CB_1step_QCV3/") 

rsconnect::deployApp('shinyApp_WT_Kcnc1_p14_CB_1step_QCV3')


# Generate Shiny app QCV3 with name V2 (PLI12 PLI23);  QCV3dim30kparam50res035 ; output/seurat/WT_Kcnc1_p14_CB_1step.sct_V5_numeric.rds
DefaultAssay(WT_Kcnc1_p14_CB_1step.sct) <- "RNA" # 

scConf = createConfig(WT_Kcnc1_p14_CB_1step.sct)

makeShinyApp(WT_Kcnc1_p14_CB_1step.sct, scConf, gene.mapping = TRUE,
             shiny.title = "WT_Kcnc1_p14_CB_1step_QCV3dim30kparam50res035",
             shiny.dir = "shinyApp_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam50res035/") 

rsconnect::deployApp('shinyApp_WT_Kcnc1_p14_CB_1step_QCV3dim30kparam50res035')



# Generate Shiny app QCV3 without name;  QCV3dim50kparam50res03 ; output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_numeric.rds
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "RNA" # 

scConf = createConfig(WT_Kcnc1_p35_CB_1step.sct)

makeShinyApp(WT_Kcnc1_p35_CB_1step.sct, scConf, gene.mapping = TRUE,
             shiny.title = "WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03",
             shiny.dir = "shinyApp_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03/") 

rsconnect::deployApp('shinyApp_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03')
##### --> This app as been deleted

# Generate Shiny app QCV3 with name V1;  QCV3dim50kparam50res03 ; output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam50res03.sct_V1_label.rds
DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "RNA" # 

scConf = createConfig(WT_Kcnc1_p35_CB_1step.sct)

makeShinyApp(WT_Kcnc1_p35_CB_1step.sct, scConf, gene.mapping = TRUE,
             shiny.title = "WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03",
             shiny.dir = "shinyApp_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03/") 

rsconnect::deployApp('shinyApp_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam50res03')



# Generate Shiny app QCV3 without name;  QCV3dim40kparam10res03 ; output/seurat/WT_Kcnc1_p180_CB_1step-QCV3dim40kparam10res03.sct_V1_numeric.rds
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV3dim40kparam10res03.sct_V1_numeric.rds")


DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "RNA" # 

scConf = createConfig(WT_Kcnc1_p180_CB_1step.sct)

makeShinyApp(WT_Kcnc1_p180_CB_1step.sct, scConf, gene.mapping = TRUE,
             shiny.title = "WT_Kcnc1_p180_CB_1step_QCV3dim40kparam10res03",
             shiny.dir = "shinyApp_WT_Kcnc1_p180_CB_1step_QCV3dim40kparam10res03/") 

rsconnect::deployApp('shinyApp_WT_Kcnc1_p180_CB_1step_QCV3dim40kparam10res03')


# Generate Shiny app QCV4 with name V1 (PLI12 PLI23);  QCV4dim50kparam20res02 ; output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V1_label.rds
WT_Kcnc1_p180_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p180_CB_1step-QCV4dim50kparam20res02.sct_V1_label.rds")

DefaultAssay(WT_Kcnc1_p180_CB_1step.sct) <- "RNA" # 

scConf = createConfig(WT_Kcnc1_p180_CB_1step.sct)

makeShinyApp(WT_Kcnc1_p180_CB_1step.sct, scConf, gene.mapping = TRUE,
             shiny.title = "WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02",
             shiny.dir = "shinyApp_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02/") 

rsconnect::deployApp('shinyApp_WT_Kcnc1_p180_CB_1step_QCV4dim50kparam20res02')


# Generate Shiny app QCV2 with name V2 (PLI12 PLI23);  QCV2dim50kparam20res03 ; output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V2_label_ReProcess.rds
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV2dim50kparam20res03.sct_V2_label_ReProcess.rds")

DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "RNA" # 

scConf = createConfig(WT_Kcnc1_p35_CB_1step.sct)

makeShinyApp(WT_Kcnc1_p35_CB_1step.sct, scConf, gene.mapping = TRUE,
             shiny.title = "WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03",
             shiny.dir = "shinyApp_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03/") 

rsconnect::deployApp('shinyApp_WT_Kcnc1_p35_CB_1step_QCV2dim50kparam20res03')



```



- [CB_p14](https://roulethomas.shinyapps.io/shinyapp_wt_kcnc1_p14_cb_1step_qcv3dim30kparam50res035/)
- [CB_p35](https://roulethomas.shinyapps.io/shinyapp_wt_kcnc1_p35_cb_1step_qcv3dim50kparam50res03/)
- [CB_p180](https://roulethomas.shinyapps.io/shinyapp_wt_kcnc1_p180_cb_1step_qcv4dim50kparam20res02/)




# Generate bigwig

Generate coverage bigwig files to then generate PCA with deeptools; maybe it will highlight the outliers replicate and suggest more easily that we should get rid of them?




--> Let's just generate the bigwig from the `possorted_genome_bam.bam` generated by 10X cellranger count 



```bash
conda activate deeptools

# raw
sbatch scripts/bamtobigwig_p14_CB_WT.sh # 27211952 ok
sbatch scripts/bamtobigwig_p14_CB_Kcnc1.sh # 27211957 ok

sbatch scripts/bamtobigwig_p35_CB_WT.sh # 27211960 ok
sbatch scripts/bamtobigwig_p35_CB_Kcnc1.sh # 27211962 ok

sbatch scripts/bamtobigwig_p180_CB_WT.sh # 27242635 xxx
sbatch scripts/bamtobigwig_p180_CB_Kcnc1.sh # 27242642 xxx

# BPM norm (=TPM)
sbatch scripts/bamtobigwig_BPMnorm_p14_CB_WT.sh # 27211977 ok
sbatch scripts/bamtobigwig_BPMnorm_p14_CB_Kcnc1.sh # 27211988 ok

sbatch scripts/bamtobigwig_BPMnorm_p35_CB_WT.sh # 27211986 ok
sbatch scripts/bamtobigwig_BPMnorm_p35_CB_Kcnc1.sh # 27212007 ok

sbatch scripts/bamtobigwig_BPMnorm_p180_CB_WT.sh # 27242775 xxx
sbatch scripts/bamtobigwig_BPMnorm_p180_CB_Kcnc1.sh # 27242783 xxx

```


Generate PCA plots


```bash
# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_BPMnorm_p14.sh # 27241818 ok

# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_BPMnorm_p14.npz \
    --transpose \
    --ntop 0 \
    --labels WT_p14_CB_Rep1 WT_p14_CB_Rep2 WT_p14_CB_Rep3 Kcnc1_p14_CB_Rep1 Kcnc1_p14_CB_Rep2 Kcnc1_p14_CB_Rep3 \
    --colors black black black blue blue blue \
    --markers o x s o x s \
    -o output/bigwig/multiBigwigSummary_BPMnorm_p14_plotPCA.pdf
## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_BPMnorm_p14.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_p14_CB_Rep1 WT_p14_CB_Rep2 WT_p14_CB_Rep3 Kcnc1_p14_CB_Rep1 Kcnc1_p14_CB_Rep2 Kcnc1_p14_CB_Rep3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_BPMnorm_p14_heatmap.pdf


# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_BPMnorm_p35.sh # 27241898 ok

# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_BPMnorm_p35.npz \
    --transpose \
    --ntop 0 \
    --labels WT_p35_CB_Rep1 WT_p35_CB_Rep2 WT_p35_CB_Rep3 Kcnc1_p35_CB_Rep1 Kcnc1_p35_CB_Rep2 Kcnc1_p35_CB_Rep3 \
    --colors black black black blue blue blue \
    --markers o x s o x s \
    -o output/bigwig/multiBigwigSummary_BPMnorm_p35_plotPCA.pdf
## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_BPMnorm_p35.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_p35_CB_Rep1 WT_p35_CB_Rep2 WT_p35_CB_Rep3 Kcnc1_p35_CB_Rep1 Kcnc1_p35_CB_Rep2 Kcnc1_p35_CB_Rep3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_BPMnorm_p35_heatmap.pdf






# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_raw_p14.sh # 27241962 xxx

# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_raw_p14.npz \
    --transpose \
    --ntop 0 \
    --labels WT_p14_CB_Rep1 WT_p14_CB_Rep2 WT_p14_CB_Rep3 Kcnc1_p14_CB_Rep1 Kcnc1_p14_CB_Rep2 Kcnc1_p14_CB_Rep3 \
    --colors black black black blue blue blue \
    --markers o x s o x s \
    -o output/bigwig/multiBigwigSummary_raw_p14_plotPCA.pdf
## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_raw_p14.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_p14_CB_Rep1 WT_p14_CB_Rep2 WT_p14_CB_Rep3 Kcnc1_p14_CB_Rep1 Kcnc1_p14_CB_Rep2 Kcnc1_p14_CB_Rep3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_raw_p14_heatmap.pdf


# Generate compile bigwig (.npz) files
sbatch scripts/multiBigwigSummary_raw_p35.sh # 27242034 xxx

# Plot
## PCA
plotPCA -in output/bigwig/multiBigwigSummary_raw_p35.npz \
    --transpose \
    --ntop 0 \
    --labels WT_p14_CB_Rep1 WT_p14_CB_Rep2 WT_p14_CB_Rep3 Kcnc1_p14_CB_Rep1 Kcnc1_p14_CB_Rep2 Kcnc1_p14_CB_Rep3 \
    --colors black black black blue blue blue \
    --markers o x s o x s \
    -o output/bigwig/multiBigwigSummary_raw_p35_plotPCA.pdf
## Heatmap
plotCorrelation \
    -in output/bigwig/multiBigwigSummary_raw_p35.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation" \
    --removeOutliers \
    --labels WT_p14_CB_Rep1 WT_p14_CB_Rep2 WT_p14_CB_Rep3 Kcnc1_p14_CB_Rep1 Kcnc1_p14_CB_Rep2 Kcnc1_p14_CB_Rep3 \
    --whatToPlot heatmap --colorMap bwr --plotNumbers \
    -o output/bigwig/multiBigwigSummary_raw_p35_heatmap.pdf




```



