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



##### QC filtering _ V2 ############################################

### V1 not super stringeant; mit > 15, rb > 22 and RNAfeaure 100 #############################################
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



##### QC filtering _ V3 ############################################

### V1 not super stringeant; mit > 15, rb > 10 and RNAfeaure 400 #############################################
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
write.table(qc_summary_combined, file = "output/seurat/QC_summary_V3.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

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
write.table(phase_summary_combined, file = "output/seurat/CellCyclePhase_V3.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





############ SAVE samples #################################################################
save_seurat_objects <- function(seurat_objects_list, output_dir) {
  for (sample_name in names(seurat_objects_list)) {
    file_name <- paste0(output_dir, sample_name, "_V3_numeric.rds")
    saveRDS(seurat_objects_list[[sample_name]], file = file_name)
  }
}
output_dir <- "output/seurat/"
## Call the function to save the Seurat objects
save_seurat_objects(seurat_objects, output_dir)
###############################################################################################
############# READ samples  ################################################
# Function to load Seurat objects
load_seurat_objects <- function(file_paths) {
  seurat_objects <- list()
  for (file_path in file_paths) {
    sample_name <- gsub("_V1_numeric.rds", "", basename(file_path))
    seurat_objects[[sample_name]] <- readRDS(file_path)
  }
  return(seurat_objects)
}
output_dir <- "output/seurat/"
file_paths <- list.files(output_dir, pattern = "_V1_numeric.rds$", full.names = TRUE)
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
WT_Kcnc1_p14_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p14_CB_1step.sct, reduction = "pca", k.param = 20, dims = 1:30)
WT_Kcnc1_p14_CB_1step.sct <- FindClusters(WT_Kcnc1_p14_CB_1step.sct, resolution = 0.4, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p14_CB_1step.sct$condition <- factor(WT_Kcnc1_p14_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim30kparam20res04.pdf", width=7, height=6)
pdf("output/seurat/UMAP_WT_Kcnc1-test.pdf", width=7, height=6)

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

WT_Kcnc1_p14_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p14_CB_1step.sct_V3_numeric.rds")
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
Kcnc1_p35_CB_Rep2$replicate <- "Rep2"
Kcnc1_p35_CB_Rep3$replicate <- "Rep3"

Kcnc1_p35_CB_Rep1$condition <- "Kcnc1"
Kcnc1_p35_CB_Rep2$condition <- "Kcnc1"
Kcnc1_p35_CB_Rep3$condition <- "Kcnc1"

set.seed(42)


# Replicate and Genotype integration (1 step integration)
## WT Rep

### Reg v1 _ better than Regv2
WT_p35_CB_Rep1 <- SCTransform(WT_p35_CB_Rep1, method = "glmGamPoi", ncells = 8496, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep2 <- SCTransform(WT_p35_CB_Rep2, method = "glmGamPoi", ncells = 11577, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
WT_p35_CB_Rep3 <- SCTransform(WT_p35_CB_Rep3, method = "glmGamPoi", ncells = 14582, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep1 <- SCTransform(Kcnc1_p35_CB_Rep1, method = "glmGamPoi", ncells = 11642, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep2 <- SCTransform(Kcnc1_p35_CB_Rep2, method = "glmGamPoi", ncells = 33698, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 
Kcnc1_p35_CB_Rep3 <- SCTransform(Kcnc1_p35_CB_Rep3, method = "glmGamPoi", ncells = 18205, verbose = TRUE, variable.features.n = 3000, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb")) 




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
WT_Kcnc1_p35_CB_1step.sct <- FindNeighbors(WT_Kcnc1_p35_CB_1step.sct, reduction = "pca", k.param = 10, dims = 1:50)
WT_Kcnc1_p35_CB_1step.sct <- FindClusters(WT_Kcnc1_p35_CB_1step.sct, resolution = 0.3, verbose = FALSE, algorithm = 4, method = "igraph") # method = "igraph" needed for large nb of cells


WT_Kcnc1_p35_CB_1step.sct$condition <- factor(WT_Kcnc1_p35_CB_1step.sct$condition, levels = c("WT", "Kcnc1")) # Reorder untreated 1st

pdf("output/seurat/UMAP_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim50kparam10res03.pdf", width=7, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE)
dev.off()



# genes

DefaultAssay(WT_Kcnc1_p35_CB_1step.sct) <- "SCT"


pdf("output/seurat/FeaturePlot_SCT_WT_Kcnc1_p35_CB-1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim50-List5.pdf", width=30, height=70)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, features = c("Gabra6", "Pax6", "Ntng1", "Grm5", "Sorcs3", "Ptprk", "Nxph1", "Cdh22", "Lef1", "Notum", "Apcdd1", "Zeb2", "Hepacam", "Aqp4", "Slc39a12", "Aldoc", "Cnp", "Cspg5", "Rgs6", "Kcnc2", "Mbp", "Mag", "Plp1", "Cfap54", "Ccdc153", "Cfap44", "Tmem212", "Ptgds", "Dcn", "Dlc1", "Pdgfrb", "Calb1", "Slc1a6", "Car8", "Lmx1a", "Adcy2", "Mdga1", "Eomes", "Kl", "Clic6", "Slc13a4", "Ttr"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()







#### QC metrics investigation ####################################
pdf("output/seurat/VlnPlot_QCmetrics_SCT_WT_Kcnc1_p35_CB_1step-1stepIntegrationRegressNotRepeated-QCV3dim50kparam15res04-countMtRbRegression.pdf", width=20, height=5)
VlnPlot(WT_Kcnc1_p35_CB_1step.sct,features = c("percent.mt", "percent.rb","nCount_RNA","nFeature_RNA","S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nFeature_RNA_SCT_WT_Kcnc1_p35_CB_1step-1stepIntegrationRegressNotRepeated-QCV3dim50kparam15res04-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p35_CB_1step.sct,features = c("nFeature_RNA")) +
  ylim(0,2000)
dev.off()

pdf("output/seurat/VlnPlot_QCmetrics_nCount_RNA_SCT_WT_Kcnc1_p35_CB_1step-1stepIntegrationRegressNotRepeated-QCV3dim50kparam15res04-countMtRbRegression.pdf", width=5, height=5)
VlnPlot(WT_Kcnc1_p35_CB_1step.sct,features = c("nCount_RNA")) +
  ylim(0,10000)
dev.off()


pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p35_CB_nFeature_RNA-1stepIntegrationRegressNotRepeated-QCV3dim50kparam15res04.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=FALSE, features = "nFeature_RNA")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p35_CB_percentmt-1stepIntegrationRegressNotRepeated-QCV3dim50kparam15res04.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.mt")
dev.off()  
pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_p35_CB_percentrb-1stepIntegrationRegressNotRepeated-QCV3dim50kparam15res04.pdf", width=5, height=5)
FeaturePlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=FALSE, features = "percent.rb")
dev.off()  
############################################################









pdf("output/seurat/UMAP_WT_Kcnc1_splitCondition-1stepIntegrationRegressNotRepeated-QCV2dim13kparam30res05.pdf", width=13, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "condition")
dev.off()
pdf("output/seurat/UMAP_WT_Kcnc1_splitReplicate-1stepIntegrationRegressNotRepeated-QCV2dim13kparam30res05.pdf", width=15, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, reduction = "umap", label=TRUE, split.by = "replicate")
dev.off()

pdf("output/seurat/FeaturePlot_QCmetrics_WT_Kcnc1_Phase-1stepIntegrationRegressNotRepeated-QCV2dim13kparam30res05.pdf", width=10, height=6)
DimPlot(WT_Kcnc1_p35_CB_1step.sct, group.by= "Phase") & 
  theme(plot.title = element_text(size=10))
dev.off()  


# save ##################
## saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step.sct_V1_numeric.rds") # regMtRbCount with QC_V3
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step.sct_V1_numeric.rds")
## saveRDS(WT_Kcnc1_p35_CB_1step.sct, file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam10res03.sct_V1_numeric.rds") # regMtRbCount with QC_V3



set.seed(42)
##########


## Let's work with 1step integration: 
# --> 1stepIntegrationRegressNotRepeatedregMtRbCou-QCV3dim50kparam10res03 = WT_Kcnc1_p35_CB_1step_V1
WT_Kcnc1_p35_CB_1step.sct <- readRDS(file = "output/seurat/WT_Kcnc1_p35_CB_1step-QCV3dim50kparam10res03.sct_V1_numeric.rds")


############################ EasyCellType automatic annotation ##########################################

### Find all markers 
all_markers <- FindAllMarkers(WT_Kcnc1_p35_CB_1step.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(all_markers, file = "output/seurat/srat_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam10res03_all_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)


# BiocManager::install("EasyCellType")
library("EasyCellType")
library("org.Mm.eg.db")
library("AnnotationDbi")

## load marker
all_markers <- read.delim("output/seurat/srat_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam10res03_all_markers.txt", header = TRUE, row.names = 1)



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

pdf("output/seurat/EasyCellType_dotplot_WT_Kcnc1_p35_CB_1step_QCV3dim50kparam10res03-clustermole_Brain.pdf", width=6, height=8)
plot_dot(test="GSEA", annot.GSEA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




# panglao --> NOT working... 'subscript out of bounds' error... I tried gene as ENSEMBL, entrezID and geneSymbo, human/mic, everything...


#




## check some genes


pdf("output/Signac/FeaturePlot_SCT_RNA_WT_Bap1KO-EpendymalCell-QCV3_dim40kparam35res05algo4feat2000_noCellCycleRegression
.pdf", width=15, height=15)
FeaturePlot(multiome_WT_Bap1KO_QCV2.sct, features = c(  "Rabl2", "Cfap54", "Ccdc153", "Foxj1", "Pifo", "Dynlrb2", "Rsph1", "Cfap44"), max.cutoff = 1, cols = c("grey", "red"))
dev.off()













```

--> How many dims to use? Not clear, using Elbow and [quantification](https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html) say 13dims (but look very few!!). We do not observe significant changes of clustering by changes the numb of dims. Only small cluster are affected (Serotonergic neurons)





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




```



--> https://roulethomas.shinyapps.io/shinyapp_rna_qcv2_clusterv2/ 

