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
sbatch --dependency=afterany:20250949 scripts/scrublet_Kcnc1_p180.sh # 20255974 Kcnc1_p180_CB_Rep3 fail in predicted but it output file; RAS ok


```
Doublet detection score (YAP1scRNAseq ranged between 0.1 to 42%): now range from 0 to 30% (median ~10%)
--> Table in `output/doublets/doublets_all.xlsx`


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


# integration 

XXXX bnelow nto mod



Then `conda activate scRNAseq` and go into R and filtered out **RNA contamination and start with SEURAT**.

Tuto seurat [here](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html)


```bash
conda activate scRNAseqV2
```

Let's 1st try to **perform clustering that look like the E775 time-point** (use same marker); display some genes and generate a **shiny app for Conchi** --> Shared to Conchi and discuss; then see if we put both exp together



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
load("output/soupX/out_E7mousecontrol.RData")
srat_WT_E7 <- CreateSeuratObject(counts = out, project = "WT_E7") # 32,285 features across 1,829 samples

load("output/soupX/out_E7mousecYAPKO.RData")
srat_cYAPKO_E7 <- CreateSeuratObject(counts = out, project = "cYAPKO_E7") # 32,285 features across 1,897 samples


# QUALITY CONTROL
## add mitochondrial and Ribosomal conta 
srat_WT_E7[["percent.mt"]] <- PercentageFeatureSet(srat_WT_E7, pattern = "^mt-")
srat_WT_E7[["percent.rb"]] <- PercentageFeatureSet(srat_WT_E7, pattern = "^Rp[sl]")

srat_cYAPKO_E7[["percent.mt"]] <- PercentageFeatureSet(srat_cYAPKO_E7, pattern = "^mt-")
srat_cYAPKO_E7[["percent.rb"]] <- PercentageFeatureSet(srat_cYAPKO_E7, pattern = "^Rp[sl]")

## add doublet information (scrublet)
doublets <- read.table("output/doublets/embryo_E7_control.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_WT_E7 <- AddMetaData(srat_WT_E7,doublets)
srat_WT_E7$Doublet_score <- as.numeric(srat_WT_E7$Doublet_score) # make score as numeric
head(srat_WT_E7[[]])

doublets <- read.table("output/doublets/embryo_E7_cYAPKO.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat_cYAPKO_E7 <- AddMetaData(srat_cYAPKO_E7,doublets)
srat_cYAPKO_E7$Doublet_score <- as.numeric(srat_cYAPKO_E7$Doublet_score) # make score as numeric
head(srat_cYAPKO_E7[[]])



## Plot
pdf("output/seurat/VlnPlot_QC_embryo_control_E7.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QC_embryo_cYAPKO_E7.pdf", width=10, height=6)
VlnPlot(srat_cYAPKO_E7, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAcount_mt_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAcount_RNAfeature_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
pdf("output/seurat/FeatureScatter_QC_rb_mt_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_rb_mt_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "percent.rb", feature2 = "percent.mt")
dev.off()
pdf("output/seurat/FeatureScatter_QC_mt_doublet_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_mt_doublet_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "percent.mt", feature2 = "Doublet_score")
dev.off()
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_cYAPKO_E7.pdf", width=5, height=5)
pdf("output/seurat/FeatureScatter_QC_RNAfeature_doublet_control_E7.pdf", width=5, height=5)
FeatureScatter(srat_cYAPKO_E7, feature1 = "nFeature_RNA", feature2 = "Doublet_score")
dev.off()




## After seeing the plot; add QC information in our seurat object
### V1 not super stringeant
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 500 & srat_WT_E7@meta.data$QC == 'Pass','Low_nFeature',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 500 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$percent.mt > 10 & srat_WT_E7@meta.data$QC == 'Pass','High_MT',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 500 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
table(srat_WT_E7[['QC']])
## 
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 500 & srat_cYAPKO_E7@meta.data$QC == 'Pass','Low_nFeature',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 500 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$percent.mt > 10 & srat_cYAPKO_E7@meta.data$QC == 'Pass','High_MT',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 500 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
table(srat_cYAPKO_E7[['QC']])


### V2 more stringeant
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC == 'Pass','Low_nFeature',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$percent.mt > 10 & srat_WT_E7@meta.data$QC == 'Pass','High_MT',srat_WT_E7@meta.data$QC)
srat_WT_E7[['QC']] <- ifelse(srat_WT_E7@meta.data$nFeature_RNA < 5000 & srat_WT_E7@meta.data$QC != 'Pass' & srat_WT_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_WT_E7@meta.data$QC,sep = ','),srat_WT_E7@meta.data$QC)
table(srat_WT_E7[['QC']])
## 
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$Is_doublet == 'True','Doublet','Pass')
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC == 'Pass','Low_nFeature',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$percent.mt > 10 & srat_cYAPKO_E7@meta.data$QC == 'Pass','High_MT',srat_cYAPKO_E7@meta.data$QC)
srat_cYAPKO_E7[['QC']] <- ifelse(srat_cYAPKO_E7@meta.data$nFeature_RNA < 5000 & srat_cYAPKO_E7@meta.data$QC != 'Pass' & srat_cYAPKO_E7@meta.data$QC != 'High_MT',paste('High_MT',srat_cYAPKO_E7@meta.data$QC,sep = ','),srat_cYAPKO_E7@meta.data$QC)
table(srat_cYAPKO_E7[['QC']])






# Quality check after QC filtering
pdf("output/seurat/VlnPlot_QCPass_embryo_control_E7.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QCPass_embryo_cYAPKO_E7.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QCPassV2_embryo_control_E7.pdf", width=10, height=6)
pdf("output/seurat/VlnPlot_QCPassV2_embryo_cYAPKO_E7.pdf", width=10, height=6)

VlnPlot(subset(srat_cYAPKO_E7, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()







## subset my seurat object to only analyze the cells that pass the QC
srat_WT_E7 <- subset(srat_WT_E7, subset = QC == 'Pass')
srat_cYAPKO_E7 <- subset(srat_cYAPKO_E7, subset = QC == 'Pass')
srat_WT_E7$condition <- "WT_E7"
srat_cYAPKO_E7$condition <- "cYAPKO_E7"



mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name


## NORMALIZE AND SCALE DATA BEFORE RUNNING CELLCYCLESORTING
srat_WT_E7 <- NormalizeData(srat_WT_E7, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_WT_E7)
srat_WT_E7 <- ScaleData(srat_WT_E7, features = all.genes) # zero-centres and scales it

srat_cYAPKO_E7 <- NormalizeData(srat_cYAPKO_E7, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(srat_cYAPKO_E7)
srat_cYAPKO_E7 <- ScaleData(srat_cYAPKO_E7, features = all.genes) # zero-centres and scales it

### CELLCYCLESORTING
srat_WT_E7 <- CellCycleScoring(srat_WT_E7, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_WT_E7[[]]$Phase)
srat_cYAPKO_E7 <- CellCycleScoring(srat_cYAPKO_E7, s.features = mmus_s, g2m.features = mmus_g2m)
table(srat_cYAPKO_E7[[]]$Phase)

set.seed(42)

# elbow
## srat_WT_E7_elbow = SCTransform(srat_WT_E7, method = "glmGamPoi", ncells = 1624, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","S.Score","G2M.Score"), verbose = TRUE, variable.features.n = 3000) %>% RunPCA(npcs = 50, verbose = FALSE)


pdf("output/seurat/Elbow_srat_WT_E7QCV2.pdf", width=10, height=10)
ElbowPlot(srat_WT_E7_elbow) # 10 or 15 or 19
dev.off()

```

