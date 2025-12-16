# Objectives

Next, we'll consider tissue expression: Perform simulation on specific celltypes/tissue -spe expressed genes:
- Isolate tissue/cell types -specifically expressed genes
- generate FASTA CDS
- Mutate these genes only! (take coding and non-coding!)




# Data download

Let's use the [Tabula Sapiens study](https://www.biorxiv.org/content/10.1101/2024.12.03.626516v1)

Tabula sapiens [paper](https://www.biorxiv.org/node/4837559.external-links.html); then [Chan Zuckerberg CELLxGENE Discover collection](https://tabula-sapiens.sf.czbiohub.org/whereisthedata); and then data for all cells downloaded [here](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) (42Gb)



--> Data from **cellxgene** in `.h5ad` format; follow `002*/009*` for data download and processing


```bash
# Data download

wget https://datasets.cellxgene.cziscience.com/5a495302-b7cd-4bf9-853e-95627b00bb03.h5ad -O input/TabulaSapiens.h5ad 
```


## Install BPCells

We need to install BPCells to load very large `.h5ad` file

```bash
# clone env with seuratv5 and ability to read .h5a file already
conda create --name condiments_Signac_V2 --clone condiments_Signac

# install BPCells
conda activate condiments_Signac_V2
#--> module load hdf5; does not do shit so re-install it

conda install -c anaconda hdf5




```


```R
# install BPCells
remotes::install_github("bnprks/BPCells/r")
#--> Fail with hdf5... So try module load hd5--> FAIL; so try install it through conda




```





Let's try install seurat first, and then BPCells; follow [seurat website guidelines](https://satijalab.org/seurat/articles/install.html)


```bash
conda create -n BPCells -c conda-forge -c bioconda \
  r-base=4.5.1


conda activate BPCells


conda install r::r-igraph 


```




```R
install.packages("remotes")


remotes::install_github("satijalab/seurat", "seurat5", quiet = FALSE)
#--> fail for igraph and leiden

install.packages(c('igraph','leiden'))
#--> fail, test installing with conda:  conda install r::r-igraph 




```

--> FAIL


Let's install it through [conda](https://github.com/bnprks/BPCells/issues/241)


```bash
conda create -n bpcells_r44 -c conda-forge -c bioconda \
  r-base=4.4 \
  r-bpcells \
  bioconductor-genomicranges \
  bioconductor-iranges \
  r-igraph \
  r-matrixstats \
  r-rspectra \
  macs3


conda install conda-forge::pandas # to avoid leidenbase error

conda activate bpcells_r44

conda env remove --name bpcells_r44
```


```R
install.packages('Seurat')
#--> Fail, `ERROR: compilation failed for package leidenbase` --> Seems need to install pandas

```


--> Here was able to install BPCells, but not seurat... Let;'s try installing it through conda too



```bash

conda create -n bpcells_r44 -c conda-forge -c bioconda \
  r-base=4.4 \
  r-bpcells \
  bioconductor-genomicranges \
  bioconductor-iranges \
  r-igraph \
  r-matrixstats \
  r-rspectra \
  macs3 \
  r-seurat


conda install bioconda::r-azimuth # R package needed!
conda install bioconda::r-biomartr # 

conda activate bpcells_r44

```

```R
library("Seurat")
library("BPCells")
library("Azimuth")

```
--> WORK!






- ERROR: When importing `.h5ad` file `Error in py_to_r_cpp(x) : negative length vectors are not allowed`; seems to be a memory error according to [this](https://github.com/satijalab/seurat/issues/203).
    --> Seems I can "use BPCells to load the matrix and create a seurat v5 object" according to [this discussion](https://github.com/satijalab/seurat/issues/7283)
        --> Need to install [BPCells](https://github.com/bnprks/BPCells)
            --> Installed BPCells with anaconda, then install seurat with install.apcakges(Seurat); and got another error related to ledeinbase;    
                --> leidenbase ERROR solve with installing pandas according to [this](https://github.com/satijalab/seurat/issues/3851)











## Convert .hd5a to seurat object, process data and save



```bash
conda activate bpcells_r44 # For seurat v5 with BPCells
```



```R
library("Seurat")
library("BPCells")
library("Azimuth")
library("biomaRt")

# ---- paths ----
file.dir <- "input/"
files.set <- c("TabulaSapiens.h5ad")

# Loop through h5ad files and output BPCells matrices on-disk
data.list <- c()
metadata.list <- c()

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i])
  data <- open_matrix_anndata_hdf5(path)
   write_matrix_dir(
     mat = data,
     dir = paste0(gsub(".h5ad", "", path), "_BP")
   )
  # Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(gsub(".h5ad", "", path), "_BP"))
  mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
  # Get metadata
  metadata.list[[i]] <- LoadH5ADobs(path = path)
  data.list[[i]] <- mat
}
# Name layers
names(data.list) <- c("TabulaSapiens")

# (optional) trim metadata columns before CreateSeuratObject
metadata <- metadata.list[[1]]

# create Seurat object with on-disk counts
seu <- CreateSeuratObject(counts = data.list, meta.data = metadata)
seu


# 

DefaultAssay(seu) <- "RNA"

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, nfeatures = 2000)

seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30)
seu <- RunUMAP(seu, dims = 1:30)

# relevant variables
pdf("output/seurat/TabulaSapiens-30dim-tissue.pdf", width=50, height=20)   
DimPlot(seu, reduction = "umap", group.by = "tissue", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-tissue_in_publication.pdf", width=20, height=10)   
DimPlot(seu, reduction = "umap", group.by = "tissue_in_publication", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-cell_type.pdf", width=100, height=20)   
DimPlot(seu, reduction = "umap", group.by = "cell_type", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-compartment.pdf", width=40, height=20)   
DimPlot(seu, reduction = "umap", group.by = "compartment", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-broad_cell_class.pdf", width=40, height=20)   
DimPlot(seu, reduction = "umap", group.by = "broad_cell_class", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-anatomical_position.pdf", width=40, height=20)   
DimPlot(seu, reduction = "umap", group.by = "anatomical_position", raster = TRUE)
dev.off()


## useless variables
pdf("output/seurat/TabulaSapiens-30dim-tissue_type.pdf", width=20, height=10)   
DimPlot(seu, reduction = "umap", group.by = "tissue_type", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-cell_type.pdf", width=100, height=20)   
DimPlot(seu, reduction = "umap", group.by = "cell_type", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-manually_annotated.pdf", width=40, height=20)   
DimPlot(seu, reduction = "umap", group.by = "manually_annotated", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-free_annotation.pdf", width=40, height=20)   
DimPlot(seu, reduction = "umap", group.by = "free_annotation", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-tissue_ontology_term_id.pdf", width=40, height=20)   
DimPlot(seu, reduction = "umap", group.by = "tissue_ontology_term_id", raster = TRUE)
dev.off()
pdf("output/seurat/TabulaSapiens-30dim-development_stage.pdf", width=40, height=20)   
DimPlot(seu, reduction = "umap", group.by = "development_stage", raster = TRUE)
dev.off()




# SAVE OUTPUT
#saveRDS(seu, file = "output/seurat/TabulaSapiens_allCells-dim30.rds")

```
--> .hd5a file succesfully converted to seurat





## Identify top marker genes in tissue and cell types


To isolate expressed/detected genes we will generate several list:
- **FindAllMarkersmipct025logfc025**: Here isolate gene expressed in at least 25% of the cell of the cluster, and globally more express in this cluster as compared to all over cluster (use `FindAllMarkers()` with `min.pct = 0.25, logfc.threshold = 0.25`)
- **MatrixFilter025**: Here we simply keep all genes express/detected in at least 25% of the cell of the cluster (use custom script)




```bash
conda activate bpcells_r44 # For seurat v5 with BPCells
```



```R
library("Seurat")
library("BPCells")
library("Azimuth")
library("biomaRt")
library("Matrix")
library("ggplot2")


set.seed(42)

# load rds object

TabulaSapiens_allCells <- readRDS(file = "output/seurat/TabulaSapiens_allCells-dim30.rds") # 

####################################
# tissue resolution ################
####################################

# FindAllMarkersmipct025logfc025  ###########################
## Unbiased cell type marker genes
Idents(TabulaSapiens_allCells) <- "tissue"
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(TabulaSapiens_allCells) <- "RNA"
TabulaSapiens_allCells <- NormalizeData(TabulaSapiens_allCells, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(TabulaSapiens_allCells)
TabulaSapiens_allCells <- ScaleData(TabulaSapiens_allCells, features = all.genes) # zero-centres and scales it

all_markers <- FindAllMarkers(TabulaSapiens_allCells, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_TabulaSapiens_allCells-dim30-FindAllMarkersmipct025logfc025-tissue.txt", sep = "\t", quote = FALSE, row.names = TRUE)



# FindAllMarkersmipct050logfc025  ###########################
## Unbiased cell type marker genes
Idents(TabulaSapiens_allCells) <- "tissue"
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(TabulaSapiens_allCells) <- "RNA"
TabulaSapiens_allCells <- NormalizeData(TabulaSapiens_allCells, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(TabulaSapiens_allCells)
TabulaSapiens_allCells <- ScaleData(TabulaSapiens_allCells, features = all.genes) # zero-centres and scales it

all_markers <- FindAllMarkers(TabulaSapiens_allCells, assay = "RNA", only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_TabulaSapiens_allCells-dim30-FindAllMarkersmipct050slogfc025-tissue.txt", sep = "\t", quote = FALSE, row.names = TRUE)






# MatrixFilter0XX  ###########################
Idents(TabulaSapiens_allCells) <- "tissue"  # replace "celltype" with your metadata column name
## Threshold for expression detection
threshold <- 0.75  # 25%
## Create output directory


## Get all cell types
celltypes <- levels(Idents(TabulaSapiens_allCells))

### Loop through each cell type
for (ct in celltypes) {
  cat("Processing:", ct, "\n")
  # Subset cells from this cell type
  cells <- WhichCells(TabulaSapiens_allCells, idents = ct)
  # Compute detection proportion (>0 counts)
  expr_matrix <- GetAssayData(TabulaSapiens_allCells, slot = "counts")[, cells]
  pct_expressed <- Matrix::rowMeans(expr_matrix > 0)
  # Keep genes expressed in at least 25% of cells
  expressed_genes <- names(pct_expressed[pct_expressed >= threshold])
  # Save one file per cell type
  outfile <- paste0("output/seurat/", ct, "-MatrixFilter075.txt")
  writeLines(expressed_genes, outfile)
}







####################################
# cell_type resolution ################
####################################



# RUN IN SLURM  #######################################
# FindAllMarkersmipct025logfc025  ###########################
## Unbiased cell type marker genes
Idents(TabulaSapiens_allCells) <- "cell_type"
## PRIOR Lets switch to RNA assay and normalize and scale before doing the DEGs
DefaultAssay(TabulaSapiens_allCells) <- "RNA"
TabulaSapiens_allCells <- NormalizeData(TabulaSapiens_allCells, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(TabulaSapiens_allCells)
TabulaSapiens_allCells <- ScaleData(TabulaSapiens_allCells, features = all.genes) # zero-centres and scales it

all_markers <- FindAllMarkers(TabulaSapiens_allCells, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all_markers, file = "output/seurat/srat_TabulaSapiens_allCells-dim30-FindAllMarkersmipct025logfc025-cell_type.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#################################################################



# MatrixFilter0XX  ###########################
Idents(TabulaSapiens_allCells) <- "cell_type"  # replace "celltype" with your metadata column name
## Threshold for expression detection
threshold <- 0.75  # 25%
## Create output directory


## Get all cell types
celltypes <- levels(Idents(TabulaSapiens_allCells))

### Loop through each cell type
for (ct in celltypes) {
  cat("Processing:", ct, "\n")
  # Subset cells from this cell type
  cells <- WhichCells(TabulaSapiens_allCells, idents = ct)
  # Compute detection proportion (>0 counts)
  expr_matrix <- GetAssayData(TabulaSapiens_allCells, slot = "counts")[, cells]
  pct_expressed <- Matrix::rowMeans(expr_matrix > 0)
  # Keep genes expressed in at least 25% of cells
  expressed_genes <- names(pct_expressed[pct_expressed >= threshold])
  # Save one file per cell type
  outfile <- paste0("output/seurat/", ct, "-MatrixFilter075.txt")
  writeLines(expressed_genes, outfile)
}














####################################
# Random extra analysis ################
####################################


# tissue ###############
Idents(TabulaSapiens_allCells) <- "tissue" 



cells.kidney <- WhichCells(TabulaSapiens_allCells, expression = grepl("kidney", tissue, ignore.case=TRUE))
cells.heart  <- WhichCells(TabulaSapiens_allCells, expression = grepl("heart",  tissue, ignore.case=TRUE))
cells.stomach    <- WhichCells(TabulaSapiens_allCells, expression = grepl("stomach",    tissue, ignore.case=TRUE))

pdf("output/seurat/TabulaSapiens-30dim-kidneyHeartStomach.pdf", width=6, height=5)
DimPlot(
  TabulaSapiens_allCells,
  reduction       = "umap",
  cols            = "grey90",
  cells.highlight = list(kidney= cells.kidney , heart= cells.heart , stomach= cells.stomach),
  cols.highlight  = c("#1f77b4", "#2ca02c", "#d62728"),
  sizes.highlight = 0.7, shuffle = TRUE
) 
dev.off()



cells.lung <- WhichCells(TabulaSapiens_allCells, expression = grepl("lung", tissue, ignore.case=TRUE))
cells.blood  <- WhichCells(TabulaSapiens_allCells, expression = grepl("blood",  tissue, ignore.case=TRUE))
cells.spleen    <- WhichCells(TabulaSapiens_allCells, expression = grepl("spleen",    tissue, ignore.case=TRUE))

pdf("output/seurat/TabulaSapiens-30dim-lungBloddSpleen.pdf", width=6, height=5)
DimPlot(
  TabulaSapiens_allCells,
  reduction       = "umap",
  cols            = "grey90",
  cells.highlight = list(lung= cells.lung, blood= cells.blood, spleen= cells.spleen),
  cols.highlight  = c("#1f77b4", "#d62728", "#2ca02c"),
  sizes.highlight = 0.7, shuffle = TRUE
) 
dev.off()



# cell_type ###############



Idents(TabulaSapiens_allCells) <- "cell_type" 


cells.kidneyepicell <- WhichCells(TabulaSapiens_allCells, expression = grepl("kidney epithelial cell", cell_type, ignore.case=TRUE))
cells.regularatricardmyoc  <- WhichCells(TabulaSapiens_allCells, expression = grepl("regular atrial cardiac myocyte",  cell_type, ignore.case=TRUE))
cells.tonguemuscl    <- WhichCells(TabulaSapiens_allCells, expression = grepl("tongue muscle cell",    cell_type, ignore.case=TRUE))

pdf("output/seurat/TabulaSapiens-30dim-KidneyepicellRegularatricardmyocTonguemuscl.pdf", width=6, height=5)
DimPlot(
  TabulaSapiens_allCells,
  reduction       = "umap",
  cols            = "grey90",
  cells.highlight = list(kidneyepicell= cells.kidneyepicell, regularatricardmyoc= cells.regularatricardmyoc, tonguemuscl= cells.tonguemuscl),
  cols.highlight  = c("#1f77b4", "#2ca02c", "#d62728"),
  sizes.highlight = 0.7, shuffle = TRUE
) 
dev.off()




cells.platelet <- WhichCells(TabulaSapiens_allCells, expression = grepl("platelet", cell_type, ignore.case=TRUE))
cells.basalcell  <- WhichCells(TabulaSapiens_allCells, expression = grepl("basal cell",  cell_type, ignore.case=TRUE))
cells.Tcell    <- WhichCells(TabulaSapiens_allCells, expression = grepl("T cell",    cell_type, ignore.case=TRUE))

pdf("output/seurat/TabulaSapiens-30dim-plateletBasalcellTcell.pdf", width=6, height=5)
DimPlot(
  TabulaSapiens_allCells,
  reduction       = "umap",
  cols            = "grey90",
  cells.highlight = list(platelet= cells.platelet, basalcell= cells.basalcell, Tcell= cells.Tcell),
  cols.highlight  = c("#1f77b4", "#d62728", "#2ca02c"),
  sizes.highlight = 0.7, shuffle = TRUE
) 
dev.off()







# Gene count per cell_type or tissue
Idents(TabulaSapiens_allCells) <- "tissue"
threshold <- 0.75  # expressed in at least 75% of cells in that cell type
celltypes <- levels(Idents(TabulaSapiens_allCells))
gene_counts <- data.frame(
  cell_type = character(),
  n_genes = numeric(),
  stringsAsFactors = FALSE
)
for (ct in celltypes) {
  cat("Processing:", ct, "\n")
  cells <- WhichCells(TabulaSapiens_allCells, idents = ct)
  expr_matrix <- GetAssayData(TabulaSapiens_allCells, slot = "counts")[, cells]
  pct_expressed <- Matrix::rowMeans(expr_matrix > 0)
  expressed_genes <- names(pct_expressed[pct_expressed >= threshold])
  # record how many genes passed
  gene_counts <- rbind(gene_counts, data.frame(cell_type = ct, n_genes = length(expressed_genes)))
}
# Compute median number of genes across all cell types
median_genes <- median(gene_counts$n_genes)
cat("\n✅ Median number of expressed genes across all cell types:", median_genes, "\n")
write.table(gene_counts, "output/seurat/gene_counts-cell_type-MatrixFilter075-tissue.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

```




- NOTE: `MatrixFilter025` adding 'count' filtering is not working great.. Even using `expr_matrix > 0` lead to very low number of genes...

--> `FindAllMarkers()` is better to isolate specifically expressed genes



Unbiased marker gene identification with find `FindAllMarkers()` ran in slurm below for cell_type:

```bash
conda activate bpcells_r44 # For seurat v5 with BPCells

sbatch scripts/FindAllMarkers_TabulaSapiens_allCells-FindAllMarkersmipct025logfc025.sh # 59562106 ok
sbatch scripts/FindAllMarkers_TabulaSapiens_allCells-FindAllMarkersmipct050logfc025.sh # 61873932 ok

```






## Isolate highly and specifically expressed genes in each cluster


Let's try to keep only genes marker of ONE specific cluster (ie. remove any genes identify as marker in more than 2 cluster!)


```bash
conda activate deseq2
```


```R
library("tidyverse")



###################
# > 25% TISSUE ##########
###################

# Keep genes that appear in only ONE tissue (ie. remove any genes found as marker in another tissue)
all_markers <- read.delim(
  "output/seurat/srat_TabulaSapiens_allCells-dim30-FindAllMarkersmipct025logfc025-tissue.txt",
  row.names = 1
) %>% as_tibble()


## Count number of highly and specficcaly expressed genes per cluster
all_markers %>%
  filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  summarise(exclusive_genes = list(sort(unique(gene))),
            n = n_distinct(gene)) %>%
  arrange(desc(n)) %>%
  dplyr::select(cluster, n) %>%
  print(n=75) 


## Count number of specific genes per cluster
specific_markers <- all_markers %>%
  filter(p_val_adj <0.05) %>%
  group_by(gene) %>%
  mutate(n_tissues = n_distinct(cluster)) %>%   # 'cluster' column = tissue name in FindAllMarkers output
  ungroup() %>%
  filter(n_tissues == 1)

specific_markers %>%
  group_by(cluster) %>%
  summarise(exclusive_genes = list(sort(unique(gene))),
            n = n_distinct(gene)) %>%
  arrange(desc(n)) %>%
  dplyr::select(cluster, n) %>%
  print(n=75) 


## Test different FC treshold

fc_thresholds <- c(0.25, 0.5, 1, 2)
all_tissues <- sort(unique(all_markers$cluster))  # cluster column holds tissue names here

marker_counts_by_fc <- tidyr::crossing(
    cluster = all_tissues,
    fc_threshold = fc_thresholds
  ) %>%
  left_join(
    all_markers %>%
      filter(p_val_adj < 0.05) %>%
      tidyr::crossing(fc_threshold = fc_thresholds) %>%
      filter(avg_log2FC >= fc_threshold) %>%
      group_by(cluster, fc_threshold) %>%
      summarise(n_genes = n_distinct(gene), .groups = "drop"),
    by = c("cluster", "fc_threshold")
  ) %>%
  mutate(n_genes = tidyr::replace_na(n_genes, 0),
         fc_threshold = paste0("log2FC_ge_", fc_threshold)) %>%
  pivot_wider(names_from = fc_threshold, values_from = n_genes) %>%
  arrange(cluster)




write.table(marker_counts_by_fc,
            "output/seurat/marker_counts-FindAllMarkersmipct025logfc025-tissue-FC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


###################
# > 50% TISSUE ##########
###################

# Keep genes that appear in only ONE tissue (ie. remove any genes found as marker in another tissue)
all_markers <- read.delim(
  "output/seurat/srat_TabulaSapiens_allCells-dim30-FindAllMarkersmipct050slogfc025-tissue.txt",
  row.names = 1
) %>% as_tibble()


## Count number of highly and specficcaly expressed genes per cluster
all_markers %>%
  filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  summarise(exclusive_genes = list(sort(unique(gene))),
            n = n_distinct(gene)) %>%
  arrange(desc(n)) %>%
  dplyr::select(cluster, n) %>%
  print(n=75) 


## Count number of specific genes per cluster
specific_markers <- all_markers %>%
  filter(p_val_adj <0.05) %>%
  group_by(gene) %>%
  mutate(n_tissues = n_distinct(cluster)) %>%   # 'cluster' column = tissue name in FindAllMarkers output
  ungroup() %>%
  filter(n_tissues == 1)

specific_markers %>%
  group_by(cluster) %>%
  summarise(exclusive_genes = list(sort(unique(gene))),
            n = n_distinct(gene)) %>%
  arrange(desc(n)) %>%
  dplyr::select(cluster, n) %>%
  print(n=75) 


## Test different FC treshold

fc_thresholds <- c(0.25, 0.5, 1, 2)
all_tissues <- sort(unique(all_markers$cluster))  # cluster column holds tissue names here

marker_counts_by_fc <- tidyr::crossing(
    cluster = all_tissues,
    fc_threshold = fc_thresholds
  ) %>%
  left_join(
    all_markers %>%
      filter(p_val_adj < 0.05) %>%
      tidyr::crossing(fc_threshold = fc_thresholds) %>%
      filter(avg_log2FC >= fc_threshold) %>%
      group_by(cluster, fc_threshold) %>%
      summarise(n_genes = n_distinct(gene), .groups = "drop"),
    by = c("cluster", "fc_threshold")
  ) %>%
  mutate(n_genes = tidyr::replace_na(n_genes, 0),
         fc_threshold = paste0("log2FC_ge_", fc_threshold)) %>%
  pivot_wider(names_from = fc_threshold, values_from = n_genes) %>%
  arrange(cluster)




write.table(marker_counts_by_fc,
            "output/seurat/marker_counts-FindAllMarkersmipct050logfc025-tissue-FC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)






###################
# > 25% CELL TYPE ##########
###################

# Keep genes that appear in only ONE cell_type (ie. remove any genes found as marker in another v)
all_markers <- read.delim(
  "output/seurat/srat_TabulaSapiens_allCells-dim30-FindAllMarkersmipct025logfc025-cell_type.txt",
  row.names = 1
) %>% as_tibble()


## Count number of highly and specficcaly expressed genes per cluster
all_markers %>%
  filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  summarise(exclusive_genes = list(sort(unique(gene))),
            n = n_distinct(gene)) %>%
  arrange(desc(n)) %>%
  dplyr::select(cluster, n) %>%
  print(n=180) 


## Count number of specific genes per cluster
specific_markers <- all_markers %>%
  filter(p_val_adj <0.05) %>%
  group_by(gene) %>%
  mutate(n_tissues = n_distinct(cluster)) %>%   # 'cluster' column = cell_type name in FindAllMarkers output
  ungroup() %>%
  filter(n_tissues == 1)

specific_markers %>%
  group_by(cluster) %>%
  summarise(exclusive_genes = list(sort(unique(gene))),
            n = n_distinct(gene)) %>%
  arrange(desc(n)) %>%
  dplyr::select(cluster, n) %>%
  print(n=180) 


## Test different FC treshold

fc_thresholds <- c(0.25, 0.5, 1, 2)
all_tissues <- sort(unique(all_markers$cluster))  # cluster column holds cell_type names here

marker_counts_by_fc <- tidyr::crossing(
    cluster = all_tissues,
    fc_threshold = fc_thresholds
  ) %>%
  left_join(
    all_markers %>%
      filter(p_val_adj < 0.05) %>%
      tidyr::crossing(fc_threshold = fc_thresholds) %>%
      filter(avg_log2FC >= fc_threshold) %>%
      group_by(cluster, fc_threshold) %>%
      summarise(n_genes = n_distinct(gene), .groups = "drop"),
    by = c("cluster", "fc_threshold")
  ) %>%
  mutate(n_genes = tidyr::replace_na(n_genes, 0),
         fc_threshold = paste0("log2FC_ge_", fc_threshold)) %>%
  pivot_wider(names_from = fc_threshold, values_from = n_genes) %>%
  arrange(cluster)




write.table(marker_counts_by_fc,
            "output/seurat/marker_counts-FindAllMarkersmipct025logfc025-cell_type-FC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


###################
# > 50% CELL TYPE ##########
###################


# Keep genes that appear in only ONE cell_type (ie. remove any genes found as marker in another v)
all_markers <- read.delim(
  "output/seurat/srat_TabulaSapiens_allCells-dim30-FindAllMarkersmipct050logfc025-cell_type.txt",
  row.names = 1
) %>% as_tibble()


## Count number of highly and specficcaly expressed genes per cluster
all_markers %>%
  filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  summarise(exclusive_genes = list(sort(unique(gene))),
            n = n_distinct(gene)) %>%
  arrange(desc(n)) %>%
  dplyr::select(cluster, n) %>%
  print(n=180) 


## Count number of specific genes per cluster
specific_markers <- all_markers %>%
  filter(p_val_adj <0.05) %>%
  group_by(gene) %>%
  mutate(n_tissues = n_distinct(cluster)) %>%   # 'cluster' column = cell_type name in FindAllMarkers output
  ungroup() %>%
  filter(n_tissues == 1)

specific_markers %>%
  group_by(cluster) %>%
  summarise(exclusive_genes = list(sort(unique(gene))),
            n = n_distinct(gene)) %>%
  arrange(desc(n)) %>%
  dplyr::select(cluster, n) %>%
  print(n=180) 


## Test different FC treshold

fc_thresholds <- c(0.25, 0.5, 1, 2)
all_tissues <- sort(unique(all_markers$cluster))  # cluster column holds cell_type names here

marker_counts_by_fc <- tidyr::crossing(
    cluster = all_tissues,
    fc_threshold = fc_thresholds
  ) %>%
  left_join(
    all_markers %>%
      filter(p_val_adj < 0.05) %>%
      tidyr::crossing(fc_threshold = fc_thresholds) %>%
      filter(avg_log2FC >= fc_threshold) %>%
      group_by(cluster, fc_threshold) %>%
      summarise(n_genes = n_distinct(gene), .groups = "drop"),
    by = c("cluster", "fc_threshold")
  ) %>%
  mutate(n_genes = tidyr::replace_na(n_genes, 0),
         fc_threshold = paste0("log2FC_ge_", fc_threshold)) %>%
  pivot_wider(names_from = fc_threshold, values_from = n_genes) %>%
  arrange(cluster)




write.table(marker_counts_by_fc,
            "output/seurat/marker_counts-FindAllMarkersmipct050logfc025-cell_type-FC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


```




