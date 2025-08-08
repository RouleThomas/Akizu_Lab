# Project

Jasmine making a review about PRC2 subunits and brain stuff. 
- At revision stage they would like to add scRNAseq information data
    - Compare level of expression of PRC2 subunits between drosphila and human
    - Strenghten that many more types of neurons in humans than in drosophila


# Pipeline

Download and re-analize public scRNAseq data from fly brain, and human brain neurons:
- Public data in [SCOPE](https://scope.aertslab.org/#/Fly_Brain/Fly_Brain%2FRavenscroft_et_al_2019_LarvalBrain.loom/welcome):
    - Davie_Janssens_Koldere_et_al_2018_AdultBrain --> Weird naming of cell types **FROM  Janssens et al., 2022 (Nature)**
    - Aerts_Fly_AdultBrain_Filtered_57k --> Fail in importing the data; see error
    - **Zhang_Human_Brain_30k** --> LOOKS GOOD! Clear cell type, cluster well **FROM Davie, Jannssens and Koldere et al., 2018 (Cell)**
    - Waddell_CentralBrain_10k --> No cell type information
    - Ravenscroft_et_al_2019_LarvalBrain --> Weird naming of cell types **FROM  Janssens et al., 2022 (Nature)**
    - Aerts_Fly_AdultBrain_Unfiltered_157k -->  No cell type information
    - **Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain** --> Downloaded from [here](https://flybrain.aertslab.org/?_inputs_&page=%22Downloads%22) --> LOOKS GOOD!!!
    - 5a3ba000_20180809__Davie_Janssens_Koldere_et_al_2018_Adult_Brain_HARMONY_SCENIC **FROM FlyCellAtlas [here](https://scope.aertslab.org/#/FlyCellAtlas/*/welcome)** -->  No cell type information
    - 5a3ba000_20200313__Davie_Janssens_Koldere_et_al_2018_Adult_Brain_BBKNN **FROM FlyCellAtlas [here](https://scope.aertslab.org/#/FlyCellAtlas/*/welcome)** -->  No cell type information
    - Aerts_Supplemental_Fly_AdultBrain_DropSeq_2k
    - **Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain** --> Downloaded from [here](https://flybrain.aertslab.org/?_inputs_&page=%22Downloads%22) --> LOOKS GOOD!!!
    - **Zhang_Human_Brain_60k_old** --> LOOKS GOOD! Clear cell type, cluster well **FROM Davie, Jannssens and Koldere et al., 2018 (Cell)**
    - GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX: Drosophila data https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107451 **FROM Davie, Jannssens and Koldere et al., 2018 (Cell)**




# Analysis 

## Davie_Janssens_Koldere_et_al_2018_AdultBrain

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")



# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Davie_Janssens_Koldere_et_al_2018_AdultBrain.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")

## See which metadata are included in the loom
names(s_cnct[['col_attrs']])


## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
### Get the cell IDs
n_cells <- s_cnct[['col_attrs']][['CellID']][['dims']]
cellids <- s_cnct[['col_attrs']][['CellID']][1:n_cells]
### Get all available metadata field names
meta_fields <- names(s_cnct[['col_attrs']])
### Exclude CellID (already handled)
meta_fields <- setdiff(meta_fields, "CellID")
### Build a metadata dataframe from all fields
metadata <- data.frame(cellID = cellids)
for (field in meta_fields) {
  # Extract the first n_cells for each metadata field
  metadata[[field]] <- s_cnct[['col_attrs']][[field]][1:n_cells]
}
### Set rownames
rownames(metadata) <- cellids
### get raw counts matrix
raw.cnts <- as.data.frame(t(s_cnct[['matrix']][,] ))
names(raw.cnts) <- cellids
rownames(raw.cnts) <- gns

## Create the Seurat object:
Davie_Janssens_Koldere_et_al_2018_AdultBrain <- CreateSeuratObject(counts = raw.cnts,
                            project = "Davie_Janssens_Koldere_et_al_2018_AdultBrain",
                            assay = "RNA",
                            meta.data = metadata)


### Find metadata containing the cell type
unique(Davie_Janssens_Koldere_et_al_2018_AdultBrain$cell_type)
unique(Davie_Janssens_Koldere_et_al_2018_AdultBrain$ClusterID)
```


--> Cluster names are not clear in here; maybe `cell_type` but not so clear.. Let's prefer another dataset. Or ClusterID contains 58 cluster, but cannot find to which they correspond.. Even from their [paper](https://www.cell.com/cell/fulltext/S0092-8674(18)30720-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418307207%3Fshowall%3Dtrue):
    - Seems complex dataset with many time points to filter.
    --> Let's investigate other loom dataset and come back if other dataset not relevant




## Aerts_Fly_AdultBrain_Filtered_57k

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")



# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Aerts_Fly_AdultBrain_Filtered_57k.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")


## See which metadata are included in the loom
names(s_cnct[['col_attrs']])


# Extract gene names (corresponds to columns of the matrix)
gns <- s_cnct[['row_attrs']][['Gene']][]
# Extract cell IDs (corresponds to rows of the matrix)
cellids <- s_cnct[['col_attrs']][['CellID']][]
# Extract count matrix (cells x genes)
mat <- s_cnct[['matrix']][,]      # dim: cells x genes
# Assign row and column names
rownames(mat) <- cellids
colnames(mat) <- gns
# If you want a data frame (optional, but not recommended for large matrices)
# raw.cnts <- as.data.frame(mat)
# For Seurat, you typically keep as a matrix (or convert to sparse Matrix)
library("Matrix")
raw.cnts <- Matrix(mat, sparse = TRUE)

# ADD METADATA
#--> FAIL!!!

Aerts_Fly_AdultBrain_Filtered_57k <- CreateSeuratObject(counts = raw.cnts, project = "Aerts_Fly_AdultBrain_Filtered_57k")



##############################
# CLUSTERING AND UMAPS ##########
##############################

DefaultAssay(Aerts_Fly_AdultBrain_Filtered_57k) <- "RNA"

## NORMALIZE AND SCALE DATA 
Aerts_Fly_AdultBrain_Filtered_57k <- NormalizeData(Aerts_Fly_AdultBrain_Filtered_57k, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(Aerts_Fly_AdultBrain_Filtered_57k)
Aerts_Fly_AdultBrain_Filtered_57k <- ScaleData(Aerts_Fly_AdultBrain_Filtered_57k, features = all.genes) # zero-centres and scales it


## Find variable features
Aerts_Fly_AdultBrain_Filtered_57k = FindVariableFeatures(Aerts_Fly_AdultBrain_Filtered_57k, nfeatures = 3000)

Aerts_Fly_AdultBrain_Filtered_57k <- RunPCA(Aerts_Fly_AdultBrain_Filtered_57k, verbose = FALSE, npcs = 50)
Aerts_Fly_AdultBrain_Filtered_57k <- RunUMAP(Aerts_Fly_AdultBrain_Filtered_57k, reduction = "pca", dims = 1:50, verbose = FALSE)
Aerts_Fly_AdultBrain_Filtered_57k <- FindNeighbors(Aerts_Fly_AdultBrain_Filtered_57k, reduction = "pca", k.param = 15, dims = 1:50)
Aerts_Fly_AdultBrain_Filtered_57k <- FindClusters(Aerts_Fly_AdultBrain_Filtered_57k, resolution = 0.3, verbose = FALSE, algorithm = 4)



# Plot to confirm it work

pdf("output/seurat/UMAP-Aerts_Fly_AdultBrain_Filtered_57k-dim50kparam15res03-seurat_clusters.pdf", width=10, height=6)
DimPlot(Aerts_Fly_AdultBrain_Filtered_57k, reduction = "umap", label=TRUE, group.by = "seurat_clusters")
dev.off()









```


--> Error cannot import METADATA so work without it.









## Zhang_Human_Brain_30k

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")
set.seed(42)

##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Zhang_Human_Brain_30k.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")

## See which metadata are included in the loom
names(s_cnct[['col_attrs']])

# Extract gene names (corresponds to columns of the matrix)
gns <- s_cnct[['row_attrs']][['Gene']][]
# Extract cell IDs (corresponds to rows of the matrix)
cellids <- s_cnct[['col_attrs']][['CellID']][]
# Extract count matrix (cells x genes)
mat <- s_cnct[['matrix']][,]      # dim: cells x genes
# Assign row and column names
rownames(mat) <- cellids
colnames(mat) <- gns
# If you want a data frame (optional, but not recommended for large matrices)
# raw.cnts <- as.data.frame(mat)
# For Seurat, you typically keep as a matrix (or convert to sparse Matrix)
library("Matrix")
raw.cnts <- Matrix(mat, sparse = TRUE)

# ADD METADATA
# Helper to extract and flatten metadata to 1D vector
get1d <- function(x) {
    arr <- x[]
    if (is.matrix(arr)) arr[1, ] else as.vector(arr)
}
# Get all metadata field names except CellID
meta_fields <- setdiff(names(s_cnct[['col_attrs']]), "CellID")
# Create metadata data frame, first column is cellID
metadata <- data.frame(cellID = cellids, stringsAsFactors = FALSE)
# Add all metadata fields
for (field in meta_fields) {
  vec <- get1d(s_cnct[['col_attrs']][[field]])
  if(length(vec) != length(cellids)) {
    warning(sprintf("Skipping metadata field '%s': length %d != number of cells %d", 
                    field, length(vec), length(cellids)))
    next
  }
  metadata[[field]] <- vec
}
# Set rownames to cell IDs
rownames(metadata) <- cellids
# Check result
str(metadata)

## Create the Seurat object:
Zhang_Human_Brain_30k <- CreateSeuratObject(counts = t(raw.cnts),
                            project = "Zhang_Human_Brain_30k",
                            assay = "RNA",
                            meta.data = metadata)


### Find metadata containing the cell type
unique(Zhang_Human_Brain_30k$annot1)
unique(Zhang_Human_Brain_30k$annot2)
unique(Zhang_Human_Brain_30k$lvl1class)



####################################
# USE UMAP FROM THE AUTHOR - Embeddings coordinates #########
####################################
# Extract original embedding coordinates and annotation
umap_x <- s_cnct[['col_attrs']][['Embeddings_X']][]
umap_y <- s_cnct[['col_attrs']][['Embeddings_Y']][]
annot1 <- s_cnct[['col_attrs']][['annot1']][]
lvl1class <- s_cnct[['col_attrs']][['lvl1class']][]


pdf("output/seurat/Embeddings-Zhang_Human_Brain_30k-annot1.pdf", width=10, height=6)
df <- data.frame(UMAP_1 = umap_x, UMAP_2 = umap_y, annot1 = annot1)
ggplot(df, aes(x = UMAP_1.1, y = UMAP_2.1, color = annot1)) +
  geom_point(size = 0.5, alpha = 0.8) +
  theme_bw() +
  labs(title = "Original Embedding from Loom", color = "annot1")
dev.off()

#--> Embeddings looks bad; lets prefer re-clustering it


##############################
# CLUSTERING AND UMAPS ##########
##############################

DefaultAssay(Zhang_Human_Brain_30k) <- "RNA"

## NORMALIZE AND SCALE DATA 
Zhang_Human_Brain_30k <- NormalizeData(Zhang_Human_Brain_30k, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(Zhang_Human_Brain_30k)
Zhang_Human_Brain_30k <- ScaleData(Zhang_Human_Brain_30k, features = all.genes) # zero-centres and scales it


## Find variable features
Zhang_Human_Brain_30k = FindVariableFeatures(Zhang_Human_Brain_30k, nfeatures = 3000)

Zhang_Human_Brain_30k <- RunPCA(Zhang_Human_Brain_30k, verbose = FALSE, npcs = 25)
Zhang_Human_Brain_30k <- RunUMAP(Zhang_Human_Brain_30k, reduction = "pca", dims = 1:25, verbose = FALSE)
Zhang_Human_Brain_30k <- FindNeighbors(Zhang_Human_Brain_30k, reduction = "pca", k.param = 15, dims = 1:25)
Zhang_Human_Brain_30k <- FindClusters(Zhang_Human_Brain_30k, resolution = 0.3, verbose = FALSE, algorithm = 4)



# Plot to confirm it work

pdf("output/seurat/UMAP-Zhang_Human_Brain_30k-dim25kparam15res03-annot1.pdf", width=10, height=6)
DimPlot(Zhang_Human_Brain_30k, reduction = "umap", label=TRUE, group.by = "annot1")
dev.off()

pdf("output/seurat/UMAP-Zhang_Human_Brain_30k-dim25kparam15res03-lvl1class.pdf", width=10, height=6)
DimPlot(Zhang_Human_Brain_30k, reduction = "umap", label=TRUE, group.by = "lvl1class")
dev.off()



### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(Zhang_Human_Brain_30k, file = "output/seurat/Zhang_Human_Brain_30k-dim25kparam15res03.rds")
Zhang_Human_Brain_30k <- readRDS("output/seurat/Zhang_Human_Brain_30k-dim25kparam15res03.rds")
##########################################################

# Volcano plot of expression

gene_list = c("EZH1", "EZH2", "EED", "SUZ12", "PHF1", "MTF2", "PHF19", "JARID2", "AEBP2", "RBBP4", "RBBP7","LCOR","LCORL") #
#--> Missing genes EPOP

pdf("output/seurat/VlnPlot-Zhang_Human_Brain_30k-annot1-genelist.pdf", width = 10, height = 3)
for (gene in gene_list) {
  print(paste("Generating plot for:", gene))
  p <- VlnPlot(Zhang_Human_Brain_30k, 
               features = gene, 
               group.by = "annot1",  # x-axis is your annotation
               pt.size = 0) +
    theme(plot.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(p)
}
dev.off()



```

--> Confirm with Naiara that cell types and sample is good to use, then change parameters to improve clustering







## Waddell_CentralBrain_10k

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Waddell_CentralBrain_10k.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")


## See which metadata are included in the loom
names(s_cnct[['col_attrs']])



```

--> No cell type information




## Ravenscroft_et_al_2019_LarvalBrain

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Ravenscroft_et_al_2019_LarvalBrain.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")

## See which metadata are included in the loom
names(s_cnct[['col_attrs']])

# Extract gene names (corresponds to columns of the matrix)
gns <- s_cnct[['row_attrs']][['Gene']][]
# Extract cell IDs (corresponds to rows of the matrix)
cellids <- s_cnct[['col_attrs']][['CellID']][]
# Extract count matrix (cells x genes)
mat <- s_cnct[['matrix']][,]      # dim: cells x genes
# Assign row and column names
rownames(mat) <- cellids
colnames(mat) <- gns
# If you want a data frame (optional, but not recommended for large matrices)
# raw.cnts <- as.data.frame(mat)
# For Seurat, you typically keep as a matrix (or convert to sparse Matrix)
library("Matrix")
raw.cnts <- Matrix(mat, sparse = TRUE)

# ADD METADATA
# Helper to extract and flatten metadata to 1D vector
get1d <- function(x) {
    arr <- x[]
    if (is.matrix(arr)) arr[1, ] else as.vector(arr)
}
# Get all metadata field names except CellID
meta_fields <- setdiff(names(s_cnct[['col_attrs']]), "CellID")
# Create metadata data frame, first column is cellID
metadata <- data.frame(cellID = cellids, stringsAsFactors = FALSE)
# Add all metadata fields
for (field in meta_fields) {
  vec <- get1d(s_cnct[['col_attrs']][[field]])
  if(length(vec) != length(cellids)) {
    warning(sprintf("Skipping metadata field '%s': length %d != number of cells %d", 
                    field, length(vec), length(cellids)))
    next
  }
  metadata[[field]] <- vec
}
# Set rownames to cell IDs
rownames(metadata) <- cellids
# Check result
str(metadata)

## Create the Seurat object:
Ravenscroft_et_al_2019_LarvalBrain <- CreateSeuratObject(counts = t(raw.cnts),
                            project = "Ravenscroft_et_al_2019_LarvalBrain",
                            assay = "RNA",
                            meta.data = metadata)


### Find metadata containing the cell type
unique(Ravenscroft_et_al_2019_LarvalBrain$annotation)
unique(Ravenscroft_et_al_2019_LarvalBrain$ClusterID)



```

--> Weird naming








## Aerts_Fly_AdultBrain_Unfiltered_157k

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Aerts_Fly_AdultBrain_Unfiltered_157k.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")


## See which metadata are included in the loom
names(s_cnct[['col_attrs']])

# Extract gene names (corresponds to columns of the matrix)
gns <- s_cnct[['row_attrs']][['Gene']][]
# Extract cell IDs (corresponds to rows of the matrix)
cellids <- s_cnct[['col_attrs']][['CellID']][]
# Extract count matrix (cells x genes)
mat <- s_cnct[['matrix']][,]      # dim: cells x genes
# Assign row and column names
rownames(mat) <- cellids
colnames(mat) <- gns
# If you want a data frame (optional, but not recommended for large matrices)
# raw.cnts <- as.data.frame(mat)
# For Seurat, you typically keep as a matrix (or convert to sparse Matrix)
library("Matrix")
raw.cnts <- Matrix(mat, sparse = TRUE)

# ADD METADATA
# Helper to extract and flatten metadata to 1D vector
get1d <- function(x) {
    arr <- x[]
    if (is.matrix(arr)) arr[1, ] else as.vector(arr)
}
# Get all metadata field names except CellID
meta_fields <- setdiff(names(s_cnct[['col_attrs']]), "CellID")
# Create metadata data frame, first column is cellID
metadata <- data.frame(cellID = cellids, stringsAsFactors = FALSE)
# Add all metadata fields
for (field in meta_fields) {
  vec <- get1d(s_cnct[['col_attrs']][[field]])
  if(length(vec) != length(cellids)) {
    warning(sprintf("Skipping metadata field '%s': length %d != number of cells %d", 
                    field, length(vec), length(cellids)))
    next
  }
  metadata[[field]] <- vec
}
# Set rownames to cell IDs
rownames(metadata) <- cellids
# Check result
str(metadata)

## Create the Seurat object:
Zhang_Human_Brain_30k <- CreateSeuratObject(counts = t(raw.cnts),
                            project = "Zhang_Human_Brain_30k",
                            assay = "RNA",
                            meta.data = metadata)


### Find metadata containing the cell type
unique(Zhang_Human_Brain_30k$annot1)
unique(Zhang_Human_Brain_30k$annot2)
unique(Zhang_Human_Brain_30k$lvl1class)



####################################
# USE UMAP FROM THE AUTHOR - Embeddings coordinates #########
####################################
# Extract original embedding coordinates and annotation
umap_x <- s_cnct[['col_attrs']][['Embeddings_X']][]
umap_y <- s_cnct[['col_attrs']][['Embeddings_Y']][]
annot1 <- s_cnct[['col_attrs']][['annot1']][]
lvl1class <- s_cnct[['col_attrs']][['lvl1class']][]


pdf("output/seurat/Embeddings-Zhang_Human_Brain_30k-annot1.pdf", width=10, height=6)
df <- data.frame(UMAP_1 = umap_x, UMAP_2 = umap_y, annot1 = annot1)
ggplot(df, aes(x = UMAP_1.1, y = UMAP_2.1, color = annot1)) +
  geom_point(size = 0.5, alpha = 0.8) +
  theme_bw() +
  labs(title = "Original Embedding from Loom", color = "annot1")
dev.off()

#--> Embeddings looks bad; lets prefer re-clustering it


##############################
# CLUSTERING AND UMAPS ##########
##############################

DefaultAssay(Zhang_Human_Brain_30k) <- "RNA"

## NORMALIZE AND SCALE DATA 
Zhang_Human_Brain_30k <- NormalizeData(Zhang_Human_Brain_30k, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(Zhang_Human_Brain_30k)
Zhang_Human_Brain_30k <- ScaleData(Zhang_Human_Brain_30k, features = all.genes) # zero-centres and scales it


## Find variable features
Zhang_Human_Brain_30k = FindVariableFeatures(Zhang_Human_Brain_30k, nfeatures = 3000)

Zhang_Human_Brain_30k <- RunPCA(Zhang_Human_Brain_30k, verbose = FALSE, npcs = 25)
Zhang_Human_Brain_30k <- RunUMAP(Zhang_Human_Brain_30k, reduction = "pca", dims = 1:25, verbose = FALSE)
Zhang_Human_Brain_30k <- FindNeighbors(Zhang_Human_Brain_30k, reduction = "pca", k.param = 15, dims = 1:25)
Zhang_Human_Brain_30k <- FindClusters(Zhang_Human_Brain_30k, resolution = 0.3, verbose = FALSE, algorithm = 4)



# Plot to confirm it work

pdf("output/seurat/UMAP-Zhang_Human_Brain_30k-dim25kparam15res03-annot1.pdf", width=10, height=6)
DimPlot(Zhang_Human_Brain_30k, reduction = "umap", label=TRUE, group.by = "annot1")
dev.off()

pdf("output/seurat/UMAP-Zhang_Human_Brain_30k-dim25kparam15res03-lvl1class.pdf", width=10, height=6)
DimPlot(Zhang_Human_Brain_30k, reduction = "umap", label=TRUE, group.by = "lvl1class")
dev.off()


```





## Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")
## See which metadata are included in the loom
names(s_cnct[['col_attrs']])

# Extract gene names (corresponds to columns of the matrix)
gns <- s_cnct[['row_attrs']][['Gene']][]
# Extract cell IDs (corresponds to rows of the matrix)
cellids <- s_cnct[['col_attrs']][['CellID']][]
# Extract count matrix (cells x genes)
mat <- s_cnct[['matrix']][,]      # dim: cells x genes
# Assign row and column names
rownames(mat) <- cellids
colnames(mat) <- gns
# If you want a data frame (optional, but not recommended for large matrices)
# raw.cnts <- as.data.frame(mat)
# For Seurat, you typically keep as a matrix (or convert to sparse Matrix)
library("Matrix")
raw.cnts <- Matrix(mat, sparse = TRUE)

# ADD METADATA
# Helper to extract and flatten metadata to 1D vector
get1d <- function(x) {
    arr <- x[]
    if (is.matrix(arr)) arr[1, ] else as.vector(arr)
}
# Get all metadata field names except CellID
meta_fields <- setdiff(names(s_cnct[['col_attrs']]), "CellID")
# Create metadata data frame, first column is cellID
metadata <- data.frame(cellID = cellids, stringsAsFactors = FALSE)
# Add all metadata fields
for (field in meta_fields) {
  vec <- get1d(s_cnct[['col_attrs']][[field]])
  if(length(vec) != length(cellids)) {
    warning(sprintf("Skipping metadata field '%s': length %d != number of cells %d", 
                    field, length(vec), length(cellids)))
    next
  }
  metadata[[field]] <- vec
}
# Set rownames to cell IDs
rownames(metadata) <- cellids
# Check result
str(metadata)

## Create the Seurat object:
Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain <- CreateSeuratObject(counts = t(raw.cnts),
                            project = "Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain",
                            assay = "RNA",
                            meta.data = metadata)


### Find metadata containing the cell type
unique(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain$cell_type)
unique(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain$leiden)
unique(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain$sample_id)
unique(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain$ClusterID)


##############################
# CLUSTERING AND UMAPS ##########
##############################

DefaultAssay(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain) <- "RNA"

## NORMALIZE AND SCALE DATA 
Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain <- NormalizeData(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain)
Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain <- ScaleData(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain, features = all.genes) # zero-centres and scales it


## Find variable features
Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain = FindVariableFeatures(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain, nfeatures = 3000)

Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain <- RunPCA(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain, verbose = FALSE, npcs = 25)
Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain <- RunUMAP(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain, reduction = "pca", dims = 1:25, verbose = FALSE)
Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain <- FindNeighbors(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain, reduction = "pca", k.param = 15, dims = 1:25)
Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain <- FindClusters(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain, resolution = 0.3, verbose = FALSE, algorithm = 3) # switch to algorithm 3 from algo 4 to avoid ERROR `long vectors not supported`



# Plot to confirm it work

pdf("output/seurat/UMAP-Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain-dim25kparam15res03algo3-cell_type.pdf", width=30, height=10)
DimPlot(Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain, reduction = "umap", label=TRUE, group.by = "cell_type", raster=FALSE)
dev.off()

```

--> All good!








## Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")
## See which metadata are included in the loom
names(s_cnct[['col_attrs']])

# Extract gene names (corresponds to columns of the matrix)
gns <- s_cnct[['row_attrs']][['Gene']][]
# Extract cell IDs (corresponds to rows of the matrix)
cellids <- s_cnct[['col_attrs']][['CellID']][]
# Extract count matrix (cells x genes)
mat <- s_cnct[['matrix']][,]      # dim: cells x genes
# Assign row and column names
rownames(mat) <- cellids
colnames(mat) <- gns
# If you want a data frame (optional, but not recommended for large matrices)
# raw.cnts <- as.data.frame(mat)
# For Seurat, you typically keep as a matrix (or convert to sparse Matrix)
library("Matrix")
raw.cnts <- Matrix(mat, sparse = TRUE)

# ADD METADATA
# Helper to extract and flatten metadata to 1D vector
get1d <- function(x) {
    arr <- x[]
    if (is.matrix(arr)) arr[1, ] else as.vector(arr)
}
# Get all metadata field names except CellID
meta_fields <- setdiff(names(s_cnct[['col_attrs']]), "CellID")
# Create metadata data frame, first column is cellID
metadata <- data.frame(cellID = cellids, stringsAsFactors = FALSE)
# Add all metadata fields
for (field in meta_fields) {
  vec <- get1d(s_cnct[['col_attrs']][[field]])
  if(length(vec) != length(cellids)) {
    warning(sprintf("Skipping metadata field '%s': length %d != number of cells %d", 
                    field, length(vec), length(cellids)))
    next
  }
  metadata[[field]] <- vec
}
# Set rownames to cell IDs
rownames(metadata) <- cellids
# Check result
str(metadata)

## Create the Seurat object:
Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain <- CreateSeuratObject(counts = t(raw.cnts),
                            project = "Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain",
                            assay = "RNA",
                            meta.data = metadata)


### Find metadata containing the cell type
unique(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain$ClusterID)
unique(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain$annotation)


##############################
# CLUSTERING AND UMAPS ##########
##############################

DefaultAssay(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain) <- "RNA"

## NORMALIZE AND SCALE DATA 
Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain <- NormalizeData(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain)
Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain <- ScaleData(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain, features = all.genes) # zero-centres and scales it


## Find variable features
Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain = FindVariableFeatures(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain, nfeatures = 3000)

Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain <- RunPCA(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain, verbose = FALSE, npcs = 25)
Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain <- RunUMAP(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain, reduction = "pca", dims = 1:25, verbose = FALSE)
Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain <- FindNeighbors(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain, reduction = "pca", k.param = 15, dims = 1:25)
Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain <- FindClusters(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain, resolution = 0.3, verbose = FALSE, algorithm = 3) # switch to algorithm 3 from algo 4 to avoid ERROR `long vectors not supported`



# Plot to confirm it work

pdf("output/seurat/UMAP-Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain-dim25kparam15res03algo3-annotation.pdf", width=10, height=6)
DimPlot(Ravenscroft_et_al_2019_LarvalBrain_FROMflybrain, reduction = "umap", label=TRUE, group.by = "annotation", raster=FALSE)
dev.off()

```

--> All good!






## Davie_Janssens_Koldere_et_al_2018_Adult_Brain_HARMONY_SCENIC

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/5a3ba000_20180809__Davie_Janssens_Koldere_et_al_2018_Adult_Brain.HARMONY_SCENIC.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")
## See which metadata are included in the loom
names(s_cnct[['col_attrs']])


```


--> No cell types




## Davie_Janssens_Koldere_et_al_2018_Adult_Brain_BBKNN

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/5a3ba000_20200313__Davie_Janssens_Koldere_et_al_2018_Adult_Brain.BBKNN.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")
## See which metadata are included in the loom
names(s_cnct[['col_attrs']])


```


--> No cell types


## Aerts_Supplemental_Fly_AdultBrain_DropSeq_2k

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Aerts_Supplemental_Fly_AdultBrain_DropSeq_2k.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")


## See which metadata are included in the loom
names(s_cnct[['col_attrs']])


# Extract gene names (corresponds to columns of the matrix)
gns <- s_cnct[['row_attrs']][['Gene']][]
# Extract cell IDs (corresponds to rows of the matrix)
cellids <- s_cnct[['col_attrs']][['CellID']][]
# Extract count matrix (cells x genes)
mat <- s_cnct[['matrix']][,]      # dim: cells x genes
# Assign row and column names
rownames(mat) <- cellids
colnames(mat) <- gns
# If you want a data frame (optional, but not recommended for large matrices)
# raw.cnts <- as.data.frame(mat)
# For Seurat, you typically keep as a matrix (or convert to sparse Matrix)
library("Matrix")
raw.cnts <- Matrix(mat, sparse = TRUE)

# ADD METADATA
# Helper to extract and flatten metadata to 1D vector
get1d <- function(x) {
    arr <- x[]
    if (is.matrix(arr)) arr[1, ] else as.vector(arr)
}
# Get all metadata field names except CellID
meta_fields <- setdiff(names(s_cnct[['col_attrs']]), "CellID")
# Create metadata data frame, first column is cellID
metadata <- data.frame(cellID = cellids, stringsAsFactors = FALSE)
# Add all metadata fields
for (field in meta_fields) {
  vec <- get1d(s_cnct[['col_attrs']][[field]])
  if(length(vec) != length(cellids)) {
    warning(sprintf("Skipping metadata field '%s': length %d != number of cells %d", 
                    field, length(vec), length(cellids)))
    next
  }
  metadata[[field]] <- vec
}
#--> Bug nb of dim

```



--> Cannot import file







## Zhang_Human_Brain_60k_old

To *import and convert to Seurat loom files* I follow discussion from [this](https://github.com/satijalab/seurat/issues/5124).



```bash
conda activate scRNAseqV2

# install SeuratDisk to import loom files Follow: https://github.com/mojaveazure/seurat-disk
remotes::install_github("mojaveazure/seurat-disk") 
#--> ALL GOOD!
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
library("SeuratDisk")


##############################
# DATA IMPORT ##########
##############################


# import loom files
## Extract the data from .loom corresponding to the cell IDs, gene names and raw count matrix:
loom.fn <- 'input/Zhang_Human_Brain_60k_old.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")

## See which metadata are included in the loom
names(s_cnct[['col_attrs']])

# Extract gene names (corresponds to columns of the matrix)
gns <- s_cnct[['row_attrs']][['Gene']][]
# Extract cell IDs (corresponds to rows of the matrix)
cellids <- s_cnct[['col_attrs']][['CellID']][]
# Extract count matrix (cells x genes)
mat <- s_cnct[['matrix']][,]      # dim: cells x genes
# Assign row and column names
rownames(mat) <- cellids
colnames(mat) <- gns
# If you want a data frame (optional, but not recommended for large matrices)
# raw.cnts <- as.data.frame(mat)
# For Seurat, you typically keep as a matrix (or convert to sparse Matrix)
library("Matrix")
raw.cnts <- Matrix(mat, sparse = TRUE)

# ADD METADATA
# Helper to extract and flatten metadata to 1D vector
get1d <- function(x) {
    arr <- x[]
    if (is.matrix(arr)) arr[1, ] else as.vector(arr)
}
# Get all metadata field names except CellID
meta_fields <- setdiff(names(s_cnct[['col_attrs']]), "CellID")
# Create metadata data frame, first column is cellID
metadata <- data.frame(cellID = cellids, stringsAsFactors = FALSE)
# Add all metadata fields
for (field in meta_fields) {
  vec <- get1d(s_cnct[['col_attrs']][[field]])
  if(length(vec) != length(cellids)) {
    warning(sprintf("Skipping metadata field '%s': length %d != number of cells %d", 
                    field, length(vec), length(cellids)))
    next
  }
  metadata[[field]] <- vec
}
# Set rownames to cell IDs
rownames(metadata) <- cellids
# Check result
str(metadata)

## Create the Seurat object:
Zhang_Human_Brain_60k_old <- CreateSeuratObject(counts = t(raw.cnts),
                            project = "Zhang_Human_Brain_60k_old",
                            assay = "RNA",
                            meta.data = metadata)


### Find metadata containing the cell type
unique(Zhang_Human_Brain_60k_old$Brain.region) # good
unique(Zhang_Human_Brain_60k_old$Cell.subtype..provided.by.authors.) # gene name
unique(Zhang_Human_Brain_60k_old$Cell.type..provided.by.authors.) # GOOD!
unique(Zhang_Human_Brain_60k_old$Cell.type..major.groups.) # GOOD!



##############################
# CLUSTERING AND UMAPS ##########
##############################

DefaultAssay(Zhang_Human_Brain_60k_old) <- "RNA"

## NORMALIZE AND SCALE DATA 
Zhang_Human_Brain_60k_old <- NormalizeData(Zhang_Human_Brain_60k_old, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(Zhang_Human_Brain_60k_old)
Zhang_Human_Brain_60k_old <- ScaleData(Zhang_Human_Brain_60k_old, features = all.genes) # zero-centres and scales it


## Find variable features
Zhang_Human_Brain_60k_old = FindVariableFeatures(Zhang_Human_Brain_60k_old, nfeatures = 3000)

Zhang_Human_Brain_60k_old <- RunPCA(Zhang_Human_Brain_60k_old, verbose = FALSE, npcs = 25)
Zhang_Human_Brain_60k_old <- RunUMAP(Zhang_Human_Brain_60k_old, reduction = "pca", dims = 1:25, verbose = FALSE)
Zhang_Human_Brain_60k_old <- FindNeighbors(Zhang_Human_Brain_60k_old, reduction = "pca", k.param = 15, dims = 1:25)
Zhang_Human_Brain_60k_old <- FindClusters(Zhang_Human_Brain_60k_old, resolution = 0.3, verbose = FALSE, algorithm = 4)



# Plot to confirm it work

pdf("output/seurat/UMAP-Zhang_Human_Brain_60k_old-dim25kparam15res03-Cell.type..provided.by.authors.pdf", width=10, height=6)
DimPlot(Zhang_Human_Brain_60k_old, reduction = "umap", label=TRUE, group.by = "Cell.type..provided.by.authors.")
dev.off()

pdf("output/seurat/UMAP-Zhang_Human_Brain_60k_old-dim25kparam15res03-Cell.type..major.groups.pdf", width=10, height=6)
DimPlot(Zhang_Human_Brain_60k_old, reduction = "umap", label=TRUE, group.by = "Cell.type..major.groups.")
dev.off()
pdf("output/seurat/UMAP-Zhang_Human_Brain_60k_old-dim25kparam15res03-Brain.region.pdf", width=10, height=6)
DimPlot(Zhang_Human_Brain_60k_old, reduction = "umap", label=TRUE, group.by = "Brain.region")
dev.off()



### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(Zhang_Human_Brain_60k_old, file = "output/seurat/Zhang_Human_Brain_60k_old-dim25kparam15res03.rds")
Zhang_Human_Brain_60k_old <- readRDS("output/seurat/Zhang_Human_Brain_60k_old-dim25kparam15res03.rds")
##########################################################

# Volcano plot of expression

gene_list = c("EZH1", "EZH2", "EED", "SUZ12", "PHF1", "MTF2", "PHF19", "JARID2", "AEBP2", "RBBP4", "RBBP7","LCOR","LCORL") #
#--> Missing genes EPOP

pdf("output/seurat/VlnPlot-Zhang_Human_Brain_60k_old-Cell.type..provided.by.authors-genelist.pdf", width = 10, height = 3)
for (gene in gene_list) {
  print(paste("Generating plot for:", gene))
  p <- VlnPlot(Zhang_Human_Brain_60k_old, 
               features = gene, 
               group.by = "Cell.type..provided.by.authors.",  # x-axis is your annotation
               pt.size = 0) +
    theme(plot.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(p)
}
dev.off()



```




## GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX


Dataset download from [GSE107451](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107451):
  - GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv.tar.gz 
  - GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv.gz 
  - GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata_README.txt 




```bash
conda activate scRNAseqV2


# extract data
tar -xvzf input/GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv.tar.gz

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
library("SeuratDisk")
library(data.table)


# Load files
counts <- readMM("matrix.mtx") #
genes <- read.delim("genes.tsv", header = FALSE)
barcodes <- read.delim("barcodes.tsv", header = FALSE)

rownames(counts) <- genes$V2        # Use gene symbols as rownames
colnames(counts) <- barcodes$V1     # Use cell barcodes as colnames


meta <- fread("input/GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv.gz")

# 3. Find intersecting barcodes
shared_barcodes <- intersect(colnames(counts), meta$new_barcode)
counts_sub <- counts[, shared_barcodes]
meta_sub <- meta[match(shared_barcodes, meta$new_barcode), ]
rownames(meta_sub) <- shared_barcodes
# Double-check
stopifnot(identical(rownames(meta_sub), colnames(counts_sub)))

# Convert to a pure data.frame BEFORE setting rownames
meta_sub_df <- as.data.frame(meta_sub)
rownames(meta_sub_df) <- shared_barcodes

# Check order, must be identical!
stopifnot(identical(rownames(meta_sub_df), colnames(counts_sub)))

# Create Seurat object
WholeBrain_57k <- CreateSeuratObject(counts = counts_sub, meta.data = meta_sub_df)

# Check annotation
head(WholeBrain_57k$annotation)
table(is.na(WholeBrain_57k$annotation))




### Find metadata containing the cell type
unique(WholeBrain_57k$annotation) # good
unique(WholeBrain_57k$Age) # 
unique(WholeBrain_57k$Genotype) # 

# SAVE all cluster name to .txt #########
annots <- unique(WholeBrain_57k$annotation)
# Keep only those that are NOT just numbers
named_annots <- annots[!grepl("^\\d+$", annots)]
print(named_annots)
write.table(named_annots, file = "output/seurat/WholeBrain_57k-named_clusters.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
##########################################



##############################
# CLUSTERING AND UMAPS ##########
##############################

DefaultAssay(WholeBrain_57k) <- "RNA"

## NORMALIZE AND SCALE DATA 
WholeBrain_57k <- NormalizeData(WholeBrain_57k, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(WholeBrain_57k)
WholeBrain_57k <- ScaleData(WholeBrain_57k, features = all.genes) # zero-centres and scales it


## Find variable features
WholeBrain_57k = FindVariableFeatures(WholeBrain_57k, nfeatures = 3000)

WholeBrain_57k <- RunPCA(WholeBrain_57k, verbose = FALSE, npcs = 25)
WholeBrain_57k <- RunUMAP(WholeBrain_57k, reduction = "pca", dims = 1:25, verbose = FALSE)
WholeBrain_57k <- FindNeighbors(WholeBrain_57k, reduction = "pca", k.param = 15, dims = 1:25)
WholeBrain_57k <- FindClusters(WholeBrain_57k, resolution = 0.3, verbose = FALSE, algorithm = 3) # algo 3 as ERROR `long vector not supported`



# Plot to confirm it work
pdf("output/seurat/UMAP-WholeBrain_57k-dim25kparam15res03-annotation.pdf", width=30, height=15)
DimPlot(WholeBrain_57k, reduction = "umap", label=TRUE, group.by = "annotation", label.size= 8)
dev.off()

pdf("output/seurat/UMAP-WholeBrain_57k-dim25kparam15res03-Age.pdf", width=10, height=6)
DimPlot(WholeBrain_57k, reduction = "umap", label=TRUE, group.by = "Age", label.size= 8)
dev.off()

pdf("output/seurat/UMAP-WholeBrain_57k-dim25kparam15res03-splitAge.pdf", width=50, height=6)
DimPlot(WholeBrain_57k, reduction = "umap", label=TRUE, split.by = "Age", label.size= 8)
dev.off()



### SHOW EXPRESSION OF SOME GENES #######################
#saveRDS(WholeBrain_57k, file = "output/seurat/WholeBrain_57k-dim25kparam15res03.rds")
WholeBrain_57k <- readRDS("output/seurat/WholeBrain_57k-dim25kparam15res03.rds")
##########################################################


# Filter cell types and time point - keep all ages except 0days
cell_types_of_interest <- c(
  "Astrocyte-like", "Dopaminergic", "Serotonergic", "Peptidergic", 
  "Olfactory_projection_neurons", "MBON", "Plasmatocytes", 
  "Ensheathing_glia", "Perineurial_glia", "Subperineurial_glia", 
  "Cortex_glia", "Chiasm_glia", "DCN"
)

WholeBrain_57k_filteredv1 <- subset(
  WholeBrain_57k,
  subset = annotation %in% cell_types_of_interest & Age != 0
)

table(WholeBrain_57k_filteredv1$annotation)
table(WholeBrain_57k_filteredv1$Age)



DefaultAssay(WholeBrain_57k_filteredv1) <- "RNA"

## NORMALIZE AND SCALE DATA 
WholeBrain_57k_filteredv1 <- NormalizeData(WholeBrain_57k_filteredv1, normalization.method = "LogNormalize", scale.factor = 10000) # accounts for the depth of sequencing
all.genes <- rownames(WholeBrain_57k_filteredv1)
WholeBrain_57k_filteredv1 <- ScaleData(WholeBrain_57k_filteredv1, features = all.genes) # zero-centres and scales it

## Find variable features
WholeBrain_57k_filteredv1 = FindVariableFeatures(WholeBrain_57k_filteredv1, nfeatures = 3000)

WholeBrain_57k_filteredv1 <- RunPCA(WholeBrain_57k_filteredv1, verbose = FALSE, npcs = 25)
WholeBrain_57k_filteredv1 <- RunUMAP(WholeBrain_57k_filteredv1, reduction = "pca", dims = 1:25, verbose = FALSE)
WholeBrain_57k_filteredv1 <- FindNeighbors(WholeBrain_57k_filteredv1, reduction = "pca", k.param = 15, dims = 1:25)
WholeBrain_57k_filteredv1 <- FindClusters(WholeBrain_57k_filteredv1, resolution = 0.3, verbose = FALSE, algorithm = 3) # algo 3 as ERROR `long vector not supported`


# Plot to confirm it work
pdf("output/seurat/UMAP-WholeBrain_57k_filteredv1-dim25kparam15res03-annotation.pdf", width=10, height=6)
DimPlot(WholeBrain_57k_filteredv1, reduction = "umap", label=TRUE, group.by = "annotation", label.size= 4)
dev.off()


```



