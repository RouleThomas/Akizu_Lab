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
    - Zhang_Human_Brain_30k --> LOOKS GOOD! Clear cell type, cluster well **FROM Davie, Jannssens and Koldere et al., 2018 (Cell)**
    - Waddell_CentralBrain_10k --> No cell type information
    - Ravenscroft_et_al_2019_LarvalBrain --> Weird naming of cell types **FROM  Janssens et al., 2022 (Nature)**
    - Aerts_Fly_AdultBrain_Unfiltered_157k -->  No cell type information
    - Davie_Janssens_Koldere_et_al_2018_AdultBrain_FROMflybrain --> Downloaded from [here](https://flybrain.aertslab.org/?_inputs_&page=%22Downloads%22)
    - 5a3ba000_20180809__Davie_Janssens_Koldere_et_al_2018_Adult_Brain_HARMONY_SCENIC **FROM FlyCellAtlas [here](https://scope.aertslab.org/#/FlyCellAtlas/*/welcome)** -->  No cell type information
    - 5a3ba000_20200313__Davie_Janssens_Koldere_et_al_2018_Adult_Brain_BBKNN **FROM FlyCellAtlas [here](https://scope.aertslab.org/#/FlyCellAtlas/*/welcome)** -->  No cell type information
    - Aerts_Supplemental_Fly_AdultBrain_DropSeq_2k





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
#--> BUG canot import error 1 vs 2 dim



```


--> Error cannot import







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

--> Cell types not clear








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



