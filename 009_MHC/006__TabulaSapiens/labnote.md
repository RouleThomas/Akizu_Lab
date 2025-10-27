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

XXXY HERE

```bash
conda activate bpcells_r44 # For seurat v5 with BPCells
```



```R
library("Seurat")
library("BPCells")
library("Azimuth")
library("biomaRt")


# load rds object
(seu, file = "output/seurat/TabulaSapiens_allCells-dim30.rds")





```




