# Project

Jasmine making a review about PRC2 subunits and brain stuff. 
- At revision stage they would like to add scRNAseq information data
    - Compare level of expression of PRC2 subunits between drosphila and human
    - Strenghten that many more types of neurons in humans than in drosophila


# Pipeline

Download and re-analize public scRNAseq data from fly brain, and human brain neurons:
- Public data in [SCOPE](https://scope.aertslab.org/#/Fly_Brain/Fly_Brain%2FRavenscroft_et_al_2019_LarvalBrain.loom/welcome):
    - Davie_Janssens_Koldere_et_al_2018_AdultBrain --> Seems not great, too complex, cell type not clear
    - Ravenscroft_et_al_2019_LarvalBrain --> XXX
    - Aerts_Fly_AdultBrain_Filtered_57k --> XXX



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



XXY HERE CHECK OTHER LOOM DATA


















