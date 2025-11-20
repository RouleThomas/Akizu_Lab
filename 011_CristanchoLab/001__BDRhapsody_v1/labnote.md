# Project

Test the BD Rhapsody pipeline for scRNAseq. Compare result with already processed data.

Design:
- Normoxia vs hypoxia at fetal (neonatal) stage
    - 2 Bio rep per condition (1 male, 1 female)


--> Previous marker genes used by Ana can be found at `011*/docs/Annotation.xlsx`

# Data access

Access Cristancho lab folder from CHOP computer: `/mnt/isilon/cristancho_data`

Access data from [Seven Bridge - BD Rhapsody](https://igor.sbgenomics.com/home)
--> Create [Seven Bridge](https://bd-rhapsody-bioinfo-docs.genomics.bd.com/setup/sbg/top_sbg_setup.html) account and ask data file access to Paulo 


```bash
wget -O output/seurat/AC_WTA_SMK_index_library_Seurat.rds \
'https://sbg-main.s3.amazonaws.com/adf82bfb-5651-4497-9669-8151287ea02d%2BAC_WTA_SMK_index_library_Seurat.rds?x-username=roulet&x-requestId=7e137078-74a4-45cf-a7a6-b89ea3420279&x-project=deoliveirp%2Fnextseq-hypoxia&response-content-disposition=attachment%3Bfilename%3DAC_WTA_SMK_index_library_Seurat.rds&response-content-type=application%2Foctet-stream&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20251120T164843Z&X-Amz-SignedHeaders=host&X-Amz-Expires=172800&X-Amz-Credential=AKIAJH6BPOGIWTDUABEQ%2F20251120%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=1858853638f9c40f2e4742e98be3ee51b135e5a5c0c4315a90838d68f4865d5a'
```


Available rds files (in *Seven Bridge*:  `Analysis/Tasks`):
- **Jul 16, 2025 15:40**: ATACMultiomewithST-CC-additionalSMK_Seurat.rds
- **May 12, 2025 13:54**: ATACMultiomewithST-CC_Seurat.rds
- **Apr 15, 2025 19:16**: AC_WTA_SMK_index_library_Seurat.rds


# Data analysis -V1



```bash
conda activate SignacV5
module load hdf5
```

```R
set.seed(42)

# library
library("Signac")
library("Seurat")


# load seurat object
multiome_WT_Bap1KO_QCV2vC1.sct <- readRDS(file = "output/seurat/ATACMultiomewithST-CC-additionalSMK_Seurat.rds")





```






