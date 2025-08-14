# Daily-To-Do-List
--> On Benchling

Server name: 

# Infos
## Usefull commands

- If a window is buggy, stuck; run: `killall -3 gnome-shell`

## File transfer
Only on the CHOP cluster or on a CHOP machine...
Just copy/paste with cp from cluster to local computer Google Drive!

## RStudio with ressource
- Go to https://respublica-web.research.chop.edu and interactive app --> Rstudio.
- Here I have access to my `/home/roulet/`

## IGV with ressource
- Maybe try Mobaxterm with ssh and -X

To run IGV: type `igv.sh` or `igv_hidpi.sh` depending on computer reslution


## Cluster-files folder-like interace
No idea! Dolphin I can see files but not transfer them...\
To allow cp/paste between Vitrual machine and computer: 


## Directories
- BIOINFORMATICS SPACE: Working directory in the cluster is `/scr1/users/roulet` 30TB limit (temporary space, belong to CHOP)-
- SPACE TO SEE FILE (folder-like env): `/home/roulet/`
- LOCAL COMPUTER FILE:\
**Box**:`/home/roulet/tsclient/roule/Box`\
**Google drive**: `/home/roulet/tsclient/roule/Google Drive Streaming`



## Run job

- Interactive `srun --mem=20g --pty bash -l`. Or `srun --partition=gpuq --mem=20g --gres=gpu:1 --pty bash -l`
- Interactive with multiple cores: `srun --mem=500g --cpus-per-task=10 --pty bash -l` (`nproc` to check nb of cores;)
- Can do parrallel processing; asking mutliple core with unique memories: `srun --cpus-per-task=8 --mem-per-cpu=40g --pty bash -l` this is equivalent of 8*40 = 320g memory; but spread in 8 cpus (use `library(furrr)` in R)
- Sbatch `sbatch job.sh` (edit script in VSC, then create it on the cluster with `touch script.sh` edit it with `nano script.sh` copy paste from VSC, and run it)
- List jobs: `squeue -u roulet` (`scancel [JOBID]` to cancel)


If encounteer `bash __vte_prompt_command command not found` error message. Do the following:
1. add this add the end of the ~/.bashrc file:
```bash
__vte_prompt_command() {
  local pwdmaxlen=30
  local pwdoffset=$(( ${#PWD} - pwdmaxlen ))
  local pwdir=${PWD:-$HOME}

  [ "${pwdoffset}" -gt "0" ] && pwdir="…${pwdir:${pwdoffset}}"
  printf "\033]0;%s@%s:%s\007" "${USER}" "${HOSTNAME%%.*}" "${pwdir}"
}
```
2. then `source ~/.bashrc`


# Transfer to new cluster RES-RHEL-RH9HPC and conda env:
## conda environment:

### deseq2
R/4.3.0 with **deseq2** and **ChIPseeker** notably + following libraries:
library("DESeq2")
library("tidyverse") 
library("RColorBrewer")
library("pheatmap") 
library("apeglm") 
library("factoextra") 
library("gridExtra")
library("rtracklayer") 
library("readxl") 
library("ggpubr")
library("dendextend") 
library("clusterProfiler") 
library("pathview") 
library("DOSE")
library("org.Hs.eg.db") 
library("enrichplot")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("meshes")
library("ReactomePA")
library("VennDiagram")
library("ggrepel")
library("fgsea")

library("UpSetR") # venn diagram like plot
library("furrr") # parralel processing



### deseq2V2

Same as deseq2 with `apeglm` R package updated (see error from `005_SNX14/001/RNAseq`)

### deseq2V3

Same as deseq2V2 but with:
- `ChIPQC` [R package installed](https://bioconductor.org/packages/release/bioc/html/ChIPQC.html)
(see creation in  `008*/001/ # Generate bigwig coverage files`)




### ChIPseqSpikeInFree
R/3.6.1 with **ChIPseqSpikeInFree** notably + following libraries:
library("Rsamtools")
library("GenomicAlignments")
library("ChIPseqSpikeInFree")
library("devtools")


PERL, CPAN, MCPAN installation

diffreps

### DiffBind
R/4.1.0 **DiffBind** notably + following libraries:

library("DiffBind") 
library("csaw") 


### featurecounts

featurecounts

### BedToBigwig

bedtools, genomeCoverageBed, bedGraphToBigWig, wiggletools 

### wigtobigwig

wigtobigwig



### scRNAseq
R/4.3.0 with deseq2 and ChIPseeker and Seurat notably + following libraries:
library("DESeq2")
library("tidyverse") 
library("RColorBrewer")
library("pheatmap") 
library("apeglm") 
library("factoextra") 
library("gridExtra")
library("rtracklayer") 
library("readxl") 
library("ggpubr")
library("dendextend") 
library("clusterProfiler") 
library("pathview") 
library("DOSE")
library("org.Hs.eg.db") 
library("enrichplot")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("meshes")
library("ReactomePA")
library("VennDiagram")
library("devtools")
library("Seurat")
library("AnnotationHub")
library("SingleCellExperiment")
library("ensembldb")



### binBw_v2
copy of scRNAseq; with R [PopSV](https://rdrr.io/github/jmonlong/PopSV/) installed
--> To use bin.bw R function (count bigwig reads in specific bed regions or bins)

WORK great!

### scRNAseqV2

--> Same scRNAseq but also allow Shiny app creation!
--> And with harmony
--> And SeuratDisk
--> And anndata (read .h5ad file)
--> And zellkonverter  (read/convert .h5ad file)

*NOTE: BUG after adding ability to read .h5ad file, no longer Leiden algorithm for command like this for example `WT_Kcnc1_p35_CX_1step.sct <- FindClusters(WT_Kcnc1_p35_CX_1step.sct, resolution = 0.7, verbose = FALSE, algorithm = 4, method = "igraph")` --> # algorithm = 4 is Leiden...*


### scRNAseqV3

--> Same `scRNAseqV2` but with MuSic2 and bisque and without  SeuratDisk, anndata, zellkonverter : Can use `algorithm = 4` please see note in `### scRNAseqV2`



### SignacV5

Env for scRNAseq and scATACseq integration. Installed in `002*/006*` *labnote.md*
Also R: BiocManager, tidyverse, EnsDb.Hsapiens.v86, devtools, ShinyCell, rsconnect (shiny App), ggbeeswarm (arrange point of violin plot)



### Signac_Pando
R4.3.3

--> clone from `SignacV5`; added [Pando](https://quadbio.github.io/Pando/index.html) and [GenomicScores](https://www.bioconductor.org/packages/devel/bioc/vignettes/GenomicScores/inst/doc/GenomicScores.html)


### valiDrops_v1

clone of `scRNAseq` with R `valiDrops` and `presto` and `segments` v1.6.4



### NeuronChat

clone of `scRNAseq` with R `NeuronChat`, `CellChat` v1 and `ComplexHeatmap`


### NeuronChat

clone of `scRNAseq` with R  `CellChat`v2 and `ComplexHeatmap` and `presto`



### monocle3
Only `monocle3`, `Seurat_v_5`, `BiocManager` and `glmGamPoi` and `EnhancedVolcano`
(Can run SCTransform in seurat v5)

Also added in R: `deseq2`, `edgeR`, `scater`

--> Can be used for DESEQ2 deg of scRNAseq

### monocle3_V1

--> Same but also with `SeuratWrappers`

### condiments_V1

-->  `condiments` only

### condiments_V5

--> `condiments` with `Seurat_v_4`

*--> Can delete condiments_V2, V3, V4*

### condiments_V6

--> `pheatmap` R package added


### granulator

R4.3.1 with `granulator` and `tidyverse`



### homer

-->  `homer` only; can be used for peak calling or motif... (used in `008*/001*`)

### homer_deseq2_V1

--> `homer` with `deseq2` in R v4.3.3  (used in `008*/001*`)


### condiments_Signac (Seurat v5 here)

--> condiments with Signac!!! (and Seurat)

--> Use to do pseudotime analysis with Signac-dependent seurat object

--> tradeseq installed (can do DEG pseudotime)

--> zellkonverter installed (read h5da)

### salmon1

Salmon quantification tool.



# cool random command


*Copy all files within a folder, excluding a specific pattern/string*: `for file in *; do if [[ ! $file =~ \.fq\.gz$ ]]; then cp "$file" /path/to/destination/; fi; done`

Count nb of unique genes in a gtf: `awk '$3 == "gene" {print $10}' your_file.gtf | sort | uniq | wc -l`




