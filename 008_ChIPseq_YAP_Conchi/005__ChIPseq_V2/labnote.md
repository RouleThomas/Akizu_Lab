# Project

**SE**

Additional replicate for `001__ChIPseq_V1`:

- input WT and YAPKO
- IP for QSER1, TEAD4, YAP1


- *NOTE: Files weirdly small in size... ~50-500Mb. Weird Expected >Go*


**Objectives:**
- Add aditional replicate to `001__ChIPseq_V1` and integrate all into `006__ChIPseq_V1V2`


# Pipeline
- Download data (wget)
- Rename files
- FastQC (fastqc)
- Trimming (fastp)
- Histone content (R)
- Mapping (bowtie2)
- Spike-in scaling (DiffBind)
- Bigwig generation (deepTools)
- peak calling (MACS2)
- peak assignment to gene (ChIPseeker)

--> Detail of the overall pipeline in `Meeting_20230919_draft.xlsx` 



# Download / import data

- *NOTE: data downloaded in local from basespace/dropbox and then imported into the cluster*
--> Data transferred from EZH1-related ChIPseq (`001*/012*`)




