# Project


WT vs *Bap1KO* mouse scRNA/ATACseq

BAP1 involved in NDD


# Docs


Integration of **scRNAseq/ATACseq** tutorial [here](https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/scAtacSeq/lab-sc_atac_seq.html) and [here](https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette).


--> Start with some QC about ATACseq and RNAseq and get back to you
--> Provide to me gene marker for expected cell types (which cell type we expect?). Can use automatic annotation but guided with cell marker known always better

What I could do 1st:
--> Use scRNASeq data to annotate cell types and then integrate ATACseq information
--> Identify motif from the ATACseq peak
--> DEGs WT vs KO
--> Diff. Acc. Regions WT vs KO; annotate diff. peak to genes (nearest TSS)
--> Then meeting from there what we do next

Next:
- Which promoters and enhancers become active in different cell types and conditions?
- Check if any TF binding sites are active in different cell types and conditions?
- Are some genes primed for expression, (e.g. the promoters show an open chromatin state, but the gene is not expressed yet)?



# file

- Transfer data from Hard drive to Google Drive then to HPC Cluster (transfer the aggregated cellranger counts output (`input_raw/A1B1_A2B2`) and the raw fastq (`input_raw/F001 and F002`))
- A = RNA, B = ATAC; 1,2 = genotypes



