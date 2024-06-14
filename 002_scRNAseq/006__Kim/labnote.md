# Background

BAP1 interact with some member of PRC2/1. Notably involved in histone de-ubiquitination (chromatin opening; as H2A Ub = compact).
BAP1 part of a big complex. KO is embryo lethal. 11 patients with HT BAP1 mutation = NDD / ASD syndrome --> 1.5 Fold more ubuiquitnation H2AK119Ub (as BAP1 not functional). 
Generate mouse model using CRE system to deplete BAP1. Strong phenotype in mice, hyperactive and burying assay (they do not burry anything, super strong behaviorial phenotype!) + dentate gyruse (part of the hippocampal formation in the temporal lobe of the brain) decrease along dev (CC3 expression (cell death marker) increase so could be a partial contributor of decreased size but not only that..)
During dentate gyrus development, Pax6 > Tbr2 > Prox1 gene expression switch (express + cell migrate) = marker of development (could be used as proxy for pseudotime analysis!). Found that Pax6 and Tbr2 decrease in the Mt (at e18). 

# Project

Casey Lim, MD/PhD student did the experiment.
Lab of Seon-Hee Kim 

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



