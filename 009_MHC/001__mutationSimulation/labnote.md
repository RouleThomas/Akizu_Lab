
# Projects

Collaboration with Joan from Fox Chase Cancer.

Simulate known mutation signature (SBS, DBS..) to the human reference genome/exome, generate mutated peptides and assess their binding affinity to all known MHC-I variants (ie. HLA variants)



# Backgrounds



MHC-I, can be produced by three genes HLA-A/B/C each individual has 6 alleles  Produce 6 different MHC-I complex per ind.
Each allele are prone to polymorphism! Notably HLA-A  Many different possibilities of sequence combination of HLAs. (>3k alleles variants)
Cancer mutations  abnormal proteins  broken down into peptides. MHC-I presents these peptides (neoantigens) on the cell surface for recognition and killing by CD8+ T cell 
MHC-I peptide binding recognition is variable; some people can present certain peptides/neoantigen that others cannot
Different mutagens (e.g. tobacco, UV, chemo), and cancer type, create distinct mutation signatures


--> Interesting docs:
- [Mutation signature](https://medium.com/@hylke.donker/mutational-signatures-explained-1dc435b2d7b7)

For SBS-96:
    6 mutation types
    4 possible 5’ flanking bases (A/C/G/T)
    4 possible 3’ flanking bases (A/C/G/T)
--> In total, there are 96 = 4 x 6 x 4 singlet classes when sorted by three letter — trinucleotide — motifs.



# Pipeline


- Collect ALL mutation signatures
    COSMIC Mutation Signature Database (include Alexandrov et al 2018)
    Kucab et al 2019; 41 environmental agents (6 DBS, 8 ID)
    Pich et al 2019 (therapies induced mutation)
- Simulate a 1000 random mutation signature profile for each known signature on human reference exome
- Confirm profile by generating mutation profile plot
- Collect mutated peptides for each signature 
- Assess binding affinity of MHC-I to these mutated peptide; for each HLA/MHC-I known variants (NetMHCPan3.0 ; PHBR score (Marty et al 2017))



# Tool


To simulate the mutation signature to genome/exome:
- [MSA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04450-8) --> More recent, test this one first!
- [SigProfilerSimulator](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03772-3); more info [here](https://osf.io/usxjz/wiki/2.%20Simulations/)

--> These tools are not good as they simulate an already exisiting cancer, or mutagens; so it is a group of mutation signature, not a mutation signature directly.

We instead gonna create our own scripts to directly simulate mutation signature. Joan provided first version.


## SigProfilerSimulator

### Install SigProfiler 

Tool library [here](https://cancer.sanger.ac.uk/signatures/tools/); adn let's follow the R version of [SigProfilerSimulatorR](https://github.com/AlexandrovLab/SigProfilerSimulatorR)


```bash
# Create a conda env with R and Python3
conda create --name SigProfiler --clone scRNAseq

pip install SigProfilerSimulator
```

```R
library("reticulate")
library("devtools")


install_github("AlexandrovLab/SigProfilerSimulatorR")

library("SigProfilerSimulatorR")
```



### Run SigProfilerSimulator


XXX


# Set up Joan scripts


Let's set up a conda environment to run these scripts. 

- Requirements: `bedtools2/2.31.0`, `samtools/1.17` (Can be loaded with `module load`)

--> For this project, keep the same organization as Joan; so wordir = `009*/001*/mutsim`


```bash
# Create conda env
conda create -n mutsim python=3.11 -y
conda activate mutsim
conda install -y -c conda-forge numpy pandas pyarrow
pip install biopython
cd mutsim

# Set up/download reference genome
## GRCh38 & GENCODE v45
cd ~/mutsim/ref
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz
#--> Not working so I downloaded it locally and transferred to HPC cluster
gunzip GRCh38.primary_assembly.genome.fa.gz

module load SAMtools/1.16.1-GCC-11.3.0
samtools faidx GRCh38.primary_assembly.genome.fa

cd gtf
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
gunzip gencode.v45.annotation.gtf.gz


# Run + generate FASTA headers (gene first):
cd gtf

XXXY HERE!!!!


scripts/make_cds_bed.sh gencode.v45.annotation.gtf > cds.bed
bedtools getfasta -fi ../ref/GRCh38.primary_assembly.genome.fa \
         -bed cds.bed -s -name+ -fo cds.fa
Header format e.g.: >ENSG00000186092.7::chr1:65564-65573(+) (gene::coord).


```













