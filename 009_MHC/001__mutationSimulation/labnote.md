
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
pip install SigProfilerPlotting
pip install SigProfilerMatrixGenerator
pip install pandas
conda install numpy=1.24
```

```R
library("reticulate")
library("devtools")


install_github("AlexandrovLab/SigProfilerSimulatorR")

library("SigProfilerSimulatorR")
```

Its very buggy try this conda version instead

```bash
conda create -n sigprofiler -c bioconda -c conda-forge sigprofilerMatrixGenerator --yes

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
conda install conda-forge::parquet-tools # Added to vizualize parquet file
conda install bioconda::pysam # Was needed for build_context_index.py
conda install anaconda::pandas # Was needed for simulate_mutations.py
pip install biopython
pip install pandas
pip install SigProfilerMatrixGenerator
pip install seaborn
conda install bioconda::tabix # for tabix bug

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

#../scripts/make_cds_bed.sh gencode.v45.annotation.gtf > cds.bed
#--> Run this in interactive

awk '$3=="CDS"{OFS="\t";
     chr=$1; start=$4-1; end=$5; strand=$7;
     gene=$10; tx=$12;
     gsub(/[\";]/,"",gene); gsub(/[\";]/,"",tx);
     print chr, start, end, gene, tx, strand
}' gencode.v45.annotation.gtf > cds.bed

module load BEDTools/2.30.0-GCC-11.3.0
bedtools getfasta -fi ../ref/GRCh38.primary_assembly.genome.fa \
         -bed cds.bed -s -name+ -fo cds.fa

# Header format e.g.: >ENSG00000186092.7::chr1:65564-65573(+) (gene::coord).

conda activate mutsim

sbatch scripts/build_parquet_array.slurm
#--> ALL GOOD

# To check file
parquet-tools show --head 5 parquet/chr1.parquet 
```

This `scripts/build_parquet_array.slurm `will:
- Parses headers with regex
- Builds per‑nucleotide rows: chr, pos, strand, ref_base, gene_id, codon_index, ref_codon, ref_aa
- Writes one Parquet per chromosome (no append) with Zstd compression.


For each bp of CDS, it will:
- Computes its genomic position and strand.
- Encodes the reference base (A, C, G, T, N → 0–4).
- Reconstructs the codon it belongs to.
- Uses the codon to determine the reference amino acid (or * if stop).
--> It output in Parquet file with columns: chr, pos, strand, ref_base, gene_id, codon_index, ref_codon, ref_aa


Then let's create **context‑index builder**; the 96 SBS possible (4 possible left bases × 6 substitutions × 4 possible right bases) context.

```bash
python scripts/build_context_index.py
```


The script read each bp of CDS, then collect its trinucleotide context (bp before and bp after); and store the information; for ex. bp at postiioon 100 it fetches 98-101=3bp sequence= ACG for ex. **Reference/middle base must be C or T**: 
- C>A, C>G, C>T
- T>A, T>C, T>G
-->  4 left bases × 6 substitutions × 4 right bases = 96 mutation contexts


## Simulate mutations

```bash
# Test
python scripts/simulate_mutations.py \
  --signatures signatures/COSMIC_v3.4_SBS_GRCh38.txt \
  --signature-name SBS5 \
  --n 1000 \
  --exome-dir parquet \
  --context-index indices/context96.pkl \
  --out parquet/test_SBS5_1k.parquet

python scripts/simulate_mutations.py \
  --signatures signatures/COSMIC_v3.4_SBS_GRCh38.txt \
  --signature-name SBS2 \
  --n 10000 \
  --exome-dir parquet \
  --context-index indices/context96.pkl \
  --out parquet/test_SBS2_10k.parquet

# Check output
parquet-tools show --head 5 parquet/test_SBS5_1k.parquet 
```



## Convert parquet to txt file for SigProfiler

Let's check that the mutation simulation worked, by generating profile mutation plot with SigProfiler:

- Convert Parquet file to a SigProfiler-compatible matrix
    - First convert .parquet into .txt file, followin example file [here](https://osf.io/s93d5/wiki/home/)



**Parquet to txt**

```bash
conda activate mutsim
python
```
```python
import pandas as pd
import numpy as np

# Load parquet
df = pd.read_parquet("parquet/test_SBS2_10k.parquet") # SBS5

# Decode ref_base integers to letters
INT2BASE = np.array(["A", "C", "G", "T", "N"])
df["ref_base_letter"] = INT2BASE[df["ref_base"]]

# Remove 'chr' prefix from chromosome names
df["chr_clean"] = df["chr"].str.replace("^chr", "", regex=True)

# Format to expected mutation text format
df_txt = pd.DataFrame({
    "Project": "Simu",
    "Sample": "SBS2_10k", # SBS5
    "ID": ".",
    "Genome": "GRCh38",
    "mut_type": "SNP",
    "chrom": df["chr_clean"],
    "pos_start": df["pos"],
    "pos_end": df["pos"],
    "ref": df["ref_base_letter"],
    "alt": df["alt_base"],
    "Type": "SOMATIC"
})

# Save to file
df_txt.to_csv("txt/input/test_SBS2_10k.txt", sep="\t", index=False) # SBS5
```




### install SigProfilerMatrixGenerator 

```bash
conda activate mutsim
python
```

```python
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', rsync=True, bash=True)
```
--> Download ref genome fail; so I found them [here](hhttps://osf.io/s93d5/wiki/1.%20Installation%20-%20Python/); need download: chromosomes, matrix, vcf_files

So I `git clone https://github.com/AlexandrovLab/SigProfilerMatrixGenerator.git` cloned the whole repository, and then genome are here: `/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/references`

This is [where](https://github.com/AlexandrovLab/SigProfilerExtractor/issues/16) to put it..

But i do not find it online!! Not present in the repo; so instead I opened the VSC code code below and ran:

```bash
pip install SigProfilerMatrixGenerator SigProfilerSimulator
```
Then type `python`
```python
from SigProfilerMatrixGenerator.install import install as genInstall
genInstall("GRCh38", rsync=True, bash=True)
# --> Work!! Then below is to find where it has install it
import os

for root, dirs, files in os.walk(os.path.expanduser("~")):
    for name in dirs:
        if name == "GRCh38":
            print("FOUND:", os.path.join(root, name))

#--> C:\Users\roule\AppData\Local\Programs\Python\Python310\Lib\site-packages\SigProfilerMatrixGenerator\references\chromosomes ; in hidden files!
```
--> But then not sure what to do with these files, so instead I follow what was recommended [here](https://github.com/AlexandrovLab/SigProfilerExtractor/issues/16)

I open Git bash on my computer and cd `009*/001*` and then:

```bash
curl -O ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerMatrixGenerator/GRCh38.tar.gz
#--> Downloaded to current working dir; that I then transfer to 

python -c "import SigProfilerMatrixGenerator as sp; print(sp.__file__)" # Say where my SigProfiler is installed
```
--> This downloaded to my current local work dir `C:\Users\roule\OneDrive\Bureau\Github\Akiz_Lab` and tehn I copy to the clsuter at `/home/roulet/anaconda3/envs/mutsim/lib/python3.11/site-packages/SigProfilerMatrixGenerator/references/chromosomes/tsb/` and unzip file with `tar -xzf`



## Generate SigProfiler matrix


```bash
conda activate mutsim

SigProfilerMatrixGenerator matrix_generator test GRCh38 txt --plot=TRUE
```


--> The msimulation does not work properly, did not reproduce SBS2 mutation... Lets write anothoer `simulate_mutation.py` script


## Fine tune simulate_mutation.py scripts

### Change simulate_mutations script


```bash
nano scripts/simulate_mutations_v2.py



python scripts/simulate_mutations_v2.py \
  --signatures signatures/COSMIC_v3.4_SBS_GRCh38.txt \
  --signature-name SBS5 \
  --n 1000 \
  --exome-dir parquet \
  --context-index indices/context96.pkl \
  --out parquet/testv2_SBS5_1k.parquet

# Check output
parquet-tools show --head 5 parquet/testv2_SBS1_1k.parquet 

parquet-tools show --head 5 parquet/test_SBS2_1k.parquet 

```


#### Parquet to txt

```bash
conda activate mutsim
python
```
```python
import pandas as pd
import numpy as np

# Load the V2-style simulated mutations
df = pd.read_parquet("parquet/testv2_SBS5_1k.parquet")  # CHANGE NAME HERE!!!!!!!!!!!!!!!!!!!!!!

# Extract ref base from context ID
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref, alts in [("C", "AGT"), ("T", "ACG")]
    for alt in alts
    for l in "ACGT"
    for r in "ACGT"
]
contexts = np.array(CONTEXTS_96)
df["ref_base"] = [c[2] for c in contexts[df["context_id"]]]  # extract REF from e.g., T[C>T]A

# Remove 'chr' prefix
df["chr_clean"] = df["chr"].str.replace("^chr", "", regex=True)

# Format correctly
df_txt = pd.DataFrame({
    "Project": "Simu",
    "Sample": "SBS5_1k_v2",  # CHANGE NAME HERE!!!!!!!!!!!!!!!!!!!!!!
    "ID": ".",
    "Genome": "GRCh38",
    "mut_type": "SNP",
    "chrom": df["chr_clean"],
    "pos_start": df["pos"],
    "pos_end": df["pos"],
    "ref": df["ref_base"],
    "alt": df["alt_base"],
    "Type": "SOMATIC"
})

# Save to tab-delimited text file
df_txt.to_csv("txt/input/testv2_SBS5_1k.txt", sep="\t", index=False)   # CHANGE NAME HERE!!!!!!!!!!!!!!!!!!!!!!

```

#### Matrix and plot

```bash
conda activate mutsim

SigProfilerMatrixGenerator matrix_generator test GRCh38 txt --plot=TRUE
```



## Simulation with annotation and plots - testing

Let's collect the annoation information from [dbNSFP](https://www.dbnsfp.org/download) --> Set up an academic account and download `dbNSFP5.1a_grch38.gz` and `dbNSFP5.1a_grch38.gz.tbi` from the GoogleDrive linke received. Then transferred these files to `/ref`


Let's **simulate+annotate** using `simulate_sbs5_array.py` with 
--> I modified a bit `scripts/simulate_sbs5_array.py` and added entire path to folders `ref/` `signature/`... --> Renamed as `scripts/simulate_array_v2.py`; I also now use `--mem=100G` memory instead 36


```bash
conda activate mutsim

sbatch scripts/submit_simulation_array_v2.sh # 44062009/44063127/44063359 tabix installation error; 44063761 cancel; 

# Test on SBS5 and SBS90
sbatch scripts/submit_simulation_array_v2-SBS90.sh # 44064135  xxx ~20min
sbatch scripts/submit_simulation_array_v2-SBS5.sh # 44064160 xxx

```

- *NOTE: Here `scripts/submit_simulation_array_v2.sh` is simply `scripts/submit_simulation_array.sh` but with `scripts/simulate_array_v2.py`*

There is 80 array jobs in the `scripts/submit_simulation_array_v2.sh` as *8 mutation counts × 10 replicates = 80 combinations → array indices 0 to 79*

Then **plots to summarise simmulation**:

--> I modified a bit `scripts/summarize_simulation_results.py`; to deal with a bug with matching file name.

```bash
# summary metrics
python scripts/summarize_simulation_results_v2.py \
 --results-dir results/SBS90 \
 --output results/SBS90/sbs90_summary_all.tsv
python scripts/summarize_simulation_results_v2.py \
 --results-dir results/SBS5 \
 --output results/SBS5/sbs5_summary_all.tsv
# plot
python scripts/plot_signature_summary.py \
  --input results/SBS90/sbs90_summary_all.tsv \
  --output results/SBS90/sbs90_summary_plots.pdf
python scripts/plot_signature_summary.py \
  --input results/SBS5/sbs5_summary_all.tsv \
  --output results/SBS5/sbs5_summary_plots.pdf
```


--> The test seems good, results are the same as the one obtained by Joan. 


## Simulation with annotation and plots - Run all mutation signature


Let's generate a custom code to run all mutation signature in a single script; script will:
- Extracts all signature names (ie. SBS1, SBS5, SBS7a, etc.) from `signatures/COSMIC_v3.4_SBS_GRCh38.txt` 
- Loops over them
- Writes individual SLURM batch scripts (one per signature) to `scripts/submit_all_signatures`
- Then run run them all `for f in scripts/submit_all_signatures/*.sh; do sbatch "$f"; done`
- Run a SLURM batch that will generate summary and plots for all signatures



```bash
conda activate mutsim

# Generate all unique scripts for each mutation automatically
bash scripts/GenerateScript-submit_simulation_array_v2.sh
#--> All scripts generated at scripts/submit_all_signatures

for f in scripts/submit_all_signatures/*.sh; do sbatch "$f"; done # XXX

XXXY BELOW NOT RUN!!!

# Generate summary metrics and plots for each signature
bash scripts/run_summary_plot_all.sh



```


