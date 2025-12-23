# Project

Simulation of Experimental and SBS signatures.

Previous version `004` got issues:
- dbNSFP, many score were missing; likely because of mismatch GENCODE use and GENCODE version used by dbNSFP;
  - now lets use: **dbNSFP v5.3 with GENCODE release 49 (Ensembl release 115, September 2025)**
- chromosme parquet not correct; phase exon junction not taken into account
  - Re-generate chromosome parquet correctly
- Work with 1-based



# Download dbNSFP5.3a

See CHOP email *dbNSFP academic user registration response* for credentials to use in [dbnsfp](https://www.dbnsfp.org/download) for already registered users.



```bash
cd ref/

wget https://dist.genos.us/academic/b2fd38/dbNSFP5.3a_grch38.gz
wget https://dist.genos.us/academic/b2fd38/dbNSFP5.3a_grch38.gz.tbi
wget https://dist.genos.us/academic/b2fd38/dbNSFP5.3a_grch38.gz.md5
```

--> All good


# Regenerate cds.fa and .parquet files using clean reference

Use (GENCODE release 49)[https://www.gencodegenes.org/human/release_49.html]:
- *GTF*: Comprehensive gene annotation, CHR, It contains the comprehensive gene annotation on the reference chromosomes only
- *FASTA*: Genome sequence, primary assembly (GRCh38), PRI,Nucleotide sequence of the GRCh38 primary genome assembly (chromosomes and scaffolds)
    The sequence region names are the same as in the GTF/GFF3 files






```bash
# Download FASTA
cd ref/
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa

# Generate bed of CDS
cd gtf/
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gunzip gencode.v49.annotation.gtf.gz

awk '$3=="CDS"{OFS="\t";
     chr=$1; start=$4-1; end=$5; strand=$7;
     gene=$10; tx=$12;
     gsub(/[\";]/,"",gene); gsub(/[\";]/,"",tx);
     print chr, start, end, gene, tx, strand
}' gencode.v49.annotation.gtf > cds.bed

# Generate FASTA of CDS
conda activate BedToBigwig

grep -v "_" cds.bed > cds.primary.bed # keep only main chr
bedtools getfasta \
    -fi ../ref/GRCh38.primary_assembly.genome.fa \
    -bed cds.primary.bed -s -name+ -fo cds.fa

# build parquet
conda activate mutsim

sbatch scripts/build_parquet_array.slurm # 61905518 NOT GOOD, phase not taken into account

sbatch scripts/build_parquet_array_version2.slurm # 62161823 ok



# To check file
parquet-tools show --head 5 parquet/chr1.parquet 



```


--> All good:
- parquet are now 1-based
- ref_AA and ref_codon are correct! Even after exon exon junction









# Build the pkl 




```bash
conda activate mutsim

# Build pkl
python scripts/build_context_index_version2.py #  



```

--> 653,843,189 total indices
--> 324,923,277 positions


## Check that number of positions match btwn parquet and pkl


Count **rows across all `*.parquet` files**


```python
import pyarrow.parquet as pq
from pathlib import Path

PARQ_DIR = Path("parquet")  

total_rows = 0
for parquet in sorted(PARQ_DIR.glob("chr*.parquet")):
    table = pq.read_table(parquet, columns=["chr"])  # read minimal data
    num_rows = table.num_rows
    total_rows += num_rows
    print(f"{parquet.name}: {num_rows:,} rows")

print(f"\n📦 Total positions in parquet: {total_rows:,}") # 324,923,277
```


Count **rows stored in `context96.pkl`** (number of unique indices)


```python
import pickle
import numpy as np
from pathlib import Path

with open("indices/context96.pkl", "rb") as f:
    context_index = pickle.load(f)

# Flatten all indices into one array
all_indices = np.concatenate(list(context_index.values()))
unique_indices = np.unique(all_indices)

print(f"🧠 Total indices stored across all 96 contexts: {len(all_indices):,}") # 495,018,673
print(f"🔍 Unique positions covered (non-redundant rows): {len(unique_indices):,}") # 324,923,277
```

--> 324,923,277 positions (= all CDS bp) in both parquet and pkl index. 




# Simulate mutation  - SBS, Exp, context



Here follow `## version4 update for missing dbNSFP5` from `009*/004*`.

- Update `scripts/simulate_array_v4.py` into `scripts/simulate_array_v4.py` to use *dbNSFP55.3* and added *stop_lost* count + work with 1-based
- Update  CADD PHRED value score also; `scripts/annotate_damage4.py` into `scripts/annotate_damage5.py`
- Update `scripts/simulate_mutations.py` into `scripts/simulate_mutations_v2.py`

```python
import gzip

with gzip.open("ref/dbNSFP5.3a_grch38.gz", "rt") as f:
    header = f.readline().strip().split("\t")

for i, col in enumerate(header, 1):
    print(f"{i}: {col}")
```

--> Double check we are good (by checking `scripts/annotate_damage3.py`):
  -  IDX_TRANSCRIPTID     = 15 - 1 --> OK
  -  IDX_VEP_CANONICAL    = 26 - 1 --> OK
  -  IDX_SIFT4G_SCORE     = 50 - 1 --> OK
  -  IDX_SIFT4G_PRED      = 52 - 1 --> OK
  -  IDX_PP2_HDIV_SCORE   = 53 - 1 --> OK
  -  IDX_PP2_HDIV_PRED    = 55 - 1 --> OK
  -  IDX_CADD_PHRED       = 150 - 1 --> NOT OK!!! Correct is 150!!!




```bash
conda activate mutsim

# Test with one simulation

python scripts/simulate_array_v5.py \
  --signature SBS90 \
  --n 4000 \
  --rep 1 \
  --seed 42 \
  --sigfile signatures/COSMIC_with_flat.txt \
  --outdir results_v5 # UPDATE PATH

python scripts/simulate_array_v5.py \
  --signature SBS1 \
  --n 4000 \
  --rep 1 \
  --seed 42 \
  --sigfile signatures/COSMIC_with_flat.txt \
  --outdir results_v5 # UPDATE PATH



XXX TO RUN BELOW!!! 

# Generate plot for all 

sbatch scripts/run_filtered_cosmic_v5.slurm # 62192890 xxx --> results_v5/
sbatch scripts/run_filtered_contexts_v5.slurm # 62192892 xxx --> results_contexts_v5/
sbatch scripts/run_filtered_experimental_v5.slurm # 62192893 xxx --> results_experimental_v5/
```







# QC simulation
## Check why so many missing in dbNSFP5 score


Let's use audit script from `009*/004*` 





```bash


# CHECK CADD PHRED - version4 0-based 1-based corrected
python scripts/audit_cadd_phred.py ref/dbNSFP5.3a_grch38.gz results_v5/SBS90/n_4000/rep_01.annot.parquet
python scripts/audit_cadd_phred.py ref/dbNSFP5.3a_grch38.gz results_v5/SBS1/n_4000/rep_01.annot.parquet

```

--> All good, I check on dbNSFP and the missing are really missing!





XXXY HERE !!!