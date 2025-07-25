# Project

Simulation of Experimental and SBS signatures.

Previous version `003` got issues, again related to the pkl generation.. Let's start with the 001 version this time, and try correct the issue related to .pkl generation and pyrimidine...


# Regenerate cds.fa and .parquet files using clean reference



```bash
# Download FASTA
rm ref/GRCh38.primary_assembly.genome.fa
cd ref/
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa

# Generate bed of CDS
cd gtf/
awk '$3=="CDS"{OFS="\t";
     chr=$1; start=$4-1; end=$5; strand=$7;
     gene=$10; tx=$12;
     gsub(/[\";]/,"",gene); gsub(/[\";]/,"",tx);
     print chr, start, end, gene, tx, strand
}' gencode.v45.annotation.gtf > cds.bed

# Generate FASTA of CDS
conda activate BedToBigwig

grep -v "_" cds.bed > cds.primary.bed # keep only main chr
bedtools getfasta \
    -fi ../ref/GRCh38.primary_assembly.genome.fa \
    -bed cds.primary.bed -s -name+ -fo cds.fa

# build parquet
conda activate mutsim

sbatch scripts/build_parquet_array.slurm # 47779545 ok
```




# Build the pkl 

Let's first build the *pkl* (label each bp that are in one of the 96 mutation context) from the *parquet* (File with chr, pos, strand, ref_base, gene_id, codon_index, ref_codon, ref_aa for each bp of CDS):
- Parquet files already generated in `001` and has been copied in `parquet/`


```bash
conda activate mutsim
# To check file
parquet-tools show --head 5 parquet/chr1.parquet 

# Build pkl
python scripts/build_context_index.py # ok
```

--> 195,291,599 total indices
--> 126,870,005 positions


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

print(f"\n📦 Total positions in parquet: {total_rows:,}") # 126,870,005
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

print(f"🧠 Total indices stored across all 96 contexts: {len(all_indices):,}") # 34,210,797
print(f"🔍 Unique positions covered (non-redundant rows): {len(unique_indices):,}") # 126,870,005
```

--> 126,870,005 positions (= all CDS bp) in both parquet and pkl index. 





# Simulate mutation


--> Update `scripts/simulate_mutations.py` based on `001__`

--> Update `scripts/simulate_array.py` based on `001__`

Three scripts `run_filtered_*.slurm`; each uses `scripts/simulate_array.py` for simulation and annotation:
- cosmic= All SBS signatures (without the artifcatual ones)
- experimental_serial= experimental signatures
- contexts=random control



```bash
conda activate mutsim

python scripts/simulate_array.py \
  --signature "SBS2" \
  --n "4000" \
  --rep "1" \
  --seed "42" \
  --sigfile signatures/COSMIC_with_flat.txt \
  --outdir results # UPDATE PATH
#--> light testing seems to be working!


# Generate plot for all
sbatch scripts/run_filtered_cosmic.slurm # 47828611 ok --> results/
sbatch scripts/run_filtered_experimental.slurm # 47875937 xxx --> results_experimental/
sbatch scripts/run_filtered_contexts.slurm # 47876020 xxx --> results_contexts/
```

--> All good all files generated with `n_*` folders in the respectrive `results*/` folders






## Generate profile plots


Let's generate profile plots of my simulation to see if these are in agreement with the known profile.

- Custom script `plot_sbs96_from_parquet.py` to make these plot from the `*.parquet`


```bash

plot_sbs96_from_parquet.py

# Check a few samples

python scripts/plot_sbs96_from_parquet.py \
  --parquet results/SBS2/n_4000/rep_01.annot.parquet \
  --fasta ref/GRCh38.primary_assembly.genome.fa \
  --sample-name SBS2-n_4000-rep_01 \
  --output-pdf plot/SBS2-n_4000-rep_01_plot.pdf
```

--> It seems to be working!!!



## Annotation summary plot

Let's generate plot with `x` = `n_mutations` and `y` = `annotation score`



```bash
conda activate mutsim

####################################
# COSMIC ##########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-cosmic.sh

# One plot with all mutations
python scripts/plot_all_signatures_combined-cosmic.py
#--> results/combined_signature_summary_errorbars.pdf

# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-cosmic-highlight_random.py
#--> results/combined_signature_summary_plots-highlight_random.pdf



####################################
# EXPERIMENTAL ##########################
####################################

XXX



####################################
# CONTEXT ##########################
####################################

XXX




```







# QC simulation


Let's double check that the simulation is working properly.
--> Check KRAS gene (ENSG00000133703 )

```python
import pandas as pd
import numpy as np

df = pd.read_parquet("results/SBS2/n_4000/rep_01.annot.parquet")
# Mapping: 0=A, 1=C, 2=G, 3=T
INT2BASE = np.array(["A", "C", "G", "T"])

# Replace if needed
df["ref_base"] = INT2BASE[df["ref_base"].values]


# inspect mutation
print(df[[
    "chr", "pos", "strand", "ref_base", "alt_base",
    "ref_codon", "mut_codon", "ref_aa", "alt_aa"
]])
```


--> I checked (on IGV):
  - ref and mut codon; all good
  - `+` `-` strand mutation correctly applied
  - mutation not always at the center of the codon











