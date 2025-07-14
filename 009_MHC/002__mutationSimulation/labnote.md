# Project

Simulation of Experimental and SBS signatures.

Previous version `001` got issues, notably we only generated mutations on the + strand pyrimidine bases (C and T):
- COSMIC, represent all mutations with Pyrimidine bases; *Even if the actual mutated base was a purin (G or A) on the - strand, it is reported as the complementary pyrimidine mutation on the + strand*
    - so actual G>T on - strand is reported as C>A on COSMIC... And in `001`, we only mutated pyrimidine
- Improve code speed
- Filter out artifactual COSMIC signatures
- New version dbNSFP5.2


# Download dbNSFP5.2

See CHOP email *dbNSFP academic user registration response* for credentials to use in [dbnsfp](https://www.dbnsfp.org/download) for already registered users.


```bash
wget https://download.genos.us/dbnsfp/academic/01f8c3/dbNSFP5.2a_grch38.gz
wget https://download.genos.us/dbnsfp/academic/01f8c3/dbNSFP5.2a_grch38.gz.tbi
wget https://download.genos.us/dbnsfp/academic/01f8c3/dbNSFP5.2a_grch38.gz.md5
```

--> All good





# Build the pkl 

Let's first build the *pkl* (label each bp that are in one of the 96 mutation context) from the *parquet* (File with chr, pos, strand, ref_base, gene_id, codon_index, ref_codon, ref_aa for each bp of CDS):
- Parquet files already generated in `001` and has been copied in `parquet/`


```bash
conda activate mutsim
# To check file
parquet-tools show --head 5 parquet/chr1.parquet 

# Build pkl
python scripts/build_context_index_3.py # ok
```

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

print(f"🧠 Total indices stored across all 96 contexts: {len(all_indices):,}") # 380,610,012
print(f"🔍 Unique positions covered (non-redundant rows): {len(unique_indices):,}") # 126,870,005
```

--> 126,870,005 positions (= all CDS bp) in both parquet and pkl index. 



# Simulate mutation

Three scripts `run_filtered_*.slurm`; each uses `simulate_sbs5_array.py` for simulation and annotation:
- cosmic= All SBS signatures (without the artifcatual ones)
- experimental_serial= experimental signatures
- contexts=random control


Let's re-name `scripts/simulate_sbs5_array.py` into `scripts/simulate_array.py`; and **update it so that it uses *dbNSFP5.2a* + path folders updated**


```bash
conda activate mutsim

sbatch scripts/run_filtered_cosmic_1.slurm # 46998992 ok  --> results/
sbatch scripts/run_filtered_experimental_1.slurm # 47058807 fail; 47174183 xxx
sbatch scripts/run_filtered_contexts_1.slurm #  47058818 ok  --> results_contexts/


```

--> XXX
