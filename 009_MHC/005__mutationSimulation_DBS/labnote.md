# Project

Simulation DBS signature; lets follow `004*` version; but many changes will be needed; we will also need a context version for DBS!

--> CDS and parquet of genome file can still be used: `../004__mutationSimulation/parquet/chr*.parquet`
    --> Then the pkl need this time to follow a 78 context size (ie. not 96 as for SBS)






# Build the pkl 

Let's first build the *pkl* (label each bp that are in one of the 78 DBS mutation context) from the *parquet* (File with chr, pos, strand, ref_base, gene_id, codon_index, ref_codon, ref_aa for each bp of CDS):
- Parquet files already generated in `004` and has been copied in `parquet/`


```bash
conda activate mutsim
# To check file
parquet-tools show --head 5 parquet/chr1.parquet 

XXXY HERE WORK ON THIS check chat

# Build pkl
python scripts/build_context78_index.py #  xxx
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




