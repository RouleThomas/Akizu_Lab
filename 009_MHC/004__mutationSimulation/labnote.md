# Project

Simulation of Experimental and SBS signatures.

Previous version `003` got issues, again related to the pkl generation.. Let's start with the 001 version this time, and try correct the issue related to .pkl generation and pyrimidine...





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
  --signature "SBS1" \
  --n "500" \
  --rep "1" \
  --seed "42" \
  --sigfile signatures/COSMIC_with_flat.txt \
  --outdir results # UPDATE PATH





sbatch scripts/run_filtered_cosmic.slurm # 





sbatch scripts/run_filtered_experimental_1.slurm # WAIT CHECK COSMIC FIRST xxx --> results_experimental/
sbatch scripts/run_filtered_contexts_1.slurm #  xxx --> results_contexts/


```

--> XXX All good all files generated with `n_*` folders in the respectrive `results*/` folders


XXXY SUFFERING HERE!




## Generate profile plots


Let's generate profile plots of my simulation to see if these are in agreement with the known profile.

- Convert simulation `.parquet` to `.txt`
- Transfer to `plot/` folder a few .txt files
- Run `SigProfilerMatrixGenerator` to generate signature profile plot 








```bash
parquet_to_sigprofiler_txt.py


python scripts/parquet_to_sigprofiler_txt.py \
  --parquet results/SBS1/n_500/rep_01.annot.parquet \
  --fasta ref/GRCh38.primary_assembly.genome.fa \
  --sample-name SBS1-n_500-rep_01 \
  --output plot/SBS1-n_500-rep_01.txt




```







```bash
SigProfilerMatrixGenerator matrix_generator test GRCh38 plot/input --plot=TRUE


```





