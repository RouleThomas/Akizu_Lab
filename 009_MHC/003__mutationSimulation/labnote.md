# Project

Simulation of Experimental and SBS signatures.

Previous version `002` got issues, again related to the pkl generation.. I realized it was trash...





# Build the pkl 

Let's first build the *pkl* (label each bp that are in one of the 96 mutation context) from the *parquet* (File with chr, pos, strand, ref_base, gene_id, codon_index, ref_codon, ref_aa for each bp of CDS):
- Parquet files already generated in `001` and has been copied in `parquet/`


```bash
conda activate mutsim
# To check file
parquet-tools show --head 5 parquet/chr1.parquet 

# Build pkl
python scripts/build_context_index_4.py # ok
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
print(f"🔍 Unique positions covered (non-redundant rows): {len(unique_indices):,}") # 126,870,004
```

--> 126,870,005 positions (= all CDS bp) in both parquet and pkl index. 





# Simulate mutation

--> Update `scripts/simulate_mutations.py` into `scripts/simulate_mutations_1.py`

--> Update `scripts/simulate_array.py` into `scripts/simulate_array_1.py`

Three scripts `run_filtered_*.slurm`; each uses `scripts/simulate_array.py` for simulation and annotation:
- cosmic= All SBS signatures (without the artifcatual ones)
- experimental_serial= experimental signatures
- contexts=random control



```bash
conda activate mutsim

sbatch scripts/run_filtered_cosmic_1.slurm # 47244234 fail --> results/
# --> NEW scripts: scripts/simulate_mutations_1.py scripts/simulate_array_1.py
sbatch scripts/run_filtered_cosmic_2.slurm # 47402117 xxx --> results/





sbatch scripts/run_filtered_experimental_1.slurm # WAIT CHECK COSMIC FIRST xxx --> results_experimental/
sbatch scripts/run_filtered_contexts_1.slurm #  xxx --> results_contexts/


```

--> XXX All good all files generated with `n_*` folders in the respectrive `results*/` folders



XXXY HERE CHECK MAIN IN CHATGPT...




## Generate profile plots


Let's generate profile plots of my simulation to see if these are in agreement with the known profile.

- Convert simulation `.parquet` to `.txt`
- Transfer to `plot/` folder a few .txt files
- Run `SigProfilerMatrixGenerator` to generate signature profile plot 








```bash
parquet_to_sigprofiler_txt.py


python scripts/parquet_to_sigprofiler_txt.py \
  --parquet results/SBS1/n_1000/rep_01.annot.parquet \
  --fasta ref/GRCh38.primary_assembly.genome.fa \
  --sample-name SBS1-n_1000-rep_01 \
  --output plot/SBS1-n_1000-rep_01.txt



```







```bash
conda activate mutsim
python
```
```python
import pandas as pd
import numpy as np
import pysam






############################################
# Load from results ######################
############################################
import pandas as pd
import pysam
import numpy as np

# === INPUTS ===
parquet_path = "results/SBS1/n_1000/rep_01.annot.parquet"  # CHANGE THIS
fasta_path = "ref/GRCh38.primary_assembly.genome.fa"       # CHANGE THIS
sample_name = "SBS1-n_1000-rep_01"                          # CHANGE THIS
output_txt = f"plot/{sample_name}.txt"

# === LOAD ===
df = pd.read_parquet(parquet_path)
fasta = pysam.FastaFile(fasta_path)

# === Normalize chromosome name (e.g. 'chr1' → '1') ===
df["chr_clean"] = df["chr"].str.replace("^chr", "", regex=True)

# === Get genome REF base (+ strand) ===
def fetch_genome_ref(row):
    try:
        return fasta.fetch(row["chr_clean"], row["pos"] - 1, row["pos"]).upper()
    except Exception:
        return "N"

df["genome_ref"] = df.apply(fetch_genome_ref, axis=1)

# === Decode simulated REF base ===
INT2BASE = {0: "A", 1: "C", 2: "G", 3: "T"}
df["simulated_ref"] = df["ref_base"].map(INT2BASE)

# === Check which entries match the genome ===
df["match"] = df["genome_ref"] == df["simulated_ref"]

# === Report ===
n_total = len(df)
n_match = df["match"].sum()
print(f"✅ Matches: {n_match} / {n_total} ({n_match/n_total:.1%})")

# === Filter only valid entries ===
df_valid = df[df["match"] & (df["genome_ref"] != "N")].copy()

# === Format final DataFrame for SigProfiler ===
df_txt = pd.DataFrame({
    "Project": "Simu",
    "Sample": sample_name,
    "ID": ".",
    "Genome": "GRCh38",
    "mut_type": "SNP",
    "chrom": df_valid["chr_clean"],
    "pos_start": df_valid["pos"],
    "pos_end": df_valid["pos"],
    "ref": df_valid["genome_ref"],
    "alt": df_valid["alt_base"],
    "Type": "SOMATIC"
})

# === Save ===
df_txt.to_csv(output_txt, sep="\t", index=False)
print(f"✅ Saved: {output_txt} with {len(df_txt)} mutations")


############################################
# Load from results_experimental ######################
############################################



```



```bash
SigProfilerMatrixGenerator matrix_generator test GRCh38 plot/input --plot=TRUE


```


