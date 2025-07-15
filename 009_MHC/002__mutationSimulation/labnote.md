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

sbatch scripts/run_filtered_cosmic_1.slurm # 46998992 ok --> results/
sbatch scripts/run_filtered_experimental_1.slurm # 47058807 fail; 47174183 ok --> results_experimental/
sbatch scripts/run_filtered_contexts_1.slurm #  47058818 ok --> results_contexts/


```

--> All good all files generated with `n_*` folders in the respectrive `results*/` folders



# Verify that the simulation were perform correctly

## Check individual mutations and dbnsfp

Import simulation in python and manually inspect whether STOP is STOP, missense is missense; using manual query in [dbnsfp](https://www.dbnsfp.org/web-query).



```python
import pandas as pd
df = pd.read_parquet("results/SBS4/n_16000/rep_01.annot.parquet")
df[['chr', 'pos', 'ref_base', 'alt_base', 'consequence']].sample(5) # check a few rows
```

--> SBS4  rep1 n 16k; first five rows:

| Row | Original             | Formatted Variant ID |
| --- | -------------------- | -------------------- |
| 1   | `chr3 184321295 2 A` | `3-184321295-G-A`    |
| 2   | `chr18 58939494 0 A` | `18-58939494-A-A`    |
| 3   | `chr14 71661412 0 G` | `14-71661412-A-G`    |
| 4   | `chr14 64404139 1 T` | `14-64404139-C-T`    |
| 5   | `chr17 50618101 2 A` | `17-50618101-G-A`    |

--> Pasted formatted variant ID in in [dbnsfp](https://www.dbnsfp.org/web-query).

--> Issue here is that these variants does not exist in dbnsfp; not annotated as they are new! **Let's filter to keep the one that already exist in dbnsfp**



```python
import pandas as pd
from pathlib import Path
import sys
sys.path.append("scripts") 
from annotate_damage3 import DBNSFP


# Load Parquet simulation
df = pd.read_parquet("results/SBS4/n_500/rep_01.annot.parquet")

# Convert numeric base to letter
BASES = ["A", "C", "G", "T"]
df["ref_nt"] = df["ref_base"].map(lambda i: BASES[i])
df["alt_nt"] = df["alt_base"].astype(str)

# Load dbNSFP
db = DBNSFP(Path("ref/dbNSFP5.2a_grch38.gz"))

# Check existence in dbNSFP
hits = []
for idx, row in df.iterrows():
    result = db.query(str(row["chr"]).replace("chr", ""), int(row["pos"]), row["ref_nt"], row["alt_nt"])
    hits.append(result is not None and result["sift4g_score"] is not None)

# Add result column
df["in_dbNSFP"] = hits

# Summary - output rows with match with dbNSFP
print(df["in_dbNSFP"].value_counts())
print(df[df["in_dbNSFP"]].head())


# check specific rows where there is a match in dbNSFP
df[['chr', 'pos', 'ref_base', 'alt_base', 'consequence','sift4g_score', 'polyphen2_hdiv_score','cadd_phred','sift4g_pred']].iloc[4]

# Inspect my local dbNSFP to check constitency with simulation
result = db.query("11", 64651633, "T", "A")
print(result)
result = db.query("9", 83970201, "T", "A")
print(result)
result = db.query("19", 6697388, "G", "A")
print(result)
result = db.query("12", 53443571, "C", "A")
print(result)
result = db.query("12", 69257848, "G", "A")
print(result)

```

--> SBS4  rep1 n 500; first five rows with dbNSFP hit:

| Row | Original             | Formatted Variant ID |
| --- | -------------------- | -------------------- |
| 1   | `chr11 64651633 3 A`  | `1-64651633-T-A`     |
| 2   | `chr9 83970201 3 A`  | `9-83970201-T-A`     |
| 3   | `chr19 6697388 2 A`  | `19-6697388-G-A`     |
| 4   | `chr12 53443571 1 A` | `12-53443571-C-A`    |
| 5   | `chr12 69257848 2 A` | `12-69257848-G-A`    |



--> Looks good, the score are the same between the known variant from dbNSFP and the one I simulated.



## Check KRAS gene




```python
import pyarrow.parquet as pq
import pandas as pd

# Adjust path
df = pq.read_table("parquet/chr12.parquet").to_pandas()  # KRAS is on chr12
df['gene_id'] = df['gene_id'].str.replace(r'\.\d+$', '', regex=True) # remove gene version
# Check how many bases are associated with KRAS
kras_rows = df[df['gene_id'] == "ENSG00000133703"]
print(kras_rows.shape)
print(kras_rows[['pos', 'strand', 'codon_index', 'ref_codon', 'ref_aa']].drop_duplicates())


# Glycine codons: GGT, GGC, GGA, GGG
gly_codon_set = {"GGT", "GGC", "GGA", "GGG"}

# Which codons code glycine in KRAS
gly_kras_codons = kras_rows[kras_rows['ref_codon'].isin(gly_codon_set)]
print(gly_kras_codons[['pos', 'strand', 'ref_codon', 'ref_aa']].drop_duplicates())
print(f"Total glycine codons in KRAS: {gly_kras_codons['codon_index'].nunique()}")
```



--> Not sure what to conclude from that part...





## Generate profile plots


Let's generate profile plots of my simulation to see if these are in agreement with the known profile.

- Convert simulation `.parquet` to `.txt`
- Transfer to `plot/` folder a few .txt files
- Run `SigProfilerMatrixGenerator` to generate signature profile plot 


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
# Load your annotated Parquet file
df = pd.read_parquet("results/SBS6/n_4000/rep_01.annot.parquet")
# Load the reference genome
fasta = pysam.FastaFile("ref/GRCh38.primary_assembly.genome.fa")  # <-- UPDATE THIS
# Clean chromosome name
df["chr_clean"] = df["chr"].str.replace("^chr", "", regex=True)

# Fetch the actual REF base from genome (+ strand)
def fetch_ref(row):
    try:
        return fasta.fetch(row["chr_clean"], row["pos"] - 1, row["pos"]).upper()
    except Exception:
        return "N"

df["ref"] = df.apply(fetch_ref, axis=1)

# Now compare to df["alt_base"] and construct the final .txt
df_txt = pd.DataFrame({
    "Project": "Simu",
    "Sample": "SBS6-n_4000-rep_01",
    "ID": ".",
    "Genome": "GRCh38",
    "mut_type": "SNP",
    "chrom": df["chr_clean"],
    "pos_start": df["pos"],
    "pos_end": df["pos"],
    "ref": df["ref"],
    "alt": df["alt_base"],
    "Type": "SOMATIC"
})

# Filter out any entries with ambiguous base
df_txt = df_txt[df_txt["ref"] != "N"]
df_txt.to_csv("plot/SBS6-n_4000-rep_01.txt", sep="\t", index=False)






# Get ALT base from simulation (integer-coded)
INT2BASE = {0: "A", 1: "C", 2: "G", 3: "T"}
df["alt"] = df["alt_base"].map(INT2BASE)
# Filter out where the reference base does not match the genome
mismatch_mask = df["ref_base"].map(INT2BASE) != df["ref"]
print(f"⚠️ Filtering out {mismatch_mask.sum()} mismatches between simulated and real genome bases.")
df = df[~mismatch_mask]
# Create SigProfiler-compatible .txt
df_txt = pd.DataFrame({
    "Project": "Simu",
    "Sample": "SBS6-n_4000-rep_01",
    "ID": ".",
    "Genome": "GRCh38",
    "mut_type": "SNP",
    "chrom": df["chr_clean"],
    "pos_start": df["pos"],
    "pos_end": df["pos"],
    "ref": df["ref"],
    "alt": df["alt"],
    "Type": "SOMATIC"
})
# Save
df_txt.to_csv("plot/SBS6-n_4000-rep_01.txt", sep="\t", index=False)







############################################
# Load from results_experimental ######################
############################################
df = pd.read_parquet("results_experimental/methanol_5f193a12cbbf/n_4000/rep_01.annot.parquet")  # CHANGE NAME HERE!!!!!!!!!!!!!!!!!!!!!!
# Extract ref base from context ID
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref, alts in [("C", "AGT"), ("T", "ACG")]
    for alt in alts
    for l in "ACGT"
    for r in "ACGT"
]
contexts = np.array(CONTEXTS_96)
df["ref_base"] = df["ref_base"].map({0: "A", 1: "C", 2: "G", 3: "T"})
# Remove 'chr' prefix
df["chr_clean"] = df["chr"].str.replace("^chr", "", regex=True)
# Format correctly
df_txt = pd.DataFrame({
    "Project": "Simu",
    "Sample": "methanol_5f193a12cbbf-n_4000-rep_01",  # CHANGE NAME HERE!!!!!!!!!!!!!!!!!!!!!!
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
df_txt.to_csv("plot/methanol_5f193a12cbbf-n_4000-rep_01.txt", sep="\t", index=False)   # CHANGE NAME HERE!!!!!!!!!!!!!!!!!!!!!!


```



```bash
SigProfilerMatrixGenerator matrix_generator test GRCh38 plot --plot=TRUE


```





--> HUGE ISSUE WHEN GENERATING THE PKL... ISSUE RELATED TO THE PYRIMIDINE....

--> Lets do a new version `003__mutationSimulation`




