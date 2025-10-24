# Project

Simulation DBS signature; lets follow `004*` version; but many changes will be needed; we will also need a context version for DBS!

--> CDS and parquet of genome file can still be used: `../004__mutationSimulation/parquet/chr*.parquet`
    --> Then the pkl need this time to follow a 78 context size (ie. not 96 as for SBS)






# Build the pkl 

Let's first build the *pkl* (label each bp that are in one of the 78 DBS mutation context) from the *parquet* (File with chr, pos, strand, ref_base, gene_id, codon_index, ref_codon, ref_aa for each bp of CDS):
- Parquet files already generated in `004` and has been copied in `parquet/`

--> Tricky part is to build the `CONTEXTS_78`; for that we need to find all 




```bash
conda activate mutsim
# To check file
parquet-tools show --head 5 parquet/chr1.parquet 


# Build pkl
python scripts/build_context78_index.py # work but issue as i force all ref to start with C or T; so fail at simulation


python scripts/build_context78_index_v2.py # OK!
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


Count **rows stored in `context78.pkl`** (number of unique indices)


```python
import pickle
import numpy as np
from pathlib import Path

with open("indices/context78.pkl", "rb") as f:
    context_index = pickle.load(f)

# Flatten all indices into one array
all_indices = np.concatenate(list(context_index.values()))
unique_indices = np.unique(all_indices)

print(f"🧠 Total indices stored across all 78 contexts: {len(all_indices):,}") # 2,967,684
print(f"🔍 Unique positions covered (non-redundant rows): {len(unique_indices):,}") # 126,870,005
```

--> 126,870,005 positions (= all CDS bp) in both parquet and pkl index. 


# Prerequiste files

- COSMIC DBS: Download data from COSMIC= `COSMIC_v3.4_DBS_GRCh38.txt`
- COSMIC DBS with Flat: Add Flat; which is for SBS `1/96=0.0104167` so for DBS it is: `1/78= 0.01282051`
    --> Simply copy COSMIC_v3.4_DBS_GRCh38.txt and add Flat with value 0.01282051

```bash
# Add Flat to COSMIC (1/78 chance to have this signature)
awk -F'\t' -v OFS='\t' '
NR==1 { print $0, "Flat"; next }   # add new column name to header
{ print $0, "0.01282051" }        # append value to each data row
' signatures/COSMIC_v3.4_DBS_GRCh38.txt > tmp && mv tmp signatures/COSMIC_v3.4_DBS_GRCh38_with_flat.txt
```

- `signatures/context_signature_list_DBS.txt`= just a list of all DBS context one per row (generate manually)
- Generate `signatures/context_sigs_fixed_DBS.txt` file which is our context version for DBS


```bash
awk '
NR==FNR {a[NR]=$1; n=NR; next}
END {
  printf "type"
  for (i=1;i<=n;i++) printf "\t" a[i]
  printf "\n"
  for (i=1;i<=n;i++) {
    printf a[i]
    for (j=1;j<=n;j++) {
      if (i==j) printf "\t1.0"; else printf "\t0.0"
    }
    printf "\n"
  }
}' signatures/context_signature_list_DBS.txt signatures/context_signature_list_DBS.txt \
> signatures/context_sigs_fixed_DBS.txt
```

- Generate `signatures/DBS_COSMIC_signatures.txt` = list of all DBS names (generate manually)


# Simulate mutation  - DBS, context

Lets follow/adapt `## version3 .json corrected and path. score corrected` from `004*`

BUT, lets ONLY do **simulation, not annotation**, as I cannot use *dbNSFP5.2a* which is designed for SNP.


--> Need to update:
- `scripts/simulate_array_v3.py` into `scripts/simulate_array_DBS.py`
- `scripts/simulate_mutations` into `scripts/simulate_mutations_DBS.py`



```bash
conda activate mutsim

# light test context
python scripts/simulate_array_DBS.py \
  --signature "AC>CA" \
  --n 500 \
  --rep 1 \
  --seed 1 \
  --sigfile signatures/context_sigs_fixed_DBS.txt \
  --contexts_list signatures/context_signature_list_DBS.txt \
  --context indices/context78.pkl \
  --exome parquet \
  --outdir results_contexts_DBS

#--> Works!

# light test DBS

python scripts/simulate_array_DBS.py \
  --signature "DBS1" \
  --n 500 \
  --rep 1 \
  --seed 1 \
  --sigfile signatures/COSMIC_v3.4_DBS_GRCh38_with_flat.txt \
  --contexts_list signatures/context_signature_list_DBS.txt \
  --context indices/context78.pkl \
  --exome parquet \
  --outdir results_DBS

#--> Works!




sbatch scripts/run_filtered_cosmic_DBS.slurm # 55679243 ok   --> results_DBS


sbatch scripts/run_filtered_contexts_DBS.slurm # 56796158 xxx   --> results_contexts_DBS



```






# Generate profile plots


Let's generate profile plots of my simulation to see if these are in agreement with the known profile.

- Custom script `scripts/plot_dbs78_from_parquet.py` adapted from `scripts/plot_sbs96_from_parquet.py` (`004*`) to make these plot from the `*.parquet`


```bash

plot_dbs78_from_parquet.py

# Check a few samples

python scripts/plot_dbs78_from_parquet.py \
  --parquet results_DBS/DBS19/n_4000/rep_05.sim.parquet \
  --fasta ref/GRCh38.primary_assembly.genome.fa \
  --sample-name DBS19-n_4000-rep_05 \
  --output-pdf plot/DBS19-n_4000-rep_05_plot.pdf

#--> Works!!!


```

--> It seems to be working!!!







# QC simulation



Check that the simulation is correct by investigating simulating parquet file


--> Generate `scripts/dbs_ref_alt.py` custom script that:
- Transform our parquet to indicate ref and alt BPs instead of using context IDs
- Generate a preview TSV (top 1000 rows)
- Generate a simple VCF for vizualization on IGV



```bash
scripts/dbs_ref_alt.py


python scripts/dbs_ref_alt.py \
  --parquet-in results_DBS/DBS19/n_4000/rep_01.sim.parquet \
  --fasta ref/GRCh38.primary_assembly.genome.fa \
  --parquet-out results_DBS/DBS19/n_4000/rep_01.with_refalt.parquet \
  --tsv-out results_DBS/DBS19/n_4000/rep_01.preview.tsv \
  --vcf-out results_DBS/DBS19/n_4000/rep_01.with_refalt.vcf









```










