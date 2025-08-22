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
## version1 .json issue with STOP

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
sbatch scripts/run_filtered_experimental.slurm # 47875937 FAIL not enough array; 49631598 ok --> results_experimental/
sbatch scripts/run_filtered_contexts.slurm # 47876020 FAIL not enough array; TO LAUNCH  --> results_contexts/
```
--> All good all files generated with `n_*` folders in the respectrive `results*/` folders

--> *NOTE: I did an error but not using enough array for experimental and context... So ran `*missing*` jobs for the missing ones.*

--> ERROR with the .json file I count the stop->stop... Start with context signature and verify that context that should not give STOP does not give STOP (`A[T>C]A`)!



# Simulate mutation 
## version2 .json corrected





```bash
conda activate mutsim

# Light testing with a context that should not have STOP `A[T>C]A`
python scripts/simulate_array_v2.py \
  --signature "A[T>C]A" \
  --n "4000" \
  --rep "1" \
  --seed "42" \
  --sigfile signatures/context_sigs_fixed.txt \
  --outdir results_contexts 
#--> STOP gained at 0 (GOOD!)

# Light testing with a context that should not have STOP `T[C>A]A`
python scripts/simulate_array_v2.py \
  --signature "T[C>A]A" \
  --n "4000" \
  --rep "1" \
  --seed "42" \
  --sigfile signatures/context_sigs_fixed.txt \
  --outdir results_contexts 
#--> STOP gained at 301 (GOOD!)





# Generate plot for all (script updated to use `simulate_array_v2.py`)
sbatch scripts/run_filtered_contexts_v2.slurm # 50364278 xxx --> results_contexts/
sbatch scripts/run_filtered_cosmic_v2.slurm # 50364847 xxx --> results/
sbatch scripts/run_filtered_experimental_v2.slurm # 50365053 xxx --> results_experimental/
```

XXX HERE!!!
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

### Version1 FAIL

Let's generate plot with `x` = `n_mutations` and `y` = `annotation score`



```bash
conda activate mutsim

####################################
# COSMIC ###########################
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
# EXPERIMENTAL #####################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-experimental.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-experimental.py
#--> results_experimental/combined_signature_summary_errorbars.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-experimental-highlight_random.py
#--> results_experimental/combined_signature_summary_plots-highlight_random.pdf

####################################
# CONTEXTS ##########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-contexts.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-contexts.py
#--> results_contexts/combined_signature_summary_errorbars.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-contexts-highlight_random.py
#--> results_contexts/combined_signature_summary_plots-highlight_random.pdf
```

--> All good, all summary annottion files and plot generated
  --> ERROR! Realized there was an error with the stop count (I count all stops; even the stop->stop). BUT, not because of these scripts, but because of the `.json`...




### Version2 - FAIL STOP, unnecessary correction

--> Same name with `*_v2*` suffix. (realized STOP issue related to `.json`, not this sumary plot part)



```bash
conda activate mutsim


####################################
# CONTEXTS ##########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot_v2-contexts.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined_v2-contexts.py
#--> results_contexts/combined_signature_summary_errorbars.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined_v2-contexts-highlight_random.py
#--> results_contexts/combined_signature_summary_plots-highlight_random.pdf



```

--> FAIL, the issue comes from the `.json`, not these scripts








# QC simulation

## Check issue with mutation center codon

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




## Check issue with stop codon - context

Let's check issue with stop codon related to context; annotation plot show none of my contexts showing 0 for stop codon, but some context should NEVER produce a stop codon, thus, be at 0 all the time, whatever the n_mutations...


STOP codon= UAA, UAG, and UGA = TAA, TAG, TGA

Correspond to contexts: 'T[C>A]A','T[C>A]G','T[C>G]A','T[T>A]A','T[T>A]G','T[T>G]A'

But then we do not know WHERE in the codon the bp is mutated! According to Joan, some context should still NEVER produce a stop codons. But I think it is wrong as mutation can be apply to first, second, or last codon, so CDS codon need to be taken into account.

```ruby
For the stop codon; maybe I 'm wrong, but I am not sure that we expect some context to always give 0 STOP. This would only happen if the mutated base was always the middle base of the codon, but here I didn’t force the mutated base to be codon position 2 (middle) in the CDS/genome.
Not sure I am clear :sweat_smile: but as a counter example: T[C>A]A lead to TAA which is a stop codon. Thus, it will only lead to a stop codon if the C>A is located at position2 in the human CDS real codon. For example if sequence is AAG.TCA.TTA --> mutation will lead to AAG.TAA.TTA=AAG.STOP.TTA (middle base pair mutated). But if codons are organized like this (3rd base pair mutated):  AA.GTC.ATT.A, then mutation is AA.GTA.ATT.A : no STOP codon in this case. In other words depending on the position of the mutated base pair any context can produce a stop codon?
```

--> No my point is not true, for example: A[T>C]A T>C will NEVER give stop codon, as stop codon do NOT have C... Let's make a plot to highlight A[T>C]A annotation; like instead of Flat, highlight this one



```bash
conda activate mutsim


####################################
# CONTEXTS ##########################
####################################

python scripts/plot_all_signatures_combined-contexts-highlight_ATCA.py
#--> results_contexts/combined_signature_summary_plots-highlight_random.pdf
```

--> A[T>C]A got a low fraction STOP, but not a value of 0... Maybe because it is a mutation that change a stop codon to a stop codon? Or other issue direclty related to how my plot is make (I think related to forward/reverse...). So probably issue with how frac_stop is computed! As should be 0 for this one!



Issue from `## Annotation summary plot` scripts; when computing the stop; I include stop_retained (stop→stop) and/or stop_lost, so even a T>C context like `A[T>C]A` ends up with non-zero "stop"...
  --> These scripts are wrong: `scripts/run_summary_plot-*.sh`
    --> NO!! These scripts were good, `.json` is not good!



