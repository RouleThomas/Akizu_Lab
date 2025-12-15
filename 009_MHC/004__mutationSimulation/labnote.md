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
sbatch scripts/run_filtered_contexts_v2.slurm # 50364278 ok --> results_contexts/
sbatch scripts/run_filtered_cosmic_v2.slurm # 50364847 ok --> results/
sbatch scripts/run_filtered_experimental_v2.slurm # 50365053 ok --> results_experimental/
```

--> All good all files generated with `n_*` folders in the respectrive `results*/` folders
  --> No file missing!




# Simulate mutation  - SBS, Exp, context
## version3 .json corrected and path. score corrected

I noticed that CADD Phred score used the raw version; not the normalized one that goes from 0-99... The wrong column might be picked in the `scripts/annotate_damage3.py`.

Below code provide cell name and ID of all dbNSFP5 db:

```python
import gzip

with gzip.open("ref/dbNSFP5.2a_grch38.gz", "rt") as f:
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
  -  IDX_CADD_PHRED       = 145 - 1 --> NOT OK!!! Correct is 144!!!


--> `scripts/annotate_damage4.py` now updated with the correct CADD_PHRED!



```bash
conda activate mutsim


# Generate plot for all (script updated to use `scripts/simulate_array_v3.py` and `scripts/annotate_damage4.py`)

sbatch scripts/run_filtered_cosmic_v3.slurm # 55881484 ok  --> results_v3/
sbatch scripts/run_filtered_contexts_v3.slurm # 55881483 ok  --> results_contexts_v3/
sbatch scripts/run_filtered_experimental_v3.slurm # 56796128 ok  --> results_experimental_v3/
```








## version4 update for missing dbNSFP5

I noticed that many dbNSFP5 scores were not retrieved, it seems to be because dbNSFP/tabix is 1-based and my simulated parquet at 0-based. So let's make them 1-based like dbNSFP5.

--> Modify `scripts/simulate_array_v3.py` into `scripts/simulate_array_v4.py`





```bash
conda activate mutsim

# Test with one simulation

python scripts/simulate_array_v4.py \
  --signature SBS90 \
  --n 4000 \
  --rep 1 \
  --seed 42 \
  --sigfile signatures/COSMIC_with_flat.txt \
  --outdir results_v4 # UPDATE PATH






# Generate plot for all (script updated to use `scripts/simulate_array_v3.py` and `scripts/annotate_damage4.py`)

sbatch scripts/run_filtered_cosmic_v4.slurm #  xxx --> results_v4/

XXX TO MODIFY 


sbatch scripts/run_filtered_contexts_v4.slurm #  xxx --> results_contexts_v4/
sbatch scripts/run_filtered_experimental_v4.slurm #  xxx --> results_experimental_v4/
```


--> Well still a lot of missing! Could be gencode version; I simualted on gencode v45 but dbNSFP5.2a use v48...






# Generate profile plots


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

python scripts/plot_sbs96_from_parquet.py \
  --parquet results/SBS4/n_4000/rep_01.annot.parquet \
  --fasta ref/GRCh38.primary_assembly.genome.fa \
  --sample-name SBS4-n_4000-rep_01 \
  --output-pdf plot/SBS4-n_4000-rep_01_plot.pdf

```

--> It seems to be working!!!



# Annotation summary plot

## Version1 FAIL

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

    --> Scripts re-launch in `### Version3`



## Version2 - FAIL STOP, unnecessary correction

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





## Version3

Let's generate plot with `x` = `n_mutations` and `y` = `annotation score`

--> Added `*_v3` at `*_summary_all_v3.tsv` and `*_summary_plots_v3.pdf`
--> Updated `scripts/summarize_simulation_results.py` into `scripts/summarize_simulation_results_v3.py`


```bash
conda activate mutsim

####################################
# CONTEXTS ##########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-contexts_v3.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-contexts_v3.py
#--> results_contexts/combined_signature_summary_errorbars.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-contexts-highlight_random_v3.py
#--> results_contexts/combined_signature_summary_plots-highlight_random.pdf

####################################
# COSMIC ###########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-cosmic_v3.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-cosmic_v3.py
#--> results/combined_signature_summary_errorbars.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-cosmic-highlight_random_v3.py
#--> results/combined_signature_summary_plots-highlight_random.pdf

####################################
# EXPERIMENTAL #####################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-experimental_v3.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-experimental_v3.py
#--> results_experimental/combined_signature_summary_errorbars.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-experimental-highlight_random_v3.py
#--> results_experimental/combined_signature_summary_plots-highlight_random.pdf

```

--> All good, all summary annottion files and plot generated
  --> ERROR! Realized there was an error with the stop count (I count all stops; even the stop->stop). BUT, not because of these scripts, but because of the `.json`...






## Version4

From the `## version3 .json corrected and path. score corrected` simulation where CADD_PHRED score corrected.. 

--> Annotate using consequences of the pathogeneicity, rather than score; so for **pathogeneicity use qualitative and not continuous variables.**
  --> Still keep same plot, just add horizontal lines for path. threshold:
    - SIFT4G: < 0.05 = damaging; > 0.05 = benign
    - PolyPhen2: 0-0.15= benign; 0.15-0.85 = possibly damaging; >0.85 = damaging
    - CADD: < 2= benign; >10 = likely damaging, >20 = damaging


--> **Exact same code as v3; I just now use the results_v3 folders!**
  --> And output `*_summary_all_v4.tsv` and `*_summary_plots_v4.pdf`


```bash
conda activate mutsim

####################################
# CONTEXTS ##########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-contexts_v4.sh

# One plot with all mutations
python scripts/plot_all_signatures_combined-contexts_v4.py
#--> results_contexts_v3/combined_signature_summary_errorbars_v4.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-contexts-highlight_random_v4.py
#--> results_contexts_v3/combined_signature_summary_plots-highlight_random_v4.pdf

####################################
# COSMIC ###########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-cosmic_v4.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-cosmic_v4.py
#--> results_v3/combined_signature_summary_errorbars_v4.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-cosmic-highlight_random_v4.py
#--> results_v3/combined_signature_summary_plots-highlight_random_v4.pdf



####################################
# EXPERIMENTAL #####################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-experimental_v4.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-experimental_v4.py
#--> results_experimental_v3/combined_signature_summary_errorbars_v4.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-experimental-highlight_random_v4.py
#--> results_experimental_v3/combined_signature_summary_plots-highlight_random_v4.pdf
```



For **SIFT4G, Polyphen2 and CADD**; lets instead use **percent of damaging/benign for the y axis** (ie. and not the score values) --> Will lead to 2 plot for SIFT4G (ie. Damaging, benign) and 3 plot for Polyphen2 (ie. Damaging, Possible damaging, benign)


- update `scripts/summarize_simulation_results_v3.py` into `scripts/summarize_simulation_results_v3_dbNSFP5.py`
- update `scripts/plot_signature_summary.py` into `scripts/plot_signature_summary_dbNSFP5.py`


```bash
conda activate mutsim

####################################
# CONTEXTS ##########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-contexts_v4_dbNSFP5.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-contexts_v4_dbNSFP5.py
#--> results_contexts_v3/combined_signature_summary_errorbars_v4_dbNSFP5.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-contexts-highlight_random_v4_dbNSFP5.py
#--> results_contexts_v3/combined_signature_summary_plots-highlight_random_v4_dbNSFP5.pdf


####################################
# COSMIC ###########################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-cosmic_v4_dbNSFP5.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-cosmic_v4_dbNSFP5.py
#--> results_v3/combined_signature_summary_errorbars_v4_dbNSFP5.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-cosmic-highlight_random_v4_dbNSFP5.py
#--> results_v3/combined_signature_summary_plots-highlight_random_v4_dbNSFP5.pdf



####################################
# EXPERIMENTAL #####################
####################################

# Plot and summary metric for each mutations
bash scripts/run_summary_plot-experimental_v4_dbNSFP5.sh
# One plot with all mutations
python scripts/plot_all_signatures_combined-experimental_v4_dbNSFP5.py
#--> results_experimental_v3/combined_signature_summary_errorbars_v4_dbNSFP5.pdf
# One plot with all mutations - Random/Flat highlighted
python scripts/plot_all_signatures_combined-experimental-highlight_random_v4_dbNSFP5.py
#--> results_experimental_v3/combined_signature_summary_plots-highlight_random_v4_dbNSFP5.pdf
```

--> All good



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




## Check why so many missing in dbNSFP5 score

Version `_v3` gave a lot of Missing mutation: mutation without a score; notably for SIFT4G, HDIV (and CADD, but far less)



```bash
conda activate mutsim

python scripts/audit_dbnsfp_hits.py ref/dbNSFP5.2a_grch38.gz results_v3/SBS90/n_4000/rep_01.annot.parquet
#--> No region hit: chr2:74523421 C > A 



# CHECK SIFT4G
python scripts/audit_sift4g.py ref/dbNSFP5.2a_grch38.gz results_v3/SBS90/n_4000/rep_01.annot.parquet

# CHECK PolyPhen2 HDIV
python scripts/audit_polyphen2_hdiv.py ref/dbNSFP5.2a_grch38.gz results_v3/SBS90/n_4000/rep_01.annot.parquet

# CHECK CADD PHRED
python scripts/audit_cadd_phred.py ref/dbNSFP5.2a_grch38.gz results_v3/SBS90/n_4000/rep_01.annot.parquet


# CHECK CADD PHRED - version4 0-based 1-based
python scripts/audit_cadd_phred.py ref/dbNSFP5.2a_grch38.gz results_v4/SBS90/n_4000/rep_01.annot.parquet

# Double check by looking at the dbNSFP5
tabix ref/dbNSFP5.2a_grch38.gz 2:74523421-74523421 #--> nothing!
```

- SIFT4G, and Polyph2 are for non-synonymous only (ie. no score for synonymous mutation)
- CADD Phred score for any kind of mutations.

--> **No score are missing for a technical reason**; score are missing, either because:
- they are not present in the dbNSFP5 file
- mutation is present but not score provided for the canonical transcript (it. but alternative transcript got a score)
- mutation is present but no score for either SIFT4G or Polyphenm (ie. but could be a score for CADD PHRED...)


--> **Testing v3 vs v4 showed that in v4, show that in v3 lot of no_region_hit; which is bad and show that position is not correct; however in v4 the missing are mostly ALT mismatch; meaning it is likely a synonymous mutation which dbNSFP5 sometime does not include**


Check more carefully why some mutations are NOT included in the dbNSFP5

```python
import pyarrow.parquet as pq
import pandas as pd

df = pq.read_table(
    "results_v3/SBS90/n_4000/rep_01.annot.parquet"
).to_pandas()


db = "ref/dbNSFP5.2a_grch38.gz"

miss = df[df["cadd_phred"].isna()].sample(20, random_state=1)

import subprocess

def tabix_has_hit(db, chrom, pos):
    reg = f"{chrom}:{pos}-{pos}"
    p = subprocess.run(["tabix", db, reg], capture_output=True, text=True)
    return bool(p.stdout.strip())

hits_pos   = sum(tabix_has_hit(db, str(r["chr"]).lstrip("chr"), int(r["pos"]))     for _,r in miss.iterrows())
hits_posp1 = sum(tabix_has_hit(db, str(r["chr"]).lstrip("chr"), int(r["pos"]) + 1) for _,r in miss.iterrows())

print("hits at pos:", hits_pos)
print("hits at pos+1:", hits_posp1)

```

--> Show that the `simulation/parquet` pos is 0-based and `dbNSFP/tabix` is 1-based --> Lets generate a new 1-based mutation as *version4*



```python
import pyarrow.parquet as pq
import pandas as pd

df = pq.read_table(
    "results_v4/SBS90/n_4000/rep_01.annot.parquet"
).to_pandas()


db = "ref/dbNSFP5.2a_grch38.gz"

miss = df[df["cadd_phred"].isna()].sample(20, random_state=1)

import subprocess

def tabix_has_hit(db, chrom, pos):
    reg = f"{chrom}:{pos}-{pos}"
    p = subprocess.run(["tabix", db, reg], capture_output=True, text=True)
    return bool(p.stdout.strip())

hits_pos   = sum(tabix_has_hit(db, str(r["chr"]).lstrip("chr"), int(r["pos"]))     for _,r in miss.iterrows())
hits_posp1 = sum(tabix_has_hit(db, str(r["chr"]).lstrip("chr"), int(r["pos"]) + 1) for _,r in miss.iterrows())

print("hits at pos:", hits_pos)
print("hits at pos+1:", hits_posp1)

```

--> Seems to contain as many missing...


So it seems that the 0 and 1-based pos is important; and also the GENCODE version may be important too... Let's try a new version 007* where I downlaod the last version of dbNSFP together with the corresponding gencode version... This should work optimally!





# Score distribution  version1 continuous values

Instead of relying on the error-bar average for each annotation score; lets check specifically their dispersion: is there some extremes or most mutations fall within the error bar?

--> we displayed the average, and the average score is low for most, but how is the distribution? Maybe there is many mutation with a high score that are mask with the average

## SBS

Let's now check for each SBS signature the profile/distribution of each scores for 4k mutations (ie. 4k is where variance stabilizes)

For each SBS signature:
- Loop through `results/*/*/` (All folders: all SBS and Flat)
- Generates one 6-panel page per replicate (6 scores annotations)

--> `*scoreUpdate` is by translating path. score to consequences:
  - SIFT4G: < 0.05 = damaging; > 0.05 = benign
  - PolyPhen2: 0-0.15= benign; 0.15-0.85 = possibly damaging; >0.85 = damaging
  - CADD: < 2= benign; >10 = likely damaging, >20 = damaging

```bash
conda activate mutsim

# light testing on one signature
python scripts/plot_distributions_per_replicate-SBS2_4k.py
#--> Works!
python scripts/plot_distributions_per_replicate-SBS2_4k-scoreUpdate.py




# Run all SBS signatures and Flat
python scripts/plot_distributions_per_replicate-SBS.py


# Run 

```
--> Plot produced at `results/*/*/*_replicate_score_distributions.pdf`


Some mutations with extreme score are masked with average.
  --> Average score can be biased and contain mutations with extreme csq



## Experimental


Same for Experimental signatures



```bash
conda activate mutsim

# Run all Experimental signatures
python scripts/plot_distributions_per_replicate-experimental.py


```
--> Plot produced at `results_experimental/*/*/*_replicate_score_distributions.pdf`




XXXY RUN THE TWO BELOW, code ready!



## Context


Same for Experimental signatures


```bash
conda activate mutsim

# Run all contexts signatures
python scripts/plot_distributions_per_replicate-contexts.py


```
--> Plot produced at `results_contexts/*/*/*_replicate_score_distributions.pdf`










# Score distribution  version2 qualitative values

Same as version1, but this time display count as begnin, damaging for the path. scores
  - SIFT4G: < 0.05 = damaging; > 0.05 = benign
  - PolyPhen2: 0-0.15= benign; 0.15-0.85 = possibly damaging; >0.85 = damaging
  - CADD: < 2= benign; >10 = likely damaging, >20 = damaging



## SBS

Let's now check for each SBS signature the profile/distribution of each scores for 4k mutations (ie. 4k is where variance stabilizes)

For each SBS signature:
- Loop through `results/*/*/` (All folders: all SBS and Flat)
- Generates one 6-panel page per replicate (6 scores annotations)

--> `*scoreUpdate` is by translating path. score to consequences:
  - SIFT4G: < 0.05 = damaging; > 0.05 = benign
  - PolyPhen2: 0-0.15= benign; 0.15-0.85 = possibly damaging; >0.85 = damaging
  - CADD: < 2= benign; >10 = likely damaging, >20 = damaging



```bash
conda activate mutsim

# light testing on one signature
python scripts/plot_distributions_per_replicate_version2-SBS2_4k.py
#--> Works!




# Run all SBS signatures and Flat
python scripts/plot_distributions_per_replicate_version2-SBS.py




```
--> Plot produced at `results_v3/*/*/*_replicate_score_distributions.pdf`


Some mutations with extreme score are masked with average.
  --> Average score can be biased and contain mutations with extreme csq




## Context


Same for Experimental signatures


```bash
conda activate mutsim

# Run all contexts signatures
python scripts/plot_distributions_per_replicate_version2-contexts.py


```
--> Plot produced at `results_contexts/*/*/*_replicate_score_distributions.pdf`




## Experimental


XXXY HERE BELOW TO RUN!!!! AFTER  sbatch scripts/run_filtered_experimental_v3.slurm # 55881533 xxx --> results_experimental_v3/

Same for Experimental signatures



```bash
conda activate mutsim

# Run all Experimental signatures
python scripts/plot_distributions_per_replicate_version2-experimental.py


```
--> Plot produced at `results_experimental/*/*/*_replicate_score_distributions.pdf`



