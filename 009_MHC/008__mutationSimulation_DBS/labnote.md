# Project

Simulation DBS signature; Let's remodify `009*/005*`; but using the correct genome parquet file, re-generated in `009*/007*`




# Build the pkl 

--> copy all parquets files from `009*/007*`

Then build the *pkl* (label each bp that are in one of the 78 DBS mutation context) from the *parquet* (File with chr, pos, strand, ref_base, gene_id, codon_index, ref_codon, ref_aa for each bp of CDS):
- Parquet files already generated in `004` and has been copied in `parquet/`

--> Tricky part is to build the `CONTEXTS_78`; for that we need to find all 




```bash
conda activate mutsim
# To check file
parquet-tools show --head 5 parquet/chr1.parquet 
#--> GOOD, and 1-based

# Build pkl
python scripts/build_context78_index_v3.py # ok
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

print(f"\n📦 Total positions in parquet: {total_rows:,}") # 324,923,277
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

print(f"🧠 Total indices stored across all 78 contexts: {len(all_indices):,}") # 126,870,005
print(f"🔍 Unique positions covered (non-redundant rows): {len(unique_indices):,}") # 324,923,277
```

--> 324,923,277 positions (= all CDS bp) in both parquet and pkl index. 


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

Follow `009*/005*` but update `scripts/simulate_array_DBS.py` into `scripts/simulate_array_DBS_v2.py`


--> Need to update:
- `scripts/simulate_array_v3.py` into `scripts/simulate_array_DBS.py` then into `scripts/simulate_array_DBS_v2.py`
- `scripts/simulate_mutations` into `scripts/simulate_mutations_DBS.py` then into `scripts/simulate_mutations_DBS_v2.py`



```bash
conda activate mutsim


# light test DBS

python scripts/simulate_array_DBS_v2.py \
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


# light test context
python scripts/simulate_array_DBS_v2.py \
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




sbatch scripts/run_filtered_cosmic_DBS_v2.slurm # 62172664 xxx --> results_DBS
sbatch scripts/run_filtered_contexts_DBS_v2.slurm # 62172724 xxx --> results_contexts_DBS

# Check it is all good
parquet-tools show --head 5 results_DBS/DBS1/n_500/rep_01.sim.parquet
parquet-tools show --head 5 results_contexts_DBS/AC\>CA/n_500/rep_01.sim.parquet

#--> Looks good!




```


--> 






# Generate profile plots


Let's generate profile plots of my simulation to see if these are in agreement with the known profile.

- Custom script `scripts/plot_dbs78_from_parquet.py` adapted from `scripts/plot_sbs96_from_parquet.py` (`004*`) to make these plot from the `*.parquet`


```bash

plot_dbs78_from_parquet.py

# Check a few samples

python scripts/plot_dbs78_from_parquet.py \
  --parquet results_DBS/DBS1/n_500/rep_01.sim.parquet \
  --fasta ref/GRCh38.primary_assembly.genome.fa \
  --sample-name DBS1-n_500-rep_01 \
  --output-pdf plot/DBS1-n_500-rep_01-plot.pdf

#--> Works!!!


```

--> All good








# Add consequence

Here more complex than SBS, as one DBS mutation can alter two AA at the same time! Here is how we will proceed:
- If mutation affect 2 AA, keep only the most deleterious one
  - STOP > missense > synonymous (keep > leave)


## Proportion of simulated mutations affected 2 AA

First lets check the **proportion of each simulated mutations affected 2 AA**. Which is interesting as higher probabilities would indicate more probability to generate neoantigens!

The below script used the codon_index (0,1,2) column from each parquet file:
- value of 0= Bases 1 and 2 of the same codon are mutated (1 AA affected)
- value of 1= Bases 2 and 3 of the same codon are mutated (1 AA affected)
- value of 2= Bases 3 of this codond, and 1 of the next codon is affected (2 AA)

--> Count in each parquet `codon_index ==2`


--> Scripts updated from `009*/005*` to account for the 1-based


```bash
conda activate mutsim


parquet-tools show --head 5 results_DBS/DBS1/n_4000/rep_01.sim.parquet



# DBS ################
## Summarise/count double_AA mutation
python scripts/summarize_doubleAA_v2.py --root results_DBS

XXXY RUN AND MODIFY BELOW AFTER ALL IS RAN

## plot all double_AA mutation
python scripts/plot_doubleAA.py --root results_DBS --out results_DBS/DBS_doubleAA_prop.pdf
python scripts/plot_doubleAA-highlight_random.py --root results_DBS --out results_DBS/DBS_doubleAA_prop-highlight_random.pdf


# context ################
## Summarise/count double_AA mutation
python scripts/summarize_doubleAA.py --root results_contexts_DBS

## plot all double_AA mutation
python scripts/plot_doubleAA.py --root results_contexts_DBS --out results_contexts_DBS/contexts_DBS_doubleAA_prop.pdf
python scripts/plot_doubleAA-highlight_random.py --root results_contexts_DBS --flat results_DBS/Flat --out results_contexts_DBS/contexts_DBS_doubleAA_prop-highlight_random.pdf


```

--> Gloablly around 1/3, but still some differences between signatures!



Let's investigate whether result make sense, notably whether the context with higher probability to affect 2 AA are sequences that are often found overlapping two codons.


```bash
conda activate mutsim

python scripts/quantify_dbs_boundary_propensity.py
```

--> Output file: `dbs_ref2_boundary_propensity.tsv`; provide proportion of each context spanning two codons

- **AT-rich pairs** (AT, TT, AC, TC) are slightly *more likely to occur at codon boundaries* (ends of codons often A/T-rich).
- **GC-rich pairs** (CG, GC, CC) are slightly *less likely to sit at boundaries*.





## Add consequence (STOP, missense, synonymous)

Let's compute consequences at the codon level (not bp level); so from each simulated parquet we will:
- build the ALT codon(s) by applying the 2-bp change,
- compare the ALT to the REF AA
- for boundary cases (two codons hits by single DBS mutation), keep only the most deleterious one (keep the other as *secondary*.)
  - When codon_index ∈ {0,1} the 2 bases sit within one codon → we build that single ALT codon and compare AA.
  - When codon_index == 2, we touch two codons: pos2 of the first and pos0 of the next. We compute both AA outcomes and pick the most deleterious as primary (and keep the other as secondary for QC).


Updated `scripts/add_csq_v1.py` into `scripts/add_csq_v2.py`



```bash
conda activate mutsim


XXXY ALL GOOD RUN THIS WHEN simu generated

python scripts/add_csq_v2.py --root results_DBS --fasta ref/GRCh38.primary_assembly.genome.fa
python scripts/add_csq_v2.py --root results_contexts_DBS --fasta ref/GRCh38.primary_assembly.genome.fa





```

--> All good, new `consequence_summary_rep*.sim.json` file generated for each replicate




## Add path. score (SIFT4G, CADD, )


Let's use the same strategy to compile path. score

The script will:
- calls dbNSFP for each variant
- Maps scores/preds to the primary AA (worst of the two AAs when a boundary codon is hit).
  - If both SNVs are in the same codon (codon_index 0/1): aggregate to one codon score (SIFT4G=min, PP2=max, CADD=max).
  - If it’s a boundary (codon_index=2): use the SNV that belongs to the primary codon (the one with the worse AA consequence).
- Writes one JSON per replicate with counts, fractions and summary stats, and (optionally) a TSV of per-variant primary scores.

--> Work at the codon level: dbNSFP contain AA_ref and AA_alt so we have information at the AA level; simply collect the corresponding score generating the same AA change.





```bash
conda activate mutsim

# Prerequisete
## Build AA index version of dbNSFP
python scripts/build_dbnsfp_aa_index.py \
  --dbnsfp ref/dbNSFP5.3a_grch38.gz \
  --out ref/dbnsfp_aa_index.parquet


# Annotate each replicate by AA change
python scripts/annotate_dbs_scores_by_aa_v2.py \
  --parquet-in results_DBS/DBS1/n_500/rep_01.sim.parquet \
  --dbnsfp-aa-index ref/dbnsfp_aa_index.parquet \
  --out-prefix results_DBS/DBS1/n_500/pathogenicity_summary \
  --write-tsv


#--> Look good!


# Run all
XXXY double check that is goos!
sbatch scripts/run_annotate_dbs_scores_by_aa_v2.slurm





```







# Annotation summary plot


XXXY  do annotatio nplot for consequences and path., score








