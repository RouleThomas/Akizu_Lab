#!/bin/bash
#SBATCH --job-name=dbs_aa_scores
#SBATCH --output=logs/dbs_scores_%A_%a.out
#SBATCH --error=logs/dbs_scores_%A_%a.err
#SBATCH --array=0-1469%50
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=defq

set -euo pipefail

# ---- env ----
source ~/.bashrc
conda activate mutsim

mkdir -p logs

# ---- params (MUST match your simulation layout) ----
SIGNATURE_LIST=( $(<signatures/DBS_COSMIC_signatures.txt) )   # e.g. DBS1 DBS2 ... and (optionally) Flat
N_LIST=(500 1000 2000 4000 8000 16000 32000)
REPS=10

TOTAL_SIGNATURES=${#SIGNATURE_LIST[@]}
TOTAL_NS=${#N_LIST[@]}
TOTAL_JOBS=$((TOTAL_SIGNATURES * TOTAL_NS * REPS))

# ---- decode array index -> (signature, n, rep) ----
IDX=${SLURM_ARRAY_TASK_ID}

SIGNATURE_IDX=$((IDX / (REPS * TOTAL_NS)))
N_IDX=$(((IDX / REPS) % TOTAL_NS))
REP=$((IDX % REPS + 1))

SIGNATURE=${SIGNATURE_LIST[$SIGNATURE_IDX]}
MUTATION_COUNT=${N_LIST[$N_IDX]}

echo "Job IDX: $IDX / $((TOTAL_JOBS-1))"
echo "Signature: $SIGNATURE"
echo "n: $MUTATION_COUNT"
echo "rep: $REP"

# ---- paths ----
IN_PARQ="results_DBS/${SIGNATURE}/n_${MUTATION_COUNT}/rep_$(printf "%02d" ${REP}).sim.parquet"
OUT_PARQ="results_DBS/${SIGNATURE}/n_${MUTATION_COUNT}/rep_$(printf "%02d" ${REP}).sim.with_aapos.parquet"

GTF="gtf/gencode.v45.annotation.gtf.gz"
CACHE="meta/cds_models_cache.pkl"
DBNSFP_AA="ref/dbnsfp_aa_index.parquet"
OUT_PREFIX="results_DBS/${SIGNATURE}/n_${MUTATION_COUNT}/pathogenicity_summary"

echo "IN_PARQ:  $IN_PARQ"
echo "OUT_PARQ: $OUT_PARQ"

# ---- sanity ----
if [[ ! -s "$IN_PARQ" ]]; then
  echo "SKIP (missing input): $IN_PARQ"
  exit 0
fi

# ---- step 1: add aapos (skip if exists) ----
if [[ ! -s "$OUT_PARQ" ]]; then
  python scripts/add_aapos_from_gtf.py \
    --parquet-in "$IN_PARQ" \
    --gtf "$GTF" \
    --parquet-out "$OUT_PARQ" \
    --cache-pkl "$CACHE"
else
  echo "SKIP add_aapos (exists)"
fi

# ---- step 2: annotate scores by AA (skip if JSON exists for this rep) ----
JSON_OUT="results_DBS/${SIGNATURE}/n_${MUTATION_COUNT}/pathogenicity_summary_rep_$(printf "%02d" ${REP}).sim.json"

if [[ ! -s "$JSON_OUT" ]]; then
  python scripts/annotate_dbs_scores_by_aa.py \
    --parquet-in "$OUT_PARQ" \
    --dbnsfp-aa-index "$DBNSFP_AA" \
    --out-prefix "$OUT_PREFIX" \
    --write-tsv
else
  echo "SKIP annotate (exists): $JSON_OUT"
fi

echo "DONE"
