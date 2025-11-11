#!/usr/bin/env python3

import sys, subprocess, pyarrow.parquet as pq
from collections import Counter, defaultdict

DB, PARQ = sys.argv[1], sys.argv[2]
df = pq.read_table(PARQ).to_pandas()
alphabet = "ACGT"

# --- tabix helpers ---
def tabix(db, region):
    p = subprocess.run(["tabix", db, region], capture_output=True, text=True)
    if p.returncode != 0 or not p.stdout.strip():
        return []
    return [l.rstrip("\n").split("\t") for l in p.stdout.splitlines()]

def rows_at(db, chrom, pos):
    c = str(chrom)
    for cand in (c.lstrip("chr"), f"chr{c.lstrip('chr')}",
                 ("MT" if c in {"M","MT","chrM"} else None),
                 ("chrM" if c in {"M","MT","chrM"} else None)):
        if not cand: continue
        rows = tabix(db, f"{cand}:{pos}-{pos}")
        if rows: return rows
    return []

# --- header-driven dbNSFP columns ---
hdr = subprocess.run(["tabix","-H",DB], capture_output=True, text=True, check=True).stdout
cols = [l for l in hdr.splitlines() if l and not l.startswith("##")][-1].lstrip("#").split("\t")
colidx = {name:i for i,name in enumerate(cols)}

F_REF, F_ALT = 2, 3
F_TX     = colidx.get("Ensembl_transcriptid")
F_CANON  = colidx.get("VEP_canonical")

# PolyPhen2 HDIV fields
F_PP2_S  = colidx.get("Polyphen2_HDIV_score")
F_PP2_P  = colidx.get("Polyphen2_HDIV_pred")

def split_list(x):
    return None if x in (None,"",".") else x.split(";")

def any_numeric(lst):
    if not lst: return False
    for v in lst:
        try:
            if v not in (None,"",".") and float(v)==float(v):
                return True
        except:
            pass
    return False

cats = Counter()
examples = {}
by_conseq = defaultdict(Counter)

# align with your plot: investigate rows where the numeric PP2 score is NA
missing_mask = df["polyphen2_hdiv_score"].isna()

for _, row in df[missing_mask].iterrows():
    chrom, pos = row["chr"], int(row["pos"])
    ref = alphabet[int(row["ref_base"])].upper()
    alt = str(row["alt_base"]).upper()
    conseq = str(row.get("consequence","unknown")).lower()

    rows = rows_at(DB, chrom, pos)
    if not rows:
        cats["no_region_hit"] += 1
        by_conseq[conseq]["no_region_hit"] += 1
        examples.setdefault("no_region_hit", (chrom,pos,ref,alt))
        continue

    # require exact REF and ALT match (ALT may be comma-separated)
    rows = [r for r in rows if len(r)>F_ALT and r[F_REF].upper()==ref and alt in r[F_ALT].upper().split(",")]
    if not rows:
        cats["alt_mismatch"] += 1
        by_conseq[conseq]["alt_mismatch"] += 1
        examples.setdefault("alt_mismatch", (chrom,pos,ref,alt))
        continue

    r = rows[0]
    pp2S = split_list(r[F_PP2_S] if (F_PP2_S is not None and F_PP2_S < len(r)) else ".")
    pp2P = split_list(r[F_PP2_P] if (F_PP2_P is not None and F_PP2_P < len(r)) else ".")
    canon = split_list(r[F_CANON]  if (F_CANON  is not None and F_CANON  < len(r)) else ".")

    # (rare) line lacks PP2 columns entirely
    if F_PP2_S is None or (pp2S is None and F_PP2_S < len(r) and r[F_PP2_S] in (".","")):
        cats["region_hit_but_no_pp2_field"] += 1
        by_conseq[conseq]["region_hit_but_no_pp2_field"] += 1
        examples.setdefault("region_hit_but_no_pp2_field", (chrom,pos,ref,alt))
        continue

    # columns exist; check if all values missing
    if not any_numeric(pp2S):
        cats["region_hit_but_no_pp2_value"] += 1
        by_conseq[conseq]["region_hit_but_no_pp2_value"] += 1
        examples.setdefault("region_hit_but_no_pp2_value", (chrom,pos,ref,alt))
        continue

    # if we reach here, dbNSFP has some PP2 number but your parquet shows NA:
    # -> canonical likely '.' while other transcripts have values
    cats["canonical_dot_but_other_tx_has_value"] += 1
    by_conseq[conseq]["canonical_dot_but_other_tx_has_value"] += 1
    examples.setdefault("canonical_dot_but_other_tx_has_value", (chrom,pos,ref,alt))

print(f"Rows with NA polyphen2_hdiv_score in your parquet: {int(missing_mask.sum())}")
print("Breakdown:", dict(cats))
print("\nBy consequence (top buckets):")
for k,v in sorted(by_conseq.items(), key=lambda kv: -sum(kv[1].values()))[:10]:
    print(f"  {k}: {dict(v)} (n={sum(v.values())})")
print("\nExamples:", examples)
