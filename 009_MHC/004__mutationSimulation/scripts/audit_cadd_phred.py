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
F_CADD       = colidx.get("CADD_phred")  # dbNSFP 5.2a GRCh38

cats = Counter()
examples = {}
by_conseq = defaultdict(Counter)

missing_mask = df["cadd_phred"].isna()

def is_num(x):
    try:
        if x in (None, "", "."): return False
        float(x)
        return True
    except:
        return False

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

    # exact REF and ALT (ALT may be comma-separated)
    rows = [r for r in rows if len(r)>F_ALT and r[F_REF].upper()==ref and alt in r[F_ALT].upper().split(",")]
    if not rows:
        cats["alt_mismatch"] += 1
        by_conseq[conseq]["alt_mismatch"] += 1
        examples.setdefault("alt_mismatch", (chrom,pos,ref,alt))
        continue

    r = rows[0]

    # CADD field presence
    if F_CADD is None or F_CADD >= len(r):
        cats["region_hit_but_no_cadd_field"] += 1
        by_conseq[conseq]["region_hit_but_no_cadd_field"] += 1
        examples.setdefault("region_hit_but_no_cadd_field", (chrom,pos,ref,alt))
        continue

    cadd_val = r[F_CADD]

    # CADD value content
    if cadd_val in (None, "", "."):
        cats["region_hit_but_no_cadd_value"] += 1
        by_conseq[conseq]["region_hit_but_no_cadd_value"] += 1
        examples.setdefault("region_hit_but_no_cadd_value", (chrom,pos,ref,alt))
        continue

    if not is_num(cadd_val):
        cats["cadd_non_numeric"] += 1
        by_conseq[conseq]["cadd_non_numeric"] += 1
        examples.setdefault("cadd_non_numeric", (chrom,pos,ref,alt))
        continue

    # If we got here, dbNSFP *has* a numeric CADD but your parquet shows NA.
    # This flags a parsing/assignment issue in the annotator (should be rare).
    cats["has_cadd_value_but_parquet_na"] += 1
    by_conseq[conseq]["has_cadd_value_but_parquet_na"] += 1
    examples.setdefault("has_cadd_value_but_parquet_na", (chrom,pos,ref,alt))

print(f"Rows with NA cadd_phred in your parquet: {int(missing_mask.sum())}")
print("Breakdown:", dict(cats))
print("\nBy consequence (top buckets):")
for k,v in sorted(by_conseq.items(), key=lambda kv: -sum(kv[1].values()))[:10]:
    print(f"  {k}: {dict(v)} (n={sum(v.values())})")
print("\nExamples:", examples)
