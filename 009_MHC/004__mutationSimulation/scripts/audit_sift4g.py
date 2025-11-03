#!/usr/bin/env python3
import sys, subprocess, pyarrow.parquet as pq
from collections import Counter, defaultdict

DB, PARQ = sys.argv[1], sys.argv[2]
df = pq.read_table(PARQ).to_pandas()
alphabet = "ACGT"

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

# header-driven dbNSFP columns
hdr = subprocess.run(["tabix","-H",DB], capture_output=True, text=True, check=True).stdout
cols = [l for l in hdr.splitlines() if l and not l.startswith("##")][-1].lstrip("#").split("\t")
colidx = {name:i for i,name in enumerate(cols)}
F_REF, F_ALT = 2, 3
F_TX      = colidx.get("Ensembl_transcriptid")
F_CANON   = colidx.get("VEP_canonical")
F_SIFT_S  = colidx.get("SIFT4G_score")
F_SIFT_P  = colidx.get("SIFT4G_pred")

def split_list(x):
    return None if x in (None,"",".") else x.split(";")

def any_numeric(lst):
    if not lst: return False
    for v in lst:
        try:
            if v not in (None,"",".") and float(v)==float(v):  # not NaN
                return True
        except: pass
    return False

cats = Counter()
examples = {}
by_conseq = defaultdict(Counter)

missing_mask = df["sift4g_score"].isna()
assert missing_mask.sum() + df["sift4g_score"].notna().sum() == len(df)

for i, row in df[missing_mask].iterrows():
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

    rows = [r for r in rows if len(r)>F_ALT and r[F_REF].upper()==ref and alt in r[F_ALT].upper().split(",")]
    if not rows:
        cats["alt_mismatch"] += 1
        by_conseq[conseq]["alt_mismatch"] += 1
        examples.setdefault("alt_mismatch", (chrom,pos,ref,alt))
        continue

    r = rows[0]
    sS = split_list(r[F_SIFT_S] if F_SIFT_S is not None and F_SIFT_S < len(r) else ".")
    cF = split_list(r[F_CANON]   if F_CANON   is not None and F_CANON   < len(r) else ".")
    # (rare) a row could lack the SIFT columns entirely in some builds
    if sS is None:
        cats["region_hit_but_no_sift_field"] += 1
        by_conseq[conseq]["region_hit_but_no_sift_field"] += 1
        examples.setdefault("region_hit_but_no_sift_field", (chrom,pos,ref,alt))
        continue

    any_val = any_numeric(sS)
    if not any_val:
        cats["region_hit_but_no_sift_value"] += 1  # SIFT present but all '.'
        by_conseq[conseq]["region_hit_but_no_sift_value"] += 1
        examples.setdefault("region_hit_but_no_sift_value", (chrom,pos,ref,alt))
        continue

    # if we get here, DB has some SIFT value(s) but your dataframe still has NA:
    # -> this means canonical was '.' and your annotator didn't fall back
    cats["canonical_dot_but_other_tx_has_value"] += 1
    by_conseq[conseq]["canonical_dot_but_other_tx_has_value"] += 1
    examples.setdefault("canonical_dot_but_other_tx_has_value", (chrom,pos,ref,alt))

print(f"Rows with NA sift4g_score in your parquet: {missing_mask.sum()}")
print("Breakdown:", dict(cats))
print("Examples:", examples)


