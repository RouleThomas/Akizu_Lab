# audit_dbnsfp_hits.py
import subprocess, sys, pyarrow.parquet as pq
from collections import Counter

db = sys.argv[1]   # e.g. ref/dbNSFP5.2a_grch38.gz
pq_in = sys.argv[2]

df = pq.read_table(pq_in).to_pandas()
alphabet = "ACGT"

def tabix_rows(db, chrom, pos):
    out = subprocess.run(["tabix", db, f"{chrom}:{pos}-{pos}"],
                         capture_output=True, text=True)
    if out.returncode != 0 or not out.stdout.strip():
        return []
    return [l.rstrip("\n").split("\t") for l in out.stdout.splitlines()]

def has_ref(db_rows, ref):
    return any(len(r)>=4 and r[2].upper()==ref for r in db_rows)

def has_alt(db_rows, ref, alt):
    return any(len(r)>=4 and r[2].upper()==ref and alt in r[3].upper().split(",") for r in db_rows)

cats = Counter()
examples = {}

for chrom, pos, r_idx, a in zip(df.chr, df.pos, df.ref_base, df.alt_base):
    ref = alphabet[int(r_idx)].upper()
    alt = str(a).upper()
    chrom_str = str(chrom)

    # try without and with "chr"
    rows = tabix_rows(db, chrom_str.lstrip("chr"), int(pos))
    if not rows:
        rows = tabix_rows(db, f"chr{chrom_str.lstrip('chr')}", int(pos))
    if not rows:
        cats["no_region_hit"] += 1
        examples.setdefault("no_region_hit", (chrom, pos, ref, alt))
        continue

    if not has_ref(rows, ref):
        # probe ±1 to flag off-by-one
        m1 = tabix_rows(db, chrom_str.lstrip("chr"), int(pos)-1)
        p1 = tabix_rows(db, chrom_str.lstrip("chr"), int(pos)+1)
        if m1 or p1:
            cats["pos_off_by_one"] += 1
            examples.setdefault("pos_off_by_one", (chrom, pos, ref, alt))
        else:
            cats["ref_mismatch"] += 1
            examples.setdefault("ref_mismatch", (chrom, pos, ref, alt))
        continue

    if not has_alt(rows, ref, alt):
        cats["alt_mismatch"] += 1
        examples.setdefault("alt_mismatch", (chrom, pos, ref, alt))
        continue

    cats["perfect_match"] += 1

print("Hit categories:", dict(cats))
print("Examples:", examples)