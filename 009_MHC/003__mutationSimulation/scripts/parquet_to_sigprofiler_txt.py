import pandas as pd, pysam

def parquet_to_txt(parquet, fasta, sample, output):
    df = pd.read_parquet(parquet)
    fasta = pysam.FastaFile(fasta)

    df["genome_ref"] = df.apply(lambda r: fasta.fetch(r.chr, r.pos-1, r.pos).upper(), axis=1)
    df["simulated_ref"] = df["ref_base"].map({0:"A",1:"C",2:"G",3:"T"})
    df = df[df["genome_ref"] == df["simulated_ref"]]

    df_txt = pd.DataFrame({
        "Project": "Simu",
        "Sample": sample,
        "ID": ".",
        "Genome": "GRCh38",
        "mut_type": "SNP",
        "chrom": df["chr"].str.replace("chr", "", regex=True),
        "pos_start": df["pos"],
        "pos_end": df["pos"],
        "ref": df["genome_ref"],
        "alt": df["alt_base"],
        "Type": "SOMATIC"
    })
    df_txt.to_csv(output, sep="\t", index=False)
    print(f"✅ Saved {len(df_txt)} mutations to {output}")
