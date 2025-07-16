import pandas as pd
import pysam

def parquet_to_txt(parquet, fasta, sample_name, output):
    df = pd.read_parquet(parquet)
    fasta = pysam.FastaFile(fasta)
    
    def get_genomic_ref(row):
        try:
            return fasta.fetch(row["chr"], row["pos"] - 1, row["pos"]).upper()
        except:
            return "N"

    df["genome_ref"] = df.apply(get_genomic_ref, axis=1)
    df["simulated_ref"] = df["ref_base"].map({0:"A", 1:"C", 2:"G", 3:"T"})

    mismatches = df[df["simulated_ref"] != df["genome_ref"]]
    matches = df[df["simulated_ref"] == df["genome_ref"]]

    print(f"\n✅ {len(matches)} / {len(df)} matches")
    print(f"❌ {len(mismatches)} mismatches\n")
    print(mismatches[["chr", "pos", "simulated_ref", "genome_ref", "alt_base", "revcomp"]].head())

    df_txt = matches.copy()
    df_txt["chr_clean"] = df_txt["chr"].str.replace("chr", "", regex=True)
    df_txt_final = pd.DataFrame({
        "Project": "Simu",
        "Sample": sample_name,
        "ID": ".",
        "Genome": "GRCh38",
        "mut_type": "SNP",
        "chrom": df_txt["chr_clean"],
        "pos_start": df_txt["pos"],
        "pos_end": df_txt["pos"],
        "ref": df_txt["genome_ref"],
        "alt": df_txt["alt_base"],
        "Type": "SOMATIC"
    })
    df_txt_final.to_csv(output, sep="\t", index=False)
    print(f"✅ Saved {len(df_txt_final)} mutations to {output}")
