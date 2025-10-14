import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
import numpy as np


# Set base directory
BASE_DIR = Path("results")

# Define score columns and their pretty names
SCORE_COLUMNS = {
    "sift4g_score": "SIFT4G Score",
    "polyphen2_hdiv_score": "PolyPhen2 HDIV Score",
    "cadd_phred": "CADD PHRED Score",
    "frac_syn": "Fraction Synonymous",
    "frac_missense": "Fraction Missense",
    "frac_stop_gained": "Fraction Stop Gained"
}

# Define how to extract binary scores from 'consequence'
def consequence_to_binary(df):
    df["frac_syn"] = (df["consequence"] == "synonymous").astype(int)
    df["frac_missense"] = (df["consequence"] == "missense").astype(int)
    df["frac_stop_gained"] = (df["consequence"] == "stop").astype(int)
    return df

# Iterate over SBS folders
for sbs_folder in sorted(BASE_DIR.glob("*/*")):
    pdf_path = sbs_folder / f"{sbs_folder.parts[-2]}_replicate_score_distributions.pdf"
    with PdfPages(pdf_path) as pdf:
        # Find all replicates
        for parquet_file in sorted(sbs_folder.glob("rep_*.annot.parquet")):
            df = pd.read_parquet(parquet_file)
            df = consequence_to_binary(df)

            fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8, 10))
            axs = axs.flatten()

            for i, (col, title) in enumerate(SCORE_COLUMNS.items()):
                ax = axs[i]
                if col not in df.columns:
                    ax.set_visible(False)
                    continue


                # Add jittered transparent points
                x_jitter = 0 + 0.1 * (np.random.rand(len(df)) - 0.5)  # small jitter around 0
                ax.scatter(x_jitter, df[col], color="black", alpha=0.2, s=6, edgecolors='none')

                ax.boxplot(df[col].dropna(), positions=[0], widths=0.5,
                           patch_artist=True, boxprops=dict(facecolor='purple', alpha=0.2))
                ax.set_title(title)
                ax.set_ylabel("Score")
                ax.set_xticks([])

                # For binary scores, add count annotations
                if col in {"frac_syn", "frac_missense", "frac_stop_gained"}:
                    count_0 = (df[col] == 0).sum()
                    count_1 = (df[col] == 1).sum()
                    ax.text(0, 0.1, f"0: {count_0}", ha='center', fontsize=8)
                    ax.text(0, 0.9, f"1: {count_1}", ha='center', fontsize=8)

            # Turn off any unused subplots
            for j in range(i + 1, len(axs)):
                axs[j].set_visible(False)

            replicate_name = parquet_file.stem
            fig.suptitle(f"{sbs_folder.parts[-2]} - {replicate_name}: Mutation Score Distributions", fontsize=12)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig(fig)
            plt.close()

    print(f"✅ Finished plotting {sbs_folder.parts[-2]}")
