import sys
import pandas as pd
import numpy as np
from scipy.io import mmread
import scrublet as scr

def run_scrublet(cellranger_path, output_path, manual_threshold=None):
    # Load the data
    matrix = mmread(f"{cellranger_path}/matrix.mtx.gz").T
    barcodes = pd.read_csv(f"{cellranger_path}/barcodes.tsv.gz", header=None, sep='\t')[0]

    # Run Scrublet
    scrub = scr.Scrublet(matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    if manual_threshold is not None:
        predicted_doublets = doublet_scores > manual_threshold

    # Create a DataFrame with the results
    results = pd.DataFrame({
        'Doublet_score': doublet_scores,
        'Is_doublet': predicted_doublets
    }, index=barcodes)

    # Write the results to a TSV file
    results.to_csv(output_path, sep='\t')

if __name__ == '__main__':
    if len(sys.argv) not in [3, 4]:
        print("Usage: python3 script.py <input_path> <output_path> [manual_threshold]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    manual_threshold = float(sys.argv[3]) if len(sys.argv) == 4 else None

    run_scrublet(input_path, output_path, manual_threshold)


