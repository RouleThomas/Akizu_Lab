import sys
import pandas as pd
import numpy as np
from scipy.io import mmread
import scrublet as scr
import matplotlib.pyplot as plt

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

    # Save the results to a TSV file
    results.to_csv(output_path, sep='\t')

    # Plot the histogram of doublet scores
    plt.hist(doublet_scores, bins=50)
    plt.xlabel('Doublet Score')
    plt.ylabel('Number of Cells')
    plt.title('Distribution of Doublet Scores')
    plt.savefig(output_path.replace('.tsv', '_doublet_scores_histogram.png'))

if __name__ == '__main__':
    if len(sys.argv) not in [3, 4]:
        print("Usage: python3 script.py <input_path> <output_path> [manual_threshold]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    manual_threshold = float(sys.argv[3]) if len(sys.argv) == 4 else None

    run_scrublet(input_path, output_path, manual_threshold)
