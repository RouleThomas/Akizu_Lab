import os
import pandas as pd
import numpy as np
from scipy.signal import argrelextrema

# Define your input list of sample names
input_list = [
  # WT
  "ESC_WT_EZH1_R1", "ESC_WT_EZH1_R2", "ESC_WT_EZH1_R3",
  "ESC_WT_EZH2_R1", "ESC_WT_EZH2_R2", "ESC_WT_EZH2_R3",
  "ESC_WT_H3K27me3_R1", "ESC_WT_H3K27me3_R2", "ESC_WT_H3K27me3_R3",
  "ESC_WT_IGG_R1", "ESC_WT_IGG_R2", "ESC_WT_IGG_R3",

  # KO
  "ESC_KO_EZH1_R1", "ESC_KO_EZH1_R2", "ESC_KO_EZH1_R3",
  "ESC_KO_EZH2_R1", "ESC_KO_EZH2_R2", "ESC_KO_EZH2_R3",
  "ESC_KO_H3K27me3_R1", "ESC_KO_H3K27me3_R2", "ESC_KO_H3K27me3_R3",

  # OEKO
  "ESC_OEKO_EZH1_R1", "ESC_OEKO_EZH1_R2", "ESC_OEKO_EZH1_R3",
  "ESC_OEKO_EZH2_R1", "ESC_OEKO_EZH2_R2", "ESC_OEKO_EZH2_R3",
  "ESC_OEKO_H3K27me3_R1", "ESC_OEKO_H3K27me3_R2", "ESC_OEKO_H3K27me3_R3"
]

# Input and output directories


input_dir = "output/bigwig"
output_dir = "output/bigwig"
os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

# Function to process each sample
def process_sample(sample_name):
    input_file = os.path.join(input_dir, f"{sample_name}.unique.dupmark.sorted.blacklist.bedGraph")
    output_file = os.path.join(output_dir, f"{sample_name}_unique_local_maxima.bed")
    
    if not os.path.exists(input_file):
        print(f"Input file for {sample_name} not found. Skipping.")
        return
    
    # Load the filtered bedGraph file
    data = pd.read_csv(input_file, sep="\t", header=None, names=["chrom", "start", "end", "score"])
    
    # Convert scores to a NumPy array
    scores = data['score'].values

    # Find local maxima
    local_maxima_indices = argrelextrema(scores, np.greater)[0]
    local_maxima = data.iloc[local_maxima_indices]

    # Save the local maxima to a new file
    local_maxima.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Processed {sample_name}: Local maxima saved to {output_file}")

# Loop through all samples
for sample in input_list:
    process_sample(sample)



