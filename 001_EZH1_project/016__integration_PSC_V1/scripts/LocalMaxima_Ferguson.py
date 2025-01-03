import os
import pandas as pd
import numpy as np
from scipy.signal import argrelextrema

# Define your input list of sample names
input_list = [
    "PSC_KOEF1aEZH1_EZH1_005R",
    "PSC_KOEF1aEZH1_EZH1_006R",
    "PSC_KOEF1aEZH1_EZH1_013R1",
    "PSC_KOEF1aEZH1_EZH2_006R",
    "PSC_KOEF1aEZH1_EZH2_013R1",
    "PSC_KOEF1aEZH1_EZH2_014R1",
    "PSC_KOEF1aEZH1_H3K27me3_005R",
    "PSC_KOEF1aEZH1_H3K27me3_006R",
    "PSC_KOEF1aEZH1_H3K27me3_013R1",
    "PSC_KOEF1aEZH1_IGG_005R",
    "PSC_KOEF1aEZH1_IGG_006R",
    "PSC_KOEF1aEZH1_IGG_013R1",
    "PSC_KOEF1aEZH1_IGG_014R1",
    "PSC_KOEF1aEZH1_SUZ12_005R",
    "PSC_KOEF1aEZH1_SUZ12_006R",
    "PSC_KOEF1aEZH1_SUZ12_013R1",
    "PSC_KO_EZH1_006R",
    "PSC_KO_EZH1_013R1",
    "PSC_KO_EZH1_014R2",
    "PSC_KO_EZH2_013R1",
    "PSC_KO_EZH2_014R1",
    "PSC_KO_EZH2_014R2",
    "PSC_KO_H3K27me3_006R",
    "PSC_KO_H3K27me3_013R1",
    "PSC_KO_H3K27me3_014R2",
    "PSC_KO_IGG_006R",
    "PSC_KO_IGG_013R1",
    "PSC_KO_IGG_014R1",
    "PSC_KO_IGG_014R2",
    "PSC_KO_SUZ12_013R1",
    "PSC_KO_SUZ12_014R1",
    "PSC_KO_SUZ12_014R2",
    "PSC_WT_EZH1_006R",
    "PSC_WT_EZH2_006R",
    "PSC_WT_EZH2_010R",
    "PSC_WT_EZH2_014R1",
    "PSC_WT_H3K27me3_006R",
    "PSC_WT_H3K27me3_010R",
    "PSC_WT_H3K27me3_013R1",
    "PSC_WT_IGG_006R",
    "PSC_WT_IGG_010R",
    "PSC_WT_IGG_013R1",
    "PSC_WT_IGG_014R1",
    "PSC_WT_SUZ12_006R",
    "PSC_WT_SUZ12_013R1",
    "PSC_WT_SUZ12_014R1"
]

# Input and output directories


input_dir = "output/bigwig"
output_dir = "output/bigwig"
os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

# Function to process each sample
def process_sample(sample_name):
    input_file = os.path.join(input_dir, f"{sample_name}.filter.blacklist.bedGraph")
    output_file = os.path.join(output_dir, f"{sample_name}_local_maxima.bed")
    
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



