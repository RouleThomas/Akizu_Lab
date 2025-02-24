import pandas as pd
import os
import numpy as np

# Define your input list of sample names
input_list = [
    "NPC_WT_H3K4me3_008",
    "NPC_WT_IGG_008",
    "NPC_WT_H3K4me3_005",
    "NPC_WT_IGG_005",
    "NPC_KO_H3K4me3_008",
    "NPC_KO_IGG_008",
    "NPC_KO_H3K4me3_005",
    "NPC_KO_IGG_005",
    "NPC_WT_H3K27me3_008",
    "NPC_WT_IGG_008",
    "NPC_WT_H3K27me3_005",
    "NPC_WT_IGG_005",
    "NPC_KO_H3K27me3_008",
    "NPC_KO_IGG_008",
    "NPC_KO_H3K27me3_005",
    "NPC_KO_IGG_005"
]

# Input and output directories


local_maxima_dir = "output/bigwig"
output_file = "output/bigwig/Ferguson_SF_99_unique.txt"


# Initialize a dictionary to store 99th percentile values
percentile_dict = {}

# Loop through each sample
for sample_name in input_list:
    input_file = os.path.join(local_maxima_dir, f"{sample_name}_unique_local_maxima.bed")
    
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}. Skipping.")
        continue
    
    # Load local maxima file
    data = pd.read_csv(input_file, sep="\t", header=None, names=["chrom", "start", "end", "score"])
    
    # Calculate the 99th percentile of the 'score' column
    percentile_99 = np.percentile(data["score"], 99)
    percentile_dict[sample_name] = percentile_99
    print(f"{sample_name}: 99th percentile = {percentile_99}")

# Save the 99th percentile values to a file
with open(output_file, "w") as f:
    for sample, percentile in percentile_dict.items():
        f.write(f"{sample}\t{percentile}\n")

print(f"99th percentile values saved to {output_file}")


