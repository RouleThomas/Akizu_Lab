import pandas as pd
import os
import numpy as np


# Input list of sample names
input_list = [
    "PSC_KOEF1aEZH1_IGG_005R",
    "PSC_KOEF1aEZH1_IGG_006R",
    "PSC_KOEF1aEZH1_IGG_013R1",
    "PSC_KOEF1aEZH1_IGG_014R1",
    "PSC_KO_IGG_006R",
    "PSC_KO_IGG_013R1",
    "PSC_KO_IGG_014R1",
    "PSC_KO_IGG_014R2",
    "PSC_WT_IGG_006R",
    "PSC_WT_IGG_010R",
    "PSC_WT_IGG_013R1",
    "PSC_WT_IGG_014R1",
]


# Input file containing 99th percentile values
percentile_file = "output/bigwig/Ferguson_SF_95.txt"
local_maxima_dir = "output/bigwig"
normalized_dir = "output/bigwig_Ferguson"
os.makedirs(normalized_dir, exist_ok=True)

# Reference sample
reference_sample = "PSC_WT_IGG_006R"

# Load the 99th percentile values from the file
percentile_dict = {}
with open(percentile_file, "r") as f:
    for line in f:
        sample, percentile = line.strip().split("\t")
        percentile_dict[sample] = float(percentile)

# Get the reference value
reference_value = percentile_dict.get(reference_sample)
if reference_value is None:
    print(f"Reference sample {reference_sample} not found in percentile_dict.")
    exit()

# Normalize only the samples in input_list
for sample_name in input_list:
    percentile_value = percentile_dict.get(sample_name)
    if percentile_value is None:
        print(f"Sample {sample_name} not found in percentile_dict. Skipping.")
        continue
    
    scaling_factor = reference_value / percentile_value
    
    input_file = os.path.join(local_maxima_dir, f"{sample_name}_local_maxima.bed")
    output_file = os.path.join(normalized_dir, f"{sample_name}_norm95.bedGraph")
    
    if os.path.exists(input_file):
        # Load local maxima file
        data = pd.read_csv(input_file, sep="\t", header=None, names=["chrom", "start", "end", "score"])
        
        # Apply scaling factor
        data["score"] *= scaling_factor
        
        # Save normalized data
        data.to_csv(output_file, sep="\t", index=False, header=False)
        print(f"Normalized {sample_name}: Scaling factor = {scaling_factor}")
    else:
        print(f"File not found: {input_file}. Skipping.")