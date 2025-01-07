import pandas as pd
import os
import numpy as np

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


local_maxima_dir = "output/bigwig"
output_file = "output/bigwig/Ferguson_SF_98.txt"


# Initialize a dictionary to store 98th percentile values
percentile_dict = {}

# Loop through each sample
for sample_name in input_list:
    input_file = os.path.join(local_maxima_dir, f"{sample_name}_local_maxima.bed")
    
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}. Skipping.")
        continue
    
    # Load local maxima file
    data = pd.read_csv(input_file, sep="\t", header=None, names=["chrom", "start", "end", "score"])
    
    # Calculate the 98th percentile of the 'score' column
    percentile_98 = np.percentile(data["score"], 98)
    percentile_dict[sample_name] = percentile_98
    print(f"{sample_name}: 98th percentile = {percentile_98}")

# Save the 98th percentile values to a file
with open(output_file, "w") as f:
    for sample, percentile in percentile_dict.items():
        f.write(f"{sample}\t{percentile}\n")

print(f"98th percentile values saved to {output_file}")


