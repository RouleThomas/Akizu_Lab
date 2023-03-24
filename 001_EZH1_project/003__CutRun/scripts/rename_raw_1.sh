#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

# Set output directory
output_dir="input"

# Function to rename files
rename_file() {
    sample_name="$1"
    input_file="$2"
    output_file="${output_dir}/${sample_name}.fq.gz"
    
    if [[ -f "$input_file" && ! -f "$output_file" ]]; then
        mv "$input_file" "$output_file"
    elif [[ ! -f "$input_file" && ! -f "$output_file" ]]; then
        echo "ERROR: Cannot Find File: ${input_file}"
        exit
    fi
}

# Rename files for HET IGG samples
rename_file "8wN_HET_IGG_R1_1" "input/R1_Het6_Igg_CKDL220021576-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_IGG_R1_2" "input/R1_Het6_Igg_CKDL220021576-1A_HF3THDSX5_L1_2.fq.gz"
rename_file "8wN_HET_IGG_R2_1" "input/R2_Het6_Igg_CKDL220021578-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_IGG_R2_2" "input/R2_Het6_Igg_CKDL220021578-1A_HF3THDSX5_L1_2.fq.gz"

# Rename files for HET H3K27me3 samples
rename_file "8wN_HET_H3K27me3_R1_1" "input/R1_Het6_K27me3_CKDL220021577-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_H3K27me3_R1_2" "input/R1_Het6_K27me3_CKDL220021577-1A_HF3THDSX5_L1_2.fq.gz"
rename_file "8wN_HET_H3K27me3_R2_1" "input/R2_Het6_K27me3_CKDL220021579-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_H3K27me3_R2_2" "input/R2_Het6_K27me3_CKDL220021579-1A_HF3THDSX5_L1_2.fq.gz"

# Rename files for KO IGG samples
rename_file "8wN_KO_IGG_R1_1" "input/R1_KO8bp_Igg_CKDL220021564-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_IGG_R1_2" "input/R1_KO8bp_Igg_CKDL220021564-1A_HF3THDSX5_L1_2.fq.gz"
rename_file "8wN_KO_IGG_R2_1" "input/R2_KO8bp_Igg_CKDL220021566-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_IGG_R2_2" "input/R2_KO8bp_Igg_CKDL220021566-1A_HF3THDSX5_L1_2.fq.gz"
