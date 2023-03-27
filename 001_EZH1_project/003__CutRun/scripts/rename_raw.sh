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

# Apply the function 

rename_file "8wN_HET_IGG_R1_1"	"input/R1_Het5_Igg_CKDL220021572-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_IGG_R1_2"	"input/R1_Het5_Igg_CKDL220021572-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_HET_H3K27me3_R1_1"	"input/R1_Het5_K27me3_CKDL220021573-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_H3K27me3_R1_2"	"input/R1_Het5_K27me3_CKDL220021573-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_HET_IGG_R2_1"	"input/R1_Het6_Igg_CKDL220021576-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_IGG_R2_2"	"input/R1_Het6_Igg_CKDL220021576-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_HET_H3K27me3_R2_1"	"input/R1_Het6_K27me3_CKDL220021577-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_H3K27me3_R2_2"	"input/R1_Het6_K27me3_CKDL220021577-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_KO_IGG_R1_1"	"input/R1_KO16bp_Igg_CKDL220021568-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_IGG_R1_2"	"input/R1_KO16bp_Igg_CKDL220021568-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_KO_H3K27me3_R1_1"	"input/R1_KO16bp_K27me3_CKDL220021569-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_H3K27me3_R1_2"	"input/R1_KO16bp_K27me3_CKDL220021569-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_KO_IGG_R2_1"	"input/R1_KO8bp_Igg_CKDL220021564-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_IGG_R2_2"	"input/R1_KO8bp_Igg_CKDL220021564-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_KO_H3K27me3_R2_1"	"input/R1_KO8bp_K27me3_CKDL220021565-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_H3K27me3_R2_2"	"input/R1_KO8bp_K27me3_CKDL220021565-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_iPSCpatient_IGG_R1_1"	"input/R1_Splc11_Igg_CKDL220021580-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_iPSCpatient_IGG_R1_2"	"input/R1_Splc11_Igg_CKDL220021580-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_iPSCpatient_H3K27me3_R1_1"	"input/R1_Splc11_K27me3_CKDL220021581-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_iPSCpatient_H3K27me3_R1_2"	"input/R1_Splc11_K27me3_CKDL220021581-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_WT_IGG_R1_1"	"input/R1_WTH9_Igg_CKDL220021556-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_WT_IGG_R1_2"	"input/R1_WTH9_Igg_CKDL220021556-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_WT_H3K27me3_R1_1"	"input/R1_WTH9_K27me3_CKDL220021557-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_WT_H3K27me3_R1_2"	"input/R1_WTH9_K27me3_CKDL220021557-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_WT_IGG_R2_1"	"input/R1_WTS2_Igg_CKDL220021560-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_WT_IGG_R2_2"	"input/R1_WTS2_Igg_CKDL220021560-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_WT_H3K27me3_R2_1"	"input/R1_WTS2_K27me3_CKDL220021561-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_WT_H3K27me3_R2_2"	"input/R1_WTS2_K27me3_CKDL220021561-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_HET_IGG_R3_1"	"input/R2_Het5_Igg_CKDL220021574-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_IGG_R3_2"	"input/R2_Het5_Igg_CKDL220021574-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_HET_H3K27me3_R3_1"	"input/R2_Het5_K27me3_CKDL220021575-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_H3K27me3_R3_2"	"input/R2_Het5_K27me3_CKDL220021575-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_HET_IGG_R4_1"	"input/R2_Het6_Igg_CKDL220021578-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_IGG_R4_2"	"input/R2_Het6_Igg_CKDL220021578-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_HET_H3K27me3_R4_1"	"input/R2_Het6_K27me3_CKDL220021579-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_HET_H3K27me3_R4_2"	"input/R2_Het6_K27me3_CKDL220021579-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_KO_IGG_R3_1"	"input/R2_KO16bp_Igg_CKDL220021570-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_IGG_R3_2"	"input/R2_KO16bp_Igg_CKDL220021570-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_KO_H3K27me3_R3_1"	"input/R2_KO16bp_K27me3_CKDL220021571-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_H3K27me3_R3_2"	"input/R2_KO16bp_K27me3_CKDL220021571-1A_HF3THDSX5_L1_2.fq.gz"

rename_file "8wN_KO_IGG_R4_1"	"input/R2_KO8bp_Igg_CKDL220021566-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_IGG_R4_2"	"input/R2_KO8bp_Igg_CKDL220021566-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_KO_H3K27me3_R4_1"	"input/R2_KO8bp_K27me3_CKDL220021567-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_KO_H3K27me3_R4_2"	"input/R2_KO8bp_K27me3_CKDL220021567-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_iPSCpatient_IGG_R2_1"	"input/R2_Splc11_Igg_CKDL220021582-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_iPSCpatient_IGG_R2_2"	"input/R2_Splc11_Igg_CKDL220021582-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_iPSCpatient_H3K27me3_R2_1"	"input/R2_Splc11_K27me3_CKDL220021583-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_iPSCpatient_H3K27me3_R2_2"	"input/R2_Splc11_K27me3_CKDL220021583-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_WT_IGG_R3_1"	"input/R2_WTH9_Igg_CKDL220021558-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_WT_IGG_R3_2"	"input/R2_WTH9_Igg_CKDL220021558-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_WT_H3K27me3_R3_1"	"input/R2_WTH9_K27me3_CKDL220021559-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_WT_H3K27me3_R3_2"	"input/R2_WTH9_K27me3_CKDL220021559-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_WT_IGG_R4_1"	"input/R2_WTS2_Igg_CKDL220021562-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_WT_IGG_R4_2"	"input/R2_WTS2_Igg_CKDL220021562-1A_HF3THDSX5_L1_2.fq.gz"


rename_file "8wN_WT_H3K27me3_R4_1"	"input/R2_WTS2_K27me3_CKDL220021563-1A_HF3THDSX5_L1_1.fq.gz"
rename_file "8wN_WT_H3K27me3_R4_2"	"input/R2_WTS2_K27me3_CKDL220021563-1A_HF3THDSX5_L1_2.fq.gz"


