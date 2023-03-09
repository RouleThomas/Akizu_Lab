#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

outdir="input"

x="NPC_WT_R1_1"
raw_f="P_WT_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_WT_R1_2"
raw_f="P_WT_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_WT_R2_1"
raw_f="P_WT_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_WT_R2_2"
raw_f="P_WT_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_WT_R3_1"
raw_f="P_WT_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_WT_R3_2"
raw_f="P_WT_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi
# KO
x="NPC_KO_R1_1"
raw_f="P_KO_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_KO_R1_2"
raw_f="P_KO_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_KO_R2_1"
raw_f="P_KO_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_KO_R2_2"
raw_f="P_KO_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_KO_R3_1"
raw_f="P_KO_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_KO_R3_2"
raw_f="P_KO_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# HET
x="NPC_HET_R1_1"
raw_f="P_HET_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_HET_R1_2"
raw_f="P_HET_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_HET_R2_1"
raw_f="P_HET_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_HET_R2_2"
raw_f="P_HET_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_HET_R3_1"
raw_f="P_HET_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_HET_R3_2"
raw_f="P_HET_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

## 2dN
x="2dN_WT_R1_1"
raw_f="N_WT_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_WT_R1_2"
raw_f="N_WT_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_WT_R2_1"
raw_f="N_WT_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_WT_R2_2"
raw_f="N_WT_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_WT_R3_1"
raw_f="N_WT_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_WT_R3_2"
raw_f="N_WT_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi
# KO
x="2dN_KO_R1_1"
raw_f="N_KO_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_KO_R1_2"
raw_f="N_KO_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_KO_R2_1"
raw_f="N_KO_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_KO_R2_2"
raw_f="N_KO_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_KO_R3_1"
raw_f="N_KO_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_KO_R3_2"
raw_f="N_KO_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# HET
x="2dN_HET_R1_1"
raw_f="N_HET_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_HET_R1_2"
raw_f="N_HET_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_HET_R2_1"
raw_f="N_HET_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_HET_R2_2"
raw_f="N_HET_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_HET_R3_1"
raw_f="N_HET_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="2dN_HET_R3_2"
raw_f="N_HET_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

## ESC
x="ESC_WT_R1_1"
raw_f="E_WT_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R1_2"
raw_f="E_WT_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R2_1"
raw_f="E_WT_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R2_2"
raw_f="E_WT_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R3_1"
raw_f="E_WT_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R3_2"
raw_f="E_WT_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi
# KO
x="ESC_KO_R1_1"
raw_f="E_KO_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R1_2"
raw_f="E_KO_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R2_1"
raw_f="E_KO_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R2_2"
raw_f="E_KO_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R3_1"
raw_f="E_KO_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R3_2"
raw_f="E_KO_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# HET
x="ESC_HET_R1_1"
raw_f="E_HET_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R1_2"
raw_f="E_HET_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R2_1"
raw_f="E_HET_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R2_2"
raw_f="E_HET_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R3_1"
raw_f="E_HET_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R3_2"
raw_f="E_HET_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

### 4-8weeks RNAseq
x="4wN_WT_R1_1"
raw_f="R1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_WT_R1_2"
raw_f="R1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_WT_R2_1"
raw_f="R2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_WT_R2_2"
raw_f="R2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_KO_R1_1"
raw_f="R3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_KO_R1_2"
raw_f="R3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_KO_R2_1"
raw_f="R4_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_KO_R2_2"
raw_f="R4_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="4wN_HET_R1_1"
raw_f="R5_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R1_2"
raw_f="R5_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R2_1"
raw_f="R6_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R2_2"
raw_f="R6_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R3_1"
raw_f="R7_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R3_2"
raw_f="R7_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R4_1"
raw_f="R8_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R4_2"
raw_f="R8_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCWT_R1_1"
raw_f="R9_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCWT_R1_2"
raw_f="R9_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCWT_R2_1"
raw_f="R10_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCWT_R2_2"
raw_f="R10_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCpatient_R1_1"
raw_f="R11_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCpatient_R1_2"
raw_f="R11_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCpatient_R2_1"
raw_f="R12_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCpatient_R2_2"
raw_f="R12_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# NR

x="8wN_WT_R1_1"
raw_f="NR1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R1_2"
raw_f="NR1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R2_1"
raw_f="NR2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R2_2"
raw_f="NR2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R3_1"
raw_f="NR3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R3_2"
raw_f="NR3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R4_1"
raw_f="NR4_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R4_2"
raw_f="NR4_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R1_1"
raw_f="NR5_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R1_2"
raw_f="NR5_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi
