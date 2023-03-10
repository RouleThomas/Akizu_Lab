#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

outdir="input"

## ESC
x="ESC_WT_R1_1"
raw_f="input/E_WT_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R1_2"
raw_f="input/E_WT_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R2_1"
raw_f="input/E_WT_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R2_2"
raw_f="input/E_WT_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R3_1"
raw_f="input/E_WT_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_WT_R3_2"
raw_f="input/E_WT_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi
# KO
x="ESC_KO_R1_1"
raw_f="input/E_KO_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R1_2"
raw_f="input/E_KO_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R2_1"
raw_f="input/E_KO_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R2_2"
raw_f="input/E_KO_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R3_1"
raw_f="input/E_KO_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_KO_R3_2"
raw_f="input/E_KO_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# HET
x="ESC_HET_R1_1"
raw_f="input/E_HET_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R1_2"
raw_f="input/E_HET_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R2_1"
raw_f="input/E_HET_2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R2_2"
raw_f="input/E_HET_2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R3_1"
raw_f="input/E_HET_3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="ESC_HET_R3_2"
raw_f="input/E_HET_3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi