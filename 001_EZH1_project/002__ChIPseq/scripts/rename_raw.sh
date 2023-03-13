#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

outdir="input"

x="NPC_WT_R1_1"
raw_f="input/P_WT_1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="NPC_WT_R1_2"
raw_f="input/P_WT_1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi