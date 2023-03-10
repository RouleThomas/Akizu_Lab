#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

outdir="input"

### 4-8weeks RNAseq
x="4wN_WT_R1_1"
raw_f="input/R1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_WT_R1_2"
raw_f="input/R1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_WT_R2_1"
raw_f="input/R2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_WT_R2_2"
raw_f="input/R2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_KO_R1_1"
raw_f="input/R3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_KO_R1_2"
raw_f="input/R3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_KO_R2_1"
raw_f="input/R4_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_KO_R2_2"
raw_f="input/R4_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="4wN_HET_R1_1"
raw_f="input/R5_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R1_2"
raw_f="input/R5_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R2_1"
raw_f="input/R6_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R2_2"
raw_f="input/R6_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R3_1"
raw_f="input/R7_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R3_2"
raw_f="input/R7_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R4_1"
raw_f="input/R8_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_HET_R4_2"
raw_f="input/R8_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCWT_R1_1"
raw_f="input/R9_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCWT_R1_2"
raw_f="input/R9_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCWT_R2_1"
raw_f="input/R10_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCWT_R2_2"
raw_f="input/R10_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCpatient_R1_1"
raw_f="input/R11_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCpatient_R1_2"
raw_f="input/R11_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCpatient_R2_1"
raw_f="input/R12_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="4wN_iPSCpatient_R2_2"
raw_f="input/R12_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# NR

x="8wN_WT_R1_1"
raw_f="input/NR1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R1_2"
raw_f="input/NR1_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R2_1"
raw_f="input/NR2_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R2_2"
raw_f="input/NR2_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R3_1"
raw_f="input/NR3_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R3_2"
raw_f="input/NR3_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R4_1"
raw_f="input/NR4_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_WT_R4_2"
raw_f="input/NR4_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R1_1"
raw_f="input/NR5_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R1_2"
raw_f="input/NR5_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R2_1"
raw_f="input/NR6_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R2_2"
raw_f="input/NR6_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R3_1"
raw_f="input/NR7_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R3_2"
raw_f="input/NR7_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R4_1"
raw_f="input/NR8_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_KO_R4_2"
raw_f="input/NR8_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_HET_R1_1"
raw_f="input/NR9_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_HET_R1_2"
raw_f="input/NR9_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_HET_R2_1"
raw_f="input/NR10_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_HET_R2_2"
raw_f="input/NR10_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="8wN_HET_R3_1"
raw_f="input/NR11_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_HET_R3_2"
raw_f="input/NR11_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi


x="8wN_HET_R4_1"
raw_f="input/NR12_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_HET_R4_2"
raw_f="input/NR12_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_iPSCpatient_R1_1"
raw_f="input/NR13_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_iPSCpatient_R1_2"
raw_f="input/NR13_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_iPSCpatient_R2_1"
raw_f="input/NR14_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="8wN_iPSCpatient_R2_2"
raw_f="input/NR14_2.fq.gz"
new_f="${outdir}/${x}.fq.gz"
if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi