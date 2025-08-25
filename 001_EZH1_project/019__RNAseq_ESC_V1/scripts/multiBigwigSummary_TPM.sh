#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b \
    output/bigwig/ESC_WT_R1.bw \
    output/bigwig/ESC_WT_R2.bw \
    output/bigwig/ESC_WT_R3.bw \
    output/bigwig/ESC_KO_R1.bw \
    output/bigwig/ESC_KO_R2.bw \
    output/bigwig/ESC_KO_R3.bw \
    output/bigwig/ESC_OEKO_R1.bw \
    output/bigwig/ESC_OEKO_R2.bw \
    output/bigwig/ESC_OEKO_R3.bw \
 -o output/bigwig/multiBigwigSummary_TPM.npz


