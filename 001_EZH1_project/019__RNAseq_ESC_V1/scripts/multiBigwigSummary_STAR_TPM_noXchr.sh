#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b \
    output/bigwig_STAR/ESC_WT_R1_noXchr.bw \
    output/bigwig_STAR/ESC_WT_R2_noXchr.bw \
    output/bigwig_STAR/ESC_WT_R3_noXchr.bw \
    output/bigwig_STAR/ESC_KO_R1_noXchr.bw \
    output/bigwig_STAR/ESC_KO_R2_noXchr.bw \
    output/bigwig_STAR/ESC_KO_R3_noXchr.bw \
    output/bigwig_STAR/ESC_OEKO_R1_noXchr.bw \
    output/bigwig_STAR/ESC_OEKO_R2_noXchr.bw \
    output/bigwig_STAR/ESC_OEKO_R3_noXchr.bw \
 -o output/bigwig_STAR/multiBigwigSummary_STAR_noXchr_TPM.npz


