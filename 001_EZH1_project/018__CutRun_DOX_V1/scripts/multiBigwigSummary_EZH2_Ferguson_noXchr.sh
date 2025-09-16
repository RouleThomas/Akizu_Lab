#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b \
    output/bigwig_Ferguson/ESC_WT_EZH2_R1_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_WT_EZH2_R2_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_WT_EZH2_R3_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_KO_EZH2_R1_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_KO_EZH2_R2_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_KO_EZH2_R3_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_noXchr_unique_norm99_initialBigwig.bw \
    output/bigwig/ESC_WT_IGG_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R3.unique.dupmark.sorted.bw \
 -o output/bigwig_Ferguson/multiBigwigSummary_EZH2_Ferguson_noXchr.npz

