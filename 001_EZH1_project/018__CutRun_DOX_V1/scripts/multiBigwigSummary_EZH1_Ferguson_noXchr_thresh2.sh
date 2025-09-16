#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b \
    output/bigwig_Ferguson/ESC_WT_EZH1_R1_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig_Ferguson/ESC_WT_EZH1_R2_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig_Ferguson/ESC_WT_EZH1_R3_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig_Ferguson/ESC_KO_EZH1_R1_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig_Ferguson/ESC_KO_EZH1_R2_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig_Ferguson/ESC_KO_EZH1_R3_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R1_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R2_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R3_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    output/bigwig/ESC_WT_IGG_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R3.unique.dupmark.sorted.bw \
 -o output/bigwig_Ferguson/multiBigwigSummary_EZH1_Ferguson_noXchr_thresh2.npz

