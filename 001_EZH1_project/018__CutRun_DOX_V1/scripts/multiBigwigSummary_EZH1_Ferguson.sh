#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b \
    output/bigwig_Ferguson/ESC_WT_EZH1_R1_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_WT_EZH1_R2_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_WT_EZH1_R3_unique_norm99_initialBigwig.bw \
    \
    output/bigwig_Ferguson/ESC_KO_EZH1_R1_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_KO_EZH1_R2_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_KO_EZH1_R3_unique_norm99_initialBigwig.bw \
    \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R1_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R2_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R3_unique_norm99_initialBigwig.bw \
    output/bigwig/ESC_WT_IGG_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R3.unique.dupmark.sorted.bw \
 -o output/bigwig_Ferguson/multiBigwigSummary_EZH1_Ferguson.npz

