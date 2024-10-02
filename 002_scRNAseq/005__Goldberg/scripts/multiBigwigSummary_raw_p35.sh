#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/WT_p35_CB_Rep1.bw \
    output/bigwig/WT_p35_CB_Rep2.bw \
    output/bigwig/WT_p35_CB_Rep3.bw \
    output/bigwig/Kcnc1_p35_CB_Rep1.bw \
    output/bigwig/Kcnc1_p35_CB_Rep2.bw \
    output/bigwig/Kcnc1_p35_CB_Rep3.bw -o output/bigwig/multiBigwigSummary_raw_p35.npz








