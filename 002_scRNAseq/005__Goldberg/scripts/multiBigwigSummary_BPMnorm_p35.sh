#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/WT_p35_CB_Rep1_BPM.bw \
    output/bigwig/WT_p35_CB_Rep2_BPM.bw \
    output/bigwig/WT_p35_CB_Rep3_BPM.bw \
    output/bigwig/Kcnc1_p35_CB_Rep1_BPM.bw \
    output/bigwig/Kcnc1_p35_CB_Rep2_BPM.bw \
    output/bigwig/Kcnc1_p35_CB_Rep3_BPM.bw -o output/bigwig/multiBigwigSummary_BPMnorm_p35.npz








