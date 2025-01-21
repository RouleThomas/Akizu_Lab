#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig/PSC_KOEF1aEZH1_EZH2_006R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KOEF1aEZH1_EZH2_013R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KOEF1aEZH1_EZH2_014R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_EZH2_013R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_EZH2_014R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_EZH2_014R2.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_EZH2_006R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_EZH2_010R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_EZH2_014R1.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_EZH2_unique_raw.npz



