#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig/PSC_KOEF1aEZH1_SUZ12_005R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KOEF1aEZH1_SUZ12_006R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KOEF1aEZH1_SUZ12_013R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_SUZ12_013R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_SUZ12_014R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_SUZ12_014R2.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_SUZ12_006R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_SUZ12_013R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_SUZ12_014R1.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_SUZ12_unique_raw.npz





