#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8





multiBigwigSummary bins -b output/bigwig/hESC_WT_DVL2_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_YAPKO_DVL2_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_EZH2_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_YAPKO_EZH2_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_EZH2_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_YAPKO_EZH2_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_QSER1_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_YAPKO_QSER1_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_QSER1_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_YAPKO_QSER1_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_input_R1.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_extendReads_hESC.npz -p max








