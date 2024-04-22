#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8



multiBigwigSummary bins -b output/bigwig/CPC_untreated_YAP1_R1.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_YAP1_R1.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_YAP1_R2.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_YAP1_R2.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_TEAD4_R1.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_TEAD4_R1.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_TEAD4_R2.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_TEAD4_R2.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_NR2F2_R1.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_NR2F2_R1.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_NR2F2_R2.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_NR2F2_R2.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_NR2F2_R3.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_NR2F2_R3.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_input_R3.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_input_R3.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_YAP1_R3.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_YAP1_R3.unique.dupmark.sorted.bw \
output/bigwig/CPC_untreated_TEAD4_R3.unique.dupmark.sorted.bw \
output/bigwig/CPC_RA_TEAD4_R3.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_extendReads_CPC.npz -p max








