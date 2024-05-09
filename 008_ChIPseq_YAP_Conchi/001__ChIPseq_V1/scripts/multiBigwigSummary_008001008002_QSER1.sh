#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8



multiBigwigSummary bins -b ../002__ChIPseq_Dixon2021/output/bigwig/hESC_WT_QSER1FLAG_R1.unique.dupmark.sorted.bw \
../002__ChIPseq_Dixon2021/output/bigwig/hESC_WT_QSER1FLAG_R2.unique.dupmark.sorted.bw \
../002__ChIPseq_Dixon2021/output/bigwig/hESC_WT_input_R1.unique.dupmark.sorted.bw \
../002__ChIPseq_Dixon2021/output/bigwig/hESC_WT_input_R2.unique.dupmark.sorted.bw \
../001__ChIPseq_V1/output/bigwig/hESC_WT_QSER1_R1.unique.dupmark.sorted.bw \
../001__ChIPseq_V1/output/bigwig/hESC_WT_QSER1_R2.unique.dupmark.sorted.bw \
../001__ChIPseq_V1/output/bigwig/hESC_WT_input_R1.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_008001008002_QSER1.npz -p max



