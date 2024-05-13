#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8



multiBigwigSummary bins -b ../002__ChIPseq_Dixon2021/output/bigwig_DiffBindTMM/hESC_WT_QSER1FLAG_R1.unique.dupmark.sorted.bw \
../002__ChIPseq_Dixon2021/output/bigwig_DiffBindTMM/hESC_WT_QSER1FLAG_R2.unique.dupmark.sorted.bw \
../002__ChIPseq_Dixon2021/output/bigwig_DiffBindTMM/hESC_WT_input_R1.unique.dupmark.sorted.bw \
../002__ChIPseq_Dixon2021/output/bigwig_DiffBindTMM/hESC_WT_input_R2.unique.dupmark.sorted.bw \
../001__ChIPseq_V1/output/bigwig_DiffBindTMM/hESC_WT_QSER1_R1.unique.dupmark.sorted.bw \
../001__ChIPseq_V1/output/bigwig_DiffBindTMM/hESC_WT_QSER1_R2.unique.dupmark.sorted.bw \
../001__ChIPseq_V1/output/bigwig_DiffBindTMM/hESC_WT_input_R1.unique.dupmark.sorted.bw -o output/bigwig_DiffBindTMM/multiBigwigSummary_008001008002_QSER1_DiffBindTMM_bin2kb.npz -p max --blackListFileName /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --chromosomesToSkip chrX chrY chrM --binSize 2000



