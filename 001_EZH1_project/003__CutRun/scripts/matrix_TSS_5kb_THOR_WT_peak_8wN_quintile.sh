#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8


computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_quintile_8wN_WT_peak_notExpress.gtf meta/ENCFF159KBI_quintile_8wN_WT_peak_quint1.gtf meta/ENCFF159KBI_quintile_8wN_WT_peak_quint2.gtf meta/ENCFF159KBI_quintile_8wN_WT_peak_quint3.gtf meta/ENCFF159KBI_quintile_8wN_WT_peak_quint4.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_THOR_WT_peak_8wN_quintile.gz \
    -p 8 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_THOR_WT_peak_8wN_quintile.bed


plotProfile -m output/deeptools/matrix_TSS_5kb_THOR_WT_peak_8wN_quintile.gz \
    -out output/deeptools/matrix_TSS_5kb_THOR_WT_peak_8wN_quintile_profile.pdf




