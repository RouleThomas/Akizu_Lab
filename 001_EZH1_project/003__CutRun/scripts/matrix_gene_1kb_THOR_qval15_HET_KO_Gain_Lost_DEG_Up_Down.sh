#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00

### HET Gain; DEG down
computeMatrix scale-regions \
    -b 2000 -a 500 \
    -R output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Down.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Gain_DEG_Down.gz \
    -p 6
plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Gain_DEG_Down.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Gain_DEG_Down_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorList blue,white,red \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15






### HET Gain; DEG down ALL GENOTYPE
computeMatrix scale-regions \
    -b 2000 -a 500 \
    -R output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Down.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Gain_DEG_Down_allGeno.gz \
    -p 6
plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Gain_DEG_Down_allGeno.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Gain_DEG_Down_allGeno_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorList blue,white,red \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15




### HET Lost; DEG down
computeMatrix scale-regions \
    -b 2000 -a 500 \
    -R output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Down.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Lost_DEG_Down_allGeno.gz \
    -p 6
plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Lost_DEG_Down_allGeno.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Lost_DEG_Down_allGeno_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorList blue,white,red \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15


### HET Lost; DEG Up
computeMatrix scale-regions \
    -b 2000 -a 500 \
    -R output/ChIPseeker/THOR_qval15_HET_Lost_DEG_Up.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Lost_DEG_Up_allGeno.gz \
    -p 6
plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Lost_DEG_Up_allGeno.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_HET_Lost_DEG_Up_allGeno_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorList blue,white,red \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15


### KO Gain; DEG Up
computeMatrix scale-regions \
    -b 2000 -a 500 \
    -R output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Up.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Gain_DEG_Up_allGeno.gz \
    -p 6
plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Gain_DEG_Up_allGeno.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Gain_DEG_Up_allGeno_heatmap.png \
    --samplesLabel "WT" "KO" "HET" \
    --colorList blue,white,red \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15



### KO Gain; DEG Down
computeMatrix scale-regions \
    -b 2000 -a 500 \
    -R output/ChIPseeker/THOR_qval15_KO_Gain_DEG_Down.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Gain_DEG_Down_allGeno.gz \
    -p 6
plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Gain_DEG_Down_allGeno.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Gain_DEG_Down_allGeno_heatmap.png \
    --samplesLabel "WT" "KO" "HET" \
    --colorList blue,white,red \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15



### KO Lost; DEG Up
computeMatrix scale-regions \
    -b 2000 -a 500 \
    -R output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Up.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Lost_DEG_Up_allGeno.gz \
    -p 6
plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Lost_DEG_Up_allGeno.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Lost_DEG_Up_allGeno_heatmap.png \
    --samplesLabel "WT" "KO" "HET" \
    --colorList blue,white,red \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15

### KO Lost; DEG Down
computeMatrix scale-regions \
    -b 2000 -a 500 \
    -R output/ChIPseeker/THOR_qval15_KO_Lost_DEG_Down.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Lost_DEG_Down_allGeno.gz \
    -p 6
plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Lost_DEG_Down_allGeno.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_KO_Lost_DEG_Down_allGeno_heatmap.png \
    --samplesLabel "WT" "KO" "HET" \
    --colorList blue,white,red \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15