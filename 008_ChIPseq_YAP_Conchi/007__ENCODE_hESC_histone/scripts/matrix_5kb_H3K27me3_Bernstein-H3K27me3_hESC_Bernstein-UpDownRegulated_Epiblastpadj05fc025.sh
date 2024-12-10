#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_H3K27me3_hESC_Bernstein-UpRegulated_Epiblastpadj05fc025.gtf meta/ENCFF159KBI_H3K27me3_hESC_Bernstein-DownRegulated_Epiblastpadj05fc025.gtf \
    -S output/ENCODE/Bernstein_H3K27me3.bigWig \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_5kb_H3K27me3_Bernstein-H3K27me3_hESC_Bernstein-UpDownRegulated_Epiblastpadj05fc025.gz \
    -p 6




plotHeatmap -m output/deeptools/matrix_5kb_H3K27me3_Bernstein-H3K27me3_hESC_Bernstein-UpDownRegulated_Epiblastpadj05fc025.gz \
    -out output/deeptools/matrix_5kb_H3K27me3_Bernstein-H3K27me3_hESC_Bernstein-UpDownRegulated_Epiblastpadj05fc025_heatmap.pdf \
    --samplesLabel "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 4



# interactive
plotHeatmap -m output/deeptools/matrix_5kb_H3K27me3_Bernstein-H3K27me3_hESC_Bernstein-UpDownRegulated_Epiblastpadj05fc025.gz \
    -out output/deeptools/matrix_5kb_H3K27me3_Bernstein-H3K27me3_hESC_Bernstein-UpDownRegulated_Epiblastpadj05fc025_heatmap1.pdf \
    --samplesLabel "H3K27me3" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 5


