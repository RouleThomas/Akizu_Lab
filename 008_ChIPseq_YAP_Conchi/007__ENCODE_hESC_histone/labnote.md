# Project and goal

Part of gastrulation paper with Conchi Estaras.

--> See `002*/003*/gastrulation paper/Figure 4` for detail: We need to know the overlap of the 199 QSER1:YAP cobound peaks (identified in `008*/001*` and `008*/003*`) with H3K27ac, H3K4me1, H3K27me3 and H3K36me3 modification marks

Pipeline:
- Identify data from ENCODE
- Download ENCODE data (bigwig and peak file)
- Follow `Figure 4` ppt

# Data download from ENCODE

Prioritize most recent datasets, or all datasets from same lab:
- **Bing Ren (GSE16256)**, UCSD_Project ROADMAP; [2013 2 Bio Rep; processed 2020](https://www.encodeproject.org/experiments/ENCSR928HYM/):
    - **H3K27me3**: *ENCFF395GVR* bigwig signal p-value; *ENCFF599KDF* bed
    - **H3K4me1**:  *ENCFF164XHJ* bigwig signal p-value; *ENCFF613QAB* bed
    - **H3K36me3**: *ENCFF483UMR* bigwig signal p-value; *ENCFF681CEO* bed
    - **H3K27ac**:  *ENCFF390JIZ* bigwig signal p-value; *ENCFF045CUG* bed
- **Bradley Bernstein (GSE29611)**, Broad_Project ENCODE; 2013
    - **H3K27me3**: *ENCFF193PKI* bigwig signal p-value; *ENCFF305KNA* bed
    - **H3K4me1**:  *ENCFF706CHK* bigwig signal p-value; *ENCFF984DGO* bed
    - **H3K36me3**: *ENCFF985CVI* bigwig signal p-value; *ENCFF504KOV* bed
    - **H3K27ac**:  *ENCFF771GNB* bigwig signal p-value; *ENCFF317QGQ* bed


--> File transfer to HPC cluster at `output/ENCODE`

--> Let's check both.


# Rename files

See in `008*/007*`: `sample_008007.xlsx` for file renaming


```bash
cd output/ENCODE

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename.txt
```

--> All good 




# Check overlap QSER1:YAP1 peaks


From the [webtool VennDiagram](https://www.bioinformatics.com.cn/plot_basic_genomic_regions_overlap_venn_diagram_026_en) output in `Gastrulation paper/output/peakOverlap/homer_peakCoordinate/noExtension` I put together the 199 (171+28) QSER1:YAP1 overlapping peaks in `QSER1YAP1_199peaks.xlsx` --> Copy in HPC cluster at `output/QSER1YAP1_199peaks.bed`


--> I used the webtool again to check for overlap I used  `output/QSER1YAP1_199peaks.bed` vs `*.bed` file from ENCODE

Let's extend of 200bp the QSER1YAP1_199peaks and check overlap again:

```bash
conda activate BedToBigwig

# extend peaks of 200bp up and downstream
bedtools slop -i output/QSER1YAP1_199peaks.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 >  output/QSER1YAP1_199peaks_extend200bp.bed
```

--> All good



Results in `Heatmap Figure 4_TR.pptx`

## Overlap using bedtools


Let's instead use bedtools, less sketchy than the app...


```bash
conda activate BedToBigwig


# QSER1:YAP1
## Direct overlap, non extension
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq | wc -l # 153
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq | wc -l # 96
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Ren_H3K27me3.bed | uniq | wc -l # 10
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Ren_H3K36me3.bed | uniq | wc -l # 1

bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq | wc -l # 115
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq | wc -l # 111
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Bernstein_H3K27me3.bed | uniq | wc -l # 11
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Bernstein_H3K36me3.bed | uniq | wc -l # 0

## 200bp extension
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq | wc -l # 167
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq | wc -l # 138
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Ren_H3K27me3.bed | uniq | wc -l # 13
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Ren_H3K36me3.bed | uniq | wc -l # 4

bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq | wc -l # 132
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq | wc -l # 159
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Bernstein_H3K27me3.bed | uniq | wc -l # 18
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Bernstein_H3K36me3.bed | uniq | wc -l # 0




# YAP1 only (3062 peaks)
## Generate YAP1 bed file
awk 'NR > 1' ../001__ChIPseq_V1/output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot.txt > output/annotation_homer_hESC_WT_YAP1_R1_annot.bed
cp ../001__ChIPseq_V1/output/ChIPseeker/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.txt output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed
## Direct overlap, non extension
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq | wc -l # 1811
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq | wc -l # 1374
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot.bed -b output/ENCODE/Ren_H3K27me3.bed | uniq | wc -l # 48
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot.bed -b output/ENCODE/Ren_H3K36me3.bed | uniq | wc -l # 54

bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq | wc -l # 1172
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq | wc -l # 1394
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot.bed -b output/ENCODE/Bernstein_H3K27me3.bed | uniq | wc -l # 62
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot.bed -b output/ENCODE/Bernstein_H3K36me3.bed | uniq | wc -l # 5

## 200bp extension
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq | wc -l # 2102
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq | wc -l # 1984
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed -b output/ENCODE/Ren_H3K27me3.bed | uniq | wc -l # 86
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed -b output/ENCODE/Ren_H3K36me3.bed | uniq | wc -l # 92

bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq | wc -l # 1350
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq | wc -l # 1900
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed -b output/ENCODE/Bernstein_H3K27me3.bed | uniq | wc -l # 95
bedtools intersect -wa -a output/annotation_homer_hESC_WT_YAP1_R1_annot_extend200bp.bed -b output/ENCODE/Bernstein_H3K36me3.bed | uniq | wc -l # 9






# QSER1 only (3062 peaks)
## Generate YAP1 bed file
awk 'NR > 1' ../001__ChIPseq_V1/output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot.txt > output/annotation_homer_hESC_WT_QSER1_pool_annot.bed # 12461
cp ../001__ChIPseq_V1/output/ChIPseeker/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.txt output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed
## Direct overlap, non extension
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq | wc -l # 4862
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq | wc -l # 3367
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Ren_H3K27me3.bed | uniq | wc -l # 1563
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Ren_H3K36me3.bed | uniq | wc -l # 58

bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq | wc -l # 3618
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq | wc -l # 4027
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Bernstein_H3K27me3.bed | uniq | wc -l # 2229
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Bernstein_H3K36me3.bed | uniq | wc -l # 14

## 200bp extension
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq | wc -l # 5746
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq | wc -l # 4726
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed -b output/ENCODE/Ren_H3K27me3.bed | uniq | wc -l # 2182
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed -b output/ENCODE/Ren_H3K36me3.bed | uniq | wc -l # 120

bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq | wc -l # 4230
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq | wc -l # 5607
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed -b output/ENCODE/Bernstein_H3K27me3.bed | uniq | wc -l # 2797
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot_extend200bp.bed -b output/ENCODE/Bernstein_H3K36me3.bed | uniq | wc -l # 33
```

--> ChIPseeker files from homer used for total peak counts (in `008*/001*` and `008*/003*`)

For emboR revision, we need the list of genes with **QSER1:YAP1 co-bound peaks (199) overlapping with H3K4me1 and H3K27ac**. Let's isolate bed file and then assign peak to genes (direct overlap AND +/-200bp)




```bash
conda activate BedToBigwig


# QSER1:YAP1
## Direct overlap, non extension
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq > output/ENCODE/QSER1YAP1_199peaks__Ren_H3K27ac.bed
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq > output/ENCODE/QSER1YAP1_199peaks__Ren_H3K4me1.bed

bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq > output/ENCODE/QSER1YAP1_199peaks__Bernstein_H3K27ac.bed
bedtools intersect -wa -a output/QSER1YAP1_199peaks.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq > output/ENCODE/QSER1YAP1_199peaks__Bernstein_H3K4me1.bed

## 200bp extension
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq > output/ENCODE/QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac.bed
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq > output/ENCODE/QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1.bed

bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq > output/ENCODE/QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac.bed
bedtools intersect -wa -a output/QSER1YAP1_199peaks_extend200bp.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq > output/ENCODE/QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1.bed
```

--> Good, next lets assign peak to genes (output in ppt: `002*/003*/gastrulation paper or QSER1 paper/Revision1 emboR/QSER1 paper_Revision_TR.pptx`)





# deeptool plots


Generate heatmap showing bigwig signal of the histone marks and EZH2 for (window of 2 / 5 / 10kb):
- YAP:QSER1 peaks
- YAP1 only peaks
- QSER1 only peaks
- Show EZH2 and H3K27me3 signal (Bernstein) in all genes and in EZH2-bound genes (from `008*/001*`) - *3D paper*


```bash
conda activate deeptools

# YAP:QSER1 peaks (H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2)
sbatch scripts/matrix_2kb_H3K4me3_QSER1YAP1peaks_Ren.sh # 29178677 ok
sbatch scripts/matrix_5kb_H3K4me3_QSER1YAP1peaks_Ren.sh # 29178682 ok
sbatch scripts/matrix_10kb_H3K4me3_QSER1YAP1peaks_Ren.sh # 29178686 ok

sbatch scripts/matrix_2kb_H3K4me3_QSER1YAP1peaks_Bernstein.sh # 29178826 ok
sbatch scripts/matrix_5kb_H3K4me3_QSER1YAP1peaks_Bernstein.sh # 29178852 ok
sbatch scripts/matrix_10kb_H3K4me3_QSER1YAP1peaks_Bernstein.sh  # 29178862 ok

# YAP:QSER1 peaks (YAP1, QSER1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2)
sbatch scripts/matrix_2kb_QSER1YAP1peaks_Ren.sh # 40532499 ok
sbatch scripts/matrix_5kb_QSER1YAP1peaks_Ren.sh # 40532511 ok

sbatch scripts/matrix_2kb_QSER1YAP1peaks_Bernstein.sh # 40532543 ok
sbatch scripts/matrix_5kb_QSER1YAP1peaks_Bernstein.sh # 40532556 ok

# YAP:QSER1 peaks (YAP1, QSER1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2, NIPBL)
sbatch scripts/matrix_2kb_QSER1YAP1peaks_Ren_2.sh #  interactive
sbatch scripts/matrix_2kb_QSER1YAP1peaks_Bernstein_2.sh #  interactive

# YAP:QSER1 peaks (YAP1, QSER1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2, NIPBL, NIBPLinput)
sbatch scripts/matrix_2kb_QSER1YAP1peaks_Ren_2input.sh #  interactive
sbatch scripts/matrix_2kb_QSER1YAP1peaks_Bernstein_2input.sh #  interactive
sbatch scripts/matrix_2kb_QSER1YAP1peaks_Ren_3input.sh #  interactive
sbatch scripts/matrix_2kb_QSER1YAP1peaks_Bernstein_3input.sh #  interactive


# 1538 YAP1 only peaks QSER1:YAP1 co bound genes (YAP1, QSER1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2, NIPBL, NIBPLinput)
sbatch scripts/matrix_2kb_H3K4me3_YAP1peaksQSER1YAP1genes_Ren.sh #  interactive
sbatch scripts/matrix_2kb_H3K4me3_YAP1peaksQSER1YAP1genes_Bernstein.sh # interactive

# 2073 QSER1 only peaks QSER1:YAP1 co bound genes (QSER1, YAP1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2, NIPBL, NIBPLinput)
sbatch scripts/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genes_Ren.sh # interactive
sbatch scripts/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genes_Bernstein.sh # interactive

# 2073 QSER1 only peaks QSER1:YAP1 co bound genes (YAP1 ,QSER1 , H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2, NIPBL, NIBPLinput)
sbatch scripts/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Ren.sh # interactive
sbatch scripts/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Bernstein.sh # interactive


# YAP1 only peaks
sbatch scripts/matrix_2kb_H3K4me3_YAP1peaks_Ren.sh # 29178940 ok
sbatch scripts/matrix_5kb_H3K4me3_YAP1peaks_Ren.sh # 29178972 ok
sbatch scripts/matrix_10kb_H3K4me3_YAP1peaks_Ren.sh # 29179019 ok

sbatch scripts/matrix_2kb_H3K4me3_YAP1peaks_Bernstein.sh # 29179075 ok
sbatch scripts/matrix_5kb_H3K4me3_YAP1peaks_Bernstein.sh # 29179088 ok
sbatch scripts/matrix_10kb_H3K4me3_YAP1peaks_Bernstein.sh # 29179106 ok

# QSER1 only peaks (H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2)
sbatch scripts/matrix_2kb_H3K4me3_QSER1peaks_Ren.sh # 29179264 ok
sbatch scripts/matrix_5kb_H3K4me3_QSER1peaks_Ren.sh # 29179285 ok
sbatch scripts/matrix_10kb_H3K4me3_QSER1peaks_Ren.sh # 29179334 ok

sbatch scripts/matrix_2kb_H3K4me3_QSER1peaks_Bernstein.sh # 29179355 ok
sbatch scripts/matrix_5kb_H3K4me3_QSER1peaks_Bernstein.sh # 29179370 ok
sbatch scripts/matrix_10kb_H3K4me3_QSER1peaks_Bernstein.sh # 29179387 ok

# QSER1 only peaks (QSER1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2)
sbatch scripts/matrix_2kb_QSER1peaks_Bernstein.sh # 40287751 ok
sbatch scripts/matrix_2kb_QSER1peaks_Ren.sh # 40292479 ok

# QSER1 only peaks (QSER1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2, NIPBL)
sbatch scripts/matrix_2kb_QSER1peaks_Bernstein_2.sh # interactive
sbatch scripts/matrix_2kb_QSER1peaks_Ren_2.sh # interactive

# QSER1 only peaks (QSER1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2, NIPBL, NIPBLinput)
sbatch scripts/matrix_2kb_QSER1peaks_Bernstein_2input.sh #  interactive
sbatch scripts/matrix_2kb_QSER1peaks_Ren_2input.sh # interactive

# EZH2 and H3K27me3 signal - All genes
sbatch scripts/matrix_5kb_H3K27me3_Bernstein_EZH2008001_allGenes.sh # 31762805 ok
sbatch scripts/matrix_10kb_H3K27me3_Bernstein_EZH2008001_allGenes.sh # 31762991 ok

# EZH2 and H3K27me3 signal - EZH2 target genes (008*/001*)
sbatch scripts/matrix_5kb_H3K27me3_Bernstein_EZH2008001_EZH2target.sh # 31763984 ok
sbatch scripts/matrix_10kb_H3K27me3_Bernstein_EZH2008001_EZH2target.sh # 31764129 ok
```

- *NOTE: naming `*_H3K4me3_*` make no sense here...; typo...*
- *NOTE: EZH2 target gene list GTF generated below in `# Heatmap H3K27me3 signal in DEG Epiblast - *3D paper*`*


--> EZH2 (`008*/001*`) and H3K27me3 (*ENCODE_Bernstein*) well correlate



# Overlap QSER1 peaks with H3K4me1 and H3K27ac - emboR revision tasks, email 4/5/2023

emboR revision tasks, email 4/5/2023:

- For 3F: How many of the QSER1 peaks are included in the heatmap? (number)
12,461 peaks (identified with homer) (`wc -l output/annotation_homer_hESC_WT_QSER1_pool_annot.bed`; file used to generate heatmaps)

- For 3F: How many of the QSER1 peaks bind H3K4me1 and H3K27ac?(number)
```bash
conda activate BedToBigwig

bedtools intersect -v -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Bernstein_H3K4me1.bed | wc -l # 8434
bedtools intersect -v -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Bernstein_H3K27ac.bed | wc -l # 8843

bedtools intersect -v -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Ren_H3K4me1.bed | wc -l # 9094
bedtools intersect -v -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Ren_H3K27ac.bed | wc -l # 7599
```


- NIPBL overlap: how many of the actual YAP and QSER1 peaks overlap with NIPBL peaks on the 1192 genes?

On the 1192 genes co-bound by QSER1 and YAP1. Here is the respective overlap for QSER1 and YAP1 with NIPBL:
(1192 is the QSER1:YAP1 cobound genes.)

Peak coordinate of QSER1 and YAP1 peaks from the 1192 QSER1:YAP1 co-bound genes as already been isolated in `008*/001*`:
- Including intergenic peaks:
    - QSER1 peak= `../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker.txt` --> convert to bed and remove header: `../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed`
    - YAP1 peak= `../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1.txt` --> convert to bed and remove header: `../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed`


```bash
conda activate BedToBigwig
# Check overlap QSER1 peaks with NIPBL
bedtools intersect -v -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader.bed | wc -l # 1548
wc -l ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed # 2073

# Check overlap YAP1 peaks with NIPBL
bedtools intersect -v -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader.bed | wc -l # 639
wc -l ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed # 1538
```

--> Among the 1192 genes cobound by QSER1 and YAP1; I exported the QSER1 and YAP1 peaks and check their overlap with NIPBL:
- Among the 1192 genes, there are 2073 QSER1 peaks, from which 1548 directly overlap with NIPBL
- Among the 1192 genes, there are 1538 YAP1 peaks, from which 639 directly overlap with NIPBL


## HUGE MISTAKE: email 20250604; I use -v instead of -wa to chekc for overlap... HERE below correction:


```bash
conda activate BedToBigwig
# Check overlap QSER1 peaks with NIPBL
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader.bed | uniq | wc -l # 525
wc -l ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed # 2073

# Check overlap YAP1 peaks with NIPBL
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader.bed | uniq | wc -l # 899
wc -l ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed # 1538
```

--> Among the 1192 genes cobound by QSER1 and YAP1; I exported the QSER1 and YAP1 peaks and check their overlap with NIPBL:
- Among the 1192 genes, there are 2073 QSER1 peaks, from which 525 directly overlap with NIPBL
- Among the 1192 genes, there are 1538 YAP1 peaks, from which 899 directly overlap with NIPBL


If we add extension to these QSER1, YAP1, NIPBL peaks:



```bash
# Add extension 
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend200bp.bed
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend500bp.bed
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend1kp.bed
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend2kp.bed
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 5000 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend5kp.bed

bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend200bp.bed
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend500bp.bed
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend1kp.bed
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend2kp.bed
bedtools slop -i ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 5000 > ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend5kp.bed



# extended files:
../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend200bp.bed
../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend500bp.bed
../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend1kp.bed
../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend2kp.bed
../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend5kp.bed
../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend200bp.bed
../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend500bp.bed
../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend1kp.bed
../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend2kp.bed
../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend200bp.bed
../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend500bp.bed
../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend1kp.bed
../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend2kp.bed




# overlap
## Check overlap QSER1 peaks with NIPBL
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend200bp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend200bp.bed | uniq | wc -l # 653
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend500bp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend500bp.bed | uniq | wc -l # 750
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend1kp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend1kp.bed | uniq | wc -l # 803
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend2kp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend2kp.bed | uniq | wc -l # 891
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader_extend5kp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend5kp.bed | uniq | wc -l # 1090

## Check overlap YAP1 peaks with NIPBL
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend200bp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend200bp.bed | uniq | wc -l # 923
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend500bp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend500bp.bed | uniq | wc -l # 941
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend1kp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend1kp.bed | uniq | wc -l # 965
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend2kp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend2kp.bed | uniq | wc -l # 1006
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_YAP1_noHeader_extend5kp.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader_extend5kp.bed | uniq | wc -l # 1082


```

--> Among the 1192 genes cobound by QSER1 and YAP1; I exported the QSER1 and YAP1 peaks and check their overlap with NIPBL:
- Among the 1192 genes, there are 2073 QSER1 peaks, from which 525 directly overlap with NIPBL (extension 200bp-> 653 peak overlap; 500bp-> 750; 1kb-> 803; 2kb-> 891; 5kb-> 1090)
- Among the 1192 genes, there are 1538 YAP1 peaks, from which 899 directly overlap with NIPBL (extension 200bp-> 923 peak overlap; 500bp-> 941; 1kb-> 965; 2kb-> 1006; 5kb-> 1082)






Same here huge mistake


- For 3F: How many of the QSER1 peaks bind H3K4me1 and H3K27ac?(number)
```bash
conda activate BedToBigwig

bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Bernstein_H3K4me1.bed | uniq | wc -l # 4027
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Bernstein_H3K27ac.bed | uniq | wc -l # 3618

bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Ren_H3K4me1.bed | uniq | wc -l # 3367
bedtools intersect -wa -a output/annotation_homer_hESC_WT_QSER1_pool_annot.bed -b output/ENCODE/Ren_H3K27ac.bed | uniq | wc -l # 4862



# Add extension 
bedtools slop -i output/ENCODE/Bernstein_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > output/ENCODE/Bernstein_H3K4me1_extend200bp.bed
bedtools slop -i output/ENCODE/Bernstein_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > output/ENCODE/Bernstein_H3K4me1_extend500bp.bed
bedtools slop -i output/ENCODE/Bernstein_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > output/ENCODE/Bernstein_H3K4me1_extend1kp.bed
bedtools slop -i output/ENCODE/Bernstein_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > output/ENCODE/Bernstein_H3K4me1_extend2kp.bed
bedtools slop -i output/ENCODE/Bernstein_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 5000 > output/ENCODE/Bernstein_H3K4me1_extend5kp.bed


bedtools slop -i output/ENCODE/Bernstein_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > output/ENCODE/Bernstein_H3K27ac_extend200bp.bed
bedtools slop -i output/ENCODE/Bernstein_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > output/ENCODE/Bernstein_H3K27ac_extend500bp.bed
bedtools slop -i output/ENCODE/Bernstein_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > output/ENCODE/Bernstein_H3K27ac_extend1kp.bed
bedtools slop -i output/ENCODE/Bernstein_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > output/ENCODE/Bernstein_H3K27ac_extend2kp.bed
bedtools slop -i output/ENCODE/Bernstein_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 5000 > output/ENCODE/Bernstein_H3K27ac_extend5kp.bed

bedtools slop -i output/ENCODE/Ren_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > output/ENCODE/Ren_H3K4me1_extend200bp.bed
bedtools slop -i output/ENCODE/Ren_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > output/ENCODE/Ren_H3K4me1_extend500bp.bed
bedtools slop -i output/ENCODE/Ren_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > output/ENCODE/Ren_H3K4me1_extend1kp.bed
bedtools slop -i output/ENCODE/Ren_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > output/ENCODE/Ren_H3K4me1_extend2kp.bed
bedtools slop -i output/ENCODE/Ren_H3K4me1.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 5000 > output/ENCODE/Ren_H3K4me1_extend5kp.bed



bedtools slop -i output/ENCODE/Ren_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 200 > output/ENCODE/Ren_H3K27ac_extend200bp.bed
bedtools slop -i output/ENCODE/Ren_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 500 > output/ENCODE/Ren_H3K27ac_extend500bp.bed
bedtools slop -i output/ENCODE/Ren_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 1000 > output/ENCODE/Ren_H3K27ac_extend1kp.bed
bedtools slop -i output/ENCODE/Ren_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 2000 > output/ENCODE/Ren_H3K27ac_extend2kp.bed
bedtools slop -i output/ENCODE/Ren_H3K27ac.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab -b 5000 > output/ENCODE/Ren_H3K27ac_extend5kp.bed




## Files extended
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend1kp.bed
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend2kp.bed
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend5kp.bed

output/ENCODE/Bernstein_H3K4me1_extend200bp.bed
output/ENCODE/Bernstein_H3K4me1_extend500bp.bed
output/ENCODE/Bernstein_H3K4me1_extend1kp.bed
output/ENCODE/Bernstein_H3K4me1_extend2kp.bed
output/ENCODE/Bernstein_H3K4me1_extend5kp.bed

output/ENCODE/Bernstein_H3K27ac_extend200bp.bed
output/ENCODE/Bernstein_H3K27ac_extend500bp.bed
output/ENCODE/Bernstein_H3K27ac_extend1kp.bed
output/ENCODE/Bernstein_H3K27ac_extend2kp.bed
output/ENCODE/Bernstein_H3K27ac_extend5kp.bed

Ren_H3K4me1_extend200bp.bed
Ren_H3K4me1_extend500bp.bed
Ren_H3K4me1_extend1kp.bed
Ren_H3K4me1_extend2kp.bed
Ren_H3K4me1_extend5kp.bed

output/ENCODE/Ren_H3K27ac_extend200bp.bed
output/ENCODE/Ren_H3K27ac_extend500bp.bed
output/ENCODE/Ren_H3K27ac_extend1kp.bed
output/ENCODE/Ren_H3K27ac_extend2kp.bed
output/ENCODE/Ren_H3K27ac_extend5kp.bed


# overlap extended
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed -b output/ENCODE/Bernstein_H3K4me1_extend200bp.bed | uniq | wc -l # 6710
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed -b output/ENCODE/Bernstein_H3K4me1_extend500bp.bed | uniq | wc -l # 8522
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend1kp.bed -b output/ENCODE/Bernstein_H3K4me1_extend1kp.bed | uniq | wc -l # 9478
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend2kp.bed -b output/ENCODE/Bernstein_H3K4me1_extend2kp.bed | uniq | wc -l # 9974
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend5kp.bed -b output/ENCODE/Bernstein_H3K4me1_extend5kp.bed | uniq | wc -l # 10721

bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed -b output/ENCODE/Bernstein_H3K27ac_extend200bp.bed | uniq | wc -l # 4420
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed -b output/ENCODE/Bernstein_H3K27ac_extend500bp.bed | uniq | wc -l # 4652
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend1kp.bed -b output/ENCODE/Bernstein_H3K27ac_extend1kp.bed | uniq | wc -l # 4820
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend2kp.bed -b output/ENCODE/Bernstein_H3K27ac_extend2kp.bed | uniq | wc -l # 5053
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend5kp.bed -b output/ENCODE/Bernstein_H3K27ac_extend5kp.bed | uniq | wc -l # 


bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed -b output/ENCODE/Ren_H3K4me1_extend200bp.bed | uniq | wc -l # 5743
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed -b output/ENCODE/Ren_H3K4me1_extend500bp.bed | uniq | wc -l # 7625
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend1kp.bed -b output/ENCODE/Ren_H3K4me1_extend1kp.bed | uniq | wc -l # 8751
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend2kp.bed -b output/ENCODE/Ren_H3K4me1_extend2kp.bed | uniq | wc -l # 9356
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend5kp.bed -b output/ENCODE/Ren_H3K4me1_extend5kp.bed | uniq | wc -l # 10343


bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend200bp.bed -b output/ENCODE/Ren_H3K27ac_extend200bp.bed | uniq | wc -l # 6005
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend500bp.bed -b output/ENCODE/Ren_H3K27ac_extend500bp.bed | uniq | wc -l # 6275
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend1kp.bed -b output/ENCODE/Ren_H3K27ac_extend1kp.bed | uniq | wc -l # 6494
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend2kp.bed -b output/ENCODE/Ren_H3K27ac_extend2kp.bed | uniq | wc -l # 6758
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks_extend5kp.bed -b output/ENCODE/Ren_H3K27ac_extend5kp.bed | uniq | wc -l # 7589


```


Let's now check the overlap for the heatmap with the QSER1 peaks with QSER1, H3K4me1, H3K27ac, H3K36me3, H3K27me3, EZH2, NIPBL

```bash
conda activate BedToBigwig

# files
../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed
../007__ENCODE_hESC_histone/output/ENCODE/Bernstein_H3K4me1.bed
../007__ENCODE_hESC_histone/output/ENCODE/Bernstein_H3K27ac.bed
../007__ENCODE_hESC_histone/output/ENCODE/Bernstein_H3K36me3.bed
../007__ENCODE_hESC_histone/output/ENCODE/Bernstein_H3K27me3.bed
../001__ChIPseq_V1/output/homer/hESC_WT_EZH2_outputPeaks.bed
../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader.bed
 


bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../007__ENCODE_hESC_histone/output/ENCODE/Bernstein_H3K4me1.bed | uniq | wc -l # 4027
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../007__ENCODE_hESC_histone/output/ENCODE/Bernstein_H3K27ac.bed | uniq | wc -l # 3618
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../007__ENCODE_hESC_histone/output/ENCODE/Bernstein_H3K36me3.bed | uniq | wc -l # 14
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../007__ENCODE_hESC_histone/output/ENCODE/Bernstein_H3K27me3.bed | uniq | wc -l # 2229
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../001__ChIPseq_V1/output/homer/hESC_WT_EZH2_outputPeaks.bed | uniq | wc -l # 445
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../008__ChIPseq_Nipbl/output/homer/Nipbl/peaks_noHeader.bed | uniq | wc -l # 2956

bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../007__ENCODE_hESC_histone/output/ENCODE/Ren_H3K4me1.bed | uniq | wc -l # 3367
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../007__ENCODE_hESC_histone/output/ENCODE/Ren_H3K27ac.bed | uniq | wc -l # 4862
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../007__ENCODE_hESC_histone/output/ENCODE/Ren_H3K36me3.bed | uniq | wc -l # 58
bedtools intersect -wa -a ../001__ChIPseq_V1/output/homer/hESC_WT_QSER1_outputPeaks.bed -b ../007__ENCODE_hESC_histone/output/ENCODE/Ren_H3K27me3.bed | uniq | wc -l # 1563


```








# ChIPseeker for histone ENCODE files





```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")

# SIMPLICATE SAMPLE #####################

# Import homer peaks
# Add header to bed
Ren_H3K27ac = as_tibble(read.table("output/ENCODE/Ren_H3K27ac.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
Ren_H3K4me1 = as_tibble(read.table("output/ENCODE/Ren_H3K4me1.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
Ren_H3K27me3 = as_tibble(read.table("output/ENCODE/Ren_H3K27me3.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
Ren_H3K36me3 = as_tibble(read.table("output/ENCODE/Ren_H3K36me3.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
Bernstein_H3K27ac = as_tibble(read.table("output/ENCODE/Bernstein_H3K27ac.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
Bernstein_H3K4me1 = as_tibble(read.table("output/ENCODE/Bernstein_H3K4me1.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
Bernstein_H3K27me3 = as_tibble(read.table("output/ENCODE/Bernstein_H3K27me3.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
Bernstein_H3K36me3 = as_tibble(read.table("output/ENCODE/Bernstein_H3K36me3.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 



## Tidy peaks 
Ren_H3K27ac_gr = makeGRangesFromDataFrame(Ren_H3K27ac,keep.extra.columns=TRUE)
Ren_H3K4me1_gr = makeGRangesFromDataFrame(Ren_H3K4me1,keep.extra.columns=TRUE)
Ren_H3K27me3_gr = makeGRangesFromDataFrame(Ren_H3K27me3,keep.extra.columns=TRUE)
Ren_H3K36me3_gr = makeGRangesFromDataFrame(Ren_H3K36me3,keep.extra.columns=TRUE)
Bernstein_H3K27ac_gr = makeGRangesFromDataFrame(Bernstein_H3K27ac,keep.extra.columns=TRUE)
Bernstein_H3K4me1_gr = makeGRangesFromDataFrame(Bernstein_H3K4me1,keep.extra.columns=TRUE)
Bernstein_H3K27me3_gr = makeGRangesFromDataFrame(Bernstein_H3K27me3,keep.extra.columns=TRUE)
Bernstein_H3K36me3_gr = makeGRangesFromDataFrame(Bernstein_H3K36me3,keep.extra.columns=TRUE)



gr_list <- list(Ren_H3K27ac=Ren_H3K27ac_gr, Ren_H3K4me1=Ren_H3K4me1_gr, Ren_H3K27me3=Ren_H3K27me3_gr,    Ren_H3K36me3 = Ren_H3K36me3_gr, Bernstein_H3K27ac = Bernstein_H3K27ac_gr, Bernstein_H3K4me1 = Bernstein_H3K4me1_gr, Bernstein_H3K27me3 = Bernstein_H3K27me3_gr, Bernstein_H3K36me3 = Bernstein_H3K36me3_gr)

## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here

                       
### Barplot
pdf("output/ChIPseeker/annotation_barplot_ENCODE.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()

## Get annotation data frame
Ren_H3K27ac_annot <- as.data.frame(peakAnnoList[["Ren_H3K27ac"]]@anno)
Ren_H3K4me1_annot <- as.data.frame(peakAnnoList[["Ren_H3K4me1"]]@anno)
Ren_H3K27me3_annot <- as.data.frame(peakAnnoList[["Ren_H3K27me3"]]@anno)
Ren_H3K36me3_annot <- as.data.frame(peakAnnoList[["Ren_H3K36me3"]]@anno)
Bernstein_H3K27ac_annot <- as.data.frame(peakAnnoList[["Bernstein_H3K27ac"]]@anno)
Bernstein_H3K4me1_annot <- as.data.frame(peakAnnoList[["Bernstein_H3K4me1"]]@anno)
Bernstein_H3K27me3_annot <- as.data.frame(peakAnnoList[["Bernstein_H3K27me3"]]@anno)
Bernstein_H3K36me3_annot <- as.data.frame(peakAnnoList[["Bernstein_H3K36me3"]]@anno)


## Convert entrez gene IDs to gene symbols
Ren_H3K27ac_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Ren_H3K27ac_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Ren_H3K27ac_annot$gene <- mapIds(org.Hs.eg.db, keys = Ren_H3K27ac_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
Ren_H3K4me1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Ren_H3K4me1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Ren_H3K4me1_annot$gene <- mapIds(org.Hs.eg.db, keys = Ren_H3K4me1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
Ren_H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Ren_H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Ren_H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = Ren_H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
Ren_H3K36me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Ren_H3K36me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Ren_H3K36me3_annot$gene <- mapIds(org.Hs.eg.db, keys = Ren_H3K36me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
Bernstein_H3K27ac_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Bernstein_H3K27ac_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Bernstein_H3K27ac_annot$gene <- mapIds(org.Hs.eg.db, keys = Bernstein_H3K27ac_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
Bernstein_H3K4me1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Bernstein_H3K4me1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Bernstein_H3K4me1_annot$gene <- mapIds(org.Hs.eg.db, keys = Bernstein_H3K4me1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
Bernstein_H3K27me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Bernstein_H3K27me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Bernstein_H3K27me3_annot$gene <- mapIds(org.Hs.eg.db, keys = Bernstein_H3K27me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
Bernstein_H3K36me3_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = Bernstein_H3K36me3_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
Bernstein_H3K36me3_annot$gene <- mapIds(org.Hs.eg.db, keys = Bernstein_H3K36me3_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



## Save output table
write.table(Ren_H3K27ac_annot, file="output/ChIPseeker/annotation_ENCODE_Ren_H3K27ac_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(Ren_H3K4me1_annot, file="output/ChIPseeker/annotation_ENCODE_Ren_H3K4me1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(Ren_H3K27me3_annot, file="output/ChIPseeker/annotation_ENCODE_Ren_H3K27me3_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(Ren_H3K36me3_annot, file="output/ChIPseeker/annotation_ENCODE_Ren_H3K36me3_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(Bernstein_H3K27ac_annot, file="output/ChIPseeker/annotation_ENCODE_Bernstein_H3K27ac_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(Bernstein_H3K4me1_annot, file="output/ChIPseeker/annotation_ENCODE_Bernstein_H3K4me1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(Bernstein_H3K27me3_annot, file="output/ChIPseeker/annotation_ENCODE_Bernstein_H3K27me3_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(Bernstein_H3K36me3_annot, file="output/ChIPseeker/annotation_ENCODE_Bernstein_H3K36me3_annot.txt", sep="\t", quote=F, row.names=F) 



## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
Ren_H3K27ac_annot_noIntergenic = tibble(Ren_H3K27ac_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
Ren_H3K4me1_annot_noIntergenic = tibble(Ren_H3K4me1_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
Ren_H3K27me3_annot_noIntergenic = tibble(Ren_H3K27me3_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
Ren_H3K36me3_annot_noIntergenic = tibble(Ren_H3K36me3_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
Bernstein_H3K27ac_annot_noIntergenic = tibble(Bernstein_H3K27ac_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
Bernstein_H3K4me1_annot_noIntergenic = tibble(Bernstein_H3K4me1_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
Bernstein_H3K27me3_annot_noIntergenic = tibble(Bernstein_H3K27me3_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
Bernstein_H3K36me3_annot_noIntergenic = tibble(Bernstein_H3K36me3_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()




write.table(Ren_H3K27ac_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_Ren_H3K27ac_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(Ren_H3K4me1_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_Ren_H3K4me1_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(Ren_H3K27me3_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_Ren_H3K27me3_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(Ren_H3K36me3_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_Ren_H3K36me3_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(Bernstein_H3K27ac_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_Bernstein_H3K27ac_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(Bernstein_H3K4me1_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_Bernstein_H3K4me1_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(Bernstein_H3K27me3_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_Bernstein_H3K27me3_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(Bernstein_H3K36me3_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_Bernstein_H3K36me3_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


```




# ChIPseeker for QSER1:YAP1 overlapping with H3K4me1, H3K27ac - emboR revision

- no extension
- QSER1:YAP1 peak extended of 200bp prior checking overlap with histone H3K4me1, H3K27ac



```bash
conda activate deseq2
```

```R
library("ChIPseeker")
library("tidyverse")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # hg 38 annot v41
library("clusterProfiler")
library("meshes")
library("ReactomePA")
library("org.Hs.eg.db")
library("VennDiagram")


# Import 
# Add header to bed
QSER1YAP1_199peaks__Ren_H3K27ac = as_tibble(read.table("output/ENCODE/QSER1YAP1_199peaks__Ren_H3K27ac.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
QSER1YAP1_199peaks__Ren_H3K4me1 = as_tibble(read.table("output/ENCODE/QSER1YAP1_199peaks__Ren_H3K4me1.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
QSER1YAP1_199peaks__Bernstein_H3K27ac = as_tibble(read.table("output/ENCODE/QSER1YAP1_199peaks__Bernstein_H3K27ac.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
QSER1YAP1_199peaks__Bernstein_H3K4me1 = as_tibble(read.table("output/ENCODE/QSER1YAP1_199peaks__Bernstein_H3K4me1.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac = as_tibble(read.table("output/ENCODE/QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1 = as_tibble(read.table("output/ENCODE/QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac = as_tibble(read.table("output/ENCODE/QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1 = as_tibble(read.table("output/ENCODE/QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1.bed")) %>%
    dplyr::rename(Chr=V1, start=V2, end=V3, name=V4) 



## Tidy peaks 
QSER1YAP1_199peaks__Ren_H3K27ac_gr = makeGRangesFromDataFrame(QSER1YAP1_199peaks__Ren_H3K27ac,keep.extra.columns=TRUE)
QSER1YAP1_199peaks__Ren_H3K4me1_gr = makeGRangesFromDataFrame(QSER1YAP1_199peaks__Ren_H3K4me1,keep.extra.columns=TRUE)
QSER1YAP1_199peaks__Bernstein_H3K27ac_gr = makeGRangesFromDataFrame(QSER1YAP1_199peaks__Bernstein_H3K27ac,keep.extra.columns=TRUE)
QSER1YAP1_199peaks__Bernstein_H3K4me1_gr = makeGRangesFromDataFrame(QSER1YAP1_199peaks__Bernstein_H3K4me1,keep.extra.columns=TRUE)
QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_gr = makeGRangesFromDataFrame(QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac,keep.extra.columns=TRUE)
QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_gr = makeGRangesFromDataFrame(QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1,keep.extra.columns=TRUE)
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_gr = makeGRangesFromDataFrame(QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac,keep.extra.columns=TRUE)
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_gr = makeGRangesFromDataFrame(QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1,keep.extra.columns=TRUE)



gr_list <- list(QSER1YAP1_199peaks__Ren_H3K27ac=QSER1YAP1_199peaks__Ren_H3K27ac_gr, QSER1YAP1_199peaks__Ren_H3K4me1=QSER1YAP1_199peaks__Ren_H3K4me1_gr, QSER1YAP1_199peaks__Bernstein_H3K27ac=QSER1YAP1_199peaks__Bernstein_H3K27ac_gr,    QSER1YAP1_199peaks__Bernstein_H3K4me1 = QSER1YAP1_199peaks__Bernstein_H3K4me1_gr, QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac = QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_gr, QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1 = QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_gr, QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac = QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_gr, QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1 = QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_gr)

## Export Gene peak assignemnt
peakAnnoList <- lapply(gr_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) # Not sure defeining the tssRegion is used here

                       
### Barplot
pdf("output/ChIPseeker/annotation_barplot_QSER1YAP1_199peaks_ENCODE_BernsteinRen.pdf", width=14, height=5)
plotAnnoBar(peakAnnoList)
dev.off()

## Get annotation data frame
QSER1YAP1_199peaks__Ren_H3K27ac_annot <- as.data.frame(peakAnnoList[["QSER1YAP1_199peaks__Ren_H3K27ac"]]@anno)
QSER1YAP1_199peaks__Ren_H3K4me1_annot <- as.data.frame(peakAnnoList[["QSER1YAP1_199peaks__Ren_H3K4me1"]]@anno)
QSER1YAP1_199peaks__Bernstein_H3K27ac_annot <- as.data.frame(peakAnnoList[["QSER1YAP1_199peaks__Bernstein_H3K27ac"]]@anno)
QSER1YAP1_199peaks__Bernstein_H3K4me1_annot <- as.data.frame(peakAnnoList[["QSER1YAP1_199peaks__Bernstein_H3K4me1"]]@anno)
QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot <- as.data.frame(peakAnnoList[["QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac"]]@anno)
QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot <- as.data.frame(peakAnnoList[["QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1"]]@anno)
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot <- as.data.frame(peakAnnoList[["QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac"]]@anno)
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot <- as.data.frame(peakAnnoList[["QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1"]]@anno)


## Convert entrez gene IDs to gene symbols
QSER1YAP1_199peaks__Ren_H3K27ac_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks__Ren_H3K27ac_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
QSER1YAP1_199peaks__Ren_H3K27ac_annot$gene <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks__Ren_H3K27ac_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
QSER1YAP1_199peaks__Ren_H3K4me1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks__Ren_H3K4me1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
QSER1YAP1_199peaks__Ren_H3K4me1_annot$gene <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks__Ren_H3K4me1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
QSER1YAP1_199peaks__Bernstein_H3K27ac_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks__Bernstein_H3K27ac_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
QSER1YAP1_199peaks__Bernstein_H3K27ac_annot$gene <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks__Bernstein_H3K27ac_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
QSER1YAP1_199peaks__Bernstein_H3K4me1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks__Bernstein_H3K4me1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
QSER1YAP1_199peaks__Bernstein_H3K4me1_annot$gene <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks__Bernstein_H3K4me1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot$gene <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot$gene <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot$gene <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot$geneSymbol <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot$geneId, column = "SYMBOL", keytype = "ENTREZID")
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot$gene <- mapIds(org.Hs.eg.db, keys = QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot$geneId, column = "ENSEMBL", keytype = "ENTREZID")



## Save output table
write.table(QSER1YAP1_199peaks__Ren_H3K27ac_annot, file="output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks__Ren_H3K27ac_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(QSER1YAP1_199peaks__Ren_H3K4me1_annot, file="output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks__Ren_H3K4me1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(QSER1YAP1_199peaks__Bernstein_H3K27ac_annot, file="output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks__Bernstein_H3K27ac_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(QSER1YAP1_199peaks__Bernstein_H3K4me1_annot, file="output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks__Bernstein_H3K4me1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot, file="output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot, file="output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot, file="output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot.txt", sep="\t", quote=F, row.names=F) 
write.table(QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot, file="output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot.txt", sep="\t", quote=F, row.names=F) 




## Keep only signals in non intergenic region ############################################# TO CHANGE IF NEEDED !!!!!!!!!!!!!!!!!!!
QSER1YAP1_199peaks__Ren_H3K27ac_annot_noIntergenic = tibble(QSER1YAP1_199peaks__Ren_H3K27ac_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
QSER1YAP1_199peaks__Ren_H3K4me1_annot_noIntergenic = tibble(QSER1YAP1_199peaks__Ren_H3K4me1_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
QSER1YAP1_199peaks__Bernstein_H3K27ac_annot_noIntergenic = tibble(QSER1YAP1_199peaks__Bernstein_H3K27ac_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
QSER1YAP1_199peaks__Bernstein_H3K4me1_annot_noIntergenic = tibble(QSER1YAP1_199peaks__Bernstein_H3K4me1_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot_noIntergenic = tibble(QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot_noIntergenic = tibble(QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot_noIntergenic = tibble(QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()
QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot_noIntergenic = tibble(QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot) %>%
    filter(annotation != c("Distal Intergenic"))%>%
    dplyr::select(geneSymbol) %>%
    unique()




write.table(QSER1YAP1_199peaks__Ren_H3K27ac_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks__Ren_H3K27ac_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(QSER1YAP1_199peaks__Ren_H3K4me1_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks__Ren_H3K4me1_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(QSER1YAP1_199peaks__Bernstein_H3K27ac_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks__Bernstein_H3K27ac_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(QSER1YAP1_199peaks__Bernstein_H3K4me1_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks__Bernstein_H3K4me1_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks_extend200bp__Ren_H3K27ac_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks_extend200bp__Ren_H3K4me1_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks_extend200bp__Bernstein_H3K27ac_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)
write.table(QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot_noIntergenic, file = "output/ChIPseeker/annotation_ENCODE_QSER1YAP1_199peaks_extend200bp__Bernstein_H3K4me1_annot_noIntergenic_geneSymbol.txt",
            quote = FALSE, 
            sep = "\t", 
            col.names = FALSE, 
            row.names = FALSE)


```





# Heatmap H3K27me3 signal in DEG Epiblast - *3D paper*


Let's generate heatmap of H3K27me3 signal (Bernstein) for DEGs (for H3K27me3- target genes, or not) in Epiblast (`002*/003*`):
- Generate gene list manually with `nano` from `DASATINIB2472hrs_response_dim30kparam15res04_allGenes.xlsx` (see file names below) + VennDiagram
- Generate gtf file of gene list
- Generate heatmap


```bash
# Generate gtf file from gene list:
../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_upregulatedq05fc025.txt # DEG up
../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_downregulatedq05fc025.txt # DEG down
../../008_ChIPseq_YAP_Conchi/007__ENCODE_hESC_histone/output/ChIPseeker/H3K27me3_hESC_Bernstein-UpRegulated_Epiblastpadj05fc025.txt # DEG up H3K27me3 target
../../008_ChIPseq_YAP_Conchi/007__ENCODE_hESC_histone/output/ChIPseeker/H3K27me3_hESC_Bernstein-DownRegulated_Epiblastpadj05fc025.txt # DEG up H3K27me3 target
../../008_ChIPseq_YAP_Conchi/001__ChIPseq_V1/output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_geneSymbol.txt # EZH2 target (from 008*/001*)



### create gtf from gene list
#### Modify the .txt file that list all genes so that it match gtf structure
sed 's/\r$//; s/.*/gene_name "&"/' ../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_upregulatedq05fc025.txt > ../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_upregulatedq05fc025_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_downregulatedq05fc025.txt > ../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_downregulatedq05fc025_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../../008_ChIPseq_YAP_Conchi/007__ENCODE_hESC_histone/output/ChIPseeker/H3K27me3_hESC_Bernstein-UpRegulated_Epiblastpadj05fc025.txt > ../../008_ChIPseq_YAP_Conchi/007__ENCODE_hESC_histone/output/ChIPseeker/H3K27me3_hESC_Bernstein-UpRegulated_Epiblastpadj05fc025_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../../008_ChIPseq_YAP_Conchi/007__ENCODE_hESC_histone/output/ChIPseeker/H3K27me3_hESC_Bernstein-DownRegulated_Epiblastpadj05fc025.txt > ../../008_ChIPseq_YAP_Conchi/007__ENCODE_hESC_histone/output/ChIPseeker/H3K27me3_hESC_Bernstein-DownRegulated_Epiblastpadj05fc025_as_gtf_geneSymbol.txt
sed 's/\r$//; s/.*/gene_name "&"/' ../../008_ChIPseq_YAP_Conchi/001__ChIPseq_V1/output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_geneSymbol.txt > ../../008_ChIPseq_YAP_Conchi/001__ChIPseq_V1/output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_as_gtf_geneSymbol.txt


## Filter the gtf
grep -Ff ../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_upregulatedq05fc025_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_upregulatedq05fc025.gtf
grep -Ff ../../002_scRNAseq/003__YAP1/output/seurat/Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_downregulatedq05fc025_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_downregulatedq05fc025.gtf
grep -Ff ../../008_ChIPseq_YAP_Conchi/007__ENCODE_hESC_histone/output/ChIPseeker/H3K27me3_hESC_Bernstein-UpRegulated_Epiblastpadj05fc025_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_H3K27me3_hESC_Bernstein-UpRegulated_Epiblastpadj05fc025.gtf
grep -Ff ../../008_ChIPseq_YAP_Conchi/007__ENCODE_hESC_histone/output/ChIPseeker/H3K27me3_hESC_Bernstein-DownRegulated_Epiblastpadj05fc025_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_H3K27me3_hESC_Bernstein-DownRegulated_Epiblastpadj05fc025.gtf
grep -Ff ../../008_ChIPseq_YAP_Conchi/001__ChIPseq_V1/output/ChIPseeker/annotation_homer_hESC_WT_EZH2_pool_annot_as_gtf_geneSymbol.txt meta/ENCFF159KBI.gtf > meta/ENCFF159KBI_homer_hESC_WT_EZH2_pool.gtf



```

Deeptools heatmap plots:

```bash
conda activate deeptools

# DEGs
sbatch scripts/matrix_5kb_H3K27me3_Bernstein-Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_upDownRegulatedq05fc025.sh # 31742879 ok
sbatch scripts/matrix_10kb_H3K27me3_Bernstein-Epiblast-DASATINIB2472hrs_response_dim30kparam15res04_upDownRegulatedq05fc025.sh # 31745396 ok

# DEGs and H3K27me3 target
sbatch scripts/matrix_5kb_H3K27me3_Bernstein-H3K27me3_hESC_Bernstein-UpDownRegulated_Epiblastpadj05fc025.sh # 31743158 ok
sbatch scripts/matrix_10kb_H3K27me3_Bernstein-H3K27me3_hESC_Bernstein-UpDownRegulated_Epiblastpadj05fc025.sh # 31745403 ok
```

--> Upregulated genes in Epiblast show higher level of H3K27me3 (good; hypothesis confirm). 













