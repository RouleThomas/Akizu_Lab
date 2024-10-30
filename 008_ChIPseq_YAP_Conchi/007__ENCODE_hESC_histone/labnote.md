# Project and goal

Part of gastrulation paper with Conchi Estaras.

--> See `002*/003*/gastrulation paper/Figure 4` for detail: We need to know the overlap of the 199 QSER1:YAP cobound peaks (identified in `008*/001*` and `008*/003*`) with H3K27ac, H3K4me1, H3K27me3 and H3K36me3 modification marks

Pipeline:
- Identify data from ENCODE
- Download ENCODE data (bigwig and peak file)
- Follow `Figure 4` ppt

# Data download from ENCODE

Prioritize most recent datasets, or all datasets from same lab:
- **Bing Ren**, UCSD_Project ROADMAP; [2013 2 Bio Rep; processed 2020](https://www.encodeproject.org/experiments/ENCSR928HYM/):
    - **H3K27me3**: *ENCFF395GVR* bigwig signal p-value; *ENCFF599KDF* bed
    - **H3K4me1**:  *ENCFF164XHJ* bigwig signal p-value; *ENCFF613QAB* bed
    - **H3K36me3**: *ENCFF483UMR* bigwig signal p-value; *ENCFF681CEO* bed
    - **H3K27ac**:  *ENCFF390JIZ* bigwig signal p-value; *ENCFF045CUG* bed
- **Bradley Bernstein**, Broad_Project ENCODE; 2013
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



