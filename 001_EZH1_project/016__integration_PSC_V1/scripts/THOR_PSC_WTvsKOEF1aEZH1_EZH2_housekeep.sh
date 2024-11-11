#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1EZH2housekeep --merge --housekeeping-genes ../002__ChIPseq/output/THOR/hk_genes_new_promoter_hg38.bed --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_housekeep --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_EZH2.config







