#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1DiffBindhistoneTMM --merge --output-dir output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_DiffBindhistoneTMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.757059179,1.219717613,0.970300273,1.425501712 output/THOR/PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1.config









