#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1DiffBindMG1655TMM --merge --output-dir output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1_DiffBindMG1655TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.708405794,1.266995959,0.895595346,1.293214709 output/THOR/PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1.config








