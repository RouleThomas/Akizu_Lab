#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1SUZ12FergusonUniqueNorm99noInput --merge --scaling-factors 1,1.333333333,0.8,0.4,1,1 --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_FergusonUniqueNorm99_noInput --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_SUZ12_noInput.config








