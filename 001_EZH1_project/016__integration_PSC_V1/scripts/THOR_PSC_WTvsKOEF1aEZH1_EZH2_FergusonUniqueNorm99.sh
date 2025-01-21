#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1EZH2FergusonUniqueNorm99 --merge --scaling-factors 1,1,1,1.25,1.25,1 --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_FergusonUniqueNorm99 --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_EZH2.config





