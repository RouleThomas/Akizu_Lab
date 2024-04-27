#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1SUZ12 --merge --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.033214861,1.075062821 output/THOR/PSC_WTvsKOEF1aEZH1_SUZ12.config




