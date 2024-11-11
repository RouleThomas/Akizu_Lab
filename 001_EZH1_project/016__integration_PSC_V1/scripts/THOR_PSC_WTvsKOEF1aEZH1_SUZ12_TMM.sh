#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1SUZ12TMM --merge --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_SUZ12.config







