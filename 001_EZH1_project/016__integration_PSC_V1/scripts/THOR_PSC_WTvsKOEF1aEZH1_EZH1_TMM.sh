#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1EZH1TMM --merge --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH1_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_EZH1.config







