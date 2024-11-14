#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1EZH2EpiCypher --merge --scaling-factors 0.35516538,3.09800838,2.071952793,0.412123093,6.532051763,2.431411381 --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_EpiCypher --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_EZH2.config








