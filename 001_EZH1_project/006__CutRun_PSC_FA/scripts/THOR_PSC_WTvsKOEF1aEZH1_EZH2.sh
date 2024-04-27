#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1EZH2 --merge --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.898318034,1.166797017 output/THOR/PSC_WTvsKOEF1aEZH1_EZH2.config




