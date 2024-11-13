#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1EZH1housekeepHOXHKnoInput --merge --housekeeping-genes meta/ENCFF159KBI_HOX_HK_genes.bed --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH1_housekeepHOXHKnoInput --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_EZH1_noInput.config







