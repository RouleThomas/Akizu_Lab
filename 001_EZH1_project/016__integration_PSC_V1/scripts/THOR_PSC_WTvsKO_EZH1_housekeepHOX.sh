#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEZH1housekeepHOX --merge --housekeeping-genes meta/ENCFF159KBI_HOX_genes.bed --output-dir output/THOR/THOR_PSC_WTvsKO_EZH1_housekeepHOX --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_EZH1.config







