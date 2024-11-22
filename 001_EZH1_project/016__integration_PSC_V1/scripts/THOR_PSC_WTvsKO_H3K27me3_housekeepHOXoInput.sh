#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOH3K27me3housekeepHOXnoInput --merge --housekeeping-genes meta/ENCFF159KBI_HOX_genes.bed --output-dir output/THOR/THOR_PSC_WTvsKO_H3K27me3_housekeepHOXnoInput --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_H3K27me3_noInput.config






