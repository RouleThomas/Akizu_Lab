#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsKOH3K27me3housekeepHOX --merge --housekeeping-genes meta/ENCFF159KBI_HOX_genes.bed --output-dir output/THOR/THOR_ESC_WTvsKO_H3K27me3_housekeepHOX --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/ESC_WTvsKO_H3K27me3.config







