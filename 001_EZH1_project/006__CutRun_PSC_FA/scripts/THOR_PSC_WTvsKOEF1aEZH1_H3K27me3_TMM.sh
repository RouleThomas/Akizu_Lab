#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1H3K27me3TMM --merge --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_H3K27me3.config




