#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCH3K27me3WTvsKOEF1aEZH1 --merge --output-dir output/THOR/THOR_PSC_H3K27me3_WTvsKOEF1aEZH1 --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_H3K27me3_WTvsKOEF1aEZH1.config







