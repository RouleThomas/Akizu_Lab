#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1H3K27me3 --merge --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.807177487,0.879636133 output/THOR/PSC_WTvsKOEF1aEZH1_H3K27me3.config




