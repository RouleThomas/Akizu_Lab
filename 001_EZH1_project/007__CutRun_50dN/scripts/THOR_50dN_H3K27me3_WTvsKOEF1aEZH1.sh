#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name 50dNH3K27me3WTvsKOEF1aEZH1 --merge --output-dir output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.991352629,0.982737431,2.093518738,1.636181604 output/THOR/50dN_H3K27me3_WTvsKOEF1aEZH1.config








