#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name 50dNH3K27me3WTvsKO --merge --output-dir output/THOR/THOR_50dN_H3K27me3_WTvsKO --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.991352629,0.982737431,0.847869994,0.736704256 output/THOR/50dN_H3K27me3_WTvsKO.config








