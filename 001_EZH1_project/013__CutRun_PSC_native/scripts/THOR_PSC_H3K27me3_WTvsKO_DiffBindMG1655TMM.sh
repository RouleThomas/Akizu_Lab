#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCH3K27me3WTvsKODiffBindMG1655TMM --merge --output-dir output/THOR/THOR_PSC_H3K27me3_WTvsKO_DiffBindMG1655TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.725517976,1.252268328,0.867830462,1.352071604 output/THOR/PSC_H3K27me3_WTvsKO.config









