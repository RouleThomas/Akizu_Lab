#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCH3K27me3WTvsKOEF1aEZH1DiffBindMG1655TMM --merge --output-dir output/THOR/THOR_PSC_H3K27me3_WTvsKOEF1aEZH1_DiffBindMG1655TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.812187029,1.417445695,0.773020592,1.139046868 output/THOR/PSC_H3K27me3_WTvsKOEF1aEZH1.config


 









