#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsHETChIPseqSpikeInFreeCorr --merge --output-dir output/THOR/THOR_ESC_WTvsHET_ChIPseqSpikeInFreeCorr --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 42.4519228,22.8861149413333,60.2108749293333,103.648971913333 output/THOR/ESC_WTvsHET_DiffBindTMM.config



