#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsKOChIPseqSpikeInFreeCorr --merge --output-dir output/THOR/THOR_ESC_WTvsKO_ChIPseqSpikeInFree_Corr --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 42.4519228,22.8861149413333,57.6873544346667,60.64140384 output/THOR/ESC_WTvsKO_DiffBindTMM.config



