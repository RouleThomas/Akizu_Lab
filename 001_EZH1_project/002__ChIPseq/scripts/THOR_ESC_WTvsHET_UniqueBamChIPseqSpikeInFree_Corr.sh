#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsHETUniqueBamChIPseqSpikeInFreeCorr --merge --output-dir output/THOR/THOR_ESC_WTvsHET_UniqueBamChIPseqSpikeInFree_Corr --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 26.7619678666667,16.6911036853333,30.8746776773333,46.6295234733333 output/THOR/ESC_WTvsHET_UniqueBamDiffBindTMM.config



