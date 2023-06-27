#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsKOUniqueBamChIPseqSpikeInFreeCorr --merge --output-dir output/THOR/THOR_ESC_WTvsKO_UniqueBamChIPseqSpikeInFree_Corr --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 26.7619678666667,16.6911036853333,18.4690520586667,24.285081256 output/THOR/ESC_WTvsKO_UniqueBamDiffBindTMM.config



