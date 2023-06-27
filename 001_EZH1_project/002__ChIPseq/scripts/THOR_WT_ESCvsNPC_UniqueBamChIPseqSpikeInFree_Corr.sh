#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCUniqueBamChIPseqSpikeInFreeCorr --merge --output-dir output/THOR/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Corr --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 26.7619678666667,16.6911036853333,5.10637761333333,5.134581652 output/THOR/WT_ESCvsNPC_UniqueBamDiffBindTMM.config



