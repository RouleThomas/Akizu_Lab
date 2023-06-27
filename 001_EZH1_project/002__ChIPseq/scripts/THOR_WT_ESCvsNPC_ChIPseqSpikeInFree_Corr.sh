#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCChIPseqSpikeInFreeCorr --merge --output-dir output/THOR/THOR_WT_ESCvsNPC_ChIPseqSpikeInFree_Corr --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 42.4519228,22.8861149413333,7.51087394666667,7.26075627866667 output/THOR/WT_ESCvsNPC_DiffBindTMM.config



