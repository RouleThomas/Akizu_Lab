#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCDiffBindTMM --merge --output-dir output/THOR/THOR_WT_ESCvsNPC_DiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.25660203,1.363597191,1.425558409,1.556845489 output/THOR/WT_ESCvsNPC_DiffBindTMM.config



