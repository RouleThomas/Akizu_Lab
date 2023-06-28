#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCrmdupTMM --merge --rmdup --output-dir output/THOR/THOR_WT_ESCvsNPC_rmdup_TMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/WT_ESCvsNPC_DiffBindTMM.config



