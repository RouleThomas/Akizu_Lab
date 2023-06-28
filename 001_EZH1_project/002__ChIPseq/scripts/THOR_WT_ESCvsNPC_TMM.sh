#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCTMM --merge --output-dir output/THOR/THOR_WT_ESCvsNPC_TMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/WT_ESCvsNPC_DiffBindTMM.config



