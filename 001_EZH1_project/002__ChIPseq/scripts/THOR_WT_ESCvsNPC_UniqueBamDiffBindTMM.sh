#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCUniqueBamDiffBindTMM --merge --output-dir output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.010075401,0.905242149,1.029216259,1.047211652 output/THOR/WT_ESCvsNPC_UniqueBamDiffBindTMM.config



