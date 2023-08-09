#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name KOESCvsNPCUniqueBamDiffBindTMM --merge --output-dir output/THOR/THOR_KO_ESCvsNPC_UniqueBamDiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 3.703952006,4.454690232,1.46400002,1.106896051 output/THOR/KO_ESCvsNPC_UniqueBamDiffBindTMM.config



