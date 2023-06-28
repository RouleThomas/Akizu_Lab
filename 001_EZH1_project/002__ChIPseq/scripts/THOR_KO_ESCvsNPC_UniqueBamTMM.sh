#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00


rgt-THOR --name KOESCvsNPCUniqueBamTMM --merge --output-dir output/THOR/THOR_KO_ESCvsNPC_UniqueBamTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/KO_ESCvsNPC_UniqueBamDiffBindTMM.config



