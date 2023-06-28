#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00


rgt-THOR --name HETESCvsNPCUniqueBamTMM --merge --output-dir output/THOR/THOR_HET_ESCvsNPC_UniqueBamTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/HET_ESCvsNPC_UniqueBamDiffBindTMM.config



