#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00


rgt-THOR --name KONPCvs2dNUniqueBamTMM --merge --output-dir output/THOR/THOR_KO_NPCvs2dN_UniqueBamTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/KO_NPCvs2dN_UniqueBamDiffBindTMM.config



