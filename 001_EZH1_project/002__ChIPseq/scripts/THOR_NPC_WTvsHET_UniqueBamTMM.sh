#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name NPCWTvsHETUniqueBamTMM --merge --output-dir output/THOR/THOR_NPC_WTvsHET_UniqueBamTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/NPC_WTvsHET_UniqueBamDiffBindTMM.config



