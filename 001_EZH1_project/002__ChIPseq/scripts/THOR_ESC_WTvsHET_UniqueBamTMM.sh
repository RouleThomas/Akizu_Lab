#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsHETUniqueBamTMM --merge --output-dir output/THOR/THOR_ESC_WTvsHET_UniqueBamTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/ESC_WTvsHET_UniqueBamDiffBindTMM.config



