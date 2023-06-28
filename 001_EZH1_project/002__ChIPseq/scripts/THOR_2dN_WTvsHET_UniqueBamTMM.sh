#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name 2dNWTvsHETUniqueBamTMM --merge --output-dir output/THOR/THOR_2dN_WTvsHET_UniqueBamTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/2dN_WTvsHET_UniqueBamDiffBindTMM.config



