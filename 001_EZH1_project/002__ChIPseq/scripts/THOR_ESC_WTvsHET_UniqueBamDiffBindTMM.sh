#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsHETUniqueBamDiffBindTMM --merge --output-dir output/THOR/THOR_ESC_WTvsHET_UniqueBamDiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.010075401,0.905242149,2.299537264,3.45113078 output/THOR/ESC_WTvsHET_UniqueBamDiffBindTMM.config



