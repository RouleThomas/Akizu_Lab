#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsKOUniqueBamDiffBindTMM --merge --output-dir output/THOR/THOR_ESC_WTvsKO_UniqueBamDiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.010075401,0.905242149,3.703952006,4.454690232 output/THOR/ESC_WTvsKO_UniqueBamDiffBindTMM.config



