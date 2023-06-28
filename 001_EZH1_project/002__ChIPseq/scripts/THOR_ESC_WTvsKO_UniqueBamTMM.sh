#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsKOUniqueBamTMM --merge --output-dir output/THOR/THOR_ESC_WTvsKO_UniqueBamTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/ESC_WTvsKO_UniqueBamDiffBindTMM.config



