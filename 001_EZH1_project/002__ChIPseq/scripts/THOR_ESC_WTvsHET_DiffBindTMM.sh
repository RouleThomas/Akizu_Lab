#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsHETDiffBindTMM --merge --output-dir output/THOR/THOR_ESC_WTvsHET_DiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.25660203,1.363597191,2.158245122,2.796412762 output/THOR/ESC_WTvsHET_DiffBindTMM.config



