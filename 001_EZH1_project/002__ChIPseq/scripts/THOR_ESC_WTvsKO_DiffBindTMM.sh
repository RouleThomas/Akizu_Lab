#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsKODiffBindTMM --merge --output-dir output/THOR/THOR_ESC_WTvsKO_DiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.25660203,1.363597191,2.237804803,3.1922803 output/THOR/ESC_WTvsKO_DiffBindTMM.config



