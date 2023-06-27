#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvs2dNUniqueBamDiffBindTMM --merge --output-dir output/THOR/THOR_WT_ESCvs2dN_UniqueBamDiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.010075401,0.905242149,1.281454092,1.073748366 output/THOR/WT_ESCvs2dN_UniqueBamDiffBindTMM.config



