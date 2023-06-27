#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvs2dNDiffBindTMM --merge --output-dir output/THOR/THOR_WT_ESCvs2dN_DiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.25660203,1.363597191,1.548264434,1.691816077 output/THOR/WT_ESCvs2dN_DiffBindTMM.config



