#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name HETESCvsNPCUniqueBamDiffBindTMM --merge --output-dir output/THOR/THOR_HET_ESCvsNPC_UniqueBamDiffBindTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 2.299537264,3.45113078,2.22126017,0.986039553 output/THOR/HET_ESCvsNPC_UniqueBamDiffBindTMM.config



