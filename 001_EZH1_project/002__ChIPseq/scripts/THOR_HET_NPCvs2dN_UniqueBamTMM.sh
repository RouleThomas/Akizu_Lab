#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00


rgt-THOR --name HETNPCvs2dNUniqueBamTMM --merge --output-dir output/THOR/THOR_HET_NPCvs2dN_UniqueBamTMM --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/HET_NPCvs2dN_UniqueBamDiffBindTMM.config



