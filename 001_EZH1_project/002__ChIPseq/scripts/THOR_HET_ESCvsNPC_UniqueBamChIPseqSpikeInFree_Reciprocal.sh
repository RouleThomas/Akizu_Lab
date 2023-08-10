#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name HETESCvsNPCUniqueBamChIPseqSpikeInFreeRec --merge --output-dir output/THOR/THOR_HET_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Reciprocal --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.032389002,0.021445641,0.458243215,0.197663149 output/THOR/HET_ESCvsNPC_UniqueBamDiffBindTMM.config



