#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCUniqueBamChIPseqSpikeInFreeReciprocal --merge --output-dir output/THOR/THOR_WT_ESCvsNPC_UniqueBamChIPseqSpikeInFree_Reciprocal --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.03736646,0.059912156,0.195833539,0.194757834 output/THOR/WT_ESCvsNPC_UniqueBamDiffBindTMM.config



