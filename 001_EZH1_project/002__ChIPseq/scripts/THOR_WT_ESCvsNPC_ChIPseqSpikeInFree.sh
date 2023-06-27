#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCChIPseqSpikeInFree --merge --output-dir output/THOR/THOR_WT_ESCvsNPC_ChIPseqSpikeInFree --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed --scaling-factors 7,4.31,1.45,1.51 output/THOR/WT_ESCvsNPC_DiffBindTMM.config



