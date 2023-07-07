#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00


rgt-THOR --name HETESCvsNPCUniqueBamhousekeep --housekeeping-genes output/THOR/hk_genes_new_promoter_hg38.bed --merge --output-dir output/THOR/THOR_HET_ESCvsNPC_UniqueBamhousekeep --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/HET_ESCvsNPC_UniqueBamDiffBindTMM.config



