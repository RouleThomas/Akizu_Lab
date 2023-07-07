#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name WTESCvsNPCuniqueBAMhousekeep --housekeeping-genes output/THOR/hk_genes_new_promoter_hg38.bed --merge --output-dir output/THOR/THOR_WT_ESCvsNPC_uniqueBAMhousekeep --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/WT_ESCvsNPC_UniqueBamDiffBindTMM.config



