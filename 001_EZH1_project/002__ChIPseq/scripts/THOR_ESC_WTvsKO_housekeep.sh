#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsKOhousekeep --housekeeping-genes output/THOR/hk_genes_new_promoter_hg38.bed --merge --output-dir output/THOR/THOR_ESC_WTvsKO_housekeep --deadzones /scr1/users/roulet/Akizu_Lab/Master/meta/hg38-blacklist.v2.bed output/THOR/ESC_WTvsKO_DiffBindTMM.config



