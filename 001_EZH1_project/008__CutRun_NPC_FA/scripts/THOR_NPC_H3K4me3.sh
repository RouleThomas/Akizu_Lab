#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCWTvsKOH3K4me3 --merge --output-dir output/THOR/THOR_NPC_WTvsKO_H3K4me3 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.695790606,1.459619623 output/THOR/NPC_NPCH3K4me3.config




