#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCWTvsKOH3K4me3 --merge --output-dir output/THOR/THOR_NPC_WTvsKO_H3K4me3 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.0857384,0.7577088,0.7403885,1.7061478  output/THOR/NPC_NPCH3K4me3.config




