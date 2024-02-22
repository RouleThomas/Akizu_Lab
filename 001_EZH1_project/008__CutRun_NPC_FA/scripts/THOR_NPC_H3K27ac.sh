#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCWTvsKOH3K27ac --merge --output-dir output/THOR/THOR_NPC_WTvsKO_H3K27ac --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.062330996,1.139698939 output/THOR/NPC_NPCH3K27ac.config




