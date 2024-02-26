#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCWTvsKOH3K27me3 --merge --output-dir output/THOR/THOR_NPC_WTvsKO_H3K27me3 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.1253542,1.1712504,0.8471138,1.0903367 output/THOR/NPC_H3K27me3.config




