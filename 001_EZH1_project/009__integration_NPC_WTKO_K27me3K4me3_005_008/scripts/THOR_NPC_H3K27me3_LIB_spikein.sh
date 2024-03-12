#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name NPCWTvsKOH3K27me3LIBspikein --merge --output-dir output/THOR/THOR_NPC_WTvsKO_H3K27me3_LIB_spikein --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.7533433,1.3660469,0.9389445,1.1421797 output/THOR/NPC_H3K27me3.config




