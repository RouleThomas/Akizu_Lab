#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCH3K4me3TMM --merge --output-dir output/THOR/THOR_NPC_H3K4me3_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/NPC_H3K4me3.config




