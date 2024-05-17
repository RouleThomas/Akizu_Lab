#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCSUZ12TMM --merge --output-dir output/THOR/THOR_NPC_SUZ12_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/NPC_SUZ12.config




