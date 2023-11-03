#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCEZH2MG1655bam --merge --output-dir output/THOR/THOR_NPC_EZH2_MG1655bam --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/NPC_EZH2_MG1655.config

