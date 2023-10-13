#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCEZH2 --merge --output-dir output/THOR/THOR_NPC_EZH2 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.266595568,0.792209099 output/THOR/NPC_EZH2.config

