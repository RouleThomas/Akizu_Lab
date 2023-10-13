#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCSUZ12 --merge --output-dir output/THOR/THOR_NPC_SUZ12 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.514593869,0.713876255 output/THOR/NPC_SUZ12.config




