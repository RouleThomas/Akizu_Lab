#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCEZH2LIBspikein --merge --output-dir output/THOR/THOR_NPC_EZH2_LIB_spikein --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.945517665,1.061144978 output/THOR/NPC_EZH2.config

