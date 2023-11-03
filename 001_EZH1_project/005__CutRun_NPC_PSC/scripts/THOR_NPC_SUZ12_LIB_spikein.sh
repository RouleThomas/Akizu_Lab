#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCSUZ12LIBspikein --merge --output-dir output/THOR/THOR_NPC_SUZ12_LIB_spikein --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.798455531,1.337645725 output/THOR/NPC_SUZ12.config




