#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCH3K4me3LIBspikein --merge --output-dir output/THOR/THOR_NPC_H3K4me3_LIB_spikein --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.930216896,1.081102236 output/THOR/NPC_H3K4me3.config




