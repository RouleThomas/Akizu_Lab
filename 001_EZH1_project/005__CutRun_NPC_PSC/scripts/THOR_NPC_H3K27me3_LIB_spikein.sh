#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCH3K27me3LIBspikein --merge --output-dir output/THOR/THOR_NPC_H3K27me3_LIB_spikein --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.901164999,1.123184975 output/THOR/NPC_H3K27me3.config




