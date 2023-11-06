#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCH3K27me3LIBspikeinNotReciprocal --merge --output-dir output/THOR/THOR_NPC_H3K27me3_LIB_spikein_NotReciprocal --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.1096747,0.8903253 output/THOR/NPC_H3K27me3.config




