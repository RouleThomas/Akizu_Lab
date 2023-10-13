#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name NPCH3K27me3 --merge --output-dir output/THOR/THOR_NPC_H3K27me3 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.170001033,0.870425955 output/THOR/NPC_H3K27me3.config




