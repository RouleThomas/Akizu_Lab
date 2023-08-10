#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsHETuniqueKeepdup --merge --output-dir output/THOR/THOR_WTvsHET_unique_Keepdup --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.37403632,1.195845537,1.30622362,1.403077904,1.749132736,1.50867655,1.980736937,1.435045328 output/THOR/WTvsHET_unique.config

