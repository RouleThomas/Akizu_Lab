#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsKOuniqueKeepdup --merge --output-dir output/THOR/THOR_WTvsKO_unique_Keepdup --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.37403632,1.195845537,1.30622362,1.403077904,1.494404427,1.504153193,0.458616468,1.902039633 output/THOR/WTvsKO_unique.config

