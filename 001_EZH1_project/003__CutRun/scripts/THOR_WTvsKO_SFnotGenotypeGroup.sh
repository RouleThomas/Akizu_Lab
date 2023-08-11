#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsKOSFnotGenotypeGroup --merge --output-dir output/THOR/THOR_WTvsKO_SFnotGenotypeGroup --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.591418337,0.341820104,0.487437186,0.323239055,0.765361475,1,0.262499154,0.885743637 output/THOR/WTvsKO.config

