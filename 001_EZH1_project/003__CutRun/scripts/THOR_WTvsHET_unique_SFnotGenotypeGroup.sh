#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsHETuniqueSFnotGenotypeGroup --merge --output-dir output/THOR/THOR_WTvsHET_unique_SFnotGenotypeGroup --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.591418337,0.341820104,0.487437186,0.323239055,0.63099691,0.53052574,0.951563458,0.349975195 output/THOR/WTvsHET_unique.config


