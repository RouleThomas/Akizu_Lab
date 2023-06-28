#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsKOrmdup --merge --rmdup --output-dir output/THOR/THOR_WTvsKO_rmdup --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 1.412253844,1.186123259,1.337020479,1.40842928,1.502057669,1.487749425,0.449980991,1.883122496  output/THOR/WTvsKO.config


