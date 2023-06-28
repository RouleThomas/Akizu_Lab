#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsKOrmdupTMM --merge --rmdup --output-dir output/THOR/THOR_WTvsKO_rmdup_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/WTvsKO.config


