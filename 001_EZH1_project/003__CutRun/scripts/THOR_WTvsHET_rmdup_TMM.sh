#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsHETrmdupTMM --merge --rmdup --output-dir output/THOR/THOR_WTvsHET_rmdup_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/WTvsHET.config



