#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsHETKeepdupTMM --merge --output-dir output/THOR/THOR_WTvsHET_Keepdup_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/WTvsHET.config



