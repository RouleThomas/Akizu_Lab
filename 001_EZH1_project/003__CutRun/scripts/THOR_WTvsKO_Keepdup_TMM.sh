#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsKOKeepdupTMM --merge --output-dir output/THOR/THOR_WTvsKO_Keepdup_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/WTvsKO.config


