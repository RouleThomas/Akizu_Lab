#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name hESCEZH2WTvsYAPKO --merge --output-dir output/THOR/THOR_hESC_EZH2_WTvsYAPKO --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/hESC_EZH2_WTvsYAPKO.config




