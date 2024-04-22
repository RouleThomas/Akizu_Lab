#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name hESCQSER1WTvsYAPKO --merge --output-dir output/THOR/THOR_hESC_QSER1_WTvsYAPKO --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/hESC_QSER1_WTvsYAPKO.config




