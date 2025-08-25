#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsOEKOH3K27me3TMM --merge --output-dir output/THOR/THOR_ESC_WTvsOEKO_H3K27me3_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/ESC_WTvsOEKO_H3K27me3.config







