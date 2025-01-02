#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOH3K27me3DiffBindTMM --merge --scaling-factors 0.761284887,0.868707086,1.166358806,2.225307487,1.391972301,0.756462554 --output-dir output/THOR/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_H3K27me3.config









