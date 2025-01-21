#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOH3K27me3FergusonUniqueNorm99 --merge --scaling-factors 1,0.86695279,2.06122449,3.482758621,2.376470588,1.091891892 --output-dir output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99 --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_H3K27me3.config



