#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99 --merge --scaling-factors 1,0.86695279,2.06122449,1.980392157,1.141242938,1.771929825 --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99 --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_H3K27me3.config











