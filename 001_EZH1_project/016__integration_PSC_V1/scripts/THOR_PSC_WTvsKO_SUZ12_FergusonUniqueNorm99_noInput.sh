#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOSUZ12FergusonUniqueNorm99noInput --merge --scaling-factors 1,1.333333333,0.8,1.333333333,1,1 --output-dir output/THOR/THOR_PSC_WTvsKO_SUZ12_FergusonUniqueNorm99_noInput --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_SUZ12_noInput.config






