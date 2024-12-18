#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1SUZ12DiffBindTMMEpiCypher --merge --scaling-factors 1.859094056,2.855253495,1.024947635,0.625249357,1.826501366,2.264176917 --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_SUZ12_DiffBindTMMEpiCypher --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_SUZ12.config











