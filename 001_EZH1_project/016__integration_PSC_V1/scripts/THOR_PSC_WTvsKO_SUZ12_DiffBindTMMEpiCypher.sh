#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOSUZ12DiffBindTMMEpiCypher --merge --scaling-factors 1.859094056,2.855253495,1.024947635,4.030675049,1.344660394,1.214055165 --output-dir output/THOR/THOR_PSC_WTvsKO_SUZ12_DiffBindTMMEpiCypher --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_SUZ12.config










