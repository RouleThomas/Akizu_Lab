#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEZH2DiffBindTMMEpiCypher --merge --scaling-factors 1.208227352,1.220149649,0.965529347,2.417191258,1.072331436,1.115060283 --output-dir output/THOR/THOR_PSC_WTvsKO_EZH2_DiffBindTMMEpiCypher --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_EZH2.config









