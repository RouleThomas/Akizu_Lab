#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1EZH2DiffBindTMMEpiCypher --merge --scaling-factors 1.208227352,1.220149649,0.965529347,1.640599725,2.329875953,0.978209407 --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_DiffBindTMMEpiCypher --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_EZH2.config










