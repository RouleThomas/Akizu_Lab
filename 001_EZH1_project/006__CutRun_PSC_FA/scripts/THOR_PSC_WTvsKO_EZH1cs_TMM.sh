#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEZH1csTMM --merge --output-dir output/THOR/THOR_PSC_WTvsKO_EZH1cs_TMM --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_EZH1cs.config




