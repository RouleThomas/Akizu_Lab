#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOH3K27me3DiffBindTMMEpiCypher --merge --scaling-factors 0.971366543,1.108432645,1.488505906,2.839471608,1.776552351,0.965213512 --output-dir output/THOR/THOR_PSC_WTvsKO_H3K27me3_DiffBindTMMEpiCypher --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_H3K27me3.config








