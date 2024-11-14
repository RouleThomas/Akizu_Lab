#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOEF1aEZH1H3K27me3DiffBindTMMEpiCypher --merge --scaling-factors 0.971366543,1.108432645,1.488505906,1.328793466,1.042801582,1.431746224 --output-dir output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_DiffBindTMMEpiCypher --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKOEF1aEZH1_H3K27me3.config








