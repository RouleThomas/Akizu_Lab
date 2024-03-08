#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOH3K27me3 --merge --output-dir output/THOR/THOR_PSC_WTvsKO_H3K27me3 --deadzones ../../Master/meta/hg38-blacklist.v2.bed --scaling-factors 0.807177487,2.725074088 output/THOR/PSC_WTvsKO_H3K27me3.config




