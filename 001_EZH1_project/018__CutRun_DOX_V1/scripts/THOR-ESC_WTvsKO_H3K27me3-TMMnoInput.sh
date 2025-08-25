#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name ESCWTvsKOH3K27me3TMMnoInput --merge --output-dir output/THOR/THOR_ESC_WTvsKO_H3K27me3_TMMnoInput --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/ESC_WTvsKO_H3K27me3_noInput.config







