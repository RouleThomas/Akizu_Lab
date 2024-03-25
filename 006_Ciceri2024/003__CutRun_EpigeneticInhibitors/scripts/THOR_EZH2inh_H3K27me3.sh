#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name EZH2inhH3K27me3 --merge --output-dir output/THOR/THOR_EZH2inh_H3K27me3 --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/EZH2inh_H3K27me3.config




