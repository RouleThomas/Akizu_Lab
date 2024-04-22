#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name CPCNR2F2untreatedvsRA --merge --output-dir output/THOR/THOR_CPC_NR2F2_untreatedvsRA --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/CPC_NR2F2_untreatedvsRA.config




