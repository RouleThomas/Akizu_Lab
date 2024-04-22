#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name CPCTEAD4untreatedvsRA --merge --output-dir output/THOR/THOR_CPC_TEAD4_untreatedvsRA --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/CPC_TEAD4_untreatedvsRA.config




