#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


rgt-THOR --name CPCYAP1untreatedvsRA --merge --output-dir output/THOR/THOR_CPC_YAP1_untreatedvsRA --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/CPC_YAP1_untreatedvsRA.config




