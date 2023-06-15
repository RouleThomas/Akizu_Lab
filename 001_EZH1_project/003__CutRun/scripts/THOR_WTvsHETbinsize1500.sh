#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsHET --merge --binsize 1500 --step 750 --output-dir output/THOR/THOR_WTvsHETbinsize1500 --report --deadzones ../../Master/meta/hg38-blacklist.v2.bed --pvalue 0.1 --scaling-factors 1.412253844,1.186123259,1.337020479,1.40842928,1.738155858,1.499589188,1.940327554,1.462426177 output/THOR/WTvsHET.config



