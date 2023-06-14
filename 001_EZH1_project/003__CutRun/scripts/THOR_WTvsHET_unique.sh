#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsHET_unique --merge --rmdup --output-dir output/THOR_WTvsHET_unique --report --deadzones ../../Master/meta/hg38-blacklist.v2.bed --pvalue 0.1 --scaling-factors 1.412253844,1.186123259,1.337020479,1.40842928,1.738155858,1.499589188,1.940327554,1.462426177 output/THOR/WTvsHET.config



