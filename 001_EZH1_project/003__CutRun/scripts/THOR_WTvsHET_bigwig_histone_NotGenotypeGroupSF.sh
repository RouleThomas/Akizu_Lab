#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


rgt-THOR --name WTvsHET --merge --output-dir output/THOR/THOR_WTvsHET_NotGenotypeGroupSF --report --deadzones ../../Master/meta/hg38-blacklist.v2.bed --pvalue 0.1 --scaling-factors 0.5914183370169948,0.3418201039555981,0.4874371859296476,0.3232390552755446,0.6309969100666774,0.5305257400697344,0.9515634580012263,0.3499751950570516 output/THOR/WTvsHET.config



