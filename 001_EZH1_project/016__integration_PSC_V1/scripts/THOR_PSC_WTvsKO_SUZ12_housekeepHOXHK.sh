#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00


rgt-THOR --name PSCWTvsKOSUZ12housekeepHOXHK --merge --housekeeping-genes meta/ENCFF159KBI_HOX_HK_genes.bed --output-dir output/THOR/THOR_PSC_WTvsKO_SUZ12_housekeepHOXHK --deadzones ../../Master/meta/hg38-blacklist.v2.bed output/THOR/PSC_WTvsKO_SUZ12.config







