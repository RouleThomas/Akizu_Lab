#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00



/scr1/users/roulet/Akizu_Lab/Master/software/cellranger-6.0.2/cellranger aggr --id=count \
                --csv=meta/aggr_corr.csv
                
                