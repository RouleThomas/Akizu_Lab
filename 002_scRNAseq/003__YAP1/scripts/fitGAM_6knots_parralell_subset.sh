#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=10




Rscript scripts/fitGAM_6knots_parralell_subset.R




