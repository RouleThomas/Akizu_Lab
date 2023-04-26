#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=500G
#SBATCH --time=100:00:00



Rscript scripts/DiffBind_all.R