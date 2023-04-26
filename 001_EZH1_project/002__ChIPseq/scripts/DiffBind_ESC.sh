#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --time=100:00:00



Rscript scripts/DiffBind_ESC.R