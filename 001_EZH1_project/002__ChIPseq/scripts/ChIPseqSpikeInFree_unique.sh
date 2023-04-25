#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=150G
#SBATCH --time=150:00:00


Rscript scripts/ChIPseqSpikeInFree_unique.R