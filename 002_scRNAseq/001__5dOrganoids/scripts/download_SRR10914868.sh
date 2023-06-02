#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00






fasterq-dump SRR10914868 --include-technical -S