#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00




fasterq-dump SRR8734991 --include-technical -S


