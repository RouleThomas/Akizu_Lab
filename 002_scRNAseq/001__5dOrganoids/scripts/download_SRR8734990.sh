#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=50:00:00




fasterq-dump SRR9169172 --include-technical -S


