#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00




cellranger aggr --id=count_50dOrga \
                --csv=meta/aggr.csv

