#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



cellranger count --id=humangastruloid_UNTREATED72hr \
                   --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/ds.9c1b2e5bdcd04a02b36fc6bc84d54461 \
                   --sample=humangastruloid_UNTREATED72hr

                   