#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



cellranger count --id=humangastruloid_DASATINIB72hr \
                   --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/ds.5ade69eb01db42bfa60310d0200a2d2d \
                   --sample=humangastruloid_DASATINIB72hr

                   