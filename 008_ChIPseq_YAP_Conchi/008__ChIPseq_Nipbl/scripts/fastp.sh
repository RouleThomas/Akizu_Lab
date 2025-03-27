#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"Nipbl"
"input"
)

for x in "${x[@]}"; do
    fastp -i input/${x}.fq.gz  \
    -o output/fastp/${x}.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




