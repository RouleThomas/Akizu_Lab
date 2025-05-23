#!/usr/bin/env bash
# ~/mutsim/scripts/make_cds_bed.sh
set -euo pipefail
gtf=$1    # gencode.v45.annotation.gtf

awk '$3=="CDS"{OFS="\t";
     chr=$1; start=$4-1; end=$5; strand=$7;
     gene=$10; tx=$12;
     gsub(/[\";]/,"",gene); gsub(/[\";]/,"",tx);
     print chr, start, end, gene, tx, strand
}' "$gtf"

