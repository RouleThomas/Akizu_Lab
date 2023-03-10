#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL



x=("4wN_WT_R1" "4wN_WT_R2" "4wN_KO_R1"
   "4wN_KO_R2" "4wN_HET_R1" "4wN_HET_R2"
   "4wN_HET_R3" "4wN_HET_R4" "4wN_iPSCWT_R1"
   "4wN_iPSCWT_R2" "4wN_iPSCpatient_R1" "4wN_iPSCpatient_R2"
   "8wN_WT_R1" "8wN_WT_R2" "8wN_KO_R1"
   "8wN_KO_R2" "8wN_HET_R1" "8wN_HET_R2"
   "8wN_HET_R3" "8wN_HET_R4" "8wN_iPSCWT_R1"
   "8wN_iPSCWT_R2" "8wN_iPSCpatient_R1" "8wN_iPSCpatient_R2")

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done