#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G


macs2_out="output/macs2"
gs_regions="../../Master/meta/hg38-blacklist.v2.bed" #10kb merge greenscreen


# average basepair q-value threshold (log5)
q=10


input_list=("EMF2_Rep1" "EMF2_Rep2" "H3K27me3_Rep1" "H3K27me3_Rep2")

for x in "${input_list[@]}"; do

macs2 callpeak -t output/bowtie2/downsample/${x}.dupmark.sorted.bam \
  -c mapped/chip/downsample/input.dupmark.sorted.bam  \
  -f BAMPE --keep-dup auto \
  --nomodel -g 373128865 \
  --outdir ${macs2_out} -n ${x} --broad
                                          
                                                 
awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/noMask_qval${q}/${x}_peaks.broadPeak                                             
    
bedtools intersect -v -wa \
  -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
  -b ${gs_regions} > ${macs2_out}/broad_gsMask_qval${q}/${x}_peaks.broadPeak                                             
                                                 
done       