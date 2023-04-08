#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=72:00:00


input_list=("NPC_WT_input_R1"
"NPC_WT_input_R2")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2_endtoend/downsample/${x}.dupmark.sorted.bam \
        --outFileName output/bigwig_downsample/${x}.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 10 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done




