#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s1-rep0.bw \
    output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s2-rep0.bw \
    output/THOR/THOR_NPC_EZH2/NPCEZH2-s1-rep0.bw  \
    output/THOR/THOR_NPC_EZH2/NPCEZH2-s2-rep0.bw \
    output/THOR/THOR_NPC_H3K27me3/NPCH3K27me3-s1-rep0.bw \
    output/THOR/THOR_NPC_H3K27me3/NPCH3K27me3-s2-rep0.bw \
    output/THOR/THOR_NPC_H3K4me3/NPCH3K4me3-s1-rep0.bw \
    output/THOR/THOR_NPC_H3K4me3/NPCH3K4me3-s2-rep0.bw \
    -o output/THOR/multiBigwigSummary_NPC_THOR.npz








