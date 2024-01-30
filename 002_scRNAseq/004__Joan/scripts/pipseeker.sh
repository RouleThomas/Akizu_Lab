#! /bin/bash
#PBS -S /bin/bash
#PBS -o logs_o1_Pipseq67
#PBS -e logs_e1_Pipseq67
#PBS -q himem_01_long
#PBS -l nodes=1:ppn=16
#PBS -d .
#PBS -N PipseqRun



module load pipseeker/2.1.4

pipseeker full --fastq JFB67/. \
--star-index-path /quobyteVolumes/Home01/fontlab/tranklb/STARINDEX/STAR2.7.10b_mm10-2020-A_tdT_EGFP \
--output-path /JFB67-results-V2.1 \
--skip-version-check 
