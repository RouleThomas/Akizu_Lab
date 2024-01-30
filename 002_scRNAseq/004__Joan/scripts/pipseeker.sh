#! /bin/bash
#PBS -S /bin/bash
#PBS -o logs_o1_Pipseq67
#PBS -e logs_e1_Pipseq67
#PBS -q himem_01_long
#PBS -l nodes=1:ppn=16
#PBS -d .
#PBS -N PipseqRun

#!/bin/bash
#SBATCH --output=logs_o1_Pipseq67   # Equivalent to #PBS -o, sets the file for standard output
#SBATCH --error=logs_e1_Pipseq67    # Equivalent to #PBS -e, sets the file for standard error
#SBATCH --partition=himem_01_long   # Equivalent to #PBS -q, specifies the partition (queue) name
#SBATCH --nodes=1                   # Equivalent to #PBS -l nodes=1, requests 1 node
#SBATCH --ntasks-per-node=16        # Equivalent to #PBS -l ppn=16, requests 16 tasks per node (cores)
#SBATCH --workdir=.                 # Equivalent to #PBS -d, sets the working directory to the current directory
#SBATCH --job-name=PipseqRun        # Equivalent to #PBS -N, sets the job name

# Additional directives for resource requests in Slurm (not directly mapped from your PBS script but potentially useful):
#SBATCH --mem=200G                  # Requests memory, equivalent to #SBATCH --mem in your usual Slurm script
#SBATCH --time=200:00:00            # Sets a time limit, equivalent to #SBATCH --time in your usual Slurm script



module load pipseeker/2.1.4

pipseeker full --fastq JFB67/. \
--star-index-path /quobyteVolumes/Home01/fontlab/tranklb/STARINDEX/STAR2.7.10b_mm10-2020-A_tdT_EGFP \
--output-path /JFB67-results-V2.1 \
--skip-version-check 
