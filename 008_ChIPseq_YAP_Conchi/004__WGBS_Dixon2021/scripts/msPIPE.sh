#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


module load singularity/3.11.1

singularity exec \
  --bind /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/input:/scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/input:ro \
  --bind /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/output_msPIPE:/scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/output_msPIPE \
  --bind /scr1/users/roulet/Akizu_Lab/Master/meta:/msPIPE/reference \
  /scr1/users/roulet/Akizu_Lab/Master/software/msPIPE/mspipe.sif \
  msPIPE.py -p /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/scripts/params_docker.conf -o /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/output_msPIPE

