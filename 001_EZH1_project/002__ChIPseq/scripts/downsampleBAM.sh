#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G
#SBATCH --time=100:00:00


in_dir="output/bowtie2_endtoend"
out_dir="output/bowtie2_endtoend/downsample"
JVARKIT_PATH="/scr1/users/roulet/Akizu_Lab/Master/software/JVARKIT"
seed=42
mkdir -p ${out_dir}

module load sam-bcf-tools/1.6
module load picard/2.26.10-Java-15

# function to find the minimum
# value in an array
minIndex(){
   arr=("$@")
   min_val=${arr[0]}
   min_idx=0
   for i in ${!arr[@]}; do
        cur_val=${arr[${i}]}
        if [[ ${cur_val} -lt ${min_val} ]]; then
                min_val=${arr[$i]}
                min_idx=${i}
        fi
   done

}


# There is no need to down-sample anything
# with a single replicate
single_rep=("")
for samp in "${single_rep[@]}"; do
  cp ${in_dir}/${samp}.dupmark.sorted.bam \
    ${out_dir}/${samp}.dupmark.sorted.bam
  cp ${in_dir}/${samp}.dupmark.sorted.bam.bai \
    ${out_dir}/${samp}.dupmark.sorted.bam.bai
done

# down-sample given two reps
pool_two=("2dN_HET_H3K27me3" "2dN_HET_input"
          "2dN_KO_H3K27me3" "2dN_KO_input" 
          "2dN_WT_H3K27me3" "2dN_WT_input"
          "ESC_HET_H3K27me3" "ESC_HET_input"
          "ESC_KO_H3K27me3" "ESC_KO_input" 
          "NPC_HET_H3K27me3" "NPC_HET_input"
          "NPC_KO_H3K27me3" "NPC_KO_input" 
          "NPC_WT_H3K27me3" "NPC_WT_input")
for samp in "${pool_two[@]}"; do
  depth1=`samtools view -c \
    ${in_dir}/${samp}_R1.dupmark.sorted.bam`
  depth2=`samtools view -c \
    ${in_dir}/${samp}_R2.dupmark.sorted.bam`
  arrName=(${depth1} ${depth2})
  minIndex "${arrName[@]}"
  let "min_idx = $min_idx + 1"
  for (( rep=1; rep<=2; rep++ )); do
    if [[ ${min_idx} -eq ${rep} ]]; then
        # copy this replicate which
        # has the smallest read depth
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam.bai \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam.bai
    else
        # downsample these replicates
        java -Xmx2g -Djava.io.tmpdir=output/tmp -jar ${JVARKIT_PATH}/jvarkit.jar biostar145820 \
            --seed ${seed} -n ${min_val} \
            -o ${out_dir}/${samp}_R${rep}.dupmark.bam \
            ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam
	samtools sort \
            -o ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam \
            ${out_dir}/${samp}_R${rep}.dupmark.bam && \
            rm ${out_dir}/${samp}_R${rep}.dupmark.bam
	samtools index ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
    fi
  done
done




# down-sample given three reps
pool_three=("ESC_WT_H3K27me3" "ESC_WT_input")
for samp in "${pool_three[@]}"; do
  depth1=`samtools view -c \
    ${in_dir}/${samp}_R1.dupmark.sorted.bam`
  depth2=`samtools view -c \
    ${in_dir}/${samp}_R2.dupmark.sorted.bam`
  depth3=`samtools view -c \
    ${in_dir}/${samp}_R3.dupmark.sorted.bam`
  arrName=(${depth1} ${depth2} ${depth3})
  minIndex "${arrName[@]}"
  let "min_idx = $min_idx + 1"
  for (( rep=1; rep<=3; rep++ )); do
    if [[ ${min_idx} -eq ${rep} ]]; then
        # copy this replicate which
        # has the smallest read depth
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam.bai \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam.bai
    else
        # downsample these replicates
        java -jar -Xmx2g -Djava.io.tmpdir=output/tmp -jar ${JVARKIT_PATH}/jvarkit.jar biostar145820 \
            --seed ${seed} -n ${min_val} \
            -o ${out_dir}/${samp}_R${rep}.dupmark.bam \
            ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam
	samtools sort \
            -o ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam \
            ${out_dir}/${samp}_R${rep}.dupmark.bam && \
            rm ${out_dir}/${samp}_R${rep}.dupmark.bam
	samtools index ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
    fi
  done
done

# down-sample given four reps
pool_four=("")
for samp in "${pool_four[@]}"; do
  depth1=`samtools view -c \
    ${in_dir}/${samp}_R1.dupmark.sorted.bam`
  depth2=`samtools view -c \
    ${in_dir}/${samp}_R2.dupmark.sorted.bam`
  depth3=`samtools view -c \
    ${in_dir}/${samp}_R3.dupmark.sorted.bam`
  depth4=`samtools view -c \
    ${in_dir}/${samp}_R4.dupmark.sorted.bam`
  arrName=(${depth1} ${depth2} ${depth3} ${depth4})
  minIndex "${arrName[@]}"
  let "min_idx = $min_idx + 1"
  for (( rep=1; rep<=4; rep++ )); do
    if [[ ${min_idx} -eq ${rep} ]]; then
        # copy this replicate which
        # has the smallest read depth
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
    else
        # downsample these replicates
        java -jar -Xmx2g -Djava.io.tmpdir=output/tmp -jar ${JVARKIT_PATH}/jvarkit.jar biostar145820 \
            --seed ${seed} -n ${min_val} \
            -o ${out_dir}/${samp}_R${rep}.dupmark.bam \
            ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam
	samtools sort \
            -o ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam \
            ${out_dir}/${samp}_R${rep}.dupmark.bam && \
            rm ${out_dir}/${samp}_R${rep}.dupmark.bam
	samtools index ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
    fi
  done
done