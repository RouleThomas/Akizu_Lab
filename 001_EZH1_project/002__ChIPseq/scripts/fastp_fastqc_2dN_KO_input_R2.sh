#!/bin/bash

fastp -i input/2dN_KO_input_R2_1.fq.gz -I input/2dN_KO_input_R2_2.fq.gz \
      -o output/fastp/2dN_KO_input_R2_1.fq.gz -O output/fastp/2dN_KO_input_R2_2.fq.gz \
	  -h output/fastp/2dN_KO_input_R2 -j output/fastp/2dN_KO_input_R2

fastqc -o output/fastqc/fastp output/fastp/2dN_KO_input_R2_1.fq.gz
fastqc -o output/fastqc/fastp output/fastp/2dN_KO_input_R2_2.fq.gz