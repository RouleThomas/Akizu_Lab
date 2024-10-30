# Project and goals 

Re-analysis of WGBS (Whole Genome Bisulfit Sequencing) dataset from (Dixon et al)[10.1126/science.abd0875].

--> 2 rep WT H1 hESC


# Download data


- Go to sra (explorer)[https://sra-explorer.info/]
- Search Bioproject  PRJNA631028 (RNAseq diff)
- Add to collections and select `Bash script for downloading FastQ files` --> copy into `scripts/download_urls.sh`

```bash
sbatch scripts/download_urls.sh # 18603395 ok
```

--> Seems like no input samples needed for WGBS (indeed it is just presence absence of mC)



## Rename files

Let's rename file with our classic nomenclature

**make sure to convert the `rename_008004.txt` into unix tab sep  format with `dos2unix`!!**

```bash
cd input

while IFS=$'\t' read -r old_name new_name
do
    mv "$old_name" "$new_name"
done < rename_008004.txt
```

--> All good 



# Data analysis

## msPIPE 

Let's try [msPIPE](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04925-2) for the analysis. [Github](https://github.com/jkimlab/msPIPE). Seems it do everyhting from preprocessing to downstream analysis! *Paper from 2022 so quite recent*


### Run through Docker


#### Install in a conda env

Below the requirements:
- R >= 4.1.3, python2 and python3
- Trim Galore
- Samtools
- cutadapt
- bowtie, bowtie2
- Bismark
- BS-Seeker2
- Other required R packages can be installed via msPIPE/bin/script/Package_install.R.


Create a conda env msPIPE

```bash
# attempt1
conda create -n msPIPE r-base=4.1.3 python

conda activate msPIPE
conda install -c conda-forge/label/cf202003 python=2.7 # fail here, impossible to install 2 python!!

conda env remove --name msPIPE

# attempt2
conda create -n msPIPE r-base=4.1.3 python=3.8
conda install trim-galore samtools cutadapt bowtie bowtie2 bismark bs-seeker2 # fail here

# attempt3 - copy bowtie2 env and install the remaining stuff

conda create --name msPIPE --clone bowtie2
conda activate msPIPE
conda install -c bioconda -c conda-forge -c anacond trim-galore # fail
conda install trim-galore 

```

--> Nightmare...




#### Install with Docker through singularity

Need to do sudo but I do not have privilege... but instead we can use `singularity` which can handle Docker container!!

--> Docker/singularity create a .sif file where all software and dependecies are installed (nothing related to a conda env).
To run the software; simply use `singularity run /DockerSoftareLocation/Docker.sif` and it will run the program!

```bash
module load singularity/3.11.1

cd /scr1/users/roulet/Akizu_Lab/Master/software

# Clone the repository if needed to examine files or use local resources
git clone https://github.com/jkimlab/msPIPE.git
cd msPIPE

# Build Singularity image from Docker Hub
singularity build mspipe.sif docker://jkimlab/mspipe:latest
# --> Fail no space left on device (space issue, not RAM related); re-direct temporary directory for installation 
export SINGULARITY_TMPDIR=/scr1/users/roulet/Akizu_Lab/Master/software/msPIPE

# Re-build, now using temporary directory for installation in /scr1
singularity build mspipe.sif docker://jkimlab/mspipe:latest

```

--> SUCCESS!!!

## Run msPIPE


```bash
srun --mem=500g --pty bash -l
cd /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021

singularity run /scr1/users/roulet/Akizu_Lab/Master/software/msPIPE/mspipe.sif

# Create meta file
nano scripts/params_docker.conf

# msPIPE docker vs singularity (! Keep `:ro`, means read only)
docker run -v /PATH/TO/INPUT/DATA:/msPIPE/data:ro -v /PATH/TO/REUSABLE/REFERENCE:/msPIPE/reference -v /PATH/TO/OUTDIR:/work_dir/ jkimlab/mspipe:latest msPIPE.py -p params_docker.conf -o result
singularity exec --bind /PATH/TO/INPUT/DATA:/msPIPE/data:ro --bind /PATH/TO/REUSABLE/REFERENCE:/msPIPE/reference --bind /PATH/TO/OUTDIR:/work_dir mspipe.sif msPIPE.py -p params_docker.conf -o result
## --> This is the same command in Docker vs singularity: Ask ChatGPT for tranlsation

# Run msPIPE
singularity exec \
  --bind /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/input:/scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/input:ro \
  --bind /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/output:/scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/output \
  /scr1/users/roulet/Akizu_Lab/Master/software/msPIPE/mspipe.sif \
  msPIPE.py -p /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/scripts/params_docker.conf -o /scr1/users/roulet/Akizu_Lab/008_ChIPseq_YAP_Conchi/004__WGBS_Dixon2021/output


```


--> Something weird with Rep2 from timGalore error at line 41188 lenght sequence and quality differ! rep1 worked great
----> Appeared that the R2 file is completely corrupted as I cannot `gunzip` it

--> File re-download in `input_R2Corr` and rename as previously; file in input replaced and msPIPE re run: 


```bash
sbatch scripts/msPIPE.sh # 18718424 fail; re-run with added --bind /scr1/users/roulet/Akizu_Lab/Master/meta:/msPIPE/reference 18748956 fail as bismark genome for CT conversion not prepared... Let's try to NOT mention fasta files to let the software prepare the genome 18757219 ok
```

--> The R2 is no more corrupted. BUT, Fail at bismark; cannot find genome file.
----> Need to specify in the command, where to look for the genome file (Even if it is indicated in `params_docker.conf` file)


--> Bug at bismark because bam file are not correctly sorted... FUCK msPIPE we cannot control each step; FUCK IT!!!



# Use pre-process file

- Collect in ENCODE other WGBS in H9 PSC/ESC. 

--> Found [here](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=WGBS&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=H1)and this [one](https://www.encodeproject.org/experiments/ENCSR744UGB/):  2 experiments with 2 bio rep; with Hg38 wig available 


Joe Ecker, Salk (ENCAN623MFB); ENCFF969PDB_R1 and ENCFF423PGW_R2
Richard Myers, HAIB (ENCAN508QJC); CpG sites coverage; ENCFF725YJG_R1 and ENCFF040LKO_R2
Joe Ecker, Salk (ENCSR744UGB); ENCFF895EIX_R1 and ENCFF459RZT_R2

- Download bigwig (`output/ENCODE`)
- Check m5C profile in regions that gain/lost EZH2 in WT vs YAPKO

--> Plot done in `008001` labnote `## EZH2 gain lost WT vs KO`


--> Ecker= (ENCAN623MFB and ENCSR744UGB) bad quality; Myers looks good






XXX HERE





# Quality control with FASTP (trim)

Run fastp
```bash
# run rep per rep
sbatch scripts/fastp_raw.sh # 18390922 ok
sbatch scripts/fastp_QSER1.sh # 18544351 ok
```
