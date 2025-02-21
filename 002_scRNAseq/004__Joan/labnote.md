Provide Pipseeker output to Joan.

- install Pipseeker 2.1.4 in the cluster
- download scRNAseq data (see Slack for code) from Novogene
- Run pipseeker


**--> TO MAKE SOME SPACE, ALL FILES AS BEEN DELETED, EXCEPT SCRIPTS and SLURMS**



# Pipseeker installation and novogene data download

Pipseeker need to be installed in unix and [previous versions](https://www.fluentbio.com/pipseeker-release-notes-archives/)

```bash
# install Pipeseeker
tar -zxvf *.gz

# --> move executbale in /Master/software
# activate executballe
chmod +x ~/PIPseeker/pipseeker

# To run it:
/scr1/users/roulet/Akizu_Lab/Master/software/pipseeker -h

# download data
cd input/
nano url.txt # copy all URL of donwload with "export link"
wget -i url.txt


```

--> Pipseeker succesfully installed

--> Data succesfully downloaded (fastq only)

*NOTE: I had to unzip manually on Windows the genome STAR file; it failed using `unzip *.zip`*; the file seems corrupted, trying to clean it with:  zip -FF broken.zip --out fixed.zip --> So I downloaded new mm10 genome from pipseeker [website](https://www.fluentbio.com/resources/pipseeker-downloads/) `tar -xvzf *.tar.gz`

# Run pipseeker


```bash
sbatch scripts/pipseeker.sh # 11850273 cancelled too long with 1 thread
sbatch scripts/pipseeker_quick.sh # 11863225 ok

sbatch scripts/pipseeker_GFP.sh # 12007387 
```

--> The `quick` version is over in less than 12hrs, the other is still running after 22hrs!

--> The custom genome is needed as it allow them to track GFP expressing cells.
----> It now work, the issue was that the .zip and unzipping; I had download file per file to avoid compression and it work!

--> I did not use in output_EGFP/ their `-min--sensitivity` and `-max--sensitivty` parameters as do not exist in my version!




