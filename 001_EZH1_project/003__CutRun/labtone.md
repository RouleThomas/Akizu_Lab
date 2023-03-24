# pipeline used from previous analyses


# Import files from Google drive to the cluster

Using `cp`; check integrity (volume) of files!! (`ls -lh`)

# File renaiming

**Move all fastq files within the input folder** (now [file] are in input/folder1/[file]) using a loop:
```bash
for file in input/*/;
    do mv "$file"* input/;
done

rmdir input/*
```
- `input/*/` this go 1 folder downstream; whatever their names
- we move file to `input/` directory 

**Rename files**
Detail about file renaiming can be found in `CutRun_infos.xlsx` in my Google drive folder.



XXX It fail let's modify _1 script wich is cleaner!



```bash
# example for 1 file:
outdir="input"

x="R1_Het5_Igg_1"
raw_f="R1_Het5_Igg_CKDL220021572-1A_HF3THDSX5_L1_1.fq.gz"
new_f="${outdir}/${x}.fq.gz"

if [[ -f "$raw_f" && ! -f "$new_f" ]]; then
	mv "$raw_f" "$new_f"
elif [[ ! -f "$raw_f" && ! -f "$new_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

# Run slurm for all:
sbatch rename_raw.sh # 11371503, not all complete; no idea why
```
--> For an unknown reason it does not want to rename some file lol! Weird as shit! I did it manually with `mv` for the remaining file...


