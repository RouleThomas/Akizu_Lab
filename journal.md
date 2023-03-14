# Daily-To-Do-List
--> On Benchling

# Infos
## Usefull commands

- If a window is buggy, stuck; run: `killall -3 gnome-shell`

## File transfer
Only on the CHOP cluster or on a CHOP machine...
Just copy/paste with cp from cluster to local computer Google Drive!

## RStudio with ressource
- Go to https://respublica-web.research.chop.edu and interactive app --> Rstudio.
- Here I have access to my `/home/roulet/`

## IGV with ressource
- Maybe try Mobaxterm with ssh and -X
XXX ??


## Cluster-files folder-like interace
No idea! Dolphin I can see files but not transfer them...\
To allow cp/paste between Vitrual machine and computer: 


## Directories
- BIOINFORMATICS SPACE: Working directory in the cluster is `/scr1/users/roulet` 30TB limit (temporary space, belong to CHOP)-
- SPACE TO SEE FILE (folder-like env): `/home/roulet/`
- LOCAL COMPUTER FILE:\
**Box**:`/home/roulet/tsclient/roule/Box`\
**Google drive**: `/home/roulet/tsclient/roule/Google Drive Streaming`



## Run job

- Interactive `srun --mem=20g --pty bash -l`.
- Sbatch `sbatch job.sh` (edit script in VSC, then create it on the cluster with `touch script.sh` edit it with `nano script.sh` copy paste from VSC, and run it)
- List jobs: `squeue -u roulet` (`scancel [JOBID]` to cancel)


If encounteer `bash __vte_prompt_command command not found` error message. Do the following:
1. add this add the end of the ~/.bashrc file:
```bash
__vte_prompt_command() {
  local pwdmaxlen=30
  local pwdoffset=$(( ${#PWD} - pwdmaxlen ))
  local pwdir=${PWD:-$HOME}

  [ "${pwdoffset}" -gt "0" ] && pwdir="…${pwdir:${pwdoffset}}"
  printf "\033]0;%s@%s:%s\007" "${USER}" "${HOSTNAME%%.*}" "${pwdir}"
}
```
2. then `source ~/.bashrc`

