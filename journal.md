# Daily-To-Do-List
--> On Benchling

# Infos
**Directories**
- BIOINFORMATICS SPACE: Working directory in the cluster is `/scr1/users/roulet` 30TB limit (temporary space, belong to CHOP)-
- SPACE TO SEE FILE (folder-like env): `/home/roulet/`
- LOCAL COMPUTER FILE: `/home/roulet/tsclient/roule/` is my C: computer, with `Box` and `GoogleDrive` (to see Google Drive/Preferences/Settings/ Stream to folder)

**Run job**

- Interactive `srun --mem=20g hostname`.

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
- Sbatch `sbatch job.sh`