#!/bin/bash

# Path to signature file (1st row contains the signature names)
SIGFILE="/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/signatures/human_sbs96_filtered_v1_0.txt"

# Where to write all generated scripts
OUTDIR="scripts/submit_all_signatures_experimental"
mkdir -p "$OUTDIR"

# Common SLURM config
COMMON_SLURM="#SBATCH --array=0-79
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=01:00:00"

# Output and error directory
LOGDIR="/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/logs"

# Mutation count list and replicates
N_LIST="(500 1000 2000 4000 8000 16000 32000 64000)"
REPS=10

# Extract signature names from the first row (after 'Type')
SIGS=$(head -n1 "$SIGFILE" | cut -f2-)

for SIG in $SIGS; do
  SCRIPT="$OUTDIR/submit_simulation_array_experimental-${SIG}.sh"
  cat > "$SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=sim_${SIG}
$COMMON_SLURM
#SBATCH --output=${LOGDIR}/${SIG}_%A_%a.out
#SBATCH --error=${LOGDIR}/${SIG}_%A_%a.err

# Parameters
SIGNATURE="${SIG}"
N_LIST=${N_LIST}
REPS=${REPS}

# Map SLURM_ARRAY_TASK_ID to mutation count and replicate
IDX=\$SLURM_ARRAY_TASK_ID
N_IDX=\$((IDX / REPS))
REP=\$((IDX % REPS + 1))
MUTATION_COUNT=\${N_LIST[\$N_IDX]}

# Run the job from scripts/ directory context
python scripts/simulate_array_experimental.py \\
  --signature \$SIGNATURE \\
  --n \$MUTATION_COUNT \\
  --rep \$REP \\
  --seed \$IDX
EOF

  chmod +x "$SCRIPT"
done

echo "✅ Generated SLURM scripts for all signatures in: $OUTDIR"
