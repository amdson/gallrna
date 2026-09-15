#!/bin/bash
#SBATCH --job-name=gallrna
#SBATCH --account=small_grains
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/pipeline-%j.out
#
# Unattended run of the Makefile pipeline (steps 0-5) for all samples.
#   sbatch scripts/pipeline.slurm                 # everything
#   sbatch scripts/pipeline.slurm citrus_sinensis # HOSTSEL subset
# Phases run in order so the cheap, high-value outputs (mapping fractions)
# land first; `make -k` keeps going past a failed sample. CONDA points at
# ~/.local/tools, which wraps fastp/python/multiqc so the job does not depend
# on the rnaseq conda env having finished installing.

set -u
# sbatch copies the script to a spool dir, so locate the repo via the submit dir
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}" || exit 1
mkdir -p logs

HOSTSEL="${1:-}"
JOBS=4
THREADS=$(( SLURM_CPUS_PER_TASK / JOBS ))
MK="make -k -j$JOBS THREADS=$THREADS CONDA=$HOME/.local/tools HOSTSEL=$HOSTSEL"

phase() {  # phase <name> <make targets...>
  local name=$1; shift
  echo "=== $(date '+%F %T')  $name: make $*"
  $MK "$@"
  echo "=== $(date '+%F %T')  $name exit $?"
}

echo "node $(hostname)  cpus $SLURM_CPUS_PER_TASK  mem ${SLURM_MEM_PER_NODE:-?}M  HOSTSEL='${HOSTSEL}'"
phase refs      refs
phase indexes   combined
phase trim      trim
phase align     fractions        # aligns everything, then the agro/plant read table
phase count     counts
phase matrices  matrices
phase qc        qc               # FastQC on all raw files + MultiQC (needs conda env)

echo
echo "=== mapping fractions"
cat logs/mapping_fractions.tsv 2>/dev/null
echo
echo "=== hisat2 overall alignment rates"
for f in 02_align/*.hisat2.log; do printf "%-16s %s\n" "$(basename "$f" .hisat2.log)" "$(grep overall "$f")"; done
echo
echo "=== count matrices"
for m in 04_matrix/*.tsv; do printf "%-40s %s genes x %s samples\n" "$m" "$(($(wc -l < "$m") - 1))" "$(($(head -1 "$m" | tr '\t' '\n' | wc -l) - 1))"; done
echo "=== done $(date '+%F %T')"
