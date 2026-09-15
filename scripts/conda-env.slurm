#!/bin/bash
#SBATCH --job-name=rnaseq-env
#SBATCH --account=small_grains
#SBATCH --partition=ceres
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=logs/conda-env-%j.out
#
# (Re)build the `rnaseq` conda env unattended. Linking ~120k files into NFS
# home is slow (>1 h) but needs almost no CPU, so this runs on 2 cores.
# Skips the rebuild if the env already looks complete.
#   sbatch scripts/conda-env.slurm
#   sbatch --dependency=afterany:<jobid> scripts/conda-env.slurm   # wait for a job first

set -u
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}" || exit 1
module load miniconda/24.7.1-2
ENV=$HOME/.conda/envs/rnaseq

if [ -x "$ENV/bin/fastp" ] && [ -x "$ENV/bin/multiqc" ] && [ -x "$ENV/bin/R" ] && [ -s "$ENV/conda-meta/history" ]; then
  echo "rnaseq env already complete at $ENV - nothing to do"; exit 0
fi
if [ -d "$ENV" ]; then
  echo "removing incomplete env $ENV"; rm -rf "$ENV"
fi
echo "=== $(date '+%F %T') mamba env create"
time mamba env create -y -f environment.yml
rc=$?
echo "=== $(date '+%F %T') exit $rc"
ls "$ENV/bin" | wc -l
"$ENV/bin/fastp" --version 2>&1; "$ENV/bin/multiqc" --version 2>&1; "$ENV/bin/python" -c 'import Bio; print("biopython", Bio.__version__)'
exit $rc
