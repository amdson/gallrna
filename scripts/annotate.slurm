#!/bin/bash
#SBATCH --job-name=annotate
#SBATCH --account=small_grains
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/annotate-%j.out
#
# Functional annotation of host genes (Makefile step 7) as a batch job: builds the
# eggNOG-mapper venv if it's missing, then runs `make annotate`. The first run also
# downloads the eggNOG database (~50 GB) into EGGNOG_DATA - keep that off home.
#   sbatch scripts/annotate.slurm                                   # citrus_sinensis
#   sbatch scripts/annotate.slurm "citrus_sinensis carica_papaya"   # other RefSeq hosts
#   EGGNOG_DATA=/90daydata/.../eggnog_data sbatch scripts/annotate.slurm

set -eu
# sbatch copies the script to a spool dir, so locate the repo via the submit dir
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}" || exit 1
mkdir -p logs

ANNOTATE="${1:-citrus_sinensis}"
echo "=== $(date '+%F %T')  node $(hostname)  cpus $SLURM_CPUS_PER_TASK  ANNOTATE='$ANNOTATE'"
echo "    EGGNOG_DATA=${EGGNOG_DATA:-<Makefile default>}"

[ -x "$HOME/.venvs/emapper/bin/emapper.py" ] || scripts/emapper-env.sh

make THREADS="$SLURM_CPUS_PER_TASK" ANNOTATE="$ANNOTATE" \
     ${EGGNOG_DATA:+EGGNOG_DATA="$EGGNOG_DATA"} annotate || {
  echo "=== $(date '+%F %T')  make annotate FAILED"
  for h in $ANNOTATE; do tail -n 30 "references/plant_host/$h/$h.emapper.log" 2>/dev/null; done
  exit 1
}

for h in $ANNOTATE; do
  t="references/plant_host/$h/$h.genes.tsv"
  echo "=== $t: $(($(wc -l < "$t") - 1)) genes"
  awk -F'\t' 'NR > 1 && $2 == "protein_coding" {n++; if ($12 != "") e++; if ($7 != "") g++;
                if ($13 != "") a++; if ($19 == 1) r++}
              END {printf "    protein-coding %d: eggNOG %d, GO %d, Arabidopsis hit %d, reciprocal best %d\n",
                   n, e, g, a, r}' "$t"
done
echo "=== done $(date '+%F %T')"
