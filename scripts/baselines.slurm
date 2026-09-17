#!/bin/bash
#SBATCH --job-name=baselines
#SBATCH --account=small_grains
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/baselines-%j.out
#
# Public healthy-tissue baselines (Makefile step 5c) as one unattended job: ENA download
# -> fastp -> HISAT2 (plant-only index) -> featureCounts -> 04_matrix/baseline_<host>.tsv.
#   sbatch scripts/baselines.slurm                # tier 1 (27 runs, ~105 GB fastq)
#   sbatch scripts/baselines.slurm BASE_TIER=2    # extra make variables pass straight through
# Phases run in order and `make -k` keeps going past a failed run; rerunning the job
# picks up where it stopped (downloads are md5-checked, so a partial one is redone).

set -u
# sbatch copies the script to a spool dir, so locate the repo via the submit dir
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}" || exit 1
mkdir -p logs

JOBS=4
THREADS=$(( SLURM_CPUS_PER_TASK / JOBS ))
MK="make -k -j$JOBS THREADS=$THREADS $*"

phase() {  # phase <name> <make targets...>
  local name=$1 rc; shift
  echo "=== $(date '+%F %T')  $name: make $*"
  $MK "$@"; rc=$?
  echo "=== $(date '+%F %T')  $name exit $rc"
}

# mkvar <VAR> [make vars...]: what the Makefile expands VAR to for this invocation
mkvar() {
  local probe="logs/.probe.$$.mk"
  printf 'include Makefile\nprint-var: ; @echo $(%s)\n' "$1" > "$probe"
  make -s -f "$probe" print-var "${@:2}"
  rm -f "$probe"
}
echo "node $(hostname)  cpus $SLURM_CPUS_PER_TASK  mem ${SLURM_MEM_PER_NODE:-?}M  args '$*'"
echo "runs: $(mkvar BASE_RUNS "$@")"

# the index build is one big single-threaded-ish step: give it all the cores
for h in $(mkvar BASE_REF "$@"); do
  phase "index $h" "references/plant_host/$h/$h.hisat2.ok" THREADS="$SLURM_CPUS_PER_TASK"
done
phase fetch   baseline-fetch
phase trim    baseline-trim
phase align   baseline-align
phase count   baseline-counts
phase matrix  baselines

echo
echo "=== baseline summary -> logs/baseline_summary.tsv"
{
  printf "run\thost_ref\ttissue\treplicate\tread_pairs\talignment_rate\tassigned_pairs\tassigned_frac\n"
  for c in 05_baseline/counts/*.txt.summary; do
    [ -f "$c" ] || continue
    r=$(basename "$c" .txt.summary)
    rate=$(grep -o '[0-9.]*% overall' "05_baseline/align/$r.hisat2.log" | cut -d% -f1)
    read -r asg tot <<< "$(awk 'NR > 1 {t += $2} $1 == "Assigned" {a = $2} END {print a, t}' "$c")"
    awk -F'\t' -v r="$r" -v rate="$rate" -v a="$asg" -v t="$tot" \
      '$1 == r {printf "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%.3f\n", r, $3, $5, $6, t, rate, a, (t ? a / t : 0)}' citrus_baselines.tsv
  done
} > logs/baseline_summary.tsv
column -t -s $'\t' logs/baseline_summary.tsv
echo
for m in 04_matrix/baseline_*.tsv; do
  [ -f "$m" ] && printf "%-40s %s genes x %s samples\n" "$m" "$(($(wc -l < "$m") - 1))" "$(($(head -1 "$m" | tr '\t' '\n' | wc -l) - 1))"
done
echo "=== done $(date '+%F %T')"
