#!/bin/bash
#SBATCH --job-name=dotplot
#SBATCH --account=small_grains
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/dotplot-%j.out
#
# Whole-genome alignment of two host assemblies for a dotplot. Default: the two parents of
# Carrizo citrange (README.md section 6.2), Poncirus ZK8 (query, y axis) against sweet orange
# DVS_A1.0 (target, x axis; the reference the CC galls are aligned to).
#   sbatch scripts/dotplot.slurm                               # citrus_sinensis poncirus_trifoliata
#   sbatch scripts/dotplot.slurm <target_host> <query_host>
# minimap2 -x asm20 (cross-species preset, up to ~5% divergence) with -c, so each block
# carries a base-level divergence (de:f). Output is a slim PAF: the 12 standard columns plus
# the tp (P primary / S secondary) and de tags, without the CIGAR strings:
#   references/plant_host/<query>/<target>__<query>.paf
# and the differences between the two assemblies from `paftools.js call` (ships with minimap2):
# slimmed to five columns: R|V, target contig, start, end, indel length (0 = substitution).
# R = target region covered by exactly one alignment, V = one difference inside such a region
#   references/plant_host/<query>/<target>__<query>.var.tsv
# Plot both with notebooks/04_parent_dotplot.ipynb.

set -euo pipefail
# sbatch copies the script to a spool dir, so locate the repo via the submit dir
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}" || exit 1
mkdir -p logs

TARGET="${1:-citrus_sinensis}"
QUERY="${2:-poncirus_trifoliata}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MOD_MINIMAP2="${MOD_MINIMAP2:-minimap2/2.24}"
tfa="references/plant_host/$TARGET/$TARGET.fasta"
qfa="references/plant_host/$QUERY/$QUERY.fasta"
out="references/plant_host/$QUERY/${TARGET}__${QUERY}.paf"

echo "=== $(date '+%F %T')  node $(hostname)  cpus $THREADS  $QUERY (query) vs $TARGET (target)"
module load "$MOD_MINIMAP2"
minimap2 --version

full="$out.full.tmp"          # with CIGAR and cs strings: large, removed at the end
minimap2 -x asm20 -c --cs -t "$THREADS" "$tfa" "$qfa" 2> "logs/dotplot-${TARGET}__${QUERY}.minimap2.log" > "$full"

awk 'BEGIN {OFS = "\t"}
         {tp = de = "";
          for (i = 13; i <= NF; i++) {if ($i ~ /^tp:A:/) tp = $i; else if ($i ~ /^de:f:/) de = $i}
          print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, tp, de}' "$full" > "$out.tmp"
mv "$out.tmp" "$out"

# variants from alignments >= 2 kb with mapq >= 5 (paftools defaults are 50 kb, for same-species assemblies)
var="${out%.paf}.var.tsv"
sort -k6,6 -k8,8n "$full" | paftools.js call -l 2000 -L 2000 - 2> "logs/dotplot-${TARGET}__${QUERY}.paftools.log" \
  | awk -F'\t' -v OFS='\t' '$1 == "R" {print "R", $2, $3, $4, 0}
         $1 == "V" {l = ($7 == "-") ? length($8) : ($8 == "-") ? length($7) : 0; print "V", $2, $3, $4, l}' > "$var.tmp"
mv "$var.tmp" "$var"
rm -f "$full"

echo "=== $out: $(wc -l < "$out") alignment blocks"
awk -F'\t' '$13 == "tp:A:P" {n++; bp += $11} END {printf "    primary blocks %d, %.1f Mb of alignment\n", n, bp / 1e6}' "$out"
awk '$1 == "R" {r += $4 - $3} $1 == "V" {v++} END {printf "=== %s: %.1f Mb covered once, %d variants\n", f, r / 1e6, v}' f="$var" "$var"
echo "=== done $(date '+%F %T')"
