#!/bin/bash
# Download one paired-end run from ENA and verify the md5 sums it publishes.
#   scripts/ena_fetch.sh <run_accession> <outdir>
# writes <outdir>/<run>_1.fastq.gz and <outdir>/<run>_2.fastq.gz. Uses the ENA
# portal filereport API (fastq_ftp, fastq_md5), so sra-tools is not needed. Runs
# that also carry an unpaired file (RUN.fastq.gz) have it ignored.
set -euo pipefail
run=$1; out=$2
mkdir -p "$out"

api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
rec=$(curl -fsS --retry 5 --retry-delay 10 "$api" | tail -n +2 | head -n 1)
[ -n "$rec" ] || { echo "ena_fetch: no ENA record for $run" >&2; exit 1; }
# the API always puts run_accession first, before the requested fields
IFS=$'\t' read -r acc ftp md5 <<< "$rec"
[ "$acc" = "$run" ] || { echo "ena_fetch: ENA returned $acc for $run" >&2; exit 1; }
IFS=';' read -ra urls <<< "$ftp"
IFS=';' read -ra sums <<< "$md5"

for mate in 1 2; do
  url=""; sum=""
  for i in "${!urls[@]}"; do
    case "${urls[$i]}" in *_${mate}.fastq.gz) url=${urls[$i]}; sum=${sums[$i]};; esac
  done
  [ -n "$url" ] || { echo "ena_fetch: $run has no _${mate}.fastq.gz at ENA: $ftp" >&2; exit 1; }
  f="$out/${run}_${mate}.fastq.gz"
  if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$sum" ]; then continue; fi
  ok=0
  for attempt in 1 2 3; do
    # -C - resumes a partial file; a corrupt one is removed and refetched
    curl -fsSL --retry 5 --retry-delay 30 -C - -o "$f.part" "https://$url" || true
    if [ -s "$f.part" ] && [ "$(md5sum "$f.part" | cut -d' ' -f1)" = "$sum" ]; then
      mv "$f.part" "$f"; ok=1; break
    fi
    echo "ena_fetch: $f attempt $attempt failed md5 check, retrying" >&2
    rm -f "$f.part"
  done
  [ $ok = 1 ] || { echo "ena_fetch: giving up on $f" >&2; exit 1; }
done
echo "ena_fetch: $run ok"
