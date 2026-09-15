#!/bin/bash
# How often does HISAT2 report a spliced alignment (N in CIGAR) on each organism?
#
# HISAT2 is a spliced aligner and the pipeline runs it with --dta over a combined
# plant+bacterium index, so the bacterial contigs get spliced alignment too even
# though bacterial transcripts have no introns. This measures the cost.
#
# Bacterial contigs are counted in full (they are tiny). The plant side is capped
# at the first PLANT_CAP alignments of one chromosome - enough for a rate, and it
# keeps the whole thing to a couple of minutes instead of an hour.
#
#   scripts/splice_rates.sh > logs/splice_rates.tsv
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
module load ${MOD_SAMTOOLS:-samtools/1.17} 2>/dev/null

PLANT_CAP=${PLANT_CAP:-2000000}
printf "sample\tagro_alignments\tagro_spliced\tplant_alignments\tplant_spliced\n"
for bam in 02_align/*.bam; do
    s=$(basename "$bam" .bam)
    mapfile -t agro < <(samtools idxstats "$bam" | awk '$1 ~ /^agro_/ {print $1}')
    [ ${#agro[@]} -eq 0 ] && continue
    read -r an ak < <(samtools view -F 0x900 "$bam" "${agro[@]}" |
        awk '{n++; if ($6 ~ /N/) k++} END {print n+0, k+0}')
    # busiest plant contig, so the cap lands on real data
    pc=$(samtools idxstats "$bam" | awk '$1 !~ /^agro_/ && $1 != "*"' | sort -k3,3nr | head -1 | cut -f1)
    read -r pn pk < <(samtools view -F 0x900 "$bam" "$pc" | head -n "$PLANT_CAP" |
        awk '{n++; if ($6 ~ /N/) k++} END {print n+0, k+0}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$an" "$ak" "$pn" "$pk"
done
