#!/bin/bash
# How much of the library is intronic? A polyA/mRNA library is almost purely
# exonic; a ribo-depleted total-RNA library retains pre-mRNA and typically runs
# 20-40% intronic. This is the read-position test for which prep was used - it
# is independent of gene-class arguments (rRNA share and the like).
#
# The pipeline counts with `-t gene`, so introns sit INSIDE the counted feature
# and Unassigned_NoFeatures is intergenic only. Re-counting the same BAM with
# `-t exon` moves intronic reads into NoFeatures, so
#     intronic = NoFeatures(exon) - NoFeatures(gene)
#
# Only hosts with gene models can be measured (citrus, tomato, papaya).
#
#   scripts/intron_fraction.sh > logs/intron_fraction.tsv
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
module load ${MOD_SUBREAD:-subread/2.0.4} 2>/dev/null

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
ANNOT=$(grep -oP '^ANNOT_HOSTS\s*:?=\s*\K.*' Makefile)

printf "sample\thost\tassigned\tintergenic\tintronic\tlocated\n"
for bam in 02_align/*.bam; do
    s=$(basename "$bam" .bam)
    host=$(grep -oP "^HOST_$s\s*=\s*\K\S+" Makefile) || continue
    [[ " $ANNOT " == *" $host "* ]] || continue          # no gene models, nothing to measure
    strain=$(grep -oP '^[0-9]+' <<<"$s")
    combo="strain_${strain}__${host}"
    exgff="$TMP/$combo.exon.gff3"

    if [ ! -s "$exgff" ]; then
        # some organellar exons carry locus_tag but no gene=; featureCounts aborts
        # on the first such record, so fill the attribute in rather than drop them
        awk -F'\t' -v OFS='\t' '$0!~/^#/ && $3=="exon" {
            if ($9 !~ /(^|;)gene=/) {
                lt="unknown"
                if (match($9, /locus_tag=[^;]+/)) lt=substr($9, RSTART+10, RLENGTH-10)
                $9 = $9 ";gene=" lt }
            print }' "references/combined/$combo.gff3" > "$exgff"
    fi

    featureCounts -p --countReadPairs -T "${THREADS:-2}" -s 0 -t exon -g gene \
        -a "$exgff" -o "$TMP/$s.exon.txt" "$bam" > "$TMP/$s.log" 2>&1
    exon_nf=$(awk '$1=="Unassigned_NoFeatures"{print $2}' "$TMP/$s.exon.txt.summary")
    gene_nf=$(awk '$1=="Unassigned_NoFeatures"{print $2}' "03_counts/$s.txt.summary")
    assigned=$(awk '$1=="Assigned"{print $2}' "03_counts/$s.txt.summary")
    [ -z "${exon_nf:-}" ] && { echo "$s: featureCounts failed, see $TMP/$s.log" >&2; continue; }
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$host" "$assigned" "$gene_nf" \
        "$((exon_nf - gene_nf))" "$((assigned + gene_nf))"
done
