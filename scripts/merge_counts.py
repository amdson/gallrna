#!/usr/bin/env python
"""Join per-sample featureCounts outputs into one matrix, filtered by contig prefix.

Usage: merge_counts.py <out.tsv> <keep|drop> <contig_prefix> <counts1.txt> [counts2.txt ...]

keep agro_  -> rows on contigs starting with the prefix (bacterial genes)
drop agro_  -> rows on all other contigs (plant genes)
Sample names are taken from the counts filenames (basename minus .txt).
"""
import os
import sys

def main():
    out, mode, prefix, files = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    matrix, order = {}, []
    samples = [os.path.basename(f).removesuffix(".txt") for f in files]
    for f in files:
        with open(f) as fh:
            for line in fh:
                if line.startswith(("#", "Geneid")):
                    continue
                cols = line.rstrip("\n").split("\t")
                gene, chrom, count = cols[0], cols[1].split(";")[0], cols[-1]
                hit = chrom.startswith(prefix)
                if (mode == "keep") != hit:
                    continue
                if gene not in matrix:
                    matrix[gene] = {}
                    order.append(gene)
                matrix[gene][f] = count
    with open(out, "w") as o:
        o.write("gene\t" + "\t".join(samples) + "\n")
        for g in order:
            o.write(g + "\t" + "\t".join(matrix[g].get(f, "0") for f in files) + "\n")
    print(f"{out}: {len(order)} genes x {len(samples)} samples")

if __name__ == "__main__":
    main()
