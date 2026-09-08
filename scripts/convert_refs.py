#!/usr/bin/env python
"""Convert SnapGene (.dna) / GenBank (.gb) replicon files into one genome FASTA + GFF3.

Usage: convert_refs.py <out_prefix> <contig_id>=<infile> [<contig_id>=<infile> ...]

Writes <out_prefix>.fasta and <out_prefix>.gff3 with contigs named exactly as given
(use the pipeline's organism prefixes, e.g. agro_C58_circ).
"""
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load(path):
    fmt = "snapgene" if path.lower().endswith(".dna") else "genbank"
    recs = list(SeqIO.parse(path, fmt))
    if len(recs) != 1:
        sys.exit(f"{path}: expected 1 record, got {len(recs)}")
    return recs[0]

def main():
    out_prefix, pairs = sys.argv[1], sys.argv[2:]
    records = []
    for pair in pairs:
        contig_id, path = pair.split("=", 1)
        r = load(path)
        records.append(SeqRecord(r.seq, id=contig_id, description=f"source={path}", features=r.features))

    SeqIO.write(records, f"{out_prefix}.fasta", "fasta")

    with open(f"{out_prefix}.gff3", "w") as gff:
        gff.write("##gff-version 3\n")
        n = 0
        for rec in records:
            gff.write(f"##sequence-region {rec.id} 1 {len(rec.seq)}\n")
            for feat in rec.features:
                if feat.type not in ("gene", "CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"):
                    continue
                q = feat.qualifiers
                name = (q.get("locus_tag") or q.get("gene") or q.get("label") or [f"feat{n}"])[0]
                attrs = f"ID={rec.id}.{feat.type}.{n};Name={name}"
                if "product" in q:
                    prod = q["product"][0].replace(";", ",").replace("=", " ")
                    attrs += f";product={prod}"
                strand = "+" if feat.location.strand != -1 else "-"
                gff.write(f"{rec.id}\tconvert_refs\t{feat.type}\t"
                          f"{int(feat.location.start) + 1}\t{int(feat.location.end)}\t.\t{strand}\t"
                          f"{'0' if feat.type == 'CDS' else '.'}\t{attrs}\n")
                n += 1
    counts = ", ".join(f"{r.id}:{len(r.seq)}bp" for r in records)
    print(f"wrote {out_prefix}.fasta / .gff3 ({n} features; {counts})")

if __name__ == "__main__":
    main()
