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
            # SnapGene/RAST GenBank files carry no `gene` features, only CDS/tRNA/rRNA.
            # featureCounts counts `-t gene`, so synthesize one gene row per feature
            # (plus the original row as its child) unless real gene features exist.
            has_gene = any(f.type == "gene" for f in rec.features)
            for feat in rec.features:
                if feat.type not in ("gene", "CDS", "tRNA", "rRNA", "ncRNA", "tmRNA"):
                    continue
                q = feat.qualifiers
                seed = [x.split(":", 1)[1] for x in q.get("db_xref", []) if x.startswith("SEED:")]
                name = (q.get("locus_tag") or q.get("gene") or q.get("label") or seed or [f"feat{n}"])[0]
                strand = "+" if feat.location.strand != -1 else "-"
                extra = ""
                if "product" in q:
                    prod = q["product"][0].replace(";", ",").replace("=", " ")
                    extra = f";product={prod}"
                gid = f"{rec.id}.gene.{n}"
                # Circular replicons: a feature crossing the origin is a compound location
                # whose outer bounds span the whole replicon. Write one row per part (same
                # ID - featureCounts merges rows sharing an ID into one meta-feature).
                parts = [(int(p.start) + 1, int(p.end)) for p in feat.location.parts]
                if len(parts) > 1:
                    print(f"note: {rec.id} {feat.type} {name!r} written as {len(parts)} parts")
                for start, end in parts:
                    if feat.type == "gene" or not has_gene:
                        gff.write(f"{rec.id}\tconvert_refs\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
                                  f"ID={gid};Name={name};feature_type={feat.type}{extra}\n")
                    if feat.type != "gene":
                        gff.write(f"{rec.id}\tconvert_refs\t{feat.type}\t{start}\t{end}\t.\t{strand}\t"
                                  f"{'0' if feat.type == 'CDS' else '.'}\t"
                                  f"ID={rec.id}.{feat.type}.{n};Parent={gid};Name={name}{extra}\n")
                n += 1
    counts = ", ".join(f"{r.id}:{len(r.seq)}bp" for r in records)
    print(f"wrote {out_prefix}.fasta / .gff3 ({n} features; {counts})")

if __name__ == "__main__":
    main()
