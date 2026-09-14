#!/usr/bin/env python3
"""One protein per gene (the longest isoform), named by gene.

  longest_proteins.py refseq  <genomic.gff3> <protein.faa> <out.faa>
      NCBI RefSeq: protein_id on each CDS -> Parent mRNA -> Parent gene. Headers become
      the gene feature ID (e.g. gene-LOC102606664), the key featureCounts -g ID writes,
      so annotation joins straight onto the count matrices.
  longest_proteins.py ensembl <pep.all.fa[.gz]> <out.faa> <out.agi.tsv>
      Ensembl Plants pep headers carry gene:, gene_symbol: and description:. Headers
      become the gene (e.g. AT1G01010) and a gene/symbol/description table is written.
"""
import gzip
import re
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    header, seq = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def attrs(col9):
    return dict(kv.split("=", 1) for kv in col9.rstrip("\n").split(";") if "=" in kv)


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def refseq(gff, faa, out):
    parent, prot_parent = {}, {}
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            if f[2] == "CDS":
                if "protein_id" in a and "Parent" in a:
                    prot_parent[a["protein_id"]] = a["Parent"].split(",")[0]
            elif "ID" in a and "Parent" in a:
                parent[a["ID"]] = a["Parent"].split(",")[0]
    best, n, unplaced = {}, 0, 0
    for header, seq in read_fasta(faa):
        n += 1
        node = prot_parent.get(header.split()[0])
        if node is None:
            unplaced += 1
            continue
        while node in parent:          # mRNA -> gene (gene rows have no Parent)
            node = parent[node]
        if len(seq) > len(best.get(node, "")):
            best[node] = seq
    write_fasta(out, best)
    print(f"{faa}: {n} proteins -> {len(best)} genes (longest isoform each); "
          f"{unplaced} proteins had no CDS in the GFF", file=sys.stderr)


def ensembl(pep, out, table):
    best, info = {}, {}
    for header, seq in read_fasta(pep):
        gene = re.search(r" gene:(\S+)", header)
        if not gene:
            continue
        gene = gene.group(1)
        symbol = re.search(r" gene_symbol:(\S+)", header)
        desc = header.split(" description:", 1)[1] if " description:" in header else ""
        desc = re.sub(r"\s*\[Source:[^\]]*\]$", "", desc)
        info[gene] = (symbol.group(1) if symbol else "", desc)
        if len(seq) > len(best.get(gene, "")):
            best[gene] = seq
    write_fasta(out, best)
    with open(table, "w") as fh:
        fh.write("agi\tsymbol\tdescription\n")
        for gene in best:
            fh.write(f"{gene}\t{info[gene][0]}\t{info[gene][1]}\n")
    print(f"{pep}: {len(best)} genes (longest isoform each)", file=sys.stderr)


if __name__ == "__main__":
    mode, args = sys.argv[1], sys.argv[2:]
    if mode == "refseq" and len(args) == 3:
        refseq(*args)
    elif mode == "ensembl" and len(args) == 3:
        ensembl(*args)
    else:
        sys.exit(__doc__)
