#!/usr/bin/env python3
"""Join per-gene functional annotation into one table keyed on count-matrix gene IDs.

  gene_annotation.py <host.gff3> <host.emapper.annotations> <host.ath.tsv> \
                     <host.ath_rev.tsv> <arabidopsis.agi.tsv> <out.genes.tsv>

One row per `gene` feature in the GFF, i.e. every feature featureCounts counts with
-t gene -g ID. Non-coding genes are included with empty protein columns.
  refseq_description  the GFF description= attribute
  em_*                eggNOG-mapper: preferred name, description, COG category, GO, EC,
                      KEGG KO/pathway, Pfam, and the taxonomic level annotations came from
  ath_*               DIAMOND best hit in Arabidopsis (Ensembl Plants TAIR10/Araport11
                      proteome); ath_qcov = aligned length / query length; ath_rbh = 1 when
                      that Arabidopsis gene's best hit back in this genome is the same gene
"""
import sys
from urllib.parse import unquote

EMAPPER_COLS = [("Preferred_name", "em_preferred_name"), ("Description", "em_description"),
                ("COG_category", "em_cog_category"), ("GOs", "em_go"), ("EC", "em_ec"),
                ("KEGG_ko", "em_kegg_ko"), ("KEGG_Pathway", "em_kegg_pathway"),
                ("PFAMs", "em_pfam"), ("max_annot_lvl", "em_annot_level")]
ATH_COLS = ["ath_agi", "ath_symbol", "ath_description", "ath_pident", "ath_qcov",
            "ath_evalue", "ath_rbh"]


def attrs(col9):
    return dict(kv.split("=", 1) for kv in col9.split(";") if "=" in kv)


def best_hits(path):
    """First hit per query in DIAMOND outfmt 6 (qseqid sseqid pident length qlen slen
    evalue bitscore); DIAMOND reports a query's hits best first."""
    hits = {}
    with open(path) as fh:
        for line in fh:
            q, s, pident, length, qlen, _slen, evalue, _bits = line.rstrip("\n").split("\t")
            if q not in hits:
                hits[q] = (s, pident, f"{int(length) / int(qlen):.2f}", evalue)
    return hits


def main(gff, emapper, fwd, rev, agi, out):
    genes, seen = [], set()
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] == "gene":
                a = attrs(f[8])
                if a["ID"] not in seen:
                    seen.add(a["ID"])
                    genes.append((a["ID"], a.get("gene_biotype", ""),
                                  unquote(a.get("description", ""))))

    em, cols = {}, None
    with open(emapper) as fh:
        for line in fh:
            if line.startswith("#query"):
                cols = line[1:].rstrip("\n").split("\t")
            elif not line.startswith("#") and cols:
                row = dict(zip(cols, line.rstrip("\n").split("\t")))
                em[row["query"]] = [("" if row.get(old, "-") == "-" else row[old])
                                    for old, _ in EMAPPER_COLS]
    if cols is None:
        sys.exit(f"{emapper}: no '#query' header line - not an eggNOG-mapper annotations file?")

    names = {}
    with open(agi) as fh:
        next(fh)
        for line in fh:
            gene, symbol, desc = line.rstrip("\n").split("\t")
            names[gene] = (symbol, desc)
    fwd_hits, rev_hits = best_hits(fwd), best_hits(rev)

    n_em = n_go = n_ath = n_rbh = 0
    with open(out, "w") as fh:
        fh.write("\t".join(["gene_id", "gene_biotype", "refseq_description"]
                           + [new for _, new in EMAPPER_COLS] + ATH_COLS) + "\n")
        for gene_id, biotype, desc in genes:
            e = em.get(gene_id, [""] * len(EMAPPER_COLS))
            n_em += gene_id in em
            n_go += bool(e[3])
            ath = [""] * len(ATH_COLS)
            if gene_id in fwd_hits:
                s, pident, qcov, evalue = fwd_hits[gene_id]
                rbh = rev_hits.get(s, ("",))[0] == gene_id
                ath = [s, *names.get(s, ("", "")), pident, qcov, evalue, str(int(rbh))]
                n_ath += 1
                n_rbh += rbh
            fh.write("\t".join([gene_id, biotype, desc] + e + ath) + "\n")
    print(f"{out}: {len(genes)} genes; eggNOG-annotated {n_em} (with GO {n_go}); "
          f"Arabidopsis hit {n_ath} (reciprocal best {n_rbh})", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    main(*sys.argv[1:])
