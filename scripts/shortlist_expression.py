#!/usr/bin/env python3
"""Shortlist method 1 (CITRUS_HOST_PLAN.md section 4 F): rank host genes by expression in
the galls, without any baseline.

For every gene: TPM per gall (featureCounts pair counts / gene length), then the gene's
within-gall percentile among the kept genes (1.0 = most expressed). A group's score is the
gene's LOWEST percentile across the group's galls, so a gene ranks high only if it is highly
expressed in every gall of the group; ties break on mean TPM. Kept genes: protein-coding and
nuclear (organellar contigs dropped). `n_same_ath_hit` counts how many citrus genes share the
same best Arabidopsis hit - a rough warning for multi-copy families whose promoter may be
hard to clone unambiguously.

Usage: shortlist_expression.py <plant matrix> <one featureCounts file (for lengths/contigs)>
                               <genes.tsv> <full table out> <top20 dir>
                               <group>=<sample>[,<sample>...] ...
"""
import csv
import os
import sys

ORGANELLE = {"NC_008334.1", "NC_037463.1"}  # C. sinensis chloroplast, mitochondrion
TOP = 20
ANNOT = ["refseq_description", "em_preferred_name", "ath_symbol", "ath_agi", "ath_rbh",
         "ath_pident", "em_cog_category", "em_pfam"]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def percentiles(values):
    """Average-rank percentile in (0, 1] for a list of numbers, ties averaged."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    pct = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean_rank = (i + j) / 2 + 1
        for k in range(i, j + 1):
            pct[order[k]] = mean_rank / len(values)
        i = j + 1
    return pct


def main():
    matrix, counts_file, genes_tsv, full_out, top_dir = sys.argv[1:6]
    groups = [(g.split("=")[0], g.split("=")[1].split(",")) for g in sys.argv[6:]]

    length, contig = {}, {}
    with open(counts_file) as fh:
        for line in fh:
            if line.startswith(("#", "Geneid")):
                continue
            c = line.rstrip("\n").split("\t")
            length[c[0]] = int(c[5])
            contig[c[0]] = c[1].split(";")[0]

    annot = {r["gene_id"]: r for r in read_tsv(genes_tsv)}
    same_hit = {}
    for r in annot.values():
        if r["ath_agi"]:
            same_hit[r["ath_agi"]] = same_hit.get(r["ath_agi"], 0) + 1

    rows = read_tsv(matrix)
    samples = [c for c in rows[0] if c != "gene"]
    # TPM over every counted gene, then keep protein-coding nuclear genes for ranking
    rpk = {s: {r["gene"]: int(r[s]) / (length[r["gene"]] / 1000) for r in rows} for s in samples}
    tpm = {s: {g: v / sum(rpk[s].values()) * 1e6 for g, v in rpk[s].items()} for s in samples}
    kept = [r["gene"] for r in rows
            if annot.get(r["gene"], {}).get("gene_biotype") == "protein_coding"
            and contig[r["gene"]] not in ORGANELLE]
    pct = {s: dict(zip(kept, percentiles([tpm[s][g] for g in kept]))) for s in samples}

    def score(g, members):
        return min(pct[s][g] for s in members)

    def mean_tpm(g, members):
        return sum(tpm[s][g] for s in members) / len(members)

    out_cols = (["gene_id"] + [f"score_{name}" for name, _ in groups]
                + [f"tpm_{s}" for s in samples] + [f"pct_{s}" for s in samples]
                + ["n_same_ath_hit"] + ANNOT)

    def record(g):
        a = annot.get(g, {})
        rec = {"gene_id": g, "n_same_ath_hit": same_hit.get(a.get("ath_agi", ""), 0)}
        for name, members in groups:
            rec[f"score_{name}"] = f"{score(g, members):.4f}"
        for s in samples:
            rec[f"tpm_{s}"] = f"{tpm[s][g]:.1f}"
            rec[f"pct_{s}"] = f"{pct[s][g]:.4f}"
        for c in ANNOT:
            rec[c] = a.get(c, "")
        return rec

    os.makedirs(os.path.dirname(full_out) or ".", exist_ok=True)
    os.makedirs(top_dir, exist_ok=True)
    with open(full_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, out_cols, delimiter="\t")
        w.writeheader()
        for g in kept:
            w.writerow(record(g))
    print(f"{full_out}: {len(kept)} protein-coding nuclear genes x {len(samples)} galls")

    for name, members in groups:
        ranked = sorted(kept, key=lambda g: (-score(g, members), -mean_tpm(g, members)))[:TOP]
        cols = (["rank", "gene_id", f"score_{name}"] + [f"tpm_{s}" for s in samples]
                + ["n_same_ath_hit"] + ANNOT)
        path = os.path.join(top_dir, f"method1_expression_{name}_top20.tsv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for i, g in enumerate(ranked, 1):
                w.writerow({"rank": i, **record(g)})
        print(f"{path}: galls {','.join(members)}")


if __name__ == "__main__":
    main()
