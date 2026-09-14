#!/usr/bin/env python3
"""One row of per-sample mapping statistics for the paper (see PIPELINE.md step 4, `make mapstats`).

Usage: mapstats.py <sample> <reference> <bam> <hisat2.log> <fastp.json> <featureCounts.txt> [threads]

Needs samtools on PATH (the Makefile loads the module). Writes a header line and one data
row to stdout. Host = every contig without the `agro_` prefix; bacterium = `agro_` contigs.

Read-level columns count primary alignments only (no secondary/supplementary), so each read
is counted once, and are expressed per total reads sequenced (2 x input pairs), not per mapped
read: unmapped reads are almost all plant, so a per-mapped-read fraction inflates the bacterial
share wherever the host maps poorly (logs/mapping_fractions.tsv has that problem).
"""
import json
import os
import re
import subprocess
import sys
import tempfile

PREFIX = "agro_"
COLUMNS = [
    "sample", "reference",
    "raw_pairs", "trimmed_pairs", "pct_pairs_kept",
    "pct_overall_alignment", "pct_concordant_unique", "pct_concordant_multi",
    "host_reads", "agro_reads", "pct_reads_host", "pct_reads_agro", "agro_reads_per_million",
    "host_mismatch_rate",
    "pct_pairs_assigned", "pct_pairs_multimapping", "pct_pairs_nofeature", "pct_pairs_ambiguous",
    "host_pairs_assigned", "agro_pairs_assigned", "host_genes_ge10", "agro_genes_ge10",
]


def samtools(*args):
    return subprocess.run(["samtools", *args], check=True, capture_output=True, text=True).stdout


def fastp_pairs(path):
    s = json.load(open(path))["summary"]
    return s["before_filtering"]["total_reads"] // 2, s["after_filtering"]["total_reads"] // 2


def hisat2_log(path):
    text = open(path).read()
    pairs = int(re.search(r"^(\d+) reads; of these:", text, re.M).group(1))
    once = int(re.search(r"(\d+) \([\d.]+%\) aligned concordantly exactly 1 time", text).group(1))
    multi = int(re.search(r"(\d+) \([\d.]+%\) aligned concordantly >1 times", text).group(1))
    overall = float(re.search(r"([\d.]+)% overall alignment rate", text).group(1))
    return pairs, once, multi, overall


def host_stats(bam, host_contigs, threads):
    """samtools stats over host contigs, primary mapped alignments only."""
    with tempfile.NamedTemporaryFile("w", suffix=".regions", delete=False) as fh:
        for name, length in host_contigs:
            fh.write(f"{name}\t1\t{length}\n")
        regions = fh.name
    try:
        out = samtools("stats", "-@", str(threads), "-F", "0x904", "-t", regions, bam)
    finally:
        os.unlink(regions)
    sn = {}
    for line in out.splitlines():
        if line.startswith("SN\t"):
            key, value = line.split("\t")[1:3]
            sn[key.rstrip(":")] = value
    return int(sn["reads mapped"]), float(sn["error rate"])


def featurecounts(path):
    summary = {}
    for line in open(path + ".summary"):
        key, value = line.rstrip("\n").split("\t")[:2]
        if key != "Status":
            summary[key] = int(value)
    host = agro = host_genes = agro_genes = 0
    for line in open(path):
        if line.startswith(("#", "Geneid")):
            continue
        cols = line.rstrip("\n").split("\t")
        count = int(cols[-1])
        if cols[1].split(";")[0].startswith(PREFIX):
            agro += count
            agro_genes += count >= 10
        else:
            host += count
            host_genes += count >= 10
    return summary, host, agro, host_genes, agro_genes


def pct(x, n):
    return f"{100.0 * x / n:.2f}" if n else "NA"


def main():
    sample, reference, bam, hlog, fjson, counts = sys.argv[1:7]
    threads = sys.argv[7] if len(sys.argv) > 7 else "1"

    raw, trimmed = fastp_pairs(fjson)
    pairs, once, multi, overall = hisat2_log(hlog)
    reads = 2 * pairs

    contigs = [l.split("\t")[:2] for l in samtools("idxstats", bam).splitlines()]
    host_contigs = [(n, ln) for n, ln in contigs if n != "*" and not n.startswith(PREFIX)]
    agro_contigs = [n for n, _ in contigs if n.startswith(PREFIX)]
    host_reads, mismatch = host_stats(bam, host_contigs, threads)
    agro_reads = int(samtools("view", "-c", "-@", threads, "-F", "0x904", bam, *agro_contigs))

    fc, host_pairs, agro_pairs, host_genes, agro_genes = featurecounts(counts)
    fragments = sum(fc.values())
    host_annotated = host_pairs > 0

    row = [
        sample, reference,
        raw, trimmed, pct(trimmed, raw),
        f"{overall:.2f}", pct(once, pairs), pct(multi, pairs),
        host_reads, agro_reads, pct(host_reads, reads), pct(agro_reads, reads),
        f"{1e6 * agro_reads / reads:.1f}",
        f"{mismatch:.5f}",
        pct(fc.get("Assigned", 0), fragments), pct(fc.get("Unassigned_MultiMapping", 0), fragments),
        pct(fc.get("Unassigned_NoFeatures", 0), fragments), pct(fc.get("Unassigned_Ambiguity", 0), fragments),
        host_pairs if host_annotated else "NA", agro_pairs,
        host_genes if host_annotated else "NA", agro_genes,
    ]
    print("\t".join(COLUMNS))
    print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main()
