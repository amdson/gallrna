# Reference genomes

Converted 2026-09-08 from uploaded SnapGene/GenBank files with `../scripts/convert_refs.py`
(re-runnable; contig names carry the pipeline's `agro_` prefix for combined-reference mapping).
Each strain = one multi-contig FASTA + GFF3 (gene/CDS/tRNA/rRNA features, Name/product attrs).

## agrobacterium/strain_1416/  — mapping reference for 1416* samples
Source: `Strain 1416/*FINAL.gb` (annotated). 5 replicons, 5.86 Mb total:
circ chr 2,837,379 | lin chr 2,043,296 | pAt1 589,645 | pAt2 188,396 | pTi 196,706.
`fragment_10294bp_LC.*`: 10,294 bp linear-chromosome fragment (8 CDS: adenylate cyclase,
RNase D, aa-tRNA synthetases, glutamate synthase, transporters) — purpose TBD; possibly
related to G-19/G-30 events. Original .gb kept alongside.

## agrobacterium/strain_29/  — mapping reference for 29* samples
Source: RAST annotation (`357.66x *.gbk`) — the "Original Download" .gb files are
sequence-only (0 CDS), same assembly, so RAST versions used. 5 replicons, 5.84 Mb total:
circ chr 2,837,203 | lin chr 2,043,404 | pAt1 589,644 | pAt2 188,391 | pTi 177,708.
NOTE: strains 1416 and 29 are near-identical in replicon sizes — close relatives;
biggest difference is pTi (196.7 vs 177.7 kb).
NOTE: 1416 annotation (SnapGene FINAL) and 29 annotation (RAST) come from different
pipelines — harmonize via orthology (pipeline step 8) before cross-strain gene comparisons.

## agrobacterium/C58/  — published reference (comparison/orthology only, no samples map to it)
Source: `C58_ALIGNED/*.dna` (SnapGene, = AE007869 / AE007870.2 / AE007872.2 / NC_003065).
4 replicons, 5.67 Mb; Atu locus tags. `genes/`: single-gene refs (6b protein, recA).

## plant_host/<host>/  — host genomes (NCBI datasets API, `make host-genomes`)
Each dir holds `<host>.fasta`, plus `<host>.gff3` when NCBI has gene models. On Ceres each dir
is a symlink into `/90daydata/small_grains/andrew.dickson/gallrna_data/references/plant_host/`.
Host per sample: `../samples.tsv`; availability notes: `plant_host/AVAILABILITY.md`.

| dir | assembly | gene models | samples mapped to it |
|---|---|---|---|
| citrus_sinensis | GCF_022201045.2 DVS_A1.0 (RefSeq release 103) | yes | 29wtHC83; the 4 Carrizo (CC) galls for now |
| poncirus_trifoliata | GCA_018350135.1 ZK8 | not on NCBI (at CGD, see CITRUS_HOST_PLAN.md §3.1) | none |
| carica_papaya | GCF_000150535.2 Papaya1.0 | yes | 29wtP46 |
| solanum_lycopersicum | GCF_036512215.1 SLM_r2.1 | yes | 1416wtT49, 29wtT29 |
| brassica_juncea | GCA_018703725.1 | no | 1416wtM26, 29wtM26 |
| euonymus_japonicus_proxy | GCA_963580455.1 (*E. europaeus*) | no | the 4 Euonymus (Eu, Geu) samples |

`make annotate` (Makefile step 7) adds, per annotated host:
`<host>.protein.faa` (RefSeq proteins), `<host>.longest.faa` (longest per gene, named by gene
ID), `<host>.emapper.*` (eggNOG-mapper), `<host>.dmnd`, `<host>.ath.tsv` / `.ath_rev.tsv`
(DIAMOND best hits vs Arabidopsis, each way) and the joined `<host>.genes.tsv`.

## plant_host/arabidopsis_thaliana/  — annotation helper, not a sample host
Ensembl Plants release 63 TAIR10 proteome: longest isoform per gene (`.longest.faa`, `.dmnd`)
and `arabidopsis_thaliana.agi.tsv` (AGI, symbol, description). Built by `make annotate`.

## eggNOG-mapper database  — `EGGNOG_DATA` (Makefile default `plant_host/eggnog_data`)
emapperdb 5.0.2 from eggnog5.embl.de, ~50 GB unpacked (eggnog.db, eggnog_proteins.dmnd,
taxonomy). On Ceres it lives at `/90daydata/small_grains/andrew.dickson/gallrna_data/references/eggnog_data`;
pass `EGGNOG_DATA` or set it in `local.mk`.

## combined/  — `make combined`
`<strain>__<host>.fasta` / `.gff3` (strain + host concatenated) and the HISAT2 index, one per
strain x host pair actually sampled. Symlink into 90daydata on Ceres.

## Also still missing
- G-19 / G-30 event/construct sequences (not in upload)
- "Strain 29 T-DNA.dna" (was in Windows listing, not uploaded)
