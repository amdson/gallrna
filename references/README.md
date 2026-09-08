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
pipelines — harmonize via orthology (pipeline step 7) before cross-strain gene comparisons.

## agrobacterium/C58/  — published reference (comparison/orthology only, no samples map to it)
Source: `C58_ALIGNED/*.dna` (SnapGene, = AE007869 / AE007870.2 / AE007872.2 / NC_003065).
4 replicons, 5.67 Mb; Atu locus tags. `genes/`: single-gene refs (6b protein, recA).

## plant_host/  — EMPTY, still needed
Host accessions (CC54, CC547, Eu635, Geu182, HC83, M26, P46, T29, T49): species unknown,
ask collaborator. Required before alignment (see confounding notes, pipeline step 3-4).

## Also still missing
- G-19 / G-30 event/construct sequences (not in upload)
- "Strain 29 T-DNA.dna" (was in Windows listing, not uploaded)
