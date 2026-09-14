# Reference: Haryono et al. 2019 pipeline

Haryono M, Cho S-T, Fang M-J, Chen A-P, Chou S-J, Lai E-M, Kuo C-H (2019).
"Differentiations in Gene Content and Expression Response to Virulence Induction Between
Two Agrobacterium Strains." Front. Microbiol. 10:1554. doi:10.3389/fmicb.2019.01554
Raw reads: NCBI SRA SRP156105. Genomes: C58 (AE007869-AE007872), 1D1609 (CP026924-CP026928).

Closest published analog to our project: comparative genomics + transcriptomics of two
Agrobacterium strains (C58 nopaline-type vs 1D1609 octopine-type), analogous to our
1416 vs 29 comparison. KEY DIFFERENCE: pure cultured bacteria (97.7% mapping rate),
not gall tissue — no host genome, no dual RNA-seq machinery anywhere in their pipeline.

## Design
- 2 strains x 2 conditions (acetosyringone-induced vs DMSO control) x 3 biological replicates = 12 samples
- Induction: 200 uM acetosyringone in AB medium pH 5.5, 16 h — mimics virulence condition in vitro
- Induction verified by qRT-PCR on virB1, virB11, virD2 (dnaE internal control)

## Wet lab
- 5 ug total RNA per sample; Ribo-Zero rRNA Removal Kit (Bacteria)
- TruSeq Stranded mRNA kit, polyA-capture step skipped (bacterial RNA)
- 101 bp single-end, HiSeq 2500, 12 libraries pooled on one flowcell
- Result: ~26M raw reads/sample; rRNA reads only 0.1% of mapped (depletion very effective)

## Bioinformatics

### Read processing
- Trim at first 5' base with Q<20; discard reads <50 bp → ~21M reads/sample retained

### Alignment
- BWA v0.7.12 to each strain's own complete genome (all replicons)
- ~97.7% mapping rate
- SAMtools v1.2 + BEDTools v2.17.0 for read counting
- Read fate: 84.1% sense CDS, 7.0% antisense CDS, 7.3% ncRNA, rest tRNA/rRNA/pseudo/intergenic

### Differential expression
- NOISeqBio (NOISeq R package): RPKM normalization (lc=1), CPM low-count filter (cpm=1)
- DE criteria: >=2-fold change in mean RPKM AND DE probability >=0.99
- Cross-checked TMM and TPKM normalization (R^2 > 0.99 vs RPKM)
- Found: 88/5,355 DE genes in C58, 155/5,630 in 1D1609
- Antisense/ncRNA deliberately not analyzed in depth (noise concerns)

### Comparative genomics
- OrthoMCL for homologous gene clusters (~4,115 shared; >1,000 unique per strain)
- MUSCLE v3.8 alignment of shared single-copy genes → PHYLIP for AA similarity (92.5% overall)
- MAUVE for whole-genome/replicon alignment
- BlastKOALA → KEGG Orthology → COG functional categories
- Per-replicon conservation gradient: circular chr (94.6% AA sim, 76% genes shared)
  > linear chr (92.0%, 59%) > pTi (71.7%, 31%) > pAt plasmids (3% genes shared)

### Cross-strain DE comparison (their Table 2 — template for our step 8)
- Each DE gene classified by homolog's status in the other strain:
  up/up (23), up/not-significant, up/absent, down/down (20), etc.
- <50% of DE genes had a homolog with the same expression pattern
  → regulatory divergence, not just gene content, explains phenotype differences

### Promoter analysis (template for our step 9)
- Input: 600 bp upstream of DE genes, 4 groups (up/down x strain)
- MEME ab initio: "-dna -nmotifs 5", -maxw swept 50 → 15
- Targeted scan: VirG-binding consensus RTTDCAWWTGHAAY, <=3 mismatches allowed
- Operon rule: genes <50 bp apart, same strand = one regulon; upstream gene's motif
  counts for downstream genes
- RESULT: VirG box found upstream of most up-regulated genes (44/52 C58, 61/74 1D1609),
  including genes off the Ti plasmid. NO robust novel motif ab initio — results were
  parameter-sensitive and unstable. Calibrates expectations for our promoter goal.

## Notable biology (relevant to interpreting our data)
- vir regulon on pTi dominates up-regulated genes; most down-regulated genes are chromosomal
- vir gene locations vary between strains (1D1609 has virA/virK/virH2 on pAt, not pTi) —
  check where our strains' vir genes sit before assuming pTi
- tzs (trans-zeatin secretion): nopaline-strain-specific, >200-fold induced in C58, absent
  in 1D1609 — check presence in 1416/29
- Iron/zinc transport genes consistently down-regulated under induction in both strains
