# gallrna: dual RNA-seq of *Agrobacterium* galls

> The one document to read. State of the code and data as of **2026-09-17**, checked against the
> Makefile, the scripts, the outputs on disk and the Slurm logs. It replaces `HANDOFF.md`,
> `PIPELINE.md`, `CITRUS_HOST_PLAN.md`, `HARYONO2019_PIPELINE.md` and `SETUP.md` (all still in
> git history up to commit `d5db1de`). Section numbers are stable: the Makefile, `SOURCES.md` and
> the scripts cite them as "README.md §N" or "README.md step N".

Contents: [1 Overview](#1-overview-and-goal) · [2 Status](#2-status-at-a-glance) ·
[3 Samples](#3-samples-and-design) · [4 Where things live](#4-where-things-live) ·
[5 Pipeline](#5-pipeline-steps-0-9) · [6 Citrus: what we know](#6-citrus-hosts-what-we-know) ·
[7 Citrus resources](#7-citrus-resources) · [8 Citrus analysis plan](#8-citrus-analysis-plan) ·
[9 Bacterial side](#9-bacterial-side-the-haryono-2019-template-steps-8-9) ·
[10 Findings](#10-findings-so-far) · [11 Decisions](#11-decisions-and-why) ·
[12 Next steps](#12-next-steps-in-order) · [13 Collaborator asks](#13-asks-for-the-collaborators) ·
[14 Setup and running](#14-setup-and-running) · [15 Gotchas](#15-gotchas-and-housekeeping) ·
[16 Other docs](#16-other-docs)

## 1. Overview and goal

Dual RNA-seq of *Agrobacterium* galls: 14 gall samples from two strains (1416 and 29) on six host
species. Gall tissue holds plant and bacterial RNA together, so each sample is aligned to one
combined reference (its strain + its host) and the reads are split by organism afterwards
(competitive mapping; method sources in `SOURCES.md` §1).

**Current focus: the 5 citrus galls**: one Hamlin sweet orange gall (29wtHC83) and four Carrizo
citrange galls (strain 1416 wild type at 54 and 547 days, plus the engineered strains G-19 and
G-30 at 547 days).

**Goal (set 2026-09-14):** give the collaborators (Shatters and Thomson labs) a **shortlist of
host promoters/genes to target in a follow-up experiment**. No gall has an uninoculated control
and every condition is n = 1, so galls are compared against public healthy-tissue RNA-seq of the
same genotype, and several simple ranking methods are compared by hand (§8.1). Everything else
in the citrus plan is supporting analysis and is optional.

Working rules: keep methods simple and established; run anything heavier than editing as a Slurm
job (§14).

## 2. Status at a glance

| step | target | status |
|---|---|---|
| 0-5 verify, QC, trim, references, align, count | `make all` | **done** for all 14 samples (2026-09-08, job 21995291). Matrices: `04_matrix/agro_strain_{1416,29}.tsv`, `plant_{citrus_sinensis,carica_papaya,solanum_lycopersicum}.tsv`. *Euonymus* and *B. juncea* have no gene models: mapping fractions only |
| 5b mapping statistics | `make mapstats` | **done** 2026-09-14 -> `logs/mapping_stats.tsv`, the per-sample table for the paper. Use it, not `logs/mapping_fractions.tsv` |
| 5c public baselines | `make baselines` | **done** 2026-09-15 (job 22046668, 3 h 25 min): 27 tier-1 runs -> `04_matrix/baseline_citrus_sinensis.tsv` (28,080 genes x 27 runs). Alignment rates came out as predicted (step 5c). Tier 2 and the 6 *Poncirus* runs not run |
| 5d shortlist method 1 | `make shortlist` | **done** 2026-09-15 -> `results/shortlist/method1_expression_{hamlin,carrizo_wt,carrizo_eng,all}_top20.tsv` |
| shortlist methods 3, 4, 6 | notebook 03 | **ready to run, not run.** The cells are written and guarded on the baseline matrix; the notebook was last executed before that matrix existed, so they still print "pending" |
| shortlist methods 2, 5; merged table; promoter sequences | — | not started (method 2 needs R/DESeq2, not installed) |
| 6 differential expression | — | not started |
| 7 host annotation | `make annotate` | sweet orange **done** 2026-09-15 -> `citrus_sinensis.genes.tsv`. Tomato: only the Arabidopsis best-hit half (job 22047373); no eggNOG run, no `genes.tsv`. Papaya not run. *Poncirus* ZK8 not supported yet |
| 8-9 bacterial orthology, promoters | — | not started |
| two-parent Carrizo reference | — | planned (§8.2 A), not started |

Unmerged work outside `main`: an exploratory tomato strain contrast sits untracked in the
`baselines` worktree (§15).

## 3. Samples and design

Raw data: `30-1348328766/00_fastq/`, Azenta/GeneWiz, Illumina NovaSeq X, 2x150 PE, 14 samples
(R1 + R2 + md5 each), 25-32 M read pairs per sample. Sample names are
`<strain><genotype><host code><gall age in days>`. Metadata: `samples.tsv`.

Libraries are **polyA mRNA** (~1% of located reads are intronic, `logs/intron_fraction.tsv`) and
**unstranded** (featureCounts `-s 1` vs `-s 2` on a 2.6 M-read slice of 29wtHC83 assigned
250,875 vs 249,061). Both were inferred from the data on 2026-09-14; the kit is still unknown
(§13). `STRAND=0` is therefore correct everywhere.

**Strain 1416 (8 samples):** wt on CC54, CC547, Eu635, Geu182, M26, T49; engineered events G-19
and G-30 on CC547. **Strain 29 (6 samples, all wt):** Eu635, Geu182, HC83, M26, P46, T29.
Shared hosts for strain comparisons: Eu635, Geu182, M26, and tomato at different ages (T49/T29).

| code | host | reference | gene models |
|---|---|---|---|
| Eu, Geu | *Euonymus japonicus* (green, golden) | *E. europaeus* proxy GCA_963580455.1 | none |
| HC | *Citrus sinensis* 'Hamlin' | DVS_A1.0 GCF_022201045.2 | RefSeq |
| CC | **Carrizo citrange** (*C. sinensis* x *P. trifoliata* F1) | DVS_A1.0 for now (§6.2) | RefSeq (sweet-orange half only) |
| M | *Brassica juncea* | GCA_018703725.1 | none |
| P | *Carica papaya* | Papaya1.0 GCF_000150535.2 | RefSeq |
| T | *Solanum lycopersicum* | SLM_r2.1 GCF_036512215.1 | RefSeq |

The Makefile spells these out by hand rather than deriving them: `SAMPLES_<host>` lists the
samples on each host, `COMBO_<sample>` names the combined reference each sample is aligned to,
`COMBOS` the nine strain x host pairs, and each count matrix lists its columns.

**Strain assemblies** (converted to FASTA/GFF3 on 2026-09-08; details in `references/README.md`):
strain 1416 and strain 29 each have a circular and a linear chromosome, pAt1, pAt2 and pTi. The
two are near-identical in replicon sizes; the biggest difference is pTi (196.7 vs 177.7 kb). 1416
is annotated from the collaborators' SnapGene `FINAL.gb` files, 29 from RAST, so gene calls must
be harmonised by orthology before cross-strain comparisons (step 8). C58 (AE007869 circular,
AE007870 linear, AE007871 pTi, AE007872 pAt) is kept for comparison only; nothing maps to it.

**Limits of the design.** No uninoculated tissue from any host; n = 1 per condition. For
wt vs G-19/G-30 on CC547, treat the result as a fold-change ranking, or treat G-19 + G-30 as two
replicates of "engineered". Still missing from the collaborators: the G-19/G-30 construct
sequences (possibly related to `Strain 1416/1416 10294 bp fragment (LC).gb`), and "Strain 29
T-DNA.dna" (listed, never uploaded).

## 4. Where things live

| what | where |
|---|---|
| code, docs, sample sheets, shortlist results | this repo, `~/gallrna`; GitHub `amdson/gallrna`, branch `main` |
| raw reads | `30-1348328766/` -> `/90daydata/small_grains/30-1348328766` (group dir, 60 GB) |
| pipeline outputs | `01_trim 02_align 03_counts 04_matrix qc references/combined references/plant_host/<host>`: symlinks into `/90daydata/small_grains/andrew.dickson/gallrna_data/` |
| public baselines | `gallrna_data/05_baseline/` on 90daydata. **No symlink in the main checkout** (§15) |
| eggNOG-mapper database | `gallrna_data/references/eggnog_data` (48 GB); set `EGGNOG_DATA` to it |
| strain references | `references/agrobacterium/` (real directory, gitignored, built by `make strain-refs`) |
| small logs and QC tables | `logs/` (gitignored). The baseline summary is in the `baselines` worktree's `logs/` |
| figures | `figures/` (gitignored; regenerated by the notebooks) |
| collaborator files (gitignored, unpublished) | `Strain 1416/`, `Strain 29/`, `C58_ALIGNED/`, `docs/` (pipeline notes, the resource-paper draft, rating tables, PDFs of Alabed 2023/2024 and the 1D159 genome paper) |

Disk use on 90daydata (2026-09-17): trim 52 GB, align 42 GB, references 66 GB, qc 1.8 GB,
**`05_baseline` 252 GB** (fastq 95, trim 96, align 62), counts and matrices under 30 MB. The
baseline fastq and trimmed reads are no longer needed once the counts are trusted; the downloads
are reproducible from `citrus_baselines.tsv`.

**90daydata deletes files not accessed for 90 days.** The first outputs date from 2026-09-08, so
untouched files become eligible around early December 2026. Nothing has been copied to
`/project/small_grains` yet. Copy the small, expensive-to-recreate files there: `04_matrix/*`,
`03_counts/`, `05_baseline/counts/`, `logs/mapping_stats.tsv`, the baseline summary and
`citrus_sinensis.genes.tsv`.

**Repo map**

| file | contents |
|---|---|
| `Makefile` | steps 0-5d and 7; `make help` lists targets. Sample, host, combo and baseline-run lists are written out literally (edit in place). `local.mk` (gitignored) holds per-cluster overrides |
| `samples.tsv` | the 14 samples: strain, genotype, host, gall age, mass, source lab |
| `citrus_baselines.tsv` | 42 public baseline runs (33 tier 1, 9 tier 2), each checked against NCBI/ENA labels on 2026-09-14 |
| `environment.yml` | conda env `rnaseq` (fastp, multiqc, biopython; also lists STAR, Salmon, MEME, DESeq2 for later steps) |
| `scripts/convert_refs.py` | SnapGene/GenBank replicons -> one FASTA + GFF3 with `agro_` contig names and `gene` rows |
| `scripts/merge_counts.py` | per-sample featureCounts files -> one matrix, keeping or dropping `agro_` contigs |
| `scripts/mapstats.py` | one row of mapping statistics per sample (step 5b) |
| `scripts/ena_fetch.sh` | one paired run from ENA with md5 check (step 5c) |
| `scripts/shortlist_expression.py` | shortlist method 1 (step 5d) |
| `scripts/longest_proteins.py`, `scripts/gene_annotation.py` | annotation helpers (step 7) |
| `scripts/intron_fraction.sh`, `scripts/splice_rates.sh` | one-off library checks -> `logs/intron_fraction.tsv`, `logs/splice_rates.tsv` |
| `scripts/pipeline.slurm`, `baselines.slurm`, `annotate.slurm`, `conda-env.slurm` | batch jobs for the heavy targets |
| `scripts/dotplot.slurm` | minimap2 whole-genome alignment of two host assemblies (default: Carrizo's parents) -> a slim PAF next to the query genome (§6.2) |
| `scripts/emapper-env.sh`, `scripts/viz-env.sh` | build the `~/.venvs/emapper` and `~/.venvs/viz` venvs |
| `scripts/claude-*` | Claude Code on compute nodes (`SLURM.md`) |
| `notebooks/01_results_overview` | everything the pipeline produced, QC to matrices (`.html` pre-rendered) |
| `notebooks/02_dataset_structure` | design matrix, reference inventory, library composition |
| `notebooks/03_citrus_shortlist` | the shortlist work off the cluster (step 5d, §8.1) |
| `notebooks/04_parent_dotplot` | dotplot, chromosome pairing and base-level divergence of Carrizo's two parents (§6.2) |
| `results/shortlist/` | the method-1 top-20 tables (in git) |

## 5. Pipeline (steps 0-9)

Steps 0-5d and 7 are automated in the Makefile. Settings: `THREADS` (default 8), `STRAND`
(default 0), `HOSTSEL="citrus_sinensis ..."` or `SAMPLES=...` to restrict to a subset. Run heavy
targets as batch jobs (§14).

### Step 0. Verify transfer (`make verify`)
`md5sum -c` on the delivered files -> `logs/md5.ok`.

### Step 1. Read QC (`make qc`)
FastQC (module) + MultiQC (conda env) -> `qc/multiqc_report.html`.

### Step 2. Trim (`make trim`)
fastp with default settings -> `01_trim/<sample>_R{1,2}.fastq.gz` plus per-sample JSON/HTML.
97-99% of pairs kept.

### Step 3. References (`make refs combined`)
- `make strain-refs`: `scripts/convert_refs.py` turns the collaborators' `.gb`/`.gbk` replicons
  into `references/agrobacterium/strain_{1416,29}/`. Strain contigs carry an `agro_` prefix
  (`agro_1416_pTi`, ...). The source files have spaces in their names, so make can't track them;
  delete the outputs to force a reconversion.
- `make host-genomes`: genome FASTA (+ GFF3 when NCBI has gene models) from the NCBI datasets API.
- `make combined`: one combined FASTA + GFF3 + HISAT2 index per strain x host pair actually
  sampled (`references/combined/<strain>__<host>`). Plant contigs keep their NCBI accessions, so
  names never collide and the organism split keys on `agro_`.
- Carrizo (CC) is aligned to sweet orange only for now (§6.2). Still to add: the G-19/G-30 event
  sequences as extra contigs for those two samples.

### Step 4. Align, competitive mapping (`make align fractions`)
`hisat2 --dta` against the sample's combined index, sorted and indexed BAM -> `02_align/`.
`samtools idxstats` by contig prefix -> `logs/mapping_fractions.tsv` (superseded by step 5b for
anything reported). HISAT2 also attempts spliced alignment on the bacterial contigs; the cost is
small (under 1% of bacterial alignments spliced in most samples, `logs/splice_rates.tsv`).

### Step 5. Count (`make counts matrices`)
`featureCounts -p --countReadPairs -s 0 -t gene -g ID` on the combined GFF3 -> `03_counts/`.
`scripts/merge_counts.py` splits by the `agro_` prefix into `04_matrix/agro_strain_<strain>.tsv`
and `04_matrix/plant_<host>.tsv` (annotated hosts only). T-DNA genes can multimap between pTi
and the transformed plant genome; flag them when interpreting bacterial counts.

### Step 5b. Mapping statistics for the paper (`make mapstats`)
**Why.** Every host is aligned to a single-copy (haploid) reference, but the plants are
heterozygous diploids (sweet orange), an F1 hybrid (Carrizo), an allotetraploid (*B. juncea*) or
a different species from the reference (*Euonymus* on the *E. europaeus* proxy). Reads from a gene
copy that differs from the reference carry extra mismatches; past HISAT2's default limit (about
5 high-quality mismatches per 150 bp read) they fail to align. This is reference bias (Stevenson
et al. 2013; `SOURCES.md` §6). Counts compared within one host genotype carry the same bias in
every sample, so it cancels and is normally not corrected. This target records, per sample, the
numbers that show the bias stayed small, as a supplementary table. Nothing downstream reads it.

**Output.** `logs/mapping_stats.tsv`, one row per sample (rows cached as
`02_align/<sample>.mapstats.tsv`). `samtools stats` reads every host alignment, so run it on a
compute node.

| columns | from | what it shows |
|---|---|---|
| `raw_pairs`, `trimmed_pairs`, `pct_pairs_kept` | fastp JSON | read retention after trimming |
| `pct_overall_alignment`, `pct_concordant_unique`, `pct_concordant_multi` | HISAT2 log | divergence from the reference lowers overall alignment; repeats or homeologs raise multi-mapping |
| `host_reads`, `agro_reads`, `pct_reads_host`, `pct_reads_agro`, `agro_reads_per_million` | samtools, primary alignments | reads per organism as a share of all reads sequenced, each read counted once |
| `host_mismatch_rate` | `samtools stats` error rate on host contigs | mismatches per aligned base: the direct measure of divergence (includes sequencing error; inbred tomato is the floor) |
| `pct_pairs_assigned`, `_multimapping`, `_nofeature`, `_ambiguous` | featureCounts summary | fate of read pairs at counting, both organisms together |
| `host_pairs_assigned`, `agro_pairs_assigned`, `host_genes_ge10`, `agro_genes_ge10` | featureCounts table | counted pairs and genes detected per organism; host columns are NA without gene models |

**How to read it** (run of 2026-09-14):

| host | overall alignment | host mismatch rate | bacterial share of reads |
|---|---|---|---|
| tomato | 96% | 0.29-0.31% | 0.14-0.17% |
| Hamlin sweet orange | 92% | 0.56% | 0.06% |
| *B. juncea* | 88-89% | 0.51-0.52% | 0.29-0.36% |
| papaya | 88% | 0.33% | 0.27% |
| Carrizo | 82-85% | 0.88-1.02% | wt 0.08-0.14%; G-19/G-30 ~0 |
| *Euonymus* (proxy) | 31-36% | 1.87-2.06% | 0.04-0.15% |

- Tomato, Hamlin, *B. juncea* and papaya align well with few multi-mapped pairs: no correction
  needed. *B. juncea*'s assembly keeps its A and B subgenomes as separate chromosomes, and reads
  still place uniquely.
- Carrizo loses reads from its *Poncirus* gene copies. Fine for CC-vs-CC contrasts; state it for
  comparisons against other hosts, or use the two-parent reference (§8.2 A).
- *Euonymus* is limited by the species difference, not ploidy: most aligned reads sit near the
  mismatch limit. Report it as a limitation of the proxy, and compare bacterial load with other
  hosts only per total reads (`pct_reads_agro`).
- **Cross-mapping bound.** G-19/G-30 galls barely express the wt bacterial T-DNA genes, so their
  `agro_reads_per_million` (2.7 and 9.8, against 638-1,426 in the wt citrus galls) is an upper
  bound on plant reads wrongly assigned to the bacterium, for citrus. Worth a sentence in the
  methods.
- 1416G-30CC547 has 16.4% concordant multi-mapping pairs against 4.7-8.1% in the other citrus
  galls (22.0 M multimapping pairs at counting vs 4.6-9.2 M). Find out what (rRNA, a repeat, the
  construct?) before trusting it in contrasts.
- `logs/mapping_fractions.tsv` counts secondary alignments too and divides by mapped reads only,
  so unmapped plant reads inflate its bacterial share (about 2x for *Euonymus*).

### Step 5c. Public healthy-tissue baselines (`make baselines`)
- Input: the runs in `citrus_baselines.tsv` (chosen in §7.3), fetched from ENA with
  `scripts/ena_fetch.sh` (portal filereport API, md5-checked; no sra-tools).
- Same settings as steps 2, 4 and 5, but against a **plant-only** HISAT2 index built next to the
  genome (`references/plant_host/<host>/<host>.*.ht2`): baselines have no bacterial reads to
  compete with, and the same NCBI GFF gives the same gene IDs as the gall matrices.
- Reference choice: Carrizo runs go to sweet orange, like the CC galls (`BASE_REF_carrizo`), so
  gall and baseline share the hybrid-mapping bias. The 6 *Poncirus* runs (PRJDB41296) are listed
  but skipped until an annotated ZK8 reference exists (`BASE_HOSTS` controls this).
- Output: `05_baseline/{fastq,trim,align,counts}/` and `04_matrix/baseline_<host>.tsv`, one
  column per run accession; join tissue and replicate labels from `citrus_baselines.tsv`.
  `scripts/baselines.slurm` also writes `logs/baseline_summary.tsv`.
- Run: `sbatch scripts/baselines.slurm` (48 cores, 4 runs at a time). `BASE_TIER=2` adds the
  tier-2 Carrizo leaf controls, `BASE_RUNS="..."` selects a subset. Rerunning resumes.

**Result (2026-09-15, 27 tier-1 runs).** The rates match the galls of the same genotype, which is
the check that gall and baseline share their mapping bias:

| set (BioProject) | runs | read pairs | alignment rate | pairs assigned to genes |
|---|---|---|---|---|
| sweet orange root (PRJNA599503) | 3 | 26-29 M | 89.0-89.7% | 80% |
| sweet orange bark (PRJNA599503) | 3 | 25-30 M | 89.7-90.2% | 78-80% |
| sweet orange young leaf (PRJNA599503) | 3 | 31-32 M | 91.8-92.4% | 80% |
| Valencia callus EV-L44, 0/15/30 d (PRJNA778304) | 9 | 25-32 M | 93.8-94.3% | 87-88% |
| Carrizo stem (PRJNA1216034) | 3 | 19-27 M | 86.0-86.2% | 81-82% |
| Carrizo young leaf (PRJNA1216034) | 3 | 21-23 M | 86.8-87.4% | 80-81% |
| Carrizo mature leaf (PRJNA1216034) | 3 | 21-25 M | 85.0-85.2% | 79% |

For comparison: Hamlin gall 92.1%, Carrizo galls 81.7-85.2%.

### Step 5d. Candidate shortlist, method 1 (`make shortlist`)
- §8.1 method 1: host genes ranked by expression in the galls themselves, no baseline.
  `scripts/shortlist_expression.py` turns the citrus gall counts into TPM (gene lengths from
  featureCounts), takes each gene's within-gall percentile among protein-coding nuclear genes
  (organellar contigs NC_008334.1 and NC_037463.1 dropped), and scores a group by the gene's
  **lowest** percentile across that group's galls, so only genes high in every gall rank. Ties
  break on mean TPM.
- Groups: `hamlin` (29wtHC83), `carrizo_wt` (1416wtCC547, 1416wtCC54), `carrizo_eng` (G-19, G-30)
  and `all`. Outputs: `results/shortlist/method1_expression_<group>_top20.tsv` (in git, with TPM
  per gall and the annotation columns) and the full 23,436-gene table
  `04_matrix/citrus_gall_expression.tsv`, for merging with methods 2-6. `n_same_ath_hit` = citrus
  genes sharing the same best Arabidopsis hit, a rough multi-copy warning.
- **Result.** The `all` top 20 is the constitutive-promoter set (ribosomal proteins,
  polyubiquitin, cyclophilin ROC1, metallothionein MT2A, NDPK1, TCTP) plus stress and wound genes
  (an MLP423-like major allergen at rank 1, LEA5/SAG21, dehydrin ERD14, GRP7, extensin). Hamlin
  adds aquaporins (PIP1, PIP2, TIP1) and extensins; G-19/G-30 add protease inhibitors, an
  endochitinase and dormancy-associated DYL1. 42 genes appear in at least one list, 7 in all
  four; 27 of the 42 carry a multi-copy or cloning warning. Which of them are gall-specific needs
  the baselines.
- The script runs in seconds. From a fresh worktree run it directly: the checkout's new file
  timestamps make `make` want to rebuild the whole gall chain first.
- **Notebook.** `notebooks/03_citrus_shortlist.ipynb` is the interactive version: filters,
  TPM/percentiles, gall concordance (G-19 and G-30 are near-replicates, Spearman 0.96), method 1
  with interpretation and figures (`figures/shortlist_*.png`), the warnings, and guarded cells
  for methods 3 (percentile shift), 4 (rank products) and 6 (tau) plus a merged table per group
  (`results/shortlist/merged_<group>.tsv`). It checks that its method-1 lists match the script's.

### Step 6. Differential expression (not started)
- R in conda env `rnaseq`: DESeq2, tximport, edgeR, pheatmap, ggplot2. Listed in
  `environment.yml`, not installed.
- Normalise each organism separately: the plant/bacterium ratio varies per sample.
- Contrasts: wt vs G-19/G-30 on CC547; strain 1416 vs 29 on shared hosts (Eu635, Geu182, M26).
  Citrus host genes: §8.2 D.

### Step 7. Host functional annotation (`make annotate`)
- Inputs: host RefSeq proteins (NCBI), the Arabidopsis proteome (Ensembl Plants release 63),
  eggNOG-mapper database emapperdb 5.0.2.
- Longest protein per gene, named by the gene ID featureCounts writes -> eggNOG-mapper 2.1.13
  (Viridiplantae scope: GO, KEGG, EC, Pfam, description) + DIAMOND best hit each way against
  Arabidopsis (AGI, symbol, reciprocal-best-hit flag `ath_rbh`).
- Output: `references/plant_host/<host>/<host>.genes.tsv`, one row per counted gene; joins onto
  `04_matrix/plant_<host>.tsv` by gene ID. The Arabidopsis columns connect to the AGI-based
  crown-gall marker lists (§7.4).
- Run: `sbatch scripts/annotate.slurm [hosts]` with `EGGNOG_DATA` pointing at 90daydata
  (environment or `local.mk`); the database would overflow the 30 GB home quota. Needs
  `scripts/emapper-env.sh` once.
- **Sweet orange (done):** 40,427 RefSeq proteins -> 23,556 protein-coding genes in a 28,080-row
  table; 97% eggNOG-annotated, 49% with GO, 94% with an Arabidopsis hit, 56% reciprocal best hits.
  The eggNOG annotation phase took 5.3 h at 32 GB; `--dbmem` with >= 64 GB would speed it up.
- **Tomato (partial):** 44,391 proteins -> 28,611 genes; `ath.tsv` and `ath_rev.tsv` exist, the
  eggNOG step and `genes.tsv` don't. Papaya not started. Not covered: *Poncirus* ZK8 (needs
  proteins from its CGD GFF) and the unannotated hosts.

### Steps 8-9. Bacterial orthology and promoters (not started)
See §9.

## 6. Citrus hosts: what we know

### 6.1 The five citrus galls

| sample | host | strain | gall age (d) | assigned pairs | genes >= 10 counts |
|---|---|---|---|---|---|
| 29wtHC83 | *C. sinensis* 'Hamlin' | 29 wt | 83 | 26.0 M | 18,239 |
| 1416wtCC547 | Carrizo citrange | 1416 wt | 547 | 21.3 M | 17,023 |
| 1416wtCC54 | Carrizo citrange | 1416 wt | 54 | 20.8 M | 17,132 |
| 1416G-19CC547 | Carrizo citrange | 1416 G-19 | 547 | 22.0 M | 17,239 |
| 1416G-30CC547 | Carrizo citrange | 1416 G-30 | 547 | 20.7 M | 16,753 |

Counts: `04_matrix/plant_citrus_sinensis.tsv`, all five aligned to and counted on sweet orange.
There is no uninoculated citrus tissue, so "what did the bacterium change" is triangulated
(§6.4) rather than read off a gall-vs-mock contrast. No citrus crown-gall transcriptome has been
published.

### 6.2 Both citrus hosts are hybrids

**CC = Carrizo citrange** (confirmed 2026-09-14): *C. sinensis* 'Washington' navel x
*P. trifoliata*, a single F1 clone (Swingle's USDA program, 1909; Carrizo and Troyer are the same
clone). It is normally raised from nucellar (clonal) seedlings, so all four CC plants should
share one genotype. There is no Carrizo assembly. Every Carrizo gene has one sweet-orange copy
and one *Poncirus* copy.

**Hamlin is a hybrid too**: sweet orange carries mandarin and pummelo ancestry, and all sweet
orange cultivars (Hamlin, Valencia, Washington navel) are somatic mutants of one original hybrid.
DVS_A1.0 (Valencia haplotype A) plus DVS_B1.0 (GCA_022201065.1, haplotype B) together describe
Hamlin almost exactly, and also the sweet-orange half of Carrizo.

How it shows in the current alignments (one sweet-orange chromosome, NC_068563.1; first 400k
primary alignments with MAPQ >= 10; `NM` = mismatches to the reference):

| sample | overall alignment | NM = 0 | NM 1-2 | NM >= 3 | MAPQ < 10 |
|---|---|---|---|---|---|
| 29wtHC83 (sweet orange) | 92.1% | 66.2% | 27.1% | 6.7% | 2.4% |
| 1416wtCC547 | 84.5% | 40.5% | 42.5% | 17.0% | 5.8% |
| 1416wtCC54 | 85.2% | 39.7% | 44.1% | 16.1% | — |
| 1416G-19CC547 | 81.7% | 45.1% | 42.3% | 12.7% | — |
| 1416G-30CC547 | 85.0% | 54.7% | 35.8% | 9.5% | — |

Consequences:
- Aligning Carrizo to sweet orange alone loses or mis-scores the *Poncirus* copies. Published
  Carrizo RNA-seq does the same: Afzal Naveed, Huguet-Tapia & Ali 2019 (*J Plant Interact*
  14:187-204, doi:10.1080/17429145.2019.1609106) mapped Carrizo roots to two *C. sinensis*
  genomes, report 55-73% transcriptome coverage, and fell back to de novo assembly. A
  cross-species test here (2026-09-08) assigned 92% as many pairs to genes as the native
  sweet-orange sample.
- All four CC galls share the hybrid genotype, so the bias is common to every CC sample and
  **cancels in CC-vs-CC contrasts**. It does **not** cancel against galls on other hosts, against
  pure sweet-orange or pure *Poncirus* baselines, or in any per-parent question. The Carrizo
  baselines carry the same bias (step 5c), which is why they are mapped the same way.
- **QC flag.** G-30 (9.5% NM >= 3) and G-19 (12.7%) look less hybrid-like than the two wt galls
  (16-17%). A different mix of expressed genes can shift this, so it isn't proof. Before relying
  on wt vs G-19/G-30, check that all four plants are Carrizo: at genes where the parents differ,
  an F1 should split about 50:50 between the parental copies (falls out of §8.2 A).

**The two parental assemblies compared** (2026-09-17; `sbatch scripts/dotplot.slurm`, 4 min, then
`notebooks/04_parent_dotplot.ipynb`; `figures/30_parent_dotplot.png`). Sweet orange DVS_A1.0 against
the NCBI copy of *Poncirus* ZK8, minimap2 `-x asm20`:
- The nine chromosomes pair one to one and are largely collinear. 58% of sweet orange's 299 Mb of
  chromosomes and 74% of *Poncirus*'s 231 Mb sit in an alignment; the gaps are mid-chromosome.
  *Poncirus* chr4, 5, 6, 7 and 9 are assembled in the opposite orientation (a convention, not
  biology).
- Off the diagonal: a reciprocal exchange between chr1 and chr9 (8.5 and 5.1 Mb of alignment) and a
  2.9 Mb block of sweet orange chr7 on *Poncirus* chr5. Real translocations or scaffolding
  differences; two assemblies can't tell.
- *Poncirus* chr8 is only 10 Mb against 31 Mb. The middle of sweet orange chr8 (3.8-24.3 Mb)
  aligns to ZK8's unplaced scaffolds instead, which carry 54 Mb of alignment in total. **The
  two-parent reference (§8.2 A) must keep the unplaced scaffolds.**
- Genome-wide the parents differ at a median 3.8% of aligned bases, but that is mostly non-coding.
  Genes are far better conserved: 85% of sweet-orange exon sequence and 87% of coding sequence
  has a *Poncirus* counterpart (scaffolds counted), with 1.5 substitutions per 100 bp in coding
  sequence and 1.8 in exons, against 2.5 in introns and 3.7 between genes (`paftools.js call`).

**What sweet-orange-only alignment costs Carrizo, measured.** Simulated 150 bp reads from the
*Poncirus* copy of each gene, scored with HISAT2's penalties: 78.5% align to sweet orange, **21.5%
are rejected**. With half the reads from each copy that is ~11% of all reads, predicting ~82%
alignment for a library that would otherwise reach 92%; the Carrizo galls align at 82-85%. So the
observed gap is explained, and it is an upper bound (soft-clipping, and pairs counted on one mate).
- Per gene the loss is uneven: the median gene loses 19% of its *Poncirus*-copy reads (-0.15 log2
  on its count); 12% of genes lose over half (count down by a quarter or more); 2% lose nearly all
  (-1 log2: the count is the sweet-orange copy alone). 18% of genes have no uniquely aligned exon
  and can't be scored. Per-gene values: `04_matrix/carrizo_poncirus_read_loss.tsv`.
- In the method-1 shortlists, 4-5 of each Carrizo top 20 lose over half of their *Poncirus*-copy
  reads (metallothionein type 2 75%, extensin-3 69%, a glycine-rich RNA-binding protein 68%, an
  extensin-like 58%, LEA5 54%). They are undercounted and still rank in the top 20, so the lists
  hold; but for these genes the two parental promoters probably differ too, which matters when
  choosing which copy's promoter to clone.
- Unaffected: every Carrizo-vs-Carrizo comparison, including gall vs Carrizo baselines (same bias
  in both). Affected: Carrizo gall vs *pure sweet-orange* tissue (the bark and callus baselines),
  biased downwards for diverged genes; *Poncirus*-only genes; anything per parental copy.

### 6.3 The bacterium's direct footprint: T-DNA expression in the galls

T-DNA transcripts are made by transformed *plant* cells from plant-active bacterial promoters, so
they are the one part of "which promoters did the bacterium influence" visible directly. Counts
from `04_matrix/agro_strain_*.tsv`:

| strain 1416 T-DNA gene | CC547 wt | CC54 wt | G-19 | G-30 |
|---|---|---|---|---|
| "D-octopine dehydrogenase" (opine synthase-like) | 11,761 | 16,287 | 12 | 6 |
| agrocinopine synthase (*acs*) | 1,184 | 2,793 | 2 | 0 |
| indoleacetamide hydrolase (*iaaH*) | 507 | 711 | 4 | 0 |
| tryptophan monooxygenase (*iaaM*) | 8 | 21 | 0 | 0 |
| E protein / gene 5 / C protein / Atu6012 | 38-953 | 588-953 | 0 | 0-1 |

| strain 29 T-DNA gene (pTi ~78.8-94.2 kb) | HC83 |
|---|---|
| hypothetical peg.5741 (likely opine synthase; confirm by BLAST) | 9,656 |
| *iaaH* | 515 |
| *ipt* (cytokinin) | 469 |
| *iaaM* | 11 |

- Both wt strains express their T-DNA in citrus. G-19/G-30 galls essentially do not, and they
  cluster apart from the wt galls (PC1, 55% of variance; `figures/10_citrus_structure.png`). The
  resource-paper draft (`docs/050136 ...docx`) describes disarming 1416 and installing GAANTRY,
  which would explain it. The G-19/G-30 construct maps are needed to be sure.
- No *ipt* is annotated inside the 1416 T-DNA window. The only 1416 *ipt* is at pTi 164.8 kb among
  the vir genes, which is where *tzs* sits. Check whether 1416's T-DNA *ipt* is missing or just
  unannotated (one of peg 5568-5570?).
- *iaaM* is barely detected in both strains while *iaaH* is well expressed. Biology, or a
  mapping/annotation artefact: worth a look in IGV.

### 6.4 Framing: finding bacterium-influenced genes without a control

Each line of evidence carries a known bias. A host gene is a good candidate when they agree:

1. **Gall vs public healthy tissue of the same genotype** (§7.3). Captures the gall program but is
   confounded with lab, tissue and age; mitigated by using several baselines and tissues.
2. **Internal contrasts** that hold host genotype constant. The best is wt CC547 vs G-19/G-30
   CC547: same host and age, with vs without wt T-DNA expression. These work on any reference.
3. **Prior crown-gall signatures** (Arabidopsis, via orthologs; §7.4), which say which way known
   gall genes should move.
4. **Parental-copy (allele-specific) expression in Carrizo.** Both copies of a gene share one
   nucleus, so a gall-induced shift in the sweet-orange : *Poncirus* ratio has to come from the
   gene's own DNA (cis), most likely its promoter. The most direct route to promoters the design
   allows (§8.2 D5, E1), and it needs the two-parent reference.

## 7. Citrus resources

### 7.1 Host references

**Sweet orange (Hamlin): keep the current reference.** NCBI RefSeq GCF_022201045.2 (DVS_A1.0,
Valencia, haplotype A), Annotation Release 103: 23,556 protein-coding genes, about 7,800 of them
"uncharacterized" and no GO/KEGG in the GFF, hence step 7. Optional later: DVS_B1.0
(GCA_022201065.1, Clemson annotation, 27,940 genes) as the second haplotype; low priority for
gene-level counts. Newer assemblies (T2T SWO GCA_046470765.1/775.1, Pera Rio, Washington Navel,
Newhall) have no NCBI annotation, so there is no reason to switch.

**Carrizo: planned two-parent reference, DVS_A1.0 + *Poncirus* ZK8.**
- *P. trifoliata* ZK8 v1.0 (HZAU; Huang et al. 2021, *Hortic Res*,
  doi:10.1038/s41438-021-00505-2, BioProject PRJNA554539; Peng et al. 2020 is a different,
  UF/JGI assembly). Same assembly as the unannotated NCBI copy already downloaded
  (GCA_018350135.1), but the Citrus Genome Database (CGD) serves gene models too:
  - `https://www.citrusgenomedb.org/jb2/data/Ptri_ZK8_v1.sorted.gff.gz` (7.5 MB; 25,680 genes,
    39,675 mRNAs, with UTRs; IDs like `Pt1g002240`)
  - `https://www.citrusgenomedb.org/jb2/data/Ptri_ZK8_v1.fasta.gz` (86 MB)
  - CGD lists annotation BUSCO 93.4%, assembly BUSCO 98.6%.
  - Contig names differ (CGD `chr1_ZK8`..., NCBI `CM031384.1`...), so **use the CGD FASTA with
    the CGD GFF** rather than renaming. These are the genome-browser data files, publicly
    readable but not a formal download page. Cite Huang et al. 2021.
- Contig names (`NC_0685xx` vs `chrN_ZK8`) and gene IDs (`gene-LOC...` vs `Pt...`) don't collide,
  so the two genomes concatenate as they are, alongside the `agro_` contigs.
- ZK8 is a *Poncirus* accession, not Carrizo's actual pollen parent. *Poncirus* is low-diversity,
  so expect a modest mismatch rate on that side.
- JGI *P. trifoliata* v1.3.1: the page returned 403 and the access terms are **(unverified)**.
  Not needed given ZK8.

**Sweet orange <-> *Poncirus* gene map** (pairs the two parental copies):
- CGD MCscan anchors `jb2/data/mcscan/Csin_DVS_A_v1.Ptri_ZK8_v1.anchors.gz`: 19,501 syntenic
  pairs. The DVS_A side uses Clemson GenBank locus tags (`KPL70_000001`), not RefSeq `LOC...`
  IDs. Same assembly, so map KPL70 -> LOC by coordinate overlap (GCA_022201045.1 GFF vs our
  RefSeq GFF).
- Cross-check: DIAMOND reciprocal best hits, RefSeq proteins vs ZK8 proteins.
- Genes without a partner are counted on their own and flagged.

### 7.2 Functional and regulatory annotation

| resource | what | IDs | use |
|---|---|---|---|
| eggNOG-mapper (step 7) | GO, KEGG, Pfam, COG | ours | primary functional annotation. **In use** |
| DIAMOND vs Arabidopsis (step 7) | best hits, RBH flag | ours | crown-gall marker lists (§7.4). **In use** |
| PlantRegMap `08-download/Citrus_sinensis/` | GO, promoter sequences, TFBS predictions, predicted TF -> target regulation | `orange1.1g...` (JGI v1) | ready-made promoter/TFBS layer; needs an orange1.1 -> LOC map (protein RBH) |
| PlantTFDB `Csi_TF_list.txt.gz` | 2,255 TFs / 1,338 loci in 58 families | `orange1.1g...` | TF families |
| JASPAR 2024 CORE plants | curated TF motifs | — | motif scanning (§8.2 E) |
| PLAZA dicots 5.0 | orthologous gene families, GO | — | blocked automated access; citrus coverage **(unverified)** |

No citrus ATAC-seq/DAP-seq promoter map was found. PRJDB41296 includes one *C. sinensis* DAP-seq
for a single TF, which isn't useful genome-wide.

### 7.3 Public healthy-tissue RNA-seq (baselines)

Run accessions for tiers 1 and 2 are in `citrus_baselines.tsv`. Tier 1 = 33 runs (9 Carrizo,
18 sweet orange, 6 *Poncirus*), tier 2 = 9 Carrizo leaf controls. The 27 Carrizo and sweet-orange
tier-1 runs are counted (step 5c).

**Carrizo citrange (genotype-matched to the CC galls)**

| tier | BioProject | tissue (n) | library | notes |
|---|---|---|---|---|
| 1 | PRJNA1216034 | wt **stems** (3, 18 d old), young leaves (3), untreated mature leaves (3) | cDNA PE HiSeq 2500, ~18-26 M | CsCOI1 study (*Hortic Res* 2025, uhaf174); paper confirms Carrizo as WT with 3 reps; they used SWO v3 + Salmon. SRA organism is just "Citrus". Ignore the cscoi1 mutant, H2O-inoculated and *Xcc*-infected sets |
| 2 | PRJNA839431 | leaf controls "CT-Car" (3) | cDNA PE NovaSeq, ~26 M | stress-combination study; SRA mislabels organism as *C. sinensis*; Cleopatra mandarin in the same project |
| 2 | PRJNA668159 | leaf "Carrizo Control" (3) | PCR PE NovaSeq, ~22 M | FT-overexpression study |
| 2 | PRJNA1053671 | leaf WT "ZC-1..3" (3) | RANDOM PCR PE NovaSeq, ~21 M | SDE transgenic study |
| 3 | PRJNA1092419 | apical leaves (6) | cDNA PE, ~25 M | volatile-exposure study; which runs are controls **(unverified)** |
| 3 | PRJNA532644 | roots, uninfected (2: 24 h, 48 h) | cDNA PE, ~12 M | *Phytophthora* study |
| 3 | PRJNA1126089 / PRJEB20758 | leaf control (1) / tissue unstated (CC T0, 2) | — | too small to use alone |

Gap: no Carrizo bark or older stem. The only stem set is 18-day-old shoots, against 54- and
547-day galls. The pure-parent bark/stem and callus sets fill the tissue gap.

***Poncirus* (one parent; not counted yet)**

| tier | BioProject | tissue (n) | library | notes |
|---|---|---|---|---|
| 1 | PRJDB41296 | ZK stem (3), ZK thorn (3) | polyA PE NovaSeq, ~23 M | same accession family as the ZK8 reference. Thorn-identity study, released 2026-06 |
| 2 | PRJNA482734 | trifoliate stem (6: "N" x3, "FN" x3), leaf | "RANDOM", PE, ~41 M | N/FN meaning unexplained; also has clementine |
| 2 | PRJNA1222399 | 'Flying Dragon' (4), 'Rich 16-6' (6) rootstock | polyA PE HiSeq 4000, ~28 M | tissue **(unverified; likely rootstock stem/bark)** |
| 3 | PRJNA487128 | ZK root, leaf, seed, fruit, ovules (3 each) | cDNA PE, ~35 M | tissue atlas |
| 3 | PRJNA816480 | trifoliate root, non-infected vs CDVd viroid | RANDOM, PE | controls only; root |

**Sweet orange (genotype-matched to Hamlin; also the other parent of Carrizo)**

| tier | BioProject | tissue (n) | library | notes |
|---|---|---|---|---|
| 1 | PRJNA599503 | **healthy bark (3)**, root (3), young leaf (3) | cDNA PE HiSeq 2500, ~28 M | cultivar and polyA vs total **(unverified)**; used in HLB meta-analyses (Peng et al. 2021). Ignore the HLB-infected counterparts |
| 1 | PRJNA778304 | Valencia embryogenic **callus**, empty-vector line EV-L44 at 0/15/30 d (3 each) | cDNA PE NovaSeq, ~27 M | proliferating-tissue comparator; ignore the MIM transgenic lines |
| 2 | PRJNA816480 | sweet orange **stem** on trifoliate rootstock, non-infected vs CDVd | RANDOM, PE | control count per tissue **(unverified)**; PMC9228058 |
| 2 | PRJNA1222399 | sweet orange scion on 'Flying Dragon' (4) / 'Rich 16-6' (6) | polyA PE HiSeq 4000, ~28 M | tissue **(unverified)** |
| 3 | PRJNA386941 | sweet orange stem on *P. trifoliata* (2) or *C. junos* (2) | SE BGISEQ-500, ~24 M x 50 bp | low replication, short SE |
| — | PRJNA703546 | Washington navel leaf | — | leaf only, skip |

Why several: any single study confounds "gall" with its own lab, cultivar and growth conditions.
A gene up in gall vs **every** healthy tissue from **two or more studies** is robust to that. A
gene up vs stem/bark but at callus level is proliferation, not a bacterium-specific effect.

### 7.4 Crown-gall expectations (Arabidopsis)

Deeken et al. 2006 *Plant Cell* 18:3617 (PMC1785400): C58 tumors vs uninfected inflorescence
stalk, 35 dpi, 4 reps, ATH1 array. Also Lee et al. 2009 *Plant Cell* (early infection) and the
Gohlke & Deeken 2014 review (PMC4006022). ArrayExpress E-GEOD-13927 holds Arabidopsis crown-gall
arrays; check whether it is the Deeken set before using it.

**Up in tumors:** GH3 auxin-responsive At2g23170 (56x), SuSy3 At3g43190, SuSy5 At5g20830, STP4
At3g19930, cell-wall invertase At3g13790, PDC1 At4g33070, ADH At1g77120, PIP2;5 At3g54820 (12x),
PIP1;3 At1g01620, SAD At1g43800, FAD3 At2g29980, At2g31360, LTP2 At2g38530 (10x), nitrate
transporters At3g45060 / At5g60780, amino-acid transporters At1g47670 / At1g25530, oligopeptide
transporters At4g21680 (17x) / At1g59740, ASA1 At5g05730, GAD1 At5g17330, PGP1 At2g36910, PGP4
At2g47000, expansins At1g69530 / At1g26770 / At2g28950, PUMP1 At3g54110.

**Down:** BCAT At3g19710 (49x), vacuolar invertases At1g12240 / At1g62660, NRT1 At2g26690 /
At3g21670, AMT2 At4g13510, AMT1;1 At2g38290, TIP2;2 At4g17340, CUT1 At1g68530, WAX2 At5g57800,
and photosynthesis broadly (light reactions, 24/27 Calvin-cycle genes).

**Processes (review):** ABA up with drought/ABA-responsive genes, suberization up and cutin down,
SA/ET markers up, hypoxia/fermentation, sink metabolism, MYB/bHLH/bZIP/AP2 TF changes, genome
hypermethylation (T-DNA oncogene promoters stay unmethylated).

When scoring these in citrus, match gene families, not just best hits: citrus GH3 genes, for
example, hit other Arabidopsis GH3 paralogs.

## 8. Citrus analysis plan

### 8.1 The candidate shortlist (the goal; decided 2026-09-14)

Don't trust one statistic: run several simple ranking methods, take each method's **top 20**, and
compare the lists **by hand** alongside `citrus_sinensis.genes.tsv`. Genes that come up under
several methods are the strongest candidates. Separate lists for Hamlin, the wt Carrizo galls,
and the G-19/G-30 Carrizo galls; the last is the closest to an engineered-strain gall.

| # | method | needs | status |
|---|---|---|---|
| 1 | **Gall expression**: within-sample percentile, consistent across the galls. Measures promoter strength | gall counts | done (step 5d) |
| 2 | **Fold change vs baseline**: DESeq2 with shrunken log2FC (`lfcShrink`) and a minimum-expression cut-off, so near-zero genes can't top the list | baselines, R | not started |
| 3 | **Percentile rank shift**: gall percentile minus the median healthy-tissue percentile. Robust to cross-study depth and library differences | baselines | notebook cell written, not run |
| 4 | **Rank products** across the galls (Breitling et al. 2004, *FEBS Lett* 573:83): rewards consistency | baselines | notebook cell written, not run |
| 5 | **RankComp / relative expression orderings** (Wang et al. 2015, *Bioinformatics* 31:62): gene pairs whose order is stable in healthy tissue and flips in the gall. Designed for one sample against a reference set from other studies | baselines | not started |
| 6 | **Tissue specificity, tau** (Yanai et al. 2005, *Bioinformatics* 21:650), gall as one tissue alongside stem, bark, leaf, root and callus | baselines | notebook cell written, not run |
| 7 | The T-DNA promoters (§6.3), listed separately as bacterial benchmarks | — | — |

Before ranking, drop genes not clearly expressed in the galls, organellar and rRNA genes, and
multi-copy families whose promoter can't be cloned unambiguously. Output: one merged table (gene,
annotation, Arabidopsis ortholog, rank under each method, expression per gall) plus the ~1.5 kb
upstream sequence of each pick. The picks then go to a reporter assay (GUS/luciferase in galls),
which is the actual measurement of promoter activity.

### 8.2 Supporting analyses (optional; explain why a candidate responds)

**A. Two-parent Carrizo reference.** Concatenate DVS_A1.0 (RefSeq FASTA + GFF) and CGD ZK8
(FASTA + GFF) into a `carrizo` host, index it with each strain as now, and realign the 4 CC galls.
Quantify per parental copy, either with Salmon on the combined transcriptome (its EM step shares
reads that fit both copies; the Carrizo CsCOI1 study used Salmon, against sweet orange only) or
with HISAT2 + featureCounts counting only confidently placed reads (MAPQ filter) and reporting
the ambiguous share. Gene-level count = sweet-orange copy + *Poncirus* copy via the §7.1 pairs.
Report each sample's genome-wide *Poncirus*-copy share (an F1 should be about 50%): that is the
§6.2 genotype check. Keep ZK8's unplaced scaffolds in the reference: much of its chr8 is only
there (§6.2). Makefile changes: add the `carrizo` host, point `HOST_*CC*` at it, add it
to `ANNOT_HOSTS`, then `make HOSTSEL=carrizo fractions counts matrices`. Keep the sweet-orange-only
CC counts for comparison. Hamlin keeps DVS_A1.0.

**B. Baselines through the same pipeline.** Done for sweet orange and Carrizo (step 5c). With a
`carrizo` host, Carrizo baselines would go to it and the *Poncirus* runs to ZK8. SE datasets
count reads, not pairs; note it in the matrix.

**C. Annotation layer.** Step 7 covers the functional part. Still to add: ZK8 proteins (from its
GFF), the Cs <-> Pt pair table (MCscan anchors + RBH), which doubles as the parental-copy map, and
the orange1.1 -> LOC map for the PlantRegMap/PlantTFDB files.

**D. Gene-level analyses.**
1. *Gall vs baselines*, per host genotype. DESeq2, baseline replicates for dispersion, `study`
   as a covariate where it isn't aliased. Run against each baseline separately, then intersect:
   call a gene only if |log2FC| > 2 and consistent across >= 2 studies. Report within-sample
   percentile ranks next to fold changes. Tissue logic: up vs bark/stem/root/leaf **and** vs
   callus = gall-specific candidate; up vs stem but callus-like = proliferation program.
2. *Signature check*: score the Deeken up/down sets (§7.4) in each gall vs baseline. If citrus
   galls don't show the canonical program (photosynthesis down, fermentation/SuSy/STP up,
   CUT1/WAX2 down), suspect the baseline before the biology.
3. *Internal contrasts* (valid on the current sweet-orange-only counts): wt CC547 vs G-19 + G-30
   CC547 (T-DNA-dependent host genes; pending the §6.2 genotype check); CC54 vs CC547 wt (gall
   age); HC83 vs CC wt (host and strain confounded: descriptive only).
4. *Cross-host agreement*: genes moving the same way in the Hamlin gall and the wt Carrizo
   galls, each relative to its own genotype's baselines. The most defensible "bacterium-driven"
   set, given the design.
5. *Parental-copy expression (Carrizo only)*: for each paired gene, the *Poncirus*-copy share in
   gall vs healthy Carrizo tissue. A gene whose share shifts responds in cis. Count GLM with a
   copy x condition interaction, or beta-binomial per gene. With n = 1 galls this is a ranking
   exercise, but the healthy baselines supply the normal ratio and its variance.

**E. Promoters.**
1. *Parental-promoter comparison* (uses D5): for genes with a gall-induced copy shift, align the
   sweet-orange and *Poncirus* promoters (1-2 kb upstream of the TSS), scan both with JASPAR
   plant motifs (FIMO), and list motifs gained or lost in the responsive copy. Same nucleus, same
   bacterium, same hormones: a sequence difference that tracks the expression difference is the
   strongest promoter evidence this design can give.
2. *Host promoters of candidate genes*: foreground = the gene sets from D, background =
   expressed genes matched on expression level. Known-motif enrichment (MEME-suite SEA/AME or
   HOMER, JASPAR 2024 plants); expect auxin (AuxRE TGTCTC / ARF), type-B ARR cytokinin (AGATHY),
   ABRE, hypoxia (HRPE), W-box if the Arabidopsis program holds. PlantRegMap's precomputed TFBS
   tables are a shortcut. Motif instances conserved between the Cs and Pt promoters are much
   stronger evidence than enrichment alone. With n = 1, promoter results are hypotheses.
3. *The bacterium's own plant-active promoters*: the intergenic regions upstream of each
   expressed T-DNA gene in 1416 and 29 (§6.3); annotate TATA and ocs/as-1 elements (TGACG, bound
   by host TGA/bZIP factors) and compare between strains. These are also the promoter parts the
   GAANTRY work would reuse.

## 9. Bacterial side: the Haryono 2019 template (steps 8-9)

Not started. The closest published analog is Haryono M, Cho S-T, Fang M-J, Chen A-P, Chou S-J,
Lai E-M, Kuo C-H (2019) "Differentiations in gene content and expression response to virulence
induction between two *Agrobacterium* strains", *Front Microbiol* 10:1554
(doi:10.3389/fmicb.2019.01554; reads SRP156105; genomes C58 AE007869-AE007872, 1D1609
CP026924-CP026928). They compare C58 (nopaline-type) with 1D1609 (octopine-type), as we compare
1416 with 29. **Key difference:** pure cultured bacteria (97.7% mapping), no host, no dual
RNA-seq.

**What they did.** 2 strains x acetosyringone-induced vs control x 3 replicates; Ribo-Zero,
stranded, 101 bp SE, ~26 M reads per sample. BWA to each strain's own genome; NOISeqBio with
RPKM, DE = >= 2-fold and probability >= 0.99: 88/5,355 DE genes in C58, 155/5,630 in 1D1609.
OrthoMCL gave ~4,115 shared gene clusters and > 1,000 strain-unique genes. Conservation falls by
replicon: circular chromosome (94.6% AA similarity, 76% of genes shared) > linear chromosome
(92.0%, 59%) > pTi (71.7%, 31%) > pAt (3% shared).

**Step 8. Cross-strain orthology.** Strain 1416 + strain 29 proteins (from the `.gb`
annotations; optionally C58) -> OrthoFinder (conda, to install). Classify each ortholog pair by
DE status in the two strains (up/up, up/absent, ...; their Table 2). In their data fewer than
half of DE genes had a homolog with the same pattern: regulatory divergence, not just gene
content. Analyse per replicon.

**Step 9. Bacterial promoters and motifs.** 600 bp upstream of DE genes, grouped up/down per
strain. Targeted scan first: VirG box consensus RTTDCAWWTGHAAY, <= 3 mismatches (`fimo`); genes
< 50 bp apart on the same strand count as one operon, and the upstream gene's motif counts for
the downstream genes. Then `meme`/`streme` ab initio. Our libraries are unstranded, so
coverage-based TSS or operon inference isn't possible from this data.

**Expectation check.** Haryono et al. found the VirG box upstream of most up-regulated genes
(44/52 in C58, 61/74 in 1D1609), including genes off the Ti plasmid, and **no robust novel motif
ab initio**, even with clean cultures and replicates; results were parameter-sensitive. Novel
promoter discovery likely needs cross-strain phylogenetic footprinting (align upstream regions
of orthologs) or a dedicated TSS library (dRNA-seq/Cappable-seq).

**Biology to check in 1416/29.** The vir regulon dominates up-regulated genes, and most
down-regulated genes are chromosomal. vir gene locations vary between strains (1D1609 has
virA/virK/virH2 on pAt), so check where ours sit before assuming pTi. *tzs* is
nopaline-strain-specific and > 200-fold induced in C58; check 1416/29 (§6.3 notes an *ipt*-like
gene among the 1416 vir genes). Iron/zinc transport genes go down under induction in both
strains.

## 10. Findings so far

- Libraries are polyA mRNA and unstranded (§3).
- Carrizo is an F1 hybrid aligned to sweet orange only: 82-85% alignment vs 92% for Hamlin, host
  mismatch rate 0.88-1.02% vs 0.56%. Fine for Carrizo-vs-Carrizo comparisons, biased against
  other hosts (§6.2).
- The public Carrizo and sweet-orange baselines align at the same rates as the galls of the same
  genotype (85-87% and 89-94%), so gall-vs-baseline comparisons share the mapping bias (step 5c).
- Wild-type citrus galls express the bacterium's T-DNA; the opine-synthase-like genes are
  strongest (9,700-16,300 pairs per gall). G-19/G-30 galls barely express it, cluster apart from
  the wt galls, and are near-replicates of each other (§6.3).
- G-19/G-30 also look less hybrid-like, so their plants' genotype needs confirming (§6.2). Their
  bacterial reads per million (2.7 and 9.8) bound plant-to-bacterium cross-mapping (step 5b).
- Method 1's shortlist is dominated by constitutive and stress/wound genes (step 5d); 27 of the
  42 genes in the union of lists carry a multi-copy or cloning warning.
- No citrus crown-gall transcriptome has been published.

## 11. Decisions, and why

- **HISAT2 + featureCounts stays** for alignment, QC and gene-level counts. Salmon only if the
  per-parental-copy Carrizo analysis is done (§8.2 A).
- **Carrizo stays on the sweet-orange reference for now**, galls and baselines alike, so the bias
  is shared. The two-parent reference is planned, not required for the shortlist.
- **Shortlist method:** six simple ranking methods, top 20 each, compared by hand against
  `genes.tsv` (§8.1). Simple and established over sophisticated.
- **Annotation kept simple:** eggNOG-mapper plus DIAMOND best hits against Arabidopsis.
- **Mapping bias is documented, not corrected** (`make mapstats`), because it cancels within a
  host genotype.
- **All non-trivial compute runs as Slurm jobs**, not in an interactive shell.
- **Plotting stack in its own venv** (`~/.venvs/viz`), so the pinned aligner/counter versions in
  the `rnaseq` env can't move.

## 12. Next steps, in order

1. **Re-execute `notebooks/03_citrus_shortlist.ipynb`** now that
   `04_matrix/baseline_citrus_sinensis.tsv` exists: methods 3, 4 and 6 and
   `results/shortlist/merged_<group>.tsv`. First join tissue labels from `citrus_baselines.tsv`
   and sanity-check the baselines (§8.2 D2: photosynthesis genes should separate leaf from
   bark/root/callus).
2. **Back up to `/project/small_grains`** (§4) and decide whether to delete the 191 GB of
   baseline fastq + trimmed reads.
3. **Method 5 (RankComp)** in the notebook; **method 2** if R/DESeq2 gets installed.
4. **Merged top-20 table, manual review, ~1.5 kb promoter sequences** for the picks (§8.1).
5. **Ask the collaborators** (§13).
6. Housekeeping (§15): the two missing symlinks, the worktrees, the tomato strain contrast.

Optional: eggNOG `--dbmem` with >= 64 GB; finish tomato and run papaya annotation; tier-2
baselines; the two-parent Carrizo reference (§8.2 A), which unlocks the *Poncirus* baselines, the
genotype check and the parental-copy analyses; step 6 DE; steps 8-9.

## 13. Asks for the collaborators

- **Any uninoculated tissue from the same plants/greenhouse**: wounded-mock or adjacent healthy
  bark from the CC547 and HC83 plants, even 2-3 samples. Worth more than all the public baselines
  combined; it would turn §8.2 D1 from exploratory into a real contrast.
- G-19/G-30 construct maps: what replaced or modified the 1416 T-DNA.
- The CC plants: confirm Carrizo, seed source, whether they are nucellar seedlings (not grafted),
  and whether the G-19/G-30 plants came from the same stock as the wt ones (§6.2 QC flag).
- RNA extraction and library kit, to confirm the polyA/unstranded inference, and the tissue
  boundaries of the gall sample (gall only, or gall + adjacent stem?).
- From `SOURCES.md` §12: strain 29's designation and publication; the tomato cultivar (T29, T49);
  the *Euonymus* cultivars; the 1416 pAt1 size (589,645 bp in our conversion vs 519,735 bp in
  Alabed et al. 2023).

## 14. Setup and running

**On Ceres (current).** Partition `ceres`, account `small_grains`; home quota 30 GB, so everything
large is symlinked into 90daydata (§4). `local.mk` (gitignored, auto-included) sets the Ceres
module versions (`MOD_SAMTOOLS = samtools/1.17`, `MOD_SUBREAD = subread/2.0.4`; the Makefile
defaults are the Atlas ones) and is the place for `EGGNOG_DATA`. Submit from a login node, not
`ceres-dtn`. `SLURM.md` covers interactive and detached sessions on compute nodes.

| job | command |
|---|---|
| whole gall pipeline (48 cpu / 128 G / 12 h) | `sbatch scripts/pipeline.slurm [host ...]` |
| baselines (48 cpu / 128 G / 24 h) | `sbatch scripts/baselines.slurm` |
| annotation (16 cpu / 32 G / 8 h) | `sbatch scripts/annotate.slurm [hosts]` |
| a single target | `srun --cpus-per-task=16 --mem=64G make -j2 align` |
| mapping statistics | `srun --cpus-per-task=16 --mem=16G make -j4 mapstats THREADS=4` |
| notebooks | `scripts/viz-env.sh` once, then `~/.venvs/viz/bin/jupyter nbconvert --to notebook --execute --inplace notebooks/<nb>.ipynb` (`notebooks/README.md`) |

**On a new cluster.**
1. Put the repo where there is room for ~420 GB of outputs (or symlink the output dirs as on
   Ceres, including `05_baseline` and `references/plant_host/arabidopsis_thaliana`), and place
   or symlink the raw reads at `30-1348328766/00_fastq/`. Copy the gitignored collaborator
   folders (`Strain 1416/`, `Strain 29/`, `C58_ALIGNED/`, `docs/`) alongside.
2. `conda env create -f environment.yml` (`scripts/conda-env.slurm` does it as a job). Until it
   exists, `make strain-refs CONDA=$(dirname $(which python3))` works after `module load python_3`.
   On Ceres, conda fails while linking files into home, which is why the notebook and
   eggNOG-mapper environments are plain venvs.
3. `module avail hisat2 samtools subread fastqc diamond`; put differing names in `local.mk`.
4. `make verify qc`, `make refs` (~3.6 GB of host genomes), then the batch jobs above.

Not bundled: the collaborators' Snippy variant calls and Bakta annotations for strain 29
(~5.9 GB, still on Atlas under `rnaseq/Strain 29/`). Prior comparative work, not needed here.

## 15. Gotchas and housekeeping

- **Two symlinks are missing from the main checkout.** `05_baseline` and
  `references/plant_host/arabidopsis_thaliana` exist on 90daydata and in the `baselines`
  worktree, but not in `~/gallrna`. Running `make baselines` or `make annotate` from the main
  checkout would write into home (30 GB quota) and re-download. Create both symlinks first.
  `.gitignore` lists `05_baseline/` with a trailing slash, which doesn't match a symlink: add
  plain `05_baseline` (as was done for `qc`).
- **`logs/baseline_summary.tsv` and the baselines Slurm logs live in
  `.claude/worktrees/baselines/logs/`**, not in the main `logs/`. Copy them before removing that
  worktree.
- **Unmerged work in the `baselines` worktree** (untracked, in no branch):
  `scripts/strain_contrast.py` (two galls of different strains on one host: TPM, percentile and
  log2 ratio per gene; a ranking, not a test) and its output
  `results/shortlist/tomato_strain_contrast.tsv` (29wtT29 vs 1416wtT49, 28,485 genes), plus a
  rendered `notebooks/03_citrus_shortlist.html`. Commit it to `main` or discard it deliberately.
- **Leftover worktrees** under `.claude/worktrees/`: `annotate-target`, `baselines`,
  `carrizo-makefile-comment`, `carrizo-plan-merge`, `citrus-host-plan`, `handoff-doc`,
  `sources-file` (on branch `mapstats`). Every branch is fully merged into `main` except local
  `worktree-sources-file` (one commit, `b0860c7`, which looks superseded by `4840030`; check
  before deleting). `carrizo-makefile-comment` has an uncommitted Makefile edit that `main`
  supersedes. The worktrees still hold copies of the old docs this file replaces.
- **Notebooks 01/02 will break on Carrizo.** `notebooks/gallrna_viz.py:32` `HOST_ORDER` still
  lists "Poncirus trifoliata", but `samples.tsv` calls the CC host "Citrus sinensis x Poncirus
  trifoliata", so the colour lookup fails. Rename that entry.
- **`notebooks/03_citrus_shortlist.ipynb` has an uncommitted re-execution** from the main
  checkout (outputs only; it predates the baseline matrix). Step 1 of §12 supersedes it.
- **Download hosts.** `eggnogdb.embl.de` no longer resolves, so the Makefile fetches the eggNOG
  database from `eggnog5.embl.de`. A compute node once failed to resolve
  `ftp.ensemblgenomes.ebi.ac.uk`; fetch the Arabidopsis proteome on the login node if that recurs
  (Makefile comment).
- **Timestamps.** NCBI zips can carry future timestamps (the protein rule touches the file). If
  make unexpectedly wants to rerun eggNOG-mapper, check that `protein.faa` < `longest.faa` <
  `emapper.annotations`. Likewise a fresh checkout's script timestamps can make `make shortlist`
  want to rebuild the gall chain; run the script directly.
- The `arabidopsis_thaliana` directory on 90daydata was deleted once while still empty (around
  19:09 on 2026-09-14, cause unknown), which failed a job. It now holds the proteome.
- **`environment.yml` lists more than is installed or used** (STAR, bowtie2, StringTie, Salmon,
  MEME, DESeq2/edgeR). R/DESeq2 is not installed on Ceres.
- **Open discrepancies** (`SOURCES.md` §12): the 1416 pAt1 size, how González-Mula et al. 2018
  split reads, and the collaborator unknowns (§13).

## 16. Other docs

| file | contents |
|---|---|
| `Makefile.annotated.md` | reading copy of the Makefile: the make syntax it uses, and input -> output per block. Not executed; goes stale when the Makefile changes |
| `SOURCES.md` | citation and justification for every dataset and tool; §12 lists open discrepancies |
| `SLURM.md` | running Claude Code and the pipeline on Ceres compute nodes |
| `references/README.md` | what's in `references/`, replicon by replicon and host by host |
| `references/plant_host/AVAILABILITY.md` | genome availability for the collaborators' full 31-species host panel (2026-09-08) |
| `notebooks/README.md` | running the notebooks, and why they use a separate venv |

For Claude sessions: memory notes in `~/.claude/projects/-home-andrew-dickson-gallrna/memory/`
(run compute through Slurm; keep methods simple and aimed at the shortlist).
