# Handoff: gallrna (2026-09-15)

> State of the repo and data on 2026-09-15, for whoever picks this up next (person or Claude
> session). Start here, then follow the pointers.

## The project

Dual RNA-seq of *Agrobacterium* galls: 14 gall samples from two strains (1416, 29) on six host
species. Each sample is aligned to one combined reference (its strain + its host), and reads are
split by organism afterwards. Current focus: the **5 citrus galls**, one Hamlin sweet orange
(29wtHC83) and four Carrizo citrange (1416 wt at 54 and 547 days, plus engineered strains G-19
and G-30).

**Goal (set 2026-09-14):** give the collaborator a **shortlist of promoters/genes to target in a
follow-up experiment** (CITRUS_HOST_PLAN.md §4 F). No gall has an uninoculated control, so galls
are compared against public healthy-tissue RNA-seq of the same genotype.

## Where things live

| what | where |
|---|---|
| code and docs | this repo; GitHub `amdson/gallrna`, branch `main` |
| raw reads | `30-1348328766/` -> `/90daydata/small_grains/30-1348328766` (group dir, 60 GB) |
| pipeline outputs | `01_trim 02_align 03_counts 04_matrix 05_baseline qc references/combined references/plant_host/<host>`, all symlinks into `/90daydata/small_grains/andrew.dickson/gallrna_data/` (~160 GB: trim 52, align 42, references 66 incl. the 48 GB eggNOG database, qc 1.8; plus the public baselines under `05_baseline`, ~105 GB raw + trimmed + BAM once step 5c finishes) |
| eggNOG-mapper database | `/90daydata/small_grains/andrew.dickson/gallrna_data/references/eggnog_data` |
| collaborator files (gitignored, unpublished) | `Strain 1416/`, `Strain 29/`, `C58_ALIGNED/`, `docs/` |

**90daydata deletes files not accessed for 90 days.** Copy anything worth keeping (count
matrices, `logs/mapping_stats.tsv`, the `*.genes.tsv` tables) to `/project`.

Cluster: Ceres (partition `ceres`, account `small_grains`). `local.mk` (per-cluster, not in git)
sets the Ceres module versions. Home has a 30 GB quota.

## Which doc answers what

| file | contents |
|---|---|
| `PIPELINE.md` | steps 0-9: what each does, status, how to run; 5b explains the mapping-QC table, 5c the public baselines |
| `CITRUS_HOST_PLAN.md` | the citrus analysis: data, hybrid hosts, baselines, the §4 F shortlist method, next steps (§5), collaborator asks (§6) |
| `SOURCES.md` | citation and justification for every dataset and tool; §12 lists open discrepancies |
| `references/README.md` | what's in `references/`, host by host |
| `samples.tsv` | the 14 samples; CC = Carrizo citrange, HC = Hamlin |
| `citrus_baselines.tsv` | 42 public baseline runs (tier 1 and 2), each checked against NCBI/ENA labels |
| `notebooks/03_citrus_shortlist.ipynb` | the shortlist work off the cluster: method 1 with figures and interpretation, methods 3/4/6 once the baselines exist (viz venv, `notebooks/README.md`) |
| `HARYONO2019_PIPELINE.md` | template paper for the bacterial steps 8-9 |
| `SETUP.md` | new-cluster setup (the committed copy is older than the local edits; see below) |

## Status

| step | status |
|---|---|
| 0-5 verify -> count | done for all 14 samples (2026-09-08). Matrices: `04_matrix/agro_strain_{1416,29}.tsv`, `plant_{citrus_sinensis,carica_papaya,solanum_lycopersicum}.tsv`. *Euonymus* and *B. juncea* have no gene models: mapping fractions only |
| 5b `make mapstats` | done 2026-09-14 -> `logs/mapping_stats.tsv`, the per-sample table for the paper. Use it, not `logs/mapping_fractions.tsv` |
| 5c `make baselines` | **running**: job 22046668 submitted 2026-09-15 from the `worktree-baselines` worktree (ENA download -> fastp -> HISAT2 plant-only -> featureCounts for the 27 tier-1 sweet-orange/Carrizo runs, ~105 GB). Output: `04_matrix/baseline_citrus_sinensis.tsv` (on 90daydata, visible from any checkout) and `logs/baseline_summary.tsv` plus the Slurm log `logs/baselines-22046668.out` in that worktree's `logs/`. Check with `squeue -u $USER`. *Poncirus* runs skipped until ZK8 is annotated |
| 6 differential expression | not started; R/DESeq2 not installed |
| 7 `make annotate` | sweet orange done 2026-09-15 -> `references/plant_host/citrus_sinensis/citrus_sinensis.genes.tsv`: 28,080 genes; of 23,556 protein-coding, 97% eggNOG-annotated, 49% with GO, 94% with an Arabidopsis hit, 56% reciprocal best hits. Papaya and tomato not run; *Poncirus* ZK8 not supported yet |
| 8-9 bacterial orthology, promoters | not started |
| citrus shortlist (plan §4 F) | method 1 (gall expression) done 2026-09-15 -> `results/shortlist/method1_expression_*_top20.tsv`; methods 2-6 wait for the baselines |

## What we've learned

- Libraries are **polyA mRNA** (~1% intronic reads) and **unstranded**, so `STRAND=0` is correct.
- **Carrizo is an F1 hybrid** (Washington navel x *Poncirus*), aligned to sweet orange only:
  82-85% overall alignment vs 92% for Hamlin; host mismatch rate 0.88-1.02% vs 0.56%. That's
  fine for Carrizo-vs-Carrizo comparisons and biased against other hosts.
- **G-19/G-30 galls barely express the wild-type T-DNA** and cluster apart from the wt galls.
  They also look less hybrid-like (fewer mismatched reads), so confirm their plants' genotype.
  Their bacterial reads per million (2.7 and 9.8) put an upper bound on plant reads wrongly
  assigned to the bacterium.
- Wild-type citrus galls express the bacterium's T-DNA; the opine-synthase-like genes are the
  strongest (9,700-16,300 reads per gall).
- No citrus crown-gall transcriptome has been published.
- When scoring the Arabidopsis crown-gall markers, match gene families, not just best hits:
  citrus GH3 genes, for example, hit other Arabidopsis GH3 paralogs.

## Decisions, and why

- **HISAT2 + featureCounts stays** for alignment, QC and gene-level counts. Salmon only if the
  per-parental-copy Carrizo analysis is done (plan §4 A).
- **Carrizo stays on the sweet-orange reference for now.** A two-parent (sweet orange + *Poncirus*
  ZK8) `carrizo` reference is planned (plan §1.1, §4 A).
- **Shortlist method:** six simple ranking methods (gall expression, shrunken fold change,
  percentile rank shift, rank products, RankComp, tau), top 20 from each, compared by hand
  against `genes.tsv` (plan §4 F).
- **Annotation kept simple:** eggNOG-mapper plus DIAMOND best hits against Arabidopsis.
- **All non-trivial compute runs as Slurm jobs**, not in an interactive shell.

## Next steps, in order

1. **Shortlist method 1: done 2026-09-15** (`make shortlist`, PIPELINE.md 5d) ->
   `results/shortlist/method1_expression_{hamlin,carrizo_wt,carrizo_eng,all}_top20.tsv` and the
   full table `04_matrix/citrus_gall_expression.tsv`. Top of the lists = constitutive promoters
   (ribosomal, ubiquitin, cyclophilin, metallothionein) plus stress/wound genes (MLP423-like,
   LEA5, dehydrin, extensins, aquaporins in Hamlin). Review by hand once methods 2-6 exist.
2. **`make baselines`** exists (PIPELINE.md 5c) and the tier-1 job is running (status table).
   When it finishes: read `logs/baseline_summary.tsv` (alignment rate and assigned fraction per
   run; the Carrizo runs should sit near the CC galls' 82-85%, the sweet-orange runs near 92%),
   then copy `04_matrix/baseline_citrus_sinensis.tsv` to `/project`. Add `05_baseline` as a
   symlink into 90daydata in the main checkout (the job ran from the worktree, which has one).
3. **Methods 2-6**, then the merged top-20 table, manual review, and ~1.5 kb promoter sequences
   for the picks.
4. **Ask the collaborators** (plan §6; SOURCES.md §12 item 5): mock tissue, G-19/G-30 construct
   maps, the CC plants' source, the library kit, strain 29's designation.

Optional: add eggNOG-mapper `--dbmem` and >=64 GB to `make annotate` (its annotation phase took
5.3 h at 32 GB); annotate papaya and tomato; build the two-parent Carrizo reference.

Run commands:
- `sbatch scripts/annotate.slurm [hosts]` with `EGGNOG_DATA` pointing at the 90daydata path
  (environment or `local.mk`)
- `make mapstats` and other heavy targets: on a compute node (`srun ... make`) or a batch job

## Gotchas and open issues

- **Uncommitted work in the main checkout.** `scripts/convert_refs.py` writes the `gene` rows the
  bacterial GFF3s need, and the existing bacterial counts depend on it. Also uncommitted:
  `.gitignore` (ignores the symlinked output dirs, `local.mk`, `figures`), `SETUP.md`,
  `environment.yml` (+biopython). Untracked, but referenced by the docs: `SLURM.md`,
  `notebooks/`, `scripts/pipeline.slurm`, `scripts/conda-env.slurm`, `scripts/claude-*`,
  `scripts/intron_fraction.sh`, `scripts/splice_rates.sh`, `scripts/viz-env.sh`. Commit them
  when ready.
- **Notebooks will break on Carrizo.** `notebooks/gallrna_viz.py` `HOST_ORDER` still lists
  "Poncirus trifoliata", but `samples.tsv` now calls the CC host "Citrus sinensis x Poncirus
  trifoliata", so the colour lookup fails. Rename that entry.
- **Download hosts.** `eggnogdb.embl.de` no longer resolves, so the Makefile fetches the eggNOG
  database from `eggnog5.embl.de`. A compute node failed to resolve
  `ftp.ensemblgenomes.ebi.ac.uk`; fetch the Arabidopsis proteome on the login node if that
  recurs (Makefile comment).
- **Timestamps.** NCBI zips can carry future timestamps (the protein rule now touches the file).
  If make unexpectedly wants to rerun eggNOG-mapper, check that `protein.faa` < `longest.faa` <
  `emapper.annotations`.
- The empty `.../plant_host/arabidopsis_thaliana` dir on 90daydata was deleted once (around
  19:09 on 2026-09-14, cause unknown), which failed a job. It now holds the proteome.
- **Leftover worktrees** under `.claude/worktrees/`: `annotate-target`, `carrizo-plan-merge`,
  `citrus-host-plan`, `sources-file` (branch `mapstats`), `handoff-doc`, and
  `carrizo-makefile-comment` (another session's; its uncommitted Makefile edit is superseded by
  `main`). All their work is in `main`. **Remote branches** `mapstats`, `worktree-annotate-target`,
  `worktree-carrizo-plan-merge` and `worktree-citrus-host-plan` are merged or superseded;
  `worktree-sources-file` (`b0860c7`) looks superseded by `4840030`, so check before deleting.
- **SOURCES.md §12** still open: the 1416 pAt1 size, the C58 accession list, and the collaborator
  unknowns.

## For Claude sessions

Memory notes in `~/.claude/projects/-home-andrew-dickson-gallrna/memory/`: run compute through
Slurm, and keep methods simple and established, aimed at the shortlist.
