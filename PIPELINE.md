# Gall Dual RNA-seq Pipeline

> Steps 0-5 and 7 are automated in `Makefile` (`make help`); sample metadata in `samples.tsv`.
> Run heavy targets on a compute node, e.g. `srun --cpus-per-task=16 --mem=64G make -j2 align`,
> or as a batch job (`sbatch scripts/annotate.slurm` for step 7). Keep the login node for editing.
> Citrus-host analysis plan (baselines, hybrids, promoters): CITRUS_HOST_PLAN.md.

## Existing data

`30-1348328766/00_fastq/` — Azenta/GeneWiz, Illumina NovaSeq X, 2x150 PE, 14 samples (R1+R2 + md5 each).
Sample naming: `<strain><genotype><host code><gall age in days>`. Gall tissue: mixed plant host +
Agrobacterium reads. Libraries are **polyA mRNA** (~1% intronic reads) and **unstranded**
(sense:antisense 50:50), both confirmed 2026-09-14.

**Strain 1416 (8 samples)**
- wt on hosts: CC54, CC547, Eu635, Geu182, M26, T49
- engineered events G-19, G-30 on host CC547 (wt/G-19/G-30 trio; n=1 per genotype — limits DE inference).
  G-19/G-30 galls express almost none of the wt 1416 T-DNA (CITRUS_HOST_PLAN.md §1.2)

**Strain 29 (6 samples, all wt)**
- hosts: Eu635, Geu182, HC83, M26, P46, T29

**Host codes** (species per sample in `samples.tsv`; reference genome per host in the Makefile):

| code | host | reference | gene models |
|---|---|---|---|
| Eu, Geu | *Euonymus japonicus* (green, golden) | *E. europaeus* proxy GCA_963580455.1 | none |
| HC | *Citrus sinensis* 'Hamlin' | DVS_A1.0 GCF_022201045.2 | RefSeq |
| CC | **Carrizo citrange** (*C. sinensis* x *P. trifoliata* F1) | DVS_A1.0 for now (see step 3) | RefSeq (sweet-orange half only) |
| M | *Brassica juncea* | GCA_018703725.1 | none |
| P | *Carica papaya* | Papaya1.0 GCF_000150535.2 | RefSeq |
| T | *Solanum lycopersicum* | SLM_r2.1 GCF_036512215.1 | RefSeq |

**Strain assemblies**: converted to FASTA/GFF3 on 2026-09-08 (`references/README.md`):
- Strain 1416: circular + linear chromosome, pAT1, pAT2, pTi (all FINAL.gb)
- Strain 29: circular + linear chromosome, pAT1, pAT2, pTi (RAST annotation)
- C58 published reference (AE007869/AE007870/AE007872) — for comparison only

**Still needed**
- Uninoculated (mock/healthy) tissue from the same plants — no host has a within-experiment
  control; citrus uses public baselines instead (CITRUS_HOST_PLAN.md §3.3)
- G-19 / G-30 construct or event sequences (check "1416 10294 bp fragment (LC).gb")
- Gene models for the Euonymus proxy and *B. juncea*: these samples have mapping fractions but
  no plant count matrix

## Pipeline

### 0. Verify transfer
- Data: fastq.gz + .md5 files
- Tools: `md5sum -c`

### 1. Read QC
- Data: raw fastq.gz
- Tools: `fastqc` (module), `multiqc` (conda env `rnaseq`)
- Check: adapter content, rRNA overrepresentation, per-base quality

### 2. Trim/filter
- Data: raw fastq.gz
- Tools: `fastp` (conda env `rnaseq`)
- Output: trimmed fastq + per-sample QC json/html

### 3. Build combined references
- Data: plant host genome FASTA + strain FASTA (chromosomes + plasmids), merged GFF3
- Tools: Biopython (`scripts/convert_refs.py`) for .gb -> FASTA/GFF3, `hisat2-build`
- Strain contigs carry an `agro_` prefix; plant contigs keep their NCBI accessions, so there are
  no name collisions and the organism split downstream keys on `agro_`
- One combined index per strain x host pair actually sampled (`references/combined/<strain>__<host>`)
- Carrizo (CC) is aligned to sweet orange only for now: fine for CC-vs-CC contrasts, biased
  against other hosts. Planned: a two-parent `carrizo` host (DVS_A1.0 + *Poncirus* ZK8), see
  CITRUS_HOST_PLAN.md §1.1 and §4 A
- Still to add: G-19/G-30 event sequences as extra contigs for those samples

### 4. Align (competitive mapping)
- Data: trimmed fastq + combined index for the sample's strain x host
- Tools: `hisat2 --dta` (module); sorted, indexed BAM (`samtools`, module)
- QC: `samtools idxstats` aggregated by contig prefix -> plant vs agro fraction per sample
  (`logs/mapping_fractions.tsv`)

### 5. Count
- Data: BAMs + merged GFF3
- Tools: `featureCounts -p --countReadPairs -s 0 -t gene -g ID` (subread module). `-s 0` is
  correct: libraries are unstranded
- `scripts/merge_counts.py` splits by the `agro_` prefix -> `04_matrix/agro_strain_<strain>.tsv`
  and `04_matrix/plant_<host>.tsv` (annotated hosts only)
- Flag T-DNA genes (multimap between pTi and transformed plant genome)

### 6. Differential expression
- Data: split count matrices + sample metadata (strain, genotype, host)
- Tools: R in conda env `rnaseq`: DESeq2, tximport, edgeR, pheatmap, ggplot2 (not installed yet)
- Normalize each organism SEPARATELY (own size factors — plant/agro ratio varies per sample)
- Key contrasts: wt vs G-19/G-30 on host CC547 (n=1 each — fold-change ranking only, or treat
  G-19+G-30 as 2 reps of "engineered"); strain 1416 vs 29 on shared hosts (Eu635, Geu182, M26)
- Citrus host genes: gall vs public healthy-tissue baselines of the same genotype, internal
  contrasts, crown-gall signatures and parental-copy expression — CITRUS_HOST_PLAN.md §4 D.
  Baseline run list: `citrus_baselines.tsv`

### 7. Host functional annotation (`make annotate`)
- Data: host RefSeq proteins (NCBI PROT_FASTA), Arabidopsis proteome (Ensembl Plants release 63),
  eggNOG-mapper database emapperdb 5.0.2 (~50 GB, from eggnog5.embl.de)
- Tools: eggNOG-mapper 2.1.13 (venv from `scripts/emapper-env.sh`), `diamond` (module),
  `scripts/longest_proteins.py`, `scripts/gene_annotation.py`
- Longest protein per gene, named by the gene ID featureCounts writes -> eggNOG-mapper
  (Viridiplantae scope: GO, KEGG, EC, Pfam, description) + DIAMOND best hit each way vs
  Arabidopsis (AGI, symbol, reciprocal-best-hit flag)
- Output: `references/plant_host/<host>/<host>.genes.tsv`, one row per counted gene; joins onto
  `04_matrix/plant_<host>.tsv` by gene ID. The Arabidopsis columns connect to the AGI-based
  crown-gall marker lists (CITRUS_HOST_PLAN.md §3.4)
- Run: `sbatch scripts/annotate.slurm` (citrus_sinensis by default; pass other RefSeq hosts as
  the first argument). Point `EGGNOG_DATA` at 90daydata (environment or `local.mk`) — the
  database would overflow the 30 GB home quota. Sweet orange: 40,427 RefSeq proteins ->
  23,556 genes
- Not covered yet: *Poncirus* ZK8 (needs proteins from its CGD GFF) and the unannotated hosts

### 8. Cross-strain orthology (after Haryono et al. 2019 — see HARYONO2019_PIPELINE.md)
- Data: strain 1416 + strain 29 protein sequences (from .gb annotations); optionally C58
- Tools: `orthofinder` (conda, to install when needed)
- Cluster homologs; classify each ortholog pair by DE status (up/up, up/absent, etc.)
- Analyze per replicon — expect circular chr most conserved, linear chr + plasmids divergent

### 9. Promoter / motif analysis (bacterial genes)
- Data: 600 bp upstream regions of DE genes (grouped up/down per strain). Libraries are
  unstranded, so coverage-based TSS/operon inference is not possible from this data
- Tools: MEME suite (conda env `rnaseq`): `fimo` for known motifs, `meme`/`streme` for ab initio
- Targeted scan first: VirG box consensus RTTDCAWWTGHAAY, <=3 mismatches; group genes <50 bp
  apart on same strand as one operon (Haryono et al. protocol)
- EXPECTATION CHECK: Haryono et al. found NO robust novel motif ab initio even with clean
  cultures + replicates — only the known VirG box. Novel promoter discovery likely needs
  cross-strain phylogenetic footprinting (align upstream regions of orthologs) or a dedicated
  TSS library (dRNA-seq/Cappable-seq)
- Host-side promoters (citrus genes the gall changes, and the T-DNA's plant-active promoters):
  CITRUS_HOST_PLAN.md §4 E
