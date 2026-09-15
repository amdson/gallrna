# Citrus host response: which host genes/promoters did the gall bacterium change?

> Scoping doc, 2026-09-14; revised the same day once the CC host was confirmed as **Carrizo
> citrange** (merges the reasoning from the CC comment in the Makefile). Covers the 5
> citrus-host samples only. Automated so far: host functional annotation (`make annotate`);
> §5 lists the other proposed Makefile additions. Public accessions were checked against
> NCBI/ENA/CGD on this date; anything not verified is marked **(unverified)**.
>
> **Goal:** the collaborator wants a shortlist of promoters/genes to target in a follow-up
> experiment. §4 F is the path to that list; the rest of §4 is supporting analysis for
> explaining why a candidate responds, and is optional.

## 1. What we have

| sample | host | strain | gall age (d) | assigned pairs | genes >=10 counts |
|---|---|---|---|---|---|
| 29wtHC83 | *Citrus sinensis* 'Hamlin' | 29 wt | 83 | 26.0 M | 18,239 |
| 1416wtCC547 | Carrizo citrange | 1416 wt | 547 | 21.3 M | 17,023 |
| 1416wtCC54 | Carrizo citrange | 1416 wt | 54 | 20.8 M | 17,132 |
| 1416G-19CC547 | Carrizo citrange | 1416 G-19 | 547 | 22.0 M | 17,239 |
| 1416G-30CC547 | Carrizo citrange | 1416 G-30 | 547 | 20.7 M | 16,753 |

Counts: `04_matrix/plant_citrus_sinensis.tsv` (all 5 currently aligned to and counted on sweet
orange; see §1.1 for why that is acceptable for some contrasts and not others).

Library facts established 2026-09-14 (needed to match public data):
- **polyA mRNA**: intronic reads are ~1% of located reads (`logs/intron_fraction.tsv`).
- **Unstranded**: featureCounts `-s 1` vs `-s 2` on a 2.6 M-read slice of 29wtHC83 gave
  250,875 vs 249,061 assigned. `STRAND=0` is correct; keep it for everything below.
- 1416G-30CC547 has 22.0 M multimapping pairs vs 4.6-9.2 M in the other four — check what
  (rRNA? repeat? construct?) before trusting it in contrasts.

**There is no uninoculated citrus tissue.** Every sample is a gall and every condition is n=1.
"What did the bacterium change" therefore has to be triangulated (§4) rather than read off a
gall-vs-mock contrast.

### 1.1 Both citrus hosts are hybrids

**CC = Carrizo citrange** (confirmed 2026-09-14): *C. sinensis* 'Washington' navel x
*P. trifoliata*, a single F1 clone (made 1909 by Swingle's USDA program; Carrizo and Troyer are
the same clone). Carrizo is normally raised from nucellar (clonal) seedlings, so all four CC
plants should share one genotype. There's no Carrizo genome assembly. Every Carrizo gene has
one sweet-orange copy and one *Poncirus* copy.

**Hamlin is a hybrid too**: sweet orange carries mandarin and pummelo ancestry, and all sweet
orange cultivars (Hamlin, Valencia, Washington navel) are somatic mutants of one original
hybrid. So DVS_A1.0 (Valencia haplotype A) plus DVS_B1.0 (GCA_022201065.1, haplotype B)
together describe Hamlin almost exactly. They also describe the sweet-orange half of Carrizo,
which is a recombined set of those two haplotypes.

How it shows up in the current alignments (one sweet-orange chromosome, NC_068563.1; first
400k primary alignments with MAPQ >= 10; `NM` = mismatches to the reference):

| sample | overall alignment | NM = 0 | NM 1-2 | NM >= 3 | MAPQ < 10 |
|---|---|---|---|---|---|
| 29wtHC83 (sweet orange) | 92.1% | 66.2% | 27.1% | 6.7% | 2.4% |
| 1416wtCC547 | 84.5% | 40.5% | 42.5% | 17.0% | 5.8% |
| 1416wtCC54 | 85.2% | 39.7% | 44.1% | 16.1% | — |
| 1416G-19CC547 | 81.7% | 45.1% | 42.3% | 12.7% | — |
| 1416G-30CC547 | 85.0% | 54.7% | 35.8% | 9.5% | — |

Consequences:
- Aligning Carrizo to sweet orange alone loses or mis-scores the *Poncirus* copies. This is what
  published Carrizo RNA-seq does. Afzal Naveed, Huguet-Tapia & Ali 2019, *J Plant Interact*
  14:187-204 (doi:10.1080/17429145.2019.1609106) mapped Carrizo roots to two *C. sinensis*
  genomes. They report only 55-73% transcriptome coverage and fell back to de novo assembly.
  Here the CC libraries align at 82-85% vs 92% for the native sweet-orange gall. A
  cross-species test (2026-09-08) assigned 92% as many pairs to genes as that native sample.
- All four CC galls share the hybrid genotype, so that mapping bias is common to every CC sample
  and **cancels in CC-vs-CC contrasts**. It does **not** cancel against galls on other hosts,
  against pure sweet-orange or pure *Poncirus* baselines, or in any per-parent question.
- **QC flag**: G-30 (9.5% NM >= 3) and G-19 (12.7%) look less hybrid-like than the two wt galls
  (16-17%). A different mix of expressed genes can shift this, so it isn't proof. Before relying
  on wt-vs-G-19/G-30, check that all four plants are Carrizo: at genes where sweet orange and
  *Poncirus* differ, an F1 should split ~50:50 between the parental copies (falls out of §4 A).

### 1.2 The bacterium's direct footprint: T-DNA expression in the galls
T-DNA transcripts are made by transformed *plant* cells from plant-active bacterial promoters,
so they are the one part of "which promoters did the bacterium influence" we can see directly.
Counts from `04_matrix/agro_strain_*.tsv`:

| strain 1416 T-DNA gene | CC547 wt | CC54 wt | G-19 | G-30 |
|---|---|---|---|---|
| "D-octopine dehydrogenase" (opine synthase-like) | 11,761 | 16,287 | 12 | 6 |
| agrocinopine synthase (*acs*) | 1,184 | 2,793 | 2 | 0 |
| indoleacetamide hydrolase (*iaaH*) | 507 | 711 | 4 | 0 |
| tryptophan monooxygenase (*iaaM*) | 8 | 21 | 0 | 0 |
| E protein / gene 5 / C protein / Atu6012 | 38-953 | 588-953 | 0 | 0-1 |

| strain 29 T-DNA gene (pTi ~78.8-94.2 kb) | HC83 |
|---|---|
| hypothetical peg.5741 (likely opine synthase — confirm by BLAST) | 9,656 |
| *iaaH* | 515 |
| *ipt* (cytokinin) | 469 |
| *iaaM* | 11 |

Observations:
- Both wt strains express their T-DNA in citrus. G-19/G-30 galls essentially do not, and they
  cluster apart from the wt galls (PC1, 55% of variance; `figures/10_citrus_structure.png`).
  The resource-paper draft (`docs/050136 ...docx`) describes disarming 1416 and installing
  GAANTRY, which would explain it — **need the G-19/G-30 construct maps** to be sure.
- No *ipt* is annotated inside the 1416 T-DNA window. The only 1416 *ipt* is at pTi
  164.8 kb among the vir genes, which is where *tzs* sits. Check whether 1416's T-DNA *ipt* is
  missing or just unannotated (one of peg 5568-5570?).
- *iaaM* is barely detected in both strains while *iaaH* is well expressed. Could be biology,
  or a mapping/annotation artefact — worth a look in IGV.

## 2. Framing

Without mock tissue, independent lines of evidence each carry a known bias. A host gene is a
good "influenced by the bacterium" candidate when they agree:

1. **Gall vs public healthy tissue of the same genotype** (§3.3): Carrizo baselines for the CC
   galls, sweet-orange baselines for Hamlin. Captures the gall program but is confounded with
   lab, tissue and age. Mitigated by using several baselines and several tissues.
2. **Internal contrasts** that hold host genotype constant. The best one is wt CC547 vs
   G-19/G-30 CC547: same host and age, with vs without wt T-DNA expression. These work on any
   reference (§1.1).
3. **Prior crown-gall signatures** (Arabidopsis, via orthologs), which say which direction
   known gall genes should move.
4. **Parental-copy (allele-specific) expression in Carrizo.** Both copies of each gene share
   one nucleus, so a gall-induced shift in the sweet-orange : *Poncirus* ratio has to come from
   the gene's own DNA (cis), most likely its promoter. This is the most direct route to
   promoters the design allows (§4 D5, E1).

## 3. Data available

### 3.1 Host references and annotation

**Sweet orange (Hamlin) — keep the current reference.**
- NCBI RefSeq GCF_022201045.2 (DVS_A1.0, Valencia, haplotype A), Annotation Release 103.
  Already in `references/plant_host/citrus_sinensis/`. 23,556 protein-coding genes, but ~7,800
  genes are "uncharacterized"/"hypothetical" and there is no GO/KEGG in the GFF, so functional
  annotation has to be added (§3.2).
- Optional later: add DVS_B1.0 (GCA_022201065.1, Clemson annotation, 27,940 genes) as the second
  sweet-orange haplotype. Low priority for gene-level counts.
- Newer assemblies (T2T SWO GCA_046470765.1/775.1, Pera Rio, Washington Navel, Newhall, ...)
  have no NCBI annotation, so there's no reason to switch.

**Carrizo — build a two-parent reference: sweet orange DVS_A1.0 + *Poncirus* ZK8.**
This follows the "concatenate both parents' genomes" note in
`references/plant_host/AVAILABILITY.md` and the "add *P. trifoliata* as a second haplotype"
option in the Makefile CC comment. Until it exists, the CC galls stay on sweet orange (valid for
CC-vs-CC contrasts only, §1.1).
- *P. trifoliata* ZK8 v1.0 (HZAU; Huang et al. 2021, *Hortic Res*, doi:10.1038/s41438-021-00505-2,
  BioProject PRJNA554539 — the paper CGD cites; Peng et al. 2020 is a different, UF/JGI assembly).
  Same assembly as the unannotated NCBI copy we already downloaded (GCA_018350135.1), but the
  Citrus Genome Database (CGD) serves the gene models too:
  - `https://www.citrusgenomedb.org/jb2/data/Ptri_ZK8_v1.sorted.gff.gz` (7.5 MB; 25,680 genes,
    39,675 mRNAs, with UTRs; IDs like `Pt1g002240`)
  - `https://www.citrusgenomedb.org/jb2/data/Ptri_ZK8_v1.fasta.gz` (86 MB)
  - CGD lists annotation BUSCO 93.4%, assembly BUSCO 98.6%.
  - Contig names differ: CGD uses `chr1_ZK8` ... `chrUn_ZK8`, NCBI uses `CM031384.1` ... plus
    `VKKW01*` scaffolds. **Use the CGD FASTA with the CGD GFF** rather than renaming. (The
    Makefile comment's `P.trifoliata_ZK8_v1.scaffolds.fa` wasn't found on CGD; the file listed
    above is the one CGD's genome browser serves.)
  - These are the genome-browser data files, not a formal download page (CGD's `/download` page
    and the Analysis page did not link them). They are publicly readable. Cite Huang et al. 2021.
- Contig names (`NC_0685xx` vs `chrN_ZK8`) and gene IDs (`gene-LOC...` vs `Pt...`) don't
  collide, so the two can be concatenated as they are, alongside the `agro_` strain contigs.
- ZK8 is a specific *Poncirus* accession, not Carrizo's actual pollen parent. *Poncirus* is
  low-diversity, so expect a modest mismatch rate on the *Poncirus* side.
- JGI *P. trifoliata* v1.3.1 (Phytozome/CGD `bio_data/1312872`): the page returned 403 and the
  access terms are **(unverified)**. Not needed given ZK8.

**Sweet orange <-> Poncirus gene map (now essential: it pairs the two parental copies)**
- CGD MCscan anchors `jb2/data/mcscan/Csin_DVS_A_v1.Ptri_ZK8_v1.anchors.gz`: 19,501 syntenic
  pairs. The DVS_A side uses the Clemson GenBank locus tags (`KPL70_000001`), not RefSeq
  `LOC...` IDs. Same assembly, so map KPL70 -> LOC by coordinate overlap
  (GCA_022201045.1 GFF vs our RefSeq GFF).
- Alternative or cross-check: DIAMOND reciprocal best hits, RefSeq proteins vs ZK8 proteins.
- Genes without a partner (lineage-specific, or lost in one parent) are counted on their own
  and flagged. "*Poncirus*-only genes" was one reason the Makefile comment gave for going
  two-parent.

### 3.2 Functional and regulatory annotation

| resource | what | IDs | use |
|---|---|---|---|
| eggNOG-mapper (run locally) | GO, KEGG, Pfam, COG for any protein set | ours | primary functional annotation for both genomes |
| DIAMOND vs TAIR10/Araport11 | Arabidopsis best hits | ours | needed to use the crown-gall marker lists (§3.4) and MapMan/Mercator4 bins |
| PlantRegMap `08-download/Citrus_sinensis/` | GO, **promoter sequences, genome-wide and in-promoter TFBS predictions, predicted TF->target regulation** | `orange1.1g...` (JGI v1 annotation) | ready-made promoter/TFBS layer; needs orange1.1 -> LOC map (protein RBH) |
| PlantTFDB `Csi_TF_list.txt.gz` | 2,255 TFs / 1,338 loci in 58 families | `orange1.1g...` | TF families for the DE lists |
| JASPAR 2024 CORE plants | curated TF motifs | — | motif scanning/enrichment (§4 E) |
| PLAZA dicots 5.0 | orthologous gene families, GO | — | site blocked automated access; citrus coverage **(unverified)** |

No citrus ATAC-seq/DAP-seq promoter map was found in this pass. PRJDB41296 includes one
*C. sinensis* DAP-seq for a single TF, which isn't useful genome-wide.

### 3.3 Public healthy-tissue RNA-seq (baselines)

No citrus crown-gall transcriptome was found; the only crown-gall host transcriptomes are in
Arabidopsis (§3.4). Candidates, best first. Run accessions for the tiered sets are in
`citrus_baselines.tsv`.

**Carrizo citrange (genotype-matched to the CC galls)**

| tier | BioProject | tissue (n) | library | notes |
|---|---|---|---|---|
| 1 | PRJNA1216034 | wt **stems** (3, 18 d old), young leaves (3), untreated mature leaves (3) | cDNA PE HiSeq 2500, ~18-26 M | CsCOI1 study (*Hortic Res* 2025, uhaf174); paper confirms Carrizo as WT and 3 reps; they used SWO v3 + Salmon. SRA organism is just "Citrus". Ignore the cscoi1 mutant, H2O-inoculated and *Xcc*-infected sets |
| 2 | PRJNA839431 | leaf controls "CT-Car" (3) | cDNA PE NovaSeq, ~26 M | stress-combination study; SRA mislabels organism as *C. sinensis*; Cleopatra mandarin samples in the same project |
| 2 | PRJNA668159 | leaf "Carrizo Control" (3) | PCR PE NovaSeq, ~22 M | FT-overexpression study |
| 2 | PRJNA1053671 | leaf WT "ZC-1..3" (3) | RANDOM PCR PE NovaSeq, ~21 M | SDE transgenic study |
| 3 | PRJNA1092419 | apical leaves (6) | cDNA PE, ~25 M | volatile-exposure study; which runs are controls **(unverified)** |
| 3 | PRJNA532644 | roots, uninfected (2: 24 h, 48 h) | cDNA PE, ~12 M | *Phytophthora* study |
| 3 | PRJNA1126089 / PRJEB20758 | leaf control (1) / tissue unstated (CC T0, 2) | — | too small to use alone |

Gap: no Carrizo bark or older stem. The only stem set is 18-day-old shoots vs 54- and 547-day
galls. The pure-parent bark/stem and callus sets below fill the tissue gap, per parental copy.

**Poncirus (one parent; tissue references and *Poncirus*-copy comparisons)**

| tier | BioProject | tissue (n) | library | notes |
|---|---|---|---|---|
| 1 | PRJDB41296 | ZK stem (3), ZK thorn (3) | polyA (poly-T beads), PE NovaSeq, ~23 M reads | **Same accession family as the ZK8 reference.** Thorn-identity study (SHI/STY TF; also PRJCA046320), released 2026-06 |
| 2 | PRJNA482734 | trifoliate stem (6: "N" x3, "FN" x3), leaf | "RANDOM" selection, PE, ~41 M | N/FN meaning unexplained; also has clementine stem/leaf |
| 2 | PRJNA1222399 | 'Flying Dragon' (4) and 'Rich 16-6' (6) rootstock | polyA PE HiSeq 4000, ~28 M | tissue not in the metadata **(unverified; likely rootstock stem/bark)**; released 2026-01 |
| 3 | PRJNA487128 | ZK root, leaf, seed, fruit, ovules (3 each) | cDNA PE, ~35 M | tissue atlas, useful for "is this gene just root-like?" |
| 3 | PRJNA816480 | trifoliate root, non-infected vs CDVd viroid | RANDOM, PE | controls only; root |

**Sweet orange (genotype-matched to Hamlin; also the other parent of Carrizo)**

| tier | BioProject | tissue (n) | library | notes |
|---|---|---|---|---|
| 1 | PRJNA599503 | **healthy bark (3)**, root (3), young leaf (3); HLB-infected counterparts | cDNA PE HiSeq 2500, ~28 M | cultivar and polyA vs total **(unverified)**; used in HLB meta-analyses (Peng et al. 2021) |
| 1 | PRJNA778304 | Valencia embryogenic **callus**, empty-vector line EV-L44 at 0/15/30 d (3 each) | cDNA PE NovaSeq, ~27 M | proliferating-tissue comparator; ignore the MIM transgenic lines |
| 2 | PRJNA816480 | sweet orange **stem** on trifoliate rootstock, non-infected vs CDVd | RANDOM, PE | control count per tissue **(unverified)**; paper PMC9228058 (*Microorganisms* 2022) |
| 2 | PRJNA1222399 | sweet orange scion on 'Flying Dragon' (4) / 'Rich 16-6' (6) | polyA PE HiSeq 4000, ~28 M | tissue **(unverified)** |
| 3 | PRJNA386941 | sweet orange stem on *P. trifoliata* (2) or *C. junos* (2) rootstock | SE BGISEQ-500, ~24 M x 50 bp | low replication, short SE |
| — | PRJNA703546 | Washington navel **leaf** | — | leaf only, skip |

Why several: any single study confounds "gall" with its own lab, cultivar and growth
conditions. A gene that is up in gall vs **every** healthy tissue from **two or more studies**
is robust to that. A gene that is up vs stem/bark but at callus level is proliferation, not a
bacterium-specific effect.

### 3.4 Crown-gall expectations (Arabidopsis)

Deeken et al. 2006 *Plant Cell* 18:3617 (PMC1785400): C58 tumors vs uninfected inflorescence
stalk, 35 dpi, 4 reps, ATH1 array. Also Lee et al. 2009 *Plant Cell* (early infection), and the
Gohlke & Deeken 2014 review (PMC4006022). Arabidopsis mature crown-gall arrays are in
ArrayExpress E-GEOD-13927 (check whether that is the Deeken set before using it).

Up in tumors (Deeken 2006): GH3 auxin-responsive At2g23170 (56x), SuSy3 At3g43190,
SuSy5 At5g20830, STP4 At3g19930, cell-wall invertase At3g13790, PDC1 At4g33070,
ADH At1g77120, PIP2;5 At3g54820 (12x), PIP1;3 At1g01620, SAD At1g43800, FAD3 At2g29980,
At2g31360, LTP2 At2g38530 (10x), nitrate transporters At3g45060 / At5g60780,
AA transporters At1g47670 / At1g25530, oligopeptide transporters At4g21680 (17x) / At1g59740,
ASA1 At5g05730, GAD1 At5g17330, PGP1 At2g36910, PGP4 At2g47000, expansins At1g69530 /
At1g26770 / At2g28950, PUMP1 At3g54110.

Down: BCAT At3g19710 (49x), vacuolar invertases At1g12240 / At1g62660, NRT1 At2g26690 /
At3g21670, AMT2 At4g13510, AMT1;1 At2g38290, TIP2;2 At4g17340, CUT1 At1g68530, WAX2 At5g57800,
and photosynthesis broadly (light reactions, 24/27 Calvin-cycle genes).

Processes (review): ABA up with drought/ABA-responsive genes, suberization up and cutin down,
SA/ET markers up, hypoxia/fermentation, sink metabolism, MYB/bHLH/bZIP/AP2 TF changes,
genome hypermethylation (T-DNA oncogene promoters stay unmethylated).

## 4. Process

**A. Two-parent Carrizo reference.**
Concatenate DVS_A1.0 (RefSeq FASTA+GFF) and CGD ZK8 (FASTA+GFF) into a `carrizo` host, index
it with each strain as now, and realign the 4 CC galls. Quantify per parental copy. Options:
- Salmon on the combined transcriptome: its EM step shares reads that fit both copies. This is
  what the Carrizo CsCOI1 study used, against sweet orange only.
- HISAT2 + featureCounts, counting only confidently placed reads (MAPQ filter) per copy and
  reporting the ambiguous share.
Gene-level count = sweet-orange copy + *Poncirus* copy (via the §3.1 pairs). Per-copy counts
are kept for D5/E1. Report each sample's genome-wide *Poncirus*-copy share (an F1 should be
~50%); that's the §1.1 genotype check. Hamlin keeps DVS_A1.0 (optionally + DVS_B1.0).

**B. Process baselines through the same pipeline.**
fastp -> HISAT2 against the plant-only index of the matching host -> featureCounts `-s 0`,
same GFF, same `-t gene -g ID` (or the same Salmon index). Pull fastq from ENA (`fastq_ftp`
field) so sra-tools isn't needed. Carrizo baselines -> `carrizo`; sweet orange -> DVS_A1.0;
*Poncirus* -> ZK8. SE datasets count reads, not pairs; note it in the matrix. Keep one matrix per
host reference: gall + baseline columns, with a `batch`/`study` column in the sample sheet.

**C. Annotation layer (once per genome).**
`make annotate` (Makefile step 7; run it as `sbatch scripts/annotate.slurm`) does the
functional part for RefSeq-annotated hosts. It keeps the longest protein per gene, runs
eggNOG-mapper (GO, KEGG, EC, Pfam, description) plus a DIAMOND best hit each way against the
Arabidopsis proteome (Ensembl Plants release 63), and writes
`references/plant_host/<host>/<host>.genes.tsv`: one row per counted gene, keyed on the same
gene IDs as `04_matrix/plant_<host>.tsv`, with a reciprocal-best-hit flag on the Arabidopsis
match. ZK8 needs a protein FASTA first (from its GFF). Still to add: the Cs<->Pt pair table
(MCscan anchors + RBH), which doubles as the Carrizo parental-copy map, and the orange1.1->LOC
map for the PlantRegMap/PlantTFDB files (TF family, TFBS).

**D. Gene-level analyses.**
1. *Gall vs baselines*, per host genotype: Carrizo galls vs Carrizo baselines, Hamlin vs sweet
   orange. DESeq2, using baseline replicates for dispersion and `study` as a covariate where it
   isn't aliased. Run it against each baseline separately, then intersect. Call a gene only if
   |log2FC| > 2 and it's consistent across >=2 studies. Report within-sample percentile ranks
   next to fold changes; ranks are more robust to cross-study batch.
   Tissue logic: up vs bark/stem/root/leaf **and** vs callus = gall-specific candidate; up vs
   stem but callus-like = proliferation program. Carrizo has no bark or callus baseline, so use
   the sweet-orange bark/callus and ZK stem sets per parental copy for that check.
2. *Signature check*: score the Deeken up/down sets (via Arabidopsis orthologs) in each gall
   vs baseline. If citrus galls don't show the canonical program (photosynthesis down,
   fermentation/SuSy/STP up, CUT1/WAX2 down), suspect the baseline before the biology.
3. *Internal contrasts* (host genotype held constant; valid even on the current sweet-orange-only
   counts):
   - wt CC547 vs G-19 + G-30 CC547: oncogene/T-DNA-dependent host genes (1 vs 2; fold-change
     ranking, or DESeq2 with dispersion borrowed from the Carrizo baselines). Pending the
     §1.1 genotype check.
   - CC54 vs CC547 wt: gall age
   - HC83 vs CC wt: host and strain confounded (29 vs 1416) — descriptive only
4. *Cross-host agreement*: genes moving the same way in the Hamlin gall and the wt Carrizo
   galls, each relative to its own genotype's baselines. Compare Hamlin both with Carrizo's
   gene-level totals and with its sweet-orange copies specifically. That's the most defensible
   "bacterium-driven" set, given the design.
5. *Parental-copy expression (Carrizo only)*: for each paired gene, the *Poncirus*-copy share in
   gall vs in healthy Carrizo tissue. A gene whose share shifts in the gall responds in cis:
   one parent's copy is induced or repressed more than the other. Model it as a count GLM with
   copy x condition interaction (DESeq2-style, sample as a factor) or beta-binomial per gene.
   Also compare wt vs G-19/G-30 for oncogene-dependent shifts. With n=1 galls this is a ranking
   exercise, but the healthy baselines supply the "normal" ratio and its variance.

**E. Promoters.**
1. *Parental-promoter comparison* (new; uses D5): for genes with a gall-induced copy shift,
   align the sweet-orange and *Poncirus* promoters (1-2 kb upstream of TSS), scan both with
   JASPAR plant motifs (FIMO), and list motifs gained or lost in the responsive copy. Same
   nucleus, same bacterium, same hormones, so a sequence difference that tracks the expression
   difference is the strongest promoter evidence this design can give.
2. *Host promoters of candidate genes*: 1-2 kb upstream of the annotated TSS (both annotations
   have UTRs). Foreground = the gene sets from D (up and down separately); background =
   expressed genes matched on expression level.
   - Known-motif enrichment (MEME-suite SEA/AME or HOMER, JASPAR 2024 plants). Expect
     auxin (AuxRE TGTCTC / ARF), type-B ARR cytokinin (AGATHY), ABRE/ABA, hypoxia (HRPE),
     W-box/SA if the Arabidopsis program holds.
   - PlantRegMap's precomputed *C. sinensis* TFBS and regulation tables: shortcut to "which
     TFs could drive this set".
   - Phylogenetic footprinting: for genes that respond in both parental copies (and in Hamlin),
     motif instances conserved between the Cs and Pt promoters are much stronger evidence than
     enrichment alone.
   - Set expectations: with n=1, promoter results are hypotheses for follow-up. Haryono et al.
     (HARYONO2019_PIPELINE.md) found ab initio motifs unstable even with replicates.
3. *The bacterium's own plant-active promoters*: extract the intergenic regions upstream of
   each expressed T-DNA gene in strains 1416 and 29 (§1.2 table), annotate TATA and ocs/as-1
   elements (TGACG, bound by host TGA/bZIP factors), and compare between strains. These drive
   the only transcripts the bacterium puts in the plant genome directly. They're also the
   promoter parts the GAANTRY/"Symbionts" work would reuse.

**F. Candidate shortlist for the follow-up experiment (decided 2026-09-14).**
Don't trust one statistic: run several simple ranking methods, take each method's **top 20**,
and compare the lists **by hand** alongside the annotation table
(`references/plant_host/<host>/<host>.genes.tsv`, from `make annotate`). Genes that come up
under several methods are the strongest candidates. Make separate lists for Hamlin, the wt
Carrizo galls, and the G-19/G-30 Carrizo galls; the last is the closest to an engineered-strain
gall.
1. **Gall expression**: within-sample percentile (or TPM), consistent across the galls. Measures
   promoter strength, and needs only the existing gall counts.
2. **Fold change vs baseline**: DESeq2 with shrunken log2FC (`lfcShrink`) and a minimum-
   expression cut-off, so near-zero genes can't top the list.
3. **Percentile rank shift**: gall percentile minus the median healthy-tissue percentile. Robust
   to cross-study depth and library differences.
4. **Rank products** across the galls (Breitling et al. 2004, *FEBS Lett* 573:83): genes ranked
   near the top of the per-gall change lists in every gall. Rewards consistency.
5. **RankComp / relative expression orderings** (Wang et al. 2015, *Bioinformatics* 31:62): gene
   pairs whose order is stable in healthy tissue and flips in the gall. Designed for one sample
   against a reference set from other studies.
6. **Tissue specificity, tau** (Yanai et al. 2005, *Bioinformatics* 21:650), with gall as one
   tissue alongside stem, bark, leaf, root and callus: gall-specific genes.
7. Listed separately: the T-DNA promoters (§1.2), as bacterial benchmarks.

Before ranking, drop genes not clearly expressed in the galls, organellar and rRNA genes, and
multi-copy families whose promoter can't be cloned unambiguously. Methods 2-6 need healthy-
tissue counts: the tier-1 runs in `citrus_baselines.tsv` are enough (Carrizo stem and leaves,
PRJNA1216034; sweet-orange bark, root, leaf and callus, PRJNA599503 / PRJNA778304). Output: one
merged table (gene, annotation, Arabidopsis ortholog, rank under each method, expression per
gall) plus the ~1.5 kb upstream sequence of each pick. The picks then go to a reporter assay
(GUS/luciferase in galls), which is the actual measurement of promoter activity.

## 5. Proposed next steps

0. **Shortlist (§4 F), first priority.** Once `make annotate` has written the sweet-orange gene
   table: method 1 from the existing gall counts, then align/count the tier-1 baselines as a
   Slurm job for methods 2-6, then the merged top-20 table for manual review.
1. Makefile: add a `carrizo` host (DVS_A1.0 FASTA+GFF + CGD ZK8 FASTA+GFF from the §3.1 URLs),
   point `HOST_*CC*` at it, add it to `ANNOT_HOSTS`, then
   `make HOSTSEL=carrizo fractions counts matrices`. Keep the current sweet-orange-only CC
   counts for comparison until the two-parent numbers are in.
2. `make baselines` target reading `citrus_baselines.tsv` (ENA fastq download -> fastp ->
   align -> count into `04_matrix/baseline_<host>.tsv`). Tier 1: 33 runs, ~0.81 billion read
   pairs (9 Carrizo stem/leaf, 6 *Poncirus* stem/thorn, 9 sweet orange bark/root/leaf,
   9 sweet orange callus). Tier 2: 9 Carrizo leaf controls, ~0.21 billion. All run accessions
   were matched to their sample labels on NCBI/ENA on 2026-09-14.
3. Annotation layer (§4 C): sweet orange via `sbatch scripts/annotate.slurm` (eggNOG-mapper in
   a venv from `scripts/emapper-env.sh`, DIAMOND from the cluster module). Next: ZK8 proteins,
   the Cs<->Pt pair table, and Salmon if it's chosen for §4 A.
4. DE notebook `notebooks/03_citrus_host_de.ipynb` for §4 D. R/DESeq2 isn't installed yet
   (PIPELINE.md step 6).

## 6. Asks for collaborators (Shatters / Thomson labs)

- **Any uninoculated tissue from the same plants/greenhouse**: wounded-mock or adjacent healthy
  bark from the CC547 and HC83 plants, even 2-3 samples. This would be worth more than all the
  public baselines combined, and would turn §4 D1 from exploratory into a real contrast.
- G-19 / G-30 construct maps: what replaced or modified the 1416 T-DNA.
- The CC plants: confirm Carrizo, seed source, whether they are nucellar seedlings (not grafted),
  and whether the G-19/G-30 plants came from the same stock as the wt ones (§1.1 QC flag).
- RNA extraction and library kit, to confirm the polyA inference and the tissue boundaries of the
  gall sample (gall only, or gall + adjacent stem?).
