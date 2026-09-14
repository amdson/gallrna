# Citrus host response: which host genes/promoters did the gall bacterium change?

> Scoping doc, 2026-09-14. Covers the 5 citrus-host samples only. Nothing here is automated yet;
> §5 lists proposed Makefile additions. Public accessions were checked against NCBI/ENA/CGD
> on this date; anything not verified is marked **(unverified)**.

## 1. What we have

| sample | host (samples.tsv) | strain | gall age (d) | assigned pairs | genes >=10 counts |
|---|---|---|---|---|---|
| 29wtHC83 | *Citrus sinensis* 'Hamlin' | 29 wt | 83 | 26.0 M | 18,239 |
| 1416wtCC547 | *Poncirus trifoliata* | 1416 wt | 547 | 21.3 M | 17,023 |
| 1416wtCC54 | *Poncirus trifoliata* | 1416 wt | 54 | 20.8 M | 17,132 |
| 1416G-19CC547 | *Poncirus trifoliata* | 1416 G-19 | 547 | 22.0 M | 17,239 |
| 1416G-30CC547 | *Poncirus trifoliata* | 1416 G-30 | 547 | 20.7 M | 16,753 |

Counts: `04_matrix/plant_citrus_sinensis.tsv` (all 5 currently mapped to sweet orange, see §3).

Library facts established 2026-09-14 (needed to match public data):
- **polyA mRNA**: intronic reads are ~1% of located reads (`logs/intron_fraction.tsv`).
- **Unstranded**: featureCounts `-s 1` vs `-s 2` on a 2.6 M-read slice of 29wtHC83 gave
  250,875 vs 249,061 assigned. `STRAND=0` is correct; keep it for everything below.
- 1416G-30CC547 has 22.0 M multimapping pairs vs 4.6-9.2 M in the other four — check what
  (rRNA? repeat? construct?) before trusting it in contrasts.

**There is no uninoculated citrus tissue.** Every sample is a gall and every condition is n=1.
"What did the bacterium change" therefore has to be triangulated (§4) rather than read off a
gall-vs-mock contrast.

### The bacterium's direct footprint: T-DNA expression in the galls
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

Without mock tissue, three independent lines of evidence each carry a known bias. A host gene
is a good "influenced by the bacterium" candidate when they agree:

1. **Gall vs public healthy tissue** (§3.3). Captures the gall program but is confounded with
   lab, cultivar, tissue and age. Mitigated by using several baselines and several tissues.
2. **Internal contrasts** that hold host accession constant. The best one is wt CC547 vs
   G-19/G-30 CC547: same host accession and same age, with vs without wt T-DNA expression.
3. **Prior crown-gall signatures** (Arabidopsis, via orthologs), which say which direction
   known gall genes should move.

Promoter analysis (§4 step E) comes after the gene sets exist.

## 3. Data available

### 3.1 Host references and annotation

**Sweet orange — keep the current reference.**
- NCBI RefSeq GCF_022201045.2 (DVS_A1.0, Valencia, haplotype A), Annotation Release 103.
  Already in `references/plant_host/citrus_sinensis/`. 23,566 protein-coding genes, but ~7,800
  genes are "uncharacterized"/"hypothetical" and there is no GO/KEGG in the GFF, so functional
  annotation has to be added (§3.2).
- Newer assemblies (T2T SWO GCA_046470765.1/775.1, Pera Rio, Washington Navel, Newhall, ...)
  have no NCBI annotation, so there's no reason to switch.

**Poncirus — an annotated reference exists; the Liftoff workaround in the Makefile is no longer needed.**
- *P. trifoliata* ZK8 v1.0 (HZAU; Peng et al. 2020, *Plant J* doi:10.1111/tpj.14993).
  Same assembly as the unannotated NCBI copy we already downloaded (GCA_018350135.1), but the
  Citrus Genome Database (CGD) serves the gene models too:
  - `https://www.citrusgenomedb.org/jb2/data/Ptri_ZK8_v1.sorted.gff.gz` (7.5 MB; 25,680 genes,
    39,675 mRNAs, with UTRs; IDs like `Pt1g002240`)
  - `https://www.citrusgenomedb.org/jb2/data/Ptri_ZK8_v1.fasta.gz` (86 MB)
  - CGD lists annotation BUSCO 93.4%, assembly BUSCO 98.6%.
  - Contig names differ: CGD uses `chr1_ZK8` ... `chrUn_ZK8`, NCBI uses `CM031384.1` ... plus
    `VKKW01*` scaffolds. **Use the CGD FASTA with the CGD GFF** rather than renaming.
  - These are the genome-browser data files, not a formal download page (CGD's `/download` page
    and the Analysis page did not link them). They are publicly readable. Cite Peng et al. 2020.
- JGI *P. trifoliata* v1.3.1 (Phytozome/CGD `bio_data/1312872`): the page returned 403 and the
  access terms are **(unverified)**. Not needed given ZK8.
- The citrus samples' host code is "CC / Citrus (trifoliate)". If the actual plants are
  a named cultivar or a citrange, tell us: ZK8 is a specific accession.

**Sweet orange <-> Poncirus gene map**
- CGD MCscan anchors `jb2/data/mcscan/Csin_DVS_A_v1.Ptri_ZK8_v1.anchors.gz`: 19,501 syntenic
  pairs. The DVS_A side uses the Clemson GenBank locus tags (`KPL70_000001`), not RefSeq
  `LOC...` IDs. Same assembly, so map KPL70 -> LOC by coordinate overlap
  (GCA_022201045.1 GFF vs our RefSeq GFF).
- Alternative or cross-check: DIAMOND reciprocal best hits, RefSeq proteins vs ZK8 proteins.

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
Arabidopsis (§3.4). Candidates, best first. Run accessions for tier 1 are in
`citrus_baselines.tsv`.

**Poncirus**

| tier | BioProject | tissue (n) | library | notes |
|---|---|---|---|---|
| 1 | PRJDB41296 | ZK stem (3), ZK thorn (3) | polyA (poly-T beads), PE NovaSeq, ~23 M reads | **Same accession family as the ZK8 reference.** Thorn-identity study (SHI/STY TF; also PRJCA046320), released 2026-06 |
| 2 | PRJNA482734 | trifoliate stem (6: "N" x3, "FN" x3), leaf | "RANDOM" selection, PE, ~41 M | N/FN meaning unexplained; also has clementine stem/leaf |
| 2 | PRJNA1222399 | 'Flying Dragon' (4) and 'Rich 16-6' (6) rootstock | polyA PE HiSeq 4000, ~28 M | tissue not in the metadata **(unverified; likely rootstock stem/bark)**; released 2026-01 |
| 3 | PRJNA487128 | ZK root, leaf, seed, fruit, ovules (3 each) | cDNA PE, ~35 M | tissue atlas, useful for "is this gene just root-like?" |
| 3 | PRJNA816480 | trifoliate root, non-infected vs CDVd viroid | RANDOM, PE | controls only; root |

**Sweet orange**

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

**A. Fix the Poncirus reference.**
Switch the 4 CC samples from `citrus_sinensis` to a CGD-ZK8-based `poncirus_trifoliata`,
realign, recount. Compare assignment vs the current cross-species numbers (Hamlin keeps
DVS_A1.0).

**B. Process baselines through the same pipeline.**
fastp -> HISAT2 against the plant-only index of the matching host -> featureCounts `-s 0`,
same GFF, same `-t gene -g ID`. Pull fastq from ENA (`fastq_ftp` field) so sra-tools isn't
needed. SE datasets count reads, not pairs; note it in the matrix. Keep one matrix per host
reference: gall + baseline columns, with a `batch`/`study` column in the sample sheet.

**C. Annotation layer (once per genome).**
eggNOG-mapper and DIAMOND->TAIR10 on RefSeq and ZK8 proteins. Cs<->Pt ortholog table (MCscan
anchors + RBH). orange1.1->LOC map to bring in the PlantRegMap/PlantTFDB files. Output: one
gene table per genome (ID, description, GO, KEGG, Arabidopsis hit, TF family, ortholog in
the other species).

**D. Gene-level analyses.**
1. *Gall vs baselines*: DESeq2 per species, using baseline replicates for dispersion and
   `study` as a covariate where it isn't aliased. Run it against each baseline separately, then
   intersect. Call a gene only if |log2FC| > 2 and it's consistent across >=2 studies. Report
   within-sample percentile ranks next to fold changes; ranks are more robust to cross-study
   batch.
   Tissue logic: up vs bark/stem/root/leaf **and** vs callus = gall-specific candidate; up vs
   stem but callus-like = proliferation program.
2. *Signature check*: score the Deeken up/down sets (via Arabidopsis orthologs) in each gall
   vs baseline. If citrus galls don't show the canonical program (photosynthesis down,
   fermentation/SuSy/STP up, CUT1/WAX2 down), suspect the baseline before the biology.
3. *Internal contrasts* (host accession held constant):
   - wt CC547 vs G-19 + G-30 CC547: oncogene/T-DNA-dependent host genes (1 vs 2; fold-change
     ranking, or DESeq2 with dispersion borrowed from the Poncirus baselines)
   - CC54 vs CC547 wt: gall age
   - HC83 vs CC wt: host species, but confounded with strain (29 vs 1416) — descriptive only
4. *Cross-species agreement*: genes moving the same way in the Hamlin gall and the wt Poncirus
   galls (via orthologs), each relative to its own species' baselines. That's the most
   defensible "bacterium-driven" set, given the design.

**E. Promoters.**
1. *Host promoters of candidate genes*: 1-2 kb upstream of the annotated TSS (both annotations
   have UTRs). Foreground = the gene sets from D (up and down separately); background =
   expressed genes matched on expression level.
   - Known-motif enrichment (MEME-suite SEA/AME or HOMER, JASPAR 2024 plants). Expect
     auxin (AuxRE TGTCTC / ARF), type-B ARR cytokinin (AGATHY), ABRE/ABA, hypoxia (HRPE),
     W-box/SA if the Arabidopsis program holds.
   - PlantRegMap's precomputed *C. sinensis* TFBS and regulation tables: shortcut to "which
     TFs could drive this set".
   - Phylogenetic footprinting: align Cs and Pt promoters of orthologs that respond in both
     species; conserved motif instances are much stronger evidence than enrichment alone.
   - Set expectations: with n=1, promoter results are hypotheses for follow-up. Haryono et al.
     (HARYONO2019_PIPELINE.md) found ab initio motifs unstable even with replicates.
2. *The bacterium's own plant-active promoters*: extract the intergenic regions upstream of
   each expressed T-DNA gene in strains 1416 and 29 (§1 table), annotate TATA and ocs/as-1
   elements (TGACG, bound by host TGA/bZIP factors), and compare between strains. These drive
   the only transcripts the bacterium puts in the plant genome directly. They're also the
   promoter parts the GAANTRY/"Symbionts" work would reuse.

## 5. Proposed next steps

1. Makefile: add CGD ZK8 as the `poncirus_trifoliata` source (FASTA+GFF from the CGD URLs in
   §3.1), move `HOST_*CC*` back to `poncirus_trifoliata`, add it to `ANNOT_HOSTS`, then
   `make HOSTSEL=poncirus_trifoliata fractions counts matrices`.
2. `make baselines` target reading `citrus_baselines.tsv` (ENA fastq download -> fastp ->
   align -> count into `04_matrix/baseline_<host>.tsv`). Start with the tier-1 runs:
   24 runs, ~0.62 billion read pairs (6 Poncirus stem/thorn, 9 sweet orange bark/root/leaf,
   9 sweet orange callus). All 24 run accessions were matched to their sample labels on
   NCBI/ENA on 2026-09-14.
3. Annotation layer (§4 C): needs eggNOG-mapper + DIAMOND in the conda env or as modules.
4. DE notebook `notebooks/03_citrus_host_de.ipynb` for §4 D. R/DESeq2 isn't installed yet
   (PIPELINE.md step 6).

## 6. Asks for collaborators (Shatters / Thomson labs)

- **Any uninoculated tissue from the same plants/greenhouse**: wounded-mock or adjacent healthy
  bark from the CC547 and HC83 plants, even 2-3 samples. This would be worth more than all the
  public baselines combined, and would turn §4 D1 from exploratory into a real contrast.
- G-19 / G-30 construct maps: what replaced or modified the 1416 T-DNA.
- Exact Poncirus germplasm for the CC plants (seedling? cultivar? rootstock-grafted?).
- RNA extraction and library kit, to confirm the polyA inference and the tissue boundaries of the
  gall sample (gall only, or gall + adjacent stem?).
