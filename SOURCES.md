# Sources

> Every external source, dataset and tool the current work depends on, and why it is used.
> Compiled 2026-09-14. DOIs were checked against Crossref / Europe PMC and genome accessions
> against the NCBI Datasets API on that date; anything not checked is marked **(unverified)**.
> When a Makefile target, reference or tool changes, amend its entry here.
> Some files referenced below (SLURM.md, notebooks/, several scripts/) were untracked in the
> working checkout when this was written.

1. [Method: dual RNA-seq by competitive mapping](#1-method-dual-rna-seq-by-competitive-mapping-to-a-concatenated-reference)
2. [Agrobacterium in galls](#2-agrobacterium-in-galls)
3. [Sequencing data and sample metadata](#3-sequencing-data-and-sample-metadata)
4. [Bacterial references](#4-bacterial-references)
5. [Plant host references](#5-plant-host-references)
6. [Pipeline tools (Makefile steps 0-5)](#6-pipeline-tools-makefile-steps-0-5)
7. [Checks and helper scripts](#7-checks-and-helper-scripts)
8. [Notebooks and figures](#8-notebooks-and-figures)
9. [Compute and environments](#9-compute-and-environments)
10. [Installed or planned, not yet used](#10-installed-or-planned-not-yet-used)
11. [Internal and collaborator documents](#11-internal-and-collaborator-documents)
12. [Discrepancies found while compiling](#12-discrepancies-found-while-compiling)

## 1. Method: dual RNA-seq by competitive mapping to a concatenated reference

What we do (Makefile §3c-5): each sample's strain genome (contigs prefixed `agro_`) is
concatenated with its host genome and indexed once with `hisat2-build`. Reads are aligned to
that single index, so each read goes to whichever organism it fits best. The organisms are split
afterwards by the `agro_` prefix (`logs/mapping_fractions.tsv`, `scripts/merge_counts.py`) and
counted against the concatenated GFF3.

**The same design, published**

| study | system | aligner | relevance |
|---|---|---|---|
| Wang Q, Shakoor N, Boyher A, Veley KM, Berry JC, Mockler TC, Bart RS (2021) Escalation in the host-pathogen arms race: a host resistance response corresponds to a heightened bacterial virulence response. *PLoS Pathog* 17:e1009175. [doi:10.1371/journal.ppat.1009175](https://doi.org/10.1371/journal.ppat.1009175) | sorghum + *Xanthomonas vasicola* pv. *holcicola* | HISAT2 2.0.6 | **Primary citation for our design.** Plant + bacterium, HISAT2, one index: "aligned the RNA-seq reads against a concatenated genome comprised of the sorghum nuclear genome …, the sorghum chloroplast genome …, the sorghum mitochondrial genome … and the *Xvh* genome to allow for the best alignment of each read." |
| Chakraborty S, Sharma R, Bhat A, et al. (2025) Partners in root nodule symbiosis respond uniquely to heavy metal stresses in a host genotype-dependent manner. *Sci Rep* 15:33518. [doi:10.1038/s41598-025-17827-z](https://doi.org/10.1038/s41598-025-17827-z) | *Medicago truncatula* nodules + *Sinorhizobium meliloti* | STAR | Plant + bacterium; host genome and annotation concatenated with the rhizobium genome before mapping. |
| Human MP, Berger DK, Crampton BG (2020) Time-course RNAseq reveals *Exserohilum turcicum* effectors and pathogenicity determinants. *Front Microbiol* 11:360. [doi:10.3389/fmicb.2020.00360](https://doi.org/10.3389/fmicb.2020.00360) | maize + *E. turcicum* (fungus) | TopHat 2.0.13 | Host and pathogen genomes concatenated. |

**Benchmarks: concatenated vs sequential mapping**

- Espindula E, Sperb ER, Bach E, Passaglia LMP (2019) The combined analysis as the best strategy
  for Dual RNA-Seq mapping. *Genet Mol Biol* 42(4):e20190215.
  [doi:10.1590/1678-4685-gmb-2019-0215](https://doi.org/10.1590/1678-4685-gmb-2019-0215).
  Simulated maize + *Herbaspirillum seropedicae*; real soybean + *Bradyrhizobium elkanii* and
  maize + *Fusarium verticillioides*; CLC Genomics Workbench. Host-first sequential mapping put
  ~13 M bacterial reads on maize; the combined reference cut that to ~0.53 M (accuracy 0.99 vs
  0.79). They then re-mapped each sorted read set to its own genome; we count from the combined
  alignment directly, as Wang et al. do.
- Fruggiero C, Aufiero G, D'Angelo D, Pasolli E, D'Agostino N (2024) Refining dual RNA-seq
  mapping: sequential and combined approaches in host-parasitic plant dynamics.
  *Front Plant Sci* 15:1483717. [doi:10.3389/fpls.2024.1483717](https://doi.org/10.3389/fpls.2024.1483717).
  *Arabidopsis* and tomato + *Cuscuta campestris*, STAR. ~90 % mapped and ~1 % cross-mapping
  either way; combined was slightly more accurate and took half the time. Two plants are far
  harder to separate than a plant and a bacterium, so this is a conservative test for us.
- Sollero BP, de Jesus Costa RM, Togawa RC, de Castro Costa É, Miller RNG (2026) Genome
  concatenation enables accurate dual RNA-seq mapping for lignocellulolytic fungi in coculture.
  *Sci Rep* 16:16898. [doi:10.1038/s41598-026-43548-y](https://doi.org/10.1038/s41598-026-43548-y).
  Two fungi, STAR: the same accuracy (>96 % mapped), <2 % cross-mapping, ~93 % less compute.

**The same idea elsewhere**

- Callari M, Batra AS, Batra RN, et al. (2018) Computational approach to discriminate human and
  mouse sequences in patient-derived tumour xenografts. *BMC Genomics* 19:19.
  [doi:10.1186/s12864-017-4414-y](https://doi.org/10.1186/s12864-017-4414-y). Combined
  human-mouse reference, STAR; 99.8 % of human and 98.9 % of mouse RNA-seq reads assigned
  correctly (MAPQ > 0).
- nf-core/dualrnaseq (B. Mika-Gospodorz, R. Hayward; Barquist lab): a community dual RNA-seq
  pipeline whose STAR mode builds "a chimeric genome index (combining host and pathogen
  genomes)". <https://nf-co.re/dualrnaseq>. No methods paper; cite it as software.
- Background: Westermann AJ, Gorski SA, Vogel J (2012) Dual RNA-seq of pathogen and host.
  *Nat Rev Microbiol* 10:618-630. [doi:10.1038/nrmicro2852](https://doi.org/10.1038/nrmicro2852).

**Counterpoint: sequential mapping**

- Infanta Saleth Teresa Eden M, Vetrivel U (2025) Optimal dual RNA-seq mapping for accurate
  pathogen detection in complex eukaryotic hosts. *Bio-protocol* 15:e5182.
  [doi:10.21769/BioProtoc.5182](https://doi.org/10.21769/BioProtoc.5182). Human macrophages +
  *M. tuberculosis*; recommend pathogen-first, then host.
- Liao Z-X, Ni Z, Wei X-L, et al. (2019) Dual RNA-seq of *Xanthomonas oryzae* pv. *oryzicola*
  infecting rice reveals novel insights into bacterial-plant interaction. *PLoS ONE* 14:e0215039.
  [doi:10.1371/journal.pone.0215039](https://doi.org/10.1371/journal.pone.0215039). Host first
  (TopHat), unmapped reads to the pathogen.

Neither compares itself against a combined reference; the three benchmarks above do.

## 2. Agrobacterium in galls

- González-Mula A, Lang J, Grandclément C, et al. (2018) Lifestyle of the biotroph
  *Agrobacterium tumefaciens* in the ecological niche constructed on its host plant.
  *New Phytol* 219:350-362. [doi:10.1111/nph.15164](https://doi.org/10.1111/nph.15164).
  First transcriptome of *A. tumefaciens* (C58) living in plant tumors (*Arabidopsis*): a
  motile-to-sessile switch, cell-surface remodelling, and use of opines, amino acids, sugars,
  organic acids, phosphate and iron. The closest precedent for the bacterial half of our data.
  How they separated plant from bacterial reads is **(unverified)**: the full text and the
  thesis behind it (theses.hal.science tel-03505903) were not reachable.
- Gonzalez-Mula A, Lachat J, Mathias L, et al. (2019) The biotroph *Agrobacterium tumefaciens*
  thrives in tumors by exploiting a wide spectrum of plant host metabolites. *New Phytol*
  222:455-467. [doi:10.1111/nph.15598](https://doi.org/10.1111/nph.15598). Tomato tumors vs
  healthy stems: metabolomics (sucrose, glutamate; enriched opines, GABA, GHB, pyruvate), Tn-seq,
  and transcriptomics of C58 grown on GHB or sucrose. Tomato is one of our hosts (T29, T49).
- Gonzalez-Mula A, Torres M, Faure D (2019) Integrative and deconvolution omics approaches to
  uncover the *Agrobacterium tumefaciens* lifestyle in plant tumors. *Plant Signal Behav*.
  [doi:10.1080/15592324.2019.1581562](https://doi.org/10.1080/15592324.2019.1581562). Short
  addendum summarising the two papers above.
- Haryono M, Cho S-T, Fang M-J, Chen A-P, Chou S-J, Lai E-M, Kuo C-H (2019) Differentiations in
  gene content and expression response to virulence induction between two *Agrobacterium*
  strains. *Front Microbiol* 10:1554.
  [doi:10.3389/fmicb.2019.01554](https://doi.org/10.3389/fmicb.2019.01554). Template for
  cross-strain orthology and promoter analysis (PIPELINE.md steps 8-9; HARYONO2019_PIPELINE.md).
- Deeken R, Engelmann J, Efetova M, Czirjak T, et al. (2006) An integrated view of gene
  expression and solute profiles of *Arabidopsis* tumors: a genome-wide approach. *Plant Cell*
  18:3617-3634. [doi:10.1105/tpc.106.044743](https://doi.org/10.1105/tpc.106.044743). Host
  crown-gall expression signature (CITRUS_HOST_PLAN.md §3.4).
- Gohlke J, Deeken R (2014) Plant responses to *Agrobacterium tumefaciens* and crown gall
  development. *Front Plant Sci* 5:155.
  [doi:10.3389/fpls.2014.00155](https://doi.org/10.3389/fpls.2014.00155). Review.
- Lee C-W, et al. (2009) *Agrobacterium tumefaciens* promotes tumor induction by modulating
  pathogen defense in *Arabidopsis*. *Plant Cell*. PMC2768927. Early infection (cited in the
  plan; volume and pages **(unverified)**).

## 3. Sequencing data and sample metadata

- **Reads:** `30-1348328766/00_fastq/`: Azenta/GENEWIZ project 30-1348328766, Illumina NovaSeq X,
  2 × 150 bp paired-end, 14 gall samples, R1 + R2 with md5 files (checked by `make verify` →
  `logs/md5.ok`). Not public.
- **Samples:** `samples.tsv`: strain, genotype, host, gall age and mass, source lab (Shatters,
  Thomson), inoculation and collection dates. Hosts were confirmed with the collaborators on
  2026-09-14 (CC = Carrizo citrange, HC = Hamlin sweet orange).
- **Library type, read back from the data** (the kit is still to be confirmed with the lab):
  - polyA mRNA: intronic reads are ~1 % of located reads (`scripts/intron_fraction.sh` →
    `logs/intron_fraction.tsv`). A ribo-depleted library would carry far more intronic signal:
    Zhao W, He X, Hoadley KA, Parker JS, Hayes DN, Perou CM (2014) Comparison of RNA-Seq by
    poly(A) capture, ribosomal RNA depletion, and DNA microarray for expression profiling.
    *BMC Genomics* 15:419. [doi:10.1186/1471-2164-15-419](https://doi.org/10.1186/1471-2164-15-419).
  - Unstranded: featureCounts `-s 1` vs `-s 2` on a 2.6 M-read slice of 29wtHC83 assigned
    250,875 vs 249,061 pairs, hence `STRAND=0` (CITRUS_HOST_PLAN.md §1).

## 4. Bacterial references

Built by `make strain-refs` with `scripts/convert_refs.py` (Biopython); details in
`references/README.md`.

**strain_1416** (1416\* samples): `Strain 1416/1416 * FINAL.gb`, collaborator files, annotated
in SnapGene. 5 replicons, 5.86 Mb.
- Alabed D, Huo N, Gu Y, Thomson JG (2023) Complete genome of *Agrobacterium fabrum* strain
  1D1416. *Microbiol Resour Announc* 12:e00264-23.
  [doi:10.1128/mra.00264-23](https://doi.org/10.1128/mra.00264-23). The strain was isolated from
  *Euonymus japonicus*, one of our hosts, and forms galls on citrus.
- Alabed D, Tibebu R, Ariyaratne M, Shao M, Milner MJ, Thomson JG (2024) Novel *Agrobacterium
  fabrum* str. 1D1416 for citrus transformation. *Microorganisms* 12:1999.
  [doi:10.3390/microorganisms12101999](https://doi.org/10.3390/microorganisms12101999).

**strain_29** (29\* samples): `Strain 29/.../Xtra RAST download/357.66*.gbk`, collaborator files.
The "Original Download" .gb files are sequence-only, so the RAST-annotated versions are used.
5 replicons, 5.84 Mb.
- No publication or strain designation in hand **(ask the collaborators)**.
- Annotation: Aziz RK, Bartels D, Best AA, et al. (2008) The RAST Server: rapid annotations using
  subsystems technology. *BMC Genomics* 9:75.
  [doi:10.1186/1471-2164-9-75](https://doi.org/10.1186/1471-2164-9-75); Brettin T, Davis JJ,
  Disz T, et al. (2015) RASTtk. *Sci Rep* 5:8365.
  [doi:10.1038/srep08365](https://doi.org/10.1038/srep08365).
- `docs/` also holds Huo N, Gu Y, McCue KF, Alabed D, Thomson JG (2019) Complete genome sequence
  of *Agrobacterium fabrum* strain 1D159. *Microbiol Resour Announc* 8:e00207-19.
  [doi:10.1128/MRA.00207-19](https://doi.org/10.1128/MRA.00207-19) (1D159 = ATCC 27912). Its
  relation to strain 29 is **(unverified)**.
- The 1416 (SnapGene) and 29 (RAST) annotations come from different pipelines. Harmonise them
  by orthology before comparing genes across strains (PIPELINE.md step 8).

**C58** (comparison and orthology only; no samples map to it): GenBank AE007869 (circular
chromosome, 2,841,580 bp), AE007870 (linear chromosome, 2,075,577 bp), AE007871 (pTi,
214,233 bp), AE007872 (pAt, 542,868 bp); local SnapGene copies in `C58_ALIGNED/`.
- Wood DW, Setubal JC, Kaul R, et al. (2001) The genome of the natural genetic engineer
  *Agrobacterium tumefaciens* C58. *Science* 294:2317-2323.
  [doi:10.1126/science.1066804](https://doi.org/10.1126/science.1066804)
- Goodner B, Hinkle G, Gattung S, et al. (2001) Genome sequence of the plant pathogen and
  biotechnology agent *Agrobacterium tumefaciens* C58. *Science* 294:2323-2328.
  [doi:10.1126/science.1066803](https://doi.org/10.1126/science.1066803)

Not used by the pipeline: the collaborator's Snippy variant calls and Bakta annotations of
strain 29 (still on Atlas; SETUP.md). Bakta: Schwengers O, Jelonek L, Dieckmann MA, et al.
(2021) *Microb Genom* 7:000685. [doi:10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685).
Snippy: <https://github.com/tseemann/snippy>.

## 5. Plant host references

All are downloaded by `make host-genomes` from the NCBI Datasets v2 API
(`https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<ACC>/download`). Cite:
O'Leary NA, Cox E, Holmes JB, et al. (2024) Exploring and retrieving sequence and metadata for
species across the tree of life with NCBI Datasets. *Sci Data* 11:732.
[doi:10.1038/s41597-024-03571-y](https://doi.org/10.1038/s41597-024-03571-y).
Assembly metadata below comes from the Datasets report (2026-09-14). The host panel and the
fallback choices are in `references/plant_host/AVAILABILITY.md`.

| Makefile key | samples | accession | assembly (cultivar), submitter | level | NCBI annotation |
|---|---|---|---|---|---|
| `euonymus_japonicus_proxy` | 1416wtEu635, 29wtEu635, 29wtGeu182, 1416wtGeu182 | [GCA_963580455.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_963580455.1/) | drEuoEuro1.1, *Euonymus europaeus*; Wellcome Sanger Institute | chromosome | none |
| `citrus_sinensis` | 29wtHC83; also the four CC galls (below) | [GCF_022201045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_022201045.2/) | DVS_A1.0 (Valencia, haplotype A); Clemson University | complete | RefSeq Annotation Release 103, 23,556 protein-coding genes |
| `brassica_juncea` | 1416wtM26, 29wtM26 | [GCA_018703725.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_018703725.1/) | ASM1870372v1 ('Sichuan Huangzi'); Hunan Agricultural University | chromosome | none |
| `carica_papaya` | 29wtP46 | [GCF_000150535.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000150535.2/) | Papaya1.0 ('SunUp'); Papaya Genome Sequencing Consortium | scaffold | RefSeq Annotation Release 100, 18,126 protein-coding genes |
| `poncirus_trifoliata` | none at present | [GCA_018350135.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_018350135.1/) | ASM1835013v1 ('ZK8'); Huazhong Agricultural University | chromosome | none (gene models at CGD, below) |
| `solanum_lycopersicum` | 1416wtT49, 29wtT29 | [GCF_036512215.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036512215.1/) | SLM_r2.1 ('Micro-Tom'); Kazusa DNA Research Institute | chromosome | RefSeq GCF_036512215.1-RS_2024_09, 28,611 protein-coding genes |

**Publications and why each assembly was chosen**

- ***Euonymus* (proxy).** Christenhusz MJM, Leitch IJ, Fay MF, et al. (2026) The genome sequence
  of *Euonymus europaeus* L., 1753 (Celastrales: Celastraceae). *Wellcome Open Res* 11:313.
  [doi:10.12688/wellcomeopenres.26648.1](https://doi.org/10.12688/wellcomeopenres.26648.1).
  NCBI has no *E. japonicus* assembly; this is the only one in the genus (plus its alternate
  haplotype). With no annotation, these samples give mapping fractions but no host count matrix.
- **Sweet orange.** Wu B, Yu Q, Deng Z, Duan Y, Luo F, Gmitter F (2023; online 2022) A
  chromosome-level phased genome enabling allele-level studies in sweet orange: a case study on
  citrus Huanglongbing tolerance. *Hortic Res* 10:uhac247.
  [doi:10.1093/hr/uhac247](https://doi.org/10.1093/hr/uhac247). This is the DVS assembly; its
  haplotype B is GCA_022201065.1 (DVS_B1.0, Clemson annotation, 27,940 genes). It is RefSeq
  annotated, and all sweet orange cultivars are somatic mutants of one hybrid, so Valencia
  stands in for Hamlin (CITRUS_HOST_PLAN.md §1.1, §3.1).
- **Carrizo citrange mapped to sweet orange** (Makefile comment on the `HOST_*CC*` lines). No
  Carrizo genome exists. Published Carrizo RNA-seq maps to *C. sinensis* too: Afzal Naveed Z,
  Huguet-Tapia J, Ali G (2019) Transcriptome profile of Carrizo citrange roots in response to
  *Phytophthora parasitica* infection. *J Plant Interact* 14:187-204.
  [doi:10.1080/17429145.2019.1609106](https://doi.org/10.1080/17429145.2019.1609106).
  In-house test (2026-09-08): the CC libraries align at 82-85 %, against 92 % for the native
  sweet-orange gall, and assign 92 % as many pairs to genes. Valid for CC-vs-CC contrasts only;
  a two-parent `carrizo` reference is planned (CITRUS_HOST_PLAN.md §4A).
- ***Brassica juncea*.** Kang L, Qian L, Zheng M, et al. (2021) Genomic insights into the origin,
  domestication and diversification of *Brassica juncea*. *Nat Genet* 53:1392-1402.
  [doi:10.1038/s41588-021-00922-y](https://doi.org/10.1038/s41588-021-00922-y). A yellow-seeded
  assembly from the same submitter ("Huangzi" = yellow seed). None of the 5 *B. juncea*
  assemblies in NCBI carries an NCBI annotation, so M26 gives mapping fractions only (it isn't
  in `ANNOT_HOSTS`).
- **Papaya.** Ming R, Hou S, Feng Y, et al. (2008) The draft genome of the transgenic tropical
  fruit tree papaya (*Carica papaya* Linnaeus). *Nature* 452:991-996.
  [doi:10.1038/nature06856](https://doi.org/10.1038/nature06856). RefSeq annotated.
- ***Poncirus* ZK8.** Huang Y, Xu Y, Jiang X, et al. (2021) Genome of a citrus rootstock and
  global DNA demethylation caused by heterografting. *Hortic Res*.
  [doi:10.1038/s41438-021-00505-2](https://doi.org/10.1038/s41438-021-00505-2) (BioProject
  PRJNA554539; the paper cited for this assembly on the
  [CGD ZK8 v1.0 page](https://www.citrusgenomedb.org/Analysis/1829479)). It is the only
  *Poncirus* assembly in NCBI. Gene models: CGD `jb2/data/Ptri_ZK8_v1.sorted.gff.gz`, used with
  CGD's own `Ptri_ZK8_v1.fasta.gz` (URLs in CITRUS_HOST_PLAN.md §3.1). A different assembly,
  the UF/JGI genome, is Peng Z, Bredeson JV, Wu GA, et al. (2020) *Plant J* 104:1215-1232,
  [doi:10.1111/tpj.14993](https://doi.org/10.1111/tpj.14993) (see §12).
- **Tomato.** Shirasawa K, Ariizumi T (2024) Near-complete genome assembly of tomato (*Solanum
  lycopersicum*) cultivar Micro-Tom. *Plant Biotechnol* 41:367-374.
  [doi:10.5511/plantbiotechnology.24.0522a](https://doi.org/10.5511/plantbiotechnology.24.0522a).
  RefSeq annotated. The cultivar of the gall plants isn't recorded **(ask)**.

## 6. Pipeline tools (Makefile steps 0-5)

Module versions are the SCINet Atlas defaults in the Makefile; on Ceres, `local.mk` overrides
samtools and subread. Conda versions come from `environment.yml`.

| tool | version | used for | citation | link |
|---|---|---|---|---|
| md5sum (GNU coreutils) | system | `make verify` | — | <https://www.gnu.org/software/coreutils/> |
| FastQC | 0.12.1 (module) | `make qc` | Andrews S (2010) FastQC: a quality control tool for high throughput sequence data | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/> |
| MultiQC | 1.35 (conda) | `make qc` | Ewels P, Magnusson M, Lundin S, Käller M (2016) *Bioinformatics* 32:3047-3048. [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) | <https://multiqc.info> |
| fastp | 1.3.6 (conda), default settings | `make trim` | Chen S (2025) fastp 1.0. *iMeta* 4:e70078. [doi:10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078); Chen S, Zhou Y, Chen Y, Gu J (2018) *Bioinformatics* 34:i884-i890. [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560) | <https://github.com/OpenGene/fastp> |
| Biopython | conda | `scripts/convert_refs.py` (.gb → FASTA + GFF3) | Cock PJA, Antao T, Chang JT, et al. (2009) *Bioinformatics* 25:1422-1423. [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) | <https://biopython.org> |
| NCBI Datasets API v2 | — | `make host-genomes` | O'Leary et al. 2024 (§5) | <https://www.ncbi.nlm.nih.gov/datasets/> |
| HISAT2 | 2.2.1 (module) | `make combined` (`hisat2-build`), `make align` (`--dta`) | Kim D, Paggi JM, Park C, Bennett C, Salzberg SL (2019) *Nat Biotechnol* 37:907-915. [doi:10.1038/s41587-019-0201-4](https://doi.org/10.1038/s41587-019-0201-4); Kim D, Langmead B, Salzberg SL (2015) *Nat Methods* 12:357-360. [doi:10.1038/nmeth.3317](https://doi.org/10.1038/nmeth.3317) | <https://daehwankimlab.github.io/hisat2/manual/> |
| SAMtools | 1.21 Atlas / 1.17 Ceres | sort, index, `idxstats` (organism fractions) | Danecek P, Bonfield JK, Liddle J, et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008. [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) | <https://www.htslib.org> |
| featureCounts (Subread) | 2.0.6 Atlas / 2.0.4 Ceres | `make counts` | Liao Y, Smyth GK, Shi W (2014) *Bioinformatics* 30:923-930. [doi:10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656) | <https://subread.sourceforge.net> |
| GNU Make | system | pipeline driver | — | <https://www.gnu.org/software/make/> |

**Why these settings**

- **HISAT2 rather than STAR:** HISAT2 needs little memory (Kim et al. 2015), which suits large
  plant hosts on a standard node (PIPELINE.md step 4). Wang et al. 2021 (§1) used HISAT2 with a
  concatenated plant + bacterial genome.
- **`--dta`:** reports alignments tailored for transcript assemblers such as StringTie, and
  requires longer anchors for novel splice sites (HISAT2 manual). Spliced alignment stays on
  for the bacterial contigs; `scripts/splice_rates.sh` → `logs/splice_rates.tsv` measures
  what that costs.
- **Single-copy host references, no allele-aware mapping.** Reads from gene copies that differ
  from the reference align less often (reference bias). The bias grows with divergence, and
  allowing more mismatches reduces it: Stevenson KR, Coolon JD, Wittkopp PJ (2013) Sources of
  bias in measures of allele-specific expression derived from RNA-seq data aligned to a single
  reference genome. *BMC Genomics* 14:536.
  [doi:10.1186/1471-2164-14-536](https://doi.org/10.1186/1471-2164-14-536); Degner JF, Marioni
  JC, Pai AA, et al. (2009) Effect of read-mapping biases on detecting allele-specific expression
  from RNA-sequencing data. *Bioinformatics* 25:3207-3212.
  [doi:10.1093/bioinformatics/btp579](https://doi.org/10.1093/bioinformatics/btp579).
  Gene-level counts compared within one host genotype share the bias, so it is not corrected.
  `make mapstats` records its size per sample (PIPELINE.md step 5b). The remedies below are for
  questions about which gene copy is expressed, and aren't used yet:
  - variant-aware index: `hisat2-build --snp` (Kim et al. 2019);
  - remap-and-filter: van de Geijn B, McVicker G, Gilad Y, Pritchard JK (2015) WASP.
    *Nat Methods* 12:1061-1063. [doi:10.1038/nmeth.3582](https://doi.org/10.1038/nmeth.3582);
    Castel SE, Levy-Moonshine A, Mohammadi P, et al. (2015) *Genome Biol* 16:195.
    [doi:10.1186/s13059-015-0762-6](https://doi.org/10.1186/s13059-015-0762-6);
  - references holding both gene copies: the DVS phased assembly (Wu et al., §5); the planned
    two-parent Carrizo reference; Rozowsky J, Abyzov A, Wang J, et al. (2011) AlleleSeq.
    *Mol Syst Biol* 7:522. [doi:10.1038/msb.2011.54](https://doi.org/10.1038/msb.2011.54);
    Kaminow B, Ballouz S, Gillis J, Dobin A (2022) *Genome Res* 32:738-749.
    [doi:10.1101/gr.275613.121](https://doi.org/10.1101/gr.275613.121);
  - assigning reads to subgenomes in polyploids (*B. juncea*): Akama S, Shimizu-Inatsugi R,
    Shimizu KK, Sese J (2014) HomeoRoq. *Nucleic Acids Res* 42:e46.
    [doi:10.1093/nar/gkt1376](https://doi.org/10.1093/nar/gkt1376); Page JT, Gingle AR,
    Udall JA (2013) PolyCat. *G3* 3:517-525.
    [doi:10.1534/g3.112.005298](https://doi.org/10.1534/g3.112.005298); Kuo T, Hatakeyama M,
    Tameshige T, et al. Homeolog expression quantification methods for allopolyploids.
    *Brief Bioinform* 21:395-407. [doi:10.1093/bib/bby121](https://doi.org/10.1093/bib/bby121).
- **featureCounts `-p --countReadPairs -s 0 -t gene -g ID`:** counts fragments rather than reads,
  unstranded (§3), at gene level. Multi-mapping reads are not counted (the featureCounts
  default); PIPELINE.md step 5 flags that T-DNA genes can be affected.
- **Gene rows for the bacterial GFF3:** SnapGene and RAST GenBank files carry CDS, tRNA and rRNA
  features but no `gene` features, so `convert_refs.py` writes a `gene` row for each one. A
  feature crossing the origin of a circular replicon is written as one row per part with a
  shared ID, and featureCounts merges rows sharing an ID into one meta-feature (Subread User's
  Guide). This change was uncommitted in the working checkout on 2026-09-14.
- **Normalise each organism separately** when DE starts: the plant : bacterium ratio varies by
  sample (PIPELINE.md step 6).

**Host functional annotation (Makefile step 7, `make annotate`; added 2026-09-15)**

| tool / data | version | citation | link |
|---|---|---|---|
| eggNOG-mapper | 2.1.13 (venv, `scripts/emapper-env.sh`) | Cantalapiedra CP, Hernandez-Plaza A, Letunic I, Bork P, Huerta-Cepas J (2021) eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol* msab293. [doi:10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293) | <https://github.com/eggnogdb/eggnog-mapper> |
| eggNOG 5.0 (emapperdb 5.0.2) | downloaded 2026-09-14 | Huerta-Cepas J, Szklarczyk D, Heller D, et al. (2019) eggNOG 5.0. *Nucleic Acids Res* 47:D309-D314. [doi:10.1093/nar/gky1085](https://doi.org/10.1093/nar/gky1085) | <http://eggnog5.embl.de/download/emapperdb-5.0.2/> |
| DIAMOND | 2.1.24 (module) | Buchfink B, Reuter K, Drost HG (2021) Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods* 18:366-368. [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x) | <https://github.com/bbuchfink/diamond> |
| Arabidopsis proteome | Ensembl Plants release 63, TAIR10 `pep.all` | citation to add | <https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/pep/> |
| host proteins | NCBI RefSeq `PROT_FASTA` | O'Leary et al. 2024 (§5) | NCBI Datasets API |

Settings: longest protein per gene (`scripts/longest_proteins.py`); eggNOG-mapper in DIAMOND mode,
`--tax_scope Viridiplantae`, default GO evidence filter; DIAMOND `blastp --more-sensitive -k 1
-e 1e-5` in both directions, with reciprocal best hits flagged (`scripts/gene_annotation.py`).
The eggNOG database is fetched from eggnog5.embl.de because emapper's own download host,
eggnogdb.embl.de, no longer resolves.

## 7. Checks and helper scripts

| script | what it does | output |
|---|---|---|
| `scripts/convert_refs.py` | GenBank → FASTA + GFF3 with `agro_` contig names (Biopython) | `references/agrobacterium/` |
| `scripts/merge_counts.py` | merges per-sample featureCounts tables, keeping or dropping `agro_` rows | `04_matrix/*.tsv` |
| `scripts/intron_fraction.sh` | re-counts BAMs with `-t exon` vs `-t gene`; intronic = difference in NoFeatures (polyA test, §3) | `logs/intron_fraction.tsv` |
| `scripts/splice_rates.sh` | share of alignments with `N` in the CIGAR, per organism (`samtools view`, `idxstats`) | `logs/splice_rates.tsv` |
| `scripts/mapstats.py` (`make mapstats`) | per-sample mapping QC for the paper: fastp retention, HISAT2 rates, primary reads per organism over all reads, host mismatch rate (`samtools stats`), featureCounts assignment | `logs/mapping_stats.tsv` |

## 8. Notebooks and figures

`notebooks/01_results_overview.ipynb` and `02_dataset_structure.ipynb` (rendered `.html` beside
them; figures go to `figures/`). They run in a separate venv (`scripts/viz-env.sh`, pins in
`notebooks/requirements.txt`) so the pipeline's pinned tool versions stay untouched.

| package | version | citation |
|---|---|---|
| JupyterLab / nbconvert / ipykernel | 4.5.10 / 7.17.1 / 6.31.0 | Kluyver T, Ragan-Kelley B, Pérez F, et al. (2016) Jupyter Notebooks – a publishing format for reproducible computational workflows. In *Positioning and Power in Academic Publishing*, 87-90. [doi:10.3233/978-1-61499-649-1-87](https://doi.org/10.3233/978-1-61499-649-1-87) |
| pandas | 2.3.3 | McKinney W (2010) Data structures for statistical computing in Python. *Proc 9th Python in Science Conf*, 56-61. [doi:10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a) |
| NumPy | 2.0.2 | Harris CR, Millman KJ, van der Walt SJ, et al. (2020) *Nature* 585:357-362. [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2) |
| SciPy | 1.13.1 | Virtanen P, Gommers R, Oliphant TE, et al. (2020) *Nat Methods* 17:261-272. [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2) |
| Matplotlib | 3.9.4 | Hunter JD (2007) *Comput Sci Eng* 9:90-95. [doi:10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55) |
| seaborn | 0.13.2 | Waskom ML (2021) *J Open Source Softw* 6:3021. [doi:10.21105/joss.03021](https://doi.org/10.21105/joss.03021) |

## 9. Compute and environments

- **Clusters:** USDA-ARS SCINet: Atlas (module names in the Makefile) and Ceres (partition `ceres`,
  account `small_grains`; `local.mk`, SLURM.md). Slurm jobs: `scripts/pipeline.slurm`
  (unattended `make` run), `scripts/conda-env.slurm` (builds the env).
  <https://scinet.usda.gov>. Acknowledgement required in manuscripts
  ([SCINet FAQ](https://scinet.usda.gov/support/faq)): "This research used resources provided by
  the SCINet project and/or the AI Center of Excellence of the USDA Agricultural Research
  Service, ARS project numbers 0201-88888-003-000D and 0201-88888-002-000D."
- **Conda env `rnaseq`** (`environment.yml`, channels conda-forge + bioconda): Grüning B, Dale R,
  Sjödin A, et al. (2018) Bioconda: sustainable and comprehensive software distribution for the
  life sciences. *Nat Methods* 15:475-476.
  [doi:10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7).
- **Analysis sessions:** Claude Code (Anthropic) on compute nodes (`scripts/claude-*.sh`,
  `scripts/claude-*.slurm`, SLURM.md).

## 10. Installed or planned, not yet used

Already in `environment.yml`, for later steps:

| tool | intended use | citation |
|---|---|---|
| STAR 2.7.11b | alternative aligner | Dobin A, Davis CA, Schlesinger F, et al. (2013) *Bioinformatics* 29:15-21. [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635) |
| Bowtie 2 2.5.4 | — | Langmead B, Salzberg SL (2012) *Nat Methods* 9:357-359. [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923) |
| SeqKit 2.10.0 | contig renaming (PIPELINE.md step 3; the Makefile uses convert_refs.py instead) | Shen W, Sipos B, Zhao L (2024) SeqKit2. *iMeta* 3:e191. [doi:10.1002/imt2.191](https://doi.org/10.1002/imt2.191); Shen W, Le S, Li Y, Hu F (2016) *PLoS ONE* 11:e0163962. [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962) |
| BEDTools | interval work | Quinlan AR, Hall IM (2010) *Bioinformatics* 26:841-842. [doi:10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033) |
| StringTie 3.0.0 | transcript assembly (why `--dta` is set) | Pertea M, Pertea GM, Antonescu CM, et al. (2015) *Nat Biotechnol* 33:290-295. [doi:10.1038/nbt.3122](https://doi.org/10.1038/nbt.3122); Kovaka S, Zimin AV, Pertea GM, et al. (2019) StringTie2. *Genome Biol* 20:278. [doi:10.1186/s13059-019-1910-1](https://doi.org/10.1186/s13059-019-1910-1) |
| Salmon 1.10.3 | per-parental-copy quantification option (CITRUS_HOST_PLAN.md §4A) | Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017) *Nat Methods* 14:417-419. [doi:10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197) |
| MEME Suite | promoter motifs (PIPELINE.md step 9) | Bailey TL, Johnson J, Grant CE, Noble WS (2015) *Nucleic Acids Res* 43:W39-W49. [doi:10.1093/nar/gkv416](https://doi.org/10.1093/nar/gkv416) |
| DESeq2 1.50.2 | differential expression (step 6) | Love MI, Huber W, Anders S (2014) *Genome Biol* 15:550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8) |
| tximport 1.38.2 | Salmon → gene counts | Soneson C, Love MI, Robinson MD (2015) *F1000Research* 4:1521. [doi:10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2) |
| edgeR 4.8.2 | differential expression | Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2025) edgeR v4. *Nucleic Acids Res* 53:gkaf018. [doi:10.1093/nar/gkaf018](https://doi.org/10.1093/nar/gkaf018); Robinson MD, McCarthy DJ, Smyth GK (2010) *Bioinformatics* 26:139-140. [doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616) |
| pheatmap, ggplot2, lftp | plotting; file transfer | CRAN / <https://lftp.yar.ru> |

Planned in PIPELINE.md and CITRUS_HOST_PLAN.md (links only until adopted):

- OrthoFinder (step 8): Emms DM, Kelly S (2019) *Genome Biol* 20:238.
  [doi:10.1186/s13059-019-1832-y](https://doi.org/10.1186/s13059-019-1832-y)
- JASPAR 2024 CORE plants (motif scans): Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R,
  et al. (2024) *Nucleic Acids Res* 52:D174-D182.
  [doi:10.1093/nar/gkad1059](https://doi.org/10.1093/nar/gkad1059)
- Liftoff (the old fallback for annotating *Poncirus*; superseded by the CGD ZK8 gene models):
  Shumate A, Salzberg SL (2021) *Bioinformatics* 37:1639-1643.
  [doi:10.1093/bioinformatics/btaa1016](https://doi.org/10.1093/bioinformatics/btaa1016)
- Citrus Genome Database (ZK8 FASTA/GFF, MCscan anchors `Csin_DVS_A_v1.Ptri_ZK8_v1.anchors.gz`):
  <https://www.citrusgenomedb.org>
- PlantRegMap / PlantTFDB (`Citrus_sinensis`), PLAZA dicots 5.0: see CITRUS_HOST_PLAN.md §3.2.
  Citations to add when used. (eggNOG-mapper and DIAMOND are now in use: §6, step 7.)
- Public healthy-tissue baselines, run lists in `citrus_baselines.tsv` (downloads from ENA,
  <https://www.ebi.ac.uk/ena/browser/view/ACCESSION>):
  - Carrizo: [PRJNA1216034](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1216034),
    [PRJNA668159](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA668159),
    [PRJNA839431](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA839431),
    [PRJNA1053671](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1053671)
  - *Poncirus*: [PRJDB41296](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB41296)
  - Sweet orange: [PRJNA599503](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA599503),
    [PRJNA778304](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA778304)

  Study notes and tiers are in CITRUS_HOST_PLAN.md §3.3.

## 11. Internal and collaborator documents

These are unpublished and gitignored; don't redistribute. Use them as context, not as
citations.

- `docs/050136 Agro resource paper-RT JT-rgs.docx`: resource-paper draft (strain 1416
  disarming, GAANTRY; CITRUS_HOST_PLAN.md §1.2).
- `docs/012226 Shatter Lab Agrobacterim_analysis_pipeline.docx`: collaborator pipeline notes
  (not reviewed for this file).
- `docs/090326 Agrobacteria Rating Table for Publication.xlsx`, `docs/Agrobacteria Rating Table.xlsx`:
  the 31-species host panel behind `references/plant_host/AVAILABILITY.md`.
- `docs/111725 Gall Visual Rating System.pdf`: gall rating scale (0 no visible nodule, 1 dead or
  not actively growing, 2 live and actively growing).
- `docs/Agro RNAseq DATA set`: not reviewed.
- `Strain 1416/`, `Strain 29/`, `C58_ALIGNED/`: strain source files (§4).
- PDFs of published papers cited above: Alabed 2023 and 2024, Huo 2019 (filed as "Hou"), RAST.

## 12. Discrepancies found while compiling

1. **Wrong *Poncirus* ZK8 citation in CITRUS_HOST_PLAN.md §3.1.** It says to cite Peng et al.
   2020 for ZK8. That paper (UF / JGI authors) describes a different assembly, probably the
   JGI v1.3.1 genome also listed on CGD. ZK8 (GCA_018350135.1, PRJNA554539, HZAU) is Huang et
   al. 2021 *Hortic Res* ([doi:10.1038/s41438-021-00505-2](https://doi.org/10.1038/s41438-021-00505-2)),
   the paper CGD cites. **Fixed in the plan 2026-09-15.**
2. **1416 pAt1 size.** Our `FINAL.gb` conversion gives 589,645 bp (references/README.md), but
   Alabed et al. 2023 say 519,735 bp. The same 519,735 bp figure appears as the AT plasmid of
   1D159 in Huo et al. 2019, which suggests a copy error in the paper. Confirm with the
   collaborators before publishing replicon statistics.
3. **C58 accessions in PIPELINE.md** list AE007869 / AE007870 / AE007872 and leave out AE007871
   (the Ti plasmid). references/README.md lists NC_003065 for it.
4. **DVS_A1.0 gene count:** NCBI reports 23,556 protein-coding genes; CITRUS_HOST_PLAN.md §3.1
   said 23,566. **Fixed 2026-09-15.**
5. **Unknowns to ask the collaborators about:** strain 29's designation and publication; the
   tomato cultivar (T29, T49); the *Euonymus* cultivars (Eu, and golden Geu) against the
   *E. europaeus* proxy; the library kit (polyA and unstranded are inferred from the data).
6. **González-Mula et al. 2018:** how they split plant and bacterial reads is unverified (§2).
