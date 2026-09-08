# Gall Dual RNA-seq Pipeline

> Steps 0-5 are automated in `Makefile` (`make help`); sample metadata in `samples.tsv`.
> Run heavy targets on a compute node: `srun --cpus-per-task=16 --mem=64G make -j2 align`

## Existing data

`30-1348328766/00_fastq/` — Azenta/GeneWiz, Illumina NovaSeq X, 2x150 PE, 14 samples (R1+R2 + md5 each).
Sample naming: `<strain><genotype><host accession>`. Gall tissue: mixed plant host + Agrobacterium reads.

**Strain 1416 (8 samples)**
- wt on hosts: CC54, CC547, Eu635, Geu182, M26, T49
- engineered events G-19, G-30 on host CC547 (wt/G-19/G-30 trio; n=1 per genotype — limits DE inference)

**Strain 29 (6 samples, all wt)**
- hosts: Eu635, Geu182, HC83, M26, P46, T29

9 distinct plant host accessions: CC54, CC547, Eu635, Geu182, HC83, M26, P46, T29, T49
(Eu635, Geu182, M26 shared between strains). Species TBD.

**Strain assemblies in hand** (on local Windows machine, GenBank format w/ annotation — need transfer to cluster + .gb → FASTA/GFF3 conversion):
- Strain 1416: circular + linear chromosome, pAT1, pAT2, pTi (all FINAL.gb)
- Strain 29: circular + linear chromosome, pAT1, pAT2, pTi, T-DNA (.gb)
- C58 published reference (AE007869/AE007870/AE007872) — for comparison only

**Still needed**
- Plant host genome(s) + annotation — species/accession identity of the 9 host codes unknown; ask collaborator
- G-19 / G-30 construct or event sequences (check "1416 10294 bp fragment (LC).gb")

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
- Data: plant host genome FASTA + strain FASTA (chromosomes + plasmids), merged GFF/GTF
- Tools: Biopython for .gb → FASTA/GFF3, `seqkit` (module) for contig renaming, `hisat2-build` or `STAR --runMode genomeGenerate`
- Prefix contigs by organism (e.g. `plant_`, `agro_`) — no name collisions
- One combined index per host-genome x strain pair (2 if all hosts share one reference genome)
- Add G-19/G-30 event sequences as extra contigs for those samples

### 4. Align (competitive mapping)
- Data: trimmed fastq + combined index for the sample's strain
- Tools: `hisat2` (module; STAR needs big-mem node for wheat-scale host)
- Output: sorted, indexed BAM (`samtools`, module)
- QC: `samtools idxstats` aggregated by contig prefix → plant vs agro fraction per sample

### 5. Count
- Data: BAMs + merged GFF/GTF
- Tools: `featureCounts` (subread module)
- Split count matrix by organism prefix → plant counts, agro counts
- Flag T-DNA genes (multimap between pTi and transformed plant genome)

### 6. Differential expression
- Data: split count matrices + sample metadata (strain, genotype, host)
- Tools: R in conda env `rnaseq`: DESeq2, tximport, edgeR, pheatmap, ggplot2
- Normalize each organism SEPARATELY (own size factors — plant/agro ratio varies per sample)
- Key contrasts: wt vs G-19/G-30 on host CC547 (n=1 each — fold-change ranking only, or treat
  G-19+G-30 as 2 reps of "engineered"); strain 1416 vs 29 on shared hosts (Eu635, Geu182, M26)

### 7. Cross-strain orthology (after Haryono et al. 2019 — see HARYONO2019_PIPELINE.md)
- Data: strain 1416 + strain 29 protein sequences (from .gb annotations); optionally C58
- Tools: `orthofinder` (conda, to install when needed)
- Cluster homologs; classify each ortholog pair by DE status (up/up, up/absent, etc.)
- Analyze per replicon — expect circular chr most conserved, linear chr + plasmids divergent

### 8. Promoter / motif analysis
- Data: 600 bp upstream regions of DE genes (grouped up/down per strain); stranded coverage
- Tools: MEME suite (conda env `rnaseq`): `fimo` for known motifs, `meme`/`streme` for ab initio;
  optionally Rockhopper/ANNOgesic for coverage-based TSS/operon inference
- Targeted scan first: VirG box consensus RTTDCAWWTGHAAY, <=3 mismatches; group genes <50 bp
  apart on same strand as one operon (Haryono et al. protocol)
- EXPECTATION CHECK: Haryono et al. found NO robust novel motif ab initio even with clean
  cultures + replicates — only the known VirG box. Novel promoter discovery likely needs
  cross-strain phylogenetic footprinting (align upstream regions of orthologs) or a dedicated
  TSS library (dRNA-seq/Cappable-seq)
