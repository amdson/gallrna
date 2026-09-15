# ============================================================================
# Gall dual RNA-seq pipeline (see PIPELINE.md for rationale)
# Run from this directory. Heavy targets belong on a compute node, e.g.:
#     srun --cpus-per-task=16 --mem=64G make -j2 align
# `make help` lists targets. THREADS and STRAND can be overridden: make align THREADS=16
# HOSTSEL="poncirus_trifoliata citrus_sinensis" (or SAMPLES=...) restricts to a subset of samples.
#
# ---------------------------- INPUTS ----------------------------------------
# FILES (must exist locally):
#   30-1348328766/00_fastq/<sample>_R{1,2}_001.fastq.gz (+ .md5)   raw reads, GeneWiz 30-1348328766
#   samples.tsv                                    sample metadata (strain/genotype/host/age)
#   "Strain 1416/1416 * FINAL.gb"                  strain 1416 assembly, annotated (collaborator USB)
#   "Strain 29/.../Xtra RAST download/357.66*.gbk" strain 29 assembly, RAST-annotated (collaborator USB)
#   scripts/convert_refs.py, scripts/merge_counts.py
# LINKS (downloaded automatically by `make host-genomes`, NCBI datasets API):
#   https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<ACC>/download
#     euonymus_japonicus_proxy  GCA_963580455.1  (E. europaeus - congener proxy, no annotation)
#     citrus_sinensis           GCF_022201045.2  (DVS_A1.0, annotated)
#     brassica_juncea           GCA_018703725.1  (no annotation)
#     carica_papaya             GCF_000150535.2  (Papaya1.0, annotated)
#     poncirus_trifoliata       GCA_018350135.1  (no annotation)
#     solanum_lycopersicum      GCF_036512215.1  (SLM_r2.1, annotated)
#   (downloaded by `make annotate`, step 7)
#   same NCBI API, include_annotation_type=PROT_FASTA   host RefSeq proteins
#   https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/.../Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
#   http://eggnog5.embl.de/download/emapperdb-5.0.2/    eggNOG-mapper database (~50 GB)
# TOOLS: cluster modules hisat2/2.2.1 samtools/1.21 subread/2.0.6 fastqc/0.12.1 diamond;
#        conda env `rnaseq` (fastp, multiqc, python+biopython) - see environment.yml;
#        venv ~/.venvs/emapper (eggNOG-mapper) - see scripts/emapper-env.sh
# ============================================================================

SHELL := /bin/bash
.SHELLFLAGS := -lc
.SECONDEXPANSION:
.SECONDARY:
.DEFAULT_GOAL := help

THREADS ?= 8
STRAND  ?= 0            # featureCounts -s: 0 until strandedness confirmed; dUTP kits are usually 2
CONDA   ?= $(HOME)/.conda/envs/rnaseq/bin
# module names as found on SCINet Atlas - override for other clusters (e.g. Ceres):
#   make align MOD_HISAT2=hisat2 MOD_SAMTOOLS=samtools ...
MOD_HISAT2   ?= hisat2/2.2.1
MOD_SAMTOOLS ?= samtools/1.21
MOD_SUBREAD  ?= subread/2.0.6
MOD_FASTQC   ?= fastqc/0.12.1
MOD_DIAMOND  ?= diamond
# per-cluster overrides (gitignored), e.g. local.mk on Ceres sets MOD_SAMTOOLS/MOD_SUBREAD
-include local.mk
RAW     := 30-1348328766/00_fastq
ALL_SAMPLES := $(shell tail -n +2 samples.tsv | cut -f1)

# --- sample -> host genome (from samples.tsv / gall tissue inventory) -------
HOST_1416wtEu635   = euonymus_japonicus_proxy
HOST_29wtEu635     = euonymus_japonicus_proxy
HOST_29wtGeu182    = euonymus_japonicus_proxy
HOST_1416wtGeu182  = euonymus_japonicus_proxy
HOST_29wtHC83      = citrus_sinensis
HOST_1416wtM26     = brassica_juncea
HOST_29wtM26       = brassica_juncea
HOST_29wtP46       = carica_papaya
# CC = Carrizo citrange, an F1 hybrid of C. sinensis x P. trifoliata (confirmed
# 2026-09-14; samples.tsv updated to match). No Carrizo genome exists, so
# the CC galls are aligned to one parent, C. sinensis, and counted on its gene
# models, as in published Carrizo RNA-seq, e.g. Afzal Naveed, Huguet-Tapia & Ali
# 2019, J Plant Interact 14:187-204, doi:10.1080/17429145.2019.1609106 (roots
# mapped to two C. sinensis genomes; they report only 55-73% transcriptome
# coverage and fell back to de novo assembly). Here the CC libraries align at
# 82-85% vs 92% for the native sweet-orange gall 29wtHC83, and a cross-species
# test (2026-09-08) assigned 92% as many pairs to genes as that native sample.
# All four CC galls share the hybrid genotype, so the bias from P. trifoliata
# alleles mapping with mismatches is common to every CC sample and cancels in
# CC-vs-CC contrasts; it does not cancel against galls on other hosts.
# If parent of origin or Poncirus-only genes matter, add P. trifoliata as a
# second haplotype instead: the ZK8 assembly under poncirus_trifoliata IS
# annotated at the Citrus Genome Database (25,680 genes;
# citrusgenomedb.org/jb2/data/Ptri_ZK8_v1.sorted.gff.gz), but its contigs are
# named chr1_ZK8.., so pair it with CGD's jb2/data/Ptri_ZK8_v1.fasta.gz, not
# NCBI's GCA_018350135.1 (CM031384.1..) used now.
# CITRUS_HOST_PLAN.md (sections 1.1, 3.1, 4A) plans that two-parent `carrizo` host:
# DVS_A1.0 + ZK8 concatenated, counted per parental copy, summed per gene pair.
HOST_1416wtCC547   = citrus_sinensis
HOST_1416wtCC54    = citrus_sinensis
HOST_1416G-19CC547 = citrus_sinensis
HOST_1416G-30CC547 = citrus_sinensis
HOST_1416wtT49     = solanum_lycopersicum
HOST_29wtT29       = solanum_lycopersicum

strainof = $(if $(filter 1416%,$1),strain_1416,strain_29)
comboof  = $(call strainof,$1)__$(HOST_$1)

HOSTS       := euonymus_japonicus_proxy citrus_sinensis brassica_juncea carica_papaya poncirus_trifoliata solanum_lycopersicum
ANNOT_HOSTS := citrus_sinensis carica_papaya solanum_lycopersicum
hostsamples = $(strip $(foreach s,$(ALL_SAMPLES),$(if $(filter $(HOST_$s),$1),$s)))
# restrict the run to some hosts (or pass SAMPLES=... directly), e.g.
#   make fractions HOSTSEL="poncirus_trifoliata citrus_sinensis"
HOSTSEL ?=
SAMPLES := $(if $(HOSTSEL),$(foreach h,$(HOSTSEL),$(call hostsamples,$h)),$(ALL_SAMPLES))
COMBOS      := $(sort $(foreach s,$(SAMPLES),$(call comboof,$s)))

ACC_euonymus_japonicus_proxy := GCA_963580455.1
ACC_citrus_sinensis          := GCF_022201045.2
ACC_brassica_juncea          := GCA_018703725.1
ACC_carica_papaya            := GCF_000150535.2
ACC_poncirus_trifoliata      := GCA_018350135.1
ACC_solanum_lycopersicum     := GCF_036512215.1

# ---------------------------- phony targets ---------------------------------
.PHONY: help all verify qc trim strain-refs host-genomes refs combined align fractions counts matrices mapstats annotate

help:
	@echo "targets: verify qc trim refs (strain-refs host-genomes) combined align fractions counts matrices mapstats all annotate"

all: verify qc fractions matrices

verify: logs/md5.ok
qc: qc/multiqc_report.html
trim: $(foreach s,$(SAMPLES),01_trim/$(s)_R1.fastq.gz)
strain-refs: references/agrobacterium/strain_1416/strain_1416.fasta references/agrobacterium/strain_29/strain_29.fasta
host-genomes: $(foreach h,$(HOSTS),references/plant_host/$(h)/$(h).fasta)
refs: strain-refs host-genomes
combined: $(foreach c,$(COMBOS),references/combined/$(c).hisat2.ok)
align: $(foreach s,$(SAMPLES),02_align/$(s).bam)
fractions: logs/mapping_fractions.tsv
counts: $(foreach s,$(SAMPLES),03_counts/$(s).txt)
matrices: 04_matrix/agro_strain_1416.tsv 04_matrix/agro_strain_29.tsv \
          $(foreach h,$(ANNOT_HOSTS),04_matrix/plant_$(h).tsv)
mapstats: logs/mapping_stats.tsv

# ---------------------------- 0. verify transfer ----------------------------
logs/md5.ok:
	mkdir -p logs
	cd $(RAW) && md5sum -c --quiet *.md5
	touch $@

# ---------------------------- 1. raw QC -------------------------------------
qc/multiqc_report.html: $(wildcard $(RAW)/*.fastq.gz)
	mkdir -p qc/fastqc
	module load $(MOD_FASTQC) && fastqc -t $(THREADS) -o qc/fastqc $(RAW)/*.fastq.gz
	$(CONDA)/multiqc -f -o qc qc/fastqc

# ---------------------------- 2. trim ---------------------------------------
01_trim/%_R1.fastq.gz 01_trim/%_R2.fastq.gz: $(RAW)/%_R1_001.fastq.gz $(RAW)/%_R2_001.fastq.gz
	mkdir -p 01_trim
	$(CONDA)/fastp -w $(THREADS) -i $(RAW)/$*_R1_001.fastq.gz -I $(RAW)/$*_R2_001.fastq.gz \
	  -o 01_trim/$*_R1.fastq.gz -O 01_trim/$*_R2.fastq.gz \
	  -j 01_trim/$*.fastp.json -h 01_trim/$*.fastp.html 2> 01_trim/$*.fastp.log

# ---------------------------- 3a. strain references -------------------------
# NOTE: .gb sources contain spaces so make cannot track them as prerequisites;
# delete the outputs (or touch the .gb files' dir) to force reconversion.
references/agrobacterium/strain_1416/strain_1416.fasta: scripts/convert_refs.py
	mkdir -p $(dir $@)
	$(CONDA)/python scripts/convert_refs.py references/agrobacterium/strain_1416/strain_1416 \
	  "agro_1416_circ=Strain 1416/1416 Circular Chromosome FINAL.gb" \
	  "agro_1416_lin=Strain 1416/1416 Linear Chromosome FINAL.gb" \
	  "agro_1416_pAt1=Strain 1416/1416 pAT1 plasmid FINAL.gb" \
	  "agro_1416_pAt2=Strain 1416/1416 pAT2 plasmid FINAL.gb" \
	  "agro_1416_pTi=Strain 1416/1416 Ti plasmid FINAL.gb"

references/agrobacterium/strain_29/strain_29.fasta: scripts/convert_refs.py
	mkdir -p $(dir $@)
	$(CONDA)/python scripts/convert_refs.py references/agrobacterium/strain_29/strain_29 \
	  "agro_29_circ=Strain 29/Agrobacterium Strain 29 Raw data/Xtra RAST download/357.660 CC.gbk" \
	  "agro_29_lin=Strain 29/Agrobacterium Strain 29 Raw data/Xtra RAST download/357.661 LC.gbk" \
	  "agro_29_pAt1=Strain 29/Agrobacterium Strain 29 Raw data/Xtra RAST download/357.662 pAT1.gbk" \
	  "agro_29_pAt2=Strain 29/Agrobacterium Strain 29 Raw data/Xtra RAST download/357.663 pAT2.gbk" \
	  "agro_29_pTi=Strain 29/Agrobacterium Strain 29 Raw data/Xtra RAST download/357.664 pTi.gbk"

# ---------------------------- 3b. host genomes (NCBI) -----------------------
references/plant_host/%.fasta:
	$(eval NAME := $(notdir $*))
	mkdir -p references/plant_host/$(NAME)
	curl -sL "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$(ACC_$(NAME))/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF" \
	  -o references/plant_host/$(NAME)/$(NAME).zip
	unzip -oq references/plant_host/$(NAME)/$(NAME).zip -d references/plant_host/$(NAME)/tmp
	find references/plant_host/$(NAME)/tmp -name "*.fna" -exec mv {} references/plant_host/$(NAME)/$(NAME).fasta \;
	-find references/plant_host/$(NAME)/tmp -name "genomic.gff" -exec mv {} references/plant_host/$(NAME)/$(NAME).gff3 \;
	rm -rf references/plant_host/$(NAME)/tmp references/plant_host/$(NAME)/$(NAME).zip

# ---------------------------- 3c. combined refs + index ---------------------
# combo name = <strain>__<host>; agro contigs carry the agro_ prefix, plant
# contigs keep their NCBI accessions (organism split downstream keys on agro_).
references/combined/%.fasta: references/agrobacterium/$$(word 1,$$(subst __, ,$$*))/$$(word 1,$$(subst __, ,$$*)).fasta \
                             references/plant_host/$$(word 2,$$(subst __, ,$$*))/$$(word 2,$$(subst __, ,$$*)).fasta
	mkdir -p references/combined
	cat $^ > $@
	cat references/agrobacterium/$(word 1,$(subst __, ,$*))/$(word 1,$(subst __, ,$*)).gff3 > references/combined/$*.gff3
	@hostgff=references/plant_host/$(word 2,$(subst __, ,$*))/$(word 2,$(subst __, ,$*)).gff3; \
	  if [ -f "$$hostgff" ]; then grep -v "^#" "$$hostgff" >> references/combined/$*.gff3; \
	  else echo "NOTE: no host annotation for $* - combined GFF is agro-only"; fi

references/combined/%.hisat2.ok: references/combined/%.fasta
	module load $(MOD_HISAT2) && hisat2-build -p $(THREADS) $< references/combined/$* > references/combined/$*.build.log 2>&1
	touch $@

# ---------------------------- 4. align (competitive) ------------------------
02_align/%.bam: 01_trim/%_R1.fastq.gz 01_trim/%_R2.fastq.gz references/combined/$$(call comboof,$$*).hisat2.ok
	mkdir -p 02_align
	module load $(MOD_HISAT2) $(MOD_SAMTOOLS) && \
	hisat2 -p $(THREADS) --dta -x references/combined/$(call comboof,$*) \
	  -1 01_trim/$*_R1.fastq.gz -2 01_trim/$*_R2.fastq.gz 2> 02_align/$*.hisat2.log \
	  | samtools sort -@ $(THREADS) -o $@ - && samtools index $@

02_align/%.idxstats: 02_align/%.bam
	module load $(MOD_SAMTOOLS) && samtools idxstats $< > $@

logs/mapping_fractions.tsv: $(foreach s,$(SAMPLES),02_align/$(s).idxstats)
	mkdir -p logs
	@echo -e "sample\tagro_reads\tplant_reads\tagro_frac" > $@
	@for s in $(SAMPLES); do \
	  awk -v s=$$s '{if ($$1 ~ /^agro_/) a += $$3; else if ($$1 != "*") p += $$3} \
	    END {printf "%s\t%d\t%d\t%.4f\n", s, a, p, a/(a+p+1e-9)}' 02_align/$$s.idxstats; \
	done >> $@
	@cat $@

# ---------------------------- 5. count --------------------------------------
03_counts/%.txt: 02_align/%.bam
	mkdir -p 03_counts
	module load $(MOD_SUBREAD) && \
	featureCounts -p --countReadPairs -T $(THREADS) -s $(STRAND) -t gene -g ID \
	  -a references/combined/$(call comboof,$*).gff3 -o $@ $< 2> 03_counts/$*.log

04_matrix/agro_%.tsv: $$(foreach s,$$(SAMPLES),$$(if $$(filter $$(call strainof,$$s),$$*),03_counts/$$s.txt))
	mkdir -p 04_matrix
	$(CONDA)/python scripts/merge_counts.py $@ keep agro_ $^

04_matrix/plant_%.tsv: $$(foreach s,$$(call hostsamples,$$*),03_counts/$$s.txt)
	mkdir -p 04_matrix
	$(CONDA)/python scripts/merge_counts.py $@ drop agro_ $^

# ---------------------------- 5b. mapping statistics ------------------------
# Per-sample mapping QC for the paper's supplementary table: read retention, alignment
# (unique/multi), host vs bacterial reads (primary alignments, as a share of all reads
# sequenced), host mismatch rate and featureCounts assignment. It records that the
# single-copy host references gave reasonable results even though the hosts are
# heterozygous, hybrid or polyploid, and it bounds plant->bacterium cross-mapping.
# Column guide and how to read it: PIPELINE.md step 5b. samtools stats reads every host
# alignment, so run it on a compute node, e.g.
#   srun --cpus-per-task=16 --mem=16G make -j4 mapstats THREADS=4
02_align/%.mapstats.tsv: 02_align/%.bam 03_counts/%.txt scripts/mapstats.py
	module load $(MOD_SAMTOOLS) && python3 scripts/mapstats.py $* $(call comboof,$*) $< \
	  02_align/$*.hisat2.log 01_trim/$*.fastp.json 03_counts/$*.txt $(THREADS) > $@.tmp
	mv $@.tmp $@

logs/mapping_stats.tsv: $(foreach s,$(SAMPLES),02_align/$(s).mapstats.tsv)
	mkdir -p logs
	awk 'FNR > 1 || NR == 1' $^ > $@
	@echo "$@: $$(($$(wc -l < $@) - 1)) samples"

# 6. differential expression: not automated yet - needs strandedness + design
# decisions (see PIPELINE.md step 6). Matrices in 04_matrix/ are DESeq2-ready.

# ---------------------------- 7. functional annotation ----------------------
# One protein per gene (longest RefSeq isoform, named by the gene ID featureCounts
# writes) -> eggNOG-mapper (GO, KEGG, EC, Pfam, description) + DIAMOND best hit in
# Arabidopsis (AGI + symbol; the crown-gall marker lists in CITRUS_HOST_PLAN.md are
# AGIs). Output: references/plant_host/<host>/<host>.genes.tsv, one row per counted
# gene, so it joins straight onto 04_matrix/plant_<host>.tsv.
# Needs scripts/emapper-env.sh once. The eggNOG database (~50 GB) downloads on first
# run into EGGNOG_DATA - on Ceres keep it off home (30 GB quota): make
# references/plant_host/eggnog_data a symlink into 90daydata like the host dirs, or
# set EGGNOG_DATA in local.mk. Works for any RefSeq-annotated host, e.g.
#   make annotate ANNOTATE="citrus_sinensis carica_papaya" THREADS=16
ANNOTATE    ?= citrus_sinensis
EGGNOG_DATA ?= references/plant_host/eggnog_data
EMAPPER     ?= $(HOME)/.venvs/emapper/bin
ATH         := references/plant_host/arabidopsis_thaliana/arabidopsis_thaliana
ATH_PEP     := https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
DIAMOND_BEST = diamond blastp --more-sensitive -k 1 -e 1e-5 --threads $(THREADS) --quiet \
               --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore

annotate: $(foreach h,$(ANNOTATE),references/plant_host/$(h)/$(h).genes.tsv)

# Fetched directly rather than with emapper's download_eggnog_data.py: that script
# points at eggnogdb.embl.de, which no longer resolves (2026-09-14), and it reports
# "Finished" after wget fails. Same files, from the host that still serves them.
EGGNOG_URL := http://eggnog5.embl.de/download/emapperdb-5.0.2
$(EGGNOG_DATA)/eggnog.db:
	mkdir -p $(EGGNOG_DATA)
	cd $(EGGNOG_DATA) && for f in eggnog_proteins.dmnd.gz eggnog.taxa.tar.gz eggnog.db.gz; do \
	  curl -fsSL --retry 3 -o $$f.part $(EGGNOG_URL)/$$f && mv $$f.part $$f || exit 1; done
	cd $(EGGNOG_DATA) && gunzip -f eggnog_proteins.dmnd.gz && tar -xzf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz
	cd $(EGGNOG_DATA) && gunzip -f eggnog.db.gz
	test -s $(EGGNOG_DATA)/eggnog_proteins.dmnd && test -s $@

references/plant_host/%.protein.faa: references/plant_host/%.fasta
	curl -fsSL --retry 3 "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$(ACC_$(notdir $*))/download?include_annotation_type=PROT_FASTA" \
	  -o $(@D)/prot.zip
	unzip -oq $(@D)/prot.zip -d $(@D)/prot_tmp
	find $(@D)/prot_tmp -name "protein.faa" -exec mv {} $@ \;
	rm -rf $(@D)/prot_tmp $(@D)/prot.zip
	test -s $@ && touch $@    # unzip keeps NCBI's timestamp, which can sit in the future

references/plant_host/%.longest.faa: references/plant_host/%.protein.faa scripts/longest_proteins.py
	python3 scripts/longest_proteins.py refseq references/plant_host/$*.gff3 $< $@

# Separate target so an existing copy is reused: on 2026-09-15 a compute node could not
# resolve ftp.ensemblgenomes.ebi.ac.uk. If a job fails here, run
#   make references/plant_host/arabidopsis_thaliana/arabidopsis_thaliana.pep.all.fa.gz
# on the login node first (9.7 MB).
$(ATH).pep.all.fa.gz:
	mkdir -p $(@D)
	curl -fsSL --retry 3 -o $@.part $(ATH_PEP) && mv $@.part $@

$(ATH).longest.faa: $(ATH).pep.all.fa.gz scripts/longest_proteins.py
	python3 scripts/longest_proteins.py ensembl $< $@ $(ATH).agi.tsv

%.dmnd: %.longest.faa
	module load $(MOD_DIAMOND) && diamond makedb --in $< -d $* --threads $(THREADS) --quiet

# best hit each way; the reverse search marks reciprocal best hits (ath_rbh)
references/plant_host/%.ath.tsv: references/plant_host/%.longest.faa $(ATH).dmnd
	module load $(MOD_DIAMOND) && $(DIAMOND_BEST) -q $< -d $(ATH).dmnd -o $@

references/plant_host/%.ath_rev.tsv: $(ATH).longest.faa references/plant_host/%.dmnd
	module load $(MOD_DIAMOND) && $(DIAMOND_BEST) -q $< -d references/plant_host/$*.dmnd -o $@

references/plant_host/%.emapper.annotations: references/plant_host/%.longest.faa $(EGGNOG_DATA)/eggnog.db
	module load $(MOD_DIAMOND) && $(EMAPPER)/emapper.py -i $< --itype proteins -m diamond \
	  --tax_scope Viridiplantae --data_dir $(EGGNOG_DATA) --cpu $(THREADS) --override \
	  -o $(notdir $*) --output_dir $(@D) --temp_dir $(@D) > $(@D)/$(notdir $*).emapper.log 2>&1

references/plant_host/%.genes.tsv: references/plant_host/%.emapper.annotations \
                                   references/plant_host/%.ath.tsv references/plant_host/%.ath_rev.tsv \
                                   scripts/gene_annotation.py
	python3 scripts/gene_annotation.py references/plant_host/$*.gff3 $(wordlist 1,3,$^) $(ATH).agi.tsv $@
