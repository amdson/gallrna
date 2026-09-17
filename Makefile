# ============================================================================
# Gall dual RNA-seq pipeline (see README.md section 5 for rationale)
# Run from this directory. Heavy targets belong on a compute node, e.g.:
#     srun --cpus-per-task=16 --mem=64G make -j2 align
# `make help` lists targets. THREADS and STRAND can be overridden: make align THREADS=16
# SAMPLES="29wtHC83 1416wtCC547" (or HOSTSEL="citrus_sinensis carica_papaya") restricts a run.
# Every list below is written out by hand, on purpose: read it top to bottom, edit it in place.
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
STRAND  ?= 0            # featureCounts -s: libraries are unstranded (confirmed 2026-09-14, README.md section 3)
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

# --- the 14 samples (samples.tsv), by host --------------------------------
# <strain><genotype><host code><gall age in days>. Order = samples.tsv = matrix column order.
ALL_SAMPLES := 1416wtEu635 29wtEu635 29wtGeu182 1416wtGeu182 29wtHC83 1416wtM26 29wtM26 29wtP46 \
               1416wtCC547 1416wtCC54 1416G-19CC547 1416G-30CC547 1416wtT49 29wtT29

SAMPLES_euonymus_japonicus_proxy := 1416wtEu635 29wtEu635 29wtGeu182 1416wtGeu182
SAMPLES_citrus_sinensis          := 29wtHC83 1416wtCC547 1416wtCC54 1416G-19CC547 1416G-30CC547
SAMPLES_brassica_juncea          := 1416wtM26 29wtM26
SAMPLES_carica_papaya            := 29wtP46
SAMPLES_solanum_lycopersicum     := 1416wtT49 29wtT29
# poncirus_trifoliata is downloaded for the Carrizo plan (README.md 8.2 A); nothing maps to it yet
SAMPLES_poncirus_trifoliata      :=

# restrict a run: make align SAMPLES="29wtHC83 1416wtCC547", or by host:
#   make fractions HOSTSEL="citrus_sinensis carica_papaya"
HOSTSEL ?=
ifdef HOSTSEL
SAMPLES := $(foreach h,$(HOSTSEL),$(SAMPLES_$(h)))
else
SAMPLES ?= $(ALL_SAMPLES)
endif

# --- the combined reference each sample is aligned to: <strain>__<host> ----
# CC = Carrizo citrange, an F1 hybrid of C. sinensis x P. trifoliata (confirmed
# 2026-09-14; samples.tsv updated to match). No Carrizo genome exists, so
# the CC galls are aligned to one parent, C. sinensis, and counted on its gene
# models, as in published Carrizo RNA-seq, e.g. Afzal Naveed, Huguet-Tapia & Ali
# 2019, J Plant Interact 14:187-204, doi:10.1080/17429145.2019.1609106 (roots
# mapped to two C. sinensis genomes; they report only 55-73% transcriptome
# coverage and fell back to de novo assembly). Here the CC libraries align at
# 82-85% vs 92% for the native sweet-orange gall 29wtHC83; the measured cost is
# about a fifth of the Poncirus-copy reads, a tenth of all reads (README.md 6.2).
# All four CC galls share the hybrid genotype, so the bias from P. trifoliata
# alleles mapping with mismatches is common to every CC sample and cancels in
# CC-vs-CC contrasts; it does not cancel against galls on other hosts.
# If parent of origin or Poncirus-only genes matter, add P. trifoliata as a
# second haplotype instead: the ZK8 assembly under poncirus_trifoliata IS
# annotated at the Citrus Genome Database (25,680 genes;
# citrusgenomedb.org/jb2/data/Ptri_ZK8_v1.sorted.gff.gz), but its contigs are
# named chr1_ZK8.., so pair it with CGD's jb2/data/Ptri_ZK8_v1.fasta.gz, not
# NCBI's GCA_018350135.1 (CM031384.1..) used now.
# README.md (sections 6.2, 7.1, 8.2 A) plans that two-parent `carrizo` host:
# DVS_A1.0 + ZK8 concatenated, counted per parental copy, summed per gene pair.
COMBO_1416wtEu635   := strain_1416__euonymus_japonicus_proxy
COMBO_29wtEu635     := strain_29__euonymus_japonicus_proxy
COMBO_29wtGeu182    := strain_29__euonymus_japonicus_proxy
COMBO_1416wtGeu182  := strain_1416__euonymus_japonicus_proxy
COMBO_29wtHC83      := strain_29__citrus_sinensis
COMBO_1416wtM26     := strain_1416__brassica_juncea
COMBO_29wtM26       := strain_29__brassica_juncea
COMBO_29wtP46       := strain_29__carica_papaya
COMBO_1416wtCC547   := strain_1416__citrus_sinensis
COMBO_1416wtCC54    := strain_1416__citrus_sinensis
COMBO_1416G-19CC547 := strain_1416__citrus_sinensis
COMBO_1416G-30CC547 := strain_1416__citrus_sinensis
COMBO_1416wtT49     := strain_1416__solanum_lycopersicum
COMBO_29wtT29       := strain_29__solanum_lycopersicum

# the nine strain x host pairs actually sampled (one HISAT2 index each, step 3c)
COMBOS := strain_1416__euonymus_japonicus_proxy strain_29__euonymus_japonicus_proxy \
          strain_29__citrus_sinensis           strain_1416__citrus_sinensis \
          strain_1416__brassica_juncea         strain_29__brassica_juncea \
          strain_29__carica_papaya \
          strain_1416__solanum_lycopersicum    strain_29__solanum_lycopersicum

# --- host genomes (NCBI accessions in the header) ---------------------------
HOSTS       := euonymus_japonicus_proxy citrus_sinensis brassica_juncea carica_papaya poncirus_trifoliata solanum_lycopersicum
# ANNOT_HOSTS = the hosts with NCBI gene models, hence a plant count matrix and `make annotate`
ANNOT_HOSTS := citrus_sinensis carica_papaya solanum_lycopersicum

ACC_euonymus_japonicus_proxy := GCA_963580455.1
ACC_citrus_sinensis          := GCF_022201045.2
ACC_brassica_juncea          := GCA_018703725.1
ACC_carica_papaya            := GCF_000150535.2
ACC_poncirus_trifoliata      := GCA_018350135.1
ACC_solanum_lycopersicum     := GCF_036512215.1

# ---------------------------- phony targets ---------------------------------
.PHONY: help all verify qc trim strain-refs host-genomes refs combined align fractions counts matrices mapstats annotate

help:
	@echo "targets: verify qc trim refs (strain-refs host-genomes) combined align fractions counts matrices mapstats all annotate baselines shortlist"

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
# One line per combo names the strain FASTA and the host FASTA that go into it
# (each GFF3 sits next to its FASTA); the pattern rule below is the shared recipe.
S1416 := references/agrobacterium/strain_1416/strain_1416.fasta
S29   := references/agrobacterium/strain_29/strain_29.fasta
references/combined/strain_1416__euonymus_japonicus_proxy.fasta: $(S1416) references/plant_host/euonymus_japonicus_proxy/euonymus_japonicus_proxy.fasta
references/combined/strain_29__euonymus_japonicus_proxy.fasta:   $(S29)   references/plant_host/euonymus_japonicus_proxy/euonymus_japonicus_proxy.fasta
references/combined/strain_29__citrus_sinensis.fasta:            $(S29)   references/plant_host/citrus_sinensis/citrus_sinensis.fasta
references/combined/strain_1416__citrus_sinensis.fasta:          $(S1416) references/plant_host/citrus_sinensis/citrus_sinensis.fasta
references/combined/strain_1416__brassica_juncea.fasta:          $(S1416) references/plant_host/brassica_juncea/brassica_juncea.fasta
references/combined/strain_29__brassica_juncea.fasta:            $(S29)   references/plant_host/brassica_juncea/brassica_juncea.fasta
references/combined/strain_29__carica_papaya.fasta:              $(S29)   references/plant_host/carica_papaya/carica_papaya.fasta
references/combined/strain_1416__solanum_lycopersicum.fasta:     $(S1416) references/plant_host/solanum_lycopersicum/solanum_lycopersicum.fasta
references/combined/strain_29__solanum_lycopersicum.fasta:       $(S29)   references/plant_host/solanum_lycopersicum/solanum_lycopersicum.fasta
# $< is the strain FASTA, $(word 2,$^) the host FASTA
references/combined/%.fasta:
	mkdir -p references/combined
	cat $^ > $@
	cat $(<:.fasta=.gff3) > references/combined/$*.gff3
	@hostgff=$(word 2,$(^:.fasta=.gff3)); \
	  if [ -f "$$hostgff" ]; then grep -v "^#" "$$hostgff" >> references/combined/$*.gff3; \
	  else echo "NOTE: no host annotation for $* - combined GFF is agro-only"; fi

references/combined/%.hisat2.ok: references/combined/%.fasta
	module load $(MOD_HISAT2) && hisat2-build -p $(THREADS) $< references/combined/$* > references/combined/$*.build.log 2>&1
	touch $@

# ---------------------------- 4. align (competitive) ------------------------
02_align/%.bam: 01_trim/%_R1.fastq.gz 01_trim/%_R2.fastq.gz references/combined/$$(COMBO_$$*).hisat2.ok
	mkdir -p 02_align
	module load $(MOD_HISAT2) $(MOD_SAMTOOLS) && \
	hisat2 -p $(THREADS) --dta -x references/combined/$(COMBO_$*) \
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
	  -a references/combined/$(COMBO_$*).gff3 -o $@ $< 2> 03_counts/$*.log

# One matrix per strain (its bacterial genes, all its samples) and one per annotated host
# (plant genes, all samples on that host); columns in this order. A matrix always takes
# every sample, whatever SAMPLES says: a partial matrix would silently replace a full one.
04_matrix/agro_strain_1416.tsv: 03_counts/1416wtEu635.txt 03_counts/1416wtGeu182.txt 03_counts/1416wtM26.txt \
                                03_counts/1416wtCC547.txt 03_counts/1416wtCC54.txt 03_counts/1416G-19CC547.txt \
                                03_counts/1416G-30CC547.txt 03_counts/1416wtT49.txt
	mkdir -p 04_matrix
	$(CONDA)/python scripts/merge_counts.py $@ keep agro_ $^

04_matrix/agro_strain_29.tsv: 03_counts/29wtEu635.txt 03_counts/29wtGeu182.txt 03_counts/29wtHC83.txt \
                              03_counts/29wtM26.txt 03_counts/29wtP46.txt 03_counts/29wtT29.txt
	mkdir -p 04_matrix
	$(CONDA)/python scripts/merge_counts.py $@ keep agro_ $^

04_matrix/plant_citrus_sinensis.tsv: 03_counts/29wtHC83.txt 03_counts/1416wtCC547.txt 03_counts/1416wtCC54.txt \
                                     03_counts/1416G-19CC547.txt 03_counts/1416G-30CC547.txt
	mkdir -p 04_matrix
	$(CONDA)/python scripts/merge_counts.py $@ drop agro_ $^

04_matrix/plant_carica_papaya.tsv: 03_counts/29wtP46.txt
	mkdir -p 04_matrix
	$(CONDA)/python scripts/merge_counts.py $@ drop agro_ $^

04_matrix/plant_solanum_lycopersicum.tsv: 03_counts/1416wtT49.txt 03_counts/29wtT29.txt
	mkdir -p 04_matrix
	$(CONDA)/python scripts/merge_counts.py $@ drop agro_ $^

# ---------------------------- 5b. mapping statistics ------------------------
# Per-sample mapping QC for the paper's supplementary table: read retention, alignment
# (unique/multi), host vs bacterial reads (primary alignments, as a share of all reads
# sequenced), host mismatch rate and featureCounts assignment. It records that the
# single-copy host references gave reasonable results even though the hosts are
# heterozygous, hybrid or polyploid, and it bounds plant->bacterium cross-mapping.
# Column guide and how to read it: README.md step 5b. samtools stats reads every host
# alignment, so run it on a compute node, e.g.
#   srun --cpus-per-task=16 --mem=16G make -j4 mapstats THREADS=4
02_align/%.mapstats.tsv: 02_align/%.bam 03_counts/%.txt scripts/mapstats.py
	module load $(MOD_SAMTOOLS) && python3 scripts/mapstats.py $* $(COMBO_$*) $< \
	  02_align/$*.hisat2.log 01_trim/$*.fastp.json 03_counts/$*.txt $(THREADS) > $@.tmp
	mv $@.tmp $@

logs/mapping_stats.tsv: $(foreach s,$(SAMPLES),02_align/$(s).mapstats.tsv)
	mkdir -p logs
	awk 'FNR > 1 || NR == 1' $^ > $@
	@echo "$@: $$(($$(wc -l < $@) - 1)) samples"

# ---------------------------- 5c. public baselines --------------------------
# Healthy-tissue RNA-seq from ENA (run list: citrus_baselines.tsv, chosen in
# README.md section 7.3), put through the same fastp / HISAT2 / featureCounts
# settings as the galls, but against the plant-only reference: there are no bacterial
# reads to compete with. Carrizo runs go to the sweet-orange reference, exactly like the
# CC galls, until the two-parent `carrizo` host exists (README.md section 8.2 A). The Poncirus
# runs are listed but in no tier's run list until ZK8 has an annotated reference.
# Output: 04_matrix/baseline_citrus_sinensis.tsv, one column per run accession in the order
# below (tissue and replicate labels are in citrus_baselines.tsv). Tier 1 = 27 runs, ~250 GB
# under 05_baseline/ (symlink it into 90daydata like the other output dirs).
# Run it as a batch job: sbatch scripts/baselines.slurm
#   make baselines BASE_TIER=2                          # add the tier-2 Carrizo leaf controls
#   make baselines BASE_RUNS="SRR32263188 SRR32263215"  # a subset
BASE_DIR := 05_baseline
# every run below is aligned to sweet orange, the Carrizo ones exactly like the CC galls
BASE_REF := citrus_sinensis

# tier 1: sweet orange, PRJNA599503 (healthy bark CK_P1-3, root CK_R1-3, young leaf CK_L1-3)
BASE_SO_BARK   := SRR10848795 SRR10848794 SRR10848793
BASE_SO_ROOT   := SRR10848792 SRR10848791 SRR10848790
BASE_SO_LEAF   := SRR10848806 SRR10848805 SRR10848796
# tier 1: Valencia embryogenic callus, empty-vector line EV-L44 at 0, 15 and 30 d, PRJNA778304
BASE_SO_CALLUS := SRR16874257 SRR16874256 SRR16874247 SRR16874246 SRR16874245 SRR16874244 SRR16874243 SRR16874242 SRR16874241
# tier 1: Carrizo citrange wt, PRJNA1216034 (18 d stems, young leaves, untreated mature leaves)
BASE_CC_STEM   := SRR32263188 SRR32263215 SRR32263214
BASE_CC_YLEAF  := SRR32263194 SRR32263193 SRR32263192
BASE_CC_MLEAF  := SRR32263197 SRR32263196 SRR32263195
# tier 2: Carrizo leaf controls, PRJNA839431 (CT-Car), PRJNA668159 (Carrizo Control), PRJNA1053671 (ZC-1..3)
BASE_CC_LEAF2  := SRR19277106 SRR19277105 SRR19277104 SRR12799421 SRR12799422 SRR12799423 SRR27236378 SRR27236377 SRR27236376
# tier 1 but not run: Poncirus ZK stem and thorn, PRJDB41296 - needs an annotated ZK8 reference (README.md 8.2 A)
BASE_PT_STEM   := DRR1031157 DRR1031158 DRR1031159
BASE_PT_THORN  := DRR1031160 DRR1031161 DRR1031162

BASE_TIER1 := $(BASE_SO_BARK) $(BASE_SO_ROOT) $(BASE_SO_LEAF) $(BASE_SO_CALLUS) $(BASE_CC_STEM) $(BASE_CC_YLEAF) $(BASE_CC_MLEAF)
BASE_TIER  ?= 1
ifeq ($(BASE_TIER),2)
BASE_RUNS  ?= $(BASE_TIER1) $(BASE_CC_LEAF2)
else
BASE_RUNS  ?= $(BASE_TIER1)
endif

.PHONY: baselines baseline-fetch baseline-trim baseline-align baseline-counts
baselines:       04_matrix/baseline_$(BASE_REF).tsv
baseline-fetch:  $(foreach r,$(BASE_RUNS),$(BASE_DIR)/fastq/$(r)_1.fastq.gz)
baseline-trim:   $(foreach r,$(BASE_RUNS),$(BASE_DIR)/trim/$(r)_R1.fastq.gz)
baseline-align:  $(foreach r,$(BASE_RUNS),$(BASE_DIR)/align/$(r).bam)
baseline-counts: $(foreach r,$(BASE_RUNS),$(BASE_DIR)/counts/$(r).txt)

$(BASE_DIR)/fastq/%_1.fastq.gz $(BASE_DIR)/fastq/%_2.fastq.gz:
	scripts/ena_fetch.sh $* $(BASE_DIR)/fastq

$(BASE_DIR)/trim/%_R1.fastq.gz $(BASE_DIR)/trim/%_R2.fastq.gz: $(BASE_DIR)/fastq/%_1.fastq.gz $(BASE_DIR)/fastq/%_2.fastq.gz
	mkdir -p $(BASE_DIR)/trim
	$(CONDA)/fastp -w $(THREADS) -i $(BASE_DIR)/fastq/$*_1.fastq.gz -I $(BASE_DIR)/fastq/$*_2.fastq.gz \
	  -o $(BASE_DIR)/trim/$*_R1.fastq.gz -O $(BASE_DIR)/trim/$*_R2.fastq.gz \
	  -j $(BASE_DIR)/trim/$*.fastp.json -h $(BASE_DIR)/trim/$*.fastp.html 2> $(BASE_DIR)/trim/$*.fastp.log

# plant-only HISAT2 index next to the genome: references/plant_host/<host>/<host>.*.ht2
references/plant_host/%.hisat2.ok: references/plant_host/%.fasta
	module load $(MOD_HISAT2) && hisat2-build -p $(THREADS) $< references/plant_host/$* > references/plant_host/$*.build.log 2>&1
	touch $@

$(BASE_DIR)/align/%.bam: $(BASE_DIR)/trim/%_R1.fastq.gz $(BASE_DIR)/trim/%_R2.fastq.gz \
                         references/plant_host/$(BASE_REF)/$(BASE_REF).hisat2.ok
	mkdir -p $(BASE_DIR)/align
	module load $(MOD_HISAT2) $(MOD_SAMTOOLS) && \
	hisat2 -p $(THREADS) --dta -x references/plant_host/$(BASE_REF)/$(BASE_REF) \
	  -1 $(BASE_DIR)/trim/$*_R1.fastq.gz -2 $(BASE_DIR)/trim/$*_R2.fastq.gz 2> $(BASE_DIR)/align/$*.hisat2.log \
	  | samtools sort -@ $(THREADS) -o $@ - && samtools index $@

$(BASE_DIR)/counts/%.txt: $(BASE_DIR)/align/%.bam
	mkdir -p $(BASE_DIR)/counts
	module load $(MOD_SUBREAD) && \
	featureCounts -p --countReadPairs -T $(THREADS) -s $(STRAND) -t gene -g ID \
	  -a references/plant_host/$(BASE_REF)/$(BASE_REF).gff3 -o $@ $< 2> $(BASE_DIR)/counts/$*.log

04_matrix/baseline_$(BASE_REF).tsv: $(foreach r,$(BASE_RUNS),$(BASE_DIR)/counts/$(r).txt)
	mkdir -p 04_matrix
	$(CONDA)/python scripts/merge_counts.py $@ drop agro_ $^

# ---------------------------- 5d. candidate shortlist, method 1 -------------
# README.md section 8.1, method 1: rank host genes by expression in the galls
# themselves (TPM -> within-gall percentile; a group's score is the gene's lowest percentile
# across the group's galls). Needs only the gall counts and genes.tsv, no baseline. Writes
# the full table to 04_matrix/ and one top-20 per group to results/shortlist/ (in git).
# Methods 2-6 need the baselines (5c) and are not automated yet.
SHORTLIST_GROUPS := hamlin=29wtHC83 carrizo_wt=1416wtCC547,1416wtCC54 \
                    carrizo_eng=1416G-19CC547,1416G-30CC547 \
                    all=29wtHC83,1416wtCC547,1416wtCC54,1416G-19CC547,1416G-30CC547
.PHONY: shortlist
shortlist: results/shortlist/method1_expression_all_top20.tsv
results/shortlist/method1_expression_all_top20.tsv: 04_matrix/plant_citrus_sinensis.tsv \
        references/plant_host/citrus_sinensis/citrus_sinensis.genes.tsv scripts/shortlist_expression.py
	python3 scripts/shortlist_expression.py $< 03_counts/29wtHC83.txt $(word 2,$^) \
	  04_matrix/citrus_gall_expression.tsv results/shortlist $(SHORTLIST_GROUPS)

# 6. differential expression: not automated yet - needs strandedness + design
# decisions (see README.md step 6). Matrices in 04_matrix/ are DESeq2-ready.

# ---------------------------- 7. functional annotation ----------------------
# One protein per gene (longest RefSeq isoform, named by the gene ID featureCounts
# writes) -> eggNOG-mapper (GO, KEGG, EC, Pfam, description) + DIAMOND best hit in
# Arabidopsis (AGI + symbol; the crown-gall marker lists in README.md section 7.4 are
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
