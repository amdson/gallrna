# Makefile, annotated

A reading copy of `Makefile` (2026-09-17). Not executed; when the Makefile changes, this goes stale.
Each block: the make lines, what they mean, input → output.

## Make syntax used in this file

| syntax | meaning |
|---|---|
| `target: prerequisites` + tab-indented lines | a **rule**: if `target` is missing or older than any prerequisite, run the indented shell lines (the **recipe**). Prerequisites are built first, by their own rules |
| `VAR := value` | set once, now. `VAR ?= value`: set only if not already set (so `make VAR=x` or the environment wins). `VAR = value`: re-expanded at every use |
| `$(VAR)`, `$(VAR_$(x))` | value of VAR; the second builds the name first (`$(COMBO_$*)` → `$(COMBO_29wtHC83)`) |
| `%` in a target | a **pattern rule**: `03_counts/%.txt` matches `03_counts/29wtHC83.txt`, and `%` = `29wtHC83`. A rule with `%` in two targets builds both in one run |
| `$@` `$<` `$^` `$*` | in a recipe: the target; the first prerequisite; all prerequisites; the `%` match. `$(@D)` = the target's directory |
| `$$` | a literal `$` for the shell (make eats one) |
| `@cmd` / `-cmd` | don't echo the command / ignore its failure |
| `$(foreach v,list,text)` | `text` once per list item, `v` substituted |
| `$(word 2,$^)`, `$(wordlist 1,3,$^)`, `$(notdir x)`, `$(dir x)`, `$(wildcard g)` | 2nd word; words 1-3; basename; directory; files matching a glob |
| `$(<:.fasta=.gff3)` | `$<` with the suffix swapped |
| `$$(COMBO_$$*)` in a prerequisite | with `.SECONDEXPANSION:`, `$$` is expanded a second time after `%` is known, so prerequisites can depend on the matched name |
| `.PHONY: name` | `name` is a label, not a file: always runs its prerequisites |
| `.SECONDARY:` | never delete intermediate files |
| `ifdef` / `ifeq` / `else` / `endif` | evaluated when the Makefile is read |
| `-include local.mk` | read that file if it exists |
| `module load x && cmd` | the cluster's environment modules; `&&` chains so a failed load aborts |

## Settings (lines 33-51)

```make
SHELL := /bin/bash
.SHELLFLAGS := -lc          # login shell, so `module` exists
THREADS ?= 8                # per command; make -j2 runs two commands at once, so 2 x THREADS cores
STRAND  ?= 0                # featureCounts -s 0: unstranded
CONDA   ?= $(HOME)/.conda/envs/rnaseq/bin
MOD_HISAT2 ?= hisat2/2.2.1  # ... module names; local.mk overrides them on Ceres
-include local.mk
RAW     := 30-1348328766/00_fastq
```

## Sample lists (lines 53-128)

```make
ALL_SAMPLES := 1416wtEu635 ... 29wtT29                  # the 14 samples, samples.tsv order
SAMPLES_citrus_sinensis := 29wtHC83 1416wtCC547 ...     # one list per host
HOSTSEL ?=
ifdef HOSTSEL
SAMPLES := $(foreach h,$(HOSTSEL),$(SAMPLES_$(h)))      # HOSTSEL="citrus_sinensis carica_papaya" -> those hosts' samples
else
SAMPLES ?= $(ALL_SAMPLES)                               # or SAMPLES="29wtHC83 ..." on the command line, or all 14
endif
COMBO_29wtHC83 := strain_29__citrus_sinensis            # per sample: which combined reference (step 3c) it aligns to
COMBOS := strain_1416__euonymus_japonicus_proxy ...     # the 9 strain x host pairs that need an index
HOSTS := ...                                            # the 6 host genomes to download
ANNOT_HOSTS := citrus_sinensis carica_papaya solanum_lycopersicum   # hosts with gene models
ACC_citrus_sinensis := GCF_022201045.2                  # NCBI accession per host
```

`SAMPLES` drives `trim`, `align`, `counts`, `fractions` and `mapstats`. The matrices ignore it (all samples, always).

## Phony targets (lines 130-150)

```make
all: verify qc fractions matrices      # `make all` = these four
trim:    $(foreach s,$(SAMPLES),01_trim/$(s)_R1.fastq.gz)        # one file per sample; its rule makes R2 too
align:   $(foreach s,$(SAMPLES),02_align/$(s).bam)
counts:  $(foreach s,$(SAMPLES),03_counts/$(s).txt)
combined: $(foreach c,$(COMBOS),references/combined/$(c).hisat2.ok)
matrices: 04_matrix/agro_strain_1416.tsv 04_matrix/agro_strain_29.tsv 04_matrix/plant_<host>.tsv (x3)
```

Names for groups of files. `make align` asks for 14 BAMs; the pattern rules below make each one, pulling in trimmed reads and indexes as needed.

## Step 0: verify (line 153)

```make
logs/md5.ok:
	cd $(RAW) && md5sum -c --quiet *.md5
	touch $@
```
In: raw fastq + `.md5`. Out: an empty stamp file, so it runs once.

## Step 1: raw QC (line 159)

```make
qc/multiqc_report.html: $(wildcard $(RAW)/*.fastq.gz)
	fastqc ... ; multiqc ...
```
In: all raw fastq. Out: `qc/fastqc/*` and one MultiQC report.

## Step 2: trim (line 165)

```make
01_trim/%_R1.fastq.gz 01_trim/%_R2.fastq.gz: $(RAW)/%_R1_001.fastq.gz $(RAW)/%_R2_001.fastq.gz
	fastp -i ... -I ... -o ... -O ... -j 01_trim/$*.fastp.json
```
Per sample. In: raw R1/R2. Out: trimmed R1/R2 + `fastp.json/html/log`. Two targets in one rule = one fastp run per sample.

## Step 3a: strain references (lines 174-190)

```make
references/agrobacterium/strain_1416/strain_1416.fasta: scripts/convert_refs.py
	python scripts/convert_refs.py <out prefix> "agro_1416_circ=Strain 1416/1416 Circular Chromosome FINAL.gb" ...
```
In: the five `.gb`/`.gbk` replicon files (not listed as prerequisites: their names contain spaces). Out: `strain_1416.fasta` + `.gff3`, contigs named `agro_1416_circ` etc. Same for strain 29. Delete the outputs to force a rerun.

## Step 3b: host genomes (line 193)

```make
references/plant_host/%.fasta:
	$(eval NAME := $(notdir $*))        # % = citrus_sinensis/citrus_sinensis -> NAME = citrus_sinensis
	curl <NCBI datasets API>/$(ACC_$(NAME))/download?...GENOME_FASTA...GENOME_GFF -o ....zip
	unzip; mv *.fna -> $(NAME).fasta; mv genomic.gff -> $(NAME).gff3 (ignored if absent); rm tmp
```
In: nothing local (download by accession). Out: `references/plant_host/<host>/<host>.fasta` and, for annotated hosts, `.gff3`.

## Step 3c: combined references + index (lines 208-230)

```make
S1416 := references/agrobacterium/strain_1416/strain_1416.fasta
references/combined/strain_29__citrus_sinensis.fasta: $(S29) references/plant_host/citrus_sinensis/citrus_sinensis.fasta
...                                                   # one such line per combo: prerequisites only, no recipe
references/combined/%.fasta:                          # the shared recipe for all of them
	cat $^ > $@                                       # strain FASTA + host FASTA
	cat $(<:.fasta=.gff3) > references/combined/$*.gff3           # strain GFF3
	@hostgff=$(word 2,$(^:.fasta=.gff3)); if [ -f "$$hostgff" ]; then grep -v "^#" "$$hostgff" >> ...gff3; fi   # + host GFF3 if it exists

references/combined/%.hisat2.ok: references/combined/%.fasta
	hisat2-build -p $(THREADS) $< references/combined/$*
	touch $@
```
In: one strain FASTA/GFF3 + one host FASTA/GFF3. Out: `references/combined/<strain>__<host>.fasta`, `.gff3`, the HISAT2 index files `*.ht2`, and a `.hisat2.ok` stamp (the index is 8 files; the stamp is the one make tracks).

## Step 4: align (lines 233-250)

```make
02_align/%.bam: 01_trim/%_R1.fastq.gz 01_trim/%_R2.fastq.gz references/combined/$$(COMBO_$$*).hisat2.ok
	hisat2 -p $(THREADS) --dta -x references/combined/$(COMBO_$*) -1 ... -2 ... 2> 02_align/$*.hisat2.log \
	  | samtools sort -o $@ - && samtools index $@
```
Per sample. In: trimmed reads + that sample's combined index (looked up through `COMBO_<sample>`). Out: sorted `.bam`, `.bam.bai`, `.hisat2.log` (alignment rates).

```make
02_align/%.idxstats: 02_align/%.bam
	samtools idxstats $< > $@
logs/mapping_fractions.tsv: $(foreach s,$(SAMPLES),02_align/$(s).idxstats)
	for each sample: awk sums column 3 (mapped reads) over agro_* contigs vs the rest
```
Out: per-sample reads per contig, then one table: sample, agro reads, plant reads, agro fraction. Superseded by step 5b for reporting.

## Step 5: count (lines 253-284)

```make
03_counts/%.txt: 02_align/%.bam
	featureCounts -p --countReadPairs -T $(THREADS) -s $(STRAND) -t gene -g ID -a references/combined/$(COMBO_$*).gff3 -o $@ $<
```
Per sample. In: BAM + the same combined GFF3. Out: `03_counts/<sample>.txt` (gene, contig, start, end, strand, length, count) and `.txt.summary` (assigned / unassigned reasons).

```make
04_matrix/agro_strain_1416.tsv: 03_counts/1416wtEu635.txt ... 03_counts/1416wtT49.txt      # its 8 samples, in this column order
	python scripts/merge_counts.py $@ keep agro_ $^                                          # keep rows on agro_* contigs
04_matrix/plant_citrus_sinensis.tsv: 03_counts/29wtHC83.txt ... 03_counts/1416G-30CC547.txt
	python scripts/merge_counts.py $@ drop agro_ $^                                          # drop them: plant genes
```
Five explicit rules. In: the listed count files. Out: one gene x sample matrix each.

## Step 5b: mapping statistics (lines 295-303)

```make
02_align/%.mapstats.tsv: 02_align/%.bam 03_counts/%.txt scripts/mapstats.py
	python3 scripts/mapstats.py <sample> $(COMBO_$*) <bam> <hisat2.log> <fastp.json> <counts> $(THREADS) > $@.tmp
	mv $@.tmp $@                       # write-then-rename: no half-written file if it dies
logs/mapping_stats.tsv: $(foreach s,$(SAMPLES),02_align/$(s).mapstats.tsv)
	awk 'FNR > 1 || NR == 1' $^ > $@   # concatenate, keeping only the first header
```
In: per sample, the BAM, HISAT2 log, fastp JSON and counts. Out: one row per sample, then the joined table. Depends on the script, so editing `mapstats.py` marks all rows stale.

## Step 5c: public baselines (lines 318-383)

```make
BASE_DIR := 05_baseline
BASE_REF := citrus_sinensis                       # all runs align to sweet orange
BASE_SO_BARK := SRR10848795 SRR10848794 SRR10848793   # ... one list per tissue
BASE_TIER1 := $(BASE_SO_BARK) $(BASE_SO_ROOT) ... $(BASE_CC_MLEAF)    # 27 runs
BASE_TIER ?= 1
ifeq ($(BASE_TIER),2)
BASE_RUNS ?= $(BASE_TIER1) $(BASE_CC_LEAF2)       # make baselines BASE_TIER=2 adds 9 leaf runs
else
BASE_RUNS ?= $(BASE_TIER1)                        # or BASE_RUNS="SRR... SRR..." for a subset
endif
baselines: 04_matrix/baseline_$(BASE_REF).tsv
baseline-fetch / -trim / -align / -counts: $(foreach r,$(BASE_RUNS),<that stage's file for r>)   # phase targets for the batch script
```

```make
$(BASE_DIR)/fastq/%_1.fastq.gz $(BASE_DIR)/fastq/%_2.fastq.gz:
	scripts/ena_fetch.sh $* $(BASE_DIR)/fastq                       # download run % from ENA, md5-checked
$(BASE_DIR)/trim/%_R1.fastq.gz $(BASE_DIR)/trim/%_R2.fastq.gz: <the two fastq>
	fastp ...                                                       # as step 2
references/plant_host/%.hisat2.ok: references/plant_host/%.fasta
	hisat2-build $< references/plant_host/$*                        # plant-only index, next to the genome
$(BASE_DIR)/align/%.bam: <trimmed pair> references/plant_host/$(BASE_REF)/$(BASE_REF).hisat2.ok
	hisat2 ... -x references/plant_host/citrus_sinensis/citrus_sinensis | samtools sort   # as step 4, plant-only index
$(BASE_DIR)/counts/%.txt: $(BASE_DIR)/align/%.bam
	featureCounts ... -a references/plant_host/citrus_sinensis/citrus_sinensis.gff3    # as step 5, host GFF3 only
04_matrix/baseline_citrus_sinensis.tsv: $(foreach r,$(BASE_RUNS),$(BASE_DIR)/counts/$(r).txt)
	python scripts/merge_counts.py $@ drop agro_ $^
```
Per run: ENA fastq → trimmed → BAM → counts, same settings as the galls. Out: `04_matrix/baseline_citrus_sinensis.tsv`, one column per run in `BASE_RUNS` order. Same gene IDs as `plant_citrus_sinensis.tsv`.

## Step 5d: shortlist, method 1 (lines 391-399)

```make
SHORTLIST_GROUPS := hamlin=29wtHC83 carrizo_wt=1416wtCC547,1416wtCC54 carrizo_eng=... all=...
shortlist: results/shortlist/method1_expression_all_top20.tsv
results/shortlist/method1_expression_all_top20.tsv: 04_matrix/plant_citrus_sinensis.tsv \
        references/plant_host/citrus_sinensis/citrus_sinensis.genes.tsv scripts/shortlist_expression.py
	python3 scripts/shortlist_expression.py <matrix> 03_counts/29wtHC83.txt <genes.tsv> \
	  04_matrix/citrus_gall_expression.tsv results/shortlist $(SHORTLIST_GROUPS)
```
In: citrus gall matrix, one counts file (for gene lengths), the annotation table. Out: `04_matrix/citrus_gall_expression.tsv` (all genes) and `results/shortlist/method1_expression_<group>_top20.tsv` x 4. Make tracks only the `all` file.

## Step 7: functional annotation (lines 415-477)

```make
ANNOTATE    ?= citrus_sinensis                     # make annotate ANNOTATE="citrus_sinensis carica_papaya"
EGGNOG_DATA ?= references/plant_host/eggnog_data   # set to the 90daydata path in local.mk
ATH := references/plant_host/arabidopsis_thaliana/arabidopsis_thaliana     # path prefix for Arabidopsis files
DIAMOND_BEST = diamond blastp --more-sensitive -k 1 -e 1e-5 ...            # one best hit per query
annotate: $(foreach h,$(ANNOTATE),references/plant_host/$(h)/$(h).genes.tsv)
```

The chain for host `H` (`%` = `H/H`, e.g. `citrus_sinensis/citrus_sinensis`):

| rule | in → out |
|---|---|
| `$(EGGNOG_DATA)/eggnog.db:` | curl three archives from eggnog5.embl.de → `eggnog.db`, `eggnog_proteins.dmnd`, taxonomy (~50 GB, once) |
| `%.protein.faa: %.fasta` | NCBI API `PROT_FASTA` → all RefSeq proteins of H (`touch` fixes a future timestamp) |
| `%.longest.faa: %.protein.faa` | `longest_proteins.py refseq <gff3>` → one protein per gene, named by gene ID |
| `$(ATH).pep.all.fa.gz:` | curl Ensembl Plants → Arabidopsis proteome |
| `$(ATH).longest.faa: $(ATH).pep.all.fa.gz` | `longest_proteins.py ensembl` → one per gene + `agi.tsv` (AGI, symbol, description) |
| `%.dmnd: %.longest.faa` | `diamond makedb` → search database (H and Arabidopsis) |
| `%.ath.tsv: %.longest.faa $(ATH).dmnd` | H proteins searched against Arabidopsis → best hit per H gene |
| `%.ath_rev.tsv: $(ATH).longest.faa %.dmnd` | the reverse → marks reciprocal best hits |
| `%.emapper.annotations: %.longest.faa eggnog.db` | `emapper.py --tax_scope Viridiplantae` → GO, KEGG, EC, Pfam, description (hours) |
| `%.genes.tsv: the three above + gene_annotation.py` | `gene_annotation.py <gff3> <emapper> <ath> <ath_rev> <agi.tsv>` → one row per counted gene |

Out: `references/plant_host/H/H.genes.tsv`, joinable to `04_matrix/plant_H.tsv` on gene ID.
