# Setting up on a new cluster (e.g. Ceres)

1. Place this directory somewhere with room for ~50 GB of outputs, and put (or
   symlink) the raw reads at `30-1348328766/00_fastq/` (fastq.gz + .md5 files).
2. Create the tool environment:
       conda env create -f environment.yml && conda activate rnaseq
3. Check module names: `module avail hisat2 samtools subread fastqc`. If they
   differ from the Atlas defaults in the Makefile header, override, e.g.:
       make align MOD_HISAT2=hisat2 MOD_SAMTOOLS=samtools/1.17
4. Build everything up to alignment stats:
       make verify qc              # md5 + FastQC/MultiQC (fast)
       make refs                   # strain .gb conversion + host genome downloads (~3.6 GB)
       srun --cpus-per-task=16 --mem=64G make -j2 fractions   # indexes + alignment
   `make help` lists all targets; PIPELINE.md explains each step.

5. For heavy work, run Claude and the pipeline on a compute node rather than
   the login node - see SLURM.md and `scripts/claude-node.sh`.

6. Visualise what came out:
       scripts/viz-env.sh                   # once, builds ~/.venvs/viz
       ~/.venvs/viz/bin/jupyter lab         # notebooks/01_results_overview.ipynb
   `notebooks/01_results_overview.html` is the same thing pre-rendered.

Ceres layout (2026-09-08)
- Repo lives in `~/gallrna` (home quota is 30 GB); everything large is on
  `/90daydata/small_grains/andrew.dickson/gallrna_data/` via symlinks:
  `01_trim 02_align 03_counts 04_matrix qc references/combined references/plant_host/<genome>`.
  Raw reads: `30-1348328766 -> /90daydata/small_grains/30-1348328766` (group dir, 60 GB).
  90daydata is purged after 90 days without access - copy anything you want to keep to /project.
- Ceres module names differ from Atlas; `local.mk` (gitignored, auto-included) sets
  `MOD_SAMTOOLS = samtools/1.17` and `MOD_SUBREAD = subread/2.0.4`.
- `HOSTSEL="poncirus_trifoliata citrus_sinensis"` restricts any target to those hosts' samples.
- Until the conda env exists, `make strain-refs CONDA=$(dirname $(which python3))` works
  after `module load python_3` (it has Biopython).

6. Visualise what came out:
       scripts/viz-env.sh                   # once, builds ~/.venvs/viz
       ~/.venvs/viz/bin/jupyter lab         # notebooks/01_results_overview.ipynb
   `notebooks/01_results_overview.html` is the same thing pre-rendered.

Notes
- `references/agrobacterium/` and `references/plant_host/<genome>/` are built by
  `make refs` from the Strain */C58 folders and NCBI links in the Makefile header.
- STRAND=0 (unstranded counting) until library strandedness is confirmed from
  the first alignments; then rerun counts with `make counts STRAND=2 -B`.
- The .gitignore keeps raw data, outputs, and unpublished collaborator files
  (Strain */, C58_ALIGNED/, docs/) out of git - only code, docs, and metadata
  are committed. Adjust if/when the strain genomes are published.
- Not bundled (still on Atlas under rnaseq/Strain 29/): the collaborator's
  Snippy variant-calling outputs and Bakta annotations (~5.9 GB) - prior
  comparative analysis, not needed by this pipeline.
