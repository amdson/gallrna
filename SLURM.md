# Running Claude on compute nodes

The login node is for editing and submitting; anything that actually burns CPU
(alignment, counting, a Claude session that drives those) belongs on a compute
node. Scripts in `scripts/` cover the usual cases.

Verified on Ceres 2026-09-08: the partition is `ceres` (there is no `short`),
the account is `small_grains`, nodes are 96-core / 752G, default walltime is 2h
and the cap is 60 days. Compute nodes accept `ssh` from the login node and can
reach `api.anthropic.com` directly, so no proxy is needed.

**Submit from a login node, not the DTN.** `ceres-dtn.scinet.usda.gov` has no
Slurm client at all - `ssh ceres.scinet.usda.gov` first. `$HOME` is shared, so
Claude's auth and settings follow you to whichever node you land on; nothing to
copy per job.

## 0. Bare template: hold a node and ssh in

The simplest thing that works - `scripts/claude-session.slurm` parks a node and
tells you how to get on it.

```bash
mkdir -p logs
sbatch scripts/claude-session.slurm
cat logs/claude-<jobid>.out
```

```
node      : ceres20-compute-75
Get a shell on it:
    ssh ceres20-compute-75
Then:
    cd /home/andrew.dickson/gallrna && claude
Release the node when done:  scancel <jobid>
```

Nothing clever: it prints that, then `sleep infinity` until the walltime ends or
you `scancel`. Edit the `#SBATCH` lines at the top to change cpu/mem/time.

## 1. Interactive session (the common case)

```bash
tmux                              # so a dropped ssh doesn't kill the job
cd ~/gallrna
scripts/claude-node.sh            # 16 cpu, 64G, 8h
```

You get Claude running on the node, in this directory. `make -j2 fractions`
and friends then run inside the allocation with no `srun` prefix needed.

```bash
scripts/claude-node.sh -c 32 -m 128G -t 24:00:00
scripts/claude-node.sh --shell          # plain shell instead of Claude
scripts/claude-node.sh -- --continue    # args after -- go to claude
```

## 2. Detached session (log off, come back later)

```bash
scripts/claude-node.sh --detach -c 32 -m 128G -t 48:00:00
# -> Submitted job 1234567 ... attach with:
srun --jobid=1234567 --overlap --pty tmux attach -t claude
```

Claude runs in tmux on the compute node; ctrl-b d detaches and leaves it
working. `scancel 1234567` ends it. Transcript lands in `logs/claude-1234567.out`.

## 3. Unattended pipeline run (no Claude)

```bash
sbatch scripts/pipeline.slurm                    # all samples, 48 cpu / 128G / 12h
sbatch scripts/pipeline.slurm citrus_sinensis    # HOSTSEL subset
tail -f logs/pipeline-<jobid>.out
```

Runs `make` phase by phase (refs, indexes, trim, fractions, counts, matrices, qc)
with `-k`, so one bad sample does not stop the rest, and prints the mapping
fractions and matrix sizes at the end. Uses `~/.local/tools` wrappers for
fastp/python/multiqc so it works even while the conda env is still installing.
First run 2026-09-08: job 21995291.

## 4. Unattended headless Claude run

```bash
mkdir -p logs
sbatch scripts/claude-job.slurm "run make fractions, then summarise the mapping rates in logs/mapping_fractions.tsv"
sbatch -c 32 --mem=128G -t 24:00:00 scripts/claude-job.slurm --prompt-file task.md
```

One prompt, no interaction, output to `logs/claude-<jobid>.out`. It runs with
`--permission-mode acceptEdits`, so Claude will edit files on its own but stops
at anything riskier - in a batch job that means the run simply ends there. Set
`CLAUDE_PERMISSION_MODE=--dangerously-skip-permissions` only for a prompt you
have read and are happy to let run any command unattended.

## Notes

- **Submit from a login node.** `ceres-dtn` (the data transfer node) has no
  Slurm client, so none of this works there.
- **Account.** `small_grains` is your default and is baked into the scripts;
  pass `-A <account>` or export `SLURM_ACCOUNT` to charge a different one.
- **Walltime.** The partition default is only 2h - always pass `-t`/`--time`
  for a long session. The `long` QOS is available if you need past 21 days.
- **Outbound HTTPS.** Confirmed working from compute nodes, but the scripts
  still check before starting Claude and print proxy instructions if a node
  ever turns out to be walled off (see `scripts/claude-env.sh`).
