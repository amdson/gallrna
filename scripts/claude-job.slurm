#!/bin/bash
# Claude Code as a batch job. Two modes:
#
#   headless (default) - run one prompt unattended, log the transcript:
#       sbatch scripts/claude-job.slurm "run make fractions and summarise the mapping rates"
#       sbatch scripts/claude-job.slurm --prompt-file notes/task.md
#
#   hold (CLAUDE_MODE=hold) - park an interactive Claude in tmux on the node
#   and keep the allocation alive so you can attach later. This is what
#   `claude-node.sh --detach` submits; you rarely invoke it by hand.
#
# Override resources at submit time, e.g.
#       sbatch -c 32 --mem=128G -t 24:00:00 scripts/claude-job.slurm "..."
#
#SBATCH --job-name=claude
#SBATCH --account=small_grains
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/claude-%j.out
#SBATCH --error=logs/claude-%j.out
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$SLURM_SUBMIT_DIR}"
cd "$PROJECT_DIR"
# shellcheck source=claude-env.sh
. scripts/claude-env.sh

echo "=== claude job ${SLURM_JOB_ID:-?} on $(hostname -s), $(date) ==="
echo "    ${SLURM_CPUS_PER_TASK:-?} cpu, ${SLURM_MEM_PER_NODE:-?}M, dir $PROJECT_DIR"
claude_net_check || exit 1

# ---------------------------------------------------------------- hold mode --
if [ "${CLAUDE_MODE:-headless}" = hold ]; then
    tmux new-session -d -s claude -c "$PROJECT_DIR" \
        "claude ${*:-}; echo '--- claude exited; shell follows ---'; exec bash -l"
    echo "tmux session 'claude' started; attach with:"
    echo "    srun --jobid=${SLURM_JOB_ID} --overlap --pty tmux attach -t claude"
    # Hold the allocation until the tmux server goes away or we hit the walltime.
    while tmux has-session -t claude 2>/dev/null; do sleep 30; done
    echo "tmux session ended, releasing allocation"
    exit 0
fi

# ------------------------------------------------------------ headless mode --
PROMPT=""
if [ "${1:-}" = "--prompt-file" ]; then
    PROMPT=$(cat "$2"); shift 2
elif [ $# -gt 0 ]; then
    PROMPT="$*"
else
    echo "no prompt given (pass one as an argument, or --prompt-file FILE)" >&2
    exit 2
fi

# acceptEdits lets Claude write files and answer its own edit prompts, but it
# still stops at anything it considers risky - which in a batch job means the
# run just ends there. Swap in --dangerously-skip-permissions only if you have
# read the prompt and are willing to let it run any command unattended.
PERM="${CLAUDE_PERMISSION_MODE:---permission-mode acceptEdits}"

echo "--- prompt ---"
echo "$PROMPT"
echo "--- output ---"
# shellcheck disable=SC2086
exec claude -p "$PROMPT" $PERM --output-format text
