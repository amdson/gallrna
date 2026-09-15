#!/bin/bash
# Run a Claude Code session on a Ceres compute node instead of the login node.
#
#   scripts/claude-node.sh                    # interactive session, 16 cpu / 64G / 8h
#   scripts/claude-node.sh -c 32 -m 128G -t 24:00:00
#   scripts/claude-node.sh --shell            # just a shell on the node
#   scripts/claude-node.sh --detach           # session survives losing your ssh
#   scripts/claude-node.sh -- --continue      # extra args go to claude
#
# Run this from a login node (ceres.scinet.usda.gov); the DTN has no Slurm.
set -euo pipefail

CPUS=16
MEM=64G
TIME=08:00:00
PARTITION="${CLAUDE_PARTITION:-ceres}"
ACCOUNT="${SLURM_ACCOUNT:-small_grains}"
DETACH=0
SHELL_ONLY=0
FORCE=0

usage() { awk 'NR>1 { if (!/^#/) exit; sub(/^# ?/, ""); print }' "$0"; exit "${1:-0}"; }

while [ $# -gt 0 ]; do
    case "$1" in
        -c|--cpus)      CPUS=$2; shift 2 ;;
        -m|--mem)       MEM=$2; shift 2 ;;
        -t|--time)      TIME=$2; shift 2 ;;
        -p|--partition) PARTITION=$2; shift 2 ;;
        -A|--account)   ACCOUNT=$2; shift 2 ;;
        --detach)       DETACH=1; shift ;;
        --shell)        SHELL_ONLY=1; shift ;;
        -f|--force)     FORCE=1; shift ;;
        -h|--help)      usage 0 ;;
        --)             shift; break ;;
        *)              echo "unknown option: $1" >&2; usage 1 ;;
    esac
done
CLAUDE_ARGS=("$@")

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=claude-env.sh
. "$HERE/claude-env.sh"

if ! command -v salloc >/dev/null 2>&1; then
    cat >&2 <<MSG
No Slurm client on $(hostname -s). You are probably on the data transfer node.
Submit from a login node instead:
    ssh $USER@ceres.scinet.usda.gov
    cd $PROJECT_DIR && scripts/claude-node.sh
MSG
    exit 1
fi

SLURM_OPTS=(--job-name=claude --partition="$PARTITION"
            --nodes=1 --ntasks=1 --cpus-per-task="$CPUS"
            --mem="$MEM" --time="$TIME")
[ -n "$ACCOUNT" ] && SLURM_OPTS+=(--account="$ACCOUNT")

# ---------------------------------------------------------------- detached --
if [ "$DETACH" = 1 ]; then
    mkdir -p "$PROJECT_DIR/logs"
    jobid=$(sbatch --parsable "${SLURM_OPTS[@]}" \
                --output="$PROJECT_DIR/logs/claude-%j.out" \
                --export=ALL,CLAUDE_MODE=hold,PROJECT_DIR="$PROJECT_DIR" \
                "$HERE/claude-job.slurm" "${CLAUDE_ARGS[@]}")
    cat <<MSG
Submitted job $jobid ($CPUS cpu, $MEM, $TIME on $PARTITION).
Claude is running in a tmux session on the node. Attach with:

    srun --jobid=$jobid --overlap --pty tmux attach -t claude

Detach again with ctrl-b d; the job keeps running. End it with:
    scancel $jobid          (log: logs/claude-$jobid.out)
MSG
    exit 0
fi

# ------------------------------------------------------------- interactive --
# salloc dies with your ssh connection - tmux on the login node prevents that.
if [ "$FORCE" = 0 ] && [ "$DETACH" = 0 ] && [ -z "${TMUX:-}" ] && [ -z "${STY:-}" ]; then
    cat >&2 <<'MSG'
Not inside tmux/screen: if your ssh drops, the allocation is cancelled and the
session dies. Start `tmux` here first, or use --detach, or pass --force.
MSG
    exit 1
fi

if [ "$SHELL_ONLY" = 1 ]; then
    remote="cd '$PROJECT_DIR' && . scripts/claude-env.sh && claude_net_check; exec bash -l"
else
    args=""
    for a in ${CLAUDE_ARGS[@]+"${CLAUDE_ARGS[@]}"}; do args+=" $(printf '%q' "$a")"; done
    remote="cd '$PROJECT_DIR' && . scripts/claude-env.sh && claude_net_check && exec claude$args"
fi

echo "Allocating $CPUS cpu / $MEM / $TIME on $PARTITION ..." >&2
exec salloc "${SLURM_OPTS[@]}" srun --pty --cpu-bind=none bash -lc "$remote"
