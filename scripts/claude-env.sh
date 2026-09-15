# Shared setup for Claude Code sessions on compute nodes.
# Sourced by claude-node.sh and claude-job.slurm; not meant to be run directly.

# Repo root, regardless of where the caller invoked us from.
PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PROJECT_DIR

# claude is a native binary in ~/.local/bin, and its auth + config live under
# $HOME, which the compute nodes mount - so there is nothing to copy per job.
export PATH="$HOME/.local/bin:$PATH"

# The Makefile loads its own modules (.SHELLFLAGS := -lc) and reaches into
# $(CONDA) by absolute path, so a session only needs the conda env on PATH for
# ad-hoc commands. Harmless if the env isn't built yet.
[ -d "$HOME/.conda/envs/rnaseq/bin" ] && export PATH="$HOME/.conda/envs/rnaseq/bin:$PATH"

# Compute nodes on some clusters have no route off-site; Claude needs one.
# Any HTTP response proves TCP+TLS reached the API.
claude_net_check() {
    if curl -sS -m 15 -o /dev/null https://api.anthropic.com/ 2>/dev/null; then
        return 0
    fi
    cat >&2 <<'MSG'
!! No route to api.anthropic.com from this compute node.
   Claude cannot run here until that is fixed. Options:
     - SOCKS proxy back through the login node (needs passwordless ssh
       compute -> login, and Claude must honour the proxy var):
           ssh -f -N -D 8888 ceres-login.scinet.usda.gov
           export HTTPS_PROXY=socks5h://127.0.0.1:8888
     - Ask VRSC (scinet_vrsc@usda.gov) which partitions have outbound HTTPS.
MSG
    return 1
}
