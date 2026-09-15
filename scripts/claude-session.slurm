#!/bin/bash
# Hold a compute node so you can ssh in and work there.
#
#   mkdir -p logs
#   sbatch scripts/claude-session.slurm
#   cat logs/claude-<jobid>.out     # tells you which node, and how to get in
#
# Change the resources below to taste; scancel <jobid> when you're done.

#SBATCH --job-name=claude
#SBATCH --account=small_grains
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/claude-%j.out

cat <<MSG
node      : $(hostname -s)
job       : $SLURM_JOB_ID
resources : $SLURM_CPUS_PER_TASK cpu, $(( SLURM_MEM_PER_NODE / 1024 ))G, until $(squeue -h -j "$SLURM_JOB_ID" -o %e)

Get a shell on it:
    ssh $(hostname -s)
  or, if ssh to compute nodes is blocked:
    srun --jobid=$SLURM_JOB_ID --overlap --pty bash -l

Then:
    cd $SLURM_SUBMIT_DIR && claude

Release the node when done:  scancel $SLURM_JOB_ID
MSG

# Keep the allocation alive until the walltime runs out or you scancel it.
sleep infinity
