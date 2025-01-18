#!/bin/bash

# First open tmux. We'll want to run metamap in the background via tmux (or nohup) given that we'll run snakemake on the cluster (via SLURM).
tmux new -s magscreen_run3

# Load modules
module load miniconda/3

# Activate mamba and snakemake environments. Use "conda activate" here because for some reason, using "source activate" doesn't allow source activating snakemake later.
conda activate /home/hl636/rds/hpc-work/conda_mamba
source activate snakemake

# Number of jobs
njobs=42 #sapphire
#njobs=35 #icelake
#njobs=100 #icelake

# Run genofan snakemake pipeline. For reference, see: https://github.com/alexmsalmeida/genofan, https://snakemake.readthedocs.io/en/stable/executing/cluster.html, https://snakemake.readthedocs.io/en/stable/executing/cli.html
# Run locally (interactively). Here -j refers to number of cores. -k, --keep-going goes on with independent jobs if a job fails.
snakemake --use-conda -k -j 1
# Run in parallel on cluster (via SLURM). Here -j refers to number of jobs. -k, --keep-going goes on with independent jobs if a job fails.
snakemake --use-conda -k -j ${njobs} --cluster-config cluster.yml --cluster 'sbatch -A JACOBS-SL3-CPU -p sapphire --time=12:00:00 --ntasks={cluster.nCPU} --mem={cluster.mem} -o {cluster.output}'

# To deactivate
conda deactivate
#mamba deactivate

# List running tmux sessions. For reference, see: https://gist.github.com/MohamedAlaa/2961058
#tmux ls
# Attach to named tmux session
#tmux a -t magscreen_run3
# Killed named tmux session
#tmux kill-session -t magscreen_run1

