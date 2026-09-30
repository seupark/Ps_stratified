#!/bin/bash
#SBATCH --job-name=ps_strat_smk
#SBATCH --partition=slim16
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err

set -euo pipefail

# --- activate the conda environment ---
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ps_stratified

# --- move to the workflow root (where config/ lives) ---
cd ~/Ps_stratified

# --- run ---
# --rerun-incomplete: safely resume files left half-written by an earlier crash
# -p: print shell commands as they run (useful in the log)
snakemake --cores "${SLURM_CPUS_PER_TASK}" --rerun-incomplete -p
