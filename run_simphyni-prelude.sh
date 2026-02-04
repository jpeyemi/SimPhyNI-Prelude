#!/bin/bash

# Edit for your HPC
# Increase time for raxml tree building
#SBATCH -p sched_any
#SBATCH --job-name=simphyni
#SBATCH --output=simphyni_out_%j.txt
#SBATCH --error=simphyni_err_%j.txt
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000

snakemake --snakefile Snakefile.py \
	--rerun-incomplete \
	--printshellcmds \
	--rerun-triggers mtime \
	--profile cluster_profile \
	--dry-run \
	--config tree_method=poppunk

echo Done!!!
