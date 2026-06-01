#!/bin/bash
#SBATCH --job-name=nf_scRNA
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

nextflow run main.nf -profile singularity -resume
