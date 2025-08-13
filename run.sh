#!/bin/bash
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -J amplicons
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=64
#SBATCH -p large
#SBATCH --time=5-00:00:00
#SBATCH --mem=1000000

source /groups/jmd20jns/miniforge/etc/profile.d/conda.sh
source /groups/jmd20jns/miniforge/etc/profile.d/mamba.sh

mamba activate fo_amplicons 

Rscript -e "targets::tar_make()"

#apptainer exec \
#    --bind ${cmdstan_tmpdir}/:/home/user/.cmdstan/ \
#    --bind /groups/jmd20jns/amplicons/:/work/ \
#    --pwd /work fo_amplicons_latest.sif Rscript -e "targets::tar_make()"
