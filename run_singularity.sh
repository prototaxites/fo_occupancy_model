#!/bin/bash
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -J amplicons
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=64
#SBATCH -p large
#SBATCH --time=5-00:00:00
#SBATCH --mem=1000G

if [ ! -d tmp/ ]; then
    mkdir tmp
fi

if [ ! -f tmp/overlay.img ]; then
    apptainer overlay create --sparse --size 2048 tmp/overlay.img
fi

apptainer exec \
    --userns --cleanenv \
    --overlay tmp/overlay.img \
    --bind /groups/jmd20jns/amplicons/:/work/ \
    --pwd /work fo_amplicons_latest.sif Rscript -e "targets::tar_make()"
