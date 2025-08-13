#!/bin/bash
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -J sing
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=20
#SBATCH -p large
#SBATCH --time=5-00:00:00
#SBATCH --mem=100G

apptainer pull -F docker://cenococcum/fo_amplicons:latest
