#!/bin/bash
#PBS -q qprod
#PBS -N all_samples
#PBS -l select=1:ncpus=36
#PBS -A OPEN-25-50
#PBS -M p22016@student.osu.cz
#PBS -m aea
#PBS -j oe
#PBS -o run_log/

cd /scratch/project/open-25-50/ExtraMyelomaTrans_DTU

source activate drimseq

snakemake --unlock --use-conda -s workflows/Snakefile
snakemake --cores 36 --use-conda -s workflows/Snakefile
