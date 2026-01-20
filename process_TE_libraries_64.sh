#!/bin/bash 
#SBATCH --job-name=process
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err 
#SBATCH --partition=nocona 
#SBATCH --nodes=1
#SBATCH --ntasks=64


cd /lustre/scratch/daray/bat1k_TE_analyses/hite

. ~/miniforge3/etc/profile.d/conda.sh
conda activate
mkdir process_64

python ../curation_templates/process_TE_libraries.py -i ../primary_assemblies_all.txt -l ../mammals.plus.covid_bats2.04072022.fa -o process_64
