#!/bin/bash
#SBATCH --job-name=props
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00

. ~/miniforge3/etc/profile.d/conda.sh
conda activate 

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/plotting
cd $WORKDIR

python ../curation_templates/te_proportion_calculator.py \
	-d primary_outs/ \
	-m ../species_mapping_full.tsv \
	-t 50 \
	-o te_proportions.tsv


