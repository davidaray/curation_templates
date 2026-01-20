#!/bin/sh
#SBATCH --job-name=prep_pantera
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1


. ~/miniforge3/etc/profile.d/conda.sh
conda activate 

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/pantera 
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies 
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates 

cd WORKDIR

python $CURATIONPATH/extract_for_pantera.py \
	-l species_pairs.txt \
	-d $GZPATH \
	-o $WORKDIR/chromosomes
