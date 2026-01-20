#!/bin/sh
#SBATCH --job-name=prep_pantera
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=36

## This script assumes that assembly.faidx.sh to prep .faidx files for the comparisons to be performed.

## Utilizes extract_for_pantera_v2.py

. ~/miniforge3/etc/profile.d/conda.sh
conda activate pantera

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies 
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates 

cd $WORKDIR/pantera2

python $CURATIONPATH/extract_for_pantera_v2.py \
	-l $WORKDIR/species_pairs.txt \
	-d $GZPATH \
	-o $WORKDIR/pantera2/chromosomes \
	-t 36 \
	--minsize 10000000 \
	-s 0.65 \
	--sizediff 0.1


