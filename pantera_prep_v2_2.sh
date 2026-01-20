#!/bin/sh
#SBATCH --job-name=prep_pantera
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=36


. ~/miniforge3/etc/profile.d/conda.sh
conda activate pantera

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/pantera
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies 
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates 

cd $WORKDIR

python $CURATIONPATH/extract_for_pantera_v2.py \
	-l species_pairs_2.txt \
	-d $GZPATH \
	-o $WORKDIR/chromosomes2 \
	-t 36 \
	--minsize 10000000 \
	-s 0.65 \
	--sizediff 0.1


