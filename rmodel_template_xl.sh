#!/bin/bash 
#SBATCH --job-name=<NAME>_rmod
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=xlquanah
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --time=50:00:00
#SBATCH --account=xlquanah

. ~/miniforge3/etc/profile.d/conda.sh
conda activate

TAXON=<NAME>

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses
cd $WORKDIR/te_curations

mkdir $WORKDIR/te_curations/${TAXON}/

python ../curation_templates/rmodel.py \
	-g /lustre/scratch/daray/bat1kdatafreeze/final_assemblies/${TAXON}.fa.gz \
	-p 36 \
	-b 50 \
	-w $WORKDIR/te_curations/${TAXON}/
