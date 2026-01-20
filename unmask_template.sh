#!/bin/sh
#SBATCH --job-name=<NAME>_unmask
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1


. ~/miniforge3/etc/profile.d/conda.sh
conda activate 

ASSEMBLY=<NAME>
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/te_curations/$ASSEMBLY
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies

python /lustre/scratch/daray/bat1k_TE_analyses/curation_templates/unmask.py \
	-i $GZPATH/${ASSEMBLY}.fa.gz \
	-o $ASSEMBLY/${ASSEMBLY}.fa.gz
