#!/bin/sh
#SBATCH --job-name=pre_pantera
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1

. ~/miniforge3/etc/profile.d/conda.sh
conda activate

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/repeatmasker

cd $WORKDIR

cat ../primary_assemblies_all.txt | while read i; do 
	cd ${i}_RM
	sbatch doLift.sh
done

	
