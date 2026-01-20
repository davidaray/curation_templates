#!/bin/bash 
#SBATCH --job-name=<NAME>_mLTR
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=12

. ~/miniforge3/etc/profile.d/conda.sh
conda activate MegaLTR

TAXON=<NAME>
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/te_curations/$TAXON
mkdir -p $WORKDIR/megaltr
cd $WORKDIR/megaltr
SOFTWARE=/lustre/work/daray/software

#Run megaltr
bash $SOFTWARE/MegaLTR/MegaLTR.sh \
	-A 1 \
	-t 12 \
	-P ${TAXON}_results \
	-T ${TAXON}.tRNA.fa \
	-F $WORKDIR/assemblies_dir/${TAXON}.fa \
	-I $WORKDIR/trnascan \
	-B rexdb-metazoa \
	-E rexdb-metazoa 


