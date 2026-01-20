#!/bin/sh
#SBATCH --job-name=pca
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=xlquanah
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=1:00:00
#SBATCH --account=xlquanah


. ~/miniforge3/etc/profile.d/conda.sh
conda activate plotting

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/pca
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates 

cd $WORKDIR

python generate_pcas4.py \
	--out_dir ../bed_files \
	--mapping_file ../species_mapping_full.txt \
	--output_dir . \
	--num_proc 12 
	
