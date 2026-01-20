#!/bin/bash 
#SBATCH --job-name=RM2bed
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --mem-per-cpu=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1

. ~/miniforge3/etc/profile.d/conda.sh
conda activate biopython

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/bed_files
CURATE=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates
GITPATH=/home/daray/gitrepositories/bioinfo_tools/
RMPATH=/lustre/scratch/daray/bat1k_TE_analyses/repeatmasker

cd $WORKDIR
for FILE in *.fa.out.gz; do
	NAME=$(basename $FILE .fa.out.gz)
	python ${GITPATH}/RM2bed_hubley.py --out_dir . --split class --out_prefix ${NAME} --ovlp lower_divergence ${RMPATH}/${NAME}_RM/${NAME}.fa.align.gz
	sort -k1,1 -k2,2n ${NAME}_rm.bed >${NAME}_rm_sorted.bed
done


