#!/bin/sh
#SBATCH --job-name=<NAME>_plot
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:20:00


. ~/miniforge3/etc/profile.d/conda.sh
conda activate plotting

ASSEMBLY=<NAME>

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/repeatmasker/${ASSEMBLY}_RM
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates 

cd $WORKDIR

python $CURATIONPATH/te_plotting2.py \
	-r ${ASSEMBLY}.fa.out.gz \
	-s ${ASSEMBLY}.summary.gz \
	-op ${ASSEMBLY} \
	-c LINE,SINE,LTR,DNA,RC \
	-proc 12
	
