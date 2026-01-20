#!/bin/bash 
#SBATCH --job-name=mHypMon1.pri_rmod
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=64


. ~/miniforge3/etc/profile.d/conda.sh
conda activate

TAXON=mHypMon1.pri

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses
cd $WORKDIR/te_curations

mkdir $WORKDIR/te_curations/${TAXON}/

python ../curation_templates/rmodel.py \
	-g /lustre/scratch/daray/bat1kdatafreeze/final_assemblies/${TAXON}.fa.gz \
	-p 64 \
	-b 50 \
	-w $WORKDIR/te_curations/${TAXON}/
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
