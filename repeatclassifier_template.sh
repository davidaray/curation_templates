#!/bin/bash 
#SBATCH --job-name=<NAME>_classify
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=60G

module --ignore-cache load gcc/10.1.0 r/4.0.2
. ~/miniforge3/etc/profile.d/conda.sh
conda activate repeatmodeler

TAXON=<NAME>
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/te_curations/$TAXON
mkdir -p $WORKDIR/repeatclassifier
cd $WORKDIR/extensions/final_consensuses

cat ${TAXON}*.fa >$WORKDIR/repeatclassifier/${TAXON}_extended_rep.fa

#Run RepeatClassifier
cd $WORKDIR/repeatclassifier
RepeatClassifier -consensi ${TAXON}_extended_rep.fa


