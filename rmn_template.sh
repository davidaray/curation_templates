#!/bin/sh
#SBATCH --job-name=<NAME>_RMN
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=12


GENOME=<NAME>

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/repeatmasker_nextflow/${GENOME}_RMN
RMPATH=/lustre/work/daray/software/RepeatMasker/RepeatMasker
SOFTWARE=/lustre/work/daray/software
RMNEXTFLOWPATH=/lustre/work/daray/software/RepeatMasker_Nextflow/
LIBPATH=/lustre/scratch/daray/bat1k_TE_analyses/bat1k_named_20241126.fa 
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies

mkdir -p $WORKDIR

cd $WORKDIR

$SOFTWARE/nextflow run $RMNEXTFLOWPATH/RepeatMasker_Nextflow.nf \
  --inputSequence $GZPATH/${GENOME}.fa.gz \
  --inputLibrary $LIBPATH \
  --cluster nocona \
  --outputDir $WORKDIR