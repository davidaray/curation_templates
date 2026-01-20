#!/bin/bash 
#SBATCH --job-name=<NAME>_hite 
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err 
#SBATCH --partition=nocona 
#SBATCH --mem-per-cpu=10G 
#SBATCH --nodes=1
#SBATCH --ntasks=10


. ~/miniforge3/etc/profile.d/conda.sh
conda activate

# Variables used in this script
GENOME=<NAME>
RUNTYPE=${GENOME}_RM
DIR=/lustre/scratch/daray/bat1k_TE_analyses/hite/$GENOME
LIBRARYPATH=/lustre/scratch/daray/bat1k_TE_analyses/mammals.plus.covid_bats2.04072022.fa
HITELIBPATH=$DIR/HiTE.out
GITPATH=/home/daray/gitrepositories/bioinfo_tools
RMASKERPATH=/lustre/work/daray/software/RepeatMasker
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates
USEARCHPATH=/lustre/work/daray/software

cd $DIR/te-aid
cat $DIR/te-aid/NOHIT/*-_rep.fa >${GENOME}_nohit_collection-_rep.fa

if [ -e ${GENOME}_usearch_60_hits.tsv ] 
  then rm ${GENOME}_usearch_60_hits.tsv
  echo ${GENOME}_usearch_60_hits.tsv" exists, removing." 
fi

$USEARCHPATH/usearch11.0.667_i86linux32 \
	-usearch_global ${GENOME}_nohit_collection-_rep.fa \
	-db /lustre/scratch/daray/bat1k_TE_analyses/mammals.plus.covid_bats2.04072022.fa \
	-strand both \
	-threads 1 \
	-id 0.60 \
	-minsl 0.80 \
	-maxsl 1.2 \
	-maxaccepts 1 \
	-maxrejects 128 \
	-userfields query+target+id+ql+tl \
	-userout ${GENOME}_usearch_60_hits.tsv

conda activate pantera
python $CURATIONPATH/hite_usearch_60.py -u ${GENOME}_usearch_60_hits.tsv -c ${GENOME}_nohit_collection-_rep.fa -m $LIBRARYPATH

cd $DIR
cat te-aid/${GENOME}_nohit_TE.cons.renamed.fa \
	te-aid/${GENOME}_class_usearch.fa \
	te-aid/${GENOME}_family_usearch.fa \
	${GENOME}_family_usearch.fa \
	>${GENOME}_preliminary_library.fa

python $CURATIONPATH/duplicate_check.py -i ${GENOME}_preliminary_library.fa
