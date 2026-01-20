#!/bin/bash 
#SBATCH --job-name=<NAME>_hite 
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err 
#SBATCH --partition=quanah 
#SBATCH --mem-per-cpu=10G 
#SBATCH --nodes=1
#SBATCH --ntasks=10

module --ignore-cache load gnu7/7.3.0 R/3.6.1-gnu7
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

cd $DIR
pwd

cp confident_TE.cons.fa ${GENOME}_confident_TE.cons.fa

#$RMASKERPATH/RepeatMasker -lib $LIBRARYPATH -nolow -pa 10 ${GENOME}_confident_TE.cons.fa

if [ -e ${GENOME}_usearch_80_hits.tsv ] 
  then rm ${GENOME}_usearch_80_hits.tsv
  echo ${GENOME}_usearch_80_hits.tsv" exists, removing." 
fi

$USEARCHPATH/usearch11.0.667_i86linux32 \
	-usearch_global ${GENOME}_confident_TE.cons.fa \
	-db ../../mammals.plus.covid_bats2.04072022.fa \
	-strand both \
	-threads 1 \
	-id 0.80 \
	-minsl 0.80 \
	-maxsl 1.2 \
	-maxaccepts 1 \
	-maxrejects 128 \
	-userfields query+target+id+ql+tl \
	-userout ${GENOME}_usearch_80_hits.tsv

conda activate pantera
python ../../curation_templates/hite_usearch.py -u ${GENOME}_usearch_80_hits.tsv -c ${GENOME}_confident_TE.cons.fa -m ../mammals.plus.covid_bats2.04072022.fa

sed "s/<TAXON>/$GENOME/g" $CURATIONPATH/TEcurate_template_hite_q.sh >${GENOME}_TEcurate_hite_q.sh

conda activate curate
sh ${GENOME}_TEcurate_hite_q.sh

pwd

rm -rf $DIR/RM*/

