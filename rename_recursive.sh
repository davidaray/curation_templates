#!/bin/bash 
#SBATCH --job-name=mChoMin1.1.hap1_RM2bed
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --mem-per-cpu=10G
#SBATCH --nodes=1
#SBATCH --ntasks=10

## Usage:
## You will need to replace aJam with the genome assembly identifier as appropriate.
## DIR may need modifying.
## sed "s|DIRPATH|path to your directory|g" template_RM_rm2bed.sh ><new file name>
## You will need to include the path to your genome assembly, GENOMEPATH.
## Genome assembly should be named aJam.fa
## If running lots of genomes, here's a quick way to do it using a 'list.txt' file that has all of your taxon names in it.
## cat list.txt | while read i; do sed "s/aJam/$i/g" <new file name> >${i}_RM2bed.sh; done
## Then, submit the scripts using a similar loop.

#This script uses a python program. You need to make sure python is loaded and ready to go by activating conda.
#Your conda installation must have biopython, pandas, and pyfaidx installed.
. ~/miniforge3/etc/profile.d/conda.sh
conda activate pigz

#makes variables used in this script and in the other
OLDNAME=mCenSen1.hap1
NEWNAME=mCenSen1.1.hap1
DIR=/lustre/scratch/daray/bat1k_TE_analyses/repeatmasker
GITPATH=/home/daray/gitrepositories/bioinfo_tools/
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates

#Creates a working directory and then goes into it to make all the files
#mkdir /lustre/scratch/acorreap/zoonomia_rm/$RUNTYPE
cd $DIR

gunzip $NEWNAME/*.out.gz 
gunzip $NEWNAME/*.align.gz
gunzip $NEWNAME/*.summary.gz
gunzip $NEWNAME/*.divsum.gz

python rename_recursive.py . $OLDNAME $NEWNAME 

pigz -p 10 $NEWNAME/*.out
pigz $NEWNAME/*.align
pigz $NEWNAME/*.summary
pigz $NEWNAME/*.divsum
