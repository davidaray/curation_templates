#!/bin/bash 
#SBATCH --job-name=aKon_tRNA
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1

. ~/miniforge3/etc/profile.d/conda.sh
conda activate trnascan

TAXON=aKon
WORKDIR=/lustre/scratch/daray/protists/te_curations/$TAXON
mkdir -p $WORKDIR/trnascan
cd $WORKDIR/trnascan
SOFTWARE=/lustre/work/daray/software

#Run trnascan
tRNAscan-SE -o ${TAXON}.out $WORKDIR/assemblies_dir/${TAXON}.fa

sed -e 1,3d ${TAXON}.out | awk '{print $1 "\t" $3 "\t" $4 "\t" "aKon"_$1"_"$2"_"$5}' >${TAXON}.tRNA.tmp 
awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' ${TAXON}.tRNA.tmp > ${TAXON}.tRNA.bed 
rm ${TAXON}.tRNA.tmp

conda activate bedtools

bedtools getfasta -fo ${TAXON}.tRNA.fa -fi $WORKDIR/assemblies_dir/${TAXON}.fa -name -bed ${TAXON}.tRNA.bed 

