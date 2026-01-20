#!/bin/bash 
#SBATCH --job-name=<NAME>_hite 
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err 
#SBATCH --partition=nocona 
#SBATCH --mem-per-cpu=10G 
#SBATCH --nodes=1
#SBATCH --ntasks=1


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

cd $DIR/te-aid || { echo "Directory $DIR/te-aid does not exist. Exiting."; exit 1; }

LIST="DNA NOHIT LINE LTR hite_LTRs RC"
for TE in $LIST; do
    cd $DIR/te-aid
    if [ -z "$(ls -A $DIR/te-aid/$TE 2>/dev/null)" ]; then
        echo "te-aid/$TE is empty. Move to next TE type."
        continue
    else
        echo "te-aid/$TE is not empty. Proceed."

        cat $DIR/te-aid/$TE/*-_rep.fa >${GENOME}_${TE}_collection-_rep.fa

        if [ -e ${GENOME}_${TE}_usearch_60_hits.tsv ]; then
            rm ${GENOME}_${TE}_usearch_60_hits.tsv
            echo "${GENOME}_${TE}_usearch_60_hits.tsv exists, removing."
        fi

        echo "Run usearch on ${GENOME}_${TE}_collection-_rep.fa"
        $USEARCHPATH/usearch11.0.667_i86linux32 \
            -usearch_global ${GENOME}_${TE}_collection-_rep.fa \
            -db $LIBRARYPATH \
            -strand both \
            -threads 1 \
            -id 0.60 \
            -minsl 0.80 \
            -maxsl 1.2 \
            -maxaccepts 1 \
            -maxrejects 128 \
            -userfields query+target+id+ql+tl \
            -userout ${GENOME}_${TE}_usearch_60_hits.tsv

        conda activate pantera

        echo "Run hite_usearch_60.py on ${GENOME}_${TE}_collection-_rep.fa"
        python $CURATIONPATH/hite_usearch_60.py \
            -u ${GENOME}_${TE}_usearch_60_hits.tsv \
            -c ${GENOME}_${TE}_collection-_rep.fa \
            -m $LIBRARYPATH \
            -t ${TE}

        cat ${GENOME}_${TE}_cons.renamed.fa \
            ${GENOME}_${TE}_class_usearch.fa \
            ${GENOME}_${TE}_family_usearch.fa \
            ${GENOME}_${TE}_same_usearch.fa \
            >${GENOME}_${TE}_preliminary.fa
    fi
	cd $DIR/
	cat te-aid/${GENOME}_*_preliminary.fa >${GENOME}_TEcurate_preliminary.fa
done

cd $DIR/

rm ${GENOME}_full_preliminary_library.fa ${GENOME}_final_library.fa ${GENOME}_duplicates.fa

cat ${GENOME}_family_usearch.fa ${GENOME}_TEcurate_preliminary.fa > ${GENOME}_full_preliminary_library.fa
			
python $CURATIONPATH/duplicate_check.py -i ${GENOME}_full_preliminary_library.fa
