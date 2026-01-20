#!/bin/bash
#SBATCH --job-name=final_process
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=3

# Load conda environment
. ~/miniforge3/etc/profile.d/conda.sh
conda activate pantera

# Variables used in this script
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/hite
MAMMALPATH=/lustre/scratch/daray/bat1k_TE_analyses/mammal_04072022_no_duplicates.fa
GITPATH=/home/daray/gitrepositories/bioinfo_tools
USEARCHPATH=/lustre/work/daray/software
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates

cd $WORKDIR

LIST="../primary_assemblies_all.txt" 
FIRST_ITEM=true
PREVIOUS_ASSEMBLY=""

while IFS= read -r LINE; do
    echo "Processing assembly: $LINE"
    echo ""
    cd $WORKDIR/$LINE
    
    if [ "$FIRST_ITEM" = true ]; then
        echo "No previous assembly" 
        echo ""
        DB_PATH="$MAMMALPATH"
    else
        echo "Previous assembly is $PREVIOUS_ASSEMBLY"
        echo ""
        KNOWN_LIBRARY="../concatenated_library_${PREVIOUS_ASSEMBLY}.fa"
		DB_PATH=$KNOWN_LIBRARY
        
        if [ ! -e "$KNOWN_LIBRARY" ]; then
            echo "Error: $KNOWN_LIBRARY does not exist."
            echo ""
            exit 1
        fi
        
        DB_PATH="$KNOWN_LIBRARY"
    fi

    if [ ! -e "${LINE}_final_library.fa" ] || [ ! -s "${LINE}_final_library.fa" ]; then
        echo "${LINE}_final_library.fa does not exist or is empty."
        echo ""
        continue
    fi

    echo "${LINE}_final_library.fa is present and not empty. Proceed."
    
    if [ -e "${LINE}_vs_mammals_usearch_60_hits.tsv" ]; then
        rm "${LINE}_vs_mammals_usearch_60_hits.tsv"
        echo "${LINE}_vs_mammals_usearch_60_hits.tsv exists, removing."
        echo ""
    fi

    echo "Run usearch on ${LINE}_final_library.fa"
    echo ""
    $USEARCHPATH/usearch11.0.667_i86linux32 \
        -usearch_global "${LINE}_final_library.fa" \
        -db "$DB_PATH" \
        -strand both \
        -threads 3 \
        -id 0.60 \
        -minsl 0.80 \
        -maxsl 1.2 \
        -maxaccepts 1 \
        -maxrejects 128 \
        -userfields query+target+id+ql+tl \
        -userout "${LINE}_vs_mammals_usearch_60_hits.tsv"
    
    if [ "$FIRST_ITEM" = true ]; then
        python "$CURATIONPATH/hite_usearch_final.py" \
            -t "${LINE}_vs_mammals_usearch_60_hits.tsv" \
            -n "${LINE}_final_library.fa" \
            -m "$DB_PATH"
		
		python "$CURATIONPATH/remove_duplicates.py" \
			-i "$WORKDIR/concatenated_library_${LINE}.fa" \
			-o "$WORKDIR"

        FIRST_ITEM=false
    else
        python "$CURATIONPATH/hite_usearch_final.py" \
            -t "${LINE}_vs_mammals_usearch_60_hits.tsv" \
            -n "${LINE}_final_library.fa" \
            -m "$DB_PATH"

		python "$CURATIONPATH/remove_duplicates.py" \
			-i "$WORKDIR/concatenated_library_${LINE}.fa" \
			-o "$WORKDIR"
    fi
    
    PREVIOUS_ASSEMBLY="$LINE"

done < "$LIST"

