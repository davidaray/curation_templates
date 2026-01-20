#!/bin/bash 
#SBATCH --job-name=<NAME>_hite 
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err 
#SBATCH --partition=nocona 
#SBATCH --mem-per-cpu=10G 
#SBATCH --nodes=1
#SBATCH --ntasks=3


. ~/miniforge3/etc/profile.d/conda.sh
conda activate

# Variables used in this script
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/hite
LIBRARYPATH=/lustre/scratch/daray/bat1k_TE_analyses/mammals.plus.covid_bats2.04072022.fa
GITPATH=/home/daray/gitrepositories/bioinfo_tools
USEARCHPATH=/lustre/work/daray/software

cd $DIR/te-aid || { echo "Directory $DIR/te-aid does not exist. Exiting."; exit 1; }

cd $WORKDIR

LIST="../primary_assemblies_all.txt" 
# Initialize a flag to track the first item
FIRST_ITEM=true

while IFS= read -r LINE; do
    if [ "$FIRST_ITEM" = true ]; then
		echo "Current assembly is $LINE. No previous assembly" 
		cd $WORKDIR/$LINE
		# Check if ${LINE}_final_library.fa file exists and is not empty
		if [ ! -e "${LINE}_final_library.fa" ] || [ ! -s "${LINE}_final_library.fa" ]; then
			echo "${LINE}_final_library.fa does not exist or is not empty."
		else
			echo "${LINE}_final_library.fa is present and not empty. Proceed."
			if [ -e ${LINE}_vs_mammals_usearch_60_hits.tsv ]; then
				rm ${LINE}_vs_mammals_usearch_60_hits.tsv
				echo "${LINE}_vs_mammals_usearch_60_hits.tsv exists, removing."
			fi

			echo "Run usearch on ${LINE}_final_library.fa"
			$USEARCHPATH/usearch11.0.667_i86linux32 \
				-usearch_global ${LINE}_final_library.fa \
				-db $LIBRARYPATH \
				-strand both \
				-threads 3 \
				-id 0.60 \
				-minsl 0.80 \
				-maxsl 1.2 \
				-maxaccepts 1 \
				-maxrejects 128 \
				-userfields query+target+id+ql+tl \
				-userout ${GENOME}_vs_mammals_usearch_60_hits.tsv

			conda activate pantera

			echo "Run hite_usearch_final.py on ${LINE}_final_library.fa"
			python $CURATIONPATH/hite_usearch_final.py \
				-t ${GENOME}_vs_mammals_usearch_60_hits.tsv \
				-n ${LINE}_final_library.fa \
				-m $LIBRARYPATH 
		fi
	else 
		TARGET_ASSEMBLY="$LINE"
		PREVIOUS_ASSEMBLY=""
		while read -r CURRENT_ITEM; do
			if [ "$CURRENT_ITEM" == "$TARGET_ITEM" ]; then
				echo "Current assembly is $CURRENT_ASSEMBLY"
				echo "Previous assembly is $PREVIOUS_ASSEMBLY"
				cd $WORKDIR/$LINE
				# Check if ${LINE}_final_library.fa file exists and is not empty
				if [ ! -e "${LINE}_final_library.fa" ] || [ ! -s "${LINE}_final_library.fa" ]; then
					echo "${LINE}_final_library.fa does not exist or is not empty."
				else
					echo "${LINE}_final_library.fa is present and not empty. Proceed."
					if [ -e ${LINE}_vs_mammals_usearch_60_hits.tsv ]; then
						rm ${LINE}_vs_mammals_usearch_60_hits.tsv
						echo "${LINE}_vs_mammals_usearch_60_hits.tsv exists, removing."
					fi
				echo "Run usearch on ${LINE}_final_library.fa"
				$USEARCHPATH/usearch11.0.667_i86linux32 \
					-usearch_global ${LINE}_final_library.fa \
					-db $LIBRARYPATH \
					-strand both \
					-threads 3 \
					-id 0.60 \
					-minsl 0.80 \
					-maxsl 1.2 \
					-maxaccepts 1 \
					-maxrejects 128 \
					-userfields query+target+id+ql+tl \
					-userout ${GENOME}_vs_mammals_usearch_60_hits.tsv

				conda activate pantera

				echo "Run hite_usearch_60.py on ${GENOME}_${TE}_collection-_rep.fa"
				python $CURATIONPATH/hite_usearch_final.py \
					-t ${GENOME}_vs_mammals_usearch_60_hits.tsv \
					-n ${LINE}_final_library.fa \
					-m ../concatenated_library_${PREVIOUS_ASSEMBLY}.fa
				PREVIOUS_ASSEMBLY="$CURRENT_ASSEMBLY"
			fi
		done < "$LIST"

done < "$LIST"


