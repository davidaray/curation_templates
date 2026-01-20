#!/bin/bash
#SBATCH --job-name=teclass2_cuda
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=matador
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1 
#SBATCH --time=2:00:00

#####
#https://github.com/IOB-Muenster/TEclass2
#
#Training model prep --> tec;ass_models_and run.sh
#####

. ~/miniforge3/etc/profile.d/conda.sh
conda activate TEclass2_cuda
#conda activate TEclass2

CUDA_LAUNCH_BLOCKING=1

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/teclass2_nancy
mkdir -p $WORKDIR/data/databases
mkdir -p $WORKDIR/models
cd $WORKDIR

CURATION_DIR=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates
BIN=/home/daray/TEclass2
LIBRARY=/lustre/scratch/daray/bat1k_TE_analyses/final_hite_library_20240919.fa

#: << 'EOF'
DATE=$(date)
echo -e "Prep start time: $DATE\n"

#List of TE models to use
TE_LIST=(
    "LTR/Copia"
    "LTR/ERVK_ltr"
    "LTR/ERVK_internal"
    "LTR/ERV1_ltr"
    "LTR/ERV1_internal"
    "LTR/ERVK"
    "LTR/ERV1"
    "LTR/ERV"
    "LTR/Gypsy"
    "LTR/ERVL"
    "LTR/LTR"
    "DNA/PiggyBac"
    "DNA/piggyBac"
    "DNA/hAT"
    "RC/Helitron"
    "DNA/Maverick"
    "DNA/Crypton"
    "DNA/Merlin"
    "DNA/PIF"
    "DNA/TcMar"
    "DNA/Transib"
    "LINE/Jockey"
    "LINE/L1"
    "LINE/L2"
    "LINE/HAL"
    "LINE/I"
    "LINE/I-Jockey"
    "LINE/Penelope"
    "LINE/Proto2"
    "LINE/R1"
    "LINE/R2"
    "LINE/RTE"
    "SINE/Ves"
    "SINE/tRNA"
    "SINE/Rhin"
    "SINE/Meg"
    "SINE/Alu"
    "SINE/B4"
    "SINE/ID"
    "SINE/5S"
    "LTR/Pao"
)


[[ -e $WORKDIR/data/bats_v1.fa ]] && rm $WORKDIR/data/bats_v1.fa

for TE in "${TE_LIST[@]}"; do 
	
	NEWLABEL=$(echo "$TE" | sed "s|/|-|g") 
	CLASS=$(echo "$NEWLABEL" | cut -d'-' -f1)
	FAMILY=$(echo "$NEWLABEL" | cut -d'-' -f2)	
	
	echo "NEWLABEL = ${NEWLABEL}"
	echo "CLASS = ${CLASS}"
	echo "FAMILY = ${FAMILY}"
	
	python $CURATION_DIR/pull_using_header_text2.py \
		-i "$LIBRARY" \
		-o "${NEWLABEL}_pulled.fa" \
		-s "$TE" || exit 1
	
	python $CURATION_DIR/replace_fasta_headers.py \
		-i "${NEWLABEL}_pulled.fa" \
		-o "${NEWLABEL}_model.fa" \
		-n "$CLASS $FAMILY" || exit 1

	cat "${NEWLABEL}_model.fa" >> $WORKDIR/data/bats_v1.fa
	
	rm "${NEWLABEL}_model.fa"
	rm "${NEWLABEL}_pulled.fa"

done
#EOF
#: << 'EOF'

DATE=$(date)
echo -e "Prep end time: $DATE\n"

echo -e "Create_config start:  $DATE\n"

cp -r $BIN/*/ .
cp $BIN/config.yml original_config.yml
sed "s/model_name: 'clust_cats_16'/model_name: 'bats_v1'/g" original_config.yml >config.yml
sed -i "s|model_save_path: 'models/'|model_save_path: '${WORKDIR}/models'|g" config.yml
sed -i "s|data/clust_cats_16.fa|${WORKDIR}/data/bats_v1.fa|g" config.yml
sed -i "s|dataset_path: 'data/databases/clust_cats_16'|dataset_path: '${WORKDIR}/data/databases/bats_v1'|g" config.yml

cp $WORKDIR/data/bats_v1.fa $WORKDIR/models

python $BIN/TEclass2.py --database -c config.yml >database_output.txt
: << 'EOF'

# Alter config.yml to elements in database with at least 500 entries. 
# New config file = config2.yml
python update_config.py
EOF

DATE=$(date)
echo -e "Create_config end:  $DATE\n"

# Bash version of update_config.py.
# Extracting relevant entries with values > 500
te_keywords=$(awk '$NF > 500 {print $1}' database_output.txt | grep -v 'total_len\|Dataset' | tr '\n' ',' | sed 's/,$//')

# Preparing the te_keywords line for config.yml
te_keywords="te_keywords: [\"$(echo $te_keywords | sed 's/,/","/g')\"]"

# Replace the te_keywords line in config.yml and save it as config1.yml
awk -v new_keywords="$te_keywords" '
/^te_keywords:/ {print new_keywords; next} 
{print}
' config.yml > config1.yml

echo "config.yml updated and saved as config1.yml"

python $BIN/TEclass2.py --train -c config1.yml

DATE=$(date)
echo -e "Training end:  $DATE\n"

#: << 'EOF'
echo -e "Classify start:  $DATE\n"

python $BIN/TEclass2.py --classify -c config1.yml \
	-o $WORKDIR \
	-f $LIBRARY

DATE=$(date)
echo -e "Classify end:  $DATE\n"

python $CURATION_DIR/filter_teclass.py \
	-if output_from_webserver/ \
	-f 0.7 \
	-o teclass_webserver_filtered_0.7.txt

python $CURATION_DIR/pull_using_header_text.py \
	-i $LIBRARY \
	-o hite_unknowns.fa \
	-s $SEARCH

python $CURATION_DIR/clean_sequences.py \
	-i hite_unknowns.fa \
	-o hite_clean_library.fa \
	-d invalid_seqs

python $CURATION_DIR/replace_fasta_headers_teclass2.py \
	-f hite_clean_library.fa \
	-r teclass_webserver_filtered_0.7.txt \
	-o hite_clean_teclass.fa

python $CURATION_DIR/replace_original_headers.py \
	-f1 $LIBRARY \
	-f2 hite_clean_teclass.fa \
	-o hite_teclass_20240919.fa
#EOF



