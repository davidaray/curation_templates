#!/bin/bash
#SBATCH --job-name=prep_model
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=matador
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=2 

#####
#https://github.com/IOB-Muenster/TEclass2
#
#Training model prep --> prep_teclass_models.sh
#####

DATE=$(date)
echo -e "Prep start time: $DATE\n"

. ~/miniforge3/etc/profile.d/conda.sh
conda activate TEclass2

WORK_DIR=/lustre/scratch/daray/bat1k_TE_analyses/teclass2
mkdir -p $WORK_DIR/data/databases
cd $WORK_DIR

CURATION_DIR=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates
BIN=/home/daray/TEclass2
LIBRARY=final_hite_library_20240919.fa

: << 'EOF'

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


[[ -e $WORK_DIR/data/teclass_models.fa ]] && rm $WORK_DIR/data/teclass_models.fa

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

	cat "${NEWLABEL}_model.fa" >> $WORK_DIR/data/teclass_models.fa
	
	rm "${NEWLABEL}_model.fa"
	rm "${NEWLABEL}_pulled.fa"

done
EOF
DATE=$(date)
echo -e "Prep end time: $DATE\n"

echo -e "Create_config start:  $DATE\n"

#cp -r $BIN/*/ .
#cp $BIN/config.yml original_config.yml
#sed "s/model_name: 'clust_cats_16'/model_name: 'bats_v1'/g" original_config.yml >config.yml
#sed -i "s|model_save_path: 'models/'|model_save_path: '$WORK_DIR/models'|g" config.yml
#sed -i "s|data/clust_cats_16.fa|$WORK_DIR/data/teclass_models.fa|g" config.yml
#sed -i "s|dataset_path: 'data/databases/clust_cats_16'|dataset_path: '$WORK_DIR/data/databases/bats_v1'|g" config.yml

#python $BIN/TEclass2.py --database -c config.yml 

DATE=$(date)
echo -e "Create_config end:  $DATE\n"

echo -e "Training start:  $DATE\n"

python $BIN/TEclass2.py --train -c config.yml

DATE=$(date)
echo -e "Training end:  $DATE\n"

echo -e "Classify start:  $DATE\n"

python $BIN/TEclass2.py --classify -c config.yml \
	-o $WORK_DIR \
	-f $LIBRARY

DATE=$(date)
echo -e "Classify end:  $DATE\n"



