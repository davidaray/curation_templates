#!/bin/bash 
#SBATCH --job-name=<TAXON>_TEcurate
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=50G

#### usage: sbatch TEcurate.sh
# Will analyze input from processed input from various programs to analyze and 
# categorize putative TEs. Output tables in $WORKDIR/prioritize. Additional output 
# files in $WORKDIR/te-aid. 
## Required prior steps:
# 1. Submit genome assembly to RepeatModeler analysis --> generate -families.fa file. 
#     Run rmodel.py using rmodel.sh. Both are available in my github repository.
#     Post-processing of this file should be to remove all identifying Class/Family info.
#     cut -d'#' -f1 aKon-families.fa >aKon-families_mod.fa
# 2. Submit output from .classified file to RepeatAfterMe (RAM) analysis using 
#     template_extend_align.sh 
#     (https://github.com/davidaray/bioinfo_tools/blob/master/template_extend_align.sh). 
#     Output from this should be in the $EXTENSIONSDIR defined below.
# 3. Collapse RAM output via cd-hit-est with our parameters. 
#     Example cd-hit-est run: cd-hit-est -i ID_newlib_4cdhit.fa -o ID_newlib.cdhit90_sS09 -c 0.90 -d 70 -aS 0.9 -n 9 -M 2200 -l 100
#     Process ID_newlib.cdhit90_sS09.clstr and ID_newlib.cdhit90_sS09.fa as necessary to
#     generate a final set of putative consensus sequences --> aKon_extended_rep.fa
# 4. Submit final set of putative consensus sequences to a new RepeatClassifier analysis using
#     repeatclassifier.sh 
#     (https://github.com/davidaray/bioinfo_tools/blob/master/repeatclassifer.sh)
#     Output from this analysis should be in $WORKDIR/repeatclassifier
#     Output should be called aKon_extended_rep.fa
## Required conda environments for this script and for steps above:
# /home/daray/extend_env.env.txt
# /home/daray/repeatmodeler.env.txt
# /home/daray/curate.env.txt
# Recreate on HPCC at TTU using: conda create --name [name of environment] --[environment file name]
## Must define all paths and NAME in the PATHS block below

##### To do?:
# 1. Incorporate tandem repeat finder to identify tandemly repeated Penelopes and LTR elements
#     common in some genomes. Need to avoid simple repeats but catch longer ones.
# 2. Incorporate the LTR identification package from Jessica Storer that I used with I. scapularis.
#     It will automatically find most LTRs and separate them from the internal portions, generating
#     a file ready to submit to Dfam.   


module --ignore-cache load gcc/10.1.0 r/4.0.2
. ~/miniforge3/etc/profile.d/conda.sh
conda activate curate

##### PATHS block
NAME=<TAXON>
MINORF=500
MINCOPY=10
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/hite/${NAME}
TARGET=${NAME}_confident_TE.cons.trimmed.fa
EXTENSIONSDIR=$WORKDIR/extensions

PFAM=/lustre/work/daray/software/pfam_db
AIDPATH=/lustre/work/daray/software/TE-Aid
ASSEMBLIESDIR=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies
AIDOUT=$WORKDIR/te-aid
mkdir -p $AIDOUT
mkdir -p $WORKDIR/prioritize
cd $WORKDIR/prioritize

#<<COMMENT
##Extract headers and subdivide names for later concatenation.
##Some hits are missing /Family. Correct for that.
echo -e "Extract headers and subdivide names for later concatenation.\n"
echo -e "Some hits are missing /Family. Correct for that."
grep ">" $WORKDIR/$TARGET | sed "s/#/-#/g" | sed "s/>//g" | sed "s|#Unknown|#Unknown/Unknown|g" | sed "s|#Satellite|#Satellite/Satellite|g" | sed "s|#LTR |#LTR/Unknown|g" | sed "s|#DNA |#DNA/Unknown|g" | sed "s|#tRNA |#tRNA/Nothing|g" | sed "s|#LINE |#LINE/Unknown|g" | sed "s|#|\\t|g" | sed "s|/|\\t|g" >${NAME}_name_class_family.txt
grep ">" $WORKDIR/$TARGET | sed "s/#/-#/g" | sed "s/>//g" | cut -d"#" -f1 >${NAME}_name.txt
grep ">" $WORKDIR/$TARGET | sed "s/#/-#/g" | sed "s/>//g" >${NAME}_original_headers.txt
DATE=$(date)
echo -e "Complete $DATE\n" 

echo -e "Concatenate table.txt."
paste ${NAME}_original_headers.txt ${NAME}_name_class_family.txt > table.txt
echo -e "Complete.\n"

##Get open reading frames for later blastp search
echo -e "Get open reading frames for later blastp search."
if [ ! -f ${NAME}_extended_rep_getorf.fa ]
	then 
	getorf $WORKDIR/$TARGET ${TARGET}_getorf.fa -minsize $MINORF
fi
DATE=$(date)
echo -e "Complete $DATE\n" 
#COMMENT

#<<COMMENT
##Test for presence of TE peptide library. Download if necessary and run blastp.
echo -e "Test for presence of TE peptide library. Download if necessary and run blastp."
if [ -e db/RepeatPeps.lib.phr ] 
  then 
    blastp -query ${TARGET}_getorf.fa  -db db/RepeatPeps.lib -outfmt 6 -evalue 1e-15 | sort -k1,1 -k12,12nr | sort -u -k1,1 | sed 's/#/--/g' > ${TARGET}_rep_blastp.out
  else
    mkdir -p $WORKDIR/prioritize/db
    cd db 
    wget https://raw.githubusercontent.com/rmhubley/RepeatMasker/master/Libraries/RepeatPeps.lib
    makeblastdb -in RepeatPeps.lib -out RepeatPeps.lib -dbtype prot &>/dev/null
    cd ..
	blastp -query ${TARGET}_getorf.fa  -num_threads 10 -db db/RepeatPeps.lib -outfmt 6 -evalue 1e-15 | sort -k1,1 -k12,12nr | sort -u -k1,1 | sed 's/#/--/g' > ${TARGET}_rep_blastp.out
fi
DATE=$(date)
echo -e "Complete $DATE\n" 
#COMMENT

#<<COMMENT
##Pull results of blastp, sort by longest hit, convert to columns, and add to growing list for concatenation.
echo -e "Pull results of blastp, sort by longest hit, convert to columns, and add to growing list for concatenation."
if [ -e ${NAME}_typelist.txt ] 
  then rm ${NAME}_typelist.txt
fi
while read -r I; do
  grep "$I" ${TARGET}_rep_blastp.out | cut -d$'\t' -f2,4 | sed "s|--|#|g" | cut -d"#" -f2,3 >rows.tmp
  COUNT=$(wc -l rows.tmp | cut -d" " -f1)
  if (( COUNT == 0 ))
    then 
      echo $COUNT "NOHIT" >>${NAME}_typelist.txt
    else 
      sort -n -k2 -r -o rows.tmp rows.tmp
      uniq rows.tmp >tetype.tmp
      TETYPES=$(tr '\n' ' ' < tetype.tmp)
      echo $COUNT $TETYPES >>${NAME}_typelist.txt
    rm tetype.tmp
  fi
done < ${NAME}_name.txt 
sed -i 's/  */\t/g' ${NAME}_typelist.txt
##Remove temporary files
rm tetype.tmp
rm rows.tmp
DATE=$(date)
echo -e "Complete $DATE\n" 
#COMMENT


#<<COMMENT
##Get TE consensus sequence lengths for later concatenation
echo -e "Get TE consensus sequence lengths for later concatenation."
seqkit fx2tab --length --name --header-line $WORKDIR/$TARGET | cut -d$'\t' -f2 >${NAME}_sizes.txt
sed -i '1d' ${NAME}_sizes.txt
DATE=$(date)
echo -e "Complete $DATE\n" 

##Build the initial table with results.
echo -e "Build the final table with results."
dos2unix ${NAME}_original_headers.txt ${NAME}_name_class_family.txt ${NAME}_sizes.txt ${NAME}_typelist.txt 
paste ${NAME}_original_headers.txt ${NAME}_name_class_family.txt ${NAME}_sizes.txt ${NAME}_typelist.txt > ${NAME}_table.txt
echo -e "Complete.\n"
#COMMENT

##Create sorting directories
echo -e "Create sorting directories."
TELIST="LINE SINE LTR RC DNA NOHIT"
for TENAME in $TELIST; do 
	mkdir $AIDOUT/$TENAME
done
DATE=$(date)
echo -e "Complete $DATE\n" 
#COMMENT

#Split consensus file to simulate having run extend_align
# Create splits directory
mkdir -p $WORKDIR/prioritize/splits
# Copy the target file to the splits directory
cp $WORKDIR/$TARGET $WORKDIR/prioritize/splits/$TARGET
# Replace '#' with '-#'
sed -i "s/#/-#/g" $WORKDIR/prioritize/splits/$TARGET
# Extract the first part of each line up to '#'
cut -d"#" -f1 $WORKDIR/prioritize/splits/$TARGET > $WORKDIR/prioritize/splits/${TARGET}.tmp
# Split the file using seqkit
seqkit split -i -O $WORKDIR/prioritize/splits $WORKDIR/prioritize/splits/${TARGET}.tmp
# Function to rename files
RENAME_FILES() {
  FILE=$1
  RENAME=$(echo $FILE | sed "s/confident_TE.cons.trimmed.fa.id_//g" | sed "s/.tmp/_rep.fa/g")
  mv $FILE $RENAME
}
export -f RENAME_FILES
# Use parallel to rename the files with up to 10 jobs
find $WORKDIR/prioritize/splits -name "*tmp" | parallel --jobs 10 RENAME_FILES
# Remove temporary files
rm $WORKDIR/prioritize/splits/${TARGET}.tmp $WORKDIR/prioritize/splits/${TARGET}_rep.fa

#<<COMMENT
##Copy created files to sorting directories.
echo -e "Copy created files to sorting directories."
TELIST="LINE SINE LTR RC DNA NOHIT"
for TENAME in $TELIST; do 
	awk '{print $2 "\t" $7}' ${NAME}_table.txt | sed "s|/|\t|g" >${NAME}_table.tmp 
	awk -v TENAME="$TENAME" '$2 == TENAME' ${NAME}_table.tmp  > ${NAME}_${TENAME}s.txt
	cut -d' ' -f1 ${NAME}_${TENAME}s.txt >${NAME}_${TENAME}s.tmp
	while read -r I; do
		CONSNAME=$(echo $I | awk '{print $1}')
		cp $WORKDIR/prioritize/splits/${NAME}_${CONSNAME}_rep.fa $AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa
	done < ${NAME}_${TENAME}s.tmp
done
rm ${NAME}_*.tmp
DATE=$(date)
echo -e "Complete $DATE\n" 
#COMMENT

##Change the header to shortened version
#echo -e "Change the header to shortened version."
#TELIST="LINE SINE LTR RC DNA"	
#for TENAME in $TELIST; do 
#	while read -r I; do
#		CONSNAME=$(echo $I | awk '{print $1}')
#		CONSNAMEMOD=${CONSNAME/-rnd-/.}
#		CONSNAMEMOD=${CONSNAMEMOD/_family-/.}
#		CLASS=$(echo $I | awk '{print $2}') 
#		FAMILY=$(echo $I | awk '{print $3}') 
#		HEADER=${CONSNAMEMOD}#${CLASS}/${FAMILY}
#		sed "s|${CONSNAME::-1}|$HEADER|g" $AIDOUT/$TENAME/${CONSNAME}_rep.fa > $AIDOUT/$TENAME/${CONSNAME}_rep_mod.fa
#	done < ${NAME}_${TENAME}s.txt 
#done
#while read -r I; do
#	CONSNAME=$(echo $I | awk '{print $1}')
#	CONSNAMEMOD=${CONSNAME/-rnd-/.}
#	CONSNAMEMOD=${CONSNAMEMOD/_family-/.}
#	HEADER=${CONSNAMEMOD}#Unknown/Unknown
#	sed "s|${CONSNAME::-1}|$HEADER|g" $AIDOUT/NOHIT/${CONSNAME}_rep.fa >$AIDOUT/NOHIT/${CONSNAME}_rep_mod.fa
#done < ${NAME}_NOHITs.txt
#DATE=$(date)
#echo -e "Complete $DATE\n" 

# Check orientation of ORF-containing hits and reverse complement if necessary
echo -e "Check orientation of ORF-containing hits and reverse complement if necessary."
TELIST="LINE SINE LTR RC DNA"
for TENAME in $TELIST; do
	mkdir -p $AIDOUT/check_orientation/$TENAME
done
# If no_blastx_hit.txt exists, erase it and start over.
if [ -f no_blastx_hit.txt ]; then
	rm no_blastx_hit.txt
fi
# Function to process each entry
PROCESS_ENTRY() {
    TENAME=$1
    I=$2
    AIDOUT=$3
    NAME=$4
    CONSNAME=$(echo $I | awk '{print $1}')
    FILE="$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa"
    echo "Analyzing $FILE"
    blastx -query "$FILE" -db db/RepeatPeps.lib -outfmt 6 -evalue 1e-15 | sort -k1,1 -k12,12nr | sort -u -k1,1 | sed 's/#/--/g' > "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_extended_rep_blastx.out"

    if [[ -s "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_extended_rep_blastx.out" ]]; then
        START=$(head -1 "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_extended_rep_blastx.out" | awk '{print $7}')
        echo "start = $START"
        END=$(head -1 "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_extended_rep_blastx.out" | awk '{print $8}')
        echo "end = $END"
        if (( START > END )); then
            echo "start > end. Hit is reversed."
            echo -e "Reverse complementing $FILE\n"
            # Reverse complement rep file
            seqkit seq -r -p -t DNA "$FILE" > "$AIDOUT/$TENAME/${NAME}_${CONSNAME}-rep_rc.tmp"
            mv "$AIDOUT/$TENAME/${NAME}_${CONSNAME}-rep_rc.tmp" "$FILE"
            # Reverse complement rep_mod file
#            seqkit seq -r -p -t DNA "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep_mod.fa" > "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep_mod.tmp"
#            mv "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep_mod.tmp" "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep_mod.fa"
        else
            echo "start < end. Hit is in correct orientation."
            echo -e "Not reverse complementing $FILE\n"
        fi
    else
        echo -e "No blastx hits for $FILE\n"
        echo "$FILE" >> "${NAME}_no_blastx_hit_${TENAME}.txt"
        mv "$AIDOUT/$TENAME/${NAME}_${CONSNAME}_extended_rep_blastx.out" "$AIDOUT/check_orientation/$TENAME"
    fi
}
export -f PROCESS_ENTRY  # Export the function to be used by parallel
# Process each TENAME in parallel
for TENAME in $TELIST; do
    cat "${NAME}_${TENAME}s.txt" | parallel --jobs 10 PROCESS_ENTRY $TENAME {} $AIDOUT $NAME
done
DATE=$(date)
echo -e "Complete $DATE\n"

##Run TE-Aid on files in each category
echo -e "Run TE-Aid on files in each category."
TELIST="LINE SINE LTR RC DNA NOHIT"	
if [ -f ${NAME}_final_table.txt ] 
	then rm ${NAME}_final_table.txt
fi
#printf "%-45s \t %-30s \t %-10s \t %-10s \t %-17s \t %-20s \t %-8s \t %-10s \t %-14s \t %-10s \t %-14s \t %-10s \t %-14s \t %-10s \t %-14s\n" "RM_ID" "Short_ID"  "Class" "Family" "Consensus_length" "90percent_consensus" "N_ORFS" "ORF1_type" "ORF1_length" "ORF2_type" "ORF2_length"	 "ORF3_type" "ORF3_length" >${NAME}_final_table.txt
echo -e "RM_ID \t Short_ID \t Class \t Family \t Consensus_length \t 90percent_consensus \t N_ORFS \t ORF1_type \t ORF1_length \t ORF2_type \t ORF2_length \t ORF3_type \t ORF3_length" >${NAME}_final_table.txt

mkdir ../assembly
gunzip -c $ASSEMBLIESDIR/${NAME}.fa.gz >../assembly/${NAME}.fa

for TENAME in $TELIST; do 
	echo "TE-Aid processing of files in "$NAMESFILE
	while read -r I; do
		CONSNAME=$(echo $I | awk '{print $1}')
		FILE=$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa
		echo "TE-Aid processing "$FILE
		#Generate reverse complement files for identifying TIRs
		echo "Generate reverse complement files for identifying TIRs."
		seqkit seq $FILE -r -p -t DNA >$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep_rc.fa
		cat $FILE $AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep_rc.fa >$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep_RC.fa
		rm $AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep_rc.fa
		#Run TE-Aid
		echo "Run TE-Aid"
		$AIDPATH/TE-Aid -q $FILE -g ../assembly/${NAME}.fa -T -o $AIDOUT/$TENAME
		mv ${FILE}.c2g.pdf $AIDOUT/$TENAME/${NAME}_${CONSNAME}.c2g.pdf
		#Gather information for final table
		tail -n +2 $AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa.genome.blastn.out | awk '{print $1 "\t" $7 "\t" $8}' | awk 'BEGIN { OFS = "\t" } { $4 = $3 - $2 + 1 } 1' > $AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa.genome.blastn.tmp
		CONSSIZE=$(seqkit fx2tab --length --name $FILE | awk '{print $2}')
		MINCONSSIZE=$(awk "BEGIN { print $CONSSIZE * 0.9 }")
		BLASTTMP=$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa.genome.blastn.tmp
		FULLCOUNT=$(awk -v MINCONSSIZE="$MINCONSSIZE" -v BLASTOUT="$BLASTTMP" '$4 > MINCONSSIZE' $BLASTTMP | wc -l)
		ROW=$(grep $CONSNAME ${NAME}_table.txt | awk -v FULLCOUNT="$FULLCOUNT" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" FULLCOUNT "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}')
		#Generate final table
		echo $ROW >> ${NAME}_final_table.txt
	done < ${NAME}_${TENAME}s.txt
done	
TELIST="LINE SINE LTR RC DNA"	
for TENAME in $TELIST; do 
	echo "TE-Aid processing of files in "$NAMESFILE
	while read -r I; do
		CONSNAME=$(echo $I | awk '{print $1}')
		FILE=$AIDOUT/check_orientation/$TENAME/${NAME}_${CONSNAME}_rep.fa
		echo "TE-Aid processing "$FILE
		#Generate reverse complement files for identifying TIRs
		seqkit seq $FILE -r -p -t DNA >$AIDOUT/check_orientation/$TENAME/${NAME}_${CONSNAME}_rep_rc.fa
		cat $FILE $AIDOUT/check_orientation/$TENAME/${NAME}_${CONSNAME}_rep_rc.fa >$AIDOUT/check_orientation/$TENAME/${NAME}_${CONSNAME}_rep_RC.fa
		rm $AIDOUT/check_orientation/$TENAME/${NAME}_${CONSNAME}_rep_rc.fa
		#Run TE-Aid
		$AIDPATH/TE-Aid -q $FILE -g ../assembly/${NAME}.fa -T -o $AIDOUT/check_orientation/$TENAME/
		mv ${FILE}.c2g.pdf $AIDOUT/check_orientation/$TENAME/${NAME}_${CONSNAME}.c2g.pdf
		#Gather information for final table
		tail -n +2 $AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa.genome.blastn.out | awk '{print $1 "\t" $7 "\t" $8}' | awk 'BEGIN { OFS = "\t" } { $4 = $3 - $2 + 1 } 1' > $AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa.genome.blastn.tmp
		CONSSIZE=$(seqkit fx2tab --length --name $FILE | awk '{print $2}')
		MINCONSSIZE=$(awk "BEGIN { print $CONSSIZE * 0.9 }")
		BLASTTMP=$AIDOUT/$TENAME/${NAME}_${CONSNAME}_rep.fa.genome.blastn.tmp
		FULLCOUNT=$(awk -v MINCONSSIZE="$MINCONSSIZE" -v BLASTOUT="$BLASTTMP" '$4 > MINCONSSIZE' $BLASTTMP | wc -l)
		ROW=$(grep $CONSNAME ${NAME}_table.txt | awk -v FULLCOUNT="$FULLCOUNT" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" FULLCOUNT "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}')
		#Generate final table
		echo $ROW >> ${NAME}_final_table.txt
	done < ${NAME}_no_blastx_hit_${TENAME}.txt
done	
rm -rf ../assembly
DATE=$(date)
echo -e "Complete $DATE\n" 

# Loop through all *_LTR_*_LTR-rep.fa files in the source directory
mkdir $AIDOUT/hite_LTRs/
# Define the directories to search for files to remove
SEARCH_DIRS=("$AIDOUT/LTR" "$AIDOUT/DNA" "$AIDOUT/LINE" "$AIDOUT/NOHIT")  
for LTR_FILE in $WORKDIR/prioritize/splits/*_LTR_*_LTR-_rep.fa; do
	# Extract the base name without the suffix
	BASE_NAME=$(basename "$LTR_FILE" LTR-_rep.fa)
	echo "Searcing for "$BASE_NAME" files"
	# Construct the corresponding *_LTR_*_INT-rep.fa file name
	INT_FILE="$WORKDIR/prioritize/splits/${BASE_NAME}INT-_rep.fa"

	# Check if the *_LTR_0_INT-rep.fa file exists
	if [ -e "$INT_FILE" ]; then
		echo $INT_FILE" exists"
		# Copy both files to the destination directory
		cp "$LTR_FILE" "$AIDOUT/hite_LTRs/"
		cp "$INT_FILE" "$AIDOUT/hite_LTRs/"
		for DIR in "${SEARCH_DIRS[@]}"; do
			find "$DIR" -type f -name "*${BASE_NAME}*" -exec mv -f {} $AIDOUT/hite_LTRs/ \; -exec echo "Moved files associated with {} from $DIR to $AIDOUT/hite_LTRs/" \;
		done
		echo "Copied pair: $LTR_FILE and $INT_FILE to $AIDOUT/hite_LTRs/"
	fi
done

## Filter hits with fewer than 10 copies >90% of full-length
echo -e "Filtering hits for $NAME with fewer than 10 copies >90% of full-length"
mkdir -p "$AIDOUT/fewhits"
sed '1d' "${NAME}_final_table.txt" | awk '{print $2 "\t" $6}' | awk -v MINCOPY="$MINCOPY" '$2 < MINCOPY' >"${NAME}_filtered_for_low_count.txt"
# Remove all blank lines from the filter file.
sed -i ':a;/^\s*$/{$d;N;ba}' ${NAME}_filtered_for_low_count.txt
SEARCH_DIRS=("$AIDOUT/LTR" "$AIDOUT/DNA" "$AIDOUT/LINE" "$AIDOUT/RC" "$AIDOUT/NOHIT")
MAX_JOBS=10  # Maximum number of parallel jobs
JOB_COUNT=0
MOVE_FILES() {
    local CONSNAME="$1"
    for DIR in "${SEARCH_DIRS[@]}"; do
        if find "$DIR" -type f -name "*${CONSNAME}*" | grep -q .; then
            echo "Found files associated with ${CONSNAME}. Moving from $DIR to fewhits."
            find "$DIR" -type f -name "*${CONSNAME}*" -exec mv -f {} "$AIDOUT/fewhits/" \; 
        fi
    done
}
while read -r I; do
    CONSNAME=$(echo "$I" | awk '{print $1}')
    CONSNAME="${NAME}_${CONSNAME}"
    MOVE_FILES "$CONSNAME" &
    ((JOB_COUNT++))

    # Control the number of concurrent jobs
    if (( JOB_COUNT >= MAX_JOBS )); then
        wait -n
        ((JOB_COUNT--))
    fi
done < "${NAME}_filtered_for_low_count.txt"
# Wait for all background jobs to finish
wait
DATE=$(date)
echo -e "Complete $DATE\n"


##Filter hite_LTR hits with fewer than 10 copies >90% of full-length
echo -e "Filtering hite_LTR hits for $NAME with fewer than 10 copies >90% of full-length"
mkdir -p "$AIDOUT/hite_LTRs/fewhits"
while read -r I; do
    CONSNAME=$(echo "$I" | awk '{print $1}')
    CONSNAME="${NAME}_${CONSNAME}"
    if find "$AIDOUT/hite_LTRs" -type f -name "*${CONSNAME}*" | grep -q .; then
		echo "Found files associated with ${CONSNAME}. Moving from hite_LTRs to hite_LTRs/fewhits."
		mv "$AIDOUT/hite_LTRs/${CONSNAME}"* "$AIDOUT/hite_LTRs/fewhits/"
	fi
    ((JOB_COUNT++))
    # Control the number of concurrent jobs
    if (( JOB_COUNT >= MAX_JOBS )); then
        wait -n
        ((JOB_COUNT--))
    fi
done < "${NAME}_filtered_for_low_count.txt"
# Wait for all background jobs to finish
wait
DATE=$(date)
echo -e "Complete $DATE\n"


##Prepare files for download and manual inspection as necessary
echo "Prepare files for download and manual inspection as necessary"
TELIST="LINE SINE LTR RC DNA NOHIT"	
for TENAME in $TELIST; do 
	mkdir -p $WORKDIR/fordownload/$TENAME
	cp $AIDOUT/$TENAME/*.pdf $AIDOUT/$TENAME/*.fa $AIDOUT/$TENAME/*self-blast.pairs.txt $WORKDIR/fordownload/$TENAME
	tar -zcf $WORKDIR/fordownload/fordownload_${TENAME}.tgz $WORKDIR/fordownload/$TENAME
done
cp -r $AIDOUT/check_orientation $WORKDIR/fordownload/check_orientation
tar -zcf $WORKDIR/fordownload/fordownload_check_orientation.tgz $WORKDIR/fordownload/check_orientation
cp -r $AIDOUT/zerohits $WORKDIR/fordownload/zerohits	
tar -zcf $WORKDIR/fordownload/fordownload_zerohits.tgz $WORKDIR/fordownload/zerohits
cp -r $AIDOUT/fewhits $WORKDIR/fordownload/fewhits
tar -zcf $WORKDIR/fordownload/fordownload_fewhits.tgz $WORKDIR/fordownload/fewhits
DATE=$(date)
echo -e "Complete $DATE\n" 


