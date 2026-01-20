#!/bin/bash 
#SBATCH --job-name=<NAME>_TEcurate
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=50G

module --ignore-cache load gcc/10.1.0 r/4.0.2
. ~/miniforge3/etc/profile.d/conda.sh
conda activate curate

##### PATHS block
NAME=<NAME>
MINORF=500
MINCOPY=10
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/hite/${NAME}
TARGET=${NAME}_confident_TE.cons.trimmed.fa
EXTENSIONSDIR=$WORKDIR/extensions

PFAM=/lustre/work/daray/software/pfam_db
AIDPATH=/lustre/work/daray/software/TE-Aid
ASSEMBLIESDIR=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies
AIDOUT=$WORKDIR/te-aid
cd $WORKDIR/prioritize

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


