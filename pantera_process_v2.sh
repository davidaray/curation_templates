#!/bin/sh
#SBATCH --job-name=pre_pantera
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=12

. ~/miniforge3/etc/profile.d/conda.sh
conda activate samtools

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/pantera2
CHRPATH=$WORKDIR/chromosomes
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates 
ARCHIVE=$WORKDIR/archive
TABLES=$WORKDIR/tables
RUNS=$WORKDIR/runs
# Define the path to the Singularity image
SINGULARITY_IMAGE=$WORKDIR/pggb_container/pggb_latest.sif
BIN=$WORKDIR/bin
PAIRS=/lustre/scratch/daray/bat1k_TE_analyses/species_pairs.txt

# Check if $ARCHIVE exists, then create it if it doesn't
if [ ! -d "$ARCHIVE" ]; then
    echo "Creating directory: $ARCHIVE"
    mkdir -p "$ARCHIVE"
else
    echo "Directory $ARCHIVE already exists."
fi
# Check if $TABLES exists, then create it if it doesn't
if [ ! -d "$TABLES" ]; then
    echo "Creating directory: $TABLES"
    mkdir -p "$TABLES"
else
    echo "Directory $TABLES already exists."
fi
# Check if $RUNS exists, then create it if it doesn't
if [ ! -d "$RUNS" ]; then
    echo "Creating directory: $RUNS"
    mkdir -p "$RUNS"
else
    echo "Directory $RUNS already exists."
fi
# Check if $BIN exists, then create it if it doesn't
if [ ! -d "$BIN" ]; then
    echo "Creating directory: $BIN"
    mkdir -p "$BIN"
else
    echo "Directory $BIN already exists."
fi

## Generate equivalency tables for chromosome tracking and process with bgzip for downstream analyses

## This script assumes pantera_prep_v2.sh has been run to generate paired chromosome files. 

cd $WORKDIR

cat $PAIRS | while read i; do
    cd $WORKDIR
	ID=$(echo $i | cut -d"." -f1)
	mkdir -p $RUNS/$ID/pggb
    INC=1
    echo -e "scaffold_names\tpantera_ID" > $TABLES/${ID}_scaffold_table.tsv

    for FILE in "$CHRPATH/${ID}"*; do
        # Replace spaces in fasta headers with underscores
        sed -i 's/^>.* /&_/; s/^>.* /_/g' "$FILE"

        # Add "A" to duplicate headers
        awk '/^>/ {header=$0; if (seen[header]++) $0=header "A"} 1' "$FILE" > "$CHRPATH/tmp.fa"
        mv "$CHRPATH/tmp.fa" "$FILE"
        FILENAME=$(basename "$FILE" .fa)
        CHROMID="${ID}_chr${INC}"

        # Update scaffold table
        echo -e "${FILENAME}\t${CHROMID}" >> "$TABLES/${ID}_scaffold_table.tsv"

        # Copy and compress with a unique name
        cp "$FILE" "$CHRPATH/${CHROMID}.fa"
        bgzip -@ 10 "$CHRPATH/${CHROMID}.fa" -c > "$CHRPATH/${CHROMID}.fa.bz"
        samtools faidx "$CHRPATH/${CHROMID}.fa.bz"
		rm "$CHRPATH/${CHROMID}.fa"

        # Move the original file to the archive
        mv "$FILE" "$ARCHIVE"
		# Generate pggb submission script
		sed "s/<NAME>/$CHROMID/g" $CURATIONPATH/pggb_template.sh >$BIN/${CHROMID}_pggb.sh
		sed -i "s|PGGB=|PGGB=$RUNS/$ID/pggb|g" $BIN/${CHROMID}_pggb.sh
		chmod +x $BIN/${CHROMID}_pggb.sh
		
		# Generate pantera submission script
		sed "s/<NAME>/$CHROMID/g" $CURATIONPATH/pantera_template_v1.sh >$BIN/${CHROMID}_pantera.sh
		sed -i "s|PGGB=|PGGB=$RUNS/$ID/pggb|g" $BIN/${CHROMID}_pantera.sh
		sed -i "s|PANTERAOUT=|PANTERAOUT=$RUNS/$ID/pantera|g" $BIN/${CHROMID}_pantera.sh
		chmod +x $BIN/${CHROMID}_pantera.sh
		
		# Submit pggb script
		cd $BIN
		sbatch ${CHROMID}_pggb.sh

        # Increment counter
        INC=$((INC + 1))
    done
done

