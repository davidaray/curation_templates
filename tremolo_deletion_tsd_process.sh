#!/bin/bash 

. ~/miniforge3/etc/profile.d/conda.sh
conda activate bedtools

# Species being examined
SPECIES=<GENOME>

# Paths
TREMOLODIR=/lustre/work/daray/software/TrEMOLO
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/tremolo2
CONTAINER=/lustre/scratch/daray/bat1k_TE_analyses/tremolo2/container
PAIRSLIST=/lustre/scratch/daray/bat1k_TE_analyses/species_pairs.txt
BAT1K=/lustre/scratch/daray/bat1k_TE_analyses
LIBRARY=/lustre/scratch/daray/bat1k_TE_analyses/bat1k_named_20241126.fa
OUTPUT=$WORKDIR/output
DATAFREEZE=/lustre/scratch/daray/bat1kdatafreeze
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies
CURATIONDIR=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates

# Navigate to output directory
cd $OUTPUT/$SPECIES

# Identify genome fasta names and library
PAIR=$(grep "$SPECIES" "$PAIRSLIST")
REF=$(echo "$PAIR" | cut -d" " -f1)
GEN=$(echo "$PAIR" | cut -d" " -f2)
LIB=$(basename "$LIBRARY")

# Output information
echo "Reference assembly = $REF"
echo "Genome to query = $GEN"
echo "TE library = $LIB"

# Unzip input assemblies
mkdir -p $WORKDIR/unzipped_assemblies
#gunzip -c $GZPATH/${REF}.fa.gz >$WORKDIR/unzipped_assemblies/${REF}.fa
gunzip -c $GZPATH/${GEN}.fa.gz >$WORKDIR/unzipped_assemblies/${GEN}.fa

SIZE_FLANK=15
GENOME=$WORKDIR/unzipped_assemblies/${GEN}.fa

awk -v sizeFlank="$SIZE_FLANK" '
            OFS="\t" {{
                print $1, $2-sizeFlank, $2, ($3-$2)":"$4":FK_L"; 
                print $1, $3, $3+sizeFlank, ($3-$2)":"$4":FK_R";
            }}' DELETION_TE.bed > FK_INS.bed

echo "FORMAT FK SEQ TE.."
bedtools getfasta -tab -fi $GENOME -bed FK_INS.bed -name+ | \
    awk -v size_flank="$SIZE_FLANK" '
        NR%2==1 {{
            split($1, sp, ":"); 
            info=sp[2];

            split(info, sp2, "|"); 
            TE=sp2[1];
            ID=sp2[2];

            seqL=$2;
            infos=sp[1]":"sp[5]":"sp[6]":"sp[2];
        }}
                
        #FK_R
        NR%2==0 && OFS="\t" {{
            print infos, TE, ID, substr(seqL, length(seqL)-size_flank, length(seqL)), substr($2, 1, size_flank), ".", ".";
        }}

        ' > FK_INS_FT.bed;

# Clean up
#rm $WORKDIR/unzipped_assemblies/${REF}.fa
rm $WORKDIR/unzipped_assemblies/${GEN}.fa

echo "GET TSD INSIDER..."
singularity exec \
	-B $TREMOLODIR:$TREMOLODIR,$WORKDIR:$WORKDIR,$BAT1K:$BAT1K,$DATAFREEZE:$DATAFREEZE \
	$CONTAINER/TrEMOLO.simg \
	python3 $TREMOLODIR/lib/python/TSD/getTSDgenome.py -t 12 \
    FK_INS_FT.bed \
    TSD_TE.tsv

# Remove existing output file before starting
if [[ -f "deletion_infos.bed" ]]; then
rm deletion_infos.bed
fi

# Process each line in TSD_TE.tsv
while read -r line; do
# Extract CHR, TSD, and COORD from the current line
CHR=$(echo "$line" | cut -d":" -f2)
TSD=$(echo "$line" | awk '{print $2}')
COORD=$(echo "$line" | cut -d":" -f3 | cut -d"-" -f2)
#COORD=$((COORD - 15))
#echo -e "$CHR\t$COORD\t$TSD"

# Find matching deletion line
DELETION_LINE=$(awk -v chr="$CHR" -v coord="$COORD" -F'\t' '$1 == chr && $2 == coord' DELETION_TE.bed)
#echo $DELETION_LINE

# Skip if no match is found
[[ -z "$DELETION_LINE" ]] && continue

# Extract fields from deletion line
read -r CHROM START STOP TE TYPE ID <<< "$DELETION_LINE"

# Append output to deletion_infos.bed
echo -e "$CHROM\t$START\t$STOP\t$TE\t$ID\t$TSD\t?\t?\t?\t?\tNONE\tINSIDER\t?\t?\tDELETION" >> deletion_infos.bed

done < TSD_TE.tsv

sed -i "s/|/\t/g" deletion_infos.bed

cat REPORT/${SPECIES}_TE_INFOS_mod.bed deletion_infos.bed >../${SPECIES}-report/${SPECIES}_all_TEs_info.tsv
	