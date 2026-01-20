#!/bin/bash
# Input arguments
REFERENCE=$1  # Reference genome (Haplotype 1)
GENOME=$2  # Haplotype 2
BEDFILE=$3  # BED file with TE insertion coordinates on genome
BEDFILEONREF=$4  # BED file with TE insertion coordinates on reference
WORKDIR=$5  # Directory to save output
N=$6 # Number of random loci to interrogate

. ~/miniforge3/etc/profile.d/conda.sh
conda activate

# 1. Create output directory
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# 2. Clean up BED files
sed -e 's/#/\t/g' -e 's/|/\t/g' "${BEDFILE}" >"${BEDFILE}_processloci"
BEDFILE=${BEDFILE}_processloci
echo $BEDFILE
sed -e 's/#/\t/g' -e 's/|/\t/g' "${BEDFILEONREF}" >"${BEDFILEONREF}_processloci"
BEDFILEONREF=${BEDFILEONREF}_processloci
echo $BEDFILEONREF

# 3. Unzip genome assemblies
REFERENCENAME=$(basename "${REFERENCE%.gz}")
if [ ! -f "$REFERENCENAME" ]; then
    gunzip -c "$REFERENCE" > "$WORKDIR/$REFERENCENAME"
fi
GENOMENAME=$(basename "${GENOME%.gz}")
if [ ! -f "$GENOMENAME" ]; then
    gunzip -c "$GENOME" > "$WORKDIR/$GENOMENAME"
fi

# 4. Select N random loci
shuf -n "$N" "$BEDFILE" | awk '{OFS="\t"; print $0}' > "genome_selected_loci.tsv"

# 5. Buffer coordinates by 500bp
awk -v OFS="\t" '{start=($2-500<0)?0:$2-500; print $1, start, $3+500, $4, $5, $6, "buff_"$NF}' "genome_selected_loci.tsv" > "genome_buffered_loci.tsv"

# 6. Get corresponding loci in reference
> reference_selected_loci.tsv
while IFS=$'\t' read -r CHROM START END ID TE_TYPE ASS_ID LOCUS; do
MATCHING_LINE=$(awk -v val="$ASS_ID" -F'\t' '$6 == val' "$BEDFILEONREF")
if [[ -n "$MATCHING_LINE" ]]; then
echo "$MATCHING_LINE" >> reference_selected_loci.tsv
fi
done < genome_selected_loci.tsv

# 7. Buffer reference coordinates
awk -v OFS="\t" '{start=($2-500<0)?0:$2-500; print $1, start, $3+500, $4, $5, $6, "buff_"$NF}' "reference_selected_loci.tsv" > "reference_buffered_loci.tsv"

# 8. Extract sequences from genome
conda activate bedtools
while IFS=$'\t' read -r CHROM START END ID TE_TYPE ASS_ID; do
    echo -e "$CHROM\t$START\t$END\t$ASS_ID" > temp.bed
    bedtools getfasta -fi "$WORKDIR/$GENOMENAME" -bed temp.bed -fo "${ASS_ID}_genome.fa"
done < "genome_buffered_loci.tsv"
rm temp.bed

# 9. Extract sequences from reference
while IFS=$'\t' read -r CHROM START END ID TE_TYPE ASS_ID; do
    echo -e "$CHROM\t$START\t$END\t$ASS_ID" > temp.bed
	bedtools getfasta -fi "$WORKDIR/$REFERENCENAME" -bed temp.bed -fo "${ASS_ID}_reference.fa"
done < "reference_buffered_loci.tsv"
rm temp.bed

conda activate muscle
# 9. Align sequences using MUSCLE
while IFS=$'\t' read -r CHROM START END ID TE_TYPE ASS_ID; do
cat "${ASS_ID}_genome.fa" "${ASS_ID}_reference.fa" >"${ASS_ID}_combined.fa"
muscle -align "${ASS_ID}_combined.fa" -output "${ASS_ID}_aligned.fa"
done < "reference_buffered_loci.tsv"

# 10. Clean up files
rm $WORKDIR/$GENOMENAME $WORKDIR/${GENOMENAME}.fai
rm $WORKDIR/$REFERENCENAME $WORKDIR/${REFERENCENAME}.fai
rm "reference_selected_loci.tsv"
rm "genome_selected_loci.tsv"
rm *_combined.fa *_reference.fa *_genome.fa

echo "Processing completed. Outputs are in $WORKDIR"
