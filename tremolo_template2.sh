#!/bin/sh
#SBATCH --job-name=<NAME>_TrEMOLO
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=10:00:00

. ~/miniforge3/etc/profile.d/conda.sh
conda activate 

# Paths
TREMOLODIR=/lustre/work/daray/software/TrEMOLO
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/tremolo2
CONTAINER=$WORKDIR/container
PAIRSLIST=/lustre/scratch/daray/bat1k_TE_analyses/species_pairs.txt
BAT1K=/lustre/scratch/daray/bat1k_TE_analyses
LIBRARY=/lustre/scratch/daray/bat1k_TE_analyses/bat1k_named_20250224.fa
OUTPUT=$WORKDIR/output
DATAFREEZE=/lustre/scratch/daray/bat1kdatafreeze
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies
CURATIONDIR=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates

# Species being examined
SPECIES=<NAME>

# Create output directory
mkdir -p $OUTPUT/$SPECIES
cd $OUTPUT/$SPECIES

#: << 'EOF'

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
gunzip -c $GZPATH/${REF}.fa.gz >$WORKDIR/unzipped_assemblies/${REF}.fa
gunzip -c $GZPATH/${GEN}.fa.gz >$WORKDIR/unzipped_assemblies/${GEN}.fa


# Create config file for this species
sed "s|/path/to/reference_file.fasta|$WORKDIR/unzipped_assemblies/${REF}.fa|g" $TREMOLODIR/config_INSIDER.yaml >$WORKDIR/bin/${SPECIES}_insider_config.yaml
sed -i "s|/path/to/genome_file.fasta|$WORKDIR/unzipped_assemblies/${GEN}.fa|g" $WORKDIR/bin/${SPECIES}_insider_config.yaml
sed -i "s|TrEMOLO_OUTPUT|$OUTPUT/$SPECIES|g" $WORKDIR/bin/${SPECIES}_insider_config.yaml
sed -i "s|/path/to/database_TE.fasta|$LIBRARY|g" $WORKDIR/bin/${SPECIES}_insider_config.yaml
sed -i "s|THREADS: 8|THREADS: 12|g" $WORKDIR/bin/${SPECIES}_insider_config.yaml

# Run TrEMOLO
singularity exec \
	-B $TREMOLODIR:$TREMOLODIR,$WORKDIR:$WORKDIR,$BAT1K:$BAT1K,$DATAFREEZE:$DATAFREEZE \
	$CONTAINER/TrEMOLO.simg snakemake \
	--snakefile $TREMOLODIR/run.snk \
	--configfile $WORKDIR/bin/${SPECIES}_insider_config.yaml
	
#Run Module Analysis 
# https://github.com/DrosophilaGenomeEvolution/TrEMOLO/blob/master/modules/2-MODULE_TE_BLAST/README.md
#singularity exec \
#	-B $TREMOLODIR:$TREMOLODIR,$WORKDIR:$WORKDIR,$BAT1K:$BAT1K,$DATAFREEZE:$DATAFREEZE \
#	$CONTAINER/TrEMOLO.simg $TREMOLODIR/modules/2-MODULE_TE_BLAST/server/scripts/buildData.sh $OUTPUT/$SPECIES

# Clean up
rm $WORKDIR/unzipped_assemblies/${REF}.fa 
rm $WORKDIR/unzipped_assemblies/${REF}.fa.fai
rm $WORKDIR/unzipped_assemblies/${GEN}.fa 
rm $WORKDIR/unzipped_assemblies/${GEN}.fa.fai
#EOF

# Process the TE_INFOS file to create a more readable report
python $CURATIONDIR/process_tremolo_bed.py -i $OUTPUT/$SPECIES/TE_INFOS.bed -o $OUTPUT/$SPECIES/REPORT/${SPECIES}_TE_INFOS_mod.bed
mkdir -p $OUTPUT/${SPECIES}-report/
cp -r $OUTPUT/${SPECIES}/REPORT/* $OUTPUT/${SPECIES}-report/
mv $OUTPUT/${SPECIES}-report/report.html $OUTPUT/${SPECIES}-report/${SPECIES}-report.html

# Process DELETION_TE.bed file to find polymorphisms on the alternate haplotype and concatenate with the TE_INFOS.bed file
sed "s/<GENOME>/$SPECIES/g" $CURATIONDIR/tremolo_deletion_tsd_process.sh >$WORKDIR/bin/${SPECIES}_tsd_process.sh 
sh $WORKDIR/bin/${SPECIES}_tsd_process.sh

# Correct for a few issues with the formatting
sed -E "s|#chrom\tstart\tend\tTE_annotation\tAssemblytics_value\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\tNEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE|#chrom\tstart\tend\tTE_annotation\tCLASS/FAMILY\tAssemblytics_value\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\tNEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE|g" "$OUTPUT/${SPECIES}-report/${SPECIES}_all_TEs_info.tsv" >$OUTPUT/${SPECIES}-report/tmp.tsv
sed -e "2,\$s/#/\t/g" "$OUTPUT/${SPECIES}-report/tmp.tsv" > "$OUTPUT/${SPECIES}-report/${SPECIES}_all_TEs_info.tsv"
rm "$OUTPUT/${SPECIES}-report/tmp.tsv"
#EOF

# Process random loci to determine if polymorphisms are apparent.
# '5' indicates the number of alignments to generate.
sh $CURATIONDIR/process_tremolo_loci_v2.sh \
	$GZPATH/${SPECIES}.hap1.fa.gz \
	$GZPATH/${SPECIES}.hap2.fa.gz \
	$OUTPUT/$SPECIES/INSIDER/TE_DETECTION/DELETION_TE.bed \
	$OUTPUT/$SPECIES/INSIDER/TE_DETECTION/DELETION_TE_ON_REF.bed \
	$OUTPUT/$SPECIES/polymorphism_check \
	5


