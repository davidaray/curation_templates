#!/bin/bash 
#SBATCH --job-name=<NAME>_postRM    # Uses the genome name for job name
#SBATCH --output=%x.%j.out             # Job name and job ID in output file
#SBATCH --error=%x.%j.err              # Job name and job ID in error file
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=12

# Activate the required conda environment
. ~/miniforge3/etc/profile.d/conda.sh
conda activate plotting

# Set necessary variables
GENOME=<NAME>
RUNTYPE=${GENOME}_RM
DIR=/lustre/scratch/daray/bat1k_TE_analyses/repeatmasker/$RUNTYPE
GITPATH=/home/daray/gitrepositories/bioinfo_tools/
CURATIONDIR=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates

# Change to the working directory
cd $DIR

# Check if the RepeatMasker is still running by looking for files containing 'machine'
if ls ${GENOME}_RM/*machine* >/dev/null 2>&1 && [ -d "${GENOME}_RM/RMPart" ]; then
    echo "Files with 'machine' in the filename found for $GENOME. RMPart folder is present. Still running RepeatMasker. Exiting."
    exit 1 
fi

# If no 'machine' files are found but RMPart folder is present, submit doLift.sh
if ! ls ${GENOME}_RM/*machine* >/dev/null 2>&1 && [ -d "${GENOME}_RM/RMPart" ]; then
    echo "No files with 'machine' in the filename found for $GENOME, but 'RMPart' folder is present. Submitting doLift.sh and exiting."
    sbatch doLift.sh
    exit 1 
fi

# Check if doLift has been run and if the RMPart directory still exists
if [[ ! -f ${GENOME}.fa.out.gz ]] && [ -d "${GENOME}_RM/RMPart" ]; then
    echo "Run doLift for $GENOME"
    sbatch doLift.sh
    exit 1 
else
    echo "doLift is finished for $GENOME. Removing leftover files and continuing."
    rm -rf RMPart
    rm ${GENOME}.fa ${GENOME}.2bit *.err *.out
fi

# Run RM2bed if the .bed file does not exist. 
if [[ ! -f ${GENOME}_rm.bed ]]; then
    echo "Running rm2bed for $GENOME"
    python $GITPATH/RM2bed.py -d . -sp class -p ${GENOME} -o lower_divergence ${GENOME}.fa.align.gz
else    
    echo "rm2bed is finished for $GENOME. Continuing."
fi

# Generate plots if they haven't been created
if [[ ! -f ${GENOME}.plot_out_pie.png ]]; then
    # Generate stacked bar and pie plots including all TE classes
    echo "Generating stacked bar and pie for $GENOME, all categories."
    python $CURATIONDIR/te_plot_stacked_pie.py \
        -r ${GENOME}.fa.out.gz \
        -s ${GENOME}.summary.gz \
        -o ${GENOME}.plot \
        -b ${GENOME}_rm.bed \
        --plot_type both \
        -proc 12 \
        -bin 1 

    # Generate stacked bar and pie plots excluding Satellite TE class
    echo "Generating stacked bar and pie for $GENOME, selected categories."
    python $CURATIONDIR/te_plot_stacked_pie.py \
        -r ${GENOME}.fa.out.gz \
        -s ${GENOME}.summary.gz \
        -o ${GENOME}.plot.nosat \
        -b ${GENOME}_rm.bed \
        --plot_type both \
        -proc 12 \
        -bin 1 \
        --classes DNA,DIRS,LINE,LTR,RC,SINE
fi
