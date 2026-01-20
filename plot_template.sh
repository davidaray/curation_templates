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

# Generate plots if they haven't been created
    echo "Generating stacked bar and pie for $GENOME, all categories."
    python $CURATIONDIR/te_plots.py \
        -f out \
		-i ${GENOME}.fa.out.gz \
        -s ${GENOME}.summary.gz \
        -o ${GENOME}.plot \
        --plot_type both \
        --num_proc 12 \
        -bin 1

    echo "Generating stacked bar and pie for $GENOME, selected categories."
    python $CURATIONDIR/te_plots.py \
        -f out \
		-i ${GENOME}.fa.out.gz \
        -s ${GENOME}.summary.gz \
        -o ${GENOME}.plot.nosat \
        --plot_type both \
        --num_proc 12 \
        -bin 1 \
        --selected_classes DNA DIRS LINE LTR RC SINE
		
		    echo "Generating stacked bar and pie for $GENOME, selected categories."
    python $CURATIONDIR/te_plots.py \
        -f out \
		-i ${GENOME}.fa.out.gz \
        -s ${GENOME}.summary.gz \
        -o ${GENOME}.plot.nounknown \
        --plot_type both \
        --num_proc 12 \
        -bin 1 \
        --selected_classes DNA DIRS LINE LTR RC SINE Satellite
