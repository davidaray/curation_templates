#!/bin/sh
#SBATCH --job-name=<NAME>_pggb
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=12

. ~/miniforge3/etc/profile.d/conda.sh
# Activate the pantera environment
conda activate pantera

CHROMID=<NAME>

WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/pantera2
CHRPATH=$WORKDIR/chromosomes
CURATIONPATH=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates 
ARCHIVE=$WORKDIR/archive
TABLES=$WORKDIR/tables
RUNS=$WORKDIR/runs
# Define the path to the Singularity image
SINGULARITY_IMAGE=$WORKDIR/pggb_container/pggb_latest.sif
BIN=$WORKDIR/bin
PGGB=

cd $WORKDIR

# Check if $PGGB exists, then create it if it doesn't
if [ ! -d "$PGGB" ]; then
    echo "Creating directory: $PGGB"
    mkdir -p "$PGGB"
else
    echo "Directory $PGGB already exists."
fi

# Define the command to run inside the Singularity container
COMMAND="pggb -i $CHRPATH/${CHROMID}.fa.bz -o $PGGB -n 2 -t 10"

# Print the command for debugging
echo "$COMMAND"

# Run pggb
singularity exec --bind "$WORKDIR:$WORKDIR" "$SINGULARITY_IMAGE" $COMMAND
		
cd $BIN
dos2unix ${CHROMID}_pantera.sh
sbatch ${CHROMID}_pantera.sh




