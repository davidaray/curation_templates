#!/bin/bash
#SBATCH --job-name=<NAME>_pantera
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

. ~/miniforge3/etc/profile.d/conda.sh
conda activate pantera

CHROMID=<NAME>

# Navigate to your working directory
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/pantera2
cd $WORKDIR

RCONTAINER=$WORKDIR/r_container
PANTERA_SOFTWARE_DIR=/lustre/work/daray/software/pantera
R_LIBS=/lustre/work/daray/software/pantera/libraries
PANTERAOUT=$WORKDIR/runs/$GENOME
R_CONTAINER=$WORKDIR/r_container
RUNS=/lustre/scratch/daray/bat1k_TE_analyses/pantera/runs
PGGB=
PANTERAOUT=

# Check if $PANTERAOUT exists, then create it if it doesn't
if [ ! -d "$PANTERAOUT" ]; then
    echo "Creating directory: $PANTERAOUT"
    mkdir -p "$PANTERAOUT"
else
    echo "Directory $PANTERAOUT already exists."
fi

# If the gfa list exists, remove it
[ -e ${CHROMID}.gfa.list.txt ] && rm ${CHROMID}.gfa.list.txt

# Generate the list of gfa files from pggb
for FILE in $PGGB/*.gfa; do
	echo $FILE >>$PGGB/${CHROMID}.gfa.list.txt
done

# Run Pantera
singularity exec --bind $PANTERA_SOFTWARE_DIR:$PANTERA_SOFTWARE_DIR,$PANTERAOUT:$PANTERAOUT,$PGGB:$PGGB,$WORKDIR:$WORKDIR \
	$RCONTAINER/r_container.sif Rscript $PANTERA_SOFTWARE_DIR/pantera.R \
	-c 12 \
	-o $PANTERAOUT \
	-g $PGGB/${CHROMID}.gfa.list.txt
