#!/bin/bash 
#SBATCH --job-name=<NAME>_hite 
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err 
#SBATCH --partition=nocona 
#SBATCH --mem-per-cpu=10G 
#SBATCH --nodes=1
#SBATCH --ntasks=36


WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/hite
TAXON=<NAME>
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies
mkdir $WORKDIR/$TAXON
CONTAINERPATH=/lustre/work/daray/software/HiTE
OUTPUTDIR=$WORKDIR/$TAXON
mkdir -p $OUTPUTDIR/assembly
gunzip -c $GZPATH/${TAXON}.fa.gz >$OUTPUTDIR/assembly/${TAXON}.fa

singularity run -B $OUTPUTDIR:/mnt --pwd /HiTE $CONTAINERPATH/HiTE.sif python main.py \
 --genome /mnt/assembly/${TAXON}.fa \
 --thread 36 \
 --outdir /mnt \
 --plant 0 \
 --miu 2.36e-09 \
 --annotate 1 \
 --domain 1 \
 --is_wicker 0
# [other parameters]

rm $OUTPUTDIR/assembly/${TAXON}.fa
rm $OUTPUTDIR/genome.rename.fa
 
 # (1) The options "--genome" and "--outdir" need to be specified as absolute paths.
 # (2) The option "-B" is used to specify directories to be mounted.
 #     It is recommended to set ${host_path} and ${container_path} to your user directory, and ensure 
 #     that all input and output files are located within the user directory.
 # (3) "--pwd /HiTE" and "python main.py" do not need to be changed.
 
 # e.g., my command: singularity run -B /home/hukang:/home/hukang --pwd /HiTE /home/hukang/HiTE.sif python main.py 
 # --genome /home/hukang/HiTE/demo/genome.fa 
 # --thread 40 
 # --outdir /home/hukang/HiTE/demo/test/
