#!/bin/bash 
#SBATCH --job-name=<NAME>_hite 
#SBATCH --output=%x.%j.out 
#SBATCH --error=%x.%j.err 
#SBATCH --partition=xlquanah 
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=10:00:00
#SBATCH --account=xlquanah

. ~/miniforge3/etc/profile.d/conda.sh
conda activate

echo "Hello world." >hello.txt
