#!/bin/sh
#SBATCH --job-name=<NAME>_compleasm
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --time=003:00:00

## Documentation
: '
MAIN MODULES
run             Run compleasm including miniprot alignment and completeness evaluation
analyze         Evaluate genome completeness from provided miniprot alignment
download        Download specified BUSCO lineage
list            List local or remote BUSCO lineages
miniprot        Run miniprot alignment

USAGE
python compleasm.py run [-h] -a ASSEMBLY_PATH -o OUTPUT_DIR [-t THREADS] 
                        [-l LINEAGE] [-L LIBRARY_PATH] [-m {lite,busco}] [--specified_contigs SPECIFIED_CONTIGS [SPECIFIED_CONTIGS ...]] 
                        [--miniprot_execute_path MINIPROT_EXECUTE_PATH] [--hmmsearch_execute_path HMMSEARCH_EXECUTE_PATH] 
                        [--autolineage] [--sepp_execute_path SEPP_EXECUTE_PATH] 
                        [--min_diff MIN_DIFF] [--min_identity MIN_IDENTITY] [--min_length_percent MIN_LENGTH_PERCENT] 
                        [--min_complete MIN_COMPLETE] [--min_rise MIN_RISE]

PARAMETERS
  -a, --assembly_path        Input genome file in FASTA format
  -o, --output_dir           The output folder
  -t, --threads              Number of threads to use
  -l, --lineage              Specify the name of the BUSCO lineage to be used. (e.g. eukaryota, primates, saccharomycetes etc.)
  -L, --library_path         Folder path to download lineages or already downloaded lineages. 
                             If not specified, a folder named "mb_downloads" will be created on the current running path by default to store the downloaded lineage files.
  -m, --mode                 The mode of evaluation. Default mode is busco. 
                             lite:  Without using hmmsearch to filtering protein alignment.
                             busco: Using hmmsearch on all candidate predicted proteins to purify the miniprot alignment to improve accuracy.
  --specified_contigs        Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.
  --outs                     output if score at least FLOAT*bestScore [0.95]
  --miniprot_execute_path    Path to miniprot executable file. 
                             If not specified, compleasm will search for miniprot in the directory where compleasm.py is located, the current execution directory, and system environment variables.
  --hmmsearch_execute_path   Path to hmmsearch executable file.
                             If not specified, compleasm will search for hmmsearch in the directory where compleasm.py is located, the current execution directory, and system environment variables.
  --autolineage              Automatically search for the best matching lineage without specifying lineage file.
  --sepp_execute_path        Path to sepp executable file. This is required if you want to use the autolineage mode.


  --min_diff               The thresholds for the best matching and second best matching. default=0.2
  --min_identity           The identity threshold for valid mapping results. default=0.4
  --min_length_percent     The fraction of protein for valid mapping results. default=0.6
  --min_complete           The length threshold for complete gene. default=0.9

EXAMPLES
# with lineage specified
python compleasm.py run -a genome.fasta -o output_dir -l eukaryota -t 8

# autolineage mode
python compleasm.py run -a genome.fasta -o output_dir -t 8 --autolineage

# with custom specified already downloaded lineage folder
python compleasm.py run -a genome.fasta -o output_dir -l eukaryota -t 8 -L /path/to/lineages_folder

# specify contigs
python compleasm.py run -a genome.fasta -o output_dir -l eukaryota -t 8 --specified_contigs chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22

LINEAGE DOWNLOAD
python compleasm.py download [-h] [-L LIBRARY_PATH] lineages [lineages ...]

positional arguments:  
  lineages                Specify the names of the BUSCO lineages to be downloaded. (e.g. eukaryota, primates, saccharomycetes etc.)

optional arguments:
  -L, --library_path      The destination folder to store the downloaded lineage files.
                          If not specified, a folder named "mb_downloads" will be created on the current running path by default.
'

### Begin work

. ~/miniforge3/etc/profile.d/conda.sh
conda activate compleasm

ASSEMBLY=<NAME>
WORKDIR=/lustre/scratch/daray/bat1k_TE_analyses/compleasm/$ASSEMBLY
GZPATH=/lustre/scratch/daray/bat1kdatafreeze/final_assemblies
LINEAGE=mammalia_odb10
LIBRARYPATH=/lustre/scratch/daray/bat1k_TE_analyses/compleasm/library

if [ -f $WORKDIR ]; then
	echo ${WORKDIR}" exists."
else 
	mkdir -p $WORKDIR
fi

cd $WORKDIR

# Lineage download
if [ -f $LIBRARYPATH/eukaryota_odb10.done ]; then
	echo "${LINEAGE} already downloaded."
else
    compleasm download -L $LIBRARYPATH $LINEAGE && touch $LIBRARYPATH/eukaryota_odb10.done
fi

# Compleasm run
compleasm run \
	-a ${GZPATH}/${ASSEMBLY}.fa.gz \
	-o $WORKDIR \
	-l $LINEAGE \
	-t 36 \
	-L $LIBRARYPATH

rm -rf $LINEAGE/

sed "1i ${ASSEMBLY}" summary.txt >${ASSEMBLY}_compleasm_summary.txt

cat ${ASSEMBLY}_compleasm_summary.txt >>../all_summaries.txt
