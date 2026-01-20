#!/bin/bash
#SBATCH --job-name=deepte
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1

. ~/miniforge3/etc/profile.d/conda.sh
conda activate deepte

cd /lustre/scratch/daray/bat1k_TE_analyses/deepte

DEEP=/lustre/work/daray/software/DeepTE
HMMER=/lustre/work/daray/software/hmmer-3.3.2/bin
CURATION_DIR=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates
LIBRARY=/lustre/scratch/daray/bat1k_TE_analyses/hite/naming_work/final_hite_library_20240915g.fa

DATE=$(date)
echo -e "Start time: $DATE\n"

python $CURATION_DIR/pull_using_header_text.py \
	-i $LIBRARY \
	-o hite_unknowns.fa \
	-s Unknown

python $CURATION_DIR/clean_sequences.py \
	-i hite_unknowns.fa \
	-o hite_clean_unknowns.fa
	
mkdir working
mkdir output

DATE=$(date)
echo -e "Start DeepTE_domain: $DATE\n"

python $DEEP/DeepTE_domain.py -d working -o output -i hite_clean_unknowns.fa -s $DEEP/supfile_dir --hmmscan $HMMER/hmmscan

sleep 10

DATE=$(date)
echo -e "Start DeepTE: $DATE\n"

python $DEEP/DeepTE.py -d working -o output -i hite_clean_unknowns.fa -sp M -m_dir $DEEP/Metazoans_model/ -modify output/opt_te_domain_pattern.txt

: << 'EOF'
DeepTE documentation (https://github.com/LiLabAtVT/DeepTE?tab=readme-ov-file)
usage:
**DeepTE**
DeepTE.py [-h] required: [-d working_dir][-o output_dir]
                         [-i ipt_seq][-sp sp_type]
                         ([-m model_name]|[-m_dir model_dir])
               optional: [-modify domain_file]
                         [-fam te_fam][-UNS yes][-prop_thr value]

arguments:
-h, --help        Show this help message and exit.

-d                Working directory to store intermediate files of each step. 
                  Default: ./.

-o                Output directory to store the output files. 
                  Default: ./.

-i                Input sequences that are unknown TE or DNA sequences.

-sp               P or M or F or O. P:Plants, M:Metazoans, F:Fungi, and O: Others.

-m                Provide one of model names: 
                  '-m P' or '-m M' or '-m F' or '-m O' or '-m U'.
                  This argument will directly download the model dir.
                  Users do not need to initiate '-m_dir'.
                  If users do not want to directly download model, please use '-m_dir', but users need to download model directory by themselves.

-m_dir            Provide model_dir that could be downloaded from website (optional requirements). 
                  If users set -UNS yes, please provide UNS_model directory that can be downlowed in the above link.

-fam              Provide TE family name for the input te sequence
                  ClassI: the input sequence is ClassI TEs
                  ClassII: the input sequence is ClassII subclass1 TEs
                  LTR: the input sequence is LTR TEs
                  nLTR: the input sequence is nLTR TEs
                  LINE: the input sequence is LINE TEs
                  SINE: the input sequence is SINE TEs
                  Domain: the input sequence is Class II subclass1 TEs with specified super families
                  If users do not initiate '-fam', DeepTE will regard your input sequences are unknown TEs.

-modify           If set this argument, users need to provide domain file generated from another script: DeepTE_domain.py.

-UNS              If set this argument, users need change the -i to the the DNA sequences; 
                  This function will classify the sequences into TEs, CDS, or Intergenic sequences; -sp and -fam do not need to provide.
                  Note: this model is used for plants rather than metazoans and fungi.

-prop_thr         Specify a probability threshold to annotate TE.
                  For example: a TE has a probability (0.6) to be ClassI.
                  If users set 0.7 as the threshold, 
                  this TE will be labeled as 'unknown', Default: 0.6.


**DeepTE_domain**
DeepTE_domain.py [-h] required: [-d working_dir][-o output_dir]
                                [-i ipt_seq][-s supfile_dir]
                                [--hmmscan hmmscan]
arguments:
-h, --help        Show this help message and exit.

-d                Working directory to store intermediate files of each step. 
                  Default: ./.

-o                Output directory to store the output files. 
                  Default: ./.

-i                Input sequences that are unknown TE sequences.

-s                Provide supplementary dir that contains required files.

--hmmscan         File path to hmmscan executable, Default: /usr/bin/hmmscan"
EOF
