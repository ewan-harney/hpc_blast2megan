#!/bin/bash

#SBATCH --job-name=blast_prep
#SBATCH --output=blast_prep
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH -o blast_prep_%a.out # Standard output
#SBATCH -e blast_prep_%a.err # Standard error
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uni.ac.uk

## load profile and conda environment (maybe don't need conda environment)
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/metabarcoding

USAGE="Usage: $(basename "$0") -F <relative path to fasta file> \n
The script assumes you have a relatively large number of ASVs to identify (more that 1000) 
and splits them into chunks of 100 ASV sequences.
Here is an example of how to run the script: \n
qsub b2m_scripts/01A_run_prep_for_blast.sh -F working_data/06_ASV_seqs.fasta \n
The script makes a directory called split_fasta and saves fasta chunks here.
It also creates a list of the files (required for step 01B) and  symbolic links 
to ncbi taxa databases (also required for step 01B). \n\n"

## List arguments
while getopts F: flag; do
	case "${flag}" in
		F) FASTA_PATH=${OPTARG};;
	esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if [ -z "${FASTA_PATH}" ]; then
   printf "\n\n${USAGE}" >&2; exit 1
fi

## PREP STEPS
## Define path variables
MAIN_DIR=$PWD
PREP_DIR="split_fasta"
## create the directory to put split files, and logs produced by the 01B script
mkdir -p ${MAIN_DIR}/${PREP_DIR}
mkdir -p ${MAIN_DIR}/logs
## make a symlink to ncbi taxdb files, which need to be in the current directory for the 01B script
ln -s /shared/genomicsdb2/shared/ncbi_nr/current/taxdb* .
## Change to the split_fasta directory
cd ${MAIN_DIR}/${PREP_DIR}

## split fasta file into chunks each with 100 sequences and generate list.txt file
awk 'BEGIN {n=0;} /^>/ {if(n%100==0){file=sprintf("chunk%d.fa",n);} print > file; n++; next;} { print >> file; }' < ${MAIN_DIR}/${FASTA_PATH}
NUM=$(ls chunk* | wc -l )
ls chunk* > ${MAIN_DIR}/${PREP_DIR}/${PREP_DIR}_list_of_${NUM}.txt
