#!/bin/bash

#SBATCH --job-name=blastn
#SBATCH --output=blastn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH -o blastn_%a.out # Standard output
#SBATCH -e blastn_%a.err # Standard error
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=6-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uni.ac.uk

## load profile
source ~/.bash_profile

USAGE="Usage: $(basename "$0") -F <relative path to fasta file> -B <absolute path to blast db> \n
The script assumes you have relatively few ASVs to identify (less that 1000) 
and blasts them against an ncbi database selected by the user.
Here is an example of how to run the script: \n
qsub b2m_scripts/01_run_blastn_simple.sh -F working_data/06_ASV_seqs.fasta -B /shared/genomicsdb2/shared/ncbi_nt/current/nt \n
The script makes a directory called blast_out and saves the results here. 
It also makes symbolic links to ncbi taxa databases. \n\n"

## List arguments
while getopts F:B: flag; do
	case "${flag}" in
		F) FASTA_PATH=${OPTARG};;
		B) DATABASE=${OPTARG};;
	esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if [ -z "${FASTA_PATH}" ] || [ -z "${DATABASE}" ]; then
   printf "\n\n${USAGE}" >&2; exit 1
fi

## PREP STEPS
## Define path variables
MAIN_DIR=$PWD
OUT_DIR="blast_out"
## create the output directory
mkdir -p ${MAIN_DIR}/${OUT_DIR}
## make a symlink to ncbi taxdb files, which need to be in the current directory
ln -s /shared/genomicsdb2/shared/ncbi_nr/current/taxdb* .

## run blast 
blastn -query ${MAIN_DIR}/${FASTA_PATH} -task blastn -db ${DATABASE} -out ${MAIN_DIR}/${OUT_DIR}/all_blast.out.tab -num_threads 4 -outfmt "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname scomnames sblastname sskingdoms stitle"
