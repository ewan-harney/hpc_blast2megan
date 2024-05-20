#!/bin/bash

#SBATCH --job-name=blastn_array
#SBATCH --output=blastn_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH -o logs/blastn_array_%a.out # Standard output
#SBATCH -e logs/blastn_array_%a.err # Standard error
#SBATCH --mem-per-cpu=5GB
#SBATCH --time=6-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uni.ac.uk

## load profile and conda environment (maybe don't need conda environment)
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/metabarcoding

USAGE="Usage: $(basename "$0") -B <absolute path to blast db> -N <number of files in the array> \n
This script must only be run following script 01A_run_prep_for_blast.sh. It uses 
slurm to submit the job as an array, speeding up the identification of sequences.
First you must determine the size of the array. This appears in the file name of 
split_fasta/split_fasta_list_of_*.txt. \n
This number must be entered TWO TIMES: 1) as the array limit AND 2) as a mandatory argument for -N. \n
Array job submission differs to normal batch job submission; here's an example of
how to run it if you have 24 files, i.e. the list is called split_fasta_list_of_24.txt: \n
sbatch --array=1-24 b2m_scripts/01B_run_blastn_array.sh -B /shared/genomicsdb2/shared/ncbi_nt/current/nt -N 24 \n
The script makes a directory called blast_out and saves the results here. \n\n"

## List arguments
while getopts B:N: flag; do
	case "${flag}" in
		B) DATABASE=${OPTARG};;
		N) NUM=${OPTARG};;
	esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if [ -z "${DATABASE}" ] || [ -z "${NUM}" ]; then
   printf "\n\n${USAGE}" >&2; exit 1
fi

## PREP STEPS
## Define path variables
MAIN_DIR=$PWD
PREP_DIR="split_fasta"
OUT_DIR="blast_out"
mkdir -p ${MAIN_DIR}/${OUT_DIR}
DATA=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <(cat ${PREP_DIR}/${PREP_DIR}_list_of_${NUM}.txt))
FASTA=$(echo "$DATA" | cut -f1 )

## run blast 
blastn -query ${MAIN_DIR}/${PREP_DIR}/${FASTA} -task blastn -db ${DATABASE} -out ${MAIN_DIR}/${OUT_DIR}/${FASTA}_blast.out.tab -num_threads 4 -outfmt "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname scomnames sblastname sskingdoms stitle"
