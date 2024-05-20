#!/bin/bash

#SBATCH --job-name=blast2lca
#SBATCH --output=blast2lca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH -o blast2lca_%a.out # Standard output
#SBATCH -e blast2lca_%a.err # Standard error
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=6:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uni.ac.uk

## load profile and conda environment (maybe don't need conda environment)
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/metabarcoding

USAGE="Usage: $(basename "$0") -B <blast percent ID value 0-100> -M <megan percent of reads 0-100> -D <abs path of megan nucl database> \n
The script takes the output from blast (located in the blast_out directory), and applies the following: \n
step1: Merge results if blast was run in array mode (automatically detected);
step2: Filter blast results by percentage ID (0-100). The user provides a 
       minimum percentage ID  with -B for filtering blast results. We recommend 85-95;
step3: Run the megan2lca (lowest common ancestor) algorithm to determine taxonomic likelihood of 
       ASV at all taxonomic levels;
step4: Taxonomic identification based on level at which minimum percentage of alignments match. 
       The user provides a minimum percentage of matching alignments required with -M.
       We recommend 70-100;
Here is an example of how to run the script: \n
qsub b2m_scripts/02_run_blast2lca.sh -B 90 -M 100 -D /shared/genomicsdb2/shared/megan/megan-nucl-Feb2022.db \n
The script assumes blast results (whether from simple or array mode) are located in the directory 
blast_out, and also saves intermediate and final files to the same directory. \n\n"

## List arguments
while getopts B:M:D: flag; do
	case "${flag}" in
		B) BPI=${OPTARG};;
		M) MPI=${OPTARG};;
		D) MEG_DB=${OPTARG};;
	esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if [ -z "${BPI}" ] || [ -z "${MPI}" ] || [ -z "${MEG_DB}" ]; then
   printf "\n\n${USAGE}" >&2; exit 1
fi

## Define path variables
MAIN_DIR=$PWD
OUT_DIR="blast_out"

## Step 1 : Check if run in array mode and, if so, merge the chunks to create all_blast.out.tab
if [ -f "${MAIN_DIR}/${OUT_DIR}/chunk0.fa_blast.out.tab" ]; then
  echo "Blast was run in array mode, merging chunks..."
  cat ${MAIN_DIR}/${OUT_DIR}/chunk* > ${MAIN_DIR}/${OUT_DIR}/all_blast.out.tab
fi

## Step 2: Remove additional taxa information and filter by user specified blast percentage identity (blast threshold or BPI)
cut -f1-12 ${MAIN_DIR}/${OUT_DIR}/all_blast.out.tab | awk -v var="${BPI}" '$3 >= var' > ${MAIN_DIR}/${OUT_DIR}/filtered_blast.out.tab


## Step 3: Run the Megan blast2lca (lowest common ancestor) algorithm
/usr/local/extras/Genomics/apps/megan/tools/blast2lca -i ${MAIN_DIR}/${OUT_DIR}/filtered_blast.out.tab -m BlastN -o ${MAIN_DIR}/${OUT_DIR}/megan_full_out.tsv -mdb ${MEG_DB}

## Step 4: Take the full megan output and only show lowest common ancestor at a user specified threshold (megan threshold or MPI) 
## This includes a long conditional awk one-liner with if else statements nested within further if else statements. Ideally the 
## megan output would be consistent, but the number of fields (columns) can vary. After applying the cut command, values seen 
## include:
## > 1 or 2 (implying a failure to assign taxonomy)
## > 16 (missing one taxonomic level)
## > 18 (expected taxonomic levels to species) 
## > 20 (expected taxonomic levels to subspecies) 
## We allow 14 to 22 columns (14 in case it is missing two taxonomic levels, 22 in case an additonal rank appears). 
## If there are 14-22 columns but no blast2lca taxonomic assignment above the MPI threshold it receives the "y unknown' label. 
## If there are less than 14 columns it receives the "x no_assignment" label
## If there are more than 22 columns it receives the "z check_meganfull" label
cut -d ';' -f1,3,4,17- ${MAIN_DIR}/${OUT_DIR}/megan_full_out.tsv | \
awk -v var="${MPI}" -F ';' '
{ \
    if (NF < 14) {print $1 "\t" "x__no assignment"} \
    else if (NF == 14) { \
        if ($(NF-1) >= var) {print $1 "\t" $(NF-2)} \
        else if ($(NF-3) >= var) {print $1 "\t" $(NF-4)} \
        else if ($(NF-5) >= var) {print $1 "\t"  $(NF-6)} \
        else if ($(NF-7) >= var) {print $1 "\t"  $(NF-8)} \
        else if ($(NF-9) >= var) {print $1 "\t" $(NF-10)} \
        else if ($(NF-11) >= var) {print $1 "\t" $(NF-12)} \
        else {print $1 "\t" "y__unknown"} \
    } \
    else if (NF == 16) { \
        if ($(NF-1) >= var) {print $1 "\t" $(NF-2)} \
        else if ($(NF-3) >= var) {print $1 "\t" $(NF-4)} \
        else if ($(NF-5) >= var) {print $1 "\t"  $(NF-6)} \
        else if ($(NF-7) >= var) {print $1 "\t"  $(NF-8)} \
        else if ($(NF-9) >= var) {print $1 "\t" $(NF-10)} \
        else if ($(NF-11) >= var) {print $1 "\t" $(NF-12)} \
        else if ($(NF-13) >= var) {print $1 "\t" $(NF-14)} \
        else {print $1 "\t" "y__unknown"} \
    } \
    else if (NF == 18) { \
        if ($(NF-1) >= var) {print $1 "\t" $(NF-2)} \
        else if ($(NF-3) >= var) {print $1 "\t" $(NF-4)} \
        else if ($(NF-5) >= var) {print $1 "\t"  $(NF-6)} \
        else if ($(NF-7) >= var) {print $1 "\t"  $(NF-8)} \
        else if ($(NF-9) >= var) {print $1 "\t" $(NF-10)} \
        else if ($(NF-11) >= var) {print $1 "\t" $(NF-12)} \
        else if ($(NF-13) >= var) {print $1 "\t" $(NF-14)} \
        else if ($(NF-15) >= var) {print $1 "\t" $(NF-16)} \
        else {print $1 "\t" "y__unknown"} \
    } \
    else if (NF == 20) { \
        if ($(NF-1) >= var) {print $1 "\t" $(NF-2)} \
        else if ($(NF-3) >= var) {print $1 "\t" $(NF-4)} \
        else if ($(NF-5) >= var) {print $1 "\t"  $(NF-6)} \
        else if ($(NF-7) >= var) {print $1 "\t"  $(NF-8)} \
        else if ($(NF-9) >= var) {print $1 "\t" $(NF-10)} \
        else if ($(NF-11) >= var) {print $1 "\t" $(NF-12)} \
        else if ($(NF-13) >= var) {print $1 "\t" $(NF-14)} \
        else if ($(NF-15) >= var) {print $1 "\t" $(NF-16)} \
        else if ($(NF-17) >= var) {print $1 "\t" $(NF-18)} \
        else {print $1 "\t" "y__unknown"} \
    } \
    else if (NF == 22) { \
        if ($(NF-1) >= var) {print $1 "\t" $(NF-2)} \
        else if ($(NF-3) >= var) {print $1 "\t" $(NF-4)} \
        else if ($(NF-5) >= var) {print $1 "\t"  $(NF-6)} \
        else if ($(NF-7) >= var) {print $1 "\t"  $(NF-8)} \
        else if ($(NF-9) >= var) {print $1 "\t" $(NF-10)} \
        else if ($(NF-11) >= var) {print $1 "\t" $(NF-12)} \
        else if ($(NF-13) >= var) {print $1 "\t" $(NF-14)} \
        else if ($(NF-15) >= var) {print $1 "\t" $(NF-16)} \
        else if ($(NF-17) >= var) {print $1 "\t" $(NF-18)} \
        else if ($(NF-19) >= var) {print $1 "\t" $(NF-20)} \
        else {print $1 "\t" "y__unknown"} \
    } \
    else {print $1 "\t" "z__check meganfull"}
}' | awk -v OFS="\t" -F '[\t__]' '{print $1 "_" $2 "\t" $3 "\t" $5}' | tr ' ' '_' > ${MAIN_DIR}/${OUT_DIR}/megan_sum_out.tsv
