#!/bin/bash

#SBATCH --job-name=makesumfile
#SBATCH --output=makesumfile
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH -o makesumfile_%a.out # Standard output
#SBATCH -e makesumfile_%a.err # Standard error
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=6:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uni.ac.uk

## load profile
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/metabarcoding

USAGE="Usage: $(basename "$0") \n
The script calls an R script that combines megan_sum_out.tsv (located in blast_out) with the dada2 output files 06_ASV_seqs.fasta and 06_ASV_counts.tsv (located in working_data) and writes the output to blast_out/ASV_taxa-summary_counts.tsv. It will not run without all 3 files. 
The script does not require any arguments, and can be run as so: \n
qsub b2m_scripts/03_run_make_summary_file.sh \n"

MAIN_DIR=$PWD

# Run the R script to combine MEGAN, sequence and counts into a single summary file
Rscript $PWD/b2m_scripts/03_make_summary_file.R
