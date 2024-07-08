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

## load profile. Conda environment required to run Rscript command
source ~/.bash_profile
conda activate /usr/local/extras/Genomics/apps/mambaforge/envs/metabarcoding

USAGE="Usage: $(basename "$0") \n
The script calls an R script that combines the two main output files from 02_run_blast2lca.sh with dada2 output files to make a final summary file:\n
The blast2lca files required are: \n
- blast_out/megan_taxonpath_out.tsv \n
- blast_out/megan_summary_out.tsv \n
Whilst the dada2 output files required are: \n
- working_data/06_ASV_seqs.fasta \n
- working_data/06_ASV_counts.tsv \n
These four files are combined and the output is written to:\n
- working_data/ASV_taxa_seq_counts.tsv \n
Additionally the script will output three files for phyloseq analyses:\n
- working_data/ps_taxamat.tsv \n
- working_data/ps_countmat.tsv \n
- working_data/ps_phylogeny.rds \n
The script will only work if all 4 input files are found, but does not require any arguments, and can be run as so: \n
qsub b2m_scripts/03_run_make_summary_files.sh \n"

MAIN_DIR=$PWD

# Run the R script to combine MEGAN, sequence and counts into a single summary file
Rscript $PWD/b2m_scripts/03_make_summary_files.R
