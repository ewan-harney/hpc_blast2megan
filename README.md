<br>
<br>
## Assigning taxonomy with BLASTn and MEGAN on UoS BESSEMER.
<br>
<font size="4">
<details><summary><font size="6"><b>1) About, credits, and other information</b></font></summary>
  <br>
  <br>
  This short HPC tutorial uses BLASTn to query the identity of nucleotide sequence data against an ncbi database, and applies the MEGAN LCA (lowest common ancestor) algorithm to provide a probable taxonomic assignment.

  The workflow was designed as an alternative method of taxonomic assignment to the dada2 "assign taxonomy" step featured in step 11 of Katy Maher's [dada2 pipeline] (https://github.com/khmaher/HPC_dada2), which is primarily designed for microbial data sets.

  For metabarcoding projects featuring eukaryotic data, using blastn to query the sequences against ncbi's nt database can provide a reliable means of identifying sequences with diverse taxonomic origins. By default BLAST will output many sequences that match the query sequence, ain some cases making it hard to tell which species the query sequence belongs to. THE LCA algorthm of MEGAN takes the 

  Although it was designed to follow on from the dada2 pipeline, this workflow can be applied to query any nucleotide sequence data in fasta format.

  The code has been written for use with the University of Sheffield's [BESSEMER](https://docs.hpc.shef.ac.uk/en/latest/bessemer/index.html) system but should be applicable to any GNU/Linux based HPC system once the appropriate modifications are made (your mileage may vary).

  Code which the user must run is highlighted in a code block like this:
  ```
  I am code - you must run me
  ```
  Sometimes the desired output from a command is included in the code block as a comment.
  For example:
  ```
  Running this command
  # Should produce this output
  ```

  Filepaths within normal text are within single quote marks, like this:

  '/home/user/a_file_path'
  <br><br>
  Contact: Ewan Harney //  e.harney@sheffield.ac.uk
  </details>
<br>
<details><summary><font size="6"><b>2) Getting set up.</b></font></summary>
  <br>
  <br>
  <font size="4"><b>2.1) Access the HPC</b></font>
  <br>
  This workflow assumes you have already been using BESSEMER to run the dada2 pipeline. If that's not the case and you wish to get set up on this particular HPC, please refer to sections 2.1 - 2.4 of the [dada2 pipeline] (https://github.com/khmaher/HPC_dada2)
  <br>

  <font size="4"><b>2.2) Navigate to your working directory</b></font>
  <br>
  Navigate to your project directory. If you have been running the dada2 analysis you likely have a 'my_project' directory within the '/fastdata' directory on BESSEMER. Within 'my_project' is the 'working_data' directory, which contains the sequence data in a file called '06_ASV_seqs.fasta'. If you have not been running the dada2 pipeline, you can navigate to the directory containing the sequence data or a parent directory, whichever you prefer (the relative path to the fasta file will be specified when running the script).
  <br>

  <font size="4"><b>2.3) Copy blast2megan scripts</b></font>
  <br>
  Clone (download) this github repository, and copy the b2m_scripts directory contained within in to your current location. You can then delete the github download.

  ```
  git clone "https://github.com/ewan-harney/hpc_blast2megan"
  cp -r hpc_blast2megan/b2m_scripts .
  rm -r hpc_blast2megan
  ```

  Check the contents of the b2m_scripts directory. There should be 5 files in the directory: 5 .sh files and 1 .R script

  ```
  ls b2m_scripts
  ```
  <br>
  <font size="4"><b>2.4) A note on editing scripts</b></font>
  <br>
  Unlike the scripts in the dada2 pipeline, the user does not provide their email address as a command line argument for these scripts, and by default will not receive email confirmation of job completion. However this can be easily altered through a small change to the resource request section of the .sh scripts. A script can be viewed and edited with the nano command, followed by the relative or absolute path to the script, e.g. :

  ```
  nano b2m_scripts/01_run_blastn_simple.sh
  ```

This will start nano. Notice that the first line of the script is #!/bin/bash, followed by an empty line, and then several lines commencing with #SBATCH. These #SBATCH arguments are used by slurm when the script is submitted (with qsub or sbatch) and allow the user to control certain parameters. Notice that the last #SBATCH line is:

#SBATCH --mail-user=user@uni.ac.uk. 

Using the arrow key, change user@uni.ac.uk to your own email address. In my case this argument would be changed to read:

#SBATCH --mail-user=e.harney@sheffield.ac.uk

Once you have made this change, you will need to save the changes. Notice at the bottom of the screen are lines of commands, such as '^G Get Help' and '^X Exit' etc. The '^' means holding down the control (windows) or command (mac) key. Pressing the 'X' key whilst holding down control/command will allow you to exit nano and save any changes. After pressing ^X you will be prompted to save the changes (options are y for yes, n for no and ^C for cancel). Press y. You will then be given the chance to rename the file if you want. In our case that's not necessary, so simply press enter to save the file with the same name and exit nano.

  <br>
  </details>
<br>
<details><summary><font size="6"><b>3) Blast sequence data against an ncbi database.</font></b></summary>
  <br>
  <br>
  <font size="4"><b>3.1) Determine how many sequences are in your fasta file</b></font>
  <br>
  The first step is to query sequence data contained with a fasta file against an ncbi database using the blastn (n for nucleotide) algorithm. Depending on the size and number of the sequences in the fasta file and the database being used, this step can be quite slow. If your fasta file contains thousands of sequences, the process can be sped up by slitting the fasta file into chunks and running blast in parallel using the array functionality of slurm. This workflow therefore contains 2 different options for running blast. If you have relatively few sequences (< 1000) you can run the single script '01_run_blastn_simple.sh'. However, if you have a relaitvely large number of sequences in your fasta file (>1000), we recommend splitting the file into chunks of 100 and running it in array mode. To do this you will run '01A_run_split_fasta.sh' and then '01B_run_blastn_array.sh'. To see how many sequences are in your fasta file, run the following:

  ```
  grep -c '>' working_data/06_ASV_seqs.fasta
  ```

  Although we use 1000 sequences as our cut-off, you can run an array with less than 1000 sequences. Equally, you can run the simple blast even if you have more than 1000 sequences, although the script may take a while to finish. In section 3.2 we describe how to run blast in simple mode, and in section 3.3 we describe how to run it in array mode.

  <br>

  <font size="4"><b>3.2) Running blastn in simple mode</b></font>
  <br>
  Running blastn in simple mode will create a new directory called blast_out in your current directory (unless this directory already exists), and also create symolic links to the ncbi taxadb files 'taxdb.btd' and 'taxdb.bti' in the current directory. It will then run blastn and the output will be saved as blast_out/all_blast.out.tab. 
  
  <b>To run the 01_run_blastn_simple.sh script you must supply two arguments:<b>
  - the relative path to the fasta file containing the sequence data (-F)
  - the location of an ncbi database on the HPC (-B)
  
  It is most likely that you will use the database nt, which contains all nucleotide sequences available on the GenBank DNA sequence database. However, in some cases you may have a a smaller or bespoke database available against which you can blast your sequences. In this README we will assume you are using nt. An example command if you have run the dada2 pipeline might be:
  
  ```
  qsub b2m_scripts/01_run_blastn_simple.sh -F working_data/06_ASV_seqs.fasta -B /shared/genomicsdb2/shared/ncbi_nt/current/nt
  ```

  <br>

  <font size="4"><b>3.3) Running blastn in array mode</b></font>
  <br>
  Running blastn in array mode requires running 2 scripts one after the other: first '01A_run_split_fasta.sh' and then '01B_run_blastn_array.sh'. The '01A_run_split_fasta.sh' script will create a directory called split_fasta. The input sequence fasta file will then be split into chunks each containing 100 sequences which will be written to split_fasta. Like the 01_run_blastn_simple.sh this script will create symolic links to the ncbi taxadb files 'taxdb.btd' and 'taxdb.bti', and will also create a directory called 'logs', which will be used by script 01B. As well as creating the chunk.fa files, it will also create a text file 'split_fasta_list_of_X.txt' with the names of the chunks for the next step. In your file the 'X' will be the number of chunks in split_fasta. This number is a parameter for script '01B_run_blastn_array.sh'.
  
  <b>To run the 01A_run_split_fasta.sh script you just need to provide the path to the sequence data:<b>
  - the relative path to the fasta file containing the sequence data (-F)
  
  An example command if you have run the dada2 pipeline might be:
  
  ```
  qsub b2m_scripts/01A_run_split_fasta.sh -F working_data/06_ASV_seqs.fasta
  ```
  
  The '01B_run_blastn_array.sh' script will then use an array to simultaneously blast multiple chunk.fa files against an ncbi database. This script will  create a new directory called blast_out in your current directory (unless this directory already exists) and write the output of each chunk to a seperate chunk.fa_blast.out.tab.
  
  <b>To run the 01A_run_split_fasta.sh script you must supply two arguments:<b>
  - the location of an ncbi database on the HPC (-B)
  - the number of input files to be run on the array (-N)
  
  As stated in section 3.2, it is most likely that you will use the database nt. The number -N is contained in the file name of 'split_fasta_list_of_X.txt' (in place of the 'X') which can be viewed with the following:
  
  ```
  ls split_fasta/split_fasta*
  ```
  
  Slurm job arrays allow many jobs to be submitted simultaneously and run in parallel. For more information on arrays refer to the Sheffield HPC documentation on [advanced job submission] (https://docs.hpc.shef.ac.uk/en/latest/hpc/scheduler/advanced_job_submission_and_control.html#gsc.tab=0). With the script itself are special arguments relating to the array, and the job submission itself is also different. If our original sequence.fasta file contained 2350 sequuences, it would have been split into 24 chunks, with the txt file named split_fasta_list_of_24.txt. This number, 24, will appear twice when we submit this job, which would be as follows:
  
  ```
  sbatch --array=1-24 b2m_scripts/01B_run_blastn_array.sh -B /shared/genomicsdb2/shared/ncbi_nt/current/nt -N 24
  ```
  
  Notice that we use sbatch instead of qsub, and that this is followed by array=1- and then the number specific to our data set. This number also appears at the end of the command following the -N flag. Also, an error and output log file for each job of the array will be written to the directory 'logs'
  
  <br><br>
  <font size="4"><b>3.4) Monitoring and assessing the result of blastn</b></font>
  <br>
  
  Blastn against the nt database can take a while to run, even if there are not too many sequences to assess. To follow the status of the job run the following command: 

  ```
  squeue --me
  ```
  
  For more information about the squeue output refer to the Sheffield HPC documentation on [squeue] (https://docs.hpc.shef.ac.uk/en/latest/referenceinfo/scheduler/SLURM/Common-commands/squeue.html#gsc.tab=0) squeue will show the status of the job, and in the case of an array, will how many of the subjobs have been submitted and how many are still queued.
  
  If blast was run in simple mode, blast_out should now contain a single file called all_blast.out.tab, and if it was run in array mode, it will contain several chunk.fa_blast.out.tab files. Look at the contents of one of the files with:
  
  ```
  head blast_out/all_blast.out.tab 
  ```
  
  or
  
  ```
  head blast_out/chunk0.fa_blast.out.tab
  ```
  
  The file(s) should be similar to the following image. Information about blast tabular output can be found at the [Metagenomics wiki] (https://www.metagenomics.wiki/tools/blast/blastn-output-format-6). The output contains a few differences from the default. The column headers in these files correspond to qseqid *saccver* pident length mismatch gapopen qstart qend sstart send evalue bitscore *staxid ssciname scomnames sblastname sskingdoms stitle* (italics highlight differet or additional columns).
  
  It is highly likely that all the rows displayed by head (the top 10) show results for the same sequence (ASV_1 if following the dada2 pipeline). This is not a mistake! Our query sequences are likely to match many sequences in the nt database. Sometimes the alignment will be much better for one species than any other, which allows us to confidently assign the sequence to that species. But often the sequence will align to multiple sequences in the database equally well. In this case, we need to class the sequence at a lower taxonomic level (e.g. genus or family). This is what we will do in the next step using the MEGAN blast2lca algorithm.
  
  <br>
  </details>
<br>
<details><summary><font size="6"><b>4) Run the MEGAN blast2lca algorithm.</font></b></summary>
  <br>
  <br>
  <font size="4"><b>4.1) Run blast2lca with stringent parameters  </b></font>
  <br>
  
  For this step we will uses [MEGAN] (https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/), a suite of bioinformatic algorithms developed by Daniel H Huson at the Universtiy of Tuebingen together with other collaborators. MEGAN contains various tools to help with the analysis of metagenomic and environmental DNA datasets. The specific tool we are interested in blast2lca, which calculates the LCA or [lowest common ancestor] (https://en.wikipedia.org/wiki/Lowest_common_ancestor) of the best results from blast.
  
  The 02_run_blast2lca.sh script takes the output from blast (located in the blast_out directory), and does the following:

  step 1: If blast was run in array mode chunks are merged (this is automatically detected);
  step 2: Blast results are filtered by pident (the percentage of identical positions: column 3 of the output). The user must provide a minimum pident (between 0 and 100) with the -B argument to filtering blast results.
  step 3: The megan2lca (lowest common ancestor) algorithm is used to determine taxonomic likelihood of a sequence at all taxonomic levels;
  step 4: Taxonomic identification is carried out based on the percentage of alignments matching at a taxonomic level. The user must provide a minimum percentage of matching alignments with the -M argument.
  
  This script will provide 2 intermediate files (filtered_blast.out.tab and megan_full_out.tsv) and it's final output file, megan_sum_out.tsv can be used in downstream analysis.

  <b>To run the 02_run_blast2lca.sh script you must supply three arguments:<b>
  - Minimum percentage identity (0-100) for the blast results to be considered by blast2lca (-B)
  - Minimum percentage of matching alignments (0-100) for taxonomic assignemnt in blast2lca (-M)
  - Absolute path to the megan nucleotide database (-D)
  
  If running the analysis on BESSEMER, the database should be available at '/shared/genomicsdb2/shared/megan/megan-nucl-Feb2022.db'. Otherwise you can download your own version from the [MEGAN Alternative Download Page] (https://unitc-my.sharepoint.com/personal/iijhu01_cloud_uni-tuebingen_de/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fiijhu01%5Fcloud%5Funi%2Dtuebingen%5Fde%2FDocuments%2FApps%2FMegan&ga=1). For deciding which values of -B and -M to use, we recommend initially using relatively strict (high) values for both, such as:
  
  - -B 95
  - -M 100
  
  Initially, we would thus recommend running the job like so:
  
  ```
  qsub scripts/02_run_blast2lca.sh -B 95 -M 100 -D /shared/genomicsdb2/shared/megan/megan-nucl-Feb2022.db
  ```
  
  The columns of the output file, megan_sum_out.tsv are as follows:
  ASV number / taxonomic rank / LCA taxon 

  <br><br>
  <font size="4"><b>4.2) Tweaking the parameters of the blast2lca script</b></font>
  <br>

  Once you have run blast2lca you can look at the taxa of your first few sequences with 
  
  ```
  head blast_out/megan_sum_out.tsv
  ```
  
  And to get a rough idea of how well blast and megan have worked, you can run the following code:

  ```
  cut -f2 blast_out/megan_sum_out.tsv | sort | uniq -c 
  ```
  
  This will show how many sequences have been assigned to each taxonomic rank. Hopefully the majority of your sequences are at s (species / subspecies) or g (genus) level. The x, y and z ranks (if they appear) refer to different failures:
  x (no assignment) there were no assignments at the given megan threshold (MPI)
  y (unknown) the megan2lca algorithm was unable to assign taxonomy (irrespective of the megan threshold) 
  z (check_meganfull) there may be a problem with the blast2lca output for this sequence - refer to that sequence in the megan_full_out.tsv result
  
  It is likely that there will be some non assigned and unknown taxa in your data, as well as some taxa assigned to lower levels such as class and order. Deciding what is a 'good' result will depend on many factors including the type of data, the organisms expected, the sampling strategy, the primers used (in metabarcoding) and a range of other factors. If you are not satisfied, you can try tweaking the values of -B and -M. 
  
  Reducing -B allows sequences with moderate blast percentage identity to be considered by blast2lca, and may reduce the number of sequences without assignment. Reducing -M can allow higher level (more specific) taxonomic assignemnts, with the caveat that there is a greater chance of misidentification.

  <br>
  </details>
<br>
<details><summary><font size="6"><b>5) Combine with dada2 output to create a summary file</font></b></summary>
  <br>
  <br>
  <font size="4"><b>5.1) Run the summary file script  </b></font>
  <br>
  
  If you have run this pipeline after dada2, you may wish to combine taxonomic assignment results with sequence data and ASV counts from dada2 to create a summary file with all information.
  
  Assuming you have created the megan_sum_out.tsv file in the previous step and have the files 06_ASV_seqs.fasta and 06_ASV_counts.tsv in your working_data directory you should be able to run the following script (no arguements need to be supplied):
  
  ```
  qsub scripts/03_run_make_summary_file.sh
  ```
  
  This script will call the R script 03_make_summary_file.R and write the output to blast_out/ASV_taxa-summary_counts.tsv. This file contains a summary of taxonmic, sequence and count results.
