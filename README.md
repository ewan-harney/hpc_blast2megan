## Assigning taxonomy with BLASTn and MEGAN on UoS BESSEMER.
This short HPC tutorial uses BLASTn to query the identity of nucleotide sequence data against an ncbi database, and applies the MEGAN LCA (lowest common ancestor) algorithm to provide a probable taxonomic assignment.
<br></br>
<font size="4">
<details><summary><font size="6"><b>1) About, credits, and other information</b></font></summary>
<br></br>

The workflow was designed as an alternative method of taxonomic assignment to the dada2 "assign taxonomy" step featured in step 11 of Katy Maher's [dada2 pipeline](https://github.com/khmaher/HPC_dada2), which is primarily designed for bacterial and microbial data sets.
<br></br>

For metabarcoding projects featuring non-microbial eukaryotic data, using blastn to query sequences against ncbi databases provides a reliable means of taxonomic identification. However, the ncbi database may contain multiple sequences that match the query sequence, making it hard to tell which taxa the query sequence belongs to. THE LCA (lowest common ancestor) algorthm of MEGAN assigns taxonomy based on the lowest common taxonomic ancestor of the species identified by BLAST. 
<br></br>

Although designed to follow the dada2 pipeline, this workflow can be applied to query any nucleotide sequence data in fasta format. The code has been written for use with the University of Sheffield's [BESSEMER](https://docs.hpc.shef.ac.uk/en/latest/bessemer/index.html) system but should be applicable to any GNU/Linux based HPC system once the appropriate modifications are made (your mileage may vary).
<br></br>

Code which the user must run is highlighted in a code block like this:

```
I am code - you must run me
```

Filepaths, directory names and file names within normal text are within single quotes, like this:

* '/home/user/a_file_path'
<br></br>

When a specific button should be pressed this is shown in double quotes, like so:

* press "y" then "enter"
<br></br>

Contact: Ewan Harney //  e.harney@sheffield.ac.uk
<br>
</details>
<br>

<details><summary><font size="6"><b>2) Getting set up</b></font></summary>
<br></br>

<font size="4"><b>2.1) Access the HPC</b></font>
<br></br>
This workflow assumes you have already been using BESSEMER to run the dada2 pipeline. If that's not the case and you wish to get set up on this particular HPC, please refer to sections 2.1 - 2.4 of the [dada2 pipeline](https://github.com/khmaher/HPC_dada2)
<br></br>

<font size="4"><b>2.2) Navigate to your working directory</b></font>
<br></br>
Navigate to your project directory. If you have been running the dada2 analysis you likely have a 'my_project' directory within the '/fastdata' directory on BESSEMER. Within 'my_project' is the 'working_data' directory, which contains the sequence data in a file called '06_ASV_seqs.fasta'. If you have not been running the dada2 pipeline, you can navigate to the directory containing your sequence data or a parent directory, whichever you prefer (the relative path to the fasta file will be specified when running the script).
<br></br>

<font size="4"><b>2.3) Copy blast2megan scripts</b></font>
<br></br>
Clone (download) this github repository, and copy the b2m_scripts directory contained within to your current location.
  
```
git clone "https://github.com/ewan-harney/hpc_blast2megan"
cp -r hpc_blast2megan/b2m_scripts .
```

Check the contents of the b2m_scripts directory. There should be 5 files in the directory: 6 .sh files and 1 .R script:

```
ls b2m_scripts
```

<font size="4"><b>2.4) A note on editing scripts</b></font>
<br></br>
Unlike the scripts in the dada2 pipeline, the user does not provide their email address as a command line argument, and by default will not receive email confirmation of job completion. However this can be easily altered through a small change to the resource request section of the .sh scripts. A script can be viewed and edited with the nano command and the relative or absolute path to the script, e.g. :

```
nano b2m_scripts/01_run_blastn_simple.sh
```

This will start nano. Notice that the first line of the script is #!/bin/bash, followed by an empty line, and then several lines commencing with #SBATCH. These #SBATCH arguments are used by slurm when the script is submitted (with qsub or sbatch) and allow the user to control certain parameters relating to job submission. Notice that the last #SBATCH line is:

* #SBATCH --mail-user=user@uni.ac.uk. 

Using the arrow key, go to this line and change user@uni.ac.uk to your own email address. In my case it would be:

* #SBATCH --mail-user=e.harney@sheffield.ac.uk

Once you have made this change, you will need to save it. Notice at the bottom of the screen are lines of commands, such as "^G Get Help" and "^X Exit" etc. The "^" means holding down the control (windows) or command (mac) key. Pressing the "X" key whilst holding down control/command will allow you to save (if there have been changes) and exit. After pressing "^X" you will be prompted to save the changes (options are "y" for yes, "n" for no and "^C" for cancel). Press "y". You will then be given the chance to rename the file if you want. In our case we can keep the old name, so simply press "enter" to save the file with the same name and exit nano.
<br>
</details>
<br>

<details><summary><font size="6"><b>3) Blast sequence data against an ncbi database</font></b></summary>
<br></br>

<font size="4"><b>3.1) Determine how many sequences are in your fasta file</b></font>
<br></br>

Depending on the size and number of sequences in our fasta file and the size of the database being used, this step can be quite slow. If your fasta file contains thousands of sequences, we can speed things up by slitting the fasta file into chunks and running blastn in parallel using the array functionality of slurm. This workflow therefore contains 2 different options for running blastn:

1. If you have < 1000 sequences we suggest running the single script:
   * '01_run_blastn_simple.sh'.
3. If you have > 1000 sequences we suggest splitting the file into chunks running it in array mode with:
   * '01A_run_prep_for_blast.sh'
   * '01B_run_blastn_array.sh'

<br>
To see how many sequences are in your fasta file, run the following:

```
grep -c '>' working_data/06_ASV_seqs.fasta
```

Although we suggest 1000 sequences as a threshold, you can run an array with less than 1000 sequences or the simple blast more than 1000 sequences (although this may take a while). In section 3.2 we describe how to run blast in simple mode, and in section 3.3 we describe how to run it in array mode.
<br></br>

<font size="4"><b>3.2) Running blastn in simple mode</b></font>
<br></br>

Running blastn in simple mode will create a new directory called blast_out in your current directory, as well as symbolic links to the ncbi taxadb files 'taxdb.btd' and 'taxdb.bti'. It will then run blastn and the output will be saved as blast_out/all_blast.out.tab. 
<br></br>

<b>To run 01_run_blastn_simple.sh you need to provide:</b>
* the relative path to the fasta file containing the sequence data (-F)
* the location of an ncbi database on the HPC (-B)
<br></br>

It is most likely that you will use the nt database, which contains all nucleotide sequences available on GenBank. However, you can also supply a smaller or bespoke indexed database. Here we will assume you are using nt. Here's an example for how to submit the job:
  
```
qsub b2m_scripts/01_run_blastn_simple.sh -F working_data/06_ASV_seqs.fasta -B /shared/genomicsdb2/shared/ncbi_nt/current/nt
```

<font size="4"><b>3.3) Running blastn in array mode</b></font>
<br></br>

Running blastn in array mode requires running 2 scripts one after the other: '01A_run_split_fasta.sh' then '01B_run_blastn_array.sh'.
<br></br>

The '01A_run_prep_for_blast.sh' splits the input fasta file into chunks.fa files each containing 100 sequences which are written to a new directory called 'split_fasta'. It also creates symbolic links to the ncbi taxadb files 'taxdb.btd' and 'taxdb.bti', and a directory called 'logs' used by script 01B. Finally it creates a text file 'split_fasta_list_of_X.txt' with the names of all the chunk.fa files to be used in the next step. In your file the 'X' will be the total number of chunk.fa files and is a parameter for script '01B_run_blastn_array.sh'.
<br></br>

<b>To run 01A_run_prep_for_blast.sh  you need to provide:</b>
* the relative path to the fasta file containing the sequence data (-F)
<br></br>

An example command if you have run the dada2 pipeline might be:
  
```
qsub b2m_scripts/01A_run_prep_for_blast.sh -F working_data/06_ASV_seqs.fasta
```
  
The '01B_run_blastn_array.sh' script will then use an array to simultaneously blast multiple chunk.fa files against an ncbi database. This script will create a new directory called blast_out in your current directory and writes the output of blasting each chunk against the database to a seperate chunk.fa_blast.out.tab.
<br></br>

<b>To run 01B_run_blastn_array.sh you need to provide:</b>
* the location of an ncbi database on the HPC (-B)
* the number of input files to be run on the array (-N)
<br></br>

As stated in section 3.2, it is most likely that you will use the database nt. The number -N is contained in the file name of 'split_fasta_list_of_X.txt' (in place of the 'X'). This can be viewed with the following command:
  
```
ls split_fasta/split_fasta*
```
  
Slurm job arrays allow batch jobs to be broken down into parts and run in parallel, but the script and it's submission are somewhat different. For more information on arrays refer to the Sheffield HPC documentation on [advanced job submission](https://docs.hpc.shef.ac.uk/en/latest/hpc/scheduler/advanced_job_submission_and_control.html#gsc.tab=0). 
<br></br>

If our original sequence.fasta file contained 2350 sequuences, it would have been split into 24 chunks, with the txt file named split_fasta_list_of_24.txt. This number, 24, will appear twice when we submit this job, which would be as follows:
  
```
sbatch --array=1-24%10 b2m_scripts/01B_run_blastn_array.sh -B /shared/genomicsdb2/shared/ncbi_nt/current/nt -N 24
```
  
Notice that we use sbatch instead of qsub, and that this is followed by array=1- and then the number specific to our data set. This number also appears at the end of the command following the -N flag. 
  
Another aspect of arrays to notice is the array 'throttle', indicated by the `%10` after `--array=1-24`. The throttle is the number of jobs that will be allowed to run simulataneously. It can take any numeric value but we suggest using 10; using a higher number can impact other users' ability to submit and run jobs on the HPC.
  
Error and output log files for each job of the array will be written to the directory 'logs'
<br></br>

<font size="4"><b>3.4) Monitoring and assessing the result of blastn</b></font>
<br></br>

Running blastn against the nt database can take a while. To follow the status of the job run the following command: 

```
squeue --me
```
  
For more information about the squeue output refer to the Sheffield HPC documentation on [squeue](https://docs.hpc.shef.ac.uk/en/latest/referenceinfo/scheduler/SLURM/Common-commands/squeue.html#gsc.tab=0). squeue will show the status of the job, and in the case of an array, how many of the 'subjobs' have been submitted and how many are still queued.
<br></br>  
If blast was run in simple mode, blast_out should now contain a single file called all_blast.out.tab, and if it was run in array mode, it will contain several chunk.fa_blast.out.tab files. Look at the contents of one of the files with:
  
```
head blast_out/all_blast.out.tab 
```
  
or 

```
head blast_out/chunk0.fa_blast.out.tab
```
  
Information about blast tabular output can be found at the [Metagenomics wiki](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6). The column headers in your files correspond to:

qseqid / *saccver* / pident / length / mismatch / gapopen / qstart / qend / sstart / send / evalue / bitscore / *staxid* / *ssciname* / *scomnames* / *sblastname* / *sskingdoms* / *stitle*

(italics highlight differences to the default output).
<br></br> 

All the rows displayed by head (the top 10) are likely to show results for the same sequence (ASV_1 if following the dada2 pipeline) because the query probably matches many sequences in the nt database. Sometimes the alignment will be much better for one species than any other, allowing taxonomic assignment to species level. But often the sequence will align comparably to multiple sequences. In this case, we need to class the sequence at a higher taxonomic level (e.g. genus or family). We will do in the next step using the MEGAN blast2lca algorithm.
<br>
</details>
<br>

<details><summary><font size="6"><b>4) Run the MEGAN blast2lca algorithm</font></b></summary>
<br></br>
  
<font size="4"><b>4.1) Run blast2lca  </b></font>
<br></br>

For this step we will use [MEGAN](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/), a suite of bioinformatic algorithms. MEGAN contains various tools: we are interested in blast2lca, which calculates the LCA or [lowest common ancestor](https://en.wikipedia.org/wiki/Lowest_common_ancestor) using multiple blast results.
<br></br>

The 02_run_blast2lca.sh script takes the output from blast, and does the following:
<br></br>

1. If blast was run in array mode, chunks are merged (this is automatically detected)
2. Blast results are filtered by pident (percentage of identical positions: column 3), which is provided by the user
3. The megan blast2lca (lowest common ancestor) algorithm determines taxonomic likelihood of a sequence at all taxonomic levels. The user can adjust the sensitivity of this by varying the top percent arguement.
4. Two main output files are generated, a summary file (the lowest common ancestor and its rank) and a taxon path file (all the taxonomic levels to the lowest common ancestor).
<br></br>

This script produces some temporary files (which are deleted), some intermediate files ('filtered_blast.out.tab', 'megan_full_out.tsv' and 'megan_taxonpath_withcounts.tsv') that may be useful, and two final output files: 'megan_summary_out.tsv' and megan_taxonpath_out.tsv.
<br></br>

<b>To run the 02_run_blast2lca.sh script you must provide:</b>
* Minimum percentage identity (0-100) for the blast results to be considered by blast2lca (-B)
* The Top Percent parameter (1-10) for lca calculation (-T)
* Absolute path to the megan nucleotide database (-D)
<br></br>

Blast will potentially output hundred of hits for each ASV. The minimum percentage identity threshold can be used to reduce the number of hits that are considered by MEGAN's blast2lca algorithm. We suggest using a relatively high value to start with (90 or 95), which can be reduced if high numbers of NAs appear in the summary file.
<br></br>

Blast2lca only retains blast hits that are less than the 'Top Percent' value away from the highest scroing hit (based on bit score). Top Percent can actually be set to anything from 0 to 100 (with a default of 10), but in our experience values of 1-10 are most appropriate. Setting Top Percent to a low value (1.5, 2) will retain a relatively smaller number of blast hits and likely result in better taxonomic assignment. However if the ASVs derive from organisms that have poor representation in the reference database then it may be wiser to set Top Percent to a higher value (5-10). For more information about Top Percent and the blast2lca algorithm please see the [MEGAN manual](https://software-ab.cs.uni-tuebingen.de/download/megan6/manual.pdf).
<br></br>

If running the analysis on BESSEMER, the megan nucleotide database should be available at '/shared/genomicsdb2/shared/megan/megan-nucl-Feb2022.db'. Otherwise you can download your own version from the [MEGAN Alternative Download Page](https://unitc-my.sharepoint.com/personal/iijhu01_cloud_uni-tuebingen_de/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fiijhu01%5Fcloud%5Funi%2Dtuebingen%5Fde%2FDocuments%2FApps%2FMegan&ga=1). For deciding which values of -B and -T to use, we recommend initially using relatively strict (high) values for both, such as:
* -B 95
* -T 2
<br></br>

Thus you might run the job like so:

```
qsub b2m_scripts/02_run_blast2lca.sh -B 95 -T 2 -D /shared/genomicsdb2/shared/megan/megan-nucl-Feb2022.db
```

<br></br>

<font size="4"><b>4.2) Tweaking the parameters of the blast2lca script</b></font>
<br><br>

You can look at the results for the first few sequences of the summary file with:

```
head blast_out/megan_summary_out.tsv
```
Columns corespond to: ASV ID / taxonomic rank / LCA taxon. 

To get a rough idea of how well blast and megan have worked, you can run the following code:

```
cut -f2 blast_out/megan_summary_out.tsv | sort | uniq -c 
```
  
This will show how many sequences have been assigned to each taxonomic rank. Hopefully the majority of your sequences are at 's' (species) or 'g' (genus) level. The 'x' level refers to ASVs where the MEGAN blast2lca algorithm was unable to assign taxonomy. Although subspecies or variants are considered (indicated by 'v') most sequences do not include this level of taxonomic information, so we do not expect a high number here
<br><br>

It is common to have some non assigned and unknown taxa in your data, as well as some taxa assigned to higher levels (e.g. class). Deciding what is a 'good' result will depend on many factors including experiment type, sampling strategy, primers used etc. If you are not satisfied, you can try tweaking the values of -B and -T. 
<br><br>

Reducing -B allows sequences with lower blast percentage identity to be considered by MEGAN, and may reduce the number of sequences without assignment.  
Reducing -T will decrease the number of sequences included in the blast2lca, potentially providing more specific taxonomic assignement.
<br><br>

Take a look at the results for the taxon path file with:

```
head blast_out/megan_taxonpath_out.tsv
```
Columns corespond to: ASV ID / domain / kingdom / phylum / class / order / family / genus / species / subspecies 

The taxon path is provided down to the lowest common ancestor. Lower taxonomic levels are then assigned an NA. If taxonomic information does not exist in the database then it will be referred to as an unknown example using the next highest available taxa. For example, molluscs in the family Lymnaeoidea (class gastropoda) do not have information available about their order (although infraclass, superorder and superfamily levels exist). In this case they will be assigned "unknown_gastropod_order". This does not mean taxonomic assignment has failed, just that the taxonomy of the organism may not conform with the standard 8/9 taxonomic ranks.  
<br>
</details>
<br>

<details><summary><font size="6"><b>5) OPTIONAL: Combine with dada2 output to create summary files</font></b></summary>
<br><br>
  
<font size="4"><b>5.1) Run the summary files script  </b></font>
<br><br>
If you have run this pipeline following dada2, you may wish to combine taxonomic assignment results with sequence data and ASV counts from dada2 to create a summary file with all the information. You may also wish to create files which can be used as input in downstream analysis with the popular community analysis R package [phyloseq](https://joey711.github.io/phyloseq/).
<br><br>
Assuming that you have created the 'megan_summary_out.tsv' and 'megan_taxonpath_out.tsv' files in the previous step, and that you have the files '06_ASV_seqs.fasta' and '06_ASV_counts.tsv' in your 'working_data' directory you should be able to run the following script (no arguments need to be supplied): 

```
qsub b2m_scripts/03_run_make_summary_files.sh
```
  
This script will call the R script 03_make_summary_files.R and write several outputs files to the working_data directory:
- ASV_taxa_seq_counts.tsv: complete summary of taxonomic results (lca taxon and taxon path to the lca taxon), sequence, and count results,
- ps_taxamat.tsv : ASV taxonomic results in matrix format for phyloseq, 
- ps_countmat.tsv : ASV counts results in matrix format for phyloseq, 
- ps_phylogeny.rds : phylogenetic tree prepared according to protocol of [Callahan et al. 2016](https://f1000research.com/articles/5-1492/v1), see subsection _Construct the phylogenetic tree_,

</font>
<br>
</details>
<br>
