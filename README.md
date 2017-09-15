# Dual CRISPR Screen Analysis Quick-Start Guide
Amanda Birmingham, CCBB, UCSD (abirmingham@ucsd.edu)

## Table of Contents
* Installation
* Library Definition File Set-Up
* Count Pipeline Execution
* Score Pipeline Execution
* Appendix: Configuration File Modifications

## Installation

This software is provided through conda, a cross-platform package manager that performs installation and building of software packages with all their required dependencies, anc can be installed on any linux-64 or osx-64 platform. **Windows installation is not supported at this time.**

Download and run the installation script:

	curl https://raw.githubusercontent.com/ucsd-ccbb/mali-dual-crispr-pipeline/master/install_dual_crispr.sh -o install_dual_crispr.sh

    	bash install_dual_crispr.sh
	
	source ~/.bashrc

   * This may take several minutes, as many software libraries are being installed!


## Library Definition File Set-Up

Information about the library or libraries of dual CRISPR constructs used in the screen is provided to the pipeline in library definition files--a specially formatted tab-delimited text files created by the user.  These are placed in the `~/dual_crispr/library_definitions` directory, where the software automatically discovers them.  Additional details about these files are available on the pipeline's wiki at [https://github.com/ucsd-ccbb/mali-dual-crispr-pipeline/wiki/Library-Definition-Files](https://github.com/ucsd-ccbb/mali-dual-crispr-pipeline/wiki/Library-Definition-Files) .

### Requirements
1. Information about your dual CRISPR construct library, including unique construct ids, their target and constituent probe ids, and their probe sequences
2. A text editor such as TextEdit, vim, or emacs (do not use Word!)

### Steps
1. Create a new text file containing the following text:

		# library_name = TestLib
		# min_trimmed_grna_len = 19
		# max_trimmed_grna_len = 21
		construct_id	target_a_id	probe_a_id	probe_a_seq	target_b_id	probe_b_id	probe_b_seq

2. Replace "TestLib" with a short, easy-to-type name descriptive of your construct library
3. Replace "19" with the minimum allowable length for a detected gRNA sequence after trimming the scaffold sequences off a read; any trimmed read shorter than this length will be excluded from analysis
4. Replace "21" with the maxiumum allowable length for a detected gRNA sequence after trimming the scaffold sequences off a read; any trimmed read longer than this length will be excluded from analysis
5. On the line directly beneath the one starting with "construct_id", paste in the tab-separated columns of information about your library, with each line containing information about one construct

	* Tab-separated columns are the default format of data copied from Excel, making them easy to produce
	* Specifications for the seven required columns are:

		* `construct_id`: the unique name of the construct, specified in the format probe\_a\_id\_\_probe\_b\_id (e.g., `SMARCA4_chr19_11094819__BRD4_chr19_15376361`, or `NonTargetingControlGuideForHuman0352__SETD2_chr3_47142972`)
		* `target_a_id`: the identifier for the probe's target, usually a gene symbol (e.g., `SMARCA4`, `BRD4`, `SETD2`, `NonTargetingControlGuideForHuman0412`, etc.)
		* `probe_a_id`: the identifier for a specific probe, e.g., `SMARCA4_chr19_11094819`, `BRD4_chr19_15376361`, NonTargetingControlGuideForHuman0352`, `SETD2_chr3_47142972`, etc.)
		* `probe_a_seq`: the DNA-alphabet sequence of the gRNA of probe A, e.g. `TTCAGGGGAAGTATTACAAA`, `AAActgcaTAGCAAGTTgA`, etc.  Sequences may include only canonical DNA bases (A, C, G, and T).  While both upper and lower-case letters may be included, the sequence will be converted to all upper-case before use in the pipeline.
		* The specifications for target\_b\_id, probe\_b\_id, and probe\_b\_seq mirror those of target\_a\_id, probe\_a\_id, and probe\_a\_seq given above.

6. Save this file with the name of your choice and a `.txt` extension
7. Place it in your `~/dual_crispr/library_definitions` directory

	* It will now be available for you to reference when calling the counting and scoring pipelines described below

## Count Pipeline Execution

The count pipeline takes in raw fastq or fastq.gz files from the sequencing center.  It trims and length-filters the reads, then identifies and counts which constructs they represent, and creates output count and plot files.  Users can familiarize themselves with its products by running the count pipeline on a tiny sample dataset provided with the software installation, using the command

	count_dual_crispr CountTest TestLib ~/dual_crispr/test_data/test_set_1 ~/dual_crispr/test_outputs

### Requirements

1. Your user name and password for the FTP server on which your fastq data reside
2. The full URL of the folder on the FTP server in which your fastq data reside

### Steps

1. Create a directory for your fastq data and download it there

	* In the following commands, replacing `fastq_dir_name` everywhere with the name of this run (e.g., `160817_D00611_0339_BHWTT2BCXX`)
	* Replace XXXX and YYYY with the user name and password of the FTP server, and ZZZZ with the full URL of the folder in which the fastq data reside
	* Since fastq data are often large, ensure you make your new directory on a drive with adequate space to hold them!

			mkdir fastq_dir

			cd fastq_dir

			wget --user XXXX --password YYYY -nd ftp://ZZZZ/*.fastq.gz

	* Depending on how much data you have, this may take from a few minutes to a few hours!

2. Run the count pipeline script

	* Provide an alphanumeric-only name for your dataset in place of `dataset_name`, and input the recognized library name for the library used in your screen (e.g., "CV4") in place of `library_name`. Specify the complete path to the fastq directory you created above in place of `fastq_dir_path`, and the complete path to the directory in which you want the folder of output files to be created in place of `output_dir_path`

			count_dual_crispr dataset_name library_name fastq_dir_path output_dir_path

3.  Wait for the run to complete

	* This frequently takes several hours for a typical dataset

4.  After the run is complete, find the results directory

	* `cd` to the directory you input above as `output_dir_path`
	* Look for a folder whose name starts with the value you input above as `dataset_name` and ends with a timestamp suffix matching the timeframe of your most recent run (e.g., `MyTestDataset_20160804205646`)

5. At this point, you may continue to run the score pipeline below if desired


## Score Pipeline Execution

The score pipeline takes in counts files, such as those produced by the count pipeline.  It annotates them with the information needed by the scoring code, determines abundance thresholds, and then calculates fitness and pi scores, and creates output score and plot files. Users can familiarize themselves with its products by running the score pipeline on a tiny sample dataset provided with the software installation, using the command

	score_dual_crispr ScoreTest LargerTestLib ~/dual_crispr/test_data/test_set_6a,~/dual_crispr/test_data/test_set_6b 21,28 ~/dual_crispr/test_outputs --test

### Requirements

1. A file (or multiple files) containing the counts for all the samples in the experiment you wish to score, such as the `*_combined_counts.txt` file produced by the count pipeline

### Steps

1. Run the score pipeline script

	* As in the count pipeline, provide an alphanumeric-only name for your dataset in place of `dataset_name`, and input the recognized library name for the library used in your screen (e.g., "CV4") in place of `library_name`.  Replace `counts_fp_or_dir` with the path identified above in step 3; if wish to combine multiple counts files, you may provide multiple paths separated by commas **ONLY** (no spaces!) `day_timepoints_str` is a list, separated by commas **ONLY** only, providing--in order--the days on which timepoints were collected (e.g., 3,14,21,28). Provide the complete path to the directory in which you want the folder of output files to be created in place of `output_dir_path` .

			score_dual_crispr dataset_name library_name counts_fp_or_dir day_timepoints_str output_dir_path

2.  Wait for the run to complete

	* This usually takes between 20 minutes and an hour for a typical dataset


## Appendix: Configuration File Modifications

All user-configurable settings for the software that are not passed in through the command-line arguments are specified in the `config.txt` file locationed in the `~/dual_crispr` directory, including number of processors on which to run the pipeline, number of mismatches allowed in assigning reads to constructs, number of iterations performed when estimating the pi score, and so forth.  Novice users should probably avoid modifying this, as it is very sensitive to mistakes in spelling, capitalization, spacing, etc, but expert users will find that it allows great flexibility in modifying the parameters and behavior of the pipelines.  Additional information on settings specified in the configuration file   is availabe on the project wiki at [https://github.com/ucsd-ccbb/mali-dual-crispr-pipeline/wiki/Configuration-File](https://github.com/ucsd-ccbb/mali-dual-crispr-pipeline/wiki/Configuration-File).
