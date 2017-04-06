# Dual CRISPR Screen Analysis Quick-Start Guide
Amanda Birmingham, CCBB, UCSD (abirmingham@ucsd.edu)

## Table of Contents
* Installation
* Library File Definition Set-Up
* Count Pipeline Execution
* Score Pipeline Execution

## Installation

### Requirements
1. A new, empty Amazon Linux AMI instance to which you have access
2. The key (.pem) file that gives you access to the AMI
3. The Public DNS value for the instance
4. A BitBucket account (user name and password) that has been granted access to the `ccb_ucsd/mali-dual-crispr-pipeline` repository
5. Your AWS Access Key ID
6. Your AWS Secret Access Key

### Steps
1. From the command line, log into your instance

	* An example command is shown below; of course, the path to the the pem file should be replaced with the path to your pem, and the  *.amazonaws.com should be replaced with the Public DNS value for your AMI:

		ssh -i ~/Keys/abirmingham_oregon.pem ec2-user@ec2-52-42-121-79.us-west-2.compute.amazonaws.com
		screen
	
	* Instructions from AWS are at [https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AccessingInstancesLinux.html](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AccessingInstancesLinux.html)
	* If you receive a message stating 'The authenticity of host ... can't be established' and asking 'Are you sure you want to continue connecting (yes/no)?', enter `yes`.	
	* If you encounter a `Permission denied (publickey)` error, remember that the permissions on your key (.pem) file must be set so that it is not public, e.g. by running `chmod 0400 ~/Keys/abirmingham_oregon.pem`
	* `screen` ensures you will be able to reconnect to the process if you are disconnected at any point; more details of its operation are available at [https://www.linux.com/learn/taking-command-terminal-gnu-screen](https://www.linux.com/learn/taking-command-terminal-gnu-screen).	

2. Download and install the `conda` package manager software
    
		curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda_py3.sh
    	bash miniconda_py3.sh
    	
    * Press the space bar to move through the license agreement
    * Enter `yes` when asked "Do you wish the installer to prepend the Miniconda3 install location to PATH in your /home/ec2-user/.bashrc ?"
    	
3. **Log out and then back in** to complete the `conda` install

		exit	
   		logout
   		
   	* Then re-enter the same ssh command you used in step 1 (Hint: if you press the up-arrow while in the terminal window, it should appear for you!)

4. Run the following commands to prepare the python environment, configure the necessary `conda` channels, install the pipeline software, and set up the test files

		screen
		conda update conda
		conda install python=3.5
		conda config --add channels bioconda
		conda config --add channels r
		conda config --add channels ccbbucsd
		conda install dual_crispr
		conda install bioconductor-qvalue
		set_up_dual_crispr

	* Order matters here so don't change anything :)
	* Enter `y` whenever prompted to proceed
    	* This may take a minute or two, as many packages are being installed    

5. (Optional) Set up the Jupyter Notebook server for remote access
	
	* This will allow you to view and run the pipeline's Jupyter Notebooks through your browser from the AWS instance.  This is not necessary (as the pipeline can be entirely run from the command line) but is sometimes convenient.
	
		# optional: set up jupyter server for remote login
		cd ~/dual_crispr
		bash set_up_jupyter_server.sh
		cd notebooks
		jupyter notebook
		
6. (Optional) Configure the built-in `aws` software to allow transfer of data back and forth from Amazon's `s3` data storage
    
    		aws configure
    	
    * Enter your AWS Access Key ID and AWS Secret Access Key when prompted, and hit Enter to accept the defaults for the additional prompts.    
	
7. Continue to one of the steps below, **OR** exit the instance with these commands

	    	source deactivate
	    	exit
	    	logout
	    
	    
## Count Pipeline Execution

The count pipeline takes in raw fastq or fastq.gz files from the sequencing center.  It trims and length-filters the reads, then identifies and counts which constructs they represent, and creates output count and plot files. 

### Requirements
1. An Amazon Linux AMI instance to which you have access that has been configured with the pipeline software
2. The key (.pem) file that gives you access to the AMI
3. The Public DNS value for the instance
7. Your user name and password for the FTP server on which your fastq data reside
8. The full URL of the folder on the FTP server in which your fastq data reside
9. The `s3` path to which you wish to store the results

### Steps
1. If you are not already logged into your instance, follow Set-Up step 1 above to do so
2. If you are not already in your `conda` environment for the pipeline, follow Set-Up step 7 to activate it
3. Create a directory for your fastq data and download it there

	* In the following commands, replacing `fastq_dir_name` everywhere with the name of this run (e.g., `160817_D00611_0339_BHWTT2BCXX`)
	* Replace XXXX and YYYY with the user name and password of the FTP server, and ZZZZ with the full URL of the folder in which the fastq data reside
		
			cd /data/raw
			mkdir fastq_dir_name
			cd fastq_dir_name   
			wget --user XXXX --password YYYY -nd ftp://ZZZZ/*.fastq.gz
			
	* Depending on how much data you have, this may take from a few minutes to a few hours!
		
4. Run the count pipeline script

	* Specify the name of the fastq directory you created above in place of `fastq_dir_name`.  Provide an alphanumeric-only name for your dataset in place of `dataset_name`, and input the recognized library name for the library used in your screen (e.g., "CV4") in place of `library_name`.

			cd ~/mali-dual-crispr-pipeline/src/python/
			python run_mali_counting.py fastq_dir_name dataset_name library_name

5.  Wait for the run to complete

	* This frequently takes several hours for a typical dataset

5.  After the run is complete, find the results directory

		cd /data/processed/
		
	* Look for a directory with a name whose timestamp suffix matches the timeframe of your most recent run (e.g., `20160801_LN229_CV4_19mer_1mm_py_20160804205646`)

6. Zip the results directory and upload it to `s3`

	* Replace `results_dir_name` everywhere below with the name identified in the last step and `s3_folder_path` with the `s3` path to which you wish to store the results

			zip -r results_dir_name.zip results_dir_name/
			aws s3 cp results_dir_name.zip s3://s3_folder_path/ 

	* At this point, you may continue to run the score pipeline below, or may exit the configured instance and return to it later.  

8. If you are ready to exit, follow Set-Up step 12 to exit your instance
			
			
## Score Pipeline Execution

The score pipeline takes in counts files, such as those produced by the count pipeline.  It annotates them with the information needed by the scoring code, determines abundance thresholds, and then calculates fitness and pi scores, and creates output score and plot files. 

### Requirements
1. An Amazon Linux AMI instance to which you have access that has been configured with the pipeline software
2. The key (.pem) file that gives you access to the AMI
3. The Public DNS value for the instance
4. A file containing the counts for all the samples in the experiment you wish to score, such as the `*_combined_counts.txt` file produced by the count pipeline
9. The `s3` path to which you wish to store the results

### Steps
1. If you are not already logged into your instance, follow Set-Up step 1 above to do so
2. If you are not already in your `conda` environment for the pipeline, follow Set-Up step 7 to activate it
3. If the counts file you plan to use is not already on the instance, download it from `s3` and make a note of the path to it

	* A sample `s3` download command looks like 
	
			aws s3 cp s3://path/to/my/counts_file.txt /data/raw/counts_file.txt
		
	* Alternately, if you are using the output of the count pipeline, you only need the full path to the relevant run's output folder that was identified in Counts Pipeline step 6

4. Run the score pipeline script

	* Provide an alphanumeric-only name for your dataset in place of `dataset_name`, and input the recognized library name for the library used in your screen (e.g., "CV4") in place of `library_name`.  Replace `counts_fp_or_dir` with the path identified above in step 3, and `day_timepoints_str` with a comma-separated string listing--in order--the days on which timepoints were collected (e.g., "3,14,21,28").

			cd ~/mali-dual-crispr-pipeline/src/python/
			python run_mali_scoring.py dataset_name library_name counts_fp_or_dir day_timepoints_str

5.  Wait for the run to complete

	* This usually takes between 20 minutes and an hour for a typical dataset

6. Follow Counts Pipeline steps 6 and 7 to locate, zip, and upload the results
7. Follow Set-Up step 12 to exit your instance
