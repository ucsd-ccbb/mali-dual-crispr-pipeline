#install.packages("testthat")
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue", suppressUpdates=TRUE)

# imports
library('testthat')

# functions
cleanUpOldTestState<-function(){
  # clean up global variables I created (there are a ton more that aren't being reset, at least until I refactor them out of existence)
  assign("gFdrsList", NULL, envir = .GlobalEnv)
  assign("gIrlsLogStrings", NULL, envir = .GlobalEnv)
  assign("gOutputBad1Bad2", NULL, envir = .GlobalEnv)
  assign("gSumBad1Bad2", NULL, envir = .GlobalEnv)
  
  # clean 
  testOutputFileNames = c("A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_f.txt", "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_fc.txt",
                      "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_p.txt", "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_pi.txt", 
                      "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr.pdf")
  for (currFileName in testOutputFileNames){
    if (file.exists(currFileName)) file.remove(currFileName)
  }
}

test_dfs_from_files_match<-function(gold_standard_dir, file_fp){
  gold_file_df = read.table(paste0(gold_standard_dir, file_fp),sep="\t",header=TRUE)
  output_file_df = read.table(paste0("../",file_fp),sep="\t",header=TRUE)
  expect_that(gold_file_df, equals(output_file_df))  
}

# top-level code
setwd("/Users/Birmingham/Work/Repositories/ccbb_tickets/20161100_mali_crispr_software/src/R/mali_crispr/")
cleanUpOldTestState()
source('mainScoring3.R')
test_dir('tests', reporter = 'Summary')