#install.packages("testthat")
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue", suppressUpdates=TRUE)

library('testthat')

setwd("/Users/Birmingham/Work/Repositories/ccbb_tickets/20161100_mali_crispr_software/src/R/mali_crispr/")
source('mainScoring2.R')

test_dfs_from_files_match<-function(gold_standard_dir, file_fp){
  gold_file_df = read.table(paste0(gold_standard_dir, file_fp),sep="\t",header=TRUE)
  output_file_df = read.table(paste0("../",file_fp),sep="\t",header=TRUE)
  expect_that(gold_file_df, equals(output_file_df))  
}

test_dir('tests', reporter = 'Summary')