#install.packages("testthat")
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue", suppressUpdates=TRUE)

# clear the environment
# ----------
rm(list = ls())
# ----------

# set up path variables
# ----------
gMainDir = "~/Work/Repositories/mali-dual-crispr-pipeline/dual_crispr"
gSourceDir = file.path(gMainDir, "distributed_files/R")
gKnownGoodDirPath = file.path(gMainDir, "distributed_files/test_data/test_scoring_workspaces")
gTestScriptDirPath = file.path(gMainDir, "tests")
# ----------

# set up input parameters
# ----------
input_filename = '/Users/Birmingham/dual_crispr/test_data/test_set_8/TestSet8_timepoint_counts.txt'
g_time_str = "21,28"
project = "Notebook8Test"
gUseSeed = TRUE
niter = 2
gAbundanceThreshsDf = read.table('/Users/Birmingham/dual_crispr/test_data/test_set_8/TestSet8_abundance_thresholds.txt', row.names = 1, header = TRUE)
gScoringDir = '/Users/Birmingham/dual_crispr/test_outputs/test_scoring_workspaces'
# ----------

# set up helper function for tests
# ----------
cleanUpOldTestState <- function() {
  existingFileNames <- list.files(gScoringDir)
  for (currFileName in existingFileNames) {
    currOutputFp = file.path(gScoringDir, currFileName)
    if (file.exists(currOutputFp))
      print(sprintf("removing existing file '%s'", currFileName))
      file.remove(currOutputFp)
  }
}

getPairedFilePaths <- function(fileName) {
  fp1 = file.path(gKnownGoodDirPath, fileName)
  fp2 = file.path(gScoringDir, fileName)
  return(list(fp1 = fp1, fp2 = fp2))
}

testDfsFromPairedFilesMatch <- function(fp1, fp2) {
  gold_file_df = read.table(fp1,
                            sep = "\t",
                            header = TRUE)
  output_file_df = read.table(fp2,
                              sep = "\t",
                              header = TRUE)
  expect_that(gold_file_df, equals(output_file_df))
}

testSizesOfPairedFilesMatch <- function(fp1, fp2) {
  size1 = file.size(fp1)
  size2 = file.size(fp2)
  expect_that(size1, equals(size2))
}

testAllFilesInDir <- function() {
  existingFileNames <- list.files(gScoringDir)
  for (currFileName in existingFileNames) {
    filePaths = getPairedFilePaths(currFileName)
    if (endsWith(currFileName, ".txt") | endsWith(currFileName, ".csv")) {
      testDfsFromPairedFilesMatch(filePaths$fp1, filePaths$fp2)
    } else {
      testSizesOfPairedFilesMatch(filePaths$fp1, filePaths$fp2)
    }
  }
}


# ----------

# actually do tests
# ----------
# get unit test library
library('testthat')

# remove any state from earlier test runs
cleanUpOldTestState()

# source the code under test, which runs it since it is one long script :-/
source(file.path(gSourceDir, "scoring_main.R"))

# test the code output against expectations
test_dir(gTestScriptDirPath, reporter = 'Summary')
# ----------
