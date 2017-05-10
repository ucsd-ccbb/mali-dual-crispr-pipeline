#install.packages("testthat")
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue", suppressUpdates=TRUE)

# clear the environment: necessary
# because code under test sets a large number of global variables
# ----------
rm(list = ls())
# ----------

# get unit test library
# ----------
library('testthat')
# ----------

# set up path variables
# ----------
gMainDir = "~/Work/Repositories/mali-dual-crispr-pipeline"
gSourceDir = file.path(gMainDir, "dual_crispr/distributed_files/R")
gKnownGoodDirPath = file.path(gMainDir, "dual_crispr/distributed_files/test_data/test_scoring_workspaces")
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
                            header = TRUE,  stringsAsFactors = FALSE)
  output_file_df = read.table(fp2,
                              sep = "\t",
                              header = TRUE,  stringsAsFactors = FALSE)
  expect_equal(gold_file_df, output_file_df, info = fp1)
}

testSizesOfPairedFilesMatch <- function(fp1, fp2) {
  size1 = file.size(fp1)
  size2 = file.size(fp2)
  expect_equal(size1, size2, info = fp1)
}

testAllFilesInDir <- function() {
  knownGoodFileNames <- list.files(gKnownGoodDirPath)
  for (currFileName in knownGoodFileNames) {
    # the scoring dirs are known different for the known-good and under-test outputs
    if (currFileName == "1_preimport_gScoringDir.txt") next

    filePaths = getPairedFilePaths(currFileName)
    if (file.exists(filePaths$fp2)) {
      if (endsWith(currFileName, ".txt") | endsWith(currFileName, ".csv")) {
        testDfsFromPairedFilesMatch(filePaths$fp1, filePaths$fp2)
      } else {
        if (!endsWith(currFileName, ".RData")) {
          testSizesOfPairedFilesMatch(filePaths$fp1, filePaths$fp2)
        }
      }
    } else {
      print(sprintf("Missing output file %s", filePaths$fp2))
    }
  }

  # outputFileNames <- list.files(gScoringDir)
  # if (!identical(knownGoodFileNames, outputFileNames)) {
  #   print(which(knownGoodFileNames != outputFileNames))
  # }
}
# ----------

# actually do tests
# ----------
# remove any state from earlier test runs
cleanUpOldTestState()

# source the code under test, which runs it since it is one long script :-/
source(file.path(gSourceDir, "scoring_main.R"))

# test the code output against expectations
test_dir(gTestScriptDirPath, reporter = 'Summary')
# ----------
