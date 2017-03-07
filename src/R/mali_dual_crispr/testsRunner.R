#install.packages("testthat")
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue", suppressUpdates=TRUE)

gMainWd = "~/Work/Repositories/ccbb_tickets_2017/ideker-dual-crispr-software/src/R/mali_dual_crispr"
TEST_SUBDIR = "tests"
TEST_OUTPUT_SUBDIR = "test_outputs"
KNOWN_OUTPUT_SUBDIR = "known_outputs"
KNOWN_INPUT_SUBDIR = "known_inputs"

getTestsDir = function() {
  file.path(gMainWd, TEST_SUBDIR)
}

getKnownInputsFp = function(fileName) {
  file.path(gMainWd, TEST_SUBDIR, KNOWN_INPUT_SUBDIR, fileName)
}

getKnownOutputsFp = function(fileName) {
  file.path(gMainWd, TEST_SUBDIR, KNOWN_OUTPUT_SUBDIR, fileName)
}

getTestOutputsFp = function(fileName) {
  file.path(gMainWd, TEST_SUBDIR, TEST_OUTPUT_SUBDIR, fileName)
}

test_dfs_from_files_match <- function(fileName) {
  gold_file_df = read.table(getKnownOutputsFp(fileName),
                            sep = "\t",
                            header = TRUE)
  output_file_df = read.table(getTestOutputsFp(fileName),
                              sep = "\t",
                              header = TRUE)
  expect_that(gold_file_df, equals(output_file_df))
}


cleanUpOldTestState <- function() {
  # clean up global variables I created (there are a ton more that aren't being reset, at least until I refactor them out of existence)
  assign("gIrlsLogStrings", NULL, envir = .GlobalEnv)
  assign("gOutputBad1", NULL, envir = .GlobalEnv)
  assign("gOutputBad2", NULL, envir = .GlobalEnv)
  
  # clean
  testOutputFileNames = c(
    "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_f.txt",
    "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_fc.txt",
    "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_p.txt",
    "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_pi.txt",
    "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr.pdf"
  )
  for (currFileName in testOutputFileNames) {
    currOutputFp = getTestOutputsFp(currFileName)
    if (file.exists(currOutputFp))
      file.remove(currOutputFp)
  }
}


# get unit test library
library('testthat')

# remove any state from earlier test runs
cleanUpOldTestState()

# load the code under test'
setwd(gMainWd)
source('sharedFunctions.R')
source('robustFittingFunctions.R')
source('tableOutputFunctions.R')
source('plotOutputFunctions.R')
source('mainScoring.R')

# set up input parameters
input_filename = getKnownInputsFp("A549_CV4_counts_w_everything.txt")
niter = 2 # this value would be 1000 in a real run, but is as small as possible here to enable fast testing
time = c(3, 14, 20, 28)
project = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr"
ab0 <-
  c(-19. , -18.5, -18.5, -19. , -19. , -19. , -19. , -19.) # abundance threshold choices

# run scoring, putting all outputs into the test outputs directory
setwd(getTestOutputsFp(""))
output = scoreDualCrisprScreen(input_filename, niter, time, project, ab0)
setwd(gMainWd)

# test the code output against expectations
gOutputBad1 = output$rep1LtTwoTimeptsAboveAbundanceThresholdMask
gOutputBad2 = output$rep2LtTwoTimeptsAboveAbundanceThresholdMask
gIrlsLogStrings = output$irlsLogStrings
test_dir(getTestsDir(), reporter = 'Summary')