#install.packages("testthat")
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue", suppressUpdates=TRUE)

gMainWd = "~/Work/Repositories/ccbb_tickets_2017/ideker-dual-crispr-software/src/R/mali_dual_crispr"
TEST_SUBDIR = "tests"
getTestsDir=function(){
  file.path(gMainWd, TEST_SUBDIR)
}

# get unit test library
library('testthat')

# get test helpers for this project 
setwd(gMainWd)
source(file.path(getTestsDir(), "helper_functions.R"))

# remove any state from earlier test runs
cleanUpOldTestState(gMainWd)

# load and run the code under test
source('sharedFunctions.R')
source('robustFittingFunctions.R')
source('tableOutputFunctions.R')
source('plotOutputFunctions.R')

# end sourcing
input_filename = getKnownInputsFp(gMainWd, "A549_CV4_counts_w_everything.txt")
setwd(getTestOutputsFp(gMainWd, ""))
source(file.path(gMainWd, 'mainScoring.R'))
setwd(gMainWd)

# test the code output against expectations
test_dir(getTestsDir(), reporter = 'Summary')