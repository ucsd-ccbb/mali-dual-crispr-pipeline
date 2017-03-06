#install.packages("testthat")
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue", suppressUpdates=TRUE)

getwd()
gMainWd = getwd()
TEST_SUBDIR = "tests"
getTestsDir=function(){
  file.path(gMainWd, TEST_SUBDIR)
}


# get unit test library
library('testthat')

# get test helpers for this project 
source(file.path(getTestsDir(), "helper_functions.R"))

# remove any state from earlier test runs
cleanUpOldTestState(gMainWd)

# load and run the code under test
source('sharedRfunctions3.R')
source('tableOutputFunctions.R')
source('plotOutputFunctions.R')

# end sourcing
input_filename = getKnownInputsFp(gMainWd, "A549_CV4_counts_w_everything.txt")
setwd(getTestOutputsFp(gMainWd, ""))
source(file.path(gMainWd, 'mainScoring3.R'))
setwd(gMainWd)

# test the code output against expectations
test_dir(getTestsDir(), reporter = 'Summary')