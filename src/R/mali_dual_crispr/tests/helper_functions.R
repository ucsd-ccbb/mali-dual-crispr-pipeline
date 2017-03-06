TEST_OUTPUT_SUBDIR = "test_outputs"
KNOWN_OUTPUT_SUBDIR = "known_outputs"
KNOWN_INPUT_SUBDIR = "known_inputs"

getKnownInputsFp=function(mainWd, fileName){
  file.path(mainWd, TEST_SUBDIR, KNOWN_INPUT_SUBDIR, fileName)
}

getKnownOutputsFp=function(mainWd, fileName){
  file.path(mainWd, TEST_SUBDIR, KNOWN_OUTPUT_SUBDIR, fileName)
}

getTestOutputsFp=function(mainWd, fileName){
  file.path(mainWd, TEST_SUBDIR, TEST_OUTPUT_SUBDIR, fileName)
}

test_dfs_from_files_match<-function(mainWd, fileName){
  gold_file_df = read.table(getKnownOutputsFp(mainWd, fileName),sep="\t",header=TRUE)
  output_file_df = read.table(getTestOutputsFp(mainWd, fileName),sep="\t",header=TRUE)
  expect_that(gold_file_df, equals(output_file_df))  
}


cleanUpOldTestState<-function(mainWd){
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
    currOutputFp = getTestOutputsFp(mainWd, currFileName)
    if (file.exists(currOutputFp)) file.remove(currOutputFp)
  }
}
