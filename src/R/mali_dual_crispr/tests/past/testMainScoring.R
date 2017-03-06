test_that("mainScoring.R produces expected files with expected values on test data", {
  f_fp = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_f.txt"
  fc_fp = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_fc.txt"
  p_fp = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_p.txt"
  pi_fp = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_pi.txt"
  # not testing pdf file, as hard to compare ...
  gold_standard_dir = "/Temp/"
  
  # Test the output files
  test_dfs_from_files_match(gold_standard_dir, f_fp)
  test_dfs_from_files_match(gold_standard_dir, fc_fp)  
  test_dfs_from_files_match(gold_standard_dir, p_fp)
  test_dfs_from_files_match(gold_standard_dir, pi_fp)
})

