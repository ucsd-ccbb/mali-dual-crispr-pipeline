test_that("scoring produces expected files with expected values on A549_CV4_3-14-21-28 test data", {
  f_fp = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_f.txt"
  fc_fp = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_fc.txt"
  p_fp = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_p.txt"
  pi_fp = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr_pi.txt"
  # not testing pdf file, as hard to compare ...
  
  # Test the output files
  test_dfs_from_files_match(f_fp)
  test_dfs_from_files_match(fc_fp)  
  test_dfs_from_files_match(p_fp)
  test_dfs_from_files_match(pi_fp)
})

test_that("scoring produces expected output text with expected values on A549_CV4_3-14-21-28 test data", {
  # Test the intermediate output values: known-good values are from notebook run
  expect_that(1013, equals(sum(gOutputBad1)))
  expect_that(837, equals(sum(gOutputBad2)))
  
  # expected output from (seeded) least-squared fitting
  expected_rls_output = "1 0.3567908 0.02872182 
2 0.3478716 0.0174502 
3 0.3437614 0.01307251 
4 0.3422289 0.01133606 
5 0.3421561 0.009601067 
6 0.3426542 0.007264814 
7 0.343023 0.005148175 
8 0.3429778 0.003888621 
9 0.3425708 0.003377016 
10 0.3419478 0.00324043 
11 0.3412225 0.003240698 
12 0.3404615 0.003271521 
13 0.3397018 0.003289718 
14 0.3389641 0.003263157 
15 0.3382626 0.003165731 
16 0.3376091 0.003006403 
17 0.3370091 0.002798842 
18 0.3364648 0.00255833 
19 0.3359763 0.002310294 
20 0.3355402 0.002076285 
21 0.3351513 0.001865624 
22 0.3348047 0.00168143 
23 0.3344956 0.001523152 
24 0.3342198 0.001387338 
25 0.3339734 0.001267053 
26 0.3337535 0.001158576 
27 0.3335572 0.001060624 
28 0.3333816 0.0009716184 

 1 
1 0.2817858 0.09006099 
2 0.2053455 0.04950015 
3 0.1629278 0.02785198 
4 0.1407763 0.01606945 
5 0.1298797 0.009105295 
6 0.1247385 0.005035238 
7 0.1223317 0.002758792 
8 0.121187 0.001508346 
9 0.1206306 0.0008258677 

 2 
1 0.274394 0.09145054 
2 0.1961574 0.04921628 
3 0.1532366 0.0271459 
4 0.1317579 0.01381723 
5 0.1221282 0.00661624 
6 0.1179432 0.003098632 
7 0.1160847 0.001461648 
8 0.1152247 0.0007148302 "
  
  outputIrlsLogStrings = paste(gIrlsLogStrings, collapse = '\n')
  expect_that(expected_rls_output, equals(outputIrlsLogStrings))  
  
  expect_that(267, equals(sum(gOutputBad1 & gOutputBad2)))  
})
