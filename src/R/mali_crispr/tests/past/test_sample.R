# test_that("prepData produces expected values on test data", {
#   testfile = "/Users/Birmingham/Desktop/A549_CV4_counts_w_everything.txt"
#   nt = 4
#   origProbes = prepData(testfile, nt)
#   revisedProbes = prepDataChanged(testfile, nt)
#   expect_that(origProbes, equals(revisedProbes))
#   #expect_that(results$nn, equals(24309))
#   #expect_that(results$n, equals(74))
# })