install.packages("testthat")
library('testthat')

source('sample.R')

test_dir('tests', reporter = 'Summary')
