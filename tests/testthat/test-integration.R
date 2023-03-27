suppressPackageStartupMessages(library(Seurat))
library(data.table)

source("testdata/SCIP.R")

# something is wrong with the progress bar
#test_that("integrate_tcr() works", {
#  expect_identical(
#    integrate_tcr(test_pbmc, test_tcr),
#    test_integrated_pbmc)
#})

test_that("percent_na() works", {
  expect_equal(percent_na(test_integrated_pbmc), 0)
})
