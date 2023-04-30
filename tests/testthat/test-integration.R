suppressPackageStartupMessages(library(Seurat))
library(data.table)
library(utils)

source("testdata/SCIP.R")

test_that("integrate_tcr() works", {
  expect_identical(
    integrate_tcr(test_pbmc, test_tcr, verbose = F),
    test_integrated_pbmc)
})

test_that("percent_na() works", {
  expect_equal(percent_na(test_integrated_pbmc), 0)
})