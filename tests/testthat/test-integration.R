suppressPackageStartupMessages(library(Seurat))
library(data.table)
library(utils)

source("testdata/SCIP.R")

test_that("integrate_tcr() works", {
  expect_identical(
    suppressMessages(invisible(integrate_tcr(test_pbmc, test_tcr))),
    test_integrated_pbmc)
})

test_that("integrate_tcr() works with verbose = FALSE", {
  expect_identical(
    integrate_tcr(test_pbmc, test_tcr, verbose = FALSE),
    test_integrated_pbmc)
})

test_that("percent_na() works", {
  expect_equal(percent_na(test_integrated_pbmc), 0)
  temp_pbmc <- test_integrated_pbmc
  for (i in 1:3) {temp_pbmc@meta.data[["barcode"]][[i]] <- NA}
  expect_equal(round(percent_na(temp_pbmc)), 4)
})
