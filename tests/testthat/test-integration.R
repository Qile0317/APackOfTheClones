source("testdata/SCIP.R")

test_that("integrate_tcr() works", {
  capture_output(
    integration_attempt <- suppressMessages(
      dev_integrate_tcr(test_pbmc, test_tcr, TRUE, FALSE)
    )
  )
  expect_identical(
    integration_attempt,
    test_integrated_pbmc
  )
})

test_that("integrate_tcr() works with verbose = FALSE", {
  expect_identical(
    dev_integrate_tcr(test_pbmc, test_tcr, verbose = FALSE, FALSE),
    test_integrated_pbmc
  )
})

test_that("percent_na() works", {
  expect_equal(percent_na(test_integrated_pbmc), 0)
  temp_pbmc <- test_integrated_pbmc
  for (i in 1:3) {temp_pbmc@meta.data[["barcode"]][[i]] <- NA}
  expect_equal(round(percent_na(temp_pbmc)), 4)
})
