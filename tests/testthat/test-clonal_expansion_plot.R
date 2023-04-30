source("testdata/SCIP.R")

test_that("The main API works, where retain_axis_scales = F", {
  expect_doppelganger(
    "main_plot",
    clonal_expansion_plot(test_pbmc, test_tcr, verbose = F, clone_scale_factor = 0.5, retain_axis_scales = F))
})

test_that("The API can take in an already-integrated seurat object, where retain_axis_scales = F", {
  expect_doppelganger(
    "main_plot",
    clonal_expansion_plot(test_integrated_pbmc, verbose = F, clone_scale_factor = 0.5, retain_axis_scales = F))
})

test_that("retain_scale() works with verbose = TRUE", {
  expect_doppelganger(
    "scale_retained_main_plot",
    suppressMessages(invisible(
      clonal_expansion_plot(
        test_pbmc, test_tcr, verbose = T, clone_scale_factor = 0.5, retain_axis_scales = T
        ))))
})