source("testdata/SCIP.R")

test_that("The main API works, where retain_axis_scales = F", {
  expect_doppelganger(
    "main_plot",
    clonal_expansion_plot(test_pbmc, test_tcr, progbar = F, clone_scale_factor = 0.5, retain_axis_scales = F))
})