source("testdata/SCIP.R")

test_that("insert_legend() works, assuming clonal_expansion_plot() works", {
  expect_doppelganger(
    "main_plot_with_legend",
    clonal_expansion_plot(
      test_pbmc, test_tcr, verbose = F, clone_scale_factor = 0.5, retain_axis_scales = F,
      add_size_legend = T, legend_sizes = c(1, 5, 10), legend_buffer = 0))
  
  expect_doppelganger(
    "main_plot_with_smaller_custom_legend",
    clonal_expansion_plot(
      test_pbmc, test_tcr, verbose = F, clone_scale_factor = 0.5, retain_axis_scales = F,
      add_size_legend = T, legend_sizes = c(1, 3), legend_buffer = 0.5))
})
