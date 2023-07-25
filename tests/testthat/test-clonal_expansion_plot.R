source("testdata/SCIP.R")

test_that("The main API works, where retain_axis_scales = F", {
  expect_doppelganger(
    "main_plot",
    clonal_expansion_plot(
      test_pbmc, test_tcr, verbose = FALSE, clone_scale_factor = 0.5,
      retain_axis_scales = FALSE, add_size_legend = FALSE
    )
  )
})

test_that("API can take an integrated object, retain_axis_scales = F", {
  expect_doppelganger(
    "main_plot",
    clonal_expansion_plot(
      test_integrated_pbmc, verbose = FALSE, clone_scale_factor = 0.5,
      retain_axis_scales = FALSE, add_size_legend = FALSE
    )
  )
})

test_that("API can use tsne and pca reductions, retain_axis_scales = F", {
  expect_doppelganger(
    "main_pca_plot",
    clonal_expansion_plot(
      test_integrated_pbmc, reduction = 'pca', verbose = FALSE,
      clone_scale_factor = 0.5, retain_axis_scales = FALSE,
      add_size_legend = FALSE
    )
  )
  expect_doppelganger(
    "main_tsne_plot",
    clonal_expansion_plot(
      test_integrated_pbmc, reduction = 'tsne', verbose = FALSE,
      clone_scale_factor = 0.5, retain_axis_scales = FALSE,
      add_size_legend = FALSE
    )
  )
})

test_that("retain_scale() works", {
  expect_doppelganger(
    "scale_retained_main_plot",
    suppressMessages(invisible(
      clonal_expansion_plot(
        test_integrated_pbmc, verbose = FALSE, clone_scale_factor = 0.5,
        retain_axis_scales = TRUE, add_size_legend = FALSE
      )
    ))
  )
})

test_that("repulsion and try_place works", {
  expect_doppelganger(
    "try_place_plot",
    clonal_expansion_plot(
      test_integrated_pbmc, verbose = FALSE, repulse = FALSE,
      clone_scale_factor = 0.5, retain_axis_scales = FALSE, try_place = TRUE,
      add_size_legend = TRUE, legend_buffer = 1.5, legend_spacing = 0.25,
      legend_position = "bottom_right", legend_sizes = c(1, 5)
    )
  )
  expect_doppelganger(
    "repulsed_tight_plot",
      clonal_expansion_plot(
        test_integrated_pbmc, verbose = FALSE, repulse = TRUE,
        clone_scale_factor = 0.5, retain_axis_scales = FALSE, try_place = TRUE,
        add_size_legend = TRUE, legend_buffer = 1.5, legend_spacing = 0.25
      )
  )
})

test_that("scrambling works", {
  expect_doppelganger(
    "scramble_plot",
    clonal_expansion_plot(
      test_integrated_pbmc, verbose = FALSE, repulse = FALSE,
      clone_scale_factor = 0.5, retain_axis_scales = FALSE,
      add_size_legend = FALSE, scramble = TRUE, ORDER = FALSE
    )
  )
})

# incredibly weird, try_place = TRUE works here but not on plot_API... WHAT?!?!
