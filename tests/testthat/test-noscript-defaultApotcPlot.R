quietly_test_that("the default plot for all methods", {
  data("combined_pbmc")
  .defaultApotcPlot <- "default_apotcplot"

  expect_doppelganger(.defaultApotcPlot, vizAPOTC(combined_pbmc))

  expect_doppelganger(.defaultApotcPlot, APOTCPlot(RunAPOTC(combined_pbmc)))

  # do some runs and sanity test that AdjustAPOTC did nothing
  combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "r1")
  combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "r2")

  # test that AdjustAPOTC does nothing with no argument

  expect_doppelganger(.defaultApotcPlot, APOTCPlot(AdjustAPOTC(combined_pbmc)))
  expect_doppelganger(
    .defaultApotcPlot, APOTCPlot(AdjustAPOTC(combined_pbmc), run_id = "r1"),
  )

  # test default legend functions

  expect_doppelganger(
    .defaultApotcPlot, overlayLegend(vizAPOTC(combined_pbmc))
  )

  expect_doppelganger(
    .defaultApotcPlot,
    overlayLegend(vizAPOTC(combined_pbmc, add_size_legend = FALSE))
  )

  expect_doppelganger(
    .defaultApotcPlot, overlayLegend(removeLegend(vizAPOTC(combined_pbmc)))
  )
  
})
