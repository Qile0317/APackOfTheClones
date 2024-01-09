test_that("the default plot for all methods", { quietly({
  data("combined_pbmc")

  expect_doppelganger(.defaultApotcPlot, vizAPOTC(combined_pbmc))

  expect_doppelganger(.defaultApotcPlot, APOTCPlot(RunAPOTC(combined_pbmc)))

  # do some runs and sanity test that AdjustAPOTC did nothing
  combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "r1")
  combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "r2")

  expect_doppelganger(.defaultApotcPlot, APOTCPlot(AdjustAPOTC(combined_pbmc)))
  expect_doppelganger(
    .defaultApotcPlot, APOTCPlot(AdjustAPOTC(combined_pbmc), run_id = "r1"),
  )
  
})})
