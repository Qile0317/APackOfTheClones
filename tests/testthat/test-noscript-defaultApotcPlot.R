test_that("the default plot for all methods", {

  data("combined_pbmc")
  .defaultApotcPlot <- "default_apotcplot"

  expect_doppelganger(.defaultApotcPlot, vizAPOTC(combined_pbmc))
  expect_doppelganger(.defaultApotcPlot, APOTCPlot(RunAPOTC(combined_pbmc)))

  # do some runs and sanity test that AdjustAPOTC did nothing
  combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "r1")
  combined_pbmc <- RunAPOTC(combined_pbmc, run_id = "r2")

  # test AdjustAPOTC

  expect_doppelganger(.defaultApotcPlot, APOTCPlot(AdjustAPOTC(combined_pbmc)))

  expect_doppelganger(
    .defaultApotcPlot,
    combined_pbmc %>%
      AdjustAPOTC(
        run_id = "r2",
        nudge_cluster = 1:17,
        nudge_vector = c(1, 1)
      ) %>%
      AdjustAPOTC(
        run_id = "r2",
        nudge_cluster = 1:17,
        nudge_vector = c(-1, -1)
      ) %>%
      APOTCPlot(run_id = "r2")
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
    .defaultApotcPlot,
    overlayLegend(
      removeLegend(vizAPOTC(combined_pbmc)),
      legend_position = "bottom left"
    )
  )

  # test that overlaying no shared clones works

  expect_doppelganger(
    .defaultApotcPlot, vizAPOTC(combined_pbmc, show_shared = list())
  )
 
})
